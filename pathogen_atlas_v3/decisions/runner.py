"""Resumable Jev decision stages with grouped feature calls and durable checkpoints."""

from __future__ import annotations

import hashlib
import json
import time
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from datetime import datetime, timezone
from pathlib import Path
from threading import Lock
from typing import Any

import pandas as pd

from pathogen_atlas_v3.common.io import append_jsonl, atomic_parquet
from pathogen_atlas_v3.decisions.jev import DecisionClient, OpenRouterJevClient, extract_noul_probability

GATE_ID = "human_pathogen_evidence"
SUCCESS = "completed"
_raw_lock = Lock()


def _clean(value: Any) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except (TypeError, ValueError):
        pass
    return str(value).strip()


def _hash_input(state: dict[str, Any]) -> str:
    return hashlib.sha256(json.dumps(state, sort_keys=True, ensure_ascii=False, default=str).encode()).hexdigest()


def _decision_state(row: Any) -> dict[str, Any]:
    return {"organism": {"organism_id": row.organism_id, "organism_group": row.organism_group,
                         "accepted_name": row.accepted_name, "matched_names": row.matched_names},
            "article": {"pmid": row.pmid, "title": _clean(row.title), "abstract": _clean(row.abstract)}}


def prepare_species_abstracts(abstracts: pd.DataFrame, species_matches: pd.DataFrame) -> pd.DataFrame:
    """Join selected fetch rows, deduplicate organism/PMID pairs, and hash inputs."""
    if not {"pmid", "title", "abstract"}.issubset(abstracts):
        raise ValueError("Abstracts need pmid, title, and abstract columns.")
    if not {"pmid", "organism_id", "organism_group", "selected_for_fetch"}.issubset(species_matches):
        raise ValueError("Species matches need pmid, organism_id, organism_group, selected_for_fetch.")
    matches = species_matches.loc[species_matches["selected_for_fetch"].fillna(False).astype(bool)].copy()
    matches["pmid"] = matches["pmid"].astype(str).str.strip()
    abstracts = abstracts.copy()
    abstracts["pmid"] = abstracts["pmid"].astype(str).str.strip()
    prepared = matches.merge(abstracts, on="pmid", how="inner", validate="many_to_one")
    names_col = "matched_names" if "matched_names" in prepared else None
    if names_col is None:
        prepared["matched_names"] = [[] for _ in range(len(prepared))]
    prepared["matched_names"] = prepared["matched_names"].map(lambda x: list(x) if isinstance(x, (list, tuple, set)) else ([str(x)] if _clean(x) else []))
    prepared["accepted_name"] = prepared.get("accepted_name", prepared.get("canonical_species"))
    if prepared["accepted_name"].isna().any():
        raise ValueError("Species matches must include accepted_name or canonical_species.")
    prepared = prepared.drop_duplicates(["organism_id", "pmid"], keep="first")
    prepared["input_hash"] = [
        _hash_input({"organism": {"organism_id": row["organism_id"], "organism_group": row["organism_group"],
                                   "accepted_name": row["accepted_name"], "matched_names": row["matched_names"]},
                     "article": {"pmid": row["pmid"], "title": _clean(row["title"]), "abstract": _clean(row["abstract"])}})
        for row in prepared[["organism_id", "organism_group", "accepted_name", "matched_names", "pmid", "title", "abstract"]].to_dict("records")
    ]
    return prepared.sort_values(["organism_group", "accepted_name", "pmid"]).reset_index(drop=True)


def _response_rows(response: dict[str, Any], question_map: dict[str, dict[str, Any]], base: dict[str, Any]) -> list[dict[str, Any]]:
    rows = []
    for qid, question in question_map.items():
        try:
            p_yes = extract_noul_probability(response, qid)
            status, error = SUCCESS, None
        except Exception as exc:
            p_yes, status, error = None, "failed", f"{type(exc).__name__}: {exc}"
        rows.append({**base, "feature_id": qid, "p_yes": p_yes, "status": status, "error": error,
                     "dimension_id": question.get("dimension_id"), "dimension": question.get("dimension_label"),
                     "category_id": question.get("category_id"), "category": question.get("category"),
                     "question": question.get("question"), "question_claim": question.get("question_claim"),
                     "definition": question.get("definition")})
    return rows


def _make_base(row: Any, *, kind: str, run_id: str, question_version: str, profile_version: str, client: DecisionClient, attempts: int = 0) -> dict[str, Any]:
    return {"organism_id": row.organism_id, "organism_group": row.organism_group,
            "accepted_name": row.accepted_name, "pmid": row.pmid, "decision_kind": kind,
            "input_hash": row.input_hash, "provider": getattr(client, "provider", "openrouter_jev"),
            "model": client.model, "question_version": question_version,
            "profile_version": profile_version, "run_id": str(run_id),
            "attempts": attempts, "completed_at": datetime.now(timezone.utc).isoformat()}


def _call_group(row: Any, *, kind: str, questions: dict[str, dict[str, Any]], run_id: str, question_version: str, profile_version: str, client: DecisionClient, raw_path: Path) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    base = _make_base(row, kind=kind, run_id=run_id, question_version=question_version, profile_version=profile_version, client=client)
    state = _decision_state(row)
    try:
        response = client.decide_many(state, questions=questions)
        base["attempts"] = int(getattr(client, "last_attempts", 1))
        status_rows = _response_rows(response, questions, base)
        raw = {**base, "state": state, "questions": list(questions), "response": response,
               "recorded_at": datetime.now(timezone.utc).isoformat()}
    except Exception as exc:
        base["attempts"] = int(getattr(exc, "attempts", 1))
        base["completed_at"] = datetime.now(timezone.utc).isoformat()
        status_rows = [{**base, "feature_id": qid, "p_yes": None, "status": "failed",
                        "error": f"{type(exc).__name__}: {exc}", "dimension_id": q.get("dimension_id"),
                        "dimension": q.get("dimension_label"), "category_id": q.get("category_id"),
                        "category": q.get("category"), "question": q.get("question"),
                        "question_claim": q.get("question_claim"), "definition": q.get("definition")} for qid, q in questions.items()]
        raw = {**base, "state": state, "questions": list(questions), "error": f"{type(exc).__name__}: {exc}", "recorded_at": datetime.now(timezone.utc).isoformat()}
    with _raw_lock:
        append_jsonl(raw_path, raw)
    for result in status_rows:
        result["raw_response_path"] = str(raw_path)
    return status_rows, raw


def run_decision_stage(
    abstracts: pd.DataFrame,
    species_matches: pd.DataFrame,
    *,
    kind: str,
    questions_by_domain: dict[str, dict[str, dict[str, Any]]],
    run_id: str,
    question_version: str,
    profile_version: str,
    ledger_path: str | Path,
    raw_response_path: str | Path,
    max_calls: int | None = None,
    max_workers: int = 32,
    client: DecisionClient | None = None,
    checkpoint_every: int = 25,
    balance_poll_interval_seconds: float = 60.0,
    log: bool = True,
) -> pd.DataFrame:
    """Run one question or grouped question set per organism/PMID and resume safely."""
    if max_workers < 1 or (max_calls is not None and max_calls < 0) or balance_poll_interval_seconds < 0:
        raise ValueError("max_workers must be >= 1, max_calls non-negative or None, and balance polling interval non-negative.")
    prepared = prepare_species_abstracts(abstracts, species_matches)
    if client is None:
        client = OpenRouterJevClient()
    output_path, raw_path = Path(ledger_path), Path(raw_response_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    raw_path.parent.mkdir(parents=True, exist_ok=True)
    ledger = pd.read_parquet(output_path) if output_path.exists() else pd.DataFrame()
    task_questions: dict[tuple[str, str], dict[str, dict[str, Any]]] = {}
    for domain in prepared["organism_group"].drop_duplicates():
        qmap = questions_by_domain.get(str(domain), {})
        if kind == "pathogen_evidence":
            qmap = {GATE_ID: qmap[GATE_ID]} if GATE_ID in qmap else qmap
        if not qmap:
            raise ValueError(f"No {kind} questions configured for {domain}.")
        task_questions[(str(domain), "*")] = qmap

    successful: set[tuple[str, str, str, str, str, str, str, str]] = set()
    if not ledger.empty:
        required_keys = {"organism_id", "pmid", "decision_kind", "input_hash", "model", "question_version", "profile_version", "run_id", "feature_id", "status"}
        if not required_keys.issubset(ledger):
            raise ValueError(f"Decision ledger missing key columns {sorted(required_keys-set(ledger.columns))}.")
        matches = ledger.loc[
            ledger["decision_kind"].astype(str).eq(kind)
            & ledger["run_id"].astype(str).eq(str(run_id))
            & ledger["question_version"].astype(str).eq(str(question_version))
            & ledger["profile_version"].astype(str).eq(str(profile_version))
            & ledger["status"].astype(str).eq(SUCCESS)
        ]
        for x in matches.itertuples(index=False):
            successful.add((str(x.organism_id), str(x.pmid), str(x.input_hash), str(x.model), str(x.question_version), str(x.profile_version), str(x.run_id), str(x.feature_id)))

    def questions_missing(row: Any) -> tuple[dict[str, dict[str, Any]], list[str]]:
        questions = task_questions[(str(row.organism_group), "*")]
        model = str(client.model)
        keys = [(str(row.organism_id), str(row.pmid), str(row.input_hash), model, str(question_version), str(profile_version), str(run_id), str(qid)) for qid in questions]
        return questions, [qid for qid, key in zip(questions, keys) if key not in successful]
    pending_total = sum(bool(questions_missing(row)[1]) for row in prepared.itertuples(index=False))
    total = min(pending_total, int(max_calls)) if max_calls is not None else pending_total
    def pending_tasks():
        submitted = 0
        for row in prepared.itertuples(index=False):
            questions, missing = questions_missing(row)
            if not missing:
                continue
            yield row, {qid: questions[qid] for qid in missing}
            submitted += 1
            if submitted >= total:
                break
    started = time.monotonic()
    balance_method = getattr(client, "remaining_balance", None)
    starting_balance: float | None = None
    balance_label = "remaining balance unavailable"
    latest_balance: float | None = None
    latest_balance_at: float | None = None
    last_balance_poll = 0.0
    credit_burn_per_second: float | None = None
    if log:
        print(f"{kind}: {len(prepared):,} pairs; {len(prepared)-pending_total:,} cached; {total:,} calls scheduled ({pending_total:,} pending); up to {max_workers} concurrent requests.")
        if total and callable(balance_method):
            starting_balance, balance_label = balance_method()
            latest_balance = starting_balance
            latest_balance_at = last_balance_poll = time.monotonic()
            print(f"{kind}: {balance_label} at start: ${starting_balance:.4f}." if starting_balance is not None else f"{kind}: {balance_label}.")
    new_rows: list[dict[str, Any]] = []
    if total:
        tasks = iter(pending_tasks())
        with ThreadPoolExecutor(max_workers=max_workers, thread_name_prefix=f"v3-{kind}") as executor:
            in_flight = set()
            for _ in range(min(total, max_workers * 2)):
                row, qmap = next(tasks)
                in_flight.add(executor.submit(_call_group, row, kind=kind, questions=qmap, run_id=run_id, question_version=question_version, profile_version=profile_version, client=client, raw_path=raw_path))
            completed = 0
            while in_flight:
                done, in_flight = wait(in_flight, return_when=FIRST_COMPLETED)
                for future in done:
                    rows, _raw = future.result()
                    new_rows.extend(rows)
                    completed += 1
                    try:
                        row, qmap = next(tasks)
                        in_flight.add(executor.submit(_call_group, row, kind=kind, questions=qmap, run_id=run_id, question_version=question_version, profile_version=profile_version, client=client, raw_path=raw_path))
                    except StopIteration:
                        pass
                if completed % max(1, checkpoint_every) < len(done) or completed == total:
                    new_frame = pd.DataFrame(new_rows)
                    if ledger.empty:
                        updated = new_frame
                    else:
                        columns = list(dict.fromkeys([*ledger.columns, *new_frame.columns]))
                        all_null = [column for column in columns if (column not in ledger or ledger[column].isna().all()) and (column not in new_frame or new_frame[column].isna().all())]
                        left_null = [column for column in ledger if ledger[column].isna().all()]
                        right_null = [column for column in new_frame if new_frame[column].isna().all()]
                        left = ledger.drop(columns=left_null)
                        right = new_frame.drop(columns=right_null)
                        updated = pd.concat([left, right], ignore_index=True, sort=False)
                        for column in all_null:
                            updated[column] = pd.NA
                    # New rows supersede only rows with the same full logical key.
                    key_cols = ["organism_id", "pmid", "decision_kind", "input_hash", "model", "question_version", "profile_version", "run_id", "feature_id"]
                    updated = updated.drop_duplicates(key_cols, keep="last")
                    atomic_parquet(updated, output_path)
                    ledger = updated
                    if log:
                        elapsed = time.monotonic() - started
                        rate = completed / elapsed if elapsed else 0
                        eta = (total - completed) / rate if rate else 0
                        failed = sum(r["status"] != SUCCESS for r in new_rows)
                        balance_note = ""
                        now = time.monotonic()
                        if callable(balance_method) and now - last_balance_poll >= balance_poll_interval_seconds:
                            observed_balance, observed_label = balance_method()
                            last_balance_poll = now
                            if observed_balance is not None:
                                if latest_balance is not None and latest_balance_at is not None and observed_label == balance_label:
                                    sample_elapsed = now - latest_balance_at
                                    sample_burn = (latest_balance - observed_balance) / sample_elapsed if sample_elapsed > 0 else 0
                                    if sample_burn > 0:
                                        credit_burn_per_second = sample_burn
                                    elif observed_balance > latest_balance:
                                        # A top-up invalidates the previous observed depletion rate.
                                        credit_burn_per_second = None
                                latest_balance = float(observed_balance)
                                latest_balance_at = now
                                balance_label = observed_label
                        if latest_balance is not None:
                            if balance_label == "account credits" and credit_burn_per_second and credit_burn_per_second > 0:
                                credit_eta_hours = latest_balance / credit_burn_per_second / 3600
                                balance_note = f"; credits ${latest_balance:.4f}; credit ETA {credit_eta_hours:.2f} h"
                            elif balance_label == "account credits":
                                balance_note = f"; credits ${latest_balance:.4f}; credit ETA calculating"
                            else:
                                balance_note = f"; {balance_label} ${latest_balance:.4f}"
                        print(f"{kind}: {completed:,}/{total:,} ({completed/total:.1%}); failures {failed:,}; ETA {eta/3600:.2f} h{balance_note}.")
                    new_rows.clear()
    if not output_path.exists():
        atomic_parquet(pd.DataFrame(columns=["organism_id", "organism_group", "pmid", "decision_kind", "feature_id", "input_hash", "provider", "model", "question_version", "profile_version", "run_id", "p_yes", "status", "attempts", "error", "completed_at", "raw_response_path"]), output_path)
    if log and total and callable(balance_method):
        ending_balance, _ = client.remaining_balance()
        if starting_balance is not None and ending_balance is not None:
            print(f"{kind}: remaining {balance_label}: ${ending_balance:.4f}; change this stage: ${starting_balance-ending_balance:.4f}.")
        elif ending_balance is not None:
            print(f"{kind}: remaining {balance_label}: ${ending_balance:.4f}.")
    return pd.read_parquet(output_path)


def run_pathogen_gate(abstracts: pd.DataFrame, species_matches: pd.DataFrame, prompts: dict[str, Any], *, questions_by_domain: dict[str, Any] | None = None, **kwargs: Any) -> pd.DataFrame:
    """Score disease attribution for every selected species–PMID pair."""
    questions = {domain: {GATE_ID: {"type": "noul", **prompts["domains"][domain]}} for domain in prompts["domains"]}
    return run_decision_stage(abstracts, species_matches, kind="pathogen_evidence", questions_by_domain=questions, **kwargs)


def run_domain_features(abstracts: pd.DataFrame, species_matches: pd.DataFrame, gate_decisions: pd.DataFrame, *, threshold: float, questions_by_domain: dict[str, dict[str, dict[str, Any]]], gate_run_id: str, **kwargs: Any) -> pd.DataFrame:
    """Ask all applicable features once per abstract after the pathogen gate passes."""
    prepared = prepare_species_abstracts(abstracts, species_matches)
    gate = gate_decisions.loc[(gate_decisions["run_id"].astype(str) == str(gate_run_id)) & (gate_decisions["decision_kind"].astype(str) == "pathogen_evidence") & (gate_decisions["status"].astype(str) == SUCCESS) & (pd.to_numeric(gate_decisions["p_yes"], errors="coerce") >= float(threshold))]
    key_cols = ["organism_id", "pmid", "input_hash"]
    passing = gate[key_cols].drop_duplicates()
    selected = prepared.merge(passing, on=key_cols, how="inner", validate="one_to_one")
    if selected.empty:
        return pd.DataFrame()
    selected_matches = selected[["organism_id", "organism_group", "pmid", "matched_names", "accepted_name"]].copy()
    selected_matches["selected_for_fetch"] = True
    return run_decision_stage(abstracts, selected_matches, kind="feature", questions_by_domain=questions_by_domain, **kwargs)
