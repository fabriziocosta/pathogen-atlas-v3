"""Seed V3 bacteria decisions from matching V2 records, retaining provenance."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from pathogen_atlas_v3.common.io import atomic_parquet
from pathogen_atlas_v3.decisions.runner import prepare_species_abstracts

LEDGER_KEY = [
    "organism_id", "pmid", "decision_kind", "input_hash", "model",
    "question_version", "profile_version", "run_id", "feature_id",
]
SEED_SOURCE = "v2_compatibility_seed"
V2_GATE_QUESTION_VERSION = "human-pathogen-disease-attribution-v1"
V2_FEATURE_QUESTION_VERSION = "v1-domain-expert-questions-73-v1"
V2_FEATURE_CATALOG_PATH = (
    Path(__file__).resolve().parents[2]
    / "configs/legacy/v2_domain_expert_questions.yaml"
)


def _text(value: Any) -> str:
    if value is None:
        return ""
    try:
        if pd.isna(value):
            return ""
    except (TypeError, ValueError):
        pass
    return str(value).strip()


def _name_signature(value: Any) -> tuple[str, ...]:
    if isinstance(value, (list, tuple, set)) or hasattr(value, "tolist"):
        try:
            values = list(value)
        except TypeError:
            values = [value]
    else:
        values = [value] if _text(value) else []
    return tuple(_text(item) for item in values if _text(item))


def _current_rows(
    abstracts: pd.DataFrame,
    species_matches: pd.DataFrame,
    source_abstracts_path: Path,
) -> pd.DataFrame:
    current = prepare_species_abstracts(abstracts, species_matches)
    current = current.loc[current["organism_group"].astype(str).eq("bacteria")].copy()
    current["pmid"] = current["pmid"].astype(str)
    current["name_signature"] = current["matched_names"].map(_name_signature)

    old_abstracts = pd.read_parquet(source_abstracts_path, columns=["pmid", "title", "abstract"])
    old_abstracts["pmid"] = old_abstracts["pmid"].astype(str)
    old_abstracts = old_abstracts.drop_duplicates("pmid", keep="last").set_index("pmid")
    old_text = old_abstracts.apply(
        lambda row: (_text(row["title"]), _text(row["abstract"])), axis=1
    ).to_dict()
    current["v2_text_signature"] = current["pmid"].map(old_text)
    current["current_text_signature"] = list(
        zip(current["title"].map(_text), current["abstract"].map(_text))
    )
    current["v2_input_unchanged"] = (
        current["v2_text_signature"].notna()
        & current["v2_text_signature"].eq(current["current_text_signature"])
    )
    return current


def _validated_feature_mapping(
    feature_ids: list[str],
    feature_catalog: dict[str, Any],
    legacy_catalog_path: str | Path,
) -> dict[str, str]:
    """Map V2 question IDs to V3 feature IDs after checking exact definitions."""
    legacy = yaml.safe_load(Path(legacy_catalog_path).read_text(encoding="utf-8")) or {}
    if legacy.get("version") != V2_FEATURE_QUESTION_VERSION:
        raise ValueError(
            f"Unexpected V2 feature-question version in {legacy_catalog_path}: "
            f"{legacy.get('version')!r}."
        )
    legacy_questions = legacy.get("questions", {})
    mapping: dict[str, str] = {}
    for feature_id in feature_ids:
        feature = feature_catalog.get(feature_id)
        if not feature:
            raise ValueError(f"V3 bacterial profile refers to unknown feature {feature_id!r}.")
        legacy_id = str(feature.get("legacy_v2_question_id", ""))
        old_question = legacy_questions.get(legacy_id)
        if not legacy_id or old_question is None:
            raise ValueError(f"No V2 question mapping for V3 feature {feature_id!r}.")
        if _text(old_question.get("definition")) != _text(feature.get("definition")):
            raise ValueError(
                f"V2/V3 definition mismatch for feature {feature_id!r}; refusing to seed its score."
            )
        if legacy_id in mapping:
            raise ValueError(f"V2 question ID {legacy_id!r} maps to multiple V3 features.")
        mapping[legacy_id] = str(feature_id)
    return mapping


def _as_v3_rows(
    source: pd.DataFrame,
    *,
    kind: str,
    question_version: str,
    profile_version: str,
    run_id: str,
    source_question_column: str,
    source_profile_column: str,
    source_run_column: str,
    source_path: str,
) -> pd.DataFrame:
    rows = source.copy()
    rows["source_question_version"] = rows[source_question_column].astype(str)
    rows["source_profile_version"] = rows[source_profile_column].astype(str)
    rows["source_run_id"] = rows[source_run_column].astype(str)
    if "source_input_hash" not in rows:
        rows["source_input_hash"] = rows["input_hash"].astype(str)
    rows["accepted_name"] = rows["canonical_species"].astype(str)
    rows["decision_kind"] = kind
    if kind == "pathogen_evidence":
        rows["feature_id"] = "human_pathogen_evidence"
    else:
        source_ids = rows.get("v3_feature_id", rows["question_id"])
        rows["feature_id"] = source_ids.astype(str)
    rows["question_version"] = str(question_version)
    rows["profile_version"] = str(profile_version)
    rows["run_id"] = str(run_id)
    rows["decision_source"] = SEED_SOURCE
    rows["source_project"] = "pathogen_atlas_v2"
    rows["source_record_path"] = source_path
    # The snapshot contains normalized probabilities, not the original raw
    # provider response bodies. Keep the source Parquet path separately.
    rows["raw_response_path"] = pd.NA
    rows["error"] = pd.NA
    columns = [
        "organism_id", "organism_group", "accepted_name", "pmid", "decision_kind",
        "feature_id", "input_hash", "provider", "model", "question_version",
        "profile_version", "run_id", "p_yes", "status", "attempts", "error",
        "completed_at", "raw_response_path", "decision_source", "source_project",
        "source_question_version", "source_profile_version", "source_run_id",
        "source_input_hash", "source_record_path",
    ]
    if kind == "feature":
        rows["dimension_id"] = pd.NA
        rows["dimension"] = pd.NA
        rows["category_id"] = pd.NA
        rows["category"] = pd.NA
        rows["question"] = pd.NA
        rows["question_claim"] = pd.NA
        rows["definition"] = pd.NA
        columns += ["dimension_id", "dimension", "category_id", "category", "question", "question_claim", "definition"]
    return rows[columns]


def seed_v2_bacteria_decisions_as_v3(
    abstracts: pd.DataFrame,
    species_matches: pd.DataFrame,
    *,
    source_root: str | Path,
    ledger_path: str | Path,
    model: str,
    pathogen_question_version: str,
    feature_question_version: str,
    gate_profile_version: str,
    feature_profile_version: str,
    feature_ids: list[str],
    feature_catalog: dict[str, Any],
    legacy_feature_catalog_path: str | Path = V2_FEATURE_CATALOG_PATH,
    pathogen_threshold: float = 0.75,
    log: bool = True,
) -> dict[str, int]:
    """Import V2 scores under V3 cache keys only for unchanged V3 inputs.

    The current question and profile versions are assigned to the imported
    ledger rows by an explicit project-level compatibility decision. Original
    V2 versions, runs, hashes, and source table paths are retained alongside
    them. Abstract/title text and matched aliases must be identical; changed or
    absent inputs are left for live V3 scoring.
    """
    source = Path(source_root)
    if log:
        print("V2 compatibility seed: checking source decisions and matching Phase 2 inputs (local; no Jev calls).")
    feature_id_map = _validated_feature_mapping(
        feature_ids, feature_catalog, legacy_feature_catalog_path
    )
    gate_path = source / "03_decisions/parsed/pathogen_evidence_jev.parquet"
    feature_path = source / "03_decisions/parsed/domain_expert_answers.parquet"
    source_abstracts_path = source / "02_abstracts/release/abstracts.parquet"
    if not all(path.exists() for path in (gate_path, feature_path, source_abstracts_path)):
        if log:
            print(f"V2 decision seed: source snapshot is incomplete under {source}.")
        return {"gate_rows": 0, "feature_rows": 0, "changed_input_pairs": 0}

    ledger = Path(ledger_path)
    if ledger.exists():
        try:
            existing = pd.read_parquet(ledger, columns=["decision_source", "decision_kind"])
        except (KeyError, ValueError):
            existing = pd.DataFrame()
        seeded = existing.loc[existing["decision_source"].astype(str).eq(SEED_SOURCE)] if "decision_source" in existing else pd.DataFrame()
        if not seeded.empty:
            gate_count = int(seeded["decision_kind"].astype(str).eq("pathogen_evidence").sum())
            feature_count = int(seeded["decision_kind"].astype(str).eq("feature").sum())
            if log:
                print(f"V2 decision seed: already initialized ({gate_count:,} gate rows; {feature_count:,} feature rows).")
            return {"gate_rows": gate_count, "feature_rows": feature_count, "changed_input_pairs": 0}

    current = _current_rows(abstracts, species_matches, source_abstracts_path)
    total_pairs = int(current[["organism_id", "pmid"]].drop_duplicates().shape[0])
    current = current.loc[current["v2_input_unchanged"]].copy()
    current = current[[
        "organism_id", "organism_group", "accepted_name", "pmid", "input_hash",
        "matched_names", "name_signature",
    ]].drop_duplicates(["accepted_name", "pmid"])
    current["name_signature"] = current["matched_names"].map(_name_signature)

    source_gate = pd.read_parquet(gate_path)
    source_gate["pmid"] = source_gate["pmid"].astype(str)
    source_gate["canonical_species"] = source_gate["canonical_species"].astype(str)
    source_gate["name_signature"] = source_gate["matched_names"].map(_name_signature)
    source_gate["source_profile_version"] = "not_recorded_in_v2_gate"
    source_gate = source_gate.rename(columns={"input_hash": "source_input_hash"})
    source_gate = source_gate.loc[
        source_gate["status"].astype(str).eq("completed")
        & source_gate["model"].astype(str).eq(str(model))
        & source_gate["question_version"].astype(str).eq(V2_GATE_QUESTION_VERSION)
    ]
    gate = source_gate.merge(
        current,
        left_on=["canonical_species", "pmid", "name_signature"],
        right_on=["accepted_name", "pmid", "name_signature"],
        how="inner",
        validate="one_to_one",
    )
    gate["p_yes"] = pd.to_numeric(gate["p_yes"], errors="coerce")
    gate = gate.loc[gate["p_yes"].between(0, 1)].copy()
    gate_rows = _as_v3_rows(
        gate,
        kind="pathogen_evidence",
        question_version=pathogen_question_version,
        profile_version=gate_profile_version,
        run_id="v3-gate-v1",
        source_question_column="question_version",
        source_profile_column="source_profile_version",
        source_run_column="run_id",
        source_path=str(gate_path),
    )

    passing_pairs = gate.loc[gate["p_yes"].ge(float(pathogen_threshold)), [
        "canonical_species", "pmid", "name_signature",
    ]].drop_duplicates()
    source_features = pd.read_parquet(feature_path)
    source_features["pmid"] = source_features["pmid"].astype(str)
    source_features["canonical_species"] = source_features["canonical_species"].astype(str)
    source_features["name_signature"] = source_features["matched_names"].map(_name_signature)
    source_features = source_features.rename(columns={"input_hash": "source_input_hash"})
    source_features = source_features.loc[
        source_features["status"].astype(str).eq("completed")
        & source_features["model"].astype(str).eq(str(model))
        & source_features["question_set_version"].astype(str).eq(V2_FEATURE_QUESTION_VERSION)
        & source_features["question_id"].astype(str).isin(feature_id_map)
    ]
    source_features["v3_feature_id"] = source_features["question_id"].astype(str).map(feature_id_map)
    features = source_features.merge(
        passing_pairs,
        on=["canonical_species", "pmid", "name_signature"],
        how="inner",
        validate="many_to_one",
    ).merge(
        current,
        left_on=["canonical_species", "pmid", "name_signature"],
        right_on=["accepted_name", "pmid", "name_signature"],
        how="inner",
        validate="many_to_one",
    )
    features["p_yes"] = pd.to_numeric(features["p_yes"], errors="coerce")
    features = features.loc[features["p_yes"].between(0, 1)].copy()
    feature_rows = _as_v3_rows(
        features,
        kind="feature",
        question_version=feature_question_version,
        profile_version=feature_profile_version,
        run_id="v3-features-v1",
        source_question_column="question_set_version",
        source_profile_column="question_set_version",
        source_run_column="expert_run_id",
        source_path=str(feature_path),
    )
    rows = pd.concat([gate_rows, feature_rows], ignore_index=True, sort=False)
    if not rows.empty:
        if ledger.exists():
            existing = pd.read_parquet(ledger)
            combined = pd.concat([rows, existing], ignore_index=True, sort=False)
        else:
            combined = rows
        combined = combined.drop_duplicates(LEDGER_KEY, keep="last")
        atomic_parquet(combined, ledger)

    valid_pairs = int(current[["organism_id", "pmid"]].drop_duplicates().shape[0])
    changed = max(0, total_pairs - valid_pairs)
    if log:
        print(
            f"V2 compatibility seed: initialized {len(gate_rows):,} pathogen-gate rows and "
            f"{len(feature_rows):,} feature rows as V3; {changed:,} bacteria inputs changed or "
            "were absent and will be scored with V3. Original V2 provenance is retained."
        )
    return {"gate_rows": len(gate_rows), "feature_rows": len(feature_rows), "changed_input_pairs": changed}
