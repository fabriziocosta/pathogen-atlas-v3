"""Phase 3: disease-attribution gate followed by domain-feature scoring."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import yaml

from pathogen_atlas_v3.common.io import load_yaml
from pathogen_atlas_v3.common.config import selected_organism_groups
from pathogen_atlas_v3.decisions.features import load_feature_catalog, questions_for_domain
from pathogen_atlas_v3.decisions.jev import OpenRouterJevClient
from pathogen_atlas_v3.decisions.runner import run_domain_features, run_pathogen_gate
from pathogen_atlas_v3.migration.v2_decision_seed import seed_v2_bacteria_decisions_as_v3


def _seed_v2_compatibility(config, abstracts, matches, groups, catalog, feature_profiles, *, log=True):
    migration = config.get("migration", {})
    if not migration.get("seed_v2_bacteria_decisions_as_v3", False) or "bacteria" not in groups:
        return {"gate_rows": 0, "feature_rows": 0, "changed_input_pairs": 0}
    settings, paths = config["decisions"], config["paths"]
    return seed_v2_bacteria_decisions_as_v3(
        abstracts,
        matches,
        source_root=paths["migration_source"],
        ledger_path=paths["decision_ledger"],
        model=str(settings.get("model", "typesafe/jev-1.14")),
        pathogen_question_version=str(settings["pathogen_question_version"]),
        feature_question_version=str(settings["feature_question_set_version"]),
        gate_profile_version="pathogen-gate-v3.1",
        feature_profile_version=str(catalog["version"]),
        feature_ids=list(feature_profiles.get("bacteria", [])),
        feature_catalog=catalog["features"],
        pathogen_threshold=float(settings.get("pathogen_threshold", .75)),
        log=log,
    )


def phase3_preflight(config_path: str = "configs/default.yaml") -> dict:
    """Prepare compatible local V2 seeds, then estimate pending calls and cost."""
    config = load_yaml(config_path)
    p = Path(config["paths"]["pubmed_release"])
    matches = pd.read_parquet(p / "species_abstracts.parquet")
    abstracts = pd.read_parquet(p / "abstracts.parquet")
    groups = selected_organism_groups(config)
    matches = matches.loc[matches["organism_group"].astype(str).isin(groups)]
    catalog, feature_profiles = load_feature_catalog("configs/features.yaml", "configs/feature_profiles.yaml", allow_provisional=True)
    seed_stats = _seed_v2_compatibility(config, abstracts, matches, groups, catalog, feature_profiles, log=True)
    has_text = abstracts["abstract"].fillna("").astype(str).str.strip().ne("")
    if "title" in abstracts:
        has_text |= abstracts["title"].fillna("").astype(str).str.strip().ne("")
    available = set(abstracts.loc[has_text, "pmid"].astype(str))
    selected = matches.loc[matches["selected_for_fetch"].fillna(False).astype(bool)]
    selected_pairs = selected.loc[selected["pmid"].astype(str).isin(available), ["organism_id", "pmid"]].drop_duplicates()
    gate_pairs = len(selected_pairs)
    decisions = config["decisions"]
    ledger_path = Path(config["paths"]["decision_ledger"])
    cached_gate = 0
    cached_feature = 0
    passing = None
    if ledger_path.exists():
        ledger = pd.read_parquet(ledger_path)
        selected_organisms = set(matches["organism_id"].astype(str))
        gate = ledger.loc[ledger["decision_kind"].astype(str).eq("pathogen_evidence") & ledger["run_id"].astype(str).eq("v3-gate-v1") & ledger["question_version"].astype(str).eq(str(decisions["pathogen_question_version"])) & ledger["status"].astype(str).eq("completed") & ledger["organism_id"].astype(str).isin(selected_organisms)]
        cached_gate = gate[["organism_id", "pmid"]].drop_duplicates().shape[0]
        passing = gate.loc[pd.to_numeric(gate["p_yes"], errors="coerce").ge(float(decisions.get("pathogen_threshold", .75))), ["organism_id", "pmid"]].drop_duplicates()
        feat = ledger.loc[ledger["decision_kind"].astype(str).eq("feature") & ledger["run_id"].astype(str).eq("v3-features-v1") & ledger["question_version"].astype(str).eq(str(decisions["feature_question_set_version"])) & ledger["profile_version"].astype(str).eq(str(catalog["version"])) & ledger["status"].astype(str).eq("completed") & ledger["organism_id"].astype(str).isin(selected_organisms)]
        if passing is not None and not passing.empty and not feat.empty:
            feat = feat.merge(passing, on=["organism_id", "pmid"], how="inner")
            completed_features = feat.groupby(["organism_id", "pmid"], as_index=False)["feature_id"].nunique()
            expected_by_organism = matches[["organism_id", "organism_group"]].drop_duplicates().assign(
                expected=lambda frame: frame["organism_group"].map(lambda group: len(feature_profiles.get(str(group), [])))
            )
            completed_features = completed_features.merge(expected_by_organism, on="organism_id", how="left")
            cached_feature = int(completed_features["feature_id"].ge(completed_features["expected"]).sum())
    pending_gate = max(0, gate_pairs - cached_gate)
    expected_fraction = float(decisions["cost_estimate"].get("expected_pathogen_fraction", .5))
    feature_calls = (int(len(passing)) + int(round(pending_gate * expected_fraction))) if passing is not None else int(gate_pairs * expected_fraction)
    input_tokens = pending_gate * int(decisions["cost_estimate"]["input_tokens_per_gate_call"]) + max(0, feature_calls-cached_feature) * int(decisions["cost_estimate"]["input_tokens_per_feature_call"])
    output_tokens = (pending_gate + max(0, feature_calls-cached_feature)) * int(decisions["cost_estimate"]["output_tokens_per_call"])
    input_price, output_price = decisions["cost_estimate"].get("input_usd_per_million_tokens"), decisions["cost_estimate"].get("output_usd_per_million_tokens")
    estimated_cost = None if input_price is None or output_price is None else input_tokens*float(input_price)/1_000_000 + output_tokens*float(output_price)/1_000_000
    return {"eligible_species_pmid_pairs": gate_pairs, "cached_gate_pairs": cached_gate, "pending_gate_calls": pending_gate,
            "estimated_feature_calls": feature_calls, "cached_feature_pairs": cached_feature,
            "v2_seeded_gate_rows": seed_stats["gate_rows"], "v2_seeded_feature_rows": seed_stats["feature_rows"],
            "v2_changed_input_pairs": seed_stats["changed_input_pairs"],
            "estimated_input_tokens": input_tokens, "estimated_output_tokens": output_tokens, "estimated_cost_usd": estimated_cost,
            "cost_note": "Approximate token assumptions; configure current model prices in configs/default.yaml to calculate USD." if estimated_cost is None else "Estimate uses configured token assumptions and model prices."}


def run_phase3(config_path: str = "configs/default.yaml", *, max_calls: int | None = None, client=None, log: bool = True):
    config = load_yaml(config_path)
    settings, paths = config["decisions"], config["paths"]
    abstracts = pd.read_parquet(Path(paths["pubmed_release"]) / "abstracts.parquet")
    matches = pd.read_parquet(Path(paths["pubmed_release"]) / "species_abstracts.parquet")
    groups = selected_organism_groups(config)
    matches = matches.loc[matches["organism_group"].astype(str).isin(groups)].copy()
    selected_pmids = set(matches["pmid"].astype(str))
    abstracts = abstracts.loc[abstracts["pmid"].astype(str).isin(selected_pmids)].copy()
    if log:
        print(f"Phase 3 selected organism groups: {', '.join(groups)}; {len(matches):,} selected species–PMID rows; {len(abstracts):,} matching abstracts.")
    prompt_path = Path("configs/prompts/pathogen_evidence.yaml")
    prompt = load_yaml(prompt_path)
    catalog, profiles = load_feature_catalog("configs/features.yaml", "configs/feature_profiles.yaml", allow_provisional=True)
    question_maps = {domain: questions_for_domain(catalog, profiles, domain, allow_provisional=True) for domain in config["project"]["domains"]}
    _seed_v2_compatibility(config, abstracts, matches, groups, catalog, profiles, log=log)
    if client is None:
        client = OpenRouterJevClient(model=settings.get("model", "typesafe/jev-1.14"), endpoint=settings.get("endpoint", "https://openrouter.ai/api/alpha/decisions"), timeout=float(settings.get("timeout_seconds", 90)), retries=int(settings.get("retries", 3)), backoff_seconds=float(settings.get("backoff_seconds", 1)))
    limit = max_calls if max_calls is not None else settings.get("max_calls")
    logging_settings = config.get("logging", {})
    log_every = int(logging_settings.get("progress_every", 1000))
    feature_log_every = int(logging_settings.get("feature_progress_every", 100))
    checkpoint_every = int(logging_settings.get("checkpoint_every", 100))
    common = {"run_id": "v3-gate-v1", "question_version": settings["pathogen_question_version"], "profile_version": "pathogen-gate-v3.1", "ledger_path": paths["decision_ledger"], "raw_response_path": paths["raw_responses"], "max_calls": limit, "max_workers": int(settings.get("max_workers", 48)), "checkpoint_every": checkpoint_every, "log_every": log_every, "balance_poll_interval_seconds": float(logging_settings.get("balance_poll_interval_seconds", 60)), "client": client, "log": log}
    gate = run_pathogen_gate(abstracts, matches, prompt, **common)
    feature_common = {**common, "run_id": "v3-features-v1", "question_version": settings["feature_question_set_version"], "profile_version": catalog["version"], "max_calls": limit, "log_every": feature_log_every}
    features = run_domain_features(abstracts, matches, gate, threshold=float(settings.get("pathogen_threshold", .75)), questions_by_domain=question_maps, gate_run_id="v3-gate-v1", gate_question_version=str(settings["pathogen_question_version"]), **feature_common)
    return {"gate": gate, "features": features}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument("--config", default="configs/default.yaml"); parser.add_argument("--max-calls", type=int); parser.add_argument("--preflight", action="store_true")
    args = parser.parse_args()
    if args.preflight:
        print(phase3_preflight(args.config)); return
    result = run_phase3(args.config, max_calls=args.max_calls)
    print(f"Decision ledger rows: {len(result['gate']):,}; feature rows: {len(result['features']):,}.")


if __name__ == "__main__": main()
