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


def phase3_preflight(config_path: str = "configs/default.yaml") -> dict:
    """Estimate eligible calls, feature calls, tokens and price from local inputs."""
    config = load_yaml(config_path)
    p = Path(config["paths"]["pubmed_release"])
    matches = pd.read_parquet(p / "species_abstracts.parquet")
    abstracts = pd.read_parquet(p / "abstracts.parquet", columns=["pmid", "abstract"])
    groups = selected_organism_groups(config)
    matches = matches.loc[matches["organism_group"].astype(str).isin(groups)]
    available = set(abstracts.loc[abstracts["abstract"].fillna("").astype(str).str.strip().ne(""), "pmid"].astype(str))
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
        gate = ledger.loc[ledger["decision_kind"].astype(str).eq("pathogen_evidence") & ledger["run_id"].astype(str).eq("v3-gate-v1") & ledger["status"].astype(str).eq("completed") & ledger["organism_id"].astype(str).isin(selected_organisms)]
        cached_gate = gate[["organism_id", "pmid"]].drop_duplicates().shape[0]
        passing = gate.loc[pd.to_numeric(gate["p_yes"], errors="coerce").ge(float(decisions.get("pathogen_threshold", .75))), ["organism_id", "pmid"]].drop_duplicates()
        feat = ledger.loc[ledger["decision_kind"].astype(str).eq("feature") & ledger["run_id"].astype(str).eq("v3-features-v1") & ledger["status"].astype(str).eq("completed") & ledger["organism_id"].astype(str).isin(selected_organisms)]
        cached_feature = feat[["organism_id", "pmid"]].drop_duplicates().shape[0]
    feature_calls = int(len(passing)) if passing is not None else int(gate_pairs * float(decisions["cost_estimate"].get("expected_pathogen_fraction", .5)))
    pending_gate = max(0, gate_pairs - cached_gate)
    input_tokens = pending_gate * int(decisions["cost_estimate"]["input_tokens_per_gate_call"]) + max(0, feature_calls-cached_feature) * int(decisions["cost_estimate"]["input_tokens_per_feature_call"])
    output_tokens = (pending_gate + max(0, feature_calls-cached_feature)) * int(decisions["cost_estimate"]["output_tokens_per_call"])
    input_price, output_price = decisions["cost_estimate"].get("input_usd_per_million_tokens"), decisions["cost_estimate"].get("output_usd_per_million_tokens")
    estimated_cost = None if input_price is None or output_price is None else input_tokens*float(input_price)/1_000_000 + output_tokens*float(output_price)/1_000_000
    return {"eligible_species_pmid_pairs": gate_pairs, "cached_gate_pairs": cached_gate, "pending_gate_calls": pending_gate,
            "estimated_feature_calls": feature_calls, "cached_feature_pairs": cached_feature,
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
    if client is None:
        client = OpenRouterJevClient(model=settings.get("model", "typesafe/jev-1.14"), endpoint=settings.get("endpoint", "https://openrouter.ai/api/alpha/decisions"), timeout=float(settings.get("timeout_seconds", 90)), retries=int(settings.get("retries", 3)), backoff_seconds=float(settings.get("backoff_seconds", 1)))
    limit = max_calls if max_calls is not None else settings.get("max_calls")
    progress_every = int(config.get("logging", {}).get("progress_every", 100))
    common = {"run_id": "v3-gate-v1", "question_version": settings["pathogen_question_version"], "profile_version": "pathogen-gate-v3.1", "ledger_path": paths["decision_ledger"], "raw_response_path": paths["raw_responses"], "max_calls": limit, "max_workers": int(settings.get("max_workers", 32)), "checkpoint_every": progress_every, "client": client, "log": log}
    gate = run_pathogen_gate(abstracts, matches, prompt, **common)
    feature_common = {**common, "run_id": "v3-features-v1", "question_version": settings["feature_question_set_version"], "profile_version": catalog["version"], "max_calls": limit}
    features = run_domain_features(abstracts, matches, gate, threshold=float(settings.get("pathogen_threshold", .75)), questions_by_domain=question_maps, gate_run_id="v3-gate-v1", **feature_common)
    return {"gate": gate, "features": features}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument("--config", default="configs/default.yaml"); parser.add_argument("--max-calls", type=int); parser.add_argument("--preflight", action="store_true")
    args = parser.parse_args()
    if args.preflight:
        print(phase3_preflight(args.config)); return
    result = run_phase3(args.config, max_calls=args.max_calls)
    print(f"Decision ledger rows: {len(result['gate']):,}; feature rows: {len(result['features']):,}.")


if __name__ == "__main__": main()
