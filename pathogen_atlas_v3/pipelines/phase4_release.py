"""Phase 4: deterministic evidence, feature, and Pathogenic Interest releases."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from pathogen_atlas_v3.common.io import atomic_parquet, load_yaml
from pathogen_atlas_v3.decisions.runner import prepare_species_abstracts
from pathogen_atlas_v3.decisions.features import load_feature_catalog
from pathogen_atlas_v3.tables.interest_index import build_pathogenic_interest_index


def build_phase4(config_path: str = "configs/default.yaml", *, decision_path: str | None = None) -> dict[str, int]:
    config = load_yaml(config_path)
    paths, decision = config["paths"], config["decisions"]
    taxonomy_root = Path(paths["taxonomy_release"])
    pubmed_root = Path(paths["pubmed_release"])
    species = pd.read_parquet(taxonomy_root / "species.parquet")
    matches = pd.read_parquet(pubmed_root / "species_abstracts.parquet")
    abstracts = pd.read_parquet(pubmed_root / "abstracts.parquet")
    ledger_path = Path(decision_path or paths["decision_ledger"])
    ledger = pd.read_parquet(ledger_path) if ledger_path.exists() else pd.DataFrame()
    out = Path(paths["release_root"]); out.mkdir(parents=True, exist_ok=True)

    if not matches.empty:
        search_counts = matches.groupby(["organism_id", "organism_group"], as_index=False).agg(n_pubmed_search_records=("pmid", "nunique"))
        available_pmids = set(abstracts.loc[abstracts["abstract"].fillna("").astype(str).str.strip().ne(""), "pmid"].astype(str)) if not abstracts.empty else set()
        fetched = matches.loc[matches["selected_for_fetch"].fillna(False).astype(bool) & matches["pmid"].astype(str).isin(available_pmids)]
        fetched_counts = fetched.groupby(["organism_id", "organism_group"], as_index=False).agg(n_abstracts_fetched=("pmid", "nunique"))
    else:
        search_counts = pd.DataFrame(columns=["organism_id", "organism_group", "n_pubmed_search_records"])
        fetched_counts = pd.DataFrame(columns=["organism_id", "organism_group", "n_abstracts_fetched"])

    current_pairs = prepare_species_abstracts(abstracts, matches)[["organism_id", "pmid", "input_hash"]].drop_duplicates() if not matches.empty else pd.DataFrame(columns=["organism_id", "pmid", "input_hash"])
    gate = ledger.loc[
        ledger["decision_kind"].astype(str).eq("pathogen_evidence")
        & ledger["run_id"].astype(str).eq("v3-gate-v1")
        & ledger["question_version"].astype(str).eq(str(decision.get("pathogen_question_version")))
    ].copy() if not ledger.empty else pd.DataFrame()
    if not gate.empty:
        gate = gate.merge(current_pairs, on=["organism_id", "pmid", "input_hash"], how="inner", validate="many_to_one")
        gate["p_yes"] = pd.to_numeric(gate["p_yes"], errors="coerce")
        gate["candidate_pathogen_evidence"] = gate["status"].astype(str).eq("completed") & gate["p_yes"].ge(float(decision.get("pathogen_threshold", .75)))
        evaluation_approved = bool(decision.get("release_evaluation_approved", False))
        gate["evaluation_approved"] = evaluation_approved
        gate["accepted_pathogen_evidence"] = gate["candidate_pathogen_evidence"] & evaluation_approved
        accepted = gate.loc[gate["accepted_pathogen_evidence"]].drop_duplicates(["organism_id", "pmid"], keep="last")
        pathogen_counts = accepted.groupby(["organism_id", "organism_group"], as_index=False).agg(n_pathogen_pmids=("pmid", "nunique"))
        evidence = gate
        evidence.to_parquet(out / "pathogen_decisions.parquet", index=False)
    else:
        accepted = pd.DataFrame(columns=["organism_id", "organism_group", "pmid", "p_yes", "accepted_pathogen_evidence"])
        pathogen_counts = pd.DataFrame(columns=["organism_id", "organism_group", "n_pathogen_pmids"])
        evidence = pd.DataFrame()

    summary = species.merge(search_counts, on=["organism_id", "organism_group"], how="left").merge(fetched_counts, on=["organism_id", "organism_group"], how="left").merge(pathogen_counts, on=["organism_id", "organism_group"], how="left")
    for col in ("n_pubmed_search_records", "n_abstracts_fetched", "n_pathogen_pmids"):
        summary[col] = pd.to_numeric(summary[col], errors="coerce").fillna(0).astype("int64")
    index = build_pathogenic_interest_index(summary[["organism_id", "organism_group", "n_pathogen_pmids", "n_abstracts_fetched", "n_pubmed_search_records"]])
    summary = summary.merge(index, on=["organism_id", "organism_group", "n_pathogen_pmids", "n_abstracts_fetched", "n_pubmed_search_records"], how="left")

    features = ledger.loc[
        ledger["decision_kind"].astype(str).eq("feature")
        & ledger["run_id"].astype(str).eq("v3-features-v1")
        & ledger["question_version"].astype(str).eq(str(decision.get("feature_question_set_version")))
    ].copy() if not ledger.empty else pd.DataFrame()
    if not features.empty:
        features = features.merge(current_pairs, on=["organism_id", "pmid", "input_hash"], how="inner", validate="many_to_one")
        catalog, _ = load_feature_catalog("configs/features.yaml", "configs/feature_profiles.yaml", allow_provisional=True)
        hierarchy = catalog["features"]
        for column, key in (("dimension_id", "dimension_id"), ("dimension", "dimension_label"),
                            ("category_id", "category_id"), ("category", "category"),
                            ("question", "question"), ("question_claim", "question_claim"),
                            ("definition", "definition")):
            labels = {feature_id: feature.get(key) for feature_id, feature in hierarchy.items()}
            features[column] = features["feature_id"].map(labels).where(features["feature_id"].isin(labels), features.get(column))
        features["p_yes"] = pd.to_numeric(features["p_yes"], errors="coerce")
        features["answer_yes"] = features["status"].astype(str).eq("completed") & features["p_yes"].ge(float(decision.get("feature_yes_threshold", .75)))
        features["scored"] = features["status"].astype(str).eq("completed") & features["p_yes"].notna()
        feature_summary = features.groupby(["organism_id", "organism_group", "feature_id", "dimension_id", "dimension", "category_id", "category", "question", "question_claim", "definition"], dropna=False, as_index=False).agg(n_scored=("scored", "sum"), n_yes=("answer_yes", "sum"), p_yes_mean=("p_yes", "mean"))
        feature_summary["share_yes"] = feature_summary["n_yes"] / feature_summary["n_scored"].replace(0, pd.NA)
        feature_summary["feature_profile_version"] = features["profile_version"].astype(str).iloc[0]
        feature_summary["evaluation_approved"] = bool(decision.get("release_evaluation_approved", False))
        features["evaluation_approved"] = bool(decision.get("release_evaluation_approved", False))
        matrix = feature_summary.pivot(index=["organism_id", "organism_group"], columns="feature_id", values="share_yes").reset_index()
        matrix.columns.name = None
    else:
        feature_summary = pd.DataFrame(columns=["organism_id", "organism_group", "feature_id", "dimension_id", "dimension", "category_id", "category", "question", "question_claim", "definition", "n_scored", "n_yes", "p_yes_mean", "share_yes"])
        matrix = pd.DataFrame(columns=["organism_id", "organism_group"])
    atomic_parquet(summary, out / "species_summary.parquet")
    atomic_parquet(index, out / "pathogenic_interest_index.parquet")
    atomic_parquet(feature_summary, out / "feature_summary.parquet")
    atomic_parquet(features, out / "feature_decisions.parquet")
    atomic_parquet(matrix, out / "species_feature_matrix.parquet")
    atomic_parquet(accepted, out / "accepted_pathogen_evidence.parquet")
    metadata = {"release_version": "v3.1", "pathogen_threshold": float(decision.get("pathogen_threshold", .75)), "feature_yes_threshold": float(decision.get("feature_yes_threshold", .75)), "evaluation_approved": bool(decision.get("release_evaluation_approved", False)), "candidate_species_pmid_records": int(gate["candidate_pathogen_evidence"].sum()) if not gate.empty else 0, "species": len(summary), "accepted_species_pmid_records": len(accepted), "feature_decision_rows": len(features), "index_version": "pathogenic-interest-index-v3.1", "group_normalization": "each organism_group normalized independently; not directly comparable across groups"}
    (out / "release_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    return {"species": len(summary), "accepted_evidence": len(accepted), "feature_decisions": len(features)}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument("--config", default="configs/default.yaml")
    args = parser.parse_args(); print(build_phase4(args.config))


if __name__ == "__main__": main()
