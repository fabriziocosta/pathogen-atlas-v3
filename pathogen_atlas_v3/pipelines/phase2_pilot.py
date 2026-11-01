"""Run a reproducible, bounded PubMed pilot for bacteria, fungi, and viruses."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

from pathogen_atlas_v3.abstracts.pubmed import acquire_pubmed
from pathogen_atlas_v3.common.config import selected_organism_groups
from pathogen_atlas_v3.common.io import load_yaml


CURATED_NAMES = {
    "bacteria": [
        "Escherichia coli", "Staphylococcus aureus", "Mycobacterium tuberculosis",
        "Pseudomonas aeruginosa", "Salmonella enterica", "Streptococcus pneumoniae",
        "Clostridioides difficile", "Listeria monocytogenes", "Treponema pallidum",
        "Vibrio cholerae",
    ],
    "fungi": [
        "Candida auris", "Blastomyces gilchristii", "Blastomyces helicus",
        "Blastomyces percusus", "Paracoccidioides lutzii", "Paracoccidioides americana",
        "Talaromyces marneffei", "Cryptococcus gattii", "Sporothrix brasiliensis",
        "Coccidioides posadasii",
    ],
    "viruses": [
        "Betacoronavirus pandemicum", "Alphainfluenzavirus influenzae",
        "Orthohepadnavirus hominoidei", "Lentivirus humimdef1", "Orthoebolavirus zairense",
        "Orthoflavivirus denguei", "Morbillivirus hominis", "Lyssavirus rabies",
        "Varicellovirus humanalpha3", "Simplexvirus humanalpha1",
    ],
}


def select_pilot_species(species: pd.DataFrame, aliases: pd.DataFrame, *, per_group: int = 10, seed: int = 42,
                         organism_groups: list[str] | None = None) -> tuple[pd.DataFrame, dict]:
    """Select curated known-pathogen names first, then deterministic rank-stratified taxa."""
    if per_group < 1:
        raise ValueError("per_group must be positive.")
    alias_index: dict[tuple[str, str], set[str]] = {}
    if not aliases.empty:
        mapped = aliases.loc[aliases["organism_id"].notna()]
        for row in mapped[["organism_group", "name", "organism_id"]].dropna().itertuples(index=False):
            alias_index.setdefault((str(row.organism_group), str(row.name).casefold()), set()).add(str(row.organism_id))
    selected_rows = []
    unmatched: dict[str, list[str]] = {}
    groups = organism_groups or ["bacteria", "fungi", "viruses"]
    for group in groups:
        frame = species.loc[species["organism_group"].astype(str).eq(group)].copy()
        frame = frame.loc[~frame["mapping_status"].astype(str).eq("ambiguous_accepted_name")]
        by_name = {str(row.accepted_name).casefold(): row for row in frame.itertuples(index=False)}
        chosen: dict[str, dict] = {}
        missing = []
        for requested in CURATED_NAMES[group]:
            row = by_name.get(requested.casefold())
            if row is None:
                ids = alias_index.get((group, requested.casefold()), set())
                if len(ids) == 1:
                    match = frame.loc[frame["organism_id"].astype(str).eq(next(iter(ids)))].iloc[0]
                    row = match
            if row is None:
                missing.append(requested)
                continue
            oid = str(row.organism_id)
            chosen.setdefault(oid, {"row": row, "requested_name": requested, "selection_reason": "curated_pathogen_candidate"})
            if len(chosen) == per_group:
                break
        unmatched[group] = missing

        rank_key = "family" if group == "viruses" else "phylum"
        def rank_value(path: object) -> str:
            parts = str(path or "").split(";")
            for part in parts:
                if "=" in part:
                    rank, value = part.split("=", 1)
                    if rank.casefold() == rank_key:
                        return value
            # V2 sources may supply an unlabelled semicolon path. Keep the value
            # as a fallback stratum instead of pretending its rank is known.
            return f"unlabelled:{parts[2] if len(parts) > 2 else 'unknown'}"

        remaining = frame.loc[~frame["organism_id"].astype(str).isin(chosen)].copy()
        remaining["_stratum"] = remaining["taxonomy_path"].map(rank_value)
        remaining["_order"] = remaining["organism_id"].astype(str).map(lambda oid: hashlib.sha256(f"{seed}:{oid}".encode()).hexdigest())
        strata = remaining.sort_values("_order").groupby("_stratum", sort=True, dropna=False)
        for _, candidates in strata:
            if len(chosen) >= per_group:
                break
            row = candidates.iloc[0]
            chosen[str(row.organism_id)] = {"row": row, "requested_name": None, "selection_reason": f"rank_stratified_{rank_key}"}
        if len(chosen) < per_group:
            leftovers = remaining.loc[~remaining["organism_id"].astype(str).isin(chosen)].sort_values("_order")
            for _, row in leftovers.iterrows():
                if len(chosen) >= per_group:
                    break
                chosen[str(row.organism_id)] = {"row": row, "requested_name": None, "selection_reason": "deterministic_fill"}

        for item in chosen.values():
            row = item["row"]
            selected_rows.append({"organism_id": str(row.organism_id), "organism_group": group,
                                  "accepted_name": str(row.accepted_name), "taxonomy_path": str(row.taxonomy_path),
                                  "requested_name": item["requested_name"], "selection_reason": item["selection_reason"]})
    selection = pd.DataFrame(selected_rows)
    if selection.empty or not set(groups).issubset(set(selection["organism_group"])):
        raise ValueError(f"Could not select the requested pilot groups: {groups}.")
    return species.loc[species["organism_id"].astype(str).isin(selection["organism_id"])].copy(), {
        "selection_version": "phase2-pilot-v1", "seed": seed, "requested_per_group": per_group,
        "selected_per_group": selection.groupby("organism_group").size().to_dict(),
        "unmatched_curated_candidates": unmatched, "selection": selection.to_dict("records"),
    }


def run_phase2_pilot(config_path: str = "configs/default.yaml", *, per_group: int = 10, seed: int = 42,
                     pilot_id: str = "pilot_v1",
                     max_fetch_batches: int = 8, retmax: int = 500, batch_size: int = 100,
                     max_workers: int = 5, request_delay_seconds: float | None = None,
                     email: str | None = None, api_key: str | None = None, log: bool = True) -> dict:
    """Acquire into a separate pilot cache/release, leaving production Phase 2 outputs untouched."""
    config = load_yaml(config_path)
    taxonomy = Path(config["paths"]["taxonomy_release"])
    species = pd.read_parquet(taxonomy / "species.parquet")
    aliases = pd.read_parquet(taxonomy / "aliases.parquet")
    groups = selected_organism_groups(config)
    pilot_species, manifest = select_pilot_species(species, aliases, per_group=per_group, seed=seed, organism_groups=groups)
    pilot_ids = set(pilot_species["organism_id"].astype(str))
    pilot_aliases = aliases.loc[aliases["organism_id"].astype("string").isin(pilot_ids)].copy()
    if not pilot_id.replace("_", "").replace("-", "").isalnum():
        raise ValueError("pilot_id may contain only letters, numbers, underscores, and hyphens.")
    root = Path(config["paths"]["data_root"]) / "02_pubmed" / pilot_id
    root.mkdir(parents=True, exist_ok=True)
    (root / "selection_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2, default=str) + "\n", encoding="utf-8")
    prior_report_path = root / "pilot_report.json"
    prior_report = json.loads(prior_report_path.read_text(encoding="utf-8")) if prior_report_path.exists() else {}
    report_fetch_limit = prior_report.get("parameters", {}).get("max_fetch_batches", max_fetch_batches) if max_fetch_batches == 0 else max_fetch_batches
    matches, abstracts = acquire_pubmed(
        pilot_species, pilot_aliases, email=email, api_key=api_key,
        cache_dir=root / "cache", release_dir=root / "release", max_workers=max_workers,
        request_delay_seconds=request_delay_seconds, max_calls=None,
        max_fetch_batches=max_fetch_batches, retmax=retmax, batch_size=batch_size,
        progress_every=int(config.get("logging", {}).get("progress_every", 100)), log=log,
    )
    metadata_path = root / "release" / "acquisition_metadata.json"
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    abstract_pmids = set(abstracts.loc[abstracts.get("abstract", pd.Series(dtype=str)).fillna("").astype(str).str.strip().ne(""), "pmid"].astype(str)) if not abstracts.empty else set()
    if not matches.empty:
        matches = matches.copy()
        matches["abstract_available"] = matches["pmid"].astype(str).isin(abstract_pmids)
    group_summaries = []
    for group in groups:
        selected_group = pilot_species.loc[pilot_species["organism_group"].astype(str).eq(group)]
        group_matches = matches.loc[matches["organism_group"].astype(str).eq(group)] if not matches.empty else pd.DataFrame()
        matched_ids = set(group_matches["organism_id"].astype(str)) if not group_matches.empty else set()
        abstract_rows = group_matches.loc[group_matches["abstract_available"].astype(bool)] if not group_matches.empty else pd.DataFrame()
        group_summaries.append({
            "organism_group": group, "taxa_selected": int(manifest["selected_per_group"].get(group, 0)),
            "taxa_with_pmid_hits": len(matched_ids), "taxa_without_pmid_hits": sorted(set(selected_group["organism_id"].astype(str)) - matched_ids),
            "pmid_matches": int(group_matches["pmid"].nunique()) if not group_matches.empty else 0,
            "species_pmid_rows_with_abstract_text": len(abstract_rows),
            "unique_abstract_pmids_with_text": int(abstract_rows["pmid"].nunique()) if not abstract_rows.empty else 0,
        })
    per_group_query_stats = metadata.get("query_summary_by_group", {})
    report = {
        "pilot_version": "phase2-pilot-v1", "pilot_id": pilot_id, "parameters": {"per_group": per_group, "seed": seed,
            "max_fetch_batches": report_fetch_limit, "batch_size": batch_size, "retmax_per_species_query": retmax,
            "max_workers": max_workers, "request_delay_seconds_override": request_delay_seconds},
        "selection_manifest": str(root / "selection_manifest.json"),
        "acquisition_metadata": metadata, "query_summary_by_group": per_group_query_stats,
        "species_group_summary": group_summaries,
        "total_species_pmid_rows": len(matches), "unique_abstracts_cached": len(abstracts),
        "abstracts_with_text": len(abstract_pmids), "output_root": str(root),
    }
    (root / "pilot_report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2, default=str) + "\n", encoding="utf-8")
    if log:
        print(json.dumps(report, ensure_ascii=False, indent=2, default=str))
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="configs/default.yaml")
    parser.add_argument("--per-group", type=int, default=10)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--max-fetch-batches", type=int, default=8)
    args = parser.parse_args()
    run_phase2_pilot(args.config, per_group=args.per_group, seed=args.seed, max_fetch_batches=args.max_fetch_batches)


if __name__ == "__main__":
    main()
