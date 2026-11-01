"""Phase 2: rate-limited PubMed searches and abstract retrieval."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import pandas as pd

from pathogen_atlas_v3.abstracts.pubmed import (
    acquire_pubmed,
    ambiguous_fungal_accepted_names,
    build_name_query,
    filter_pubmed_search_aliases,
)
from pathogen_atlas_v3.common.config import selected_organism_groups
from pathogen_atlas_v3.common.io import atomic_parquet, load_yaml


def phase2_preflight(config_path: str = "configs/default.yaml") -> dict:
    """Summarize selected taxa, cache hits, and lower-bound ETA without network calls."""
    config = load_yaml(config_path)
    groups = selected_organism_groups(config)
    taxonomy_root = Path(config["paths"]["taxonomy_release"])
    species = pd.read_parquet(taxonomy_root / "species.parquet")
    aliases = pd.read_parquet(taxonomy_root / "aliases.parquet")
    settings = config["pubmed"]
    sampling = settings.get("sampling", {})
    retmax = min(10000, max(1, int(settings.get("retmax", 10000))))
    delay = settings.get("request_delay_seconds")
    if delay is None:
        from pathogen_atlas_v3.abstracts.pubmed import _ncbi_key_from_keychain
        import os
        has_key = bool(os.environ.get("NCBI_API_KEY") or os.environ.get("EUTILS_API_KEY") or _ncbi_key_from_keychain())
        delay = 0.20 if has_key else 0.34
    searchable_aliases, excluded_alias_counts = filter_pubmed_search_aliases(
        aliases, excluded_ambiguous_aliases=settings.get("exclude_ambiguous_aliases", []),
    )
    alias_map = searchable_aliases.groupby("organism_id")["name"].apply(lambda values: list(dict.fromkeys(map(str, values.dropna())))).to_dict() if not searchable_aliases.empty else {}
    ambiguous_fungal_names = ambiguous_fungal_accepted_names(species)
    per_group = []
    root = Path(config["paths"].get("pubmed_domain_root", Path(config["paths"]["pubmed_release"]).parent / "domains"))
    for group in groups:
        group_species = species.loc[species["organism_group"].astype(str).eq(group)]
        cache_path = root / group / "cache" / "pubmed_search_results.jsonl"
        cached = {}
        if cache_path.exists():
            for line in cache_path.read_text(encoding="utf-8").splitlines():
                try:
                    item = json.loads(line)
                    cached[str(item["query_id"])] = item
                except (ValueError, KeyError):
                    continue
        incomplete_cached = sum(not item.get("search_complete", True) for item in cached.values())
        segment_cache_path = root / group / "cache" / "pubmed_search_overflow_segments.jsonl"
        cached_overflow_segments = 0
        if segment_cache_path.exists():
            for line in segment_cache_path.read_text(encoding="utf-8").splitlines():
                try:
                    json.loads(line)
                    cached_overflow_segments += 1
                except ValueError:
                    continue
        pending = 0
        species_without_searchable_names = 0
        for row in group_species[["organism_id", "accepted_name"]].drop_duplicates("organism_id").itertuples(index=False):
            raw_names = alias_map.get(row.organism_id, [])
            names = [name for name in raw_names if len(" ".join(str(name).split())) >= 3 and not re.fullmatch(r"\d+(?:\.\d+)?", " ".join(str(name).split()))]
            accepted_name = " ".join(str(row.accepted_name).split())
            if accepted_name.casefold() not in ambiguous_fungal_names and row.accepted_name not in names:
                names.insert(0, row.accepted_name)
            if not names:
                species_without_searchable_names += 1
                continue
            query = build_name_query(names)
            digest = hashlib.sha256(f"{query}\nretmax={retmax}".encode()).hexdigest()
            if f"q_{digest[:20]}" not in cached:
                pending += 1
        per_group.append({
            "organism_group": group,
            "species": len(group_species),
            "cached_queries": len(cached),
            "pending_queries": pending,
            "species_without_searchable_names": species_without_searchable_names,
            "incomplete_cached_queries": incomplete_cached,
            "cached_overflow_segments": cached_overflow_segments,
            "cache_path": str(cache_path),
        })
    pending_total = sum(item["pending_queries"] for item in per_group)
    eta_seconds = pending_total * float(delay)
    return {"status": "ready", "selected_groups": groups, "per_group": per_group, "pending_queries": pending_total,
            "excluded_pubmed_aliases": excluded_alias_counts,
            "request_delay_seconds": float(delay), "estimated_search_floor_minutes": round(eta_seconds / 60, 1),
            "retmax_per_query": retmax,
            "sampling": {"max_pmids_per_species": sampling.get("max_pmids_per_species", 200),
                         "oldest_pmids_per_species": int(sampling.get("oldest_pmids_per_species", 30)),
                         "random_seed": int(sampling.get("random_seed", 42))},
            "estimated_search_note": "Lower-bound new ESearch time; excludes retries, overflow expansion, throttling, and abstract retrieval."}


def run_phase2(config_path: str = "configs/default.yaml", *, max_calls: int | None = None, max_fetch_batches: int | None = None, log: bool = True):
    config = load_yaml(config_path)
    species = pd.read_parquet(config["paths"]["taxonomy_release"] + "/species.parquet")
    aliases = pd.read_parquet(config["paths"]["taxonomy_release"] + "/aliases.parquet")
    settings, paths = config["pubmed"], config["paths"]
    sampling = settings.get("sampling", {})
    groups = selected_organism_groups(config)
    root = Path(paths.get("pubmed_domain_root", Path(paths["pubmed_release"]).parent / "domains"))
    request_limit = max_calls if max_calls is not None else settings.get("max_calls")
    fetch_limit = max_fetch_batches if max_fetch_batches is not None else settings.get("max_fetch_batches")
    progress_every = int(config.get("logging", {}).get("progress_every", 100))
    results = []
    remaining_calls = None if request_limit is None else int(request_limit)
    for group in groups:
        group_species = species.loc[species["organism_group"].astype(str).eq(group)].copy()
        ids = set(group_species["organism_id"].astype(str))
        group_aliases = aliases.loc[aliases["organism_id"].astype(str).isin(ids)].copy()
        group_root = root / group
        pending = next((x["pending_queries"] for x in phase2_preflight(config_path)["per_group"] if x["organism_group"] == group), 0)
        group_limit = None if remaining_calls is None else min(remaining_calls, pending)
        if log:
            print(f"Phase 2 · {group}: {len(group_species):,} taxa; aliases {len(group_aliases):,}; cache {group_root / 'cache'}; release {group_root / 'release'}.")
        group_matches, group_abstracts = acquire_pubmed(
            group_species, group_aliases,
            cache_dir=group_root / "cache", release_dir=group_root / "release",
            max_workers=int(settings.get("max_workers", 10)), request_delay_seconds=settings.get("request_delay_seconds"),
            max_calls=group_limit, max_fetch_batches=fetch_limit,
            retmax=int(settings.get("retmax", 10000)), batch_size=int(settings.get("batch_size", 200)),
            max_pmids_per_species=sampling.get("max_pmids_per_species", 200),
            oldest_pmids_per_species=int(sampling.get("oldest_pmids_per_species", 30)),
            pmid_sample_random_seed=int(sampling.get("random_seed", 42)),
            progress_every=progress_every,
            expand_overflow=group_limit is None,
            excluded_ambiguous_aliases=settings.get("exclude_ambiguous_aliases", []),
            log=log,
        )
        results.append((group_matches, group_abstracts))
        if remaining_calls is not None:
            remaining_calls = max(0, remaining_calls - int(group_limit or 0))
    matches = pd.concat([item[0] for item in results], ignore_index=True, sort=False) if results else pd.DataFrame()
    abstracts = pd.concat([item[1] for item in results], ignore_index=True, sort=False) if results else pd.DataFrame()
    if not matches.empty:
        matches = matches.drop_duplicates(["organism_id", "pmid"], keep="last")
    if not abstracts.empty:
        abstracts = abstracts.drop_duplicates("pmid", keep="last")
    # Refresh selected groups in the combined release while preserving other
    # groups' fetched abstracts. Abstract records are PMID-keyed and do not
    # carry an organism group, so use the previous match table to find the PMIDs
    # still needed by unselected groups.
    release_root = Path(paths["pubmed_release"])
    release_root.mkdir(parents=True, exist_ok=True)
    matches_path = release_root / "species_abstracts.parquet"
    previous_matches = pd.read_parquet(matches_path) if matches_path.exists() else pd.DataFrame()
    preserved_matches = previous_matches
    preserved_pmids = None
    if not previous_matches.empty and "organism_group" in previous_matches:
        preserved_matches = previous_matches.loc[~previous_matches["organism_group"].astype(str).isin(groups)]
        abstract_matches = preserved_matches
        if "selected_for_sampling" in preserved_matches:
            abstract_matches = preserved_matches.loc[preserved_matches["selected_for_sampling"].fillna(False).astype(bool)]
        preserved_pmids = set(abstract_matches["pmid"].astype(str)) if "pmid" in abstract_matches else set()
    combined_matches = pd.concat([preserved_matches, matches], ignore_index=True, sort=False) if not matches.empty or not preserved_matches.empty else pd.DataFrame()
    if not combined_matches.empty:
        combined_matches = combined_matches.drop_duplicates(["organism_id", "pmid"], keep="last")
    atomic_parquet(combined_matches, matches_path)

    abstracts_path = release_root / "abstracts.parquet"
    previous_abstracts = pd.read_parquet(abstracts_path) if abstracts_path.exists() else pd.DataFrame()
    if preserved_pmids is not None and not previous_abstracts.empty and "pmid" in previous_abstracts:
        previous_abstracts = previous_abstracts.loc[previous_abstracts["pmid"].astype(str).isin(preserved_pmids)]
    combined_abstracts = pd.concat([previous_abstracts, abstracts], ignore_index=True, sort=False) if not abstracts.empty or not previous_abstracts.empty else pd.DataFrame()
    if not combined_abstracts.empty:
        combined_abstracts = combined_abstracts.drop_duplicates(["pmid"], keep="last")
    atomic_parquet(combined_abstracts, abstracts_path)
    return matches, abstracts


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument("--config", default="configs/default.yaml"); parser.add_argument("--max-calls", type=int)
    args = parser.parse_args(); matches, abstracts = run_phase2(args.config, max_calls=args.max_calls)
    print(f"Stored {len(matches):,} organism–PMID matches and {len(abstracts):,} abstracts.")


if __name__ == "__main__": main()
