"""Phase 2: rate-limited PubMed searches and abstract retrieval."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import pandas as pd

from pathogen_atlas_v3.abstracts.pubmed import acquire_pubmed
from pathogen_atlas_v3.abstracts.pubmed import build_name_query
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
    retmax = min(10000, max(1, int(settings.get("retmax", 10000))))
    delay = settings.get("request_delay_seconds")
    if delay is None:
        from pathogen_atlas_v3.abstracts.pubmed import _ncbi_key_from_keychain
        import os
        has_key = bool(os.environ.get("NCBI_API_KEY") or os.environ.get("EUTILS_API_KEY") or _ncbi_key_from_keychain())
        delay = 0.20 if has_key else 0.34
    alias_map = aliases.groupby("organism_id")["name"].apply(lambda values: list(dict.fromkeys(map(str, values.dropna())))).to_dict() if not aliases.empty else {}
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
        for row in group_species[["organism_id", "accepted_name"]].drop_duplicates("organism_id").itertuples(index=False):
            raw_names = alias_map.get(row.organism_id, [])
            names = [name for name in raw_names if len(" ".join(str(name).split())) >= 3 and not re.fullmatch(r"\d+(?:\.\d+)?", " ".join(str(name).split()))]
            if row.accepted_name not in names:
                names.insert(0, row.accepted_name)
            query = build_name_query(names)
            digest = hashlib.sha256(f"{query}\nretmax={retmax}".encode()).hexdigest()
            if f"q_{digest[:20]}" not in cached:
                pending += 1
        per_group.append({
            "organism_group": group,
            "species": len(group_species),
            "cached_queries": len(cached),
            "pending_queries": pending,
            "incomplete_cached_queries": incomplete_cached,
            "cached_overflow_segments": cached_overflow_segments,
            "cache_path": str(cache_path),
        })
    pending_total = sum(item["pending_queries"] for item in per_group)
    eta_seconds = pending_total * float(delay)
    return {"status": "ready", "selected_groups": groups, "per_group": per_group, "pending_queries": pending_total,
            "request_delay_seconds": float(delay), "estimated_search_floor_minutes": round(eta_seconds / 60, 1),
            "retmax_per_query": retmax,
            "estimated_search_note": "Lower-bound new ESearch time; excludes retries, overflow expansion, throttling, and abstract retrieval."}


def run_phase2(config_path: str = "configs/default.yaml", *, max_calls: int | None = None, max_fetch_batches: int | None = None, log: bool = True):
    config = load_yaml(config_path)
    species = pd.read_parquet(config["paths"]["taxonomy_release"] + "/species.parquet")
    aliases = pd.read_parquet(config["paths"]["taxonomy_release"] + "/aliases.parquet")
    settings, paths = config["pubmed"], config["paths"]
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
            progress_every=progress_every,
            expand_overflow=group_limit is None,
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
    # Refresh selected groups in the combined release while preserving all other groups.
    release_root = Path(paths["pubmed_release"])
    release_root.mkdir(parents=True, exist_ok=True)
    for filename, incoming, keys in (("species_abstracts.parquet", matches, ["organism_id", "pmid"]), ("abstracts.parquet", abstracts, ["pmid"])):
        target = release_root / filename
        previous = pd.read_parquet(target) if target.exists() else pd.DataFrame()
        if not previous.empty and "organism_group" in previous:
            previous = previous.loc[~previous["organism_group"].astype(str).isin(groups)]
        combined = pd.concat([previous, incoming], ignore_index=True, sort=False) if not incoming.empty or not previous.empty else pd.DataFrame()
        if not combined.empty:
            combined = combined.drop_duplicates(keys, keep="last")
        atomic_parquet(combined, target)
    return matches, abstracts


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument("--config", default="configs/default.yaml"); parser.add_argument("--max-calls", type=int)
    args = parser.parse_args(); matches, abstracts = run_phase2(args.config, max_calls=args.max_calls)
    print(f"Stored {len(matches):,} organism–PMID matches and {len(abstracts):,} abstracts.")


if __name__ == "__main__": main()
