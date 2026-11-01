"""Phase 2: rate-limited PubMed searches and abstract retrieval."""

from __future__ import annotations

import argparse

import pandas as pd

from pathogen_atlas_v3.abstracts.pubmed import acquire_pubmed
from pathogen_atlas_v3.common.io import load_yaml


def run_phase2(config_path: str = "configs/default.yaml", *, max_calls: int | None = None, max_fetch_batches: int | None = None, log: bool = True):
    config = load_yaml(config_path)
    species = pd.read_parquet(config["paths"]["taxonomy_release"] + "/species.parquet")
    aliases = pd.read_parquet(config["paths"]["taxonomy_release"] + "/aliases.parquet")
    settings = config["pubmed"]
    return acquire_pubmed(species, aliases, email=None, api_key=None,
                          cache_dir=config["paths"]["pubmed_cache"], release_dir=config["paths"]["pubmed_release"],
                          max_workers=int(settings.get("max_workers", 10)), request_delay_seconds=settings.get("request_delay_seconds"),
                          max_calls=max_calls if max_calls is not None else settings.get("max_calls"),
                          max_fetch_batches=max_fetch_batches if max_fetch_batches is not None else settings.get("max_fetch_batches"),
                          batch_size=int(settings.get("batch_size", 200)), log=log)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__); parser.add_argument("--config", default="configs/default.yaml"); parser.add_argument("--max-calls", type=int)
    args = parser.parse_args(); matches, abstracts = run_phase2(args.config, max_calls=args.max_calls)
    print(f"Stored {len(matches):,} organism–PMID matches and {len(abstracts):,} abstracts.")


if __name__ == "__main__": main()
