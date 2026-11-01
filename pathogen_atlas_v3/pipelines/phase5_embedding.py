"""Phase 5: fit the shared-core and per-domain IMSTE spaces."""

from __future__ import annotations

import argparse

import pandas as pd

from pathogen_atlas_v3.common.io import load_yaml
from pathogen_atlas_v3.common.config import selected_organism_groups
from pathogen_atlas_v3.embeddings.imste import fit_feature_space


def build_phase5(config_path: str = "configs/default.yaml", *, shared_feature_ids: list[str] | None = None, settings_overrides: dict | None = None):
    config = load_yaml(config_path)
    matrix = pd.read_parquet(config["paths"]["release_root"] + "/species_feature_matrix.parquet")
    groups = selected_organism_groups(config)
    matrix = matrix.loc[matrix["organism_group"].astype(str).isin(groups)].copy()
    settings = {**config.get("embedding", {}), **(settings_overrides or {})}
    ids = shared_feature_ids or ["airborne", "epidemiology", "geography", "zoonosis", "respiratory", "direct_contact", "outbreaks", "patient_age"]
    return fit_feature_space(matrix, shared_feature_ids=ids, output_dir=config["paths"]["embedding_root"], settings=settings)


def main():
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument("--config", default="configs/default.yaml")
    args=parser.parse_args(); print(build_phase5(args.config).groupby("embedding_scope").size().to_dict())


if __name__ == "__main__": main()
