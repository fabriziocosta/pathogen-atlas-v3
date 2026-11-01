"""Phase 1: assemble taxonomy and aliases from pinned domain sources."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from pathogen_atlas_v3.common.io import atomic_parquet, load_yaml
from pathogen_atlas_v3.migration.version2 import migrate_v2_snapshot
from pathogen_atlas_v3.taxonomy.adapters import ADAPTERS, load_ictv_vmr_export


def run_phase1(config_path: str | Path = "configs/taxonomy/source_files.yaml", *, migration_source: str | Path = "data/migration/version2_snapshot", output_root: str | Path = ".") -> dict[str, int]:
    output_root = Path(output_root)
    migration_source = Path(migration_source)
    counts: dict[str, int] = {}
    # The complete V2 data snapshot is the starting point for the bacterial cohort.
    manifest = migrate_v2_snapshot(migration_source, output_root)
    counts.update({f"migrated_{key}": int(value) for key, value in manifest["counts"].items()})
    config = load_yaml(config_path)
    species_frames = [pd.read_parquet(output_root / "data/01_taxonomy/release/species.parquet")]
    alias_frames = [pd.read_parquet(output_root / "data/01_taxonomy/release/aliases.parquet")]
    variant_frames: list[pd.DataFrame] = []
    for group, settings in config.get("sources", {}).items():
        if group == "bacteria":
            continue
        source_path = Path(settings["path"])
        if not source_path.exists():
            counts[f"{group}_source_pending"] = 1
            continue
        if settings.get("release", "set-before-import") == "set-before-import":
            raise ValueError(f"Pin a source release identifier for {group} in {config_path}.")
        species, aliases = ADAPTERS[group](source_path, release=str(settings["release"]))
        if group == "viruses" and settings.get("variant_metadata_path"):
            variant_path = Path(settings["variant_metadata_path"])
            if variant_path.exists():
                variants, variant_aliases = load_ictv_vmr_export(variant_path, species, release=str(settings["variant_metadata_release"]))
                variant_frames.append(variants)
                alias_frames.append(variant_aliases)
                counts["virus_variant_records"] = len(variants)
                counts["virus_variant_aliases"] = len(variant_aliases)
            else:
                counts["viruses_variant_source_pending"] = 1
        species_frames.append(species); alias_frames.append(aliases)
        counts[f"{group}_species"] = len(species)
        counts[f"{group}_aliases"] = len(aliases)
    all_species = pd.concat(species_frames, ignore_index=True, sort=False).drop_duplicates("organism_id")
    all_aliases = pd.concat(alias_frames, ignore_index=True, sort=False).drop_duplicates(["organism_id", "name", "source"])
    for column in ("organism_id", "organism_group", "accepted_name", "source_taxon_id", "taxonomy_source", "taxonomy_release", "rank", "parent_taxon_id", "taxonomy_path", "source_record_json", "mapping_status"):
        if column in all_species:
            all_species[column] = all_species[column].astype("string")
    for column in ("organism_id", "organism_group", "name", "name_type", "source", "source_record_id", "source_record_json", "mapping_status"):
        if column in all_aliases:
            all_aliases[column] = all_aliases[column].astype("string")
    collisions = all_species["accepted_name"].str.casefold().duplicated(keep=False)
    all_species.loc[collisions, "mapping_status"] = "ambiguous_accepted_name"
    atomic_parquet(all_species, output_root / "data/01_taxonomy/release/species.parquet")
    atomic_parquet(all_aliases, output_root / "data/01_taxonomy/release/aliases.parquet")
    variants_out = pd.concat(variant_frames, ignore_index=True, sort=False) if variant_frames else pd.DataFrame(columns=["organism_id", "organism_group", "source_taxon_id", "accepted_name", "taxonomy_source", "taxonomy_release", "isolate_id", "exemplar_status", "virus_name", "virus_name_abbreviation", "isolate_designation", "genbank_accession", "genome_coverage", "genome", "host_source", "source_record_json", "mapping_status"])
    atomic_parquet(variants_out, output_root / "data/01_taxonomy/release/virus_variants.parquet")
    counts.update({"species_total": len(all_species), "aliases_total": len(all_aliases), "unresolved_alias_mappings": int(all_aliases["organism_id"].isna().sum()), "unresolved_virus_variant_records": int(variants_out["organism_id"].isna().sum()), "accepted_name_collision_records": int(collisions.sum())})
    return counts


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="configs/taxonomy/source_files.yaml")
    parser.add_argument("--migration-source", default="data/migration/version2_snapshot")
    args = parser.parse_args()
    print(run_phase1(args.config, migration_source=args.migration_source))


if __name__ == "__main__":
    main()
