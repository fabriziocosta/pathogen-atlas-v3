# Local V3 data

All files in this directory are local inputs, caches, or generated outputs and are ignored by Git.

## Taxonomy source inputs

- `source_taxonomy/fungi/SFPK_Published_Names_Export_16102026.zip`: Species Fungorum snapshot archived 12 October 2026, exported 16 October 2026. Distributed by Royal Botanic Gardens, Kew under CC BY 3.0. SHA-256 is recorded in `configs/migrations/v2_snapshot.yaml`.
- `source_taxonomy/fungi/Species_Fungorum_Snapshot_Mapping.xlsx`: source documentation describing the snapshot columns.
- `source_taxonomy/viruses/ICTV_Master_Species_List_2025_MSL41.v1.xlsx`: ICTV MSL41 v1, 2025 taxonomy ratified February 2026, distributed under CC BY 4.0. SHA-256 is recorded in `configs/migrations/v2_snapshot.yaml`.
- `source_taxonomy/viruses/VMR_MSL41.v1.20260729.xlsx`: matching ICTV Virus Metadata Resource release. It supplies exemplar and additional isolate names, isolate designations, GenBank accessions, genome composition, and host source; the importer joins it to MSL41 taxon IDs and writes `data/01_taxonomy/release/virus_variants.parquet` plus searchable virus-name aliases.
- `migration/version2_snapshot/`: complete local copy of V2's `data/` directory at the source commit recorded in `configs/migrations/v2_snapshot.yaml`.

Source download locations and pinned release identifiers are in `configs/taxonomy/source_files.yaml`.
