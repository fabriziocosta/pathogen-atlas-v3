# V3 notebooks

Run `00_workflow.ipynb` for settings and orientation, then Phases 1–5 in order. Edit `run.organism_groups` in `configs/default.yaml` to select bacteria, fungi, viruses, or any combination for Phases 2–5. Phase 1 builds the full master taxonomy.

After Phases 1 and 4, `06_wikipedia_references.ipynb` and `07_gbif_references.ipynb` create resumable reference-link tables for the selected organism groups. They include only species with a completed positive pathogen assessment. The default selection is currently `viruses`; set it to `bacteria`, `fungi`, or a combination to resolve those groups too. The Explorer build reads these tables by `organism_id`, so same-name taxa in different groups remain distinct.

After Phase 4, `08_pathogenic_interest_index.ipynb` rebuilds and inspects the `pathogenic_interest_index.parquet` release table. It computes the index separately within each organism group and writes the table consumed by the Explorer data builder. Run `python sites/pathogen-atlas-explorer/scripts/build_data.py` afterwards to refresh the web app snapshot.

The notebooks follow the inspectable V2 pattern: show configured inputs and output paths, display a preflight before network work, expose limits and logging intervals, call named Python phase functions, and inspect produced tables. Reusable logic remains in `pathogen_atlas_v3/`.

Phase 2 query and fetch logs include ETA and resume from per-domain caches. Phase 3 preflight estimates the selected calls; Jev logs show concurrency, checkpoint progress, failures, ETA, and remaining OpenRouter balance. Data and generated outputs remain local and outside Git.
