# Pathogen Atlas Evidence Explorer

This is a static, client-side explorer. Build Phase 4 release tables and Phase 5 coordinate files, then build the static data from the V3 project root:

```bash
python sites/pathogen-atlas-explorer/scripts/build_data.py
python launch_explorer.py
```

Species metadata and fitted coordinates are loaded once. Supporting evidence is stored in organism-keyed JSON shards and fetched only after a species is selected. The explorer uses the local V2 Wikipedia and GBIF mapping tables where exact species names match; new fungal and viral mappings can be added to the same link-table contract.

The page does not calculate embeddings in the browser. If no V3 coordinates are available, it shows instructions to build the release and embedding phases. Plotly is loaded from its pinned CDN URL in `dist/index.html`.
