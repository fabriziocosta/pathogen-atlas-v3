# V3 workflow

V3 is an independent Python project. Run commands from the `version3/` root. The notebooks show configuration, input tables, source functions, and resulting output paths before calling a phase. Phase implementations remain in `pathogen_atlas_v3/` so a notebook does not conceal reusable logic.

## Install and configure

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[notebooks,dev]'
```

Set `NCBI_EMAIL` for live PubMed acquisition. Set `NCBI_API_KEY` to raise the E-utilities request allowance. Set `OPENROUTER_API_KEY` only when running Jev. Secrets belong in the environment or an ignored local configuration, never in tracked YAML.

## Ordered phases

1. **Taxonomy** — run `notebooks/01_taxonomy.ipynb`. It migrates the local V2 bacterial taxonomy and aliases, imports pinned Species Fungorum and ICTV MSL41 files, and joins VMR MSL41 isolate metadata to ICTV taxon IDs.
2. **PubMed** — run `notebooks/02_pubmed.ipynb`. NCBI searches are keyed by query hash and written incrementally to local JSONL caches. Existing query and PMID cache entries are reused. By default, all search results are retained for counts, while up to 200 PMIDs per species are selected for abstract fetch and scoring: the 30 lowest numeric PMIDs first, then a seeded uniform sample from the remaining results (`pubmed.sampling` in `configs/default.yaml`; set the cap to `null` for a full fetch). The random seed and per-row sample selection are recorded. `batch_size: 200` separately controls the number of PMID IDs per EFetch request. The default worker pool keeps 48 requests in flight, while a shared rate limiter paces request starts: 0.34 seconds without `NCBI_API_KEY`, and 0.125 seconds with a key (about 8 requests/second), leaving headroom under NCBI's documented limits of 3 requests/second without a key and 10 with one. `NCBI_EMAIL` identifies the tool contact; the API key, not email alone, enables the higher default rate. Search results above NCBI's 10,000 PMID ESearch response limit are expanded with cached overflow requests. The notebook exposes worker count, optional request-delay override, query-call limit, and fetch-batch limit. See the [NCBI usage guidelines](https://eutilities.github.io/site/API_Key/usageandkey/) and [ESearch limits](https://www.ncbi.nlm.nih.gov/books/NBK25499/).
   Before query construction, VMR isolate designations and configured ambiguous virus abbreviations are excluded. For fungi, synonym names that map to multiple taxon IDs and accepted names flagged as ambiguous are excluded from species-level searches. These source names remain in taxonomy provenance; unique fungal synonyms can still be searched. Preflight reports taxa without any unambiguous searchable name.
   For a first live check, run `notebooks/02a_pubmed_pilot.ipynb` (or `pathogen-atlas-v3-phase2-pilot`). It selects ten curated/reproducible taxa per group, caps each search at 500 returned PMIDs and fetching at eight batches, and writes only under `data/02_pubmed/pilot_v1/`, leaving the full Phase 2 release paths untouched. If NCBI rate-limits the API key or IP, raise `request_delay_seconds`; the pilot completed at 0.5 seconds between requests after an initial rate-limit response at 0.11 seconds.
3. **Jev decisions** — run `notebooks/03_decisions.ipynb`. The disease attribution gate is one typed question per selected species–PMID pair. Passing abstracts are sent once more with all applicable features for that organism group. Successful decisions resume only when organism, PMID, input hash, run, model, question version, and feature-profile version match. Failed and malformed outputs remain retryable. Concurrency, max calls, ETA, and remaining OpenRouter balance are visible.
4. **Release tables** — run `notebooks/04_release.ipynb`. Tables are deterministic transformations of saved decisions. The 0.75 pathogen gate is configurable; abstracts at or above it are accepted pathogen evidence. Feature thresholds are configurable and only affect feature summaries. The Pathogenic Interest Index retains its three input counts and normalizes within organism group.
5. **Embedding** — run `notebooks/05_embedding.ipynb`. Requires the optional visualization dependencies. It fits a shared-core map and separate organism-group maps. Non-applicable features are excluded from a profile; missing observations in selected features are mean-imputed within that map, never zero-imputed.

The Phase 3 notebook defaults to preflight only. Set `RUN_PHASE = True` and a bounded `MAX_CALLS` for a pilot. Set the limit to `None` only after the question set, taxonomy profile, estimated call volume, and available OpenRouter balance have been reviewed. Stop safely at any checkpoint and rerun the notebook to resume.

## Data and Git

All local data and generated outputs are ignored by V3 Git policy. The complete V2 `data/` tree is copied under `data/migration/version2_snapshot/`; the migration utility never writes to V2. Only code, configs, notebooks, docs, schemas, and migration manifests are synced. The V2 source commit and expected snapshot counts are recorded in `configs/migrations/v2_snapshot.yaml`.

## First-run limitations

The local V2 migration populates the bacterial cohort and historical records. Pinned Species Fungorum and ICTV MSL41 exports are also stored locally under `data/source_taxonomy/` and loaded by Phase 1. V2 model decisions are stored in the legacy table; V3 prompt/profile versions require fresh V3 decisions before current evidence can be released.
