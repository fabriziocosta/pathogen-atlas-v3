# V3 implementation status

## Implemented locally

- Independent V3 Git repository and phase-oriented Python package.
- Complete local copy of V2 `data/` plus the existing species Wikipedia/GBIF association tables.
- Pinned Species Fungorum (101,843 accepted species), ICTV MSL41 (17,554 virus species), and matching ICTV VMR MSL41 (19,285 isolate records) source files. Forty ambiguous fungal accepted names and 21,291 unresolved alias mappings remain traceable for review.
- Migration to V3 organism, alias, abstract-match, and legacy-decision schemas, with the old Pathogenic Interest Index retained separately.
- 73 inherited bacteria features and initial fungal/viral candidate profiles in a versioned catalogue.
- Source-file adapters for LPSN, the dated Species Fungorum snapshot, and pinned ICTV MSL41; full taxonomy source rows are retained in provenance JSON.
- Rate-limited resumable PubMed acquisition, typed Jev gate/feature calls, local JSONL raw responses, parquet decision ledger, Phase 4 tables, IMSTE coordinates, and a static Explorer build.
- Fake-client tests for all three groups, decision resume/retry and taxonomy/index contracts.

## Still required before a production multi-domain release

- Review the 40 fungal accepted-name collisions and the source-level fungal/virus mapping decisions before publishing a release.
- Review each provisional shared feature and fungal/viral feature with domain experts; curate labelled examples for the gate and feature questions.
- Run and evaluate bounded Jev pilots for each domain. Keep `release_evaluation_approved: false` until thresholds are supported by that evaluation.
- Re-score bacterial abstracts under the V3 prompt/profile versions; V2 decisions remain historical.
- Fit and inspect V3 shared and domain-specific IMSTE spaces after scored features exist.
- Extend Wikipedia and GBIF link association coverage to the new fungal/viral cohorts.

A bounded live PubMed pilot has run for ten taxa in each organism group. Its query and abstract caches are local under `data/02_pubmed/`; see `pilot_review.md` in the `pilot_v2_filtered_aliases` directory. No Jev calls have been made, and the capped search results are not full acquisition releases.
