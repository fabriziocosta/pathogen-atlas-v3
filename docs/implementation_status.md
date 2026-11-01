# V3 implementation status

## Implemented locally

- Independent V3 Git repository and phase-oriented Python package.
- Complete local copy of V2 `data/` plus the existing species Wikipedia/GBIF association tables.
- Pinned Species Fungorum (101,843 accepted species), ICTV MSL41 (17,554 virus species), and matching ICTV VMR MSL41 (19,285 isolate records) source files. Forty ambiguous fungal accepted names and 21,291 unresolved alias mappings remain traceable for review.
- Migration to V3 organism, alias, abstract-match, and legacy-decision schemas, with the old Pathogenic Interest Index retained separately.
- 73 inherited bacteria features and initial fungal/viral candidate profiles in a versioned catalogue.
- Source-file adapters for LPSN, the dated Species Fungorum snapshot, and pinned ICTV MSL41; full taxonomy source rows are retained in provenance JSON.
- Rate-limited resumable PubMed acquisition, typed Jev gate/feature calls, local JSONL raw responses, parquet decision ledger, Phase 4 tables, IMSTE coordinates, and a static Explorer build. The current viral release contains 187,325 completed pathogen-gate assessments; the Explorer exposes their PMID, title, Jev confidence, and accepted pathogen evidence at the configured 0.75 threshold.
- Fake-client tests for all three groups, decision resume/retry and taxonomy/index contracts.

## Still required before a production multi-domain release

- Review the 40 fungal accepted-name collisions and the source-level fungal/virus mapping decisions before publishing a release.
- Review each provisional shared feature and fungal/viral feature with domain experts; curate labelled examples for the gate and feature questions.
- Run and evaluate bounded Jev pilots for each domain. The configured 0.75 pathogen gate currently defines accepted evidence; feature interpretations still need expert review and evaluation.
- Initialize matching bacterial V2 gate and feature scores into the V3 ledger with provenance; re-score only changed or missing inputs.
- Fit and inspect V3 shared and domain-specific IMSTE spaces after scored features exist.
- Extend Wikipedia and GBIF link association coverage to the new fungal/viral cohorts.

V3 viral PubMed acquisition covers all 17,554 taxonomy species and has 144,811 fetched abstracts. The viral pathogen-gate run is complete for 187,325 sampled species–PMID records. Feature scoring has 46,299 completed and 303,902 failed decision rows, so virus feature profiles are incomplete; retrying those failures and evaluating the Jev labels remain production-release requirements. Raw provider responses and source caches stay local and are not included in the public Explorer package.
