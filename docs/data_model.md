# Multi-domain data model

Each record carries an `organism_group` (`bacteria`, `fungi`, or `viruses`) independently of its source-specific taxonomy ranks. A normalized species row preserves the accepted name, source taxon identifier, taxonomy source, release identifier, rank, parent identifier, and a readable taxonomy path. Alias records retain their source and mapping type. Identical names across domains remain separate identities.

The fungal snapshot supplies Species Fungorum accepted names and taxonomy together with Index Fungorum nomenclatural IDs, source links, and synonym mappings. The virus cohort also has a `virus_variants.parquet` table from ICTV's matching Virus Metadata Resource, keyed to MSL species IDs and retaining exemplar/additional isolate names, designations, accessions, genome details, host source, and the complete source row. Ambiguous accepted fungal names and unresolved alias or taxon mappings remain traceable for review rather than being silently discarded.

## Feature catalogue

`configs/features.yaml` is the versioned, editable hierarchy of features: **dimension → category → feature**. The dimension identifies what kind of information a feature measures; the category groups related concepts within that dimension; and the feature contains the actual question and scoring instructions. `configs/feature_profiles.yaml` defines the stable scoring order and applicability for each organism group. The loader flattens the hierarchy for the decision runner while preserving dimension and category metadata for the Explorer. Every feature keeps a stable ID, yes/no definition, answer type, question version, explicit applicability, evidence instructions, and review state. Reclassifying a feature does not change its scoring cache identity; changing its interpretation requires a new question version.

The top-level dimensions and their scope are documented in [Feature taxonomy decision](feature_taxonomy_decision.md). In brief, they distinguish pathogen biology, transmission, reservoirs and exposure, clinical presentation, host factors, healthcare context, epidemiology and public health, and abstract reporting. Feature answers remain statements about what the abstract reports; they are not automatically claims of biological truth or causal validation.

The 73 inherited bacteria questions preserve their legacy identifiers and source mapping. Fungal and viral candidates add organism-specific concepts. Cross-domain questions and domain-specific profiles are marked provisional until reviewed by domain experts and tested on labelled examples; provisional results are for bounded pilots. `not_scored` and `not_applicable` are not “no”.

## Decisions and evidence

The disease gate and feature assessments share one normalized ledger. Its cache identity includes organism group and stable ID, PMID, abstract/title hash, model, question version, feature-profile version, and run ID. One Jev request contains all currently applicable feature questions for one passing abstract. Raw responses are appended before normalization; malformed or failed rows remain in the ledger with errors and can be retried.

The disease gate defaults to `p_yes >= 0.75`. A sampled abstract that meets this threshold is accepted pathogen evidence and is included in the Phase 4 release. Feature scores remain separate supporting characterizations; they do not override the disease-attribution gate.

## Pathogenic Interest Index

For each organism, let:

- `n_pathogen_pmids` be the count of distinct accepted species–PMID records;
- `n_abstracts_fetched` be the count of distinct organism-associated abstracts actually present;
- `n_pubmed_search_records` be the count of distinct species-associated PMIDs returned by the literature searches.

V3 uses the V2 sampling design by default: up to 200 PMIDs per species, keeping the 30 lowest numeric PMIDs first (the V2 proxy for oldest records), then drawing the remaining slots uniformly with seed 42. These settings are configurable under `pubmed.sampling`; set `max_pmids_per_species: null` to fetch all search results. `batch_size: 200` is a separate setting for PMIDs per EFetch request. The match table retains every search result and marks sample membership with `selected_for_sampling`, its selection reason with `sample_selection`, and successful abstract availability with `selected_for_fetch`.

The estimated support is `n_pathogen_pmids / n_abstracts_fetched * n_pubmed_search_records`. `n_pubmed_search_records` counts the full search result set, while `n_abstracts_fetched` counts successfully fetched abstracts from the per-species sample. The estimate scales the sampled evidence rate to the full search count. The index is `log(estimated support)` divided by the maximum log estimate within the same organism group. It is zero when a component makes an estimate unavailable or there are no accepted records. All component counts, sample settings, and the index version are recorded. Group-wise normalization means index values are not directly comparable across bacteria, fungi, and viruses.
