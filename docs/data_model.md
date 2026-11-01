# Multi-domain data model

Each record carries an `organism_group` (`bacteria`, `fungi`, or `viruses`) independently of its source-specific taxonomy ranks. A normalized species row preserves the accepted name, source taxon identifier, taxonomy source, release identifier, rank, parent identifier, and a readable taxonomy path. Alias records retain their source and mapping type. Identical names across domains remain separate identities.

The fungal snapshot supplies Species Fungorum accepted names and taxonomy together with Index Fungorum nomenclatural IDs, source links, and synonym mappings. The virus cohort also has a `virus_variants.parquet` table from ICTV's matching Virus Metadata Resource, keyed to MSL species IDs and retaining exemplar/additional isolate names, designations, accessions, genome details, host source, and the complete source row. Ambiguous accepted fungal names and unresolved alias or taxon mappings remain traceable for review rather than being silently discarded.

## Feature catalogue

`configs/features.yaml` is the versioned, editable hierarchy of features: **dimension → category → feature**. The dimension identifies what kind of information a feature measures; the category groups related concepts within that dimension; and the feature contains the actual question and scoring instructions. `configs/feature_profiles.yaml` defines the stable scoring order and applicability for each organism group. The loader flattens the hierarchy for the decision runner while preserving dimension and category metadata for the Explorer. Every feature keeps a stable ID, yes/no definition, answer type, question version, explicit applicability, evidence instructions, and review state. Reclassifying a feature does not change its scoring cache identity; changing its interpretation requires a new question version.

The top-level dimensions and their scope are documented in [Feature taxonomy decision](feature_taxonomy_decision.md). In brief, they distinguish pathogen biology, transmission, reservoirs and exposure, clinical presentation, host factors, healthcare context, epidemiology and public health, and abstract reporting. Feature answers remain statements about what the abstract reports; they are not automatically claims of biological truth or causal validation.

The 73 inherited bacteria questions preserve their legacy identifiers and source mapping. Fungal and viral candidates add organism-specific concepts. Cross-domain questions and domain-specific profiles are marked provisional until reviewed by domain experts and tested on labelled examples; provisional results are for bounded pilots. `not_scored` and `not_applicable` are not “no”.

## Decisions and evidence

The disease gate and feature assessments share one normalized ledger. Its cache identity includes organism group and stable ID, PMID, abstract/title hash, model, question version, feature-profile version, and run ID. One Jev request contains all currently applicable feature questions for one passing abstract. Raw responses are appended before normalization; malformed or failed rows remain in the ledger with errors and can be retried.

The disease gate defaults to `p_yes >= 0.75` for feature selection. This is a pilot setting, not an approved release threshold. Phase 4 writes candidate results but does not mark them accepted while `release_evaluation_approved` is false. Do not interpret candidate probabilities or provisional features as final labels.

## Pathogenic Interest Index

For each organism, let:

- `n_pathogen_pmids` be the count of distinct accepted species–PMID records;
- `n_abstracts_fetched` be the count of distinct organism-associated abstracts actually present;
- `n_pubmed_search_records` be the count of distinct species-associated PMIDs returned by the literature searches.

The estimated support is `n_pathogen_pmids / n_abstracts_fetched * n_pubmed_search_records`. The index is `log(estimated support)` divided by the maximum log estimate within the same organism group. It is zero when a component makes an estimate unavailable or there are no accepted records. All components and the index formula version are kept in the table. Group-wise normalization means index values are not directly comparable across bacteria, fungi, and viruses.
