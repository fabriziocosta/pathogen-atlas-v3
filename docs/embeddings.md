# Species embedding method

V3 stores a fitted two-dimensional coordinate per organism and embedding scope. The **shared-core** scope uses only catalogue features declared for bacteria, fungi, and viruses. Each **within-domain** scope uses features selected for that organism group. A feature that does not apply is excluded from that scope rather than encoded as zero.

The matrix contains per-species feature prevalence: the fraction of scored pathogen-evidence abstracts that answered yes. Unscored observations remain missing. Within a fitted scope, remaining missing values are imputed with the observed mean for that feature, then the columns are standardized before IMSTE is fitted. The chosen feature IDs, seed, and settings are written beside the coordinate table. Re-running the same build with identical inputs and settings is intended to produce the same release coordinates.

IMSTE is an optional pinned dependency installed with `pip install -e '.[viz]'`. The V3 Explorer consumes saved coordinates; it never fits the embedding in the browser. Shared-core and domain-specific coordinates have distinct scopes and should not be mixed in a single map view.
