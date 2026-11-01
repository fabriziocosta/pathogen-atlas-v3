<div align="center">
  <img src="assets/pathogen-atlas-logo.png" alt="Pathogen Atlas" width="150">
</div>

# Pathogen Atlas V3

The Pathogen Atlas brings together published reports of microorganisms associated with human disease and the evidence that supports those associations. Version 3 is being built to represent bacteria, fungi, and viruses in one auditable research resource while respecting the differences in their taxonomy and biology.

The Atlas is an evidence map, not a complete census of pathogens or proof that an organism is harmless when evidence is sparse. Its organism profiles reflect the literature retrieved, the abstracts available, and the questions used to assess them. Feature scores describe what those abstracts report; they are not measures of disease prevalence or risk.

V3 keeps a shared catalogue of research features and selects relevant features for each organism group. Species identity and literature evidence retain source and version information, and model probabilities remain separate from released labels until evaluated against human-labelled examples.

## Explore

The static multi-domain Evidence Explorer source is in [`sites/pathogen-atlas-explorer`](sites/pathogen-atlas-explorer/). Build its data from the Phase 4 releases after scored evidence and embedding coordinates exist, then run `python launch_explorer.py` from this directory. The current bootstrapped release contains taxonomy and migration inputs, but no new V3 Jev decisions or fitted coordinates yet.

## Project structure

- `pathogen_atlas_v3/` contains reusable, inspectable phase code.
- `configs/` contains versioned taxonomy, query, prompt, and feature definitions.
- `notebooks/` exposes phase inputs and calls in execution order; reusable logic stays in Python modules.
- `data/` holds the local V2 migration snapshot, caches, and generated releases. All data files are excluded from Git.
- `docs/` records data contracts, scientific decisions, migration rules, and methods.

See [the project workflow](docs/workflow.md) to install and run phases, [the multi-domain data model](docs/data_model.md) for organism and feature contracts, and [the embedding method](docs/embeddings.md) for species coordinates.
