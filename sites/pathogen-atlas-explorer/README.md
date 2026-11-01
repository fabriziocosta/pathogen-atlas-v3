# Pathogen Atlas Evidence Explorer

The Explorer is published at <https://fabriziocosta.github.io/pathogen-atlas-v3/>.

It is a static site. The HTML application is tracked in this repository; its generated data and evidence shards are distributed as the `explorer-site-data` GitHub Release asset and are assembled into a GitHub Pages artifact by `.github/workflows/deploy-explorer.yml`. At deploy time, the workflow exports the current three-level feature catalogue from `configs/features.yaml` and `configs/feature_profiles.yaml`, so the interface stays aligned with the versioned YAML without checking generated data into Git.

To update the public data snapshot:

1. Build the site data from the current V3 releases with `python sites/pathogen-atlas-explorer/scripts/build_data.py`.
2. From `sites/pathogen-atlas-explorer/dist`, create `explorer-site-data.tar.gz` containing `explorer-data.json.gz` and the `evidence/` directory (when present).
3. Upload the archive to the existing `explorer-site-data` release with `gh release upload explorer-site-data explorer-site-data.tar.gz --repo fabriziocosta/pathogen-atlas-v3 --clobber`.
4. Run the `Deploy Pathogen Atlas Explorer` workflow from GitHub Actions, or push a change to the Explorer HTML or this README.

The public snapshot is reviewed project output. Do not include API credentials, local configuration, raw provider responses, or unpublished sensitive material in it.
