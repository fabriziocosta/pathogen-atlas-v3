"""Migrate V2 local releases into the V3 canonical local data layout.

The V2 snapshot is retained unchanged as provenance. Imported decision rows are
historical only: they retain their original question/model/run versions and are
not promoted to V3 decisions unless a later compatibility review approves them.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

from pathogen_atlas_v3.common.io import atomic_parquet


def _read(root: Path, relative: str) -> pd.DataFrame:
    path = root / relative
    return pd.read_parquet(path) if path.exists() else pd.DataFrame()


def _safe_id(group: str, name: Any) -> str:
    import hashlib
    canonical = " ".join(str(name).strip().split()).casefold()
    return f"{group}:{hashlib.sha256(canonical.encode()).hexdigest()[:20]}"


def migrate_v2_snapshot(source_root: str | Path, output_root: str | Path) -> dict[str, Any]:
    """Convert compatible V2 release tables while preserving source columns."""
    source = Path(source_root)
    output = Path(output_root)
    species = _read(source, "01_species/release/species_canonical.parquet")
    taxonomy = _read(source, "01_species/release/species_taxonomy.parquet")
    aliases = _read(source, "01_species/release/species_synonym_mapping.parquet")
    abstracts = _read(source, "02_abstracts/release/abstracts.parquet")
    matches = _read(source, "02_abstracts/release/abstract_species_matches.parquet")
    gate = _read(source, "03_decisions/parsed/pathogen_evidence_jev.parquet")
    features = _read(source, "03_decisions/parsed/domain_expert_answers.parquet")
    v2_index = _read(source, "04_tables/release/pathogenic_index.parquet")

    if species.empty or taxonomy.empty:
        raise FileNotFoundError("V2 species_canonical and species_taxonomy releases are required.")

    species = species.merge(taxonomy, on="canonical_species", how="left", suffixes=("", "_taxonomy"))
    species["organism_id"] = species["canonical_species"].map(lambda x: _safe_id("bacteria", x))
    species["organism_group"] = "bacteria"
    species["accepted_name"] = species["canonical_species"]
    species["source_taxon_id"] = species.get("lpsn_record_id")
    species["taxonomy_source"] = "lpsn"
    species["taxonomy_release"] = species.get("source_retrieved_at", pd.Series("v2-local-capture", index=species.index))
    species["rank"] = "species"
    species["parent_taxon_id"] = pd.NA
    ranks = [c for c in ("domain", "kingdom", "phylum", "class", "order", "family", "genus") if c in species]
    species["taxonomy_path"] = species[ranks].fillna("").agg(";".join, axis=1) if ranks else ""
    species["source_record_json"] = "{}"
    species["mapping_status"] = species.get("status", "accepted")
    species_out = species[["organism_id", "organism_group", "accepted_name", "source_taxon_id", "taxonomy_source", "taxonomy_release", "rank", "parent_taxon_id", "taxonomy_path", "source_record_json", "mapping_status"]].drop_duplicates("organism_id")

    if not aliases.empty:
        aliases = aliases.merge(species_out[["organism_id", "accepted_name"]], left_on="canonical_species", right_on="accepted_name", how="left")
        aliases["organism_group"] = "bacteria"
        aliases["name"] = aliases["input_name"]
        aliases["name_type"] = aliases.get("mapping_type", "synonym")
        aliases["source_record_id"] = aliases.get("input_lpsn_record_id")
        aliases["source_record_json"] = "{}"
        aliases["mapping_status"] = aliases.get("status", "mapped")
        aliases_out = aliases[["organism_id", "organism_group", "name", "name_type", "source", "source_record_id", "source_record_json", "mapping_status"]]
    else:
        aliases_out = pd.DataFrame(columns=["organism_id", "organism_group", "name", "name_type", "source", "source_record_id", "source_record_json", "mapping_status"])

    if not abstracts.empty:
        abstracts["pmid"] = abstracts["pmid"].astype(str)
        abstract_out = abstracts.drop_duplicates("pmid")
    else:
        abstract_out = abstracts
    if not matches.empty:
        lookup = species_out[["organism_id", "accepted_name"]]
        matches = matches.merge(lookup, left_on="canonical_species", right_on="accepted_name", how="left")
        matches["organism_group"] = "bacteria"
        matches["matched_names"] = matches.get("matched_names", pd.Series([[]] * len(matches), index=matches.index))
        matches["selected_for_fetch"] = matches.get("selected_for_fetch", True)
        match_out = matches[["organism_id", "organism_group", "accepted_name", "pmid", "matched_names", "query_id", "query_hash", "selected_for_fetch"]].dropna(subset=["organism_id"])
    else:
        match_out = matches

    # Keep all inherited decision fields, plus canonical organism identity.
    id_lookup = species_out.set_index("accepted_name")["organism_id"]
    legacy_decisions = []
    for kind, frame in (("pathogen_evidence", gate), ("feature", features)):
        if frame.empty:
            continue
        work = frame.copy()
        work["organism_id"] = work["canonical_species"].map(id_lookup)
        work["organism_group"] = "bacteria"
        work["decision_kind"] = kind
        work["feature_id"] = work["question_id"] if kind == "feature" else "human_pathogen_evidence"
        work["question_version"] = work.get("question_version", work.get("question_set_version", "v2-version-unknown"))
        work["profile_version"] = work.get("question_set_version", "v2-domain-expert-questions-73-v1")
        work["run_id"] = work.get("run_id", work.get("expert_run_id", work.get("gate_run_id", "v2-run-unknown")))
        work = work.rename(columns={"canonical_species": "legacy_canonical_species"})
        legacy_decisions.append(work)
    legacy_out = pd.concat(legacy_decisions, ignore_index=True, sort=False) if legacy_decisions else pd.DataFrame()

    if not v2_index.empty:
        v2_index["organism_id"] = v2_index["canonical_species"].map(id_lookup)
        v2_index["organism_group"] = "bacteria"
        v2_index["legacy_index_version"] = "pathogenic-interest-index-v2"
        v2_index = v2_index.rename(columns={"canonical_species": "legacy_canonical_species", "n_causation": "n_pathogen_pmids_v2", "n_articles": "n_pubmed_search_records_v2", "n_sampled": "n_abstracts_fetched_v2", "pathogenic_index": "pathogenic_interest_index_v2"})
        legacy_index = v2_index
    else:
        legacy_index = pd.DataFrame(columns=["organism_id", "organism_group", "legacy_canonical_species", "n_pathogen_pmids_v2", "n_pubmed_search_records_v2", "n_abstracts_fetched_v2", "estimated_causation", "log_estimated_causation", "pathogenic_interest_index_v2", "legacy_index_version"])

    outputs = {
        "data/01_taxonomy/release/species.parquet": species_out,
        "data/01_taxonomy/release/aliases.parquet": aliases_out,
        "data/02_pubmed/release/abstracts.parquet": abstract_out,
        "data/02_pubmed/release/species_abstracts.parquet": match_out,
        "data/03_decisions/legacy_v2/decisions.parquet": legacy_out,
        "data/03_decisions/legacy_v2/pathogenic_interest_index.parquet": legacy_index,
    }
    for relative, frame in outputs.items():
        if not frame.empty or relative.endswith("aliases.parquet") or relative.endswith("decisions.parquet"):
            atomic_parquet(frame, output / relative)

    report = {
        "migration_id": "v2-to-v3-bacteria-v1",
        "migrated_at": datetime.now(timezone.utc).isoformat(),
        "source_root": str(source.resolve()),
        "source_data_snapshot": "V2 local snapshot; no source files modified",
        "counts": {
            "species": len(species_out), "aliases": len(aliases_out),
            "abstracts": len(abstract_out), "species_abstract_matches": len(match_out),
            "legacy_decisions": len(legacy_out),
            "legacy_index_rows": len(legacy_index),
        },
        "outputs": {key: len(value) for key, value in outputs.items()},
        "legacy_decisions_promoted": False,
        "notes": "Feature/gate scores retain V2 provenance. Changed V3 question or profile versions require rescoring.",
    }
    report_path = output / "data/migration/v2_migration_manifest.json"
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", default="data/migration/version2_snapshot")
    parser.add_argument("--output", default=".")
    args = parser.parse_args()
    report = migrate_v2_snapshot(args.source, args.output)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
