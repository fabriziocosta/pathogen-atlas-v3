"""Canonical V3 data contracts shared across organism groups."""

from __future__ import annotations

from typing import Final

DOMAINS: Final[tuple[str, ...]] = ("bacteria", "fungi", "viruses")

SPECIES_COLUMNS: Final[tuple[str, ...]] = (
    "organism_id", "organism_group", "accepted_name", "source_taxon_id",
    "taxonomy_source", "taxonomy_release", "rank", "parent_taxon_id",
    "taxonomy_path", "source_record_json", "mapping_status",
)
ALIAS_COLUMNS: Final[tuple[str, ...]] = (
    "organism_id", "organism_group", "name", "name_type", "source",
    "source_record_id", "source_record_json", "mapping_status",
)
VIRUS_VARIANT_COLUMNS: Final[tuple[str, ...]] = (
    "organism_id", "organism_group", "source_taxon_id", "accepted_name",
    "taxonomy_source", "taxonomy_release", "isolate_id", "exemplar_status",
    "virus_name", "virus_name_abbreviation", "isolate_designation",
    "genbank_accession", "genome_coverage", "genome", "host_source",
    "source_record_json", "mapping_status",
)
ABSTRACT_COLUMNS: Final[tuple[str, ...]] = (
    "pmid", "title", "abstract", "publication_date", "publication_year",
    "authors", "journal", "doi", "mesh_terms", "keywords", "source",
)
SPECIES_ABSTRACT_COLUMNS: Final[tuple[str, ...]] = (
    "organism_id", "organism_group", "pmid", "matched_names",
    "query_id", "query_hash", "selected_for_fetch",
)
DECISION_COLUMNS: Final[tuple[str, ...]] = (
    "organism_id", "organism_group", "pmid", "decision_kind", "feature_id",
    "input_hash", "provider", "model", "question_version", "profile_version",
    "run_id", "p_yes", "status", "attempts", "error", "completed_at",
    "raw_response_path",
)


def validate_domain(domain: str) -> str:
    """Normalize and validate an organism-group identifier."""
    value = str(domain).strip().lower()
    if value not in DOMAINS:
        raise ValueError(f"Unsupported organism_group {domain!r}; expected one of {DOMAINS}.")
    return value
