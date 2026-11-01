"""Taxonomy file adapters for LPSN, fungal checklists, and ICTV MSL exports.

Adapters preserve the source taxon identifiers and source release. They accept
CSV, TSV, or Excel exports so a release is pinned locally before processing.
"""

from __future__ import annotations

import hashlib
import json
import re
import zipfile
from pathlib import Path

import pandas as pd

from pathogen_atlas_v3.common.schemas import DOMAINS, validate_domain


def _text(value: object) -> str:
    return "" if pd.isna(value) else str(value).strip()


def _id(group: str, source: str, taxon_id: object, name: str) -> str:
    identity = _text(taxon_id) or _text(name).casefold()
    digest = hashlib.sha256(f"{group}|{source}|{identity}".encode()).hexdigest()[:20]
    return f"{group}:{digest}"


def read_source_file(path: str | Path, *, sheet_name: str | int = 0) -> pd.DataFrame:
    """Read a pinned taxonomy release file based on its extension."""
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix in {".tsv", ".tab"}:
        return pd.read_csv(path, sep="\t", dtype=object)
    if suffix in {".xlsx", ".xls"}:
        return pd.read_excel(path, sheet_name=sheet_name, dtype=object)
    if suffix in {".csv", ".txt"}:
        return pd.read_csv(path, dtype=object)
    if suffix == ".parquet":
        return pd.read_parquet(path)
    if suffix == ".zip":
        with zipfile.ZipFile(path) as archive:
            csv_members = [name for name in archive.namelist() if name.lower().endswith(".csv")]
            if len(csv_members) != 1:
                raise ValueError(f"Expected one CSV in taxonomy archive {path}; found {len(csv_members)}.")
            with archive.open(csv_members[0]) as stream:
                return pd.read_csv(stream, dtype=object)
    raise ValueError(f"Unsupported taxonomy file type: {suffix}")


def normalize_taxonomy_records(
    records: pd.DataFrame,
    *,
    organism_group: str,
    source: str,
    release: str,
    name_column: str,
    taxon_id_column: str | None = None,
    alias_column: str | None = None,
    rank_column: str | None = None,
    parent_id_column: str | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Normalize source rows to V3 species and alias records.

    Input must contain accepted species-level names. Taxonomic classification
    can be retained in ``taxonomy_path`` and remains source-specific.
    """
    group = validate_domain(organism_group)
    if name_column not in records:
        raise ValueError(f"Missing accepted-name column {name_column!r}.")
    if taxon_id_column and taxon_id_column not in records:
        raise ValueError(f"Missing taxon ID column {taxon_id_column!r}.")
    normalized_columns = {str(column).strip().casefold(): column for column in records.columns}
    rank_order = ("realm", "subrealm", "kingdom", "subkingdom", "phylum", "subphylum", "class", "subclass", "order", "suborder", "family", "subfamily", "genus", "subgenus")
    rank_columns = [(rank, normalized_columns[rank]) for rank in rank_order if rank in normalized_columns]
    species_rows: list[dict[str, object]] = []
    alias_rows: list[dict[str, object]] = []
    for row in records.to_dict("records"):
        accepted = _text(row.get(name_column))
        if not accepted:
            continue
        taxon_id = row.get(taxon_id_column) if taxon_id_column else None
        organism_id = _id(group, source, taxon_id, accepted)
        hierarchy = ";".join(f"{rank}={_text(row.get(column))}" for rank, column in rank_columns if _text(row.get(column)))
        species_rows.append({
            "organism_id": organism_id, "organism_group": group,
            "accepted_name": accepted, "source_taxon_id": _text(taxon_id) or None,
            "taxonomy_source": source, "taxonomy_release": release,
            "rank": _text(row.get(rank_column)) if rank_column else "species",
            "parent_taxon_id": _text(row.get(parent_id_column)) or None if parent_id_column else None,
            "taxonomy_path": hierarchy,
            "source_record_json": json.dumps(row, default=str, ensure_ascii=False, sort_keys=True),
            "mapping_status": "accepted",
        })
        alias_rows.append({
            "organism_id": organism_id, "organism_group": group,
            "name": accepted, "name_type": "accepted", "source": source,
            "source_record_id": _text(taxon_id) or None,
            "source_record_json": json.dumps(row, default=str, ensure_ascii=False, sort_keys=True),
            "mapping_status": "accepted",
        })
        if alias_column and alias_column in records:
            raw_aliases = row.get(alias_column)
            if isinstance(raw_aliases, (list, tuple, set)):
                values = raw_aliases
            else:
                values = str(raw_aliases).replace("|", ";").split(";") if _text(raw_aliases) else []
            for alias in values:
                alias = _text(alias)
                if alias and alias.casefold() != accepted.casefold():
                    alias_rows.append({
                        "organism_id": organism_id, "organism_group": group,
                        "name": alias, "name_type": "synonym", "source": source,
                        "source_record_id": _text(taxon_id) or None,
                        "source_record_json": json.dumps(row, default=str, ensure_ascii=False, sort_keys=True),
                        "mapping_status": "source_alias",
                    })
    species = pd.DataFrame(species_rows)
    aliases = pd.DataFrame(alias_rows)
    if species.empty:
        raise ValueError("No accepted species records were found.")
    if species["source_taxon_id"].dropna().duplicated().any():
        raise ValueError(f"Source taxon IDs are not unique for {group} records.")
    collisions = species["accepted_name"].str.casefold().duplicated(keep=False)
    if collisions.any():
        collision_ids = set(species.loc[collisions, "organism_id"])
        species.loc[species["organism_id"].isin(collision_ids), "mapping_status"] = "ambiguous_accepted_name"
        aliases.loc[aliases["organism_id"].isin(collision_ids), "mapping_status"] = "ambiguous_accepted_name"
    return species, aliases


def load_lpsn_export(path: str | Path, *, release: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Normalize an LPSN export; V2 parquet tables are handled by migration."""
    frame = read_source_file(path)
    name = next((c for c in ("canonical_species", "canonical_name", "name") if c in frame), None)
    identifier = next((c for c in ("lpsn_record_id", "record_id", "id") if c in frame), None)
    if name is None:
        raise ValueError("LPSN input must include canonical_species, canonical_name, or name.")
    return normalize_taxonomy_records(frame, organism_group="bacteria", source="lpsn", release=release, name_column=name, taxon_id_column=identifier, alias_column="synonyms")


def load_species_fungorum_export(path: str | Path, *, release: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    frame = read_source_file(path)
    columns = {str(column).strip().casefold(): column for column in frame.columns}
    name = next((columns[key] for key in ("accepted_name", "current_name", "scientificname", "name") if key in columns), None)
    identifier = next((columns[key] for key in ("species_fungorum_id", "index fungorum id", "if_id", "record_id", "id") if key in columns), None)
    rank_col = columns.get("rank")
    status_col = columns.get("taxonomic status")
    if name is None or identifier is None:
        raise ValueError("Species Fungorum input requires an accepted-name field and source record ID.")
    if rank_col and status_col:
        is_species = frame[rank_col].fillna("").astype(str).str.casefold().eq("species")
        is_accepted = frame[status_col].fillna("").astype(str).str.casefold().eq("accepted")
        species_rows = frame.loc[is_species & is_accepted].copy()
        if species_rows.empty:
            raise ValueError("Species Fungorum export contains no accepted species rows.")
        species, accepted_aliases = normalize_taxonomy_records(species_rows, organism_group="fungi", source="species_fungorum", release=release, name_column=name, taxon_id_column=identifier)
        accepted_id_col = columns.get("accepted name index fungorum id", identifier)
        accepted_ids = {str(row.source_taxon_id): row.organism_id for row in species.itertuples(index=False)}
        canonical_by_id = {row.organism_id: row.accepted_name for row in species.itertuples(index=False)}
        ambiguous_ids = set(species.loc[species["mapping_status"].eq("ambiguous_accepted_name"), "organism_id"])
        alias_rows = accepted_aliases.to_dict("records")
        for row in frame.loc[is_species].to_dict("records"):
            organism_id = accepted_ids.get(_text(row.get(accepted_id_col)))
            alias_name = _text(row.get(name))
            if organism_id and alias_name:
                canonical = canonical_by_id[organism_id]
                if alias_name.casefold() != canonical.casefold():
                    alias_rows.append({"organism_id": organism_id, "organism_group": "fungi", "name": alias_name, "name_type": "synonym", "source": "species_fungorum", "source_record_id": _text(row.get(identifier)) or None, "source_record_json": json.dumps(row, default=str, ensure_ascii=False, sort_keys=True), "mapping_status": "ambiguous_accepted_name" if organism_id in ambiguous_ids else "source_alias"})
            elif alias_name:
                alias_rows.append({"organism_id": None, "organism_group": "fungi", "name": alias_name, "name_type": "unresolved_source_name", "source": "species_fungorum", "source_record_id": _text(row.get(identifier)) or None, "source_record_json": json.dumps(row, default=str, ensure_ascii=False, sort_keys=True), "mapping_status": "unresolved_accepted_name_target"})
        aliases = pd.DataFrame(alias_rows).drop_duplicates(["organism_id", "name", "source"])
        return species, aliases
    return normalize_taxonomy_records(frame, organism_group="fungi", source="species_fungorum", release=release, name_column=name, taxon_id_column=identifier, alias_column="synonyms")


def load_ictv_msl_export(path: str | Path, *, release: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    frame = read_source_file(path, sheet_name="MSL") if Path(path).suffix.lower() in {".xlsx", ".xls"} else read_source_file(path)
    columns = {str(column).strip().casefold(): column for column in frame.columns}
    name = next((columns[key] for key in ("species_name", "species", "virus species", "name") if key in columns), None)
    identifier = next((columns[key] for key in ("ictv_taxon_id", "species #", "taxon id", "taxon_id", "ictv_id") if key in columns), None)
    if name is None:
        raise ValueError("ICTV MSL input is missing a species-name field.")
    frame = frame.loc[frame[name].fillna("").astype(str).str.strip().ne("")].copy()
    aliases_column = next((columns[key] for key in ("synonyms", "virus names", "virus name(s)") if key in columns), None)
    species, aliases = normalize_taxonomy_records(frame, organism_group="viruses", source="ictv_msl", release=release, name_column=name, taxon_id_column=identifier, alias_column=aliases_column)
    aliases["name_type"] = aliases["name_type"].replace({"synonym": "virus_or_strain_alias"})
    return species, aliases


def load_ictv_vmr_export(path: str | Path, species: pd.DataFrame, *, release: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Normalize ICTV VMR exemplar/additional isolate metadata and name aliases."""
    frame = read_source_file(path, sheet_name="VMR MSL41") if Path(path).suffix.lower() in {".xlsx", ".xls"} else read_source_file(path)
    columns = {str(column).strip().casefold(): column for column in frame.columns}
    taxon_col = next((columns[key] for key in ("ictv_id", "ictv taxon id", "taxon id") if key in columns), None)
    if taxon_col is None:
        raise ValueError("ICTV VMR input must include the ICTV_ID taxon field.")
    ids = {str(row.source_taxon_id): (row.organism_id, row.accepted_name) for row in species.itertuples(index=False)}
    field_map = {
        "isolate_id": ("isolate id",), "exemplar_status": ("exemplar or additional isolate",),
        "virus_name": ("virus name(s)", "virus name"),
        "virus_name_abbreviation": ("virus name abbreviation(s)", "virus name abbreviation"),
        "isolate_designation": ("virus isolate designation", "isolate designation"),
        "genbank_accession": ("virus genbank accession", "genbank accession"),
        "genome_coverage": ("genome coverage",), "genome": ("genome",), "host_source": ("host source",),
    }
    resolved = {key: next((columns[name] for name in names if name in columns), None) for key, names in field_map.items()}
    variant_rows = []
    alias_rows = []
    text_alias_fields = ("virus_name", "virus_name_abbreviation", "isolate_designation")
    for row in frame.to_dict("records"):
        source_taxon_id = _text(row.get(taxon_col))
        organism_id, accepted_name = ids.get(source_taxon_id, (None, None))
        variant = {
            "organism_id": organism_id, "organism_group": "viruses",
            "source_taxon_id": source_taxon_id or None, "accepted_name": accepted_name,
            "taxonomy_source": "ictv_vmr", "taxonomy_release": release,
            "source_record_json": json.dumps(row, default=str, ensure_ascii=False, sort_keys=True),
            "mapping_status": "mapped" if organism_id else "unresolved_ictv_taxon_id",
        }
        for key, column in resolved.items():
            variant[key] = _text(row.get(column)) or None if column else None
        variant_rows.append(variant)
        for key in text_alias_fields:
            alias_value = _text(variant.get(key))
            if not alias_value:
                continue
            delimiter = r"\s*[;|,]\s*" if key == "virus_name_abbreviation" else r"\s*[;|]\s*"
            for alias_name in filter(None, (part.strip() for part in re.split(delimiter, alias_value))):
                alias_rows.append({
                    "organism_id": organism_id, "organism_group": "viruses", "name": alias_name,
                    "name_type": "virus_or_strain_alias", "source": "ictv_vmr",
                    "source_record_id": _text(variant.get("isolate_id")) or None,
                    "source_record_json": variant["source_record_json"],
                    "mapping_status": "source_alias" if organism_id else "unresolved_ictv_taxon_id",
                })
    variants = pd.DataFrame(variant_rows)
    aliases = pd.DataFrame(alias_rows)
    if not aliases.empty:
        aliases = aliases.drop_duplicates(["organism_id", "name", "source", "source_record_id"])
    return variants, aliases


ADAPTERS = {"bacteria": load_lpsn_export, "fungi": load_species_fungorum_export, "viruses": load_ictv_msl_export}
