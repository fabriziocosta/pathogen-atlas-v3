"""Build compact explorer metadata and lazy evidence shards from V3 releases."""

from __future__ import annotations

import gzip
import hashlib
import json
import shutil
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from pathogen_atlas_v3.common.io import load_yaml  # noqa: E402
from pathogen_atlas_v3.decisions.features import load_feature_catalog  # noqa: E402


def _safe_json(value):
    if isinstance(value, (list, tuple, set)): return list(value)
    if value is None: return None
    try:
        if pd.isna(value): return None
    except (TypeError, ValueError):
        pass
    if hasattr(value, "item"):
        try: return value.item()
        except ValueError: pass
    return value


def build_site_data(config_path: str = "configs/default.yaml") -> dict:
    config = load_yaml(ROOT / config_path)
    paths = config["paths"]
    release = ROOT / paths["release_root"]
    taxonomy = ROOT / paths["taxonomy_release"]
    species_path = release / "species_summary.parquet"
    feature_path = release / "feature_summary.parquet"
    coordinate_path = ROOT / paths["embedding_root"] / "coordinates.parquet"
    if not species_path.exists():
        raise FileNotFoundError("Run Phase 4 before building Explorer data.")
    species = pd.read_parquet(species_path)
    features = pd.read_parquet(feature_path) if feature_path.exists() else pd.DataFrame()
    coords = pd.read_parquet(coordinate_path) if coordinate_path.exists() else pd.DataFrame()
    if not coords.empty:
        # The shared coordinates provide one cross-domain placement. A domain-specific
        # map can be chosen by changing embedding_scope in the client UI later.
        shared = coords.loc[coords["embedding_scope"].eq("shared_core")]
        if not shared.empty:
            coords = shared.drop_duplicates("organism_id")[["organism_id", "x", "y"]]
        else:
            coords = coords.drop_duplicates("organism_id")[["organism_id", "x", "y"]]
        species = species.merge(coords, on="organism_id", how="left")
    else:
        species["x"], species["y"] = None, None
    aliases_path = taxonomy / "aliases.parquet"
    aliases = pd.read_parquet(aliases_path) if aliases_path.exists() else pd.DataFrame()
    aliases_by_id = aliases.groupby("organism_id")["name"].apply(lambda x: list(dict.fromkeys(map(str, x.dropna())))).to_dict() if not aliases.empty else {}
    legacy_link_root = ROOT / "data/migration/v2_site_snapshot"
    wikipedia_path, gbif_path = legacy_link_root / "species_wikipedia_links.csv", legacy_link_root / "species_gbif_links.csv"
    wikipedia = pd.read_csv(wikipedia_path) if wikipedia_path.exists() else pd.DataFrame()
    gbif = pd.read_csv(gbif_path) if gbif_path.exists() else pd.DataFrame()
    wikipedia_by_name = {str(r.canonical_species).casefold(): r.wikipedia_url for r in wikipedia.itertuples(index=False) if pd.notna(getattr(r, "wikipedia_url", None)) and str(r.wikipedia_url).strip()} if not wikipedia.empty else {}
    gbif_by_name = {str(r.canonical_species).casefold(): r.gbif_url for r in gbif.itertuples(index=False) if pd.notna(getattr(r, "gbif_url", None)) and str(r.gbif_url).strip()} if not gbif.empty else {}
    feature_map: dict[str, dict] = {}
    feature_definitions = {}
    if not features.empty:
        for row in features.itertuples(index=False):
            feature_map.setdefault(str(row.organism_id), {})[str(row.feature_id)] = None if pd.isna(row.share_yes) else float(row.share_yes)
            feature_definitions[str(row.feature_id)] = {"label": str(row.feature_id).replace("_", " ").upper(), "definition": "" if pd.isna(row.definition) else str(row.definition), "category": "" if pd.isna(row.category) else str(row.category)}
    catalog, profiles = load_feature_catalog(ROOT / "configs/features.yaml", ROOT / "configs/feature_profiles.yaml", allow_provisional=True)
    for fid, row in catalog["features"].items():
        feature_definitions[fid] = {
            "label": row.get("label", fid), "definition": row.get("definition", "") or "",
            "dimension_id": row["dimension_id"], "dimension": row["dimension_label"],
            "category_id": row["category_id"], "category": row["category_label"],
        }
    evidence_path = release / "accepted_pathogen_evidence.parquet"
    evidence = pd.read_parquet(evidence_path) if evidence_path.exists() else pd.DataFrame()
    feature_decisions_path = release / "feature_decisions.parquet"
    feature_decisions = pd.read_parquet(feature_decisions_path) if feature_decisions_path.exists() else pd.DataFrame()
    abstracts_path = ROOT / paths["pubmed_release"] / "abstracts.parquet"
    abstracts = pd.read_parquet(abstracts_path) if abstracts_path.exists() else pd.DataFrame()
    if not evidence.empty and not abstracts.empty:
        evidence = evidence.merge(abstracts[["pmid", "title", "publication_date", "publication_year", "doi"]].drop_duplicates("pmid"), on="pmid", how="left")
    if not evidence.empty and not feature_decisions.empty:
        yes = feature_decisions.loc[feature_decisions["status"].astype(str).eq("completed") & feature_decisions["answer_yes"].fillna(False).astype(bool)]
        supported = yes.groupby(["organism_id", "pmid"])["feature_id"].apply(lambda x: sorted(set(map(str, x)))).to_dict()
        evidence["feature_ids_yes"] = [supported.get((str(r.organism_id), str(r.pmid)), []) for r in evidence.itertuples(index=False)]
    evidence_root = ROOT / "sites/pathogen-atlas-explorer/dist/evidence"
    evidence_root.mkdir(parents=True, exist_ok=True)
    evidence_files = {}
    if not evidence.empty:
        for organism_id, group in evidence.groupby("organism_id"):
            selected = group.sort_values("pmid", key=lambda s: s.astype(str))
            rows = [{key: _safe_json(val) for key, val in row.items()} for row in selected.to_dict("records")]
            shard = hashlib.sha256(str(organism_id).encode()).hexdigest()[:16] + ".json"
            (evidence_root / shard).write_text(json.dumps(rows, ensure_ascii=False, default=str), encoding="utf-8")
            evidence_files[str(organism_id)] = f"evidence/{shard}"
    species_rows = []
    for row in species.to_dict("records"):
        # Full source rows remain in local taxonomy provenance tables; only the
        # normalized hierarchy and display fields go to the public client data.
        item = {key: _safe_json(value) for key, value in row.items() if key != "source_record_json"}
        oid = str(item["organism_id"])
        item["synonyms"] = [name for name in aliases_by_id.get(oid, []) if name.casefold() != str(item.get("accepted_name", "")).casefold()]
        item["feature_values"] = feature_map.get(oid, {})
        item["evidence_file"] = evidence_files.get(oid)
        item["wikipedia_url"] = wikipedia_by_name.get(str(item.get("accepted_name", "")).casefold())
        item["gbif_url"] = gbif_by_name.get(str(item.get("accepted_name", "")).casefold())
        item["taxonomy"] = {part.split("=", 1)[0]: part.split("=", 1)[1] for part in str(item.get("taxonomy_path") or "").split(";") if "=" in part}
        species_rows.append(item)
    stats = {"species_with_profiles": sum(any(v is not None for v in item["feature_values"].values()) for item in species_rows),
             "accepted_species_pmid_records": int(species.get("n_pathogen_pmids", pd.Series(dtype=int)).fillna(0).sum()),
             "pubmed_search_records": int(species.get("n_pubmed_search_records", pd.Series(dtype=int)).fillna(0).sum()),
             "feature_count": len(feature_definitions)}
    # Publish the same hierarchy used by the scoring config so the Explorer can
    # present feature summaries as dimension > category > feature.
    feature_hierarchy = []
    for dimension_id, dimension in catalog["dimensions"].items():
        categories = []
        for category_id, category in dimension.get("categories", {}).items():
            feature_ids = [fid for fid in category.get("features", {})]
            if feature_ids:
                categories.append({"id": category_id, "label": category["label"], "feature_ids": feature_ids})
        if categories:
            feature_hierarchy.append({"id": dimension_id, "label": dimension["label"], "categories": categories})
    payload = {"release": "Pathogen Atlas V3", "species": species_rows, "features": feature_definitions,
               "feature_hierarchy": feature_hierarchy,
               "groups": config["project"]["domains"], "coordinates_ready": bool(not coords.empty), "stats": stats}
    dist = ROOT / "sites/pathogen-atlas-explorer/dist"; dist.mkdir(parents=True, exist_ok=True)
    json_path = dist / "explorer-data.json"
    json_path.write_text(json.dumps(payload, ensure_ascii=False, default=str), encoding="utf-8")
    with json_path.open("rb") as source, gzip.open(dist / "explorer-data.json.gz", "wb", compresslevel=9) as compressed:
        shutil.copyfileobj(source, compressed)
    print(f"Explorer data: {len(species_rows):,} species; {len(feature_definitions):,} features; {len(evidence_files):,} lazy evidence shards.")
    return payload


if __name__ == "__main__": build_site_data()
