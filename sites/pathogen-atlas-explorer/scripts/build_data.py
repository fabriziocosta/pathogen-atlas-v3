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
    reference_root = ROOT / "data/06_references/release"
    wikipedia_path, gbif_path = reference_root / "species_wikipedia_links.csv", reference_root / "species_gbif_links.csv"
    wikipedia = pd.read_csv(wikipedia_path, keep_default_na=False, dtype=str) if wikipedia_path.exists() else pd.DataFrame()
    gbif = pd.read_csv(gbif_path, keep_default_na=False, dtype=str) if gbif_path.exists() else pd.DataFrame()
    wikipedia_by_id = {str(r.organism_id): r.wikipedia_url for r in wikipedia.itertuples(index=False) if "organism_id" in wikipedia and str(getattr(r, "wikipedia_url", "")).strip()} if not wikipedia.empty else {}
    gbif_by_id = {str(r.organism_id): r.gbif_url for r in gbif.itertuples(index=False) if "organism_id" in gbif and str(getattr(r, "gbif_url", "")).strip()} if not gbif.empty else {}
    # Keep V2 associations available until the matching V3 organism group has
    # been refreshed by the new reference notebooks.
    legacy_link_root = ROOT / "data/migration/v2_site_snapshot"
    legacy_wikipedia_path = legacy_link_root / "species_wikipedia_links.csv"
    legacy_gbif_path = legacy_link_root / "species_gbif_links.csv"
    legacy_wikipedia = pd.read_csv(legacy_wikipedia_path) if legacy_wikipedia_path.exists() else pd.DataFrame()
    legacy_gbif = pd.read_csv(legacy_gbif_path) if legacy_gbif_path.exists() else pd.DataFrame()
    wikipedia_by_name = {str(r.canonical_species).casefold(): r.wikipedia_url for r in legacy_wikipedia.itertuples(index=False) if pd.notna(getattr(r, "wikipedia_url", None)) and str(r.wikipedia_url).strip()} if not legacy_wikipedia.empty else {}
    gbif_by_name = {str(r.canonical_species).casefold(): r.gbif_url for r in legacy_gbif.itertuples(index=False) if pd.notna(getattr(r, "gbif_url", None)) and str(r.gbif_url).strip()} if not legacy_gbif.empty else {}
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
            "question": row.get("question", "") or "", "question_claim": row.get("question_claim", "") or "",
            "inclusion_examples": row.get("inclusion_examples", []), "exclusion_examples": row.get("exclusion_examples", []),
            "dimension_id": row["dimension_id"], "dimension": row["dimension_label"],
            "category_id": row["category_id"], "category": row["category_label"],
        }
    # Keep the full, completed Jev gate assessments available for inspection.
    # V3 labels are provisional until the held-out evaluation is approved, so
    # the public Explorer must not hide model confidence behind the accepted
    # evidence table (which is intentionally empty before that approval).
    evidence_path = release / "pathogen_decisions.parquet"
    evidence = pd.read_parquet(evidence_path) if evidence_path.exists() else pd.DataFrame()
    if not evidence.empty:
        evidence = evidence.loc[evidence["status"].astype(str).eq("completed")].copy()
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
        scored = feature_decisions.loc[feature_decisions["status"].astype(str).eq("completed") & feature_decisions["p_yes"].notna()]
        confidences = {
            (str(organism_id), str(pmid)): {
                str(feature_id): float(p_yes) for feature_id, p_yes in zip(group["feature_id"], group["p_yes"])
            }
            for (organism_id, pmid), group in scored.groupby(["organism_id", "pmid"], sort=False)
        }
        evidence["feature_confidences"] = [confidences.get((str(r.organism_id), str(r.pmid)), {}) for r in evidence.itertuples(index=False)]
    evidence_root = ROOT / "sites/pathogen-atlas-explorer/dist/evidence"
    evidence_root.mkdir(parents=True, exist_ok=True)
    evidence_files = {}
    if not evidence.empty:
        for organism_id, group in evidence.groupby("organism_id"):
            selected = group.sort_values("pmid", key=lambda s: s.astype(str))
            public_columns = ["organism_id", "organism_group", "accepted_name", "pmid", "p_yes",
                              "provider", "model", "candidate_pathogen_evidence", "evaluation_approved",
                              "title", "publication_date", "publication_year", "doi",
                              "feature_ids_yes", "feature_confidences"]
            rows = [{key: _safe_json(val) for key, val in row.items() if key in public_columns}
                    for row in selected.to_dict("records")]
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
        name_key = str(item.get("accepted_name", "")).casefold()
        legacy_name_fallback = name_key if str(item.get("organism_group", "")).casefold() == "bacteria" else None
        item["wikipedia_url"] = wikipedia_by_id.get(oid) or (wikipedia_by_name.get(legacy_name_fallback) if legacy_name_fallback else None)
        item["gbif_url"] = gbif_by_id.get(oid) or (gbif_by_name.get(legacy_name_fallback) if legacy_name_fallback else None)
        item["taxonomy"] = {part.split("=", 1)[0]: part.split("=", 1)[1] for part in str(item.get("taxonomy_path") or "").split(";") if "=" in part}
        species_rows.append(item)
    accepted_count = int(species.get("n_pathogen_pmids", pd.Series(dtype=int)).fillna(0).sum())
    candidate_counts: dict[str, int] = {}
    assessment_counts: dict[str, int] = {}
    candidate_record_count = 0
    if not evidence.empty:
        assessment_counts = evidence.groupby("organism_id")["pmid"].nunique().astype(int).to_dict()
        candidates = evidence.loc[evidence["candidate_pathogen_evidence"].fillna(False).astype(bool)]
        candidate_counts = candidates.groupby("organism_id")["pmid"].nunique().astype(int).to_dict()
        candidate_record_count = int(candidates[["organism_id", "pmid"]].drop_duplicates().shape[0])
    for item in species_rows:
        oid = str(item["organism_id"])
        item["n_jev_assessments"] = int(assessment_counts.get(oid, 0))
        item["n_jev_candidate_pmids"] = int(candidate_counts.get(oid, 0))
    feature_status_counts = feature_decisions["status"].astype(str).value_counts().to_dict() if not feature_decisions.empty else {}
    stats = {"species_with_profiles": sum(any(v is not None for v in item["feature_values"].values()) for item in species_rows),
             "accepted_species_pmid_records": accepted_count,
             "jev_completed_assessments": int(evidence[["organism_id", "pmid"]].drop_duplicates().shape[0]) if not evidence.empty else 0,
             "jev_candidate_species_pmid_records": candidate_record_count,
             "feature_assessments_completed": int(feature_status_counts.get("completed", 0)),
             "feature_assessments_failed": int(feature_status_counts.get("failed", 0)),
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
               "feature_hierarchy": feature_hierarchy, "feature_profiles": profiles,
               "groups": config["project"]["domains"], "coordinates_ready": bool(not coords.empty), "stats": stats}
    dist = ROOT / "sites/pathogen-atlas-explorer/dist"; dist.mkdir(parents=True, exist_ok=True)
    json_path = dist / "explorer-data.json"
    json_path.write_text(json.dumps(payload, ensure_ascii=False, default=str), encoding="utf-8")
    with json_path.open("rb") as source, gzip.open(dist / "explorer-data.json.gz", "wb", compresslevel=9) as compressed:
        shutil.copyfileobj(source, compressed)
    print(f"Explorer data: {len(species_rows):,} species; {len(feature_definitions):,} features; {len(evidence_files):,} lazy evidence shards.")
    return payload


if __name__ == "__main__": build_site_data()
