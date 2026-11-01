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
    index_path = release / "pathogenic_interest_index.parquet"
    feature_path = release / "feature_summary.parquet"
    coordinate_path = ROOT / paths["embedding_root"] / "coordinates.parquet"
    if not species_path.exists():
        raise FileNotFoundError("Run Phase 4 before building Explorer data.")
    species = pd.read_parquet(species_path)
    if not index_path.exists():
        raise FileNotFoundError("Run Phase 4 and notebooks/08_pathogenic_interest_index.ipynb before building Explorer data.")
    index = pd.read_parquet(index_path)
    index_columns = [column for column in (
        "organism_id", "organism_group", "n_pathogen_pmids", "n_abstracts_fetched",
        "n_pubmed_search_records", "estimated_pathogen_pmids",
        "log_estimated_pathogen_pmids", "pathogenic_interest_index", "index_version", "index_source",
    ) if column in index.columns]
    index_by_id = {
        str(row["organism_id"]): {key: _safe_json(value) for key, value in row.items() if key != "organism_id"}
        for row in index[index_columns].to_dict("records")
    }
    features = pd.read_parquet(feature_path) if feature_path.exists() else pd.DataFrame()
    coords = pd.read_parquet(coordinate_path) if coordinate_path.exists() else pd.DataFrame()
    v2_release_root = ROOT / "data/migration/version2_snapshot/04_tables/release"
    v2_species_path = v2_release_root / "species_summary.parquet"
    v2_evidence_path = v2_release_root / "evidence_records.parquet"
    v2_feature_matrix_path = v2_release_root / "species_feature_matrix.parquet"
    v2_feature_summary_path = v2_release_root / "species_feature_summary.parquet"
    v2_answers_path = v2_release_root / "domain_expert_answers.parquet"
    v2_coordinates_path = ROOT / "data/migration/version2_snapshot/05_plots/release/species_feature_space_coordinates.parquet"
    v2_species = pd.read_parquet(v2_species_path) if v2_species_path.exists() else pd.DataFrame()
    v2_evidence = pd.read_parquet(v2_evidence_path) if v2_evidence_path.exists() else pd.DataFrame()
    v2_feature_matrix = pd.read_parquet(v2_feature_matrix_path) if v2_feature_matrix_path.exists() else pd.DataFrame()
    v2_feature_summary = pd.read_parquet(v2_feature_summary_path) if v2_feature_summary_path.exists() else pd.DataFrame()
    v2_answers = pd.read_parquet(v2_answers_path) if v2_answers_path.exists() else pd.DataFrame()
    v2_coordinates = pd.read_parquet(v2_coordinates_path) if v2_coordinates_path.exists() else pd.DataFrame()
    v2_summary_by_name = {
        str(row.canonical_species).casefold(): row._asdict()
        for row in v2_species.itertuples(index=False)
    } if not v2_species.empty else {}
    bacteria_organism_by_name = {
        str(row.accepted_name).casefold(): str(row.organism_id)
        for row in species.loc[species["organism_group"].eq("bacteria")].itertuples(index=False)
    }
    v2_features_by_name = {}
    if not v2_feature_matrix.empty:
        feature_columns = [column for column in v2_feature_matrix.columns if column.startswith("feature_")]
        for row in v2_feature_matrix.itertuples(index=False):
            v2_features_by_name[str(row.canonical_species).casefold()] = {
                column.removeprefix("feature_"): (None if pd.isna(getattr(row, column)) else float(getattr(row, column)))
                for column in feature_columns
            }
    v2_feature_counts_by_name: dict[str, dict] = {}
    if not v2_feature_summary.empty:
        for row in v2_feature_summary.itertuples(index=False):
            v2_feature_counts_by_name.setdefault(str(row.canonical_species).casefold(), {})[str(row.question_id)] = {
                "n_yes": int(row.n_yes) if pd.notna(row.n_yes) else 0,
                "n_scored": int(row.n_answered) if pd.notna(row.n_answered) else 0,
            }
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
    if not v2_coordinates.empty:
        v2_coordinates["organism_id"] = v2_coordinates["canonical_species"].astype(str).str.casefold().map(bacteria_organism_by_name)
        v2_coordinates = v2_coordinates.loc[v2_coordinates["organism_id"].notna()].drop_duplicates("organism_id")
        legacy_x = v2_coordinates.set_index("organism_id")["imste_1"].to_dict()
        legacy_y = v2_coordinates.set_index("organism_id")["imste_2"].to_dict()
        species["x"] = species["x"].where(species["x"].notna(), species["organism_id"].map(legacy_x))
        species["y"] = species["y"].where(species["y"].notna(), species["organism_id"].map(legacy_y))
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
    feature_counts_map: dict[str, dict] = {}
    feature_definitions = {}
    if not features.empty:
        for row in features.itertuples(index=False):
            organism_id, feature_id = str(row.organism_id), str(row.feature_id)
            feature_map.setdefault(organism_id, {})[feature_id] = None if pd.isna(row.share_yes) else float(row.share_yes)
            feature_counts_map.setdefault(organism_id, {})[feature_id] = {
                "n_yes": int(row.n_yes) if pd.notna(row.n_yes) else 0,
                "n_scored": int(row.n_scored) if pd.notna(row.n_scored) else 0,
            }
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
    # The Explorer displays accepted human-pathogen evidence only.
    evidence_path = release / "pathogen_decisions.parquet"
    evidence = pd.read_parquet(evidence_path) if evidence_path.exists() else pd.DataFrame()
    completed_gate_count = 0
    if not evidence.empty:
        evidence = evidence.loc[evidence["status"].astype(str).eq("completed")].copy()
        completed_gate_count = int(evidence[["organism_id", "pmid"]].drop_duplicates().shape[0])
        # Keep the same pathogen-only evidence policy across release sources.
        evidence = evidence.loc[evidence["accepted_pathogen_evidence"].fillna(False).astype(bool)].copy()
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
    if not v2_evidence.empty:
        legacy_evidence = v2_evidence.loc[v2_evidence["accepted"].fillna(False).astype(bool)].copy()
        legacy_evidence["organism_id"] = legacy_evidence["canonical_species"].astype(str).str.casefold().map(bacteria_organism_by_name)
        legacy_evidence = legacy_evidence.loc[legacy_evidence["organism_id"].notna()].copy()
        legacy_evidence["organism_group"] = "bacteria"
        legacy_evidence["accepted_name"] = legacy_evidence["canonical_species"].astype(str)
        legacy_evidence["pmid"] = legacy_evidence["pmid"].astype(str)
        legacy_evidence["p_yes"] = pd.to_numeric(legacy_evidence["p_yes"], errors="coerce")
        legacy_evidence["candidate_pathogen_evidence"] = True
        legacy_evidence["accepted_pathogen_evidence"] = True
        legacy_evidence["evaluation_approved"] = True
        legacy_evidence["title"] = legacy_evidence["title"].fillna("").astype(str).str.strip()
        missing_title_pmids = set(legacy_evidence.loc[legacy_evidence["title"].eq(""), "pmid"].astype(str))
        raw_articles_path = ROOT / "data/migration/version2_snapshot/02_abstracts/raw/pubmed_article_records.jsonl"
        if missing_title_pmids and raw_articles_path.exists():
            recovered_titles = {}
            with raw_articles_path.open(encoding="utf-8") as raw_articles:
                for line in raw_articles:
                    article = json.loads(line)
                    pmid = str(article.get("PMID", ""))
                    if pmid not in missing_title_pmids:
                        continue
                    title_parts = article.get("TI") or article.get("BTI") or article.get("CTI") or []
                    if isinstance(title_parts, str):
                        title_parts = [title_parts]
                    recovered_titles[pmid] = " ".join(map(str, title_parts)).strip()
                    missing_title_pmids.discard(pmid)
                    if not missing_title_pmids:
                        break
            legacy_evidence["title"] = [title or recovered_titles.get(pmid, "")
                                         for title, pmid in zip(legacy_evidence["title"], legacy_evidence["pmid"].astype(str))]
        if not v2_answers.empty:
            scored_legacy = v2_answers.loc[v2_answers["status"].astype(str).eq("completed") & v2_answers["p_yes"].notna()].copy()
            scored_legacy["organism_id"] = scored_legacy["canonical_species"].astype(str).str.casefold().map(bacteria_organism_by_name)
            scored_legacy["pmid"] = scored_legacy["pmid"].astype(str)
            scored_legacy = scored_legacy.loc[scored_legacy["organism_id"].notna()].copy()
            legacy_yes = scored_legacy.loc[scored_legacy["answer_yes"].fillna(False).astype(bool)]
            legacy_supported = legacy_yes.groupby(["organism_id", "pmid"])["question_id"].apply(lambda x: sorted(set(map(str, x)))).to_dict()
            legacy_confidence = {
                (str(organism_id), str(pmid)): {
                    str(feature_id): float(p_yes) for feature_id, p_yes in zip(group["question_id"], group["p_yes"])
                }
                for (organism_id, pmid), group in scored_legacy.groupby(["organism_id", "pmid"], sort=False)
            }
            legacy_evidence["feature_ids_yes"] = [legacy_supported.get((str(row.organism_id), str(row.pmid)), []) for row in legacy_evidence.itertuples(index=False)]
            legacy_evidence["feature_confidences"] = [legacy_confidence.get((str(row.organism_id), str(row.pmid)), {}) for row in legacy_evidence.itertuples(index=False)]
        else:
            legacy_evidence["feature_ids_yes"] = [[] for _ in range(len(legacy_evidence))]
            legacy_evidence["feature_confidences"] = [{} for _ in range(len(legacy_evidence))]
        evidence = pd.concat([evidence, legacy_evidence], ignore_index=True, sort=False)
    evidence_root = ROOT / "sites/pathogen-atlas-explorer/dist/evidence"
    evidence_root.mkdir(parents=True, exist_ok=True)
    # The accepted pathogen cohort can shrink when its release threshold changes.
    # Remove stale shards so the release archive contains only linked evidence.
    for stale_shard in evidence_root.glob("*.json"):
        stale_shard.unlink()
    evidence_files = {}
    if not evidence.empty:
        for organism_id, group in evidence.groupby("organism_id"):
            selected = group.sort_values("pmid", key=lambda s: s.astype(str))
            public_columns = ["organism_id", "organism_group", "accepted_name", "pmid", "p_yes",
                              "provider", "model", "candidate_pathogen_evidence", "accepted_pathogen_evidence", "evaluation_approved",
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
        item["feature_counts"] = feature_counts_map.get(oid, {})
        item["evidence_file"] = evidence_files.get(oid)
        name_key = str(item.get("accepted_name", "")).casefold()
        # The released index table is the source of both the score and its
        # supporting counts. It normalizes within each organism group.
        index_row = index_by_id.get(oid)
        if index_row and str(index_row.get("organism_group")) == str(item.get("organism_group")):
            item.update(index_row)
        if item.get("organism_group") == "bacteria" and name_key in v2_summary_by_name:
            legacy = v2_summary_by_name[name_key]
            item["feature_values"] = v2_features_by_name.get(name_key, {})
            item["feature_counts"] = v2_feature_counts_by_name.get(name_key, {})
        legacy_name_fallback = name_key if str(item.get("organism_group", "")).casefold() == "bacteria" else None
        item["wikipedia_url"] = wikipedia_by_id.get(oid) or (wikipedia_by_name.get(legacy_name_fallback) if legacy_name_fallback else None)
        item["gbif_url"] = gbif_by_id.get(oid) or (gbif_by_name.get(legacy_name_fallback) if legacy_name_fallback else None)
        item["taxonomy"] = {part.split("=", 1)[0]: part.split("=", 1)[1] for part in str(item.get("taxonomy_path") or "").split(";") if "=" in part}
        species_rows.append(item)
    accepted_counts: dict[str, int] = {}
    assessment_counts: dict[str, int] = {}
    if not evidence.empty:
        assessment_counts = evidence.groupby("organism_id")["pmid"].nunique().astype(int).to_dict()
        accepted_counts = assessment_counts
    for item in species_rows:
        oid = str(item["organism_id"])
        item["n_jev_assessments"] = int(assessment_counts.get(oid, 0))
        item["n_jev_accepted_pmids"] = int(accepted_counts.get(oid, 0))
    # Keep only human-pathogen profiles: the established V2 evidence cohort for
    # bacteria and the accepted V3 evidence cohort for viruses.
    species_rows = [item for item in species_rows if int(item.get("n_pathogen_pmids", 0) or 0) > 0]
    accepted_count = sum(int(item.get("n_pathogen_pmids", 0) or 0) for item in species_rows)
    v3_feature_status = feature_decisions["status"].astype(str).value_counts().to_dict() if not feature_decisions.empty else {}
    v2_feature_status = v2_answers["status"].astype(str).value_counts().to_dict() if not v2_answers.empty else {}
    feature_status_counts = {status: int(v3_feature_status.get(status, 0)) + int(v2_feature_status.get(status, 0))
                             for status in set(v3_feature_status) | set(v2_feature_status)}
    stats = {"species_with_profiles": sum(any(v is not None for v in item["feature_values"].values()) for item in species_rows),
             "accepted_species_pmid_records": accepted_count,
             "jev_completed_assessments": completed_gate_count,
             "jev_accepted_species_pmid_records": int(evidence[["organism_id", "pmid"]].drop_duplicates().shape[0]) if not evidence.empty else 0,
             "bacteria_pathogen_species": sum(item["organism_group"] == "bacteria" for item in species_rows),
             "virus_pathogen_species": sum(item["organism_group"] == "viruses" for item in species_rows),
             "feature_assessments_completed": int(feature_status_counts.get("completed", 0)),
             "feature_assessments_failed": int(feature_status_counts.get("failed", 0)),
             "pubmed_search_records": sum(int(item.get("n_pubmed_search_records", 0) or 0) for item in species_rows),
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
    coordinates_ready = bool(species_rows) and all(pd.notna(item.get("x")) and pd.notna(item.get("y")) for item in species_rows)
    payload = {"release": "Pathogen Atlas V3", "species": species_rows, "features": feature_definitions,
               "feature_hierarchy": feature_hierarchy, "feature_profiles": profiles,
               "groups": sorted({item["organism_group"] for item in species_rows}),
               "coordinates_ready": coordinates_ready, "stats": stats}
    dist = ROOT / "sites/pathogen-atlas-explorer/dist"; dist.mkdir(parents=True, exist_ok=True)
    json_path = dist / "explorer-data.json"
    json_path.write_text(json.dumps(payload, ensure_ascii=False, default=str), encoding="utf-8")
    with json_path.open("rb") as source, gzip.open(dist / "explorer-data.json.gz", "wb", compresslevel=9) as compressed:
        shutil.copyfileobj(source, compressed)
    print(f"Explorer data: {len(species_rows):,} species; {len(feature_definitions):,} features; {len(evidence_files):,} lazy evidence shards.")
    return payload


if __name__ == "__main__": build_site_data()
