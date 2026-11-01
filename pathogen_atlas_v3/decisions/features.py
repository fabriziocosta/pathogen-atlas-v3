"""Load and validate the master V3 feature catalogue and domain profiles."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

from pathogen_atlas_v3.common.schemas import validate_domain


def load_feature_catalog(catalog_path: str | Path, profile_path: str | Path, *, allow_provisional: bool = False) -> tuple[dict[str, Any], list[str]]:
    catalog = yaml.safe_load(Path(catalog_path).read_text(encoding="utf-8")) or {}
    profiles = yaml.safe_load(Path(profile_path).read_text(encoding="utf-8")) or {}
    if catalog.get("version") is None or catalog.get("schema_version") != 2:
        raise ValueError("Feature catalogue needs a version and hierarchical schema_version: 2.")
    dimensions = catalog.get("dimensions", {})
    if not dimensions:
        raise ValueError("Feature catalogue needs at least one top-level dimension.")
    features: dict[str, dict[str, Any]] = {}
    for dimension_id, dimension in dimensions.items():
        for category_id, category in dimension.get("categories", {}).items():
            for feature_id, source in category.get("features", {}).items():
                if feature_id in features:
                    raise ValueError(f"Feature ID {feature_id!r} appears more than once in the hierarchy.")
                if source.get("feature_id", feature_id) != feature_id:
                    raise ValueError(f"Feature key {feature_id!r} does not match its feature_id.")
                feature = {
                    **source,
                    "feature_id": source.get("feature_id", feature_id),
                    "dimension_id": dimension_id,
                    "dimension_label": dimension.get("label", dimension_id),
                    "category_id": category_id,
                    "category_label": category.get("label", category_id),
                    # Retain the legacy category field for release tables and
                    # downstream readers while carrying the full hierarchy.
                    "category": category_id,
                }
                features[feature_id] = feature
    if not features:
        raise ValueError("Feature catalogue hierarchy has no features.")
    # Keep the nested source of truth and expose a normalized flat lookup to the
    # scoring pipeline, where feature IDs are used as stable cache keys.
    catalog["features"] = features
    for feature_id, feature in features.items():
        required = {"feature_id", "definition", "question", "question_claim", "inclusion_examples", "exclusion_examples", "question_version", "organism_groups", "answer_type", "instructions", "criteria", "review_status", "dimension_id", "category_id"}
        if required - set(feature):
            raise ValueError(f"Feature {feature_id!r} missing fields {sorted(required - set(feature))}.")
        if not feature["inclusion_examples"] or not feature["exclusion_examples"]:
            raise ValueError(f"Feature {feature_id!r} needs at least one inclusion and exclusion example.")
        if feature["answer_type"] != "noul":
            raise ValueError(f"Feature {feature_id!r} must use Jev noul type.")
    return catalog, profiles.get("profiles", {})


def questions_for_domain(catalog: dict[str, Any], profiles: dict[str, list[str]], domain: str, *, allow_provisional: bool = False) -> dict[str, dict[str, Any]]:
    group = validate_domain(domain)
    ordered_ids = profiles.get(group)
    if not isinstance(ordered_ids, list):
        raise ValueError(f"No ordered feature profile for {group}.")
    features = catalog["features"]
    result: dict[str, dict[str, Any]] = {}
    for feature_id in ordered_ids:
        feature = features.get(feature_id)
        if not feature or group not in feature["organism_groups"]:
            raise ValueError(f"Profile {group} contains an invalid or inapplicable feature {feature_id!r}.")
        if not allow_provisional and feature["review_status"] != "approved":
            continue
        contextual_rule = {
            "bacteria": "Assess the named bacterium only; distinguish disease attribution from isolation, carriage, and co-detection.",
            "fungi": "Assess the named fungus only; distinguish infection from colonisation, environmental contamination, allergy, and toxin-only illness.",
            "viruses": "Assess the named virus only; distinguish active disease from exposure, serology, environmental detection, and mention of a related strain or subtype.",
        }[group]
        includes = " ".join(feature["inclusion_examples"])
        excludes = " ".join(feature["exclusion_examples"])
        instruction = (
            f"{contextual_rule} {feature['instructions']} "
            f"Apply this exact question: {feature['question']} "
            f"Inclusion example: {includes} Exclusion example: {excludes} "
            "Examples are schematic and not exhaustive. Use only the title and abstract; do not infer missing evidence."
        )
        result[feature_id] = {**feature, "type": feature["answer_type"], "instructions": instruction}
    return result
