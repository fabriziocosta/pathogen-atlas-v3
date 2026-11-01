"""Load and validate the master V3 feature catalogue and domain profiles."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

from pathogen_atlas_v3.common.schemas import validate_domain


def load_feature_catalog(catalog_path: str | Path, profile_path: str | Path, *, allow_provisional: bool = False) -> tuple[dict[str, Any], list[str]]:
    catalog = yaml.safe_load(Path(catalog_path).read_text(encoding="utf-8")) or {}
    profiles = yaml.safe_load(Path(profile_path).read_text(encoding="utf-8")) or {}
    features = catalog.get("features", {})
    if not features or catalog.get("version") is None:
        raise ValueError("Feature catalogue needs a version and non-empty features map.")
    for feature_id, feature in features.items():
        required = {"feature_id", "definition", "question_version", "organism_groups", "answer_type", "instructions", "criteria", "review_status"}
        if required - set(feature):
            raise ValueError(f"Feature {feature_id!r} missing fields {sorted(required - set(feature))}.")
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
        instruction = f"{contextual_rule} {feature['instructions']}"
        result[feature_id] = {**feature, "type": feature["answer_type"], "instructions": instruction}
    return result
