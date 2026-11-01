"""Export the versioned V3 feature hierarchy for the static Explorer."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[3]


def build_feature_catalog(output: str | Path) -> dict:
    catalog = yaml.safe_load((ROOT / "configs/features.yaml").read_text(encoding="utf-8"))
    profile_data = yaml.safe_load((ROOT / "configs/feature_profiles.yaml").read_text(encoding="utf-8"))
    dimensions = catalog["dimensions"]
    features: dict[str, dict] = {}
    hierarchy: list[dict] = []

    for dimension_id, dimension in dimensions.items():
        categories = []
        for category_id, category in dimension.get("categories", {}).items():
            feature_ids = list(category.get("features", {}))
            for feature_id, feature in category.get("features", {}).items():
                features[feature_id] = {
                    "label": feature.get("label", feature_id.replace("_", " ").upper()),
                    "definition": feature.get("definition", "") or "",
                    "question": feature.get("question", "") or "",
                    "question_claim": feature.get("question_claim", "") or "",
                    "inclusion_examples": feature.get("inclusion_examples", []),
                    "exclusion_examples": feature.get("exclusion_examples", []),
                    "dimension_id": dimension_id,
                    "dimension": dimension.get("label", dimension_id),
                    "category_id": category_id,
                    "category": category.get("label", category_id),
                }
            if feature_ids:
                categories.append({
                    "id": category_id,
                    "label": category.get("label", category_id),
                    "feature_ids": feature_ids,
                })
        if categories:
            hierarchy.append({
                "id": dimension_id,
                "label": dimension.get("label", dimension_id),
                "categories": categories,
            })

    result = {
        "catalog_version": catalog["version"],
        "profile_version": profile_data["version"],
        "features": features,
        "feature_hierarchy": hierarchy,
        "feature_profiles": profile_data["profiles"],
    }
    destination = Path(output)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(result, ensure_ascii=False, separators=(",", ":")), encoding="utf-8")
    print(f"Explorer feature catalogue: {len(features):,} features across {len(hierarchy)} dimensions → {destination}")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="sites/pathogen-atlas-explorer/dist/feature-catalog.json")
    args = parser.parse_args()
    build_feature_catalog(ROOT / args.output)


if __name__ == "__main__":
    main()
