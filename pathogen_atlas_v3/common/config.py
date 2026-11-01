"""Validated helpers for active V3 run configuration."""

from __future__ import annotations

from typing import Any


def selected_organism_groups(config: dict[str, Any]) -> list[str]:
    """Return the requested organism groups, defaulting to every supported group."""
    supported = [str(value) for value in config.get("project", {}).get("domains", [])]
    selected = config.get("run", {}).get("organism_groups", supported)
    if isinstance(selected, str):
        selected = [selected]
    selected = [str(value) for value in selected]
    unknown = sorted(set(selected) - set(supported))
    if unknown:
        raise ValueError(f"Unknown organism_groups {unknown}; supported groups are {supported}.")
    if not selected:
        raise ValueError("run.organism_groups must select at least one supported group.")
    if len(selected) != len(set(selected)):
        raise ValueError("run.organism_groups must not contain duplicates.")
    return selected
