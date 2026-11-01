"""Small, explicit I/O helpers for configuration and durable artifacts."""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any

import pandas as pd
import yaml


def load_yaml(path: str | Path) -> dict[str, Any]:
    with Path(path).open(encoding="utf-8") as handle:
        value = yaml.safe_load(handle) or {}
    if not isinstance(value, dict):
        raise ValueError(f"Expected a YAML mapping in {path}.")
    return value


def atomic_parquet(frame: pd.DataFrame, path: str | Path) -> None:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_name(f".{target.name}.tmp")
    frame.to_parquet(temporary, index=False)
    os.replace(temporary, target)


def append_jsonl(path: str | Path, record: dict[str, Any]) -> None:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(record, ensure_ascii=False, sort_keys=True, default=str))
        handle.write("\n")


def read_parquet_or_empty(path: str | Path, columns: list[str] | tuple[str, ...]) -> pd.DataFrame:
    target = Path(path)
    if not target.exists():
        return pd.DataFrame(columns=list(columns))
    return pd.read_parquet(target)
