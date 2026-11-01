"""Reproducible shared-core and within-domain IMSTE coordinates."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


def fit_feature_space(feature_matrix: pd.DataFrame, *, shared_feature_ids: list[str], output_dir: str | Path, settings: dict[str, Any] | None = None) -> pd.DataFrame:
    """Fit one shared map plus one map per organism group without zero-imputation."""
    try:
        from imste import IteratedMinimumSpanningTreeEmbedder
        from sklearn.preprocessing import StandardScaler
    except ImportError as exc:
        raise RuntimeError("Install V3 viz dependencies with `pip install -e '.[viz]'` to fit IMSTE.") from exc
    settings = dict(settings or {})
    output = Path(output_dir); output.mkdir(parents=True, exist_ok=True)
    identifiers = ["organism_id", "organism_group"]
    feature_cols = [c for c in feature_matrix.columns if c not in identifiers]
    if not feature_cols:
        raise ValueError("Feature matrix has no feature columns.")
    rows = []
    scopes: list[tuple[str, pd.DataFrame, list[str]]] = []
    shared_cols = [c for c in shared_feature_ids if c in feature_cols]
    if shared_cols:
        scopes.append(("shared_core", feature_matrix, shared_cols))
    for group, frame in feature_matrix.groupby("organism_group", sort=True):
        cols = [c for c in feature_cols if pd.to_numeric(frame[c], errors="coerce").notna().any()]
        if cols:
            scopes.append((str(group), frame, cols))
    for scope, frame, cols in scopes:
        numeric = frame[cols].apply(pd.to_numeric, errors="coerce")
        # Missing observations are imputed to the feature mean, never to zero.
        numeric = numeric.dropna(axis=1, how="all")
        if numeric.shape[0] < 2 or numeric.shape[1] < 1:
            continue
        numeric = numeric.fillna(numeric.mean(axis=0))
        variable = numeric.var(axis=0, ddof=0)
        numeric = numeric.loc[:, variable.gt(0)]
        if numeric.empty:
            continue
        values = StandardScaler().fit_transform(numeric.to_numpy(dtype=float))
        count = len(frame)
        neighbors = min(max(1, int(settings.get("mst_neighbors", 15))), count - 1)
        mapper = IteratedMinimumSpanningTreeEmbedder(
            n_msts=int(settings.get("n_msts", 15)), n_epochs=int(settings.get("n_epochs", 300)),
            mst_method=str(settings.get("mst_method", "famst")), mst_neighbors=neighbors,
            mst_inter_component_edges=int(settings.get("mst_inter_component_edges", 5)),
            mst_max_neighbors=settings.get("mst_max_neighbors"), negative_ratio=int(settings.get("negative_ratio", 5)),
            lambda_rep=float(settings.get("lambda_rep", .5)), learning_rate=float(settings.get("learning_rate", .05)),
            batch_size=int(settings.get("batch_size", 32768)), epsilon=float(settings.get("epsilon", 1e-4)),
            logistic_margin=float(settings.get("logistic_margin", 1)), logistic_temperature=float(settings.get("logistic_temperature", .5)),
            mst_representatives_per_component=int(settings.get("mst_representatives_per_component", 10)),
            random_state=int(settings.get("seed", 42)), device=str(settings.get("device", "auto")),
        )
        coordinates = mapper.fit_transform(values)
        rows.extend({"organism_id": oid, "organism_group": group, "embedding_scope": scope, "x": float(x), "y": float(y), "features_used": list(numeric.columns), "embedding_method": "IMSTE"} for oid, group, x, y in zip(frame.organism_id, frame.organism_group, coordinates[:, 0], coordinates[:, 1]))
        (output / f"{scope}_metadata.json").write_text(json.dumps({"scope": scope, "feature_ids": list(numeric.columns), "imputation": "feature mean within scope; no zero imputation", "seed": int(settings.get("seed", 42))}, indent=2) + "\n", encoding="utf-8")
    result = pd.DataFrame(rows)
    result.to_parquet(output / "coordinates.parquet", index=False)
    return result
