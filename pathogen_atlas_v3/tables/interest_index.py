"""Versioned Pathogenic Interest Index calculations."""

from __future__ import annotations

import numpy as np
import pandas as pd

INDEX_VERSION = "pathogenic-interest-index-v3.2-sampled"
REQUIRED = ("organism_id", "organism_group", "n_pathogen_pmids", "n_abstracts_fetched", "n_pubmed_search_records")


def build_pathogenic_interest_index(summary: pd.DataFrame) -> pd.DataFrame:
    """Build group-normalized V3 index and retain every component count.

    Estimate accepted evidence as n_pathogen_pmids / n_abstracts_fetched *
    n_pubmed_search_records, take log, then normalize within organism_group.
    V3.2 uses the configured reproducible per-species abstract sample, matching
    the V2 sampling design; the sample size and seed are stored in Phase 2
    acquisition metadata.
    Groups are intentionally not normalized together.
    """
    missing = set(REQUIRED) - set(summary.columns)
    if missing:
        raise ValueError(f"Missing index inputs: {sorted(missing)}")
    result = summary.loc[:, REQUIRED].copy()
    if result["organism_id"].duplicated().any():
        raise ValueError("Index input must contain one row per organism_id.")
    for name in REQUIRED[2:]:
        result[name] = pd.to_numeric(result[name], errors="coerce").fillna(0)
        if result[name].lt(0).any():
            raise ValueError(f"{name} cannot be negative.")
    if result["n_pathogen_pmids"].gt(result["n_abstracts_fetched"]).any():
        raise ValueError("n_pathogen_pmids cannot exceed n_abstracts_fetched.")
    if result["n_abstracts_fetched"].gt(result["n_pubmed_search_records"]).any():
        raise ValueError("n_abstracts_fetched cannot exceed n_pubmed_search_records.")
    eligible = result["n_pathogen_pmids"].gt(0) & result["n_abstracts_fetched"].gt(0) & result["n_pubmed_search_records"].gt(0)
    result["estimated_pathogen_pmids"] = 0.0
    result.loc[eligible, "estimated_pathogen_pmids"] = (
        result.loc[eligible, "n_pathogen_pmids"] / result.loc[eligible, "n_abstracts_fetched"]
        * result.loc[eligible, "n_pubmed_search_records"]
    )
    result["log_estimated_pathogen_pmids"] = 0.0
    result.loc[eligible, "log_estimated_pathogen_pmids"] = np.log(result.loc[eligible, "estimated_pathogen_pmids"])
    result["pathogenic_interest_index"] = 0.0
    for _, indices in result.groupby("organism_group", sort=False).groups.items():
        valid = [idx for idx in indices if eligible.loc[idx]]
        max_log = float(result.loc[valid, "log_estimated_pathogen_pmids"].max()) if valid else 0.0
        if max_log > 0:
            result.loc[valid, "pathogenic_interest_index"] = result.loc[valid, "log_estimated_pathogen_pmids"] / max_log
    result["index_version"] = INDEX_VERSION
    return result.sort_values(["organism_group", "organism_id"]).reset_index(drop=True)
