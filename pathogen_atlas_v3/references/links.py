"""Resumable Wikipedia and GBIF link resolution for V3 taxonomy records."""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Callable
from urllib.error import HTTPError, URLError
from urllib.parse import quote, urlencode
from urllib.request import Request, urlopen

import pandas as pd


GROUP_TO_KINGDOM = {"bacteria": "Bacteria", "fungi": "Fungi", "viruses": "Viruses"}
WIKIPEDIA_FIELDS = ["organism_id", "organism_group", "accepted_name", "wikipedia_title", "wikipedia_url", "match_status", "match_type"]
GBIF_FIELDS = ["organism_id", "organism_group", "accepted_name", "gbif_taxon_key", "gbif_scientific_name", "gbif_url", "match_status", "match_type", "confidence"]


def _request_json(request: Request, *, timeout: int, retries: int = 7) -> object:
    for attempt in range(retries):
        try:
            with urlopen(request, timeout=timeout) as response:
                return json.load(response)
        except (HTTPError, URLError, TimeoutError) as error:
            if attempt == retries - 1:
                raise RuntimeError(f"Reference lookup request failed: {error}") from error
            retry_after = error.headers.get("Retry-After") if isinstance(error, HTTPError) else None
            delay = float(retry_after) if retry_after and retry_after.isdigit() else min(2 * (2**attempt), 60)
            time.sleep(delay)
    raise AssertionError("unreachable")


def _wikipedia_batch(names: list[str]) -> dict[str, dict[str, str]]:
    query = urlencode({
        "action": "query", "format": "json", "formatversion": "2", "redirects": "1",
        "maxlag": "5", "titles": "|".join(names),
    })
    request = Request(
        f"https://en.wikipedia.org/w/api.php?{query}",
        headers={"User-Agent": "PathogenAtlasV3/1.0 (taxonomy reference link resolver)"},
    )
    payload = _request_json(request, timeout=30)
    if not isinstance(payload, dict) or payload.get("error"):
        raise RuntimeError(f"Wikipedia API returned an error: {payload.get('error') if isinstance(payload, dict) else payload}")
    query_data = payload.get("query", {})
    title_map: dict[str, str] = {}
    for mapping_type in ("normalized", "redirects"):
        for mapping in query_data.get(mapping_type, []):
            title_map[mapping["from"].casefold()] = mapping["to"]
    found_by_title = {
        page["title"].casefold(): page["title"]
        for page in query_data.get("pages", []) if "missing" not in page
    }
    resolved = {}
    for requested in names:
        title = requested
        visited: set[str] = set()
        while title.casefold() in title_map and title.casefold() not in visited:
            visited.add(title.casefold())
            title = title_map[title.casefold()]
        canonical_title = found_by_title.get(title.casefold())
        if canonical_title:
            url = "https://en.wikipedia.org/wiki/" + quote(canonical_title.replace(" ", "_"), safe="()_,-.")
            resolved[requested] = {"wikipedia_title": canonical_title, "wikipedia_url": url, "match_status": "found", "match_type": "exact_or_redirect"}
        else:
            resolved[requested] = {"wikipedia_title": "", "wikipedia_url": "", "match_status": "not_found", "match_type": "exact_title_query"}
    return resolved


def _gbif_batch(records: list[dict[str, str]]) -> dict[str, dict[str, str | int | float]]:
    body = json.dumps([
        {"scientificName": row["accepted_name"], "kingdom": GROUP_TO_KINGDOM[row["organism_group"]]}
        for row in records
    ]).encode("utf-8")
    request = Request(
        "https://api.gbif.org/v2/species/match", data=body,
        headers={"Accept": "application/json", "Content-Type": "application/json", "User-Agent": "PathogenAtlasV3/1.0 (GBIF reference link resolver)"},
        method="POST",
    )
    payload = _request_json(request, timeout=60)
    if not isinstance(payload, list) or len(payload) != len(records):
        raise RuntimeError(f"GBIF returned {len(payload) if isinstance(payload, list) else 'an invalid number of'} results for {len(records)} names.")
    matches = {}
    for row, result in zip(records, payload):
        usage = result.get("usage") or {}
        diagnostics = result.get("diagnostics") or {}
        match_type = str(diagnostics.get("matchType") or "NONE")
        canonical_name = str(usage.get("canonicalName") or "")
        key = str(usage.get("key") or "")
        exact = bool(key) and str(usage.get("rank") or "").upper() == "SPECIES" and match_type == "EXACT" and canonical_name.casefold() == row["accepted_name"].casefold()
        matches[row["organism_id"]] = {
            "gbif_taxon_key": key if exact else "",
            "gbif_scientific_name": str(usage.get("name") or "") if exact else "",
            "gbif_url": f"https://www.gbif.org/species/{quote(key, safe='')}" if exact else "",
            "match_status": "found" if exact else "not_found",
            "match_type": match_type,
            "confidence": diagnostics.get("confidence", ""),
        }
    return matches


def _selected_records(species_path: str | Path, organism_groups: list[str], organism_ids: set[str] | None = None) -> list[dict[str, str]]:
    unknown = sorted(set(organism_groups) - set(GROUP_TO_KINGDOM))
    if unknown:
        raise ValueError(f"Unsupported organism groups: {unknown}")
    frame = pd.read_parquet(species_path, columns=["organism_id", "organism_group", "accepted_name", "rank", "mapping_status"])
    frame = frame.loc[frame.organism_group.isin(organism_groups) & frame["rank"].astype(str).str.casefold().eq("species")]
    if organism_ids is not None:
        frame = frame.loc[frame.organism_id.astype(str).isin(organism_ids)]
    frame = frame.loc[frame.accepted_name.notna() & frame.organism_id.notna()]
    frame = frame.drop_duplicates("organism_id")
    return [
        {"organism_id": str(row.organism_id), "organism_group": str(row.organism_group), "accepted_name": str(row.accepted_name).strip()}
        for row in frame.itertuples(index=False) if str(row.accepted_name).strip()
    ]


def _write(rows: dict[str, dict], path: Path, fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    output = pd.DataFrame((rows[key] for key in sorted(rows)), columns=fields)
    output.to_csv(path, index=False)


def _resolve(
    species_path: str | Path, output_path: str | Path, groups: list[str], *, force_refresh: bool,
    fields: list[str], batch_size: int, delay: float, lookup: Callable,
    organism_ids: set[str] | None = None,
) -> pd.DataFrame:
    records = _selected_records(species_path, groups, organism_ids)
    if not records:
        raise ValueError(f"No species-rank records found for {groups} in {species_path}")
    path = Path(output_path)
    rows: dict[str, dict] = {}
    selected_ids = {row["organism_id"] for row in records}
    selected_groups = set(groups)
    if path.exists():
        previous = pd.read_csv(path, keep_default_na=False, dtype=str)
        if set(fields).issubset(previous.columns):
            for item in previous.to_dict("records"):
                oid = str(item["organism_id"])
                outside_scope = item.get("organism_group") not in selected_groups
                selected_cache = oid in selected_ids and not force_refresh
                if item.get("match_status") in {"found", "not_found"} and (outside_scope or selected_cache):
                    rows[oid] = {field: item[field] for field in fields}
    pending = [record for record in records if record["organism_id"] not in rows]
    print(f"Resolving {len(pending):,} records; reusing {len(rows):,} completed rows.", flush=True)
    for start in range(0, len(pending), batch_size):
        batch = pending[start:start + batch_size]
        if fields is WIKIPEDIA_FIELDS:
            names = list(dict.fromkeys(row["accepted_name"] for row in batch))
            results = lookup(names)
            for record in batch:
                rows[record["organism_id"]] = {**record, **results[record["accepted_name"]]}
        else:
            results = lookup(batch)
            for record in batch:
                rows[record["organism_id"]] = {**record, **results[record["organism_id"]]}
        _write(rows, path, fields)
        print(f"Resolved {min(start + len(batch), len(pending)):,}/{len(pending):,} pending records.", flush=True)
        if delay:
            time.sleep(delay)
    _write(rows, path, fields)
    print(f"Saved {len(rows):,} associations to {path}", flush=True)
    return pd.read_csv(path, keep_default_na=False)


def resolve_wikipedia_links(species_path: str | Path, output_path: str | Path, organism_groups: list[str], *, force_refresh: bool = False, organism_ids: set[str] | None = None) -> pd.DataFrame:
    return _resolve(species_path, output_path, organism_groups, force_refresh=force_refresh, fields=WIKIPEDIA_FIELDS, batch_size=25, delay=1.1, lookup=_wikipedia_batch, organism_ids=organism_ids)


def resolve_gbif_links(species_path: str | Path, output_path: str | Path, organism_groups: list[str], *, force_refresh: bool = False, organism_ids: set[str] | None = None) -> pd.DataFrame:
    return _resolve(species_path, output_path, organism_groups, force_refresh=force_refresh, fields=GBIF_FIELDS, batch_size=100, delay=0.25, lookup=_gbif_batch, organism_ids=organism_ids)
