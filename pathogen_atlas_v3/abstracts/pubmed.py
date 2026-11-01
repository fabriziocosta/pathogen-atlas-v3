"""Rate-limited, concurrent, resumable PubMed acquisition via NCBI E-utilities."""

from __future__ import annotations

import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import threading
import time
import xml.etree.ElementTree as ET
from datetime import date, timedelta
from urllib.parse import urlencode
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from pathlib import Path
from typing import Any

import pandas as pd
import requests

from pathogen_atlas_v3.common.io import append_jsonl, atomic_parquet

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
_session_local = threading.local()


def resolve_ncbi_email(email: str | None = None) -> str:
    """Use NCBI_EMAIL first, then the configured local Git author contact."""
    value = (email or os.environ.get("NCBI_EMAIL", "")).strip()
    if value:
        return value
    if shutil.which("git"):
        try:
            result = subprocess.run(["git", "config", "user.email"], capture_output=True, text=True, check=False, timeout=5)
        except (OSError, subprocess.TimeoutExpired):
            return ""
        value = result.stdout.strip() if result.returncode == 0 else ""
    return value if "@" in value else ""


def _ncbi_key_from_keychain() -> str | None:
    """Read the optional NCBI key from the AMR-compatible macOS Keychain item."""
    if sys.platform != "darwin" or shutil.which("security") is None:
        return None
    try:
        result = subprocess.run(
            ["security", "find-generic-password", "-s", "NCBI", "-a", "NCBI_API_KEY", "-w"],
            capture_output=True, text=True, check=False, timeout=10,
        )
    except (OSError, subprocess.TimeoutExpired):
        return None
    key = result.stdout.strip()
    return key if result.returncode == 0 and key else None


class RateLimiter:
    """Coordinate request starts across worker threads to honor NCBI limits."""
    def __init__(self, delay: float):
        self.delay = max(0.0, float(delay))
        self.lock = threading.Lock()
        self.next_time = 0.0

    def wait(self) -> None:
        with self.lock:
            now = time.monotonic()
            sleep_for = max(0.0, self.next_time - now)
            self.next_time = max(now, self.next_time) + self.delay
        if sleep_for:
            time.sleep(sleep_for)

    def cooldown(self, seconds: float) -> None:
        """Pause all workers after the provider reports rate limiting."""
        with self.lock:
            self.next_time = max(self.next_time, time.monotonic() + max(0.0, float(seconds)))


def build_name_query(names: list[str]) -> str:
    clean = sorted({normalized for name in names if (normalized := " ".join(str(name).split())) and len(normalized) >= 3 and not re.fullmatch(r"\d+(?:\.\d+)?", normalized)}, key=str.casefold)
    if not clean:
        raise ValueError("At least one accepted or synonym name is required.")
    quoted = " OR ".join(f'"{name.replace(chr(34), "")}"[Title/Abstract]' for name in clean)
    human_context = '(humans[MeSH Terms] OR human[Title/Abstract] OR patient[Title/Abstract] OR patients[Title/Abstract] OR infection[Title/Abstract] OR disease[Title/Abstract])'
    return f"({quoted}) AND {human_context}"


def _thread_session() -> requests.Session:
    if not hasattr(_session_local, "value"):
        _session_local.value = requests.Session()
    return _session_local.value


def _request(url: str, params: dict[str, Any], *, email: str, api_key: str | None, limiter: RateLimiter, requester=None, retries: int = 3) -> requests.Response:
    query = {**params, "tool": "pathogen_atlas_v3", "email": email}
    if api_key:
        query["api_key"] = api_key
    encoded_length = len(url) + 1 + len(urlencode(query, doseq=True))
    use_post = requester is None and encoded_length > 1500
    for attempt in range(retries + 1):
        limiter.wait()
        try:
            if requester is not None:
                # Keep the injected requester interface stable for deterministic tests.
                response = requester(url, params=query, timeout=120)
            elif use_post:
                response = _thread_session().post(url, data=query, timeout=120)
            else:
                response = _thread_session().get(url, params=query, timeout=120)
            if getattr(response, "status_code", None) == 414 and not use_post and requester is None:
                # Some proxies apply a lower URL limit than the preflight bound.
                use_post = True
                response = _thread_session().post(url, data=query, timeout=120)
            response.raise_for_status()
            return response
        except requests.RequestException as exc:
            status = getattr(getattr(exc, "response", None), "status_code", None)
            if attempt < retries and (status == 429 or status is None or status >= 500):
                retry_after = None
                response_headers = getattr(getattr(exc, "response", None), "headers", {}) or {}
                try:
                    retry_after = float(response_headers.get("Retry-After"))
                except (TypeError, ValueError):
                    pass
                if status == 429:
                    # A 429 is shared rate pressure: pause every worker rather
                    # than letting concurrent retries continue at full speed.
                    limiter.cooldown(max(retry_after or 0, min(120, 10 * (2 ** attempt))))
                else:
                    time.sleep(max(retry_after or 0, min(30, 2 ** attempt)))
                continue
            detail = f"HTTP {status}" if status is not None else "network failure"
            raise RuntimeError(f"NCBI E-utilities request failed after {attempt + 1} attempt(s): {detail}; URL and credentials redacted.") from None
    raise RuntimeError("NCBI E-utilities request failed; URL and credentials redacted.")


def _bounded_results(tasks, worker, max_workers: int):
    """Yield successes while keeping a small multiple of workers in flight.

    If one worker fails, stop submitting new work but drain the remaining
    in-flight successes before raising. Callers can checkpoint those results.
    """
    iterator = iter(tasks)
    failure = None
    with ThreadPoolExecutor(max_workers=max_workers, thread_name_prefix="v3-pubmed") as pool:
        futures = set()
        for _ in range(max_workers * 2):
            try:
                futures.add(pool.submit(worker, next(iterator)))
            except StopIteration:
                break
        while futures:
            done, futures = wait(futures, return_when=FIRST_COMPLETED)
            for future in done:
                try:
                    result = future.result()
                except Exception as exc:
                    if failure is None:
                        failure = exc
                    continue
                yield result
                if failure is None:
                    try:
                        futures.add(pool.submit(worker, next(iterator)))
                    except StopIteration:
                        pass
    if failure is not None:
        raise failure


def _fetch_batches(match_frame: pd.DataFrame, article_cache: dict[str, dict[str, Any]], *, batch_size: int, max_batches: int | None) -> tuple[list[str], list[list[str]]]:
    """Select missing PMIDs fairly across domains when a bounded fetch is requested."""
    if match_frame.empty:
        return [], []
    missing = set(match_frame["pmid"].astype(str)) - set(article_cache)
    if not missing:
        return [], []
    if max_batches is None:
        selected = sorted(missing)
        return selected, [selected[i:i + batch_size] for i in range(0, len(selected), batch_size)]
    if max_batches < 0:
        raise ValueError("max_fetch_batches cannot be negative.")
    groups = sorted(match_frame["organism_group"].dropna().astype(str).unique())
    queues = {
        group: sorted(set(match_frame.loc[match_frame["organism_group"].astype(str).eq(group), "pmid"].astype(str)) & missing)
        for group in groups
    }
    offsets = {group: 0 for group in groups}
    selected: set[str] = set()
    batches: list[list[str]] = []
    while len(batches) < max_batches:
        progressed = False
        for group in groups:
            if len(batches) >= max_batches:
                break
            queue = queues[group]
            start = offsets[group]
            while start < len(queue) and queue[start] in selected:
                start += 1
            chunk = [pmid for pmid in queue[start:start + batch_size] if pmid not in selected]
            offsets[group] = min(len(queue), start + batch_size)
            if chunk:
                selected.update(chunk)
                batches.append(chunk)
                progressed = True
        if not progressed:
            break
    return sorted(selected), batches


def sample_species_pmids(
    matches: pd.DataFrame,
    *,
    max_pmids_per_species: int | None = 200,
    oldest_pmids_per_species: int = 30,
    random_seed: int = 42,
) -> pd.DataFrame:
    """Mark a reproducible per-organism PMID sample using the V2 policy.

    Retain the oldest PMID identifiers first, then draw the remaining sample
    uniformly without replacement. PMID order is the same convention used in
    V2; this is an identifier-based proxy for publication age, not a date sort.
    All input search matches remain in the returned frame for full search counts.
    """
    required = {"organism_id", "organism_group", "pmid"}
    if not required.issubset(matches.columns):
        raise ValueError(f"Sampling input is missing columns: {sorted(required - set(matches.columns))}.")
    if max_pmids_per_species is not None and int(max_pmids_per_species) < 1:
        raise ValueError("max_pmids_per_species must be positive or None.")
    if int(oldest_pmids_per_species) < 0:
        raise ValueError("oldest_pmids_per_species cannot be negative.")
    result = matches.copy().reset_index(drop=True)
    result["pmid"] = result["pmid"].astype(str)
    result["selected_for_sampling"] = False
    result["sample_selection"] = "not_sampled"
    for _, indices in result.groupby("organism_id", sort=False).groups.items():
        group = result.loc[indices]
        pmids = sorted(set(group["pmid"].astype(str)))
        if max_pmids_per_species is None or len(pmids) <= int(max_pmids_per_species):
            sampled = pmids
            oldest = pmids
            random_sample: list[str] = []
            selection_for_pmid = {pmid: "all" for pmid in sampled}
        else:
            maximum = int(max_pmids_per_species)
            oldest_count = min(int(oldest_pmids_per_species), maximum, len(pmids))
            oldest = sorted(pmids, key=lambda value: (0, int(value)) if value.isdigit() else (1, value))[:oldest_count]
            remaining = [pmid for pmid in pmids if pmid not in set(oldest)]
            budget = maximum - len(oldest)
            if len(remaining) > budget:
                random_sample = pd.Series(remaining).sample(n=budget, random_state=int(random_seed)).astype(str).tolist()
            else:
                random_sample = sorted(remaining, key=lambda value: (0, int(value)) if value.isdigit() else (1, value))
            sampled = sorted(set(oldest + random_sample), key=lambda value: (0, int(value)) if value.isdigit() else (1, value))
            selection_for_pmid = {pmid: "oldest" for pmid in oldest}
            selection_for_pmid.update({pmid: "random" for pmid in random_sample})
        sampled_set = set(sampled)
        selected_indices = group.index[group["pmid"].isin(sampled_set)]
        result.loc[selected_indices, "selected_for_sampling"] = True
        result.loc[selected_indices, "sample_selection"] = result.loc[selected_indices, "pmid"].map(selection_for_pmid)
    return result


def parse_article_xml(payload: str) -> list[dict[str, Any]]:
    root = ET.fromstring(payload)
    rows = []
    for article in root.findall(".//PubmedArticle"):
        pmid = article.findtext(".//PMID")
        medline = article.find(".//MedlineCitation/Article")
        if not pmid or medline is None:
            continue
        title_node = medline.find("ArticleTitle")
        title = "".join(title_node.itertext()) if title_node is not None else ""
        abstract = " ".join("".join(node.itertext()).strip() for node in medline.findall(".//Abstract/AbstractText"))
        year = medline.findtext(".//PubDate/Year") or ""
        pub_date = medline.findtext(".//PubDate/MedlineDate") or year
        authors = [" ".join(part for part in (node.findtext("ForeName"), node.findtext("LastName")) if part) for node in medline.findall(".//AuthorList/Author")]
        doi = next((node.text for node in article.findall(".//ELocationID") if node.attrib.get("EIdType") == "doi"), None)
        rows.append({"pmid": str(pmid), "title": title, "abstract": abstract, "publication_date": pub_date,
                     "publication_year": year, "authors": authors, "journal": medline.findtext(".//Journal/Title") or "",
                     "doi": doi, "mesh_terms": ["".join(n.itertext()) for n in article.findall(".//MeshHeading/DescriptorName")],
                     "keywords": ["".join(n.itertext()) for n in article.findall(".//Keyword")], "source": "pubmed"})
    return rows


def _expand_query_overflow(
    query_rows: list[dict[str, Any]],
    *,
    query_cache_path: Path,
    segment_cache_path: Path,
    retmax: int,
    email: str,
    api_key: str | None,
    limiter: RateLimiter,
    progress_every: int,
    requester=None,
    log: bool = True,
) -> list[dict[str, Any]]:
    """Expand over-cap PubMed searches into resumable date/UID segments.

    NCBI limits ESearch's PubMed ID list to 10,000 records. Partitioning by
    Entrez date avoids result truncation; exceptionally dense single days are
    subdivided by PubMed UID ranges.
    """
    if not query_rows:
        return query_rows
    segment_cache: dict[str, dict[str, Any]] = {}
    if segment_cache_path.exists():
        for line in segment_cache_path.read_text(encoding="utf-8").splitlines():
            try:
                item = json.loads(line)
                segment_cache[str(item["segment_id"])] = item
            except (ValueError, KeyError):
                continue
    segment_cache_path.parent.mkdir(parents=True, exist_ok=True)
    started_all = time.monotonic()
    segment_requests_all = 0
    start_date = date(1800, 1, 1)
    end_date = date.today()
    max_uid = 1_000_000_000

    for query_index, item in enumerate(query_rows, 1):
        original_pmids = {str(pmid) for pmid in item.get("pmids", [])}
        base_query = str(item["query"])
        query_id = str(item["query_id"])
        stack: list[tuple[date, date, int | None, int | None]] = [(start_date, end_date, None, None)]
        expanded_pmids: set[str] = set()
        expanded_count = 0
        segment_count = 0
        while stack:
            mindate, maxdate, uid_low, uid_high = stack.pop()
            segment_key_value = {
                "query_id": query_id,
                "mindate": mindate.isoformat(),
                "maxdate": maxdate.isoformat(),
                "uid_low": uid_low,
                "uid_high": uid_high,
                "retmax": retmax,
            }
            segment_id = hashlib.sha256(json.dumps(segment_key_value, sort_keys=True).encode()).hexdigest()
            cached = segment_cache.get(segment_id)
            if cached is None:
                term = base_query
                if uid_low is not None and uid_high is not None:
                    term = f"({base_query}) AND {uid_low}:{uid_high}[UID]"
                response = _request(
                    f"{EUTILS}/esearch.fcgi",
                    {
                        "db": "pubmed", "term": term, "retmax": retmax, "retmode": "json",
                        "mindate": mindate.strftime("%Y/%m/%d"),
                        "maxdate": maxdate.strftime("%Y/%m/%d"),
                        "datetype": "edat",
                    },
                    email=email, api_key=api_key, limiter=limiter, requester=requester,
                )
                result = response.json().get("esearchresult", {})
                cached = {
                    **segment_key_value,
                    "segment_id": segment_id,
                    "count": int(result.get("count", 0)),
                    "pmids": [str(pmid) for pmid in result.get("idlist", [])],
                    "completed_at": time.time(),
                }
                append_jsonl(segment_cache_path, cached)
                segment_cache[segment_id] = cached
                segment_requests_all += 1
            segment_count += 1
            count = int(cached.get("count", 0))
            segment_pmids = {str(pmid) for pmid in cached.get("pmids", [])}
            if count == 0 or (count <= retmax and len(segment_pmids) >= count):
                expanded_pmids.update(segment_pmids)
                expanded_count += count
            elif uid_low is not None and uid_high is not None:
                if uid_high <= uid_low:
                    raise RuntimeError(f"PubMed overflow segment {segment_id} remains truncated at a single UID.")
                midpoint = (uid_low + uid_high) // 2
                stack.append((mindate, maxdate, midpoint + 1, uid_high))
                stack.append((mindate, maxdate, uid_low, midpoint))
            elif mindate < maxdate:
                midpoint = mindate + (maxdate - mindate) // 2
                stack.append((midpoint + timedelta(days=1), maxdate, None, None))
                stack.append((mindate, midpoint, None, None))
            else:
                stack.append((mindate, maxdate, 1, max_uid))

            if log and segment_requests_all and segment_requests_all % max(1, progress_every // 10) == 0:
                elapsed = max(0.01, time.monotonic() - started_all)
                rate = segment_requests_all / elapsed
                queued_floor = len(stack)
                print(
                    f"PubMed overflow: query {query_index}/{len(query_rows)}; "
                    f"{segment_count} segments examined; {queued_floor} queued; "
                    f"ETA at least {queued_floor / rate / 60:.1f} min."
                )

        merged_pmids = sorted(original_pmids | expanded_pmids)
        verified = expanded_count == len(expanded_pmids) and original_pmids.issubset(expanded_pmids)
        updated = {
            **item,
            "pmids": merged_pmids,
            "count": len(merged_pmids),
            "search_complete": verified,
            "search_incomplete_count": 0 if verified else max(1, len(merged_pmids) - expanded_count),
            "overflow_expanded": True,
            "overflow_date_field": "edat",
            "overflow_segments": segment_count,
            "overflow_completed_at": time.time(),
        }
        append_jsonl(query_cache_path, updated)
        item.clear()
        item.update(updated)
        if log:
            elapsed = max(0.01, time.monotonic() - started_all)
            rate = query_index / elapsed
            eta = (len(query_rows) - query_index) / rate / 60 if rate else 0
            status = "complete" if verified else "verification mismatch"
            print(
                f"PubMed overflow: {query_index}/{len(query_rows)} species queries; "
                f"{len(merged_pmids):,} PMIDs; {segment_count} segments; {status}; ETA {eta:.1f} min."
            )
    return query_rows


def acquire_pubmed(
    species: pd.DataFrame,
    aliases: pd.DataFrame,
    *,
    email: str | None = None,
    api_key: str | None = None,
    cache_dir: str | Path = "data/02_pubmed/cache",
    release_dir: str | Path = "data/02_pubmed/release",
    max_workers: int = 10,
    request_delay_seconds: float | None = None,
    max_calls: int | None = None,
    max_fetch_batches: int | None = None,
    retmax: int = 10000,
    batch_size: int = 200,
    max_pmids_per_species: int | None = 200,
    oldest_pmids_per_species: int = 30,
    pmid_sample_random_seed: int = 42,
    progress_every: int = 25,
    expand_overflow: bool = False,
    requester=None,
    log: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Search each organism, fetch each PMID once, and checkpoint raw JSONL caches."""
    email = resolve_ncbi_email(email)
    api_key = (api_key or os.environ.get("NCBI_API_KEY", "") or os.environ.get("EUTILS_API_KEY", "")).strip() or _ncbi_key_from_keychain()
    if not email:
        raise ValueError("Set NCBI_EMAIL before live PubMed requests.")
    if max_workers < 1 or batch_size < 1 or progress_every < 1:
        raise ValueError("max_workers, batch_size, and progress_every must be positive.")
    # Keep a margin under NCBI's keyed ceiling; simultaneous jobs may share a
    # key or network address, so running exactly at the nominal limit is brittle.
    request_delay_seconds = (0.20 if api_key else 0.34) if request_delay_seconds is None else float(request_delay_seconds)
    retmax = min(max(1, int(retmax)), 10000)
    cache, release = Path(cache_dir), Path(release_dir)
    cache.mkdir(parents=True, exist_ok=True); release.mkdir(parents=True, exist_ok=True)
    query_path, article_path = cache / "pubmed_search_results.jsonl", cache / "pubmed_article_records.jsonl"
    query_cache: dict[str, dict[str, Any]] = {}
    if query_path.exists():
        for line in query_path.read_text(encoding="utf-8").splitlines():
            try:
                item = json.loads(line); query_cache[str(item["query_id"])] = item
            except (ValueError, KeyError):
                continue
    article_cache: dict[str, dict[str, Any]] = {}
    if article_path.exists():
        for line in article_path.read_text(encoding="utf-8").splitlines():
            try:
                item = json.loads(line); article_cache[str(item["pmid"])] = item
            except (ValueError, KeyError):
                continue
    aliases_by_id = aliases.groupby("organism_id")["name"].apply(lambda x: list(dict.fromkeys(map(str, x.dropna())))).to_dict() if not aliases.empty else {}
    tasks = []
    active_query_ids: set[str] = set()
    for row in species[["organism_id", "organism_group", "accepted_name"]].drop_duplicates("organism_id").itertuples(index=False):
        raw_names = list(aliases_by_id.get(row.organism_id, []))
        names = [name for name in raw_names if len(" ".join(str(name).split())) >= 3 and not re.fullmatch(r"\d+(?:\.\d+)?", " ".join(str(name).split()))]
        excluded_alias_count = len(raw_names) - len(names)
        if row.accepted_name not in names:
            names.insert(0, row.accepted_name)
        query = build_name_query(names)
        digest = hashlib.sha256(f"{query}\nretmax={retmax}".encode()).hexdigest()
        query_id = f"q_{digest[:20]}"
        active_query_ids.add(query_id)
        if query_id not in query_cache:
            reusable = next((item for item in query_cache.values() if item.get("query") == query and (item.get("search_complete", False) or len(item.get("pmids", [])) >= retmax)), None)
            if reusable is not None:
                pmids = list(reusable.get("pmids", []))[:retmax]
                count = int(reusable.get("count", len(pmids)))
                query_cache[query_id] = {**reusable, "query_id": query_id, "query_hash": digest,
                    "organism_id": row.organism_id, "organism_group": row.organism_group,
                    "accepted_name": row.accepted_name, "matched_names": names, "pmids": pmids,
                    "search_complete": len(pmids) >= count, "search_incomplete_count": max(0, count - len(pmids)),
                    "excluded_alias_count": excluded_alias_count}
                append_jsonl(query_path, query_cache[query_id])
            else:
                tasks.append((row, names, query, digest, query_id, excluded_alias_count))
    if max_calls is not None:
        tasks = tasks[:int(max_calls)]
    if log:
        print(f"PubMed search: {len(species):,} species rows; {len(query_cache):,} cached queries; {len(tasks):,} pending.")
    limiter = RateLimiter(request_delay_seconds)
    def search(task: tuple[Any, ...]) -> dict[str, Any]:
        row, names, query, digest, query_id, excluded_alias_count = task
        response = _request(f"{EUTILS}/esearch.fcgi", {"db": "pubmed", "term": query, "retmax": retmax, "retmode": "json"}, email=email, api_key=api_key, limiter=limiter, requester=requester)
        result = response.json().get("esearchresult", {})
        pmids = result.get("idlist", [])
        count = int(result.get("count", 0))
        return {"query_id": query_id, "query_hash": digest, "organism_id": row.organism_id, "organism_group": row.organism_group,
                "accepted_name": row.accepted_name, "matched_names": names, "query": query,
                "pmids": pmids, "count": count, "search_complete": len(pmids) >= count,
                "search_incomplete_count": max(0, count - len(pmids)), "excluded_alias_count": excluded_alias_count,
                "completed_at": time.time()}
    search_start = time.monotonic()
    for index, item in enumerate(_bounded_results(tasks, search, max_workers), 1):
        query_cache[item["query_id"]] = item; append_jsonl(query_path, item)
        if log and (index % progress_every == 0 or index == len(tasks)):
            rate = index / max(0.01, time.monotonic() - search_start)
            print(f"PubMed search: {index:,}/{len(tasks):,}; ETA {(len(tasks)-index)/rate/60:.1f} min.")
    active_queries = [query_cache[query_id] for query_id in active_query_ids if query_id in query_cache]
    if expand_overflow:
        overflow_queries = [item for item in active_queries if not item.get("search_complete", True)]
        if log and overflow_queries:
            print(f"PubMed overflow: expanding {len(overflow_queries):,} species queries beyond the {retmax:,}-PMID ESearch cap.")
        _expand_query_overflow(
            overflow_queries,
            query_cache_path=query_path,
            segment_cache_path=cache / "pubmed_search_overflow_segments.jsonl",
            retmax=retmax,
            email=email,
            api_key=api_key,
            limiter=limiter,
            progress_every=progress_every,
            requester=requester,
            log=log,
        )
    associations = []
    for item in active_queries:
        for pmid in item.get("pmids", []):
            associations.append({"organism_id": item["organism_id"], "organism_group": item["organism_group"], "accepted_name": item["accepted_name"], "pmid": str(pmid),
                                 "matched_names": item["matched_names"], "query_id": item["query_id"], "query_hash": item["query_hash"], "selected_for_fetch": True})
    match_frame = pd.DataFrame(associations, columns=["organism_id", "organism_group", "accepted_name", "pmid", "matched_names", "query_id", "query_hash", "selected_for_fetch"])
    if not match_frame.empty:
        match_frame = match_frame.drop_duplicates(["organism_id", "pmid"])
    match_frame = sample_species_pmids(
        match_frame,
        max_pmids_per_species=max_pmids_per_species,
        oldest_pmids_per_species=oldest_pmids_per_species,
        random_seed=pmid_sample_random_seed,
    )
    match_frame["selected_for_fetch"] = False
    sample_matches = match_frame.loc[match_frame["selected_for_sampling"]].copy()
    missing_pmids, batches = _fetch_batches(sample_matches, article_cache, batch_size=batch_size, max_batches=int(max_fetch_batches) if max_fetch_batches is not None else None)
    def fetch(batch: list[str]) -> list[dict[str, Any]]:
        parse_error = None
        for attempt in range(3):
            response = _request(
                f"{EUTILS}/efetch.fcgi",
                {"db": "pubmed", "id": ",".join(batch), "retmode": "xml"},
                email=email, api_key=api_key, limiter=limiter, requester=requester,
            )
            try:
                return parse_article_xml(response.text)
            except ET.ParseError as exc:
                parse_error = exc
                if attempt < 2:
                    time.sleep(2 ** attempt)
        if len(batch) > 1:
            midpoint = len(batch) // 2
            if log:
                print(f"PubMed fetch: malformed XML for a {len(batch)}-PMID batch; retrying as {midpoint} and {len(batch)-midpoint} PMID sub-batches.")
            return fetch(batch[:midpoint]) + fetch(batch[midpoint:])
        raise RuntimeError("PubMed returned malformed XML for a single-PMID fetch after retries.") from parse_error
    if log:
        sampled_rows = sample_matches[["organism_id", "pmid"]].drop_duplicates()
        sampled_count = len(sampled_rows)
        species_count = sample_matches["organism_id"].nunique()
        cap_text = "all PMIDs" if max_pmids_per_species is None else f"up to {int(max_pmids_per_species)} per species (oldest {min(int(oldest_pmids_per_species), int(max_pmids_per_species))}, then seed {int(pmid_sample_random_seed)} sample)"
        print(f"PubMed abstract sample: {sampled_count:,} species–PMID rows across {species_count:,} species; {cap_text}.")
        print(f"PubMed abstract fetch: {len(article_cache):,} cached; {len(missing_pmids):,} sampled PMIDs missing across {len(batches):,} batches.")
    fetch_start = time.monotonic()
    for index, batch_rows in enumerate(_bounded_results(batches, fetch, max_workers), 1):
        for item in batch_rows:
            article_cache[str(item["pmid"])] = item; append_jsonl(article_path, item)
        if log and (index % progress_every == 0 or index == len(batches)):
            rate = index / max(0.01, time.monotonic() - fetch_start)
            print(f"PubMed fetch: {index:,}/{len(batches):,} batches; ETA {(len(batches)-index)/rate/60:.1f} min.")
    sampled_pmids = set(sample_matches["pmid"].astype(str))
    abstracts_frame = pd.DataFrame([record for pmid, record in article_cache.items() if str(pmid) in sampled_pmids])
    if not match_frame.empty:
        if "pmid" in abstracts_frame and "abstract" in abstracts_frame:
            match_frame = match_frame.merge(abstracts_frame[["pmid", "abstract"]], on="pmid", how="left")
        else:
            match_frame["abstract"] = ""
        match_frame["selected_for_fetch"] = (
            match_frame["selected_for_sampling"]
            & match_frame["abstract"].fillna("").astype(str).str.strip().ne("")
        )
        match_frame = match_frame.drop(columns="abstract")
    atomic_parquet(match_frame, release / "species_abstracts.parquet")
    atomic_parquet(abstracts_frame, release / "abstracts.parquet")
    active_group_summaries = {}
    for group in sorted({str(item["organism_group"]) for item in active_queries}):
        group_queries = [item for item in active_queries if str(item["organism_group"]) == group]
        active_group_summaries[group] = {
            "query_count": len(group_queries), "total_pubmed_search_count": sum(int(item.get("count", 0)) for item in group_queries),
            "pmids_returned": sum(len(item.get("pmids", [])) for item in group_queries),
            "incomplete_queries": sum(not item.get("search_complete", True) for item in group_queries),
            "pmids_beyond_retmax": sum(int(item.get("search_incomplete_count", 0)) for item in group_queries),
            "excluded_short_or_numeric_aliases": sum(int(item.get("excluded_alias_count", 0)) for item in group_queries),
        }
    metadata = {
        "query_count": len(active_queries),
        "incomplete_query_count": sum(not item.get("search_complete", True) for item in active_queries),
        "search_incomplete_pmid_count": sum(int(item.get("search_incomplete_count", 0)) for item in active_queries),
        "overflow_expanded_query_count": sum(bool(item.get("overflow_expanded")) for item in active_queries),
        "query_summary_by_group": active_group_summaries,
        "species_abstract_matches": len(match_frame),
        "sampled_species_pmid_matches": int(match_frame["selected_for_sampling"].sum()) if not match_frame.empty else 0,
        "sample_cap_per_species": max_pmids_per_species,
        "oldest_pmids_per_species": int(oldest_pmids_per_species),
        "pmid_sample_random_seed": int(pmid_sample_random_seed),
        "abstract_count": len(abstracts_frame),
        "cached_at": pd.Timestamp.utcnow().isoformat(),
    }
    (release / "acquisition_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    return match_frame, abstracts_frame
