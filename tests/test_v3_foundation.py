from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest
import requests

from pathogen_atlas_v3.decisions.features import load_feature_catalog, questions_for_domain
from pathogen_atlas_v3.decisions.jev import OpenRouterJevClient, extract_noul_probability
from pathogen_atlas_v3.decisions.runner import prepare_species_abstracts, run_decision_stage, run_domain_features, run_pathogen_gate
from pathogen_atlas_v3.abstracts.pubmed import RateLimiter, _fetch_batches, _request, acquire_pubmed, build_name_query
from pathogen_atlas_v3.pipelines.phase2_pilot import select_pilot_species
from pathogen_atlas_v3.tables.interest_index import build_pathogenic_interest_index
from pathogen_atlas_v3.taxonomy.adapters import load_ictv_msl_export, load_ictv_vmr_export, load_species_fungorum_export


class FakeClient:
    model = "typesafe/jev-test"
    provider = "fake"

    def __init__(self, *, malformed=False, fail=False):
        self.calls = []
        self.malformed = malformed
        self.fail = fail
        self.last_attempts = 1

    def decide_many(self, state, *, questions):
        self.calls.append((state, questions))
        if self.fail:
            raise RuntimeError("temporary failure")
        answers = {}
        for qid in questions:
            if not self.malformed:
                answers[qid] = {"type": "noul", "noul": 0.91}
        return {"answers": answers}


class FakeHTTPResponse:
    def __init__(self, payload, status=200):
        self.payload, self.status_code = payload, status

    def raise_for_status(self):
        if self.status_code >= 400:
            error = requests.HTTPError(f"HTTP {self.status_code}")
            error.response = self
            raise error

    def json(self):
        return self.payload


class FakeSession:
    def __init__(self, responses):
        self.responses = list(responses)
        self.calls = []

    def post(self, url, **kwargs):
        self.calls.append((url, kwargs))
        response = self.responses.pop(0)
        if isinstance(response, Exception):
            raise response
        return response


class FakeGetResponse:
    def __init__(self, *, json_value=None, text="", status=200):
        self._json, self.text, self.status_code = json_value, text, status
    def raise_for_status(self):
        if self.status_code >= 400:
            error = requests.HTTPError(f"HTTP {self.status_code}")
            error.response = self
            raise error
    def json(self): return self._json


def tables(abstract="Disease in a patient."):
    abstracts = pd.DataFrame([
        {"pmid": "1", "title": "Human infection", "abstract": abstract},
        {"pmid": "2", "title": "Not selected", "abstract": "Another abstract"},
    ])
    matches = pd.DataFrame([
        {"organism_id": "fungi:1", "organism_group": "fungi", "accepted_name": "Candida test", "pmid": "1", "matched_names": ["Candida test", "Candida old"], "selected_for_fetch": True},
        {"organism_id": "fungi:1", "organism_group": "fungi", "accepted_name": "Candida test", "pmid": "1", "matched_names": ["Candida test"], "selected_for_fetch": True},
        {"organism_id": "fungi:1", "organism_group": "fungi", "accepted_name": "Candida test", "pmid": "2", "matched_names": ["Candida test"], "selected_for_fetch": False},
    ])
    return abstracts, matches


def test_prepare_selected_rows_and_duplicate_species_pmid(tmp_path):
    abstracts, matches = tables()
    result = prepare_species_abstracts(abstracts, matches)
    assert len(result) == 1
    assert result.iloc[0].matched_names == ["Candida test", "Candida old"]
    client = FakeClient()
    run_decision_stage(abstracts, matches, kind="feature", questions_by_domain={"fungi": {"a": {"type": "noul", "instructions": "q", "criteria": {}}}}, run_id="state", question_version="q1", profile_version="p1", ledger_path=tmp_path/"state.parquet", raw_response_path=tmp_path/"state.jsonl", max_calls=1, client=client, log=False)
    assert client.calls[0][0]["organism"]["organism_group"] == "fungi"


def test_runner_resume_changed_input_question_version_and_retry(tmp_path):
    abstracts, matches = tables()
    questions = {"feature_a": {"type": "noul", "definition": "Question?", "category": "shared", "instructions": "Answer from text.", "criteria": {"true": "yes", "false": "no"}}}
    client = FakeClient()
    kwargs = dict(kind="feature", questions_by_domain={"fungi": questions}, run_id="test-run", question_version="q-v1", profile_version="p-v1", ledger_path=tmp_path/"ledger.parquet", raw_response_path=tmp_path/"raw.jsonl", client=client, log=False)
    first = run_decision_stage(abstracts, matches, **kwargs)
    assert len(client.calls) == 1
    assert first.iloc[-1].status == "completed"
    run_decision_stage(abstracts, matches, **kwargs)
    assert len(client.calls) == 1
    changed_abstracts, _ = tables("Disease in a patient, revised.")
    run_decision_stage(changed_abstracts, matches, **kwargs)
    assert len(client.calls) == 2
    run_decision_stage(abstracts, matches, **{**kwargs, "question_version": "q-v2"})
    assert len(client.calls) == 3

    failed = FakeClient(malformed=True)
    retry_kwargs = {**kwargs, "ledger_path": tmp_path/"retry.parquet", "raw_response_path": tmp_path/"retry.jsonl", "client": failed}
    failed_ledger = run_decision_stage(abstracts, matches, **retry_kwargs)
    assert failed_ledger.iloc[-1].status == "failed"
    recovered = FakeClient()
    recovered_ledger = run_decision_stage(abstracts, matches, **{**retry_kwargs, "client": recovered})
    assert len(recovered.calls) == 1
    assert recovered_ledger.iloc[-1].status == "completed"


def test_bounded_fake_end_to_end_all_organism_groups(tmp_path):
    groups = [("bacteria", "Bacillus test"), ("fungi", "Candida test"), ("viruses", "Example virus")]
    abstracts = pd.DataFrame([{"pmid": str(i), "title": f"Human disease {i}", "abstract": "A human case with infection."} for i in range(1, 4)])
    matches = pd.DataFrame([{"organism_id": f"{group}:{i}", "organism_group": group, "accepted_name": name, "pmid": str(i), "matched_names": [name], "selected_for_fetch": True} for i, (group, name) in enumerate(groups, 1)])
    prompts = {"domains": {group: {"instructions": f"Assess human disease for {group}.", "criteria": {"true": "direct evidence", "false": "no evidence"}} for group, _ in groups}}
    feature = {"one_feature": {"type": "noul", "definition": "Feature?", "category": "test", "instructions": "Assess this feature.", "criteria": {"true": "yes", "false": "no"}}}
    client = FakeClient()
    common = {"run_id": "gate-run", "question_version": "gate-v1", "profile_version": "gate-profile-v1", "ledger_path": tmp_path/"ledger.parquet", "raw_response_path": tmp_path/"raw.jsonl", "max_calls": 3, "client": client, "log": False}
    gate = run_pathogen_gate(abstracts, matches, prompts, **common)
    assert gate.loc[gate.decision_kind.eq("pathogen_evidence"), "status"].eq("completed").sum() == 3
    features = run_domain_features(abstracts, matches, gate, threshold=.75, questions_by_domain={group: feature for group, _ in groups}, gate_run_id="gate-run", run_id="features-run", question_version="features-v1", profile_version="profiles-v1", ledger_path=tmp_path/"ledger.parquet", raw_response_path=tmp_path/"raw.jsonl", max_calls=3, client=client, log=False)
    selected = features.loc[features.decision_kind.eq("feature")]
    assert set(selected.organism_group) == {"bacteria", "fungi", "viruses"}
    assert len(selected) == 3


def test_phase2_pilot_selection_is_bounded_deterministic_and_multidomain():
    species = pd.DataFrame([
        {"organism_id": f"{group}:{i}", "organism_group": group,
         "accepted_name": f"{group} candidate {i}", "taxonomy_path": f"phylum=Phylum{i%3};family=Family{i%4}",
         "mapping_status": "accepted"}
        for group in ("bacteria", "fungi", "viruses") for i in range(15)
    ])
    first, manifest = select_pilot_species(species, pd.DataFrame(), per_group=5, seed=12)
    second, _ = select_pilot_species(species, pd.DataFrame(), per_group=5, seed=12)
    assert first.organism_id.tolist() == second.organism_id.tolist()
    assert first.groupby("organism_group").size().to_dict() == {"bacteria": 5, "fungi": 5, "viruses": 5}
    assert manifest["selected_per_group"] == {"bacteria": 5, "fungi": 5, "viruses": 5}


def test_bounded_abstract_fetch_allocates_batches_across_groups():
    matches = pd.DataFrame([
        {"organism_group": group, "pmid": f"{group}-{i}"}
        for group in ("bacteria", "fungi", "viruses") for i in range(300)
    ])
    pmids, batches = _fetch_batches(matches, {}, batch_size=100, max_batches=3)
    assert len(batches) == 3 and len(pmids) == 300
    assert [batch[0].split("-")[0] for batch in batches] == ["bacteria", "fungi", "viruses"]


def test_typed_payload_probability_validation_and_transient_retry():
    session = FakeSession([requests.ConnectionError("offline"), FakeHTTPResponse({"answers": {"q": {"type": "noul", "noul": .8}}})])
    client = OpenRouterJevClient(api_key="test", retries=1, backoff_seconds=0, session=session)
    result = client.decide_many({"organism": {"accepted_name": "X"}}, questions={"q": {"type": "noul", "instructions": "Question", "criteria": {"true": "yes", "false": "no"}}})
    assert client.last_attempts == 2
    assert session.calls[0][1]["json"]["model"] == "typesafe/jev-1.14"
    assert session.calls[0][1]["json"]["questions"]["q"]["type"] == "noul"
    assert extract_noul_probability(result, "q") == .8
    with pytest.raises(ValueError):
        extract_noul_probability({"answers": {"q": {"type": "text", "text": "yes"}}}, "q")


def test_domain_profiles_are_ordered_and_exclude_inapplicable():
    root = Path(__file__).parents[1]
    catalog, profiles = load_feature_catalog(root/"configs/features.yaml", root/"configs/feature_profiles.yaml", allow_provisional=True)
    assert len(catalog["features"]) == 88
    assert "anaerobic" in questions_for_domain(catalog, profiles, "bacteria", allow_provisional=True)
    assert "anaerobic" not in questions_for_domain(catalog, profiles, "fungi", allow_provisional=True)
    assert "antifungal_resistance" in questions_for_domain(catalog, profiles, "fungi", allow_provisional=True)
    assert "viral_latency" in questions_for_domain(catalog, profiles, "viruses", allow_provisional=True)
    assert "zoonosis" in questions_for_domain(catalog, profiles, "fungi", allow_provisional=True)


def test_taxonomy_adapters_preserve_domain_source_ids_and_aliases(tmp_path):
    fungi = tmp_path/"fungi.csv"
    pd.DataFrame([{"current_name": "Candida test", "species_fungorum_id": "123", "synonyms": "Old Candida;Candida prior", "genus": "Candida"}]).to_csv(fungi, index=False)
    species, aliases = load_species_fungorum_export(fungi, release="SF-2026-01")
    assert species.iloc[0].organism_group == "fungi"
    assert species.iloc[0].taxonomy_release == "SF-2026-01"
    assert set(aliases.name) == {"Candida test", "Old Candida", "Candida prior"}
    viruses = tmp_path/"ictv.csv"
    pd.DataFrame([{"Species": "Example virus species", "Taxon ID": "ICTV42", "Virus names": "Example isolate"}]).to_csv(viruses, index=False)
    virus_species, virus_aliases = load_ictv_msl_export(viruses, release="MSL40")
    assert virus_species.iloc[0].source_taxon_id == "ICTV42"
    assert "virus_or_strain_alias" in set(virus_aliases.name_type)


def test_pinned_source_layouts_map_fungal_synonyms_and_ictv_hierarchy(tmp_path):
    fungi = tmp_path / "fungi.csv"
    pd.DataFrame([
        {"Index Fungorum ID": "100", "Name": "Candida accepted", "Rank": "Species", "Taxonomic Status": "Accepted", "Accepted Name Index Fungorum ID": "100", "Accepted Name": "Candida accepted", "Genus": "Candida"},
        {"Index Fungorum ID": "101", "Name": "Candida synonym", "Rank": "Species", "Taxonomic Status": "Synonym", "Accepted Name Index Fungorum ID": "100", "Accepted Name": "Candida accepted", "Genus": "Candida"},
        {"Index Fungorum ID": "102", "Name": "Unresolved fungus", "Rank": "Species", "Taxonomic Status": "Synonym", "Accepted Name Index Fungorum ID": "999", "Accepted Name": "Missing target", "Genus": "Unresolved"},
    ]).to_csv(fungi, index=False)
    fungus_species, fungus_aliases = load_species_fungorum_export(fungi, release="SF-snapshot")
    assert len(fungus_species) == 1
    assert fungus_species.iloc[0].source_taxon_id == "100"
    assert set(fungus_aliases.loc[fungus_aliases.organism_id.notna(), "name"]) == {"Candida accepted", "Candida synonym"}
    unresolved = fungus_aliases.loc[fungus_aliases.mapping_status.eq("unresolved_accepted_name_target")]
    assert len(unresolved) == 1 and pd.isna(unresolved.iloc[0].organism_id)

    virus_book = tmp_path / "ictv.xlsx"
    with pd.ExcelWriter(virus_book) as writer:
        pd.DataFrame({"File version": ["41"]}).to_excel(writer, sheet_name="Version", index=False)
        pd.DataFrame([{"Realm": "Riboviria", "Kingdom": "Orthornavirae", "Phylum": "Pisuviricota", "Class": "Pisoniviricetes", "Order": "Nidovirales", "Family": "Coronaviridae", "Genus": "Betacoronavirus", "Species": "Examplevirus humanus", "ICTV_ID": "ICTV123", "Genome": "ssRNA(+)"}]).to_excel(writer, sheet_name="MSL", index=False)
    virus_species, _ = load_ictv_msl_export(virus_book, release="MSL41-v1")
    assert virus_species.iloc[0].source_taxon_id == "ICTV123"
    assert "family=Coronaviridae" in virus_species.iloc[0].taxonomy_path
    assert '"Genome": "ssRNA(+)"' in virus_species.iloc[0].source_record_json

    vmr_book = tmp_path / "vmr.xlsx"
    with pd.ExcelWriter(vmr_book) as writer:
        pd.DataFrame([{"ICTV_ID": "ICTV123", "Isolate ID": "VMR1", "Species": "Examplevirus humanus", "Virus name(s)": "Example virus", "Virus isolate designation": "H1", "Virus GENBANK accession": "AB123", "Genome": "ssRNA(+)"}]).to_excel(writer, sheet_name="VMR MSL41", index=False)
    variants, variant_aliases = load_ictv_vmr_export(vmr_book, virus_species, release="VMR41")
    assert variants.iloc[0].organism_id == virus_species.iloc[0].organism_id
    assert variants.iloc[0].genbank_accession == "AB123"
    assert set(variant_aliases.name) == {"Example virus", "H1"}


def test_pubmed_uses_fake_requester_and_resumes_query_and_abstract_cache(tmp_path, monkeypatch):
    assert '"Candida test"[Title/Abstract]' in build_name_query(["Candida test", "Old Candida"])
    assert '"Old Candida"[Title/Abstract]' in build_name_query(["Candida test", "Old Candida"])
    assert '"17"[Title/Abstract]' not in build_name_query(["Candida test", "17"])
    calls = []
    monkeypatch.setenv("NCBI_API_KEY", "test-key")
    xml = "<PubmedArticleSet><PubmedArticle><MedlineCitation><PMID>42</PMID><Article><ArticleTitle>Example article</ArticleTitle><Abstract><AbstractText>Human disease abstract.</AbstractText></Abstract><Journal><Title>Journal</Title></Journal><AuthorList><Author><ForeName>A</ForeName><LastName>Person</LastName></Author></AuthorList></Article></MedlineCitation></PubmedArticle></PubmedArticleSet>"
    def requester(url, *, params, timeout):
        calls.append((url, params))
        if url.endswith("esearch.fcgi"):
            return FakeGetResponse(json_value={"esearchresult": {"count": "10001", "idlist": ["42"]}})
        return FakeGetResponse(text=xml)
    species = pd.DataFrame([{"organism_id": "fungi:1", "organism_group": "fungi", "accepted_name": "Candida test"}])
    aliases = pd.DataFrame([{"organism_id": "fungi:1", "name": "Candida test"}, {"organism_id": "fungi:1", "name": "Old Candida"}])
    matches, abstracts = acquire_pubmed(species, aliases, email="test@example.org", cache_dir=tmp_path/"cache", release_dir=tmp_path/"release", max_workers=1, request_delay_seconds=0, requester=requester, log=False)
    assert len(matches) == 1 and matches.iloc[0].selected_for_fetch
    assert calls[0][1]["api_key"] == "test-key"
    assert pd.read_json(tmp_path/"release"/"acquisition_metadata.json", typ="series")["incomplete_query_count"] == 1
    assert abstracts.iloc[0].title == "Example article"
    assert len(calls) == 2
    acquire_pubmed(species, aliases, email="test@example.org", cache_dir=tmp_path/"cache", release_dir=tmp_path/"release", max_workers=1, request_delay_seconds=0, retmax=500, requester=requester, log=False)
    assert len(calls) == 3  # a different result cap gets a separate query cache identity


def test_pubmed_retries_rate_limits_and_redacts_request_credentials(monkeypatch):
    monkeypatch.setattr("pathogen_atlas_v3.abstracts.pubmed.time.sleep", lambda _seconds: None)
    replies = [FakeGetResponse(json_value={}, status=429) for _ in range(4)]
    def requester(_url, *, params, timeout):
        assert params["api_key"] == "secret-key-for-test"
        return replies.pop(0)
    with pytest.raises(RuntimeError) as error:
        _request("https://example.test/query", {}, email="test@example.org", api_key="secret-key-for-test", limiter=RateLimiter(0), requester=requester)
    assert "429" in str(error.value)
    assert "secret-key-for-test" not in str(error.value)


def test_interest_index_retains_counts_and_normalizes_per_domain():
    frame = pd.DataFrame([
        {"organism_id": "b1", "organism_group": "bacteria", "n_pathogen_pmids": 2, "n_abstracts_fetched": 10, "n_pubmed_search_records": 100},
        {"organism_id": "b2", "organism_group": "bacteria", "n_pathogen_pmids": 1, "n_abstracts_fetched": 10, "n_pubmed_search_records": 10},
        {"organism_id": "f1", "organism_group": "fungi", "n_pathogen_pmids": 2, "n_abstracts_fetched": 10, "n_pubmed_search_records": 100},
    ])
    result = build_pathogenic_interest_index(frame).set_index("organism_id")
    assert result.loc["b1", "pathogenic_interest_index"] == 1.0
    assert result.loc["f1", "pathogenic_interest_index"] == 1.0
    assert result.loc["b1", "n_pubmed_search_records"] == 100
