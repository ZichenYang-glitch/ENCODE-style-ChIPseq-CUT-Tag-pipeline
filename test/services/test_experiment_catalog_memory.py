"""In-memory parity for Slice A catalog maintenance and combination queries."""

from __future__ import annotations

from datetime import datetime, timedelta, timezone
from itertools import count

import pytest

from encode_pipeline.platform.experiment_catalog import (
    CatalogBlockRequiresCatalogCommand,
    CatalogError,
)
from encode_pipeline.services.data_registry import DataRegistryService
from encode_pipeline.services.data_registry_repositories import (
    InMemoryDataRegistryRepository,
)
from encode_pipeline.services.experiment_catalog import ExperimentCatalogService
from encode_pipeline.services.experiment_catalog_repositories import (
    InMemoryExperimentCatalogRepository,
)


def _ids(prefix: str, start: int = 1):
    serial = count(start)
    return lambda: f"{prefix}{next(serial):032x}"


class _Clock:
    def __init__(self) -> None:
        self.now = datetime(2026, 11, 18, tzinfo=timezone.utc)

    def __call__(self) -> datetime:
        current = self.now
        self.now += timedelta(seconds=1)
        return current


def _stack():
    clock = _Clock()
    registry = InMemoryDataRegistryRepository(legacy_created_at=clock.now)
    data = DataRegistryService(
        repository=registry,
        project_id_factory=_ids("prj_"),
        sample_id_factory=_ids("smp_"),
        sample_revision_id_factory=_ids("smpr_"),
        now_factory=clock,
    )
    catalog = ExperimentCatalogService(
        InMemoryExperimentCatalogRepository(registry),
        now_factory=clock,
        term_id_factory=_ids("term_"),
        term_revision_id_factory=_ids("trmv_"),
        experiment_id_factory=_ids("exp_"),
        experiment_revision_id_factory=_ids("expr_"),
        sample_revision_id_factory=_ids("smpr_", start=1000),
    )
    return data, catalog


def _write(operation: str, body: dict):
    key = "".join(character if character.isalnum() else "-" for character in operation)
    return {
        "actor_user_id": "usr_" + "a" * 32,
        "idempotency_key": key[:120],
        "raw_body": __import__("json").dumps(body, sort_keys=True).encode(),
    }


def _term(catalog: ExperimentCatalogService, term_type: str, name: str) -> dict:
    body = {
        "term_type": term_type,
        "display_name": name,
        "namespace": None,
        "accession": None,
        "aliases": [],
    }
    return catalog.create_term(**_write(f"create-{term_type}-{name}", body), **body)


def _pin(term: dict) -> dict:
    return {
        "term_id": term["term_id"],
        "term_revision_id": term["revision"]["term_revision_id"],
    }


def _known(term: dict) -> dict:
    return {"status": "known", **_pin(term)}


def _multi(term: dict) -> dict:
    return {"status": "known", "values": [_pin(term)]}


def _single(status: str) -> dict:
    return {"status": status}


def _multi_status(status: str) -> dict:
    return {"status": status, "values": []}


def _catalog(**overrides: dict) -> dict:
    block = {
        "schema_version": 1,
        "species": _single("unknown"),
        "tissue": _single("unknown"),
        "cell_type": _multi_status("unknown"),
        "cell_line": _multi_status("not_applicable"),
        "treatment": _multi_status("not_applicable"),
        "tag": _multi_status("not_applicable"),
    }
    block.update(overrides)
    return block


def test_combination_query_keeps_historical_names_and_unresolved_traces() -> None:
    data, catalog = _stack()
    project = data.create_project("Treg project")
    human = _term(catalog, "species", "Homo sapiens")
    blood = _term(catalog, "tissue", "blood")
    treg = _term(catalog, "cell_type", "regulatory T cell")
    keratinocyte = _term(catalog, "cell_type", "keratinocyte")
    cuttag = _term(catalog, "assay", "CUT&Tag")
    atac = _term(catalog, "assay", "ATAC-seq")
    target = _term(catalog, "target", "H3K27ac")
    batch = _term(catalog, "batch", "batch-1")
    role = _term(catalog, "replicate_role", "biological replicate")
    ambiguous_a = _term(catalog, "assay", "CUT&Tag candidate")
    ambiguous_b = _term(catalog, "assay", "CUT&RUN candidate")

    def sample(key: str, catalog_block: dict | None):
        attributes = {"organism": "human"}
        if catalog_block is not None:
            attributes["catalog"] = catalog_block
        if catalog_block is None:
            return data.import_sample(
                project.project_id,
                stable_key=key,
                display_name=key,
                attributes=attributes,
            )
        with pytest.raises(CatalogBlockRequiresCatalogCommand):
            data.import_sample(
                project.project_id,
                stable_key=key,
                display_name=key,
                attributes=attributes,
            )
        created = data.import_sample(
            project.project_id,
            stable_key=key,
            display_name=key,
            attributes={"organism": "human"},
        )
        body = {
            "expected_revision_number": 1,
            "display_name": key,
            "catalog": catalog_block,
        }
        revised = catalog.revise_sample_catalog(
            **_write(f"catalog-{key}", body),
            sample_id=created.sample.sample_id,
            **body,
        )
        return revised

    treg_sample = sample(
        "treg",
        _catalog(species=_known(human), tissue=_known(blood), cell_type=_multi(treg)),
    )
    skin_sample = sample(
        "skin",
        _catalog(cell_type=_multi(keratinocyte)),
    )
    unknown_sample = sample("unknown-cell", _catalog())
    plain = data.import_sample(
        project.project_id,
        stable_key="plain",
        display_name="plain",
        attributes={"organism": "human"},
    )

    def experiment(key: str, assay: dict, sample_revision_id: str, **extra):
        body = {
            "stable_key": key,
            "display_name": key,
            "assay": assay,
            "target": _known(target),
            "batch": _known(batch),
            "samples": [
                {
                    "sample_revision_id": sample_revision_id,
                    "replicate_group": "group-a",
                    "replicate_role": _known(role),
                }
            ],
        }
        body.update(extra)
        return catalog.create_experiment(
            **_write(f"experiment-{key}", body),
            project_id=project.project_id,
            **body,
        )

    cut_treg = experiment("cut-treg", _known(cuttag), treg_sample["sample_revision_id"])
    atac_skin = experiment("atac-skin", _known(atac), skin_sample["sample_revision_id"])
    unknown_experiment = experiment(
        "cut-unknown",
        _known(cuttag),
        unknown_sample["sample_revision_id"],
    )
    ambiguous = experiment(
        "ambiguous-assay",
        {
            "status": "ambiguous",
            "candidates": [_pin(ambiguous_a), _pin(ambiguous_b)],
        },
        plain.revision.sample_revision_id,
    )

    listed = catalog.list_experiments(
        project_id=project.project_id,
        field_terms={"assay": [cuttag["term_id"]], "cell_type": [treg["term_id"]]},
        field_statuses={},
    )
    assert [item["experiment_id"] for item in listed["experiments"]] == [
        cut_treg["experiment_id"]
    ]
    assert "candidates" not in listed["experiments"][0]["assay"]
    assert listed["experiments"][0]["samples"][0]["cell_type"]["values"][0][
        "display_name"
    ] == "regulatory T cell"

    either_cell = catalog.list_experiments(
        project_id=project.project_id,
        field_terms={
            "cell_type": [treg["term_id"], keratinocyte["term_id"]],
            "assay": [atac["term_id"]],
        },
        field_statuses={},
    )
    assert [item["experiment_id"] for item in either_cell["experiments"]] == [
        atac_skin["experiment_id"]
    ]

    unknown_only = catalog.list_experiments(
        project_id=project.project_id,
        field_terms={},
        field_statuses={"cell_type": "unknown"},
    )
    assert [item["experiment_id"] for item in unknown_only["experiments"]] == [
        unknown_experiment["experiment_id"]
    ]

    absent = catalog.list_samples(
        project_id=project.project_id,
        field_terms={},
        field_statuses={},
        catalog_state="absent",
    )
    assert [item["sample_id"] for item in absent["samples"]] == [plain.sample.sample_id]

    detail = catalog.get_experiment_revision(
        ambiguous["experiment_id"],
        ambiguous["experiment_revision_id"],
    )
    assert [item["display_name"] for item in detail["assay"]["candidates"]] == [
        "CUT&Tag candidate",
        "CUT&RUN candidate",
    ]
    member = catalog.list_experiment_revisions(experiment_id=ambiguous["experiment_id"])
    assert "candidates" not in member["experiment_revisions"][0]["assay"]

    renamed = catalog.revise_term(
        **_write(
            "rename-treg",
            {
                "expected_revision_number": 1,
                "change_kind": "editorial",
                "display_name": "Treg",
                "namespace": None,
                "accession": None,
                "aliases": ["Treg"],
            },
        ),
        term_id=treg["term_id"],
        expected_revision_number=1,
        change_kind="editorial",
        display_name="Treg",
        namespace=None,
        accession=None,
        aliases=["Treg"],
    )
    assert renamed["revision"]["display_name"] == "Treg"
    again = catalog.list_experiments(
        project_id=project.project_id,
        field_terms={"cell_type": [treg["term_id"]]},
        field_statuses={},
    )
    assert (
        again["experiments"][0]["samples"][0]["cell_type"]["values"][0]["display_name"]
        == "regulatory T cell"
    )

    trace = catalog.list_experiment_run_results(
        project_id=project.project_id,
        experiment_revision_id=cut_treg["experiment_revision_id"],
    )
    assert trace == {
        "project_id": project.project_id,
        "experiment_revision_id": cut_treg["experiment_revision_id"],
        "trace_status": "unresolved",
        "reason": "experiment_snapshot_binding_absent",
        "results": [],
        "next_cursor": None,
    }

    copied = data.revise_sample(
        treg_sample["sample_id"],
        display_name="treg renamed",
        attributes={"organism": "human"},
    )
    assert '"display_name":"regulatory T cell"' not in copied.canonical_payload
    assert treg["term_id"] in copied.canonical_payload
    with pytest.raises(CatalogBlockRequiresCatalogCommand):
        data.revise_sample(
            treg_sample["sample_id"],
            display_name="treg renamed",
            attributes={"organism": "human", "catalog": {}},
        )


def test_same_key_replays_and_conflict_does_not_store_a_receipt() -> None:
    data, catalog = _stack()
    project = data.create_project("Pilot")
    body = {
        "term_type": "tag",
        "display_name": "spike-in",
        "namespace": None,
        "accession": None,
        "aliases": [],
    }
    raw = __import__("json").dumps(body).encode()
    created = catalog.create_term(
        actor_user_id="usr_" + "b" * 32,
        idempotency_key="term-1",
        raw_body=raw,
        **body,
    )
    replay = catalog.create_term(
        actor_user_id="usr_" + "b" * 32,
        idempotency_key="term-1",
        raw_body=raw,
        **body,
    )
    assert replay == created
    with pytest.raises(CatalogError) as conflict:
        catalog.create_term(
            actor_user_id="usr_" + "b" * 32,
            idempotency_key="term-1",
            raw_body=b'{"display_name":"other"}',
            **body,
        )
    assert conflict.value.code == "catalog_idempotency_conflict"

    sample = data.import_sample(
        project.project_id,
        stable_key="s1",
        display_name="S1",
        attributes={},
    )
    bad = {
        "expected_revision_number": 9,
        "display_name": "S1",
        "catalog": None,
    }
    with pytest.raises(CatalogError) as revision_conflict:
        catalog.revise_sample_catalog(
            actor_user_id="usr_" + "b" * 32,
            idempotency_key="sample-1",
            raw_body=__import__("json").dumps(bad).encode(),
            sample_id=sample.sample.sample_id,
            **bad,
        )
    assert revision_conflict.value.code == "catalog_revision_conflict"
    corrected = {
        "expected_revision_number": 1,
        "display_name": "S1",
        "catalog": None,
    }
    revised = catalog.revise_sample_catalog(
        actor_user_id="usr_" + "b" * 32,
        idempotency_key="sample-1",
        raw_body=__import__("json").dumps(corrected).encode(),
        sample_id=sample.sample.sample_id,
        **corrected,
    )
    assert revised["revision_number"] == 2
    assert revised["catalog_state"] == "absent"
