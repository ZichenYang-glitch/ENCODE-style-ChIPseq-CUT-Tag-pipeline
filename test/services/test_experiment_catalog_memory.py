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
    assert (
        listed["experiments"][0]["samples"][0]["cell_type"]["values"][0]["display_name"]
        == "regulatory T cell"
    )

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
        "ok": True,
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


def _plant_historical_catalog(data, revision, catalog_value):
    payload = __import__("json").loads(revision.canonical_payload)
    attributes = dict(payload["attributes"])
    attributes["catalog"] = catalog_value
    canonical = __import__(
        "encode_pipeline.platform.data_registry",
        fromlist=["canonical_sample_revision_payload"],
    ).canonical_sample_revision_payload(
        display_name=payload["display_name"],
        attributes=attributes,
    )
    digest = __import__(
        "encode_pipeline.platform.data_registry",
        fromlist=["build_sample_revision_payload_digest"],
    ).build_sample_revision_payload_digest(canonical)
    planted = __import__("dataclasses").replace(
        revision,
        canonical_payload=canonical,
        payload_digest=digest,
    )
    data._repository._revisions[revision.sample_revision_id] = planted
    return planted


def test_historical_catalog_keys_stay_free_attributes() -> None:
    data, catalog = _stack()
    project = data.create_project("Notebook")
    lookalike = {
        "schema_version": 1,
        "species": {"status": "unknown"},
        "tissue": {"status": "unknown"},
        "cell_type": {
            "status": "known",
            "values": [
                {
                    "term_id": "term_" + "a" * 32,
                    "term_revision_id": "trmv_" + "a" * 32,
                }
            ],
        },
        "cell_line": {"status": "not_applicable", "values": []},
        "treatment": {"status": "not_applicable", "values": []},
        "tag": {"status": "not_applicable", "values": []},
    }
    cell = _term(catalog, "cell_type", "regulatory T cell")
    species = _term(catalog, "species", "human")
    tissue = _term(catalog, "tissue", "blood")
    for key, value in (("notebook", "paper notebook 17"), ("lookalike", lookalike)):
        imported = data.import_sample(
            project.project_id,
            stable_key=key,
            display_name=key,
            attributes={"organism": "human"},
        )
        planted = _plant_historical_catalog(data, imported.revision, value)
        before = planted.canonical_payload
        before_digest = planted.payload_digest
        listed = catalog.list_samples(project_id=project.project_id)
        matched = [
            item
            for item in listed["samples"]
            if item["sample_id"] == imported.sample.sample_id
        ]
        assert matched[0]["catalog_state"] == "absent"
        assert "cell_type" not in matched[0]
        historical = catalog.get_sample_revision(
            imported.sample.sample_id,
            planted.sample_revision_id,
        )
        assert historical["catalog_state"] == "absent"
        renamed = data.revise_sample(
            imported.sample.sample_id,
            display_name=f"{key} renamed",
            attributes={"organism": "human"},
        )
        assert renamed.revision_number == 2
        stored = data._repository.get_sample_revision(planted.sample_revision_id)
        assert stored.canonical_payload == before
        assert stored.payload_digest == before_digest
        catalog._repository.rebuild_sample_projection(planted.sample_revision_id)
        rebuilt = data._repository.get_sample_revision(planted.sample_revision_id)
        assert rebuilt.canonical_payload == before
        assert rebuilt.payload_digest == before_digest
        assert (
            catalog.get_sample_revision(
                imported.sample.sample_id,
                planted.sample_revision_id,
            )["catalog_state"]
            == "absent"
        )
        classified = catalog.revise_sample_catalog(
            **_write(f"classify-{key}", {"catalog": "pending"}),
            sample_id=imported.sample.sample_id,
            expected_revision_number=2,
            display_name=None,
            catalog={
                "schema_version": 1,
                "species": {
                    "status": "known",
                    "term_id": species["term_id"],
                    "term_revision_id": species["revision"]["term_revision_id"],
                },
                "tissue": {
                    "status": "known",
                    "term_id": tissue["term_id"],
                    "term_revision_id": tissue["revision"]["term_revision_id"],
                },
                "cell_type": {
                    "status": "known",
                    "values": [
                        {
                            "term_id": cell["term_id"],
                            "term_revision_id": cell["revision"]["term_revision_id"],
                        }
                    ],
                },
                "cell_line": {"status": "not_applicable", "values": []},
                "treatment": {"status": "not_applicable", "values": []},
                "tag": {"status": "not_applicable", "values": []},
            },
        )
        assert classified["catalog_state"] == "present"
        assert classified["revision_number"] == 3
        assert classified["cell_type"]["values"][0]["term_id"] == cell["term_id"]
        assert (
            catalog.get_sample_revision(
                imported.sample.sample_id,
                planted.sample_revision_id,
            )["catalog_state"]
            == "absent"
        )
        assert (
            data._repository.get_sample_revision(
                planted.sample_revision_id
            ).canonical_payload
            == before
        )


def test_receipt_target_rejects_a_different_object() -> None:
    data, catalog = _stack()
    actor = "usr_" + "b" * 32
    first = _term(catalog, "tag", "alpha")
    second = _term(catalog, "tag", "beta")
    archive = {"expected_revision_number": 1}
    raw = __import__("json").dumps(archive).encode()
    archived = catalog.archive_term(
        actor_user_id=actor,
        idempotency_key="archive-once",
        raw_body=raw,
        term_id=first["term_id"],
        expected_revision_number=1,
    )
    with pytest.raises(CatalogError) as conflict:
        catalog.archive_term(
            actor_user_id=actor,
            idempotency_key="archive-once",
            raw_body=raw,
            term_id=second["term_id"],
            expected_revision_number=1,
        )
    assert conflict.value.status_code == 409
    assert conflict.value.code == "catalog_idempotency_conflict"
    assert catalog.get_term(second["term_id"])["archived_at"] is None
    replay = catalog.archive_term(
        actor_user_id=actor,
        idempotency_key="archive-once",
        raw_body=raw,
        term_id=first["term_id"],
        expected_revision_number=1,
    )
    assert replay["term_id"] == archived["term_id"]
    assert replay["archived_at"] == archived["archived_at"]

    left = data.create_project("Left")
    right = data.create_project("Right")
    sample = data.import_sample(
        left.project_id,
        stable_key="shared",
        display_name="Shared",
        attributes={},
    )
    assay = _term(catalog, "assay", "CUT&Tag")
    target = _term(catalog, "target", "H3K27ac")
    batch = _term(catalog, "batch", "batch-1")
    role = _term(catalog, "replicate_role", "biological replicate")

    def known(term: dict) -> dict:
        return {
            "status": "known",
            "term_id": term["term_id"],
            "term_revision_id": term["revision"]["term_revision_id"],
        }

    body = {
        "stable_key": "shared-exp",
        "display_name": "Shared",
        "assay": known(assay),
        "target": known(target),
        "batch": known(batch),
        "samples": [
            {
                "sample_revision_id": sample.revision.sample_revision_id,
                "replicate_role": known(role),
            }
        ],
    }
    created = catalog.create_experiment(
        actor_user_id=actor,
        idempotency_key="create-once",
        raw_body=__import__("json").dumps(body).encode(),
        project_id=left.project_id,
        **body,
    )
    with pytest.raises(CatalogError) as project_conflict:
        catalog.create_experiment(
            actor_user_id=actor,
            idempotency_key="create-once",
            raw_body=__import__("json").dumps(body).encode(),
            project_id=right.project_id,
            **body,
        )
    assert project_conflict.value.code == "catalog_idempotency_conflict"
    assert catalog.list_experiments(project_id=right.project_id)["experiments"] == []
    again = catalog.create_experiment(
        actor_user_id=actor,
        idempotency_key="create-once",
        raw_body=__import__("json").dumps(body).encode(),
        project_id=left.project_id,
        **body,
    )
    assert again["experiment_id"] == created["experiment_id"]
    assert again["experiment_revision_id"] == created["experiment_revision_id"]
    assert catalog.list_experiments(project_id=left.project_id)["ok"] is True


def test_memory_adoption_commits_one_subset_and_replays() -> None:
    data, catalog = _stack()
    project = data.create_project("Memory adoption")
    assay = _term(catalog, "assay", "CUT&Tag")
    role = _term(catalog, "replicate_role", "biological replicate")
    curated = data.import_sample(
        project.project_id,
        stable_key="curated",
        display_name="curated",
        attributes={"legacy_note": "keep-me", "organism": "human"},
    )
    catalog_body = {
        "expected_revision_number": 1,
        "display_name": "curated",
        "catalog": _catalog(),
    }
    sample = catalog.revise_sample_catalog(
        **_write("memory-sample", catalog_body),
        sample_id=curated.sample.sample_id,
        **catalog_body,
    )
    experiment_body = {
        "project_id": project.project_id,
        "stable_key": "memory-exp",
        "display_name": "memory-exp",
        "assay": _known(assay),
        "target": _single("unknown"),
        "batch": _single("unknown"),
        "samples": [
            {
                "sample_revision_id": sample["sample_revision_id"],
                "replicate_role": _known(role),
            }
        ],
    }
    experiment = catalog.create_experiment(
        **_write("memory-experiment", experiment_body),
        **experiment_body,
    )
    from encode_pipeline.platform.catalog_curation_preview import (
        M2_RULE_VERSION,
        M2_VOCABULARY_VERSION,
    )

    item = {
        "item_id": "item_1",
        "target_kind": "experiment",
        "target_id": experiment["experiment_id"],
        "base_revision_id": experiment["experiment_revision_id"],
        "base_revision_number": experiment["revision_number"],
        "allowed_fields": ["display_name"],
        "before": {"display_name": "memory-exp"},
        "after": {"display_name": "memory-exp adopted"},
        "citations": {
            "rule_version": M2_RULE_VERSION,
            "vocabulary_version": M2_VOCABULARY_VERSION,
            "reason_codes": ["MISSING_EVIDENCE_UNKNOWN"],
            "source_refs": [
                {
                    "kind": "experiment_revision",
                    "ref_id": experiment["experiment_revision_id"],
                    "evidence_token": "memory-exp",
                }
            ],
        },
    }
    register_body = {"source": "offline_suggestion", "items": [item]}
    registered = catalog.register_proposal(
        **_write("memory-register", register_body),
        project_id=project.project_id,
        source="offline_suggestion",
        items=[item],
        rule_version=None,
        vocabulary_version=None,
        projection_version=None,
    )
    preview = catalog.preview_proposal(
        actor_user_id="usr_" + "a" * 32,
        raw_body=b"{}",
        proposal_id=registered["proposal_id"],
        selected_item_ids=["item_1"],
    )
    adopt_body = {
        "preview_digest": preview["preview_digest"],
        "selected_item_ids": ["item_1"],
        "expected_revisions": preview["expected_revisions"],
    }
    adopted = catalog.adopt_proposal(
        **_write("memory-adopt", adopt_body),
        proposal_id=registered["proposal_id"],
        selected_item_ids=["item_1"],
        expected_revisions=preview["expected_revisions"],
        preview_digest=preview["preview_digest"],
    )
    assert adopted["replayed"] is False
    replay = catalog.adopt_proposal(
        **_write("memory-adopt", adopt_body),
        proposal_id=registered["proposal_id"],
        selected_item_ids=["item_1"],
        expected_revisions=preview["expected_revisions"],
        preview_digest=preview["preview_digest"],
    )
    assert replay["replayed"] is True
    assert replay["receipt_id"] == adopted["receipt_id"]
    with pytest.raises(CatalogError) as second:
        catalog.adopt_proposal(
            **_write("memory-adopt-again", adopt_body),
            proposal_id=registered["proposal_id"],
            selected_item_ids=["item_1"],
            expected_revisions=preview["expected_revisions"],
            preview_digest=preview["preview_digest"],
        )
    assert second.value.code == "catalog_proposal_invalidated"
    revisions = catalog.list_experiment_revisions(
        experiment_id=experiment["experiment_id"]
    )
    assert len(revisions["experiment_revisions"]) == 2
