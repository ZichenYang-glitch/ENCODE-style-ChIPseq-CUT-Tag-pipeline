"""Behavior shared by the in-memory and SQLite experiment catalogs."""

from datetime import datetime, timedelta, timezone
from hashlib import sha256
from itertools import count
import json
from unittest.mock import Mock

import pytest

from encode_pipeline.persistence import create_database_engine, create_session_factory
from encode_pipeline.persistence.data_registry import SqlAlchemyDataRegistryRepository
from encode_pipeline.persistence.experiment_catalog import (
    SqlAlchemyExperimentCatalogRepository,
)
from encode_pipeline.persistence.migrations import upgrade_database
from encode_pipeline.persistence.models import UserAccountRow
from encode_pipeline.platform.catalog_proposal_adoption import stable_proposal_id
from encode_pipeline.platform.experiment_catalog import CatalogError
from encode_pipeline.services.data_registry import DataRegistryService
from encode_pipeline.services.data_registry_repositories import (
    InMemoryDataRegistryRepository,
)
from encode_pipeline.services.experiment_catalog import ExperimentCatalogService
from encode_pipeline.services.experiment_catalog_repositories import (
    InMemoryExperimentCatalogRepository,
)


START = datetime(2026, 11, 18, tzinfo=timezone.utc)
ACTOR = "usr_" + "c" * 32
UNKNOWN = {"status": "unknown"}


@pytest.fixture(params=("memory", "sqlite"))
def catalog_stack(request, tmp_path):
    ticks = count()

    def now():
        return START + timedelta(seconds=next(ticks))

    if request.param == "memory":
        registry = InMemoryDataRegistryRepository(legacy_created_at=START)
        repository = InMemoryExperimentCatalogRepository(
            registry, is_administrator=lambda user_id: user_id == ACTOR
        )
        yield (
            DataRegistryService(repository=registry, now_factory=now),
            (ExperimentCatalogService(repository, now_factory=now)),
        )
        return

    url = f"sqlite:///{tmp_path / 'catalog.db'}"
    upgrade_database(url)
    engine = create_database_engine(url)
    try:
        factory = create_session_factory(engine)
        with factory.begin() as session:
            session.add(
                UserAccountRow(
                    user_id=ACTOR,
                    username="catalog-admin",
                    role="administrator",
                    status="enabled",
                    password_hash="fixture",
                    created_at=START,
                    updated_at=START,
                    password_changed_at=START,
                )
            )
        yield (
            DataRegistryService(
                repository=SqlAlchemyDataRegistryRepository(factory), now_factory=now
            ),
            ExperimentCatalogService(
                SqlAlchemyExperimentCatalogRepository(factory), now_factory=now
            ),
        )
    finally:
        engine.dispose()


def _write(command, key, **body):
    return command(
        actor_user_id=ACTOR,
        idempotency_key=key,
        raw_body=json.dumps(body, sort_keys=True).encode(),
        **body,
    )


def _term(catalog, term_type, name, **overrides):
    body = {
        "term_type": term_type,
        "display_name": name,
        "namespace": None,
        "accession": None,
        "aliases": [],
    } | overrides
    return _write(
        catalog.create_term,
        f"term-{term_type}-{sha256(name.encode()).hexdigest()}",
        **body,
    )


def _pin(term):
    return {
        "term_id": term["term_id"],
        "term_revision_id": term["revision"]["term_revision_id"],
    }


def _material(**overrides):
    return {
        "schema_version": 1,
        "species": UNKNOWN,
        "tissue": UNKNOWN,
        **{
            field: {"status": "not_applicable", "values": []}
            for field in ("cell_type", "cell_line", "treatment", "tag")
        },
    } | overrides


def _experiment(sample_revision_id, **overrides):
    return {
        "display_name": "Pilot experiment",
        "assay": UNKNOWN,
        "target": UNKNOWN,
        "batch": UNKNOWN,
        "samples": [
            {
                "sample_revision_id": sample_revision_id,
                "replicate_role": UNKNOWN,
                "replicate_group": "",
            }
        ],
    } | overrides


@pytest.mark.parametrize(
    ("change", "error"),
    [
        ("rename_term", "catalog_revision_conflict"),
        ("archive_term", "catalog_term_archived"),
        ("revise_experiment", "catalog_revision_conflict"),
        ("archive_experiment", "catalog_experiment_archived"),
    ],
)
def test_model_checks_pins_in_the_proposal_registration_transaction(
    catalog_stack, monkeypatch, change, error
):
    data, catalog = catalog_stack
    project = data.create_project("Synthetic review")
    sample = data.import_sample(
        project.project_id, stable_key="sample", display_name="Sample", attributes={}
    )
    fields = _experiment(sample.revision.sample_revision_id)
    experiment = _write(
        catalog.create_experiment,
        "experiment",
        project_id=project.project_id,
        stable_key="synthetic",
        **fields,
    )
    target = _term(catalog, "target", "H2A.Z")
    model = Mock()
    model.classify.return_value = "h2az"
    monkeypatch.setattr(catalog, "_catalog_model", model)
    body = {
        "project_id": project.project_id,
        "experiment_id": experiment["experiment_id"],
        "experiment_revision_id": experiment["experiment_revision_id"],
        "field": "target",
        "reason_code": "TARGET_TERM_MATCH",
        "clue": "H2A.Z",
        "vocabulary_term_revision_id": target["revision"]["term_revision_id"],
    }

    with pytest.raises(CatalogError) as collision:
        _write(catalog.generate_classification_suggestion, "experiment", **body)
    assert collision.value.code == "catalog_idempotency_conflict"
    model.classify.assert_not_called()

    register = catalog._repository.register_proposal

    def register_after_change(*args, **kwargs):
        term_args = {"term_id": target["term_id"], "expected_revision_number": 1}
        experiment_args = {
            "experiment_id": experiment["experiment_id"],
            "expected_revision_number": 1,
        }
        if change == "rename_term":
            _write(
                catalog.revise_term,
                "edit",
                **term_args,
                change_kind="editorial",
                display_name="Updated target",
                aliases=[],
                namespace=None,
                accession=None,
            )
        elif change == "archive_term":
            _write(catalog.archive_term, "edit", **term_args)
        elif change == "revise_experiment":
            _write(
                catalog.revise_experiment,
                "edit",
                **experiment_args,
                **(fields | {"display_name": "Updated experiment"}),
            )
        else:
            _write(catalog.archive_experiment, "edit", **experiment_args)
        return register(*args, **kwargs)

    monkeypatch.setattr(catalog._repository, "register_proposal", register_after_change)
    with pytest.raises(CatalogError) as conflict:
        _write(catalog.generate_classification_suggestion, "generate", **body)
    assert conflict.value.code == error
    model.classify.assert_called_once()
    with pytest.raises(CatalogError) as missing:
        catalog.get_proposal(
            proposal_id=stable_proposal_id(ACTOR, "generate", project.project_id),
            include_candidates=True,
        )
    assert missing.value.code == "catalog_not_found"


def test_implicit_candidates_come_from_current_target_pins_in_the_selected_project(
    catalog_stack, monkeypatch
):
    data, catalog = catalog_stack
    local = data.create_project("Selected")
    other = data.create_project("Other")

    def record(project, key, target=None):
        sample = data.import_sample(
            project.project_id, stable_key=key, display_name=key, attributes={}
        )
        return _write(
            catalog.create_experiment,
            key,
            project_id=project.project_id,
            stable_key=key,
            **_experiment(
                sample.revision.sample_revision_id,
                target={"status": "known", **_pin(target)} if target else UNKNOWN,
            ),
        )

    first = _term(catalog, "target", "H3K4me3")
    second = _term(catalog, "target", "H2A.Z")
    record(other, "foreign", second)
    record(local, "local", first)
    selected = record(local, "selected")
    model = Mock()
    model.classify.return_value = "insufficient_evidence"
    monkeypatch.setattr(catalog, "_catalog_model", model)
    body = {
        "project_id": local.project_id,
        "experiment_id": selected["experiment_id"],
        "experiment_revision_id": selected["experiment_revision_id"],
        "field": "target",
        "reason_code": "TARGET_TERM_MATCH",
        "clue": "H3K4",
        "vocabulary_term_revision_id": None,
    }
    result = _write(catalog.generate_classification_suggestion, "first", **body)
    assert result["outcome"] == "insufficient_evidence"
    assert model.classify.call_args.kwargs["candidates"] == (("h3k4me3", "H3K4me3"),)
    _write(
        catalog.revise_term,
        "edit",
        term_id=first["term_id"],
        expected_revision_number=1,
        change_kind="editorial",
        display_name="Revised label",
        aliases=[],
        namespace=None,
        accession=None,
    )
    with pytest.raises(CatalogError) as missing:
        _write(catalog.generate_classification_suggestion, "second", **body)
    assert missing.value.path == "vocabulary_term_revision_id"
    model.classify.assert_called_once()


def test_term_paging_rejects_removed_cursor_and_keeps_archived_history(catalog_stack):
    data, catalog = catalog_stack
    project = data.create_project("Pilot")
    sample = data.import_sample(
        project.project_id, stable_key="pilot", display_name="Pilot", attributes={}
    )
    original = _term(catalog, "assay", "CUT&Tag", aliases=["CUT and Tag"])
    other = _term(catalog, "assay", "ATAC-seq")
    tag = _term(catalog, "tag", "control")

    def term_ids(page):
        return [term["term_id"] for term in page["terms"]]

    query = {
        "term_type": "assay",
        "limit": 1,
        "created_from": START + timedelta(seconds=2),
    }
    first = catalog.list_terms(**query, created_before=START + timedelta(seconds=4))
    assert term_ids(first) == [original["term_id"]]
    assert term_ids(
        catalog.list_terms(
            **query,
            created_before=START + timedelta(seconds=4),
            cursor=first["next_cursor"],
        )
    ) == [other["term_id"]]
    assert catalog.list_terms(term_type="tag")["terms"][0]["term_id"] == tag["term_id"]

    experiment = _write(
        catalog.create_experiment,
        "create-pilot",
        project_id=project.project_id,
        stable_key="pilot",
        **_experiment(
            sample.revision.sample_revision_id,
            assay={"status": "known", **_pin(original)},
        ),
    )
    command = {
        "term_id": original["term_id"],
        "expected_revision_number": 1,
        "change_kind": "editorial",
        "display_name": "CUT&Tag (edited)",
        "namespace": "assays",
        "accession": "cat-1",
        "aliases": ["CUT and Tag"],
    }
    revised = _write(catalog.revise_term, "revise-assay", **command)
    assert revised["revision"]["revision_number"] == 2
    assert revised["revision"]["display_name"] == command["display_name"]
    assert catalog.get_term(original["term_id"]) == revised

    archive = {"term_id": original["term_id"], "expected_revision_number": 2}
    archived = _write(catalog.archive_term, "archive-assay", **archive)
    with pytest.raises(CatalogError, match="already archived") as failure:
        _write(catalog.archive_term, "archive-assay-again", **archive)
    assert failure.value.code == "catalog_already_archived"
    with pytest.raises(CatalogError) as failure:
        _write(catalog.revise_term, "edit-archived", **(command | archive))
    assert failure.value.code == "catalog_term_archived"

    with pytest.raises(CatalogError) as failure:
        catalog.list_terms(
            **query,
            created_before=START + timedelta(seconds=4),
            cursor=first["next_cursor"],
        )
    assert failure.value.code == "catalog_cursor_not_found"
    assert term_ids(catalog.list_terms(term_type="assay")) == [other["term_id"]]
    assert term_ids(catalog.list_terms(term_type="assay", include_archived=True)) == [
        original["term_id"],
        other["term_id"],
    ]
    assert archived["archived_at"] is not None
    historical = catalog.get_experiment_revision(
        experiment["experiment_id"], experiment["experiment_revision_id"]
    )
    assert historical["assay"] == experiment["assay"]
    assert historical["assay"]["display_name"] == "CUT&Tag"


def test_sample_filters_history_and_ownership(catalog_stack):
    data, catalog = catalog_stack
    project = data.create_project("Pilot")
    left = data.import_sample(
        project.project_id, stable_key="left", display_name="Left", attributes={}
    )
    right = data.import_sample(
        project.project_id, stable_key="right", display_name="Right", attributes={}
    )
    species = _term(catalog, "species", "Homo sapiens")
    cell = _term(catalog, "cell_type", "T cell")
    change = {
        "sample_id": left.sample.sample_id,
        "expected_revision_number": 1,
        "display_name": None,
        "catalog": _material(
            species={"status": "known", **_pin(species)},
            cell_type={"status": "ambiguous", "values": [], "candidates": [_pin(cell)]},
        ),
    }
    classified = _write(catalog.revise_sample_catalog, "classify-left", **change)
    assert (
        _write(catalog.revise_sample_catalog, "classify-left", **change) == classified
    )
    assert classified["display_name"] == "Left"
    assert classified["cell_type"]["candidates"][0]["term_id"] == cell["term_id"]

    def sample_ids(**query):
        page = catalog.list_samples(project_id=project.project_id, **query)
        return [row["sample_id"] for row in page["samples"]]

    assert sample_ids(catalog_state="present") == [left.sample.sample_id]
    assert sample_ids(catalog_state="absent") == [right.sample.sample_id]
    assert sample_ids(field_terms={"species": [species["term_id"]]}) == [
        left.sample.sample_id
    ]
    assert sample_ids(field_statuses={"cell_type": "ambiguous"}) == [
        left.sample.sample_id
    ]
    assert sample_ids(field_statuses={"cell_type": "unknown"}) == []
    assert sample_ids(
        created_from=START + timedelta(seconds=2),
        created_before=START + timedelta(seconds=3),
    ) == [right.sample.sample_id]
    assert sample_ids(created_before=START + timedelta(seconds=1)) == []
    page = catalog.list_samples(project_id=project.project_id, limit=1)
    assert [row["sample_id"] for row in page["samples"]] == [left.sample.sample_id]
    assert sample_ids(limit=1, cursor=page["next_cursor"]) == [right.sample.sample_id]

    owner = {"sample_id": left.sample.sample_id, "limit": 1}
    history = catalog.list_sample_revisions(**owner)
    next_page = catalog.list_sample_revisions(**owner, cursor=history["next_cursor"])
    assert [
        (row["revision_number"], row["catalog_state"])
        for row in (history["sample_revisions"] + next_page["sample_revisions"])
    ] == [(1, "absent"), (2, "present")]
    assert "candidates" not in next_page["sample_revisions"][0]["cell_type"]
    old_revision = left.revision.sample_revision_id
    for wrong_owner, revision_id in (
        (right.sample.sample_id, old_revision),
        (left.sample.sample_id, "smpr_" + "f" * 32),
    ):
        with pytest.raises(CatalogError) as failure:
            catalog.get_sample_revision(wrong_owner, revision_id)
        assert failure.value.code == "catalog_not_found"
    with pytest.raises(CatalogError) as failure:
        catalog.list_sample_revisions(
            sample_id=right.sample.sample_id, cursor=history["next_cursor"]
        )
    assert failure.value.code == "catalog_cursor_invalid"


def test_experiment_revision_and_archive_are_scoped_and_replayable(catalog_stack):
    data, catalog = catalog_stack
    project = data.create_project("Pilot")
    foreign = data.create_project("Other project")
    local = data.import_sample(
        project.project_id, stable_key="pilot", display_name="Pilot", attributes={}
    )
    other = data.import_sample(
        foreign.project_id, stable_key="other", display_name="Other", attributes={}
    )
    scope = {"project_id": project.project_id, "stable_key": "pilot"}
    for revision_id in (other.revision.sample_revision_id, "smpr_" + "f" * 32):
        with pytest.raises(CatalogError) as failure:
            _write(
                catalog.create_experiment,
                "create-pilot",
                **scope,
                **_experiment(revision_id),
            )
        assert failure.value.code == "catalog_payload_invalid"
    assert catalog.list_experiments(project_id=project.project_id)["experiments"] == []

    body = _experiment(local.revision.sample_revision_id)
    created = _write(catalog.create_experiment, "create-pilot", **scope, **body)
    with pytest.raises(CatalogError) as failure:
        _write(catalog.create_experiment, "duplicate-stable-key", **scope, **body)
    assert failure.value.code == "catalog_revision_conflict"
    experiment_id = created["experiment_id"]
    revision_id = created["experiment_revision_id"]
    assert created["samples"][0].get("replicate_group") in (None, "")

    revise = body | {
        "experiment_id": experiment_id,
        "expected_revision_number": 1,
        "display_name": "Edited pilot experiment",
    }
    revised = _write(catalog.revise_experiment, "revise-pilot", **revise)
    assert revised["revision_number"] == 2
    assert revised["experiment_revision_id"] != revision_id
    assert _write(catalog.revise_experiment, "revise-pilot", **revise) == revised

    unrelated_id = "exp_" + "f" * 32
    for command in (catalog.revise_experiment, catalog.archive_experiment):
        payload = (
            revise
            if command == catalog.revise_experiment
            else {"expected_revision_number": 2}
        )
        with pytest.raises(CatalogError) as failure:
            _write(
                command, "unowned-write", **(payload | {"experiment_id": unrelated_id})
            )
        assert failure.value.code == "catalog_not_found"
    with pytest.raises(CatalogError) as failure:
        catalog.get_experiment_revision(unrelated_id, revision_id)
    assert failure.value.code == "catalog_not_found"

    archive = {"experiment_id": experiment_id, "expected_revision_number": 2}
    archived = _write(catalog.archive_experiment, "archive-pilot", **archive)
    assert _write(catalog.archive_experiment, "archive-pilot", **archive) == archived
    assert archived["archived_at"] is not None
    with pytest.raises(CatalogError) as failure:
        _write(catalog.archive_experiment, "archive-pilot-again", **archive)
    assert failure.value.code == "catalog_already_archived"
    with pytest.raises(CatalogError) as failure:
        _write(catalog.revise_experiment, "edit-archived-pilot", **(revise | archive))
    assert failure.value.code == "catalog_experiment_archived"
    assert catalog.list_experiments(project_id=project.project_id)["experiments"] == []
    assert [
        (row["revision_number"], row["display_name"])
        for row in catalog.list_experiment_revisions(experiment_id=experiment_id)[
            "experiment_revisions"
        ]
    ] == [(1, body["display_name"]), (2, revise["display_name"])]
    assert (
        catalog.get_experiment_revision(experiment_id, revision_id)["display_name"]
        == body["display_name"]
    )


def test_catalog_traces_distinguish_empty_history_from_wrong_project(catalog_stack):
    data, catalog = catalog_stack
    project = data.create_project("Pilot")
    foreign = data.create_project("Other project")
    sample = data.import_sample(
        project.project_id, stable_key="pilot", display_name="Pilot", attributes={}
    )
    sample_revision_id = sample.revision.sample_revision_id
    experiment = _write(
        catalog.create_experiment,
        "create-pilot",
        project_id=project.project_id,
        stable_key="pilot",
        **_experiment(sample_revision_id),
    )
    for trace, revision in (
        (catalog.list_sample_run_results, {"sample_revision_id": sample_revision_id}),
        (
            catalog.list_experiment_run_results,
            {"experiment_revision_id": experiment["experiment_revision_id"]},
        ),
    ):
        empty = trace(project_id=project.project_id, **revision)
        assert empty["results"] == []
        assert empty["next_cursor"] is None
        with pytest.raises(CatalogError) as failure:
            trace(project_id=foreign.project_id, **revision)
        assert (failure.value.status_code, failure.value.code) == (
            404,
            "catalog_not_found",
        )

    data.archive_project(project.project_id)
    with pytest.raises(CatalogError) as failure:
        catalog.list_experiment_revisions(experiment_id=experiment["experiment_id"])
    assert failure.value.code == "catalog_not_found"


@pytest.fixture
def proposal_stack(catalog_stack):
    from encode_pipeline.platform.catalog_curation_preview import (
        M2_RULE_VERSION,
        M2_VOCABULARY_VERSION,
    )

    data, catalog = catalog_stack
    project = data.create_project("Curation")
    sample = data.import_sample(
        project.project_id,
        stable_key="pilot",
        display_name="Pilot",
        attributes={"notebook": "pilot entry"},
    )
    classified = _write(
        catalog.revise_sample_catalog,
        "classify",
        sample_id=sample.sample.sample_id,
        expected_revision_number=1,
        display_name=None,
        catalog=_material(),
    )
    experiment = _write(
        catalog.create_experiment,
        "create",
        project_id=project.project_id,
        stable_key="pilot",
        **_experiment(classified["sample_revision_id"]),
    )
    tag = _term(catalog, "tag", "curated")
    batch = _term(catalog, "batch", "batch-1")
    common = {
        "rule_version": M2_RULE_VERSION,
        "vocabulary_version": M2_VOCABULARY_VERSION,
    }
    experiment_item = {
        "item_id": "experiment",
        "target_kind": "experiment",
        "target_id": experiment["experiment_id"],
        "base_revision_id": experiment["experiment_revision_id"],
        "base_revision_number": 1,
        "allowed_fields": ["batch"],
        "before": {"batch": UNKNOWN},
        "after": {"batch": {"status": "known", **_pin(batch)}},
        "citations": common
        | {
            "reason_codes": ["BATCH_TERM_MATCH"],
            "source_refs": [
                {
                    "kind": "experiment_revision",
                    "ref_id": experiment["experiment_revision_id"],
                },
                {
                    "kind": "catalog_projection_row",
                    "ref_id": experiment["experiment_id"],
                },
                {"kind": "term_revision", "ref_id": _pin(batch)["term_revision_id"]},
            ],
        },
    }
    sample_item = {
        "item_id": "sample",
        "target_kind": "sample_catalog",
        "target_id": sample.sample.sample_id,
        "base_revision_id": classified["sample_revision_id"],
        "base_revision_number": 2,
        "allowed_fields": ["display_name", "tag"],
        "before": {
            "display_name": "Pilot",
            "tag": {"status": "not_applicable", "values": []},
        },
        "after": {
            "display_name": "Curated pilot",
            "tag": {"status": "known", "values": [_pin(tag)]},
        },
        "citations": common
        | {
            "reason_codes": ["TAG_TERM_MATCH"],
            "source_refs": [
                {"kind": "sample_revision", "ref_id": classified["sample_revision_id"]},
                {"kind": "catalog_projection_row", "ref_id": sample.sample.sample_id},
                {
                    "kind": "catalog_projection_row",
                    "ref_id": classified["sample_revision_id"],
                },
            ],
        },
    }
    return data, catalog, project.project_id, experiment_item, sample_item


def _register(catalog, project_id, *items, key="register-proposal"):
    return _write(
        catalog.register_proposal,
        key,
        project_id=project_id,
        source="offline_suggestion",
        items=list(items),
        rule_version=None,
        vocabulary_version=None,
        projection_version=None,
    )


def _preview(catalog, proposal_id, *item_ids, actor=ACTOR):
    return catalog.preview_proposal(
        actor_user_id=actor,
        raw_body=json.dumps({"selected_item_ids": item_ids}).encode(),
        proposal_id=proposal_id,
        selected_item_ids=list(item_ids),
    )


def _adopt(catalog, preview, key="adopt-proposal"):
    return _write(
        catalog.adopt_proposal,
        key,
        **{
            field: preview[field]
            for field in (
                "proposal_id",
                "selected_item_ids",
                "expected_revisions",
                "preview_digest",
            )
        },
    )


def _proposal_receipts(catalog, project_id, proposal_id=None):
    return catalog.list_adoption_receipts(
        actor_user_id=ACTOR, project_id=project_id, proposal_id=proposal_id
    )["receipts"]


def test_adoption_replays_both_writes_and_keeps_the_original_revisions(proposal_stack):
    data, catalog, project_id, experiment, sample = proposal_stack
    proposal = _register(catalog, project_id, experiment, sample)
    assert _register(catalog, project_id, experiment, sample) == proposal
    proposal_id = proposal["proposal_id"]
    preview = _preview(catalog, proposal_id, "sample", "experiment")
    assert preview["has_conflicts"] is False
    assert _proposal_receipts(catalog, project_id) == []
    assert (
        catalog.list_sample_revisions(sample_id=sample["target_id"])[
            "sample_revisions"
        ][-1]["revision_number"]
        == 2
    )

    receipt = _adopt(catalog, preview)
    assert receipt["replayed"] is False
    assert _adopt(catalog, preview) == receipt | {"replayed": True}
    assert _proposal_receipts(catalog, project_id) == [receipt]
    assert _proposal_receipts(catalog, project_id, proposal_id) == [receipt]
    assert _proposal_receipts(catalog, data.create_project("Other").project_id) == []
    assert (
        catalog.get_proposal(proposal_id=proposal_id, include_candidates=False)[
            "proposal_state"
        ]
        == "adopted"
    )
    for retry in (
        lambda: _preview(catalog, proposal_id, "sample"),
        lambda: _adopt(catalog, preview, key="fresh-adoption-key"),
    ):
        with pytest.raises(CatalogError) as failure:
            retry()
        assert failure.value.code == "catalog_proposal_invalidated"

    experiment_history = catalog.list_experiment_revisions(
        experiment_id=experiment["target_id"]
    )["experiment_revisions"]
    assert [row["revision_number"] for row in experiment_history] == [1, 2]
    old, new = experiment_history
    assert old["batch"] == UNKNOWN
    assert new["batch"]["display_name"] == "batch-1"
    assert new["samples"] == old["samples"]
    sample_history = catalog.list_sample_revisions(sample_id=sample["target_id"])[
        "sample_revisions"
    ]
    assert [row["revision_number"] for row in sample_history] == [1, 2, 3]
    assert sample_history[1]["tag"] == {"status": "not_applicable", "values": []}
    assert sample_history[2]["tag"]["values"][0]["display_name"] == "curated"
    assert sample_history[2]["display_name"] == "Curated pilot"
    stored = data._repository.get_sample_revision(
        sample_history[2]["sample_revision_id"]
    )
    assert (
        json.loads(stored.canonical_payload)["attributes"]["notebook"] == "pilot entry"
    )


def test_stale_sample_adoption_keeps_other_target_at_its_original_revision(
    proposal_stack,
):
    _, catalog, project_id, experiment, sample = proposal_stack
    proposal_id = _register(catalog, project_id, experiment, sample)["proposal_id"]
    preview = _preview(catalog, proposal_id, "experiment", "sample")
    _write(
        catalog.revise_sample_catalog,
        "manual-sample",
        sample_id=sample["target_id"],
        expected_revision_number=2,
        display_name="Manual correction",
        catalog=_material(),
    )
    with pytest.raises(CatalogError) as failure:
        _adopt(catalog, preview)
    assert failure.value.code == "catalog_revision_conflict"
    assert _proposal_receipts(catalog, project_id, proposal_id) == []
    assert (
        catalog.list_experiment_revisions(experiment_id=experiment["target_id"])[
            "experiment_revisions"
        ][0]["batch"]
        == UNKNOWN
    )
    assert (
        len(
            catalog.list_experiment_revisions(experiment_id=experiment["target_id"])[
                "experiment_revisions"
            ]
        )
        == 1
    )
    with pytest.raises(CatalogError) as failure:
        _preview(catalog, "cprop_" + "0" * 32, "experiment")
    assert failure.value.code == "catalog_not_found"


def test_proposal_preview_resolves_projection_sources_in_its_own_project(
    proposal_stack,
):
    data, catalog, project_id, experiment, sample = proposal_stack
    outside = data.create_project("Other")
    foreign = data.import_sample(
        outside.project_id, stable_key="foreign", display_name="Foreign", attributes={}
    )
    for number, (kind, ref) in enumerate(
        (
            ("catalog_projection_row", foreign.sample.sample_id),
            ("catalog_projection_row", "smp_" + "f" * 32),
            ("experiment_revision", "expr_" + "f" * 32),
            ("term_revision", "trmv_" + "f" * 32),
            ("sample_revision", "smpr_" + "f" * 32),
        )
    ):
        poisoned = experiment | {
            "citations": experiment["citations"]
            | {"source_refs": [{"kind": kind, "ref_id": ref}]}
        }
        proposal_id = _register(
            catalog, project_id, poisoned, key=f"register-invalid-{number}"
        )["proposal_id"]
        with pytest.raises(CatalogError) as failure:
            _preview(catalog, proposal_id, "experiment")
        assert failure.value.code == "catalog_proposal_invalidated"

    with pytest.raises(CatalogError) as failure:
        _register(catalog, outside.project_id, sample, key="register-foreign-target")
    assert failure.value.code == "catalog_not_found"
    assert _proposal_receipts(catalog, project_id) == []


@pytest.mark.parametrize("catalog_stack", ["sqlite"], indirect=True)
def test_sqlite_adoption_rolls_back_an_already_inserted_experiment_revision(
    proposal_stack, monkeypatch
):
    _, catalog, project_id, experiment, sample = proposal_stack
    proposal_id = _register(catalog, project_id, experiment, sample)["proposal_id"]
    preview = _preview(catalog, proposal_id, "experiment", "sample")
    assert preview["selected_item_ids"] == ["experiment", "sample"]
    initial_sample = catalog.list_sample_revisions(sample_id=sample["target_id"])[
        "sample_revisions"
    ][0]
    # The experiment is applied first. The subsequent sample insert collides
    # with the existing primary key and must roll back both changes and receipts.
    with monkeypatch.context() as faults:
        faults.setattr(
            catalog,
            "_sample_revision_id_factory",
            lambda: initial_sample["sample_revision_id"],
        )
        with pytest.raises(CatalogError) as failure:
            _adopt(catalog, preview)
        assert failure.value.code == "catalog_revision_conflict"

    history = catalog.list_experiment_revisions(experiment_id=experiment["target_id"])
    assert [
        row["experiment_revision_id"] for row in history["experiment_revisions"]
    ] == [experiment["base_revision_id"]]
    revisions = catalog.list_sample_revisions(sample_id=sample["target_id"])
    assert [row["revision_number"] for row in revisions["sample_revisions"]] == [1, 2]
    assert _proposal_receipts(catalog, project_id, proposal_id) == []

    recovered = _adopt(catalog, preview)
    assert recovered["replayed"] is False
    assert _adopt(catalog, preview) == recovered | {"replayed": True}
    assert _proposal_receipts(catalog, project_id, proposal_id) == [recovered]


def test_preview_requires_an_administrator_and_the_registered_curation_version(
    proposal_stack,
):
    from encode_pipeline.platform.catalog_curation_preview import (
        M2_PROJECTION_VERSION,
        M2_VOCABULARY_VERSION,
    )

    _, catalog, project_id, experiment, _ = proposal_stack
    proposal_id = _register(catalog, project_id, experiment)["proposal_id"]
    unprivileged = "usr_" + "f" * 32
    with pytest.raises(CatalogError) as failure:
        _preview(catalog, proposal_id, "experiment", actor=unprivileged)
    assert (failure.value.status_code, failure.value.code) == (
        403,
        "ADMINISTRATOR_REQUIRED",
    )
    with pytest.raises(CatalogError) as failure:
        catalog.list_adoption_receipts(
            actor_user_id=unprivileged, project_id=project_id, proposal_id=proposal_id
        )
    assert failure.value.status_code == 403

    prior = _preview(catalog, proposal_id, "experiment")
    catalog.set_active_curation_versions(
        rule_version="rules/2026-11-19",
        vocabulary_version=M2_VOCABULARY_VERSION,
        projection_version=M2_PROJECTION_VERSION,
    )
    for stale in (
        lambda: _preview(catalog, proposal_id, "experiment"),
        lambda: _adopt(catalog, prior),
    ):
        with pytest.raises(CatalogError) as failure:
            stale()
        assert failure.value.code == "catalog_proposal_invalidated"
    assert _proposal_receipts(catalog, project_id, proposal_id) == []
    assert (
        catalog.list_experiment_revisions(experiment_id=experiment["target_id"])[
            "experiment_revisions"
        ][0]["batch"]
        == UNKNOWN
    )


def test_rebuild_preserves_pinned_experiment_revision_and_catalogued_sample(
    proposal_stack,
):
    _, catalog, _, experiment, sample = proposal_stack
    experiment_id = experiment["target_id"]
    sample_id = sample["target_id"]
    experiment_revision_id = experiment["base_revision_id"]
    sample_revision_id = sample["base_revision_id"]
    old_experiment = catalog.get_experiment_revision(
        experiment_id, experiment_revision_id
    )
    old_sample = catalog.get_sample_revision(sample_id, sample_revision_id)

    catalog._repository.rebuild_experiment_projection(experiment_revision_id)
    catalog._repository.rebuild_sample_projection(sample_revision_id)
    assert (
        catalog.get_experiment_revision(experiment_id, experiment_revision_id)
        == old_experiment
    )
    assert catalog.get_sample_revision(sample_id, sample_revision_id) == old_sample
    for rebuild, missing in (
        (catalog._repository.rebuild_experiment_projection, "expr_" + "f" * 32),
        (catalog._repository.rebuild_sample_projection, "smpr_" + "f" * 32),
    ):
        with pytest.raises(CatalogError) as failure:
            rebuild(missing)
        assert failure.value.code == "catalog_not_found"


def test_model_shortlist_and_pinned_revision_match_across_catalogs(catalog_stack):
    data, catalog = catalog_stack
    calls = []

    class InsufficientModel:
        def classify(self, **request):
            calls.append(request)
            return "insufficient_evidence"

    catalog._catalog_model = InsufficientModel()
    project = data.create_project("Vocabulary projection")
    foreign = data.create_project("Other project")
    sample = data.import_sample(
        project.project_id, stable_key="test", display_name="Test", attributes={}
    )
    experiment = _write(
        catalog.create_experiment,
        "experiment",
        project_id=project.project_id,
        stable_key="experiment",
        **_experiment(sample.revision.sample_revision_id),
    )
    targets = [
        _term(catalog, "target", target)
        for target in ("CTCF", "H2A.Z", "H3K4me3", "H3K27ac", "POLR2A")
    ]
    for target in targets:
        key = target["term_id"]
        _write(
            catalog.create_experiment,
            key,
            project_id=project.project_id,
            stable_key=key,
            **_experiment(
                sample.revision.sample_revision_id,
                target={"status": "known", **_pin(target)},
            ),
        )
    tissue = _term(catalog, "tissue", "Brain")
    pin = _pin(targets[1])["term_revision_id"]
    fields = {
        "project_id": project.project_id,
        "experiment_id": experiment["experiment_id"],
        "experiment_revision_id": experiment["experiment_revision_id"],
        "field": "target",
        "reason_code": "TARGET_TERM_MATCH",
        "clue": "H2A.Z",
        "vocabulary_term_revision_id": None,
    }
    for i, (overrides, code) in enumerate(
        (
            ({}, "catalog_payload_invalid"),
            ({"project_id": foreign.project_id}, "catalog_not_found"),
            (
                {"experiment_revision_id": "expr_" + "9" * 32},
                "catalog_revision_conflict",
            ),
            ({"field": "batch"}, "catalog_payload_invalid"),
            (
                {"vocabulary_term_revision_id": _pin(tissue)["term_revision_id"]},
                "catalog_payload_invalid",
            ),
        )
    ):
        with pytest.raises(CatalogError) as refused:
            _write(
                catalog.generate_classification_suggestion,
                f"refused-{i}",
                **(fields | overrides),
            )
        assert refused.value.code == code
    assert not calls
    result = _write(
        catalog.generate_classification_suggestion,
        "pinned",
        **(fields | {"vocabulary_term_revision_id": pin}),
    )
    assert result["outcome"] == "insufficient_evidence" and result["proposal"] is None
    assert result["experiment_revision_id"] == experiment["experiment_revision_id"]
    assert [call["candidates"] for call in calls] == [(("h2az", "H2A.Z"),)]
