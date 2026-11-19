"""Behavior shared by the in-memory and SQLite experiment catalogs."""

from datetime import datetime, timedelta, timezone
from hashlib import sha256
from itertools import count
import json

import pytest

from encode_pipeline.persistence import create_database_engine, create_session_factory
from encode_pipeline.persistence.data_registry import SqlAlchemyDataRegistryRepository
from encode_pipeline.persistence.experiment_catalog import (
    SqlAlchemyExperimentCatalogRepository,
)
from encode_pipeline.persistence.migrations import upgrade_database
from encode_pipeline.persistence.models import UserAccountRow
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
        repository = InMemoryExperimentCatalogRepository(registry)
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
