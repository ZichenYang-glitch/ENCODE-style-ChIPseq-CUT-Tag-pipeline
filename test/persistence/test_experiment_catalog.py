"""SQLite evidence for catalog migration, query plans, and sample trace."""

from __future__ import annotations

from datetime import datetime, timezone
from hashlib import sha256

import pytest

from sqlalchemy import text
from sqlalchemy.orm import Session

from encode_pipeline.persistence import (
    create_database_engine,
    create_session_factory,
)
from encode_pipeline.persistence.migrations import (
    downgrade_database,
    upgrade_database,
)
from encode_pipeline.persistence.experiment_catalog import (
    CatalogWrite,
    SqlAlchemyExperimentCatalogRepository,
)
from encode_pipeline.persistence.data_registry import (
    SqlAlchemyDataRegistryRepository,
    _project_row,
    _sample_revision_row,
    _sample_row,
)
from encode_pipeline.persistence.models import (
    ExperimentRevisionRow,
    SampleRevisionRow,
    UserAccountRow,
)
from encode_pipeline.platform.data_registry import (
    SAMPLE_REVISION_PAYLOAD_DIGEST_SCHEME,
    Project,
    ProjectKind,
    Sample,
    SampleRevision,
    build_sample_revision_payload_digest,
    canonical_sample_revision_payload,
)
from encode_pipeline.platform.experiment_catalog import (
    CatalogError,
    CatalogFilters,
    parse_page_filters,
)
from encode_pipeline.services.data_registry import DataRegistryService


NOW = datetime(2026, 11, 18, tzinfo=timezone.utc)
ACTOR = "usr_" + "c" * 32


def _filters(**kwargs) -> CatalogFilters:
    values = {
        "revision_view": None,
        "include_archived": False,
        "term_type": None,
        "created_from": None,
        "created_before": None,
        "catalog_state": None,
        "field_terms": {},
        "field_statuses": {},
    }
    values.update(kwargs)
    return parse_page_filters(**values)


def _write(key: str, body: bytes, operation: str) -> CatalogWrite:
    return CatalogWrite(
        actor_user_id=ACTOR,
        idempotency_key=key,
        request_sha256=sha256(body).hexdigest(),
        operation=operation,
        now=NOW,
    )


def test_catalog_migration_backfills_heads_and_downgrades(tmp_path) -> None:
    database_url = f"sqlite:///{tmp_path / 'catalog.db'}"
    upgrade_database(database_url)
    registry = DataRegistryService(
        repository=__import__(
            "encode_pipeline.persistence.data_registry",
            fromlist=["SqlAlchemyDataRegistryRepository"],
        ).SqlAlchemyDataRegistryRepository(
            create_session_factory(create_database_engine(database_url))
        )
    )
    project = registry.create_project("Existing")
    imported = registry.import_sample(
        project.project_id,
        stable_key="existing",
        display_name="Existing",
        attributes={"organism": "human"},
    )
    engine = create_database_engine(database_url)
    with engine.connect() as connection:
        heads = connection.execute(
            text("SELECT count(*) FROM sample_revision_heads")
        ).scalar()
        projections = connection.execute(
            text("SELECT count(*) FROM sample_catalog_projection")
        ).scalar()
    assert heads == 1
    assert projections == 1
    engine.dispose()

    downgrade_database(database_url, "20261110_16")
    downgraded = create_database_engine(database_url)
    names = set(__import__("sqlalchemy").inspect(downgraded).get_table_names())
    assert "experiments" not in names
    assert "sample_revision_heads" not in names
    with downgraded.connect() as connection:
        assert connection.scalar(text("SELECT count(*) FROM sample_revisions")) == 1
    downgraded.dispose()

    upgrade_database(database_url)
    restored = create_database_engine(database_url)
    with restored.connect() as connection:
        head = connection.execute(
            text(
                "SELECT sample_revision_id, revision_number "
                "FROM sample_revision_heads WHERE sample_id = :sample_id"
            ),
            {"sample_id": imported.sample.sample_id},
        ).one()
        projection = connection.scalar(
            text(
                "SELECT count(*) FROM sample_catalog_projection "
                "WHERE sample_revision_id = :revision_id"
            ),
            {"revision_id": imported.revision.sample_revision_id},
        )
    assert head.sample_revision_id == imported.revision.sample_revision_id
    assert head.revision_number == 1
    # Upgrade backfill restores the head pointer and does not invent a catalog.
    assert projection == 0
    restored.dispose()


def test_experiment_query_plan_uses_the_assay_index(tmp_path) -> None:
    database_url = f"sqlite:///{tmp_path / 'plans.db'}"
    upgrade_database(database_url)
    engine = create_database_engine(database_url)
    session_factory = create_session_factory(engine)
    with Session(engine) as session:
        session.add(
            UserAccountRow(
                user_id=ACTOR,
                username="catalog-admin",
                role="administrator",
                status="enabled",
                password_hash="hash",
                created_at=NOW,
                updated_at=NOW,
                password_changed_at=NOW,
            )
        )
        session.commit()
    registry = DataRegistryService(
        repository=__import__(
            "encode_pipeline.persistence.data_registry",
            fromlist=["SqlAlchemyDataRegistryRepository"],
        ).SqlAlchemyDataRegistryRepository(session_factory)
    )
    catalog = SqlAlchemyExperimentCatalogRepository(session_factory)
    project = registry.create_project("Plan")
    sample = registry.import_sample(
        project.project_id,
        stable_key="sample",
        display_name="Sample",
        attributes={},
    )
    assay = catalog.create_term(
        _write(b"assay", b"assay", "catalogCreateTerm"),
        term_id="term_" + "a" * 32,
        term_revision_id="trmv_" + "a" * 32,
        term_type="assay",
        display_name="CUT&Tag",
        namespace=None,
        accession=None,
        aliases=(),
    )
    target = catalog.create_term(
        _write(b"target", b"target", "catalogCreateTerm"),
        term_id="term_" + "b" * 32,
        term_revision_id="trmv_" + "b" * 32,
        term_type="target",
        display_name="H3K27ac",
        namespace=None,
        accession=None,
        aliases=(),
    )
    batch = catalog.create_term(
        _write(b"batch", b"batch", "catalogCreateTerm"),
        term_id="term_" + "c" * 32,
        term_revision_id="trmv_" + "c" * 32,
        term_type="batch",
        display_name="batch-1",
        namespace=None,
        accession=None,
        aliases=(),
    )
    role = catalog.create_term(
        _write(b"role", b"role", "catalogCreateTerm"),
        term_id="term_" + "d" * 32,
        term_revision_id="trmv_" + "d" * 32,
        term_type="replicate_role",
        display_name="biological replicate",
        namespace=None,
        accession=None,
        aliases=(),
    )

    def known(term: dict) -> dict:
        return {
            "status": "known",
            "term_id": term["term_id"],
            "term_revision_id": term["revision"]["term_revision_id"],
        }

    for index in range(200):
        catalog.create_experiment(
            _write(
                f"exp-{index}".encode(),
                f"exp-{index}".encode(),
                "catalogCreateExperiment",
            ),
            project_id=project.project_id,
            experiment_id=f"exp_{index:032x}",
            experiment_revision_id=f"expr_{index:032x}",
            stable_key=f"plan-{index}",
            display_name=f"Plan {index}",
            assay=known(assay),
            target=known(target),
            batch=known(batch),
            samples=[
                {
                    "sample_revision_id": sample.revision.sample_revision_id,
                    "replicate_role": known(role),
                }
            ],
        )
    filters = _filters(
        kind="experiments",
        project_id=project.project_id,
        field_terms={"assay": [assay["term_id"]]},
    )
    plan = catalog.explain_experiment_query(filters)
    assert "ix_experiment_catalog_projection_assay" in plan
    assert "SCAN sample_revisions" not in plan
    assert "SCAN experiment_revisions" not in plan
    revision_id = f"expr_{0:032x}"
    with Session(engine) as db:
        before = db.get(ExperimentRevisionRow, revision_id).canonical_payload
    catalog.rebuild_experiment_projection(revision_id)
    with Session(engine) as db:
        after = db.get(ExperimentRevisionRow, revision_id).canonical_payload
    assert after == before
    with Session(engine) as db:
        account = db.get(UserAccountRow, ACTOR)
        assert account is not None
        account.role = "member"
        db.commit()
    with pytest.raises(CatalogError) as caught:
        catalog.create_term(
            _write("member-write", b"member-write", "catalogCreateTerm"),
            term_id="term_" + "e" * 32,
            term_revision_id="trmv_" + "e" * 32,
            term_type="tag",
            display_name="member tag",
            namespace=None,
            accession=None,
            aliases=(),
        )
    assert caught.value.status_code == 403
    assert caught.value.code == "ADMINISTRATOR_REQUIRED"
    engine.dispose()


def _parent_revision(display_name: str, attributes: dict, **identity) -> SampleRevision:
    canonical = canonical_sample_revision_payload(
        display_name=display_name,
        attributes=attributes,
    )
    return SampleRevision(
        canonical_payload=canonical,
        payload_digest_scheme=SAMPLE_REVISION_PAYLOAD_DIGEST_SCHEME,
        payload_digest=build_sample_revision_payload_digest(canonical),
        created_at=NOW,
        revision_number=1,
        **identity,
    )


def test_parent_database_keeps_free_catalog_attributes(tmp_path) -> None:
    database_url = f"sqlite:///{tmp_path / 'parent.db'}"
    upgrade_database(database_url, "20261110_16")
    project_id = "prj_" + "1" * 32
    lookalike = {
        "schema_version": 1,
        "species": {"status": "unknown"},
        "tissue": {"status": "unknown"},
        "cell_type": {
            "status": "known",
            "values": [
                {
                    "term_id": "term_" + "e" * 32,
                    "term_revision_id": "trmv_" + "e" * 32,
                }
            ],
        },
        "cell_line": {"status": "not_applicable", "values": []},
        "treatment": {"status": "not_applicable", "values": []},
        "tag": {"status": "not_applicable", "values": []},
    }
    records = {
        "notebook": ("paper notebook 17", "smp_" + "2" * 32, "smpr_" + "2" * 32),
        "lookalike": (lookalike, "smp_" + "3" * 32, "smpr_" + "3" * 32),
    }
    engine = create_database_engine(database_url)
    with Session(engine) as session:
        session.add(
            UserAccountRow(
                user_id=ACTOR,
                username="parent-admin",
                role="administrator",
                status="enabled",
                password_hash="hash",
                created_at=NOW,
                updated_at=NOW,
                password_changed_at=NOW,
            )
        )
        session.add(
            _project_row(
                Project(
                    project_id=project_id,
                    display_name="Parent project",
                    kind=ProjectKind.USER,
                    created_at=NOW,
                )
            )
        )
        session.flush()
        for name, (catalog_value, sample_id, revision_id) in records.items():
            session.add(
                _sample_row(
                    Sample(
                        sample_id=sample_id,
                        project_id=project_id,
                        stable_key=name,
                        created_at=NOW,
                    )
                )
            )
            session.flush()
            session.add(
                _sample_revision_row(
                    _parent_revision(
                        name,
                        {"organism": "human", "catalog": catalog_value},
                        sample_revision_id=revision_id,
                        project_id=project_id,
                        sample_id=sample_id,
                    )
                )
            )
        session.commit()
    engine.dispose()

    upgrade_database(database_url)
    engine = create_database_engine(database_url)
    factory = create_session_factory(engine)
    registry = DataRegistryService(
        repository=SqlAlchemyDataRegistryRepository(factory),
        sample_revision_id_factory=iter(
            ["smpr_" + "4" * 32, "smpr_" + "5" * 32, "smpr_" + "6" * 32]
        ).__next__,
    )
    catalog = SqlAlchemyExperimentCatalogRepository(factory)
    listed = catalog.list_samples(
        filters=_filters(kind="samples", project_id=project_id),
        cursor=None,
        limit=20,
    )
    by_sample = {item["sample_id"]: item for item in listed["samples"]}
    species = catalog.create_term(
        _write("species", b"species", "catalogCreateTerm"),
        term_id="term_" + "a" * 32,
        term_revision_id="trmv_" + "a" * 32,
        term_type="species",
        display_name="human",
        namespace=None,
        accession=None,
        aliases=(),
    )
    tissue = catalog.create_term(
        _write("tissue", b"tissue", "catalogCreateTerm"),
        term_id="term_" + "b" * 32,
        term_revision_id="trmv_" + "b" * 32,
        term_type="tissue",
        display_name="blood",
        namespace=None,
        accession=None,
        aliases=(),
    )
    cell = catalog.create_term(
        _write("cell", b"cell", "catalogCreateTerm"),
        term_id="term_" + "c" * 32,
        term_revision_id="trmv_" + "c" * 32,
        term_type="cell_type",
        display_name="regulatory T cell",
        namespace=None,
        accession=None,
        aliases=(),
    )
    for index, (name, (catalog_value, sample_id, revision_id)) in enumerate(
        records.items()
    ):
        with Session(engine) as session:
            original = session.get(SampleRevisionRow, revision_id)
            before = original.canonical_payload
            before_digest = original.payload_digest
        assert by_sample[sample_id]["catalog_state"] == "absent"
        assert "cell_type" not in by_sample[sample_id]
        historical = catalog.get_sample_revision(sample_id, revision_id)
        assert historical["catalog_state"] == "absent"
        renamed = registry.revise_sample(
            sample_id,
            display_name=f"{name} renamed",
            attributes={"organism": "human"},
        )
        assert renamed.revision_number == 2
        catalog.rebuild_sample_projection(revision_id)
        with Session(engine) as session:
            stored = session.get(SampleRevisionRow, revision_id)
            assert stored.canonical_payload == before
            assert stored.payload_digest == before_digest
            assert "paper notebook 17" in before or "term_" + "e" * 32 in before
        assert catalog.get_sample_revision(sample_id, revision_id)["catalog_state"] == (
            "absent"
        )
        classified = catalog.revise_sample_catalog(
            _write(f"classify-{index}", f"classify-{index}".encode(), "sample"),
            sample_id=sample_id,
            sample_revision_id="smpr_" + str(7 + index) * 32,
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
        assert catalog.get_sample_revision(sample_id, revision_id)["catalog_state"] == (
            "absent"
        )
        with Session(engine) as session:
            stored = session.get(SampleRevisionRow, revision_id)
            assert stored.canonical_payload == before
            assert stored.payload_digest == before_digest
        del catalog_value
    engine.dispose()


def test_receipt_target_rejects_a_different_object(tmp_path) -> None:
    database_url = f"sqlite:///{tmp_path / 'receipts.db'}"
    upgrade_database(database_url)
    engine = create_database_engine(database_url)
    factory = create_session_factory(engine)
    with Session(engine) as session:
        session.add(
            UserAccountRow(
                user_id=ACTOR,
                username="receipt-admin",
                role="administrator",
                status="enabled",
                password_hash="hash",
                created_at=NOW,
                updated_at=NOW,
                password_changed_at=NOW,
            )
        )
        session.commit()
    registry = DataRegistryService(repository=SqlAlchemyDataRegistryRepository(factory))
    catalog = SqlAlchemyExperimentCatalogRepository(factory)
    left = registry.create_project("Left")
    right = registry.create_project("Right")
    sample = registry.import_sample(
        left.project_id,
        stable_key="shared",
        display_name="Shared",
        attributes={},
    )

    def make_term(label: str, term_type: str, name: str) -> dict:
        return catalog.create_term(
            _write(label, label.encode(), "catalogCreateTerm"),
            term_id="term_" + label.encode().hex()[:32].ljust(32, "0"),
            term_revision_id="trmv_" + label.encode().hex()[:32].ljust(32, "1"),
            term_type=term_type,
            display_name=name,
            namespace=None,
            accession=None,
            aliases=(),
        )

    alpha = make_term("alpha", "tag", "alpha")
    beta = make_term("beta", "tag", "beta")
    archive = b'{"expected_revision_number":1}'
    archived = catalog.archive_term(
        _write("archive-once", archive, "catalogArchiveTerm"),
        term_id=alpha["term_id"],
        expected_revision_number=1,
    )
    with pytest.raises(CatalogError) as conflict:
        catalog.archive_term(
            _write("archive-once", archive, "catalogArchiveTerm"),
            term_id=beta["term_id"],
            expected_revision_number=1,
        )
    assert conflict.value.status_code == 409
    assert conflict.value.code == "catalog_idempotency_conflict"
    assert catalog.get_term(beta["term_id"])["archived_at"] is None
    replay = catalog.archive_term(
        _write("archive-once", archive, "catalogArchiveTerm"),
        term_id=alpha["term_id"],
        expected_revision_number=1,
    )
    assert replay["term_id"] == archived["term_id"]
    with engine.connect() as connection:
        assert (
            connection.scalar(
                text(
                    "SELECT count(*) FROM catalog_write_receipts "
                    "WHERE idempotency_key = 'archive-once'"
                )
            )
            == 1
        )
    assay = make_term("assay", "assay", "CUT&Tag")
    target = make_term("target", "target", "H3K27ac")
    batch = make_term("batch", "batch", "batch-1")
    role = make_term("role", "replicate_role", "biological replicate")

    def known(term: dict) -> dict:
        return {
            "status": "known",
            "term_id": term["term_id"],
            "term_revision_id": term["revision"]["term_revision_id"],
        }

    body = b'{"experiment":"shared"}'
    created = catalog.create_experiment(
        _write("create-once", body, "catalogCreateExperiment"),
        project_id=left.project_id,
        experiment_id="exp_" + "8" * 32,
        experiment_revision_id="expr_" + "8" * 32,
        stable_key="shared-exp",
        display_name="Shared",
        assay=known(assay),
        target=known(target),
        batch=known(batch),
        samples=[
            {
                "sample_revision_id": sample.revision.sample_revision_id,
                "replicate_role": known(role),
            }
        ],
    )
    with pytest.raises(CatalogError) as project_conflict:
        catalog.create_experiment(
            _write("create-once", body, "catalogCreateExperiment"),
            project_id=right.project_id,
            experiment_id="exp_" + "9" * 32,
            experiment_revision_id="expr_" + "9" * 32,
            stable_key="shared-exp",
            display_name="Shared",
            assay=known(assay),
            target=known(target),
            batch=known(batch),
            samples=[
                {
                    "sample_revision_id": sample.revision.sample_revision_id,
                    "replicate_role": known(role),
                }
            ],
        )
    assert project_conflict.value.code == "catalog_idempotency_conflict"
    empty = catalog.list_experiments(
        filters=_filters(kind="experiments", project_id=right.project_id),
        cursor=None,
        limit=20,
    )
    assert empty["ok"] is True
    assert empty["experiments"] == []
    with engine.connect() as connection:
        assert (
            connection.scalar(
                text("SELECT count(*) FROM experiments WHERE project_id = :project"),
                {"project": right.project_id},
            )
            == 0
        )
    again = catalog.create_experiment(
        _write("create-once", body, "catalogCreateExperiment"),
        project_id=left.project_id,
        experiment_id="exp_" + "a" * 32,
        experiment_revision_id="expr_" + "a" * 32,
        stable_key="shared-exp",
        display_name="Shared",
        assay=known(assay),
        target=known(target),
        batch=known(batch),
        samples=[
            {
                "sample_revision_id": sample.revision.sample_revision_id,
                "replicate_role": known(role),
            }
        ],
    )
    assert again["experiment_id"] == created["experiment_id"]
    archive_body = b'{"expected_revision_number":1}'
    archived = catalog.archive_experiment(
        _write("archive-experiment", archive_body, "catalogArchiveExperiment"),
        experiment_id=created["experiment_id"],
        expected_revision_number=1,
    )
    assert archived["experiment_id"] == created["experiment_id"]
    assert archived["revision_number"] == 1
    with pytest.raises(CatalogError) as archived_again:
        catalog.archive_experiment(
            _write(
                "archive-experiment-again", archive_body, "catalogArchiveExperiment"
            ),
            experiment_id=created["experiment_id"],
            expected_revision_number=1,
        )
    assert archived_again.value.code == "catalog_already_archived"
    engine.dispose()
