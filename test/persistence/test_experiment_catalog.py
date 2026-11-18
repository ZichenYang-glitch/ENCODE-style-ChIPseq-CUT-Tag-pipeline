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
from encode_pipeline.persistence.models import ExperimentRevisionRow, UserAccountRow
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
