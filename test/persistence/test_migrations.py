"""Tests for the initial Alembic platform schema."""

from __future__ import annotations

from sqlalchemy import inspect, text

from encode_pipeline.persistence import create_database_engine
from encode_pipeline.persistence.migrations import (
    downgrade_database,
    upgrade_database,
)


EXPECTED_TABLES = {
    "alembic_version",
    "run_artifacts",
    "run_events",
    "run_execution_assignments",
    "run_logs",
    "run_workflow_build_identities",
    "runs",
}


def test_initial_migration_creates_versioned_run_schema(tmp_path):
    database_url = f"sqlite:///{tmp_path / 'platform.db'}"

    upgrade_database(database_url)
    engine = create_database_engine(database_url)
    inspector = inspect(engine)

    assert set(inspector.get_table_names()) == EXPECTED_TABLES
    with engine.connect() as connection:
        assert connection.scalar(text("SELECT version_num FROM alembic_version")) == (
            "20260712_03"
        )
        assert connection.scalar(text("PRAGMA foreign_keys")) == 1
        assert connection.scalar(text("PRAGMA journal_mode")) == "wal"

    event_constraints = {
        constraint["name"]
        for constraint in inspector.get_unique_constraints("run_events")
    }
    log_constraints = {
        constraint["name"]
        for constraint in inspector.get_unique_constraints("run_logs")
    }
    assignment_constraints = {
        constraint["name"]
        for constraint in inspector.get_unique_constraints("run_execution_assignments")
    }
    assert "uq_run_events_run_sequence" in event_constraints
    assert "uq_run_logs_run_stream_sequence" in log_constraints
    assert "uq_run_execution_assignments_job_id" in assignment_constraints
    assert {
        constraint["name"]
        for constraint in inspector.get_check_constraints("run_execution_assignments")
    } == {"ck_run_execution_assignments_claim_requires_dispatch"}
    assert inspector.get_foreign_keys("run_events")[0]["options"] == {
        "ondelete": "CASCADE"
    }
    assignment_foreign_key = inspector.get_foreign_keys("run_execution_assignments")[0]
    assert assignment_foreign_key["referred_table"] == "runs"
    assert assignment_foreign_key["options"] == {"ondelete": "CASCADE"}
    assert inspector.get_pk_constraint("run_execution_assignments")[
        "constrained_columns"
    ] == ["run_id"]
    build_foreign_key = inspector.get_foreign_keys("run_workflow_build_identities")[0]
    assert build_foreign_key["referred_table"] == "runs"
    assert build_foreign_key["options"] == {"ondelete": "CASCADE"}
    assert inspector.get_pk_constraint("run_workflow_build_identities")[
        "constrained_columns"
    ] == ["run_id"]
    engine.dispose()


def test_initial_migration_can_downgrade_and_reapply(tmp_path):
    database_url = f"sqlite:///{tmp_path / 'platform.db'}"
    upgrade_database(database_url)

    downgrade_database(database_url, "base")
    downgraded_engine = create_database_engine(database_url)
    assert inspect(downgraded_engine).get_table_names() == ["alembic_version"]
    downgraded_engine.dispose()

    upgrade_database(database_url)
    upgraded_engine = create_database_engine(database_url)
    assert set(inspect(upgraded_engine).get_table_names()) == EXPECTED_TABLES
    upgraded_engine.dispose()


def test_execution_assignment_migration_can_downgrade_independently(tmp_path):
    database_url = f"sqlite:///{tmp_path / 'platform.db'}"
    upgrade_database(database_url)

    downgrade_database(database_url, "20260711_01")
    downgraded_engine = create_database_engine(database_url)
    assert (
        "run_execution_assignments" not in inspect(downgraded_engine).get_table_names()
    )
    with downgraded_engine.connect() as connection:
        assert connection.scalar(text("SELECT version_num FROM alembic_version")) == (
            "20260711_01"
        )
    downgraded_engine.dispose()

    upgrade_database(database_url)
    upgraded_engine = create_database_engine(database_url)
    assert "run_execution_assignments" in inspect(upgraded_engine).get_table_names()
    upgraded_engine.dispose()


def test_build_identity_migration_does_not_backfill_legacy_planned_runs(tmp_path):
    database_url = f"sqlite:///{tmp_path / 'platform.db'}"
    upgrade_database(database_url, "20260711_02")
    engine = create_database_engine(database_url)
    with engine.begin() as connection:
        connection.execute(
            text(
                """
                INSERT INTO runs (
                    run_id, workflow_id, inputs, status, created_at, updated_at,
                    started_at, ended_at, current_stage, cancellation_reason,
                    error, tags
                ) VALUES (
                    :run_id, :workflow_id, :inputs, 'planned', :created_at,
                    :updated_at, NULL, NULL, 'preflight', NULL, NULL, :tags
                )
                """
            ),
            {
                "run_id": "legacy-planned",
                "workflow_id": "encode-style-chipseq-cuttag-atac-mnase",
                "inputs": '{"config":{},"samples":null,"options":{}}',
                "created_at": "2026-07-12 00:00:00",
                "updated_at": "2026-07-12 00:00:00",
                "tags": "{}",
            },
        )
    engine.dispose()

    upgrade_database(database_url)
    upgraded = create_database_engine(database_url)
    with upgraded.connect() as connection:
        assert (
            connection.scalar(
                text("SELECT status FROM runs WHERE run_id = 'legacy-planned'")
            )
            == "planned"
        )
        assert (
            connection.scalar(
                text("SELECT count(*) FROM run_workflow_build_identities")
            )
            == 0
        )
    upgraded.dispose()
