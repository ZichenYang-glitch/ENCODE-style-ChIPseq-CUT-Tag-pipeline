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
    "run_qc_metrics",
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
            "20260712_05"
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
    } == {
        "ck_run_execution_assignments_ack_requires_request",
        "ck_run_execution_assignments_claim_requires_dispatch",
        "ck_run_execution_assignments_request_reason_pair",
        "ck_run_execution_assignments_request_requires_claim",
    }
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
    qc_foreign_key = inspector.get_foreign_keys("run_qc_metrics")[0]
    assert qc_foreign_key["referred_table"] == "runs"
    assert qc_foreign_key["options"] == {"ondelete": "CASCADE"}
    assert {
        constraint["name"]
        for constraint in inspector.get_unique_constraints("run_qc_metrics")
    } == {"uq_run_qc_metrics_run_metric"}
    assert {
        constraint["name"]
        for constraint in inspector.get_check_constraints("run_qc_metrics")
    } == {
        "ck_run_qc_metrics_flag",
        "ck_run_qc_metrics_scope",
        "ck_run_qc_metrics_scope_identifiers",
        "ck_run_qc_metrics_value_text_length",
    }
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


def test_pr123_assignment_upgrades_with_empty_cancellation_evidence(tmp_path):
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
                    'legacy-running', 'workflow', '{}', 'running',
                    :timestamp, :timestamp, :timestamp, NULL, 'execution',
                    NULL, NULL, '{}'
                )
                """
            ),
            {"timestamp": "2026-07-12 00:00:00"},
        )
        connection.execute(
            text(
                """
                INSERT INTO run_execution_assignments (
                    run_id, job_id, backend, queue_name, created_at,
                    dispatched_at, claimed_at
                ) VALUES (
                    'legacy-running', 'legacy-job', 'rq', 'default',
                    :timestamp, :timestamp, :timestamp
                )
                """
            ),
            {"timestamp": "2026-07-12 00:00:00"},
        )
    engine.dispose()

    upgrade_database(database_url)
    upgraded = create_database_engine(database_url)
    with upgraded.connect() as connection:
        row = (
            connection.execute(
                text(
                    """
                SELECT job_id, dispatched_at, claimed_at,
                       cancellation_requested_at, cancellation_reason,
                       cancellation_acknowledged_at
                FROM run_execution_assignments
                WHERE run_id = 'legacy-running'
                """
                )
            )
            .mappings()
            .one()
        )
        assert row["job_id"] == "legacy-job"
        assert row["dispatched_at"] is not None
        assert row["claimed_at"] is not None
        assert row["cancellation_requested_at"] is None
        assert row["cancellation_reason"] is None
        assert row["cancellation_acknowledged_at"] is None
    upgraded.dispose()


def test_cancellation_intent_migration_downgrades_without_losing_assignment(tmp_path):
    database_url = f"sqlite:///{tmp_path / 'platform.db'}"
    upgrade_database(database_url)
    current = create_database_engine(database_url)
    with current.begin() as connection:
        connection.execute(
            text(
                """
                INSERT INTO runs (
                    run_id, workflow_id, inputs, status, created_at, updated_at,
                    started_at, ended_at, current_stage, cancellation_reason,
                    error, tags
                ) VALUES (
                    'downgrade-running', 'workflow', '{}', 'running',
                    :timestamp, :timestamp, :timestamp, NULL, 'execution',
                    NULL, NULL, '{}'
                )
                """
            ),
            {"timestamp": "2026-07-12 00:00:00"},
        )
        connection.execute(
            text(
                """
                INSERT INTO run_execution_assignments (
                    run_id, job_id, backend, queue_name, created_at,
                    dispatched_at, claimed_at, cancellation_requested_at,
                    cancellation_reason, cancellation_acknowledged_at
                ) VALUES (
                    'downgrade-running', 'downgrade-job', 'rq', 'default',
                    :timestamp, :timestamp, :timestamp, :timestamp,
                    'User requested cancellation.', :timestamp
                )
                """
            ),
            {"timestamp": "2026-07-12 00:00:00"},
        )
    current.dispose()

    downgrade_database(database_url, "20260712_03")
    engine = create_database_engine(database_url)
    inspector = inspect(engine)
    columns = {
        column["name"] for column in inspector.get_columns("run_execution_assignments")
    }
    assert "cancellation_requested_at" not in columns
    assert "cancellation_reason" not in columns
    assert "cancellation_acknowledged_at" not in columns
    with engine.connect() as connection:
        assert connection.scalar(text("SELECT version_num FROM alembic_version")) == (
            "20260712_03"
        )
        assignment = (
            connection.execute(
                text(
                    """
                    SELECT job_id, dispatched_at, claimed_at
                    FROM run_execution_assignments
                    WHERE run_id = 'downgrade-running'
                    """
                )
            )
            .mappings()
            .one()
        )
        assert assignment["job_id"] == "downgrade-job"
        assert assignment["dispatched_at"] is not None
        assert assignment["claimed_at"] is not None
    engine.dispose()


def test_qc_metric_migration_upgrades_current_main_without_changing_existing_rows(
    tmp_path,
):
    database_url = f"sqlite:///{tmp_path / 'platform.db'}"
    upgrade_database(database_url, "20260712_04")
    current = create_database_engine(database_url)
    with current.begin() as connection:
        connection.execute(
            text(
                """
                INSERT INTO runs (
                    run_id, workflow_id, inputs, status, created_at, updated_at,
                    started_at, ended_at, current_stage, cancellation_reason,
                    error, tags
                ) VALUES (
                    'existing-run', 'workflow', '{}', 'succeeded',
                    :timestamp, :timestamp, :timestamp, :timestamp,
                    'execution', NULL, NULL, '{}'
                )
                """
            ),
            {"timestamp": "2026-07-12 00:00:00"},
        )
        connection.execute(
            text(
                """
                INSERT INTO run_artifacts (
                    artifact_id, run_id, artifact_type, name, uri, mime_type,
                    produced_at, artifact_metadata
                ) VALUES (
                    'artifact-1', 'existing-run', 'file', 'summary.tsv',
                    'run://runs/existing-run/artifacts/artifact-1',
                    'text/tab-separated-values', :timestamp,
                    '{"output_type":"qc_summary"}'
                )
                """
            ),
            {"timestamp": "2026-07-12 00:00:00"},
        )
    current.dispose()

    upgrade_database(database_url)
    upgraded = create_database_engine(database_url)
    with upgraded.connect() as connection:
        assert connection.scalar(text("SELECT count(*) FROM runs")) == 1
        assert connection.scalar(text("SELECT count(*) FROM run_artifacts")) == 1
        assert connection.scalar(text("SELECT count(*) FROM run_qc_metrics")) == 0
        assert connection.scalar(text("SELECT version_num FROM alembic_version")) == (
            "20260712_05"
        )
    upgraded.dispose()


def test_qc_metric_migration_downgrades_independently(tmp_path):
    database_url = f"sqlite:///{tmp_path / 'platform.db'}"
    upgrade_database(database_url)

    downgrade_database(database_url, "20260712_04")
    downgraded = create_database_engine(database_url)
    inspector = inspect(downgraded)
    assert "run_qc_metrics" not in inspector.get_table_names()
    assert "run_artifacts" in inspector.get_table_names()
    with downgraded.connect() as connection:
        assert connection.scalar(text("SELECT version_num FROM alembic_version")) == (
            "20260712_04"
        )
    downgraded.dispose()
