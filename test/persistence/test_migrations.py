"""Tests for the initial Alembic platform schema."""

from __future__ import annotations

import json

import pytest
from sqlalchemy import inspect, text

from encode_pipeline.persistence import (
    SqlAlchemyRunRepository,
    create_database_engine,
    create_session_factory,
)
from encode_pipeline.persistence.migrations import (
    downgrade_database,
    upgrade_database,
)


pytestmark = pytest.mark.full_main


EXPECTED_TABLES = {
    "alembic_version",
    "run_artifacts",
    "run_events",
    "run_execution_assignments",
    "run_logs",
    "run_qc_metrics",
    "run_result_attempts",
    "run_result_states",
    "run_workflow_build_identities",
    "runs",
    "validated_input_snapshots",
}


def test_initial_migration_creates_versioned_run_schema(tmp_path):
    database_url = f"sqlite:///{tmp_path / 'platform.db'}"

    upgrade_database(database_url)
    engine = create_database_engine(database_url)
    inspector = inspect(engine)

    assert set(inspector.get_table_names()) == EXPECTED_TABLES
    with engine.connect() as connection:
        assert connection.scalar(text("SELECT version_num FROM alembic_version")) == (
            "20260717_08"
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
    snapshot_foreign_key = inspector.get_foreign_keys("validated_input_snapshots")[0]
    assert snapshot_foreign_key["referred_table"] == "runs"
    assert snapshot_foreign_key["options"] == {"ondelete": "RESTRICT"}
    assert {
        constraint["name"]
        for constraint in inspector.get_unique_constraints("validated_input_snapshots")
    } == {"uq_validated_input_snapshots_consumed_run_id"}
    assert {
        constraint["name"]
        for constraint in inspector.get_check_constraints("validated_input_snapshots")
    } == {
        "ck_validated_input_snapshots_consumption_pair",
        "ck_validated_input_snapshots_digest_length",
        "ck_validated_input_snapshots_expiry",
        "ck_validated_input_snapshots_success",
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
    assert inspector.get_pk_constraint("run_result_states")["constrained_columns"] == [
        "run_id"
    ]
    assert inspector.get_foreign_keys("run_result_states")[0]["options"] == {
        "ondelete": "CASCADE"
    }
    assert {
        constraint["name"]
        for constraint in inspector.get_check_constraints("run_result_states")
    } == {
        "ck_run_result_states_artifact_attempt",
        "ck_run_result_states_artifact_binding",
        "ck_run_result_states_artifact_outcome",
        "ck_run_result_states_artifact_reason",
        "ck_run_result_states_nonnegative_revisions",
        "ck_run_result_states_qc_attempt",
        "ck_run_result_states_qc_binding",
        "ck_run_result_states_qc_outcome",
        "ck_run_result_states_qc_reason",
    }
    assert inspector.get_pk_constraint("run_result_attempts")[
        "constrained_columns"
    ] == ["attempt_id"]
    assert inspector.get_foreign_keys("run_result_attempts")[0]["options"] == {
        "ondelete": "CASCADE"
    }
    assert {
        constraint["name"]
        for constraint in inspector.get_check_constraints("run_result_attempts")
    } == {"ck_run_result_attempts_binding"}
    assert {
        index["name"] for index in inspector.get_indexes("run_result_attempts")
    } == {"ix_run_result_attempts_run_id"}
    assert {column["name"] for column in inspector.get_columns("run_artifacts")} >= {
        "revision"
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
        assert (
            connection.scalar(text("SELECT count(*) FROM validated_input_snapshots"))
            == 0
        )
        assert connection.scalar(text("SELECT version_num FROM alembic_version")) == (
            "20260717_08"
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


def test_validated_snapshot_migration_upgrades_current_main_without_changing_runs(
    tmp_path,
):
    database_url = f"sqlite:///{tmp_path / 'platform.db'}"
    upgrade_database(database_url, "20260712_05")
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
                    'existing-run', 'workflow', '{}', 'created',
                    :timestamp, :timestamp, NULL, NULL, NULL, NULL, NULL, '{}'
                )
                """
            ),
            {"timestamp": "2026-07-14 00:00:00"},
        )
    current.dispose()

    upgrade_database(database_url)
    upgraded = create_database_engine(database_url)
    with upgraded.connect() as connection:
        assert connection.scalar(text("SELECT count(*) FROM runs")) == 1
        assert (
            connection.scalar(text("SELECT count(*) FROM validated_input_snapshots"))
            == 0
        )
        assert connection.scalar(text("SELECT version_num FROM alembic_version")) == (
            "20260717_08"
        )
    upgraded.dispose()


def test_run_history_index_migration_preserves_rows_and_supports_all_query_shapes(
    tmp_path,
) -> None:
    database_url = f"sqlite:///{tmp_path / 'run-history-upgrade.db'}"
    upgrade_database(database_url, "20260714_06")
    engine = create_database_engine(database_url)
    with engine.begin() as connection:
        connection.execute(
            text(
                "INSERT INTO runs "
                "(run_id, workflow_id, inputs, status, created_at, updated_at, "
                "started_at, ended_at, current_stage, cancellation_reason, "
                "error, tags) VALUES "
                "(:run_id, :workflow_id, '{}', :status, :created_at, :updated_at, "
                "NULL, NULL, NULL, NULL, NULL, '{}')"
            ),
            {
                "run_id": "run-1",
                "workflow_id": "workflow-a",
                "status": "created",
                "created_at": "2026-07-14 08:00:00",
                "updated_at": "2026-07-14 08:00:00",
            },
        )
    engine.dispose()

    upgrade_database(database_url)
    upgraded = create_database_engine(database_url)
    expected_indexes = {
        "ix_runs_created_run_id",
        "ix_runs_workflow_created_run_id",
        "ix_runs_status_created_run_id",
        "ix_runs_workflow_status_created_run_id",
    }
    assert expected_indexes <= {
        index["name"] for index in inspect(upgraded).get_indexes("runs")
    }
    with upgraded.connect() as connection:
        assert connection.scalar(text("SELECT count(*) FROM runs")) == 1
        assert connection.scalar(text("SELECT version_num FROM alembic_version")) == (
            "20260717_08"
        )
        plans = {
            "ix_runs_created_run_id": (
                "SELECT run_id FROM runs "
                "ORDER BY created_at DESC, run_id DESC LIMIT 51",
                {},
            ),
            "ix_runs_workflow_created_run_id": (
                "SELECT run_id FROM runs WHERE workflow_id = :workflow_id "
                "ORDER BY created_at DESC, run_id DESC LIMIT 51",
                {"workflow_id": "workflow-a"},
            ),
            "ix_runs_status_created_run_id": (
                "SELECT run_id FROM runs WHERE status = :status "
                "ORDER BY created_at DESC, run_id DESC LIMIT 51",
                {"status": "created"},
            ),
            "ix_runs_workflow_status_created_run_id": (
                "SELECT run_id FROM runs "
                "WHERE workflow_id = :workflow_id AND status = :status "
                "ORDER BY created_at DESC, run_id DESC LIMIT 51",
                {"workflow_id": "workflow-a", "status": "created"},
            ),
        }
        for index_name, (statement, parameters) in plans.items():
            detail = " ".join(
                str(row[3])
                for row in connection.execute(
                    text(f"EXPLAIN QUERY PLAN {statement}"), parameters
                )
            )
            assert index_name in detail
            assert "USE TEMP B-TREE FOR ORDER BY" not in detail
    upgraded.dispose()

    downgrade_database(database_url, "20260714_06")
    downgraded = create_database_engine(database_url)
    assert expected_indexes.isdisjoint(
        {index["name"] for index in inspect(downgraded).get_indexes("runs")}
    )
    with downgraded.connect() as connection:
        assert connection.scalar(text("SELECT count(*) FROM runs")) == 1
    downgraded.dispose()

    upgrade_database(database_url)


def test_result_generation_migration_keeps_legacy_results_explicitly_unbound(
    tmp_path,
) -> None:
    database_url = f"sqlite:///{tmp_path / 'legacy-results.db'}"
    upgrade_database(database_url, "20260714_07")
    engine = create_database_engine(database_url)
    with engine.begin() as connection:
        connection.execute(
            text(
                "INSERT INTO runs "
                "(run_id, workflow_id, inputs, status, created_at, updated_at, "
                "started_at, ended_at, current_stage, cancellation_reason, error, tags) "
                "VALUES ('legacy-run', 'workflow', '{}', 'succeeded', :now, :now, "
                ":now, :now, 'execution', NULL, NULL, '{}')"
            ),
            {"now": "2026-07-17 00:00:00"},
        )
        connection.execute(
            text(
                "INSERT INTO run_artifacts "
                "(artifact_id, run_id, artifact_type, name, uri, mime_type, "
                "produced_at, artifact_metadata) VALUES "
                "('artifact-1', 'legacy-run', 'file', 'summary.tsv', "
                "'run://runs/legacy-run/artifacts/artifact-1', 'text/plain', :now, "
                '\'{"output_type":"qc_summary"}\')'
            ),
            {"now": "2026-07-17 00:00:00"},
        )
        connection.execute(
            text(
                "INSERT INTO run_qc_metrics "
                "(metric_id, run_id, metric_key, display_name, value_text, unit, "
                "scope, sample_id, experiment_id, assay, qc_flag, "
                "source_artifact_id, produced_at) VALUES "
                "(:metric_id, 'legacy-run', 'sequencing.total_reads', "
                "'Total reads', '10', 'count', 'sample', 'S1', NULL, 'rnaseq', "
                "NULL, 'artifact-1', :now)"
            ),
            {
                "metric_id": "qcmetric-" + "a" * 64,
                "now": "2026-07-17 00:00:00",
            },
        )
    engine.dispose()

    upgrade_database(database_url)
    upgraded = create_database_engine(database_url)
    with upgraded.connect() as connection:
        state = (
            connection.execute(
                text(
                    "SELECT artifact_revision, artifact_generation, qc_revision, "
                    "qc_generation FROM run_result_states WHERE run_id='legacy-run'"
                )
            )
            .mappings()
            .one()
        )
        assert dict(state) == {
            "artifact_revision": 0,
            "artifact_generation": None,
            "qc_revision": 0,
            "qc_generation": None,
        }
        assert connection.scalar(text("SELECT count(*) FROM run_result_attempts")) == 0
        assert (
            connection.scalar(
                text(
                    "SELECT revision FROM run_artifacts "
                    "WHERE run_id='legacy-run' AND artifact_id='artifact-1'"
                )
            )
            is None
        )

    repository = SqlAlchemyRunRepository(create_session_factory(upgraded))
    with pytest.raises(ValueError, match="generation is unbound"):
        repository.list_artifacts("legacy-run")
    with pytest.raises(ValueError, match="generation is unbound"):
        repository.list_qc_metrics("legacy-run")
    upgraded.dispose()


def test_v030_supported_prior_schema_upgrade_preserves_complete_product_record(
    tmp_path,
) -> None:
    """Revision 07 is the supported pre-generation schema fixture."""
    database_url = f"sqlite:///{tmp_path / 'v030-upgrade.db'}"
    upgrade_database(database_url, "20260714_07")
    engine = create_database_engine(database_url)
    now = "2026-07-17 12:00:00"
    run_id = "release-upgrade-run"
    workflow_id = "bulk-rnaseq"
    artifact_id = "release-artifact"
    metric_id = "qcmetric-" + "b" * 64
    snapshot_id = "vsnap_" + "c" * 32
    digest = "d" * 64
    with engine.begin() as connection:
        connection.execute(
            text(
                """
                INSERT INTO runs (
                    run_id, workflow_id, inputs, status, created_at, updated_at,
                    started_at, ended_at, current_stage, cancellation_reason,
                    error, tags
                ) VALUES (
                    :run_id, :workflow_id, :inputs, 'succeeded', :now, :now,
                    :now, :now, 'results', NULL, NULL, :tags
                )
                """
            ),
            {
                "run_id": run_id,
                "workflow_id": workflow_id,
                "inputs": '{"config":{"standard":{}},"samples":[],"options":{}}',
                "now": now,
                "tags": '{"release_fixture":"v0.3.0"}',
            },
        )
        connection.execute(
            text(
                """
                INSERT INTO run_events (
                    event_id, run_id, sequence, event_type, timestamp, status,
                    stage, message, context, issue
                ) VALUES (
                    'event-1', :run_id, 1, 'status_changed', :now, 'succeeded',
                    'results', 'Results indexed.', :context, NULL
                )
                """
            ),
            {"run_id": run_id, "now": now, "context": '{"attempt_id":"attempt-1"}'},
        )
        connection.execute(
            text(
                """
                INSERT INTO run_logs (
                    chunk_id, run_id, stream_name, sequence, timestamp, lines
                ) VALUES (
                    'chunk-1', :run_id, 'stdout', 1, :now, :lines
                )
                """
            ),
            {"run_id": run_id, "now": now, "lines": '["workflow complete"]'},
        )
        connection.execute(
            text(
                """
                INSERT INTO run_execution_assignments (
                    run_id, job_id, backend, queue_name, created_at,
                    dispatched_at, claimed_at, cancellation_requested_at,
                    cancellation_reason, cancellation_acknowledged_at
                ) VALUES (
                    :run_id, 'release-job', 'rq', 'release-queue', :now,
                    :now, :now, NULL, NULL, NULL
                )
                """
            ),
            {"run_id": run_id, "now": now},
        )
        connection.execute(
            text(
                """
                INSERT INTO run_workflow_build_identities (
                    run_id, workflow_id, adapter_version, scheme,
                    logical_entrypoint, digest, captured_at
                ) VALUES (
                    :run_id, :workflow_id, '0.3.0', 'sha256-tree-v1',
                    'main.nf', :digest, :now
                )
                """
            ),
            {
                "run_id": run_id,
                "workflow_id": workflow_id,
                "digest": digest,
                "now": now,
            },
        )
        connection.execute(
            text(
                """
                INSERT INTO validated_input_snapshots (
                    snapshot_id, workflow_id, adapter_version, schema_version,
                    schema_dialect, canonical_payload, payload_digest_scheme,
                    payload_digest, validation_outcome, validation_issue_codes,
                    validated_at, expires_at, build_adapter_version, build_scheme,
                    build_logical_entrypoint, build_digest, build_captured_at,
                    consumed_run_id, consumed_at
                ) VALUES (
                    :snapshot_id, :workflow_id, '0.3.0', '1.0.0',
                    'https://json-schema.org/draft/2020-12/schema',
                    :payload, 'sha256-canonical-json-v1', :digest,
                    'adapter_validation_succeeded', :issues, :now, :expires,
                    '0.3.0', 'sha256-tree-v1', 'main.nf', :digest, :now,
                    :run_id, :now
                )
                """
            ),
            {
                "snapshot_id": snapshot_id,
                "workflow_id": workflow_id,
                "payload": '{"config":{"standard":{}},"samples":[],"options":{}}',
                "digest": digest,
                "issues": "[]",
                "now": now,
                "expires": "2026-07-18 12:00:00",
                "run_id": run_id,
            },
        )
        connection.execute(
            text(
                """
                INSERT INTO run_artifacts (
                    artifact_id, run_id, artifact_type, name, uri, mime_type,
                    produced_at, artifact_metadata
                ) VALUES (
                    :artifact_id, :run_id, 'file', 'summary.tsv',
                    :uri, 'text/tab-separated-values', :now, :metadata
                )
                """
            ),
            {
                "artifact_id": artifact_id,
                "run_id": run_id,
                "uri": f"run://runs/{run_id}/artifacts/{artifact_id}",
                "now": now,
                "metadata": '{"output_type":"qc_summary","sample_id":"S1"}',
            },
        )
        connection.execute(
            text(
                """
                INSERT INTO run_qc_metrics (
                    metric_id, run_id, metric_key, display_name, value_text,
                    unit, scope, sample_id, experiment_id, assay, qc_flag,
                    source_artifact_id, produced_at
                ) VALUES (
                    :metric_id, :run_id, 'sequencing.total_reads', 'Total reads',
                    '10', 'count', 'sample', 'S1', NULL, 'rnaseq', 'pass',
                    :artifact_id, :now
                )
                """
            ),
            {
                "metric_id": metric_id,
                "run_id": run_id,
                "artifact_id": artifact_id,
                "now": now,
            },
        )
    engine.dispose()

    upgrade_database(database_url)
    upgraded = create_database_engine(database_url)
    with upgraded.connect() as connection:
        assert connection.scalar(text("SELECT version_num FROM alembic_version")) == (
            "20260717_08"
        )
        run = (
            connection.execute(
                text(
                    "SELECT workflow_id, status, current_stage, tags "
                    "FROM runs WHERE run_id=:run_id"
                ),
                {"run_id": run_id},
            )
            .mappings()
            .one()
        )
        assert run["workflow_id"] == workflow_id
        assert run["status"] == "succeeded"
        assert run["current_stage"] == "results"
        assert json.loads(run["tags"]) == {"release_fixture": "v0.3.0"}
        event = (
            connection.execute(
                text(
                    "SELECT event_id, sequence, event_type, status, stage, message, "
                    "context, issue FROM run_events WHERE run_id=:run_id"
                ),
                {"run_id": run_id},
            )
            .mappings()
            .one()
        )
        assert dict(event) == {
            "event_id": "event-1",
            "sequence": 1,
            "event_type": "status_changed",
            "status": "succeeded",
            "stage": "results",
            "message": "Results indexed.",
            "context": '{"attempt_id":"attempt-1"}',
            "issue": None,
        }
        log = (
            connection.execute(
                text(
                    "SELECT chunk_id, stream_name, sequence, lines "
                    "FROM run_logs WHERE run_id=:run_id"
                ),
                {"run_id": run_id},
            )
            .mappings()
            .one()
        )
        assert dict(log) == {
            "chunk_id": "chunk-1",
            "stream_name": "stdout",
            "sequence": 1,
            "lines": '["workflow complete"]',
        }
        assignment = (
            connection.execute(
                text(
                    "SELECT job_id, backend, queue_name, dispatched_at, claimed_at, "
                    "cancellation_requested_at, cancellation_reason, "
                    "cancellation_acknowledged_at FROM run_execution_assignments "
                    "WHERE run_id=:run_id"
                ),
                {"run_id": run_id},
            )
            .mappings()
            .one()
        )
        assert assignment["job_id"] == "release-job"
        assert assignment["backend"] == "rq"
        assert assignment["queue_name"] == "release-queue"
        assert assignment["dispatched_at"] is not None
        assert assignment["claimed_at"] is not None
        assert assignment["cancellation_requested_at"] is None
        assert assignment["cancellation_reason"] is None
        assert assignment["cancellation_acknowledged_at"] is None
        build_identity = (
            connection.execute(
                text(
                    "SELECT workflow_id, adapter_version, scheme, logical_entrypoint, "
                    "digest FROM run_workflow_build_identities "
                    "WHERE run_id=:run_id"
                ),
                {"run_id": run_id},
            )
            .mappings()
            .one()
        )
        assert dict(build_identity) == {
            "workflow_id": workflow_id,
            "adapter_version": "0.3.0",
            "scheme": "sha256-tree-v1",
            "logical_entrypoint": "main.nf",
            "digest": digest,
        }
        snapshot = (
            connection.execute(
                text(
                    "SELECT workflow_id, adapter_version, schema_version, "
                    "canonical_payload, payload_digest_scheme, payload_digest, "
                    "validation_outcome, validation_issue_codes, "
                    "build_adapter_version, build_scheme, "
                    "build_logical_entrypoint, build_digest, consumed_run_id "
                    "FROM validated_input_snapshots "
                    "WHERE snapshot_id=:snapshot_id"
                ),
                {"snapshot_id": snapshot_id},
            )
            .mappings()
            .one()
        )
        assert snapshot["workflow_id"] == workflow_id
        assert snapshot["adapter_version"] == "0.3.0"
        assert snapshot["schema_version"] == "1.0.0"
        assert json.loads(snapshot["canonical_payload"]) == {
            "config": {"standard": {}},
            "samples": [],
            "options": {},
        }
        assert snapshot["payload_digest_scheme"] == "sha256-canonical-json-v1"
        assert snapshot["payload_digest"] == digest
        assert snapshot["validation_outcome"] == "adapter_validation_succeeded"
        assert json.loads(snapshot["validation_issue_codes"]) == []
        assert snapshot["build_adapter_version"] == "0.3.0"
        assert snapshot["build_scheme"] == "sha256-tree-v1"
        assert snapshot["build_logical_entrypoint"] == "main.nf"
        assert snapshot["build_digest"] == digest
        assert snapshot["consumed_run_id"] == run_id
        state = (
            connection.execute(
                text(
                    "SELECT artifact_revision, artifact_generation, qc_revision, "
                    "qc_generation FROM run_result_states WHERE run_id=:run_id"
                ),
                {"run_id": run_id},
            )
            .mappings()
            .one()
        )
        assert dict(state) == {
            "artifact_revision": 0,
            "artifact_generation": None,
            "qc_revision": 0,
            "qc_generation": None,
        }
        artifact = (
            connection.execute(
                text(
                    "SELECT artifact_type, name, uri, mime_type, "
                    "artifact_metadata, revision FROM run_artifacts "
                    "WHERE run_id=:run_id AND artifact_id=:artifact_id"
                ),
                {"run_id": run_id, "artifact_id": artifact_id},
            )
            .mappings()
            .one()
        )
        assert artifact["artifact_type"] == "file"
        assert artifact["name"] == "summary.tsv"
        assert artifact["uri"] == f"run://runs/{run_id}/artifacts/{artifact_id}"
        assert artifact["mime_type"] == "text/tab-separated-values"
        assert json.loads(artifact["artifact_metadata"]) == {
            "output_type": "qc_summary",
            "sample_id": "S1",
        }
        assert artifact["revision"] is None
        metric = (
            connection.execute(
                text(
                    "SELECT metric_key, display_name, value_text, unit, scope, "
                    "sample_id, experiment_id, assay, qc_flag, source_artifact_id "
                    "FROM run_qc_metrics WHERE run_id=:run_id"
                ),
                {"run_id": run_id},
            )
            .mappings()
            .one()
        )
        assert dict(metric) == {
            "metric_key": "sequencing.total_reads",
            "display_name": "Total reads",
            "value_text": "10",
            "unit": "count",
            "scope": "sample",
            "sample_id": "S1",
            "experiment_id": None,
            "assay": "rnaseq",
            "qc_flag": "pass",
            "source_artifact_id": artifact_id,
        }

    repository = SqlAlchemyRunRepository(create_session_factory(upgraded))
    with pytest.raises(ValueError, match="generation is unbound"):
        repository.list_artifacts(run_id)
    with pytest.raises(ValueError, match="generation is unbound"):
        repository.list_qc_metrics(run_id)
    upgraded.dispose()
