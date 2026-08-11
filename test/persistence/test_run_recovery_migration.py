"""Run-recovery marker migration integrity and reversibility tests."""

from __future__ import annotations

from contextlib import contextmanager

import pytest
from alembic import op as alembic_op
from sqlalchemy import inspect, text
from sqlalchemy.exc import IntegrityError

from encode_pipeline.persistence import create_database_engine
from encode_pipeline.persistence.migration_admission import (
    verify_migration_execution_inventory,
)
from encode_pipeline.persistence.migrations import (
    downgrade_database,
    upgrade_database,
)


PRIOR_REVISION = "20260807_12"
RECOVERY_REVISION = "20260809_13"
LEGACY_ASSIGNMENT_COLUMNS = {
    "run_id",
    "job_id",
    "backend",
    "queue_name",
    "created_at",
    "dispatched_at",
    "claimed_at",
    "cancellation_requested_at",
    "cancellation_reason",
    "cancellation_acknowledged_at",
}
RECOVERY_COLUMNS = {
    "managed_container_scope",
    "managed_container_endpoint_identity",
    "requeue_requested_at",
    "requeue_confirmed_at",
}
LEGACY_ASSIGNMENT_CONSTRAINTS = {
    "ck_run_execution_assignments_ack_requires_request",
    "ck_run_execution_assignments_claim_requires_dispatch",
    "ck_run_execution_assignments_request_reason_pair",
    "ck_run_execution_assignments_request_requires_claim",
}
RECOVERY_CONSTRAINTS = {
    "ck_run_execution_assignments_cleanup_endpoint_format",
    "ck_run_execution_assignments_cleanup_identity_pair",
    "ck_run_execution_assignments_cleanup_scope_format",
    "ck_run_execution_assignments_requeue_requires_dispatch",
    "ck_run_execution_assignments_requeue_confirm_requires_request",
    "ck_run_execution_assignments_requeue_confirmation_order",
}


def test_rev12_to_rev13_zero_backfills_and_preserves_assignment_evidence(
    tmp_path,
) -> None:
    database_url = f"sqlite:///{tmp_path / 'legacy.db'}"
    upgrade_database(database_url, PRIOR_REVISION)
    engine = create_database_engine(database_url)
    with engine.begin() as connection:
        _insert_run(connection, "queued-run", status="queued")
        _insert_run(connection, "cancelled-run", status="cancelled")
        _insert_assignment(
            connection,
            run_id="queued-run",
            job_id="queued-job",
            dispatched_at="2026-08-09 10:01:00",
        )
        _insert_assignment(
            connection,
            run_id="cancelled-run",
            job_id="cancelled-job",
            dispatched_at="2026-08-09 10:01:00",
            claimed_at="2026-08-09 10:02:00",
            cancellation_requested_at="2026-08-09 10:03:00",
            cancellation_reason="operator request",
            cancellation_acknowledged_at="2026-08-09 10:04:00",
        )
    engine.dispose()

    upgrade_database(database_url)
    engine = create_database_engine(database_url)
    with engine.connect() as connection:
        rows = (
            connection.execute(
                text(
                    "SELECT run_id, job_id, backend, queue_name, created_at, "
                    "dispatched_at, claimed_at, cancellation_requested_at, "
                    "cancellation_reason, cancellation_acknowledged_at, "
                    "managed_container_scope, managed_container_endpoint_identity, "
                    "requeue_requested_at, requeue_confirmed_at "
                    "FROM run_execution_assignments ORDER BY run_id"
                )
            )
            .mappings()
            .all()
        )
        assert connection.scalar(text("SELECT version_num FROM alembic_version")) == (
            RECOVERY_REVISION
        )

    assert [row["run_id"] for row in rows] == ["cancelled-run", "queued-run"]
    cancelled, queued = rows
    assert dict(cancelled) == {
        "run_id": "cancelled-run",
        "job_id": "cancelled-job",
        "backend": "rq",
        "queue_name": "default",
        "created_at": "2026-08-09 10:00:00",
        "dispatched_at": "2026-08-09 10:01:00",
        "claimed_at": "2026-08-09 10:02:00",
        "cancellation_requested_at": "2026-08-09 10:03:00",
        "cancellation_reason": "operator request",
        "cancellation_acknowledged_at": "2026-08-09 10:04:00",
        "managed_container_scope": None,
        "managed_container_endpoint_identity": None,
        "requeue_requested_at": None,
        "requeue_confirmed_at": None,
    }
    assert queued["job_id"] == "queued-job"
    assert queued["dispatched_at"] == "2026-08-09 10:01:00"
    assert queued["claimed_at"] is None
    assert queued["managed_container_scope"] is None
    assert queued["managed_container_endpoint_identity"] is None
    assert queued["requeue_requested_at"] is None
    assert queued["requeue_confirmed_at"] is None
    engine.dispose()


def test_rev13_is_the_sole_head_and_enforces_requeue_marker_constraints(
    tmp_path,
) -> None:
    database_url = f"sqlite:///{tmp_path / 'platform.db'}"
    upgrade_database(database_url)
    engine = create_database_engine(database_url)
    inspector = inspect(engine)

    inventory = verify_migration_execution_inventory()
    assert inventory.heads == (RECOVERY_REVISION,)
    assert {
        column["name"] for column in inspector.get_columns("run_execution_assignments")
    } == LEGACY_ASSIGNMENT_COLUMNS | RECOVERY_COLUMNS
    assert {
        column["name"]: column["nullable"]
        for column in inspector.get_columns("run_execution_assignments")
        if column["name"] in RECOVERY_COLUMNS
    } == {
        "managed_container_scope": True,
        "managed_container_endpoint_identity": True,
        "requeue_requested_at": True,
        "requeue_confirmed_at": True,
    }
    constraint_names = {
        constraint["name"]
        for constraint in inspector.get_check_constraints("run_execution_assignments")
    }
    assert constraint_names == LEGACY_ASSIGNMENT_CONSTRAINTS | RECOVERY_CONSTRAINTS

    invalid_assignments = (
        {
            "run_id": "scope-without-endpoint",
            "job_id": "job-scope-without-endpoint",
            "dispatched_at": "2026-08-09 11:00:00",
            "scope": "a" * 64,
            "endpoint": None,
            "requeue_requested_at": None,
            "requeue_confirmed_at": None,
        },
        {
            "run_id": "endpoint-without-scope",
            "job_id": "job-endpoint-without-scope",
            "dispatched_at": "2026-08-09 11:00:00",
            "scope": None,
            "endpoint": "b" * 64,
            "requeue_requested_at": None,
            "requeue_confirmed_at": None,
        },
        {
            "run_id": "short-scope",
            "job_id": "job-short-scope",
            "dispatched_at": "2026-08-09 11:00:00",
            "scope": "a" * 63,
            "endpoint": "b" * 64,
            "requeue_requested_at": None,
            "requeue_confirmed_at": None,
        },
        {
            "run_id": "uppercase-endpoint",
            "job_id": "job-uppercase-endpoint",
            "dispatched_at": "2026-08-09 11:00:00",
            "scope": "a" * 64,
            "endpoint": "B" * 64,
            "requeue_requested_at": None,
            "requeue_confirmed_at": None,
        },
        {
            "run_id": "request-without-dispatch",
            "job_id": "job-request-without-dispatch",
            "dispatched_at": None,
            "scope": None,
            "endpoint": None,
            "requeue_requested_at": "2026-08-09 11:01:00",
            "requeue_confirmed_at": None,
        },
        {
            "run_id": "confirmation-without-request",
            "job_id": "job-confirmation-without-request",
            "dispatched_at": "2026-08-09 11:00:00",
            "scope": None,
            "endpoint": None,
            "requeue_requested_at": None,
            "requeue_confirmed_at": "2026-08-09 11:01:00",
        },
        {
            "run_id": "confirmation-before-request",
            "job_id": "job-confirmation-before-request",
            "dispatched_at": "2026-08-09 11:00:00",
            "scope": None,
            "endpoint": None,
            "requeue_requested_at": "2026-08-09 11:02:00",
            "requeue_confirmed_at": "2026-08-09 11:01:00",
        },
    )
    with engine.begin() as connection:
        for assignment in invalid_assignments:
            _insert_run(connection, assignment["run_id"], status="queued")

    for assignment in invalid_assignments:
        with pytest.raises(IntegrityError):
            with engine.begin() as connection:
                connection.execute(
                    text(
                        "INSERT INTO run_execution_assignments "
                        "(run_id, job_id, backend, queue_name, created_at, "
                        "dispatched_at, managed_container_scope, "
                        "managed_container_endpoint_identity, requeue_requested_at, "
                        "requeue_confirmed_at) VALUES "
                        "(:run_id, :job_id, 'rq', 'default', "
                        "'2026-08-09 11:00:00', :dispatched_at, :scope, :endpoint, "
                        ":requeue_requested_at, :requeue_confirmed_at)"
                    ),
                    assignment,
                )

    with engine.connect() as connection:
        assert connection.scalar(text("SELECT version_num FROM alembic_version")) == (
            RECOVERY_REVISION
        )
        assert (
            connection.scalar(text("SELECT count(*) FROM run_execution_assignments"))
            == 0
        )
    engine.dispose()


def test_rev13_downgrade_removes_only_recovery_fields(tmp_path) -> None:
    database_url = f"sqlite:///{tmp_path / 'platform.db'}"
    upgrade_database(database_url)
    engine = create_database_engine(database_url)
    with engine.begin() as connection:
        _insert_run(connection, "queued-run", status="queued")
        connection.execute(
            text(
                "INSERT INTO run_execution_assignments "
                "(run_id, job_id, backend, queue_name, created_at, dispatched_at, "
                "requeue_requested_at, requeue_confirmed_at) VALUES "
                "('queued-run', 'queued-job', 'rq', 'default', "
                "'2026-08-09 12:00:00', '2026-08-09 12:01:00', "
                "'2026-08-09 12:02:00', '2026-08-09 12:03:00')"
            )
        )
    engine.dispose()

    downgrade_database(database_url, PRIOR_REVISION)
    engine = create_database_engine(database_url)
    inspector = inspect(engine)
    assert "_alembic_tmp_run_execution_assignments" not in inspector.get_table_names()
    assert {
        column["name"] for column in inspector.get_columns("run_execution_assignments")
    } == LEGACY_ASSIGNMENT_COLUMNS
    assert not (
        RECOVERY_CONSTRAINTS
        & {
            constraint["name"]
            for constraint in inspector.get_check_constraints(
                "run_execution_assignments"
            )
        }
    )
    with engine.connect() as connection:
        assignment = (
            connection.execute(
                text(
                    "SELECT run_id, job_id, dispatched_at "
                    "FROM run_execution_assignments"
                )
            )
            .mappings()
            .one()
        )
        assert dict(assignment) == {
            "run_id": "queued-run",
            "job_id": "queued-job",
            "dispatched_at": "2026-08-09 12:01:00",
        }
        assert connection.scalar(text("SELECT version_num FROM alembic_version")) == (
            PRIOR_REVISION
        )
    engine.dispose()


def test_rev13_upgrade_failure_rolls_back_and_retries_only_redundant_residue(
    tmp_path,
    monkeypatch,
) -> None:
    database_url = f"sqlite:///{tmp_path / 'platform.db'}"
    upgrade_database(database_url, PRIOR_REVISION)
    engine = create_database_engine(database_url)
    with engine.begin() as connection:
        _insert_run(connection, "queued-run", status="queued")
        _insert_assignment(
            connection,
            run_id="queued-run",
            job_id="queued-job",
            dispatched_at="2026-08-09 13:01:00",
        )
    engine.dispose()

    original_batch_alter_table = alembic_op.batch_alter_table

    @contextmanager
    def fail_after_batch_rebuild(*args, **kwargs):
        with original_batch_alter_table(*args, **kwargs) as batch_op:
            yield batch_op
        raise RuntimeError("injected post-DDL migration failure")

    with monkeypatch.context() as patch:
        patch.setattr(
            alembic_op,
            "batch_alter_table",
            fail_after_batch_rebuild,
        )
        with pytest.raises(
            RuntimeError,
            match="injected post-DDL migration failure",
        ):
            upgrade_database(database_url)

    engine = create_database_engine(database_url)
    inspector = inspect(engine)
    assert "_alembic_tmp_run_execution_assignments" in inspector.get_table_names()
    assert {
        column["name"] for column in inspector.get_columns("run_execution_assignments")
    } == LEGACY_ASSIGNMENT_COLUMNS
    with engine.connect() as connection:
        assert connection.scalar(text("SELECT version_num FROM alembic_version")) == (
            PRIOR_REVISION
        )
        assert connection.execute(
            text("SELECT run_id, job_id, dispatched_at FROM run_execution_assignments")
        ).one() == ("queued-run", "queued-job", "2026-08-09 13:01:00")
    engine.dispose()

    upgrade_database(database_url)
    engine = create_database_engine(database_url)
    assert (
        "_alembic_tmp_run_execution_assignments"
        not in inspect(engine).get_table_names()
    )
    with engine.connect() as connection:
        assert connection.scalar(text("SELECT version_num FROM alembic_version")) == (
            RECOVERY_REVISION
        )
        assert connection.execute(
            text(
                "SELECT requeue_requested_at, requeue_confirmed_at "
                "FROM run_execution_assignments"
            )
        ).one() == (None, None)
    engine.dispose()


def test_rev13_refuses_preexisting_batch_table_without_deleting_unknown_data(
    tmp_path,
) -> None:
    database_url = f"sqlite:///{tmp_path / 'platform.db'}"
    upgrade_database(database_url, PRIOR_REVISION)
    engine = create_database_engine(database_url)
    with engine.begin() as connection:
        _insert_run(connection, "queued-run", status="queued")
        _insert_assignment(
            connection,
            run_id="queued-run",
            job_id="queued-job",
            dispatched_at="2026-08-09 14:01:00",
        )
        connection.execute(
            text(
                'CREATE TABLE "_alembic_tmp_run_execution_assignments" '
                "(marker TEXT PRIMARY KEY, payload TEXT NOT NULL)"
            )
        )
        connection.execute(
            text(
                'INSERT INTO "_alembic_tmp_run_execution_assignments" '
                "(marker, payload) VALUES ('operator-owned', 'must-survive')"
            )
        )
    engine.dispose()

    with pytest.raises(RuntimeError, match="ambiguous batch table"):
        upgrade_database(database_url)

    engine = create_database_engine(database_url)
    inspector = inspect(engine)
    assert {
        column["name"] for column in inspector.get_columns("run_execution_assignments")
    } == LEGACY_ASSIGNMENT_COLUMNS
    assert {
        column["name"]
        for column in inspector.get_columns("_alembic_tmp_run_execution_assignments")
    } == {"marker", "payload"}
    with engine.connect() as connection:
        assert connection.scalar(text("SELECT version_num FROM alembic_version")) == (
            PRIOR_REVISION
        )
        assert connection.execute(
            text("SELECT run_id, job_id, dispatched_at FROM run_execution_assignments")
        ).one() == ("queued-run", "queued-job", "2026-08-09 14:01:00")
        assert connection.execute(
            text('SELECT marker, payload FROM "_alembic_tmp_run_execution_assignments"')
        ).one() == ("operator-owned", "must-survive")
    engine.dispose()


def test_rev13_refuses_exact_batch_shape_when_it_contains_unique_data(
    tmp_path,
    monkeypatch,
) -> None:
    database_url = f"sqlite:///{tmp_path / 'platform.db'}"
    upgrade_database(database_url, PRIOR_REVISION)
    engine = create_database_engine(database_url)
    with engine.begin() as connection:
        _insert_run(connection, "queued-run", status="queued")
        _insert_assignment(
            connection,
            run_id="queued-run",
            job_id="queued-job",
            dispatched_at="2026-08-09 15:01:00",
        )
    engine.dispose()

    original_batch_alter_table = alembic_op.batch_alter_table

    @contextmanager
    def fail_after_batch_rebuild(*args, **kwargs):
        with original_batch_alter_table(*args, **kwargs) as batch_op:
            yield batch_op
        raise RuntimeError("injected post-DDL migration failure")

    with monkeypatch.context() as patch:
        patch.setattr(
            alembic_op,
            "batch_alter_table",
            fail_after_batch_rebuild,
        )
        with pytest.raises(RuntimeError, match="injected post-DDL"):
            upgrade_database(database_url)

    engine = create_database_engine(database_url)
    with engine.begin() as connection:
        _insert_run(connection, "temp-only-run", status="queued")
        connection.execute(
            text(
                'INSERT INTO "_alembic_tmp_run_execution_assignments" '
                "(run_id, job_id, backend, queue_name, created_at, dispatched_at, "
                "requeue_requested_at, requeue_confirmed_at) VALUES "
                "('temp-only-run', 'temp-only-job', 'rq', 'default', "
                "'2026-08-09 15:00:00', '2026-08-09 15:01:00', NULL, NULL)"
            )
        )
    engine.dispose()

    with pytest.raises(RuntimeError, match="ambiguous batch table"):
        upgrade_database(database_url)

    engine = create_database_engine(database_url)
    with engine.connect() as connection:
        assert connection.scalar(text("SELECT version_num FROM alembic_version")) == (
            PRIOR_REVISION
        )
        assert connection.execute(
            text("SELECT run_id, job_id FROM run_execution_assignments")
        ).one() == ("queued-run", "queued-job")
        assert connection.execute(
            text('SELECT run_id, job_id FROM "_alembic_tmp_run_execution_assignments"')
        ).one() == ("temp-only-run", "temp-only-job")
    engine.dispose()


def test_rev13_refuses_batch_residue_with_cleanup_identity_evidence(
    tmp_path,
    monkeypatch,
) -> None:
    database_url = f"sqlite:///{tmp_path / 'platform.db'}"
    upgrade_database(database_url, PRIOR_REVISION)
    engine = create_database_engine(database_url)
    with engine.begin() as connection:
        _insert_run(connection, "queued-run", status="queued")
        _insert_assignment(
            connection,
            run_id="queued-run",
            job_id="queued-job",
            dispatched_at="2026-08-09 16:01:00",
        )
    engine.dispose()

    original_batch_alter_table = alembic_op.batch_alter_table

    @contextmanager
    def fail_after_batch_rebuild(*args, **kwargs):
        with original_batch_alter_table(*args, **kwargs) as batch_op:
            yield batch_op
        raise RuntimeError("injected post-DDL migration failure")

    with monkeypatch.context() as patch:
        patch.setattr(
            alembic_op,
            "batch_alter_table",
            fail_after_batch_rebuild,
        )
        with pytest.raises(RuntimeError, match="injected post-DDL"):
            upgrade_database(database_url)

    engine = create_database_engine(database_url)
    with engine.begin() as connection:
        connection.execute(
            text(
                'INSERT INTO "_alembic_tmp_run_execution_assignments" '
                "(run_id, job_id, backend, queue_name, created_at, dispatched_at, "
                "managed_container_scope, managed_container_endpoint_identity) "
                "VALUES ('queued-run', 'queued-job', 'rq', 'default', "
                "'2026-08-09 10:00:00', '2026-08-09 16:01:00', "
                ":scope, :endpoint)"
            ),
            {"scope": "a" * 64, "endpoint": "b" * 64},
        )
        assert connection.execute(
            text(
                "SELECT managed_container_scope, "
                "managed_container_endpoint_identity "
                'FROM "_alembic_tmp_run_execution_assignments" '
                "WHERE run_id = 'queued-run'"
            )
        ).one() == ("a" * 64, "b" * 64)
    engine.dispose()

    with pytest.raises(RuntimeError, match="ambiguous batch table"):
        upgrade_database(database_url)

    engine = create_database_engine(database_url)
    with engine.connect() as connection:
        assert connection.scalar(text("SELECT version_num FROM alembic_version")) == (
            PRIOR_REVISION
        )
        staged = connection.execute(
            text(
                "SELECT managed_container_scope, "
                "managed_container_endpoint_identity "
                'FROM "_alembic_tmp_run_execution_assignments" '
                "WHERE run_id = 'queued-run'"
            )
        ).one()
        assert staged == ("a" * 64, "b" * 64)
    engine.dispose()


def _insert_run(connection, run_id: str, *, status: str) -> None:
    connection.execute(
        text(
            "INSERT INTO runs "
            "(run_id, workflow_id, inputs, status, created_at, updated_at, tags) "
            "VALUES (:run_id, 'workflow-a', '{}', :status, "
            "'2026-08-09 10:00:00', '2026-08-09 10:00:00', '{}')"
        ),
        {"run_id": run_id, "status": status},
    )


def _insert_assignment(
    connection,
    *,
    run_id: str,
    job_id: str,
    dispatched_at: str | None,
    claimed_at: str | None = None,
    cancellation_requested_at: str | None = None,
    cancellation_reason: str | None = None,
    cancellation_acknowledged_at: str | None = None,
) -> None:
    connection.execute(
        text(
            "INSERT INTO run_execution_assignments "
            "(run_id, job_id, backend, queue_name, created_at, dispatched_at, "
            "claimed_at, cancellation_requested_at, cancellation_reason, "
            "cancellation_acknowledged_at) VALUES "
            "(:run_id, :job_id, 'rq', 'default', '2026-08-09 10:00:00', "
            ":dispatched_at, :claimed_at, :cancellation_requested_at, "
            ":cancellation_reason, :cancellation_acknowledged_at)"
        ),
        {
            "run_id": run_id,
            "job_id": job_id,
            "dispatched_at": dispatched_at,
            "claimed_at": claimed_at,
            "cancellation_requested_at": cancellation_requested_at,
            "cancellation_reason": cancellation_reason,
            "cancellation_acknowledged_at": cancellation_acknowledged_at,
        },
    )
