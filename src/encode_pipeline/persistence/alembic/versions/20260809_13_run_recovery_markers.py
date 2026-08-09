"""Add one-shot administrator requeue request and confirmation markers.

Revision ID: 20260809_13
Revises: 20260807_12
"""

from __future__ import annotations

from alembic import op
import sqlalchemy as sa


revision = "20260809_13"
down_revision = "20260807_12"
branch_labels = None
depends_on = None


_ASSIGNMENT_TABLE = "run_execution_assignments"
_BATCH_TEMP_TABLE = "_alembic_tmp_run_execution_assignments"
_LEGACY_COLUMNS = (
    ("run_id", "VARCHAR(128)", False, None, 1),
    ("job_id", "VARCHAR(255)", False, None, 0),
    ("backend", "VARCHAR(64)", False, None, 0),
    ("queue_name", "VARCHAR(128)", False, None, 0),
    ("created_at", "DATETIME", False, None, 0),
    ("dispatched_at", "DATETIME", True, None, 0),
    ("claimed_at", "DATETIME", True, None, 0),
    ("cancellation_requested_at", "DATETIME", True, None, 0),
    ("cancellation_reason", "TEXT", True, None, 0),
    ("cancellation_acknowledged_at", "DATETIME", True, None, 0),
)
_RECOVERY_COLUMNS = _LEGACY_COLUMNS + (
    ("requeue_requested_at", "DATETIME", True, None, 0),
    ("requeue_confirmed_at", "DATETIME", True, None, 0),
)
_LEGACY_CHECKS = {
    (
        "ck_run_execution_assignments_ack_requires_request",
        "cancellation_acknowledged_at IS NULL OR cancellation_requested_at IS NOT NULL",
    ),
    (
        "ck_run_execution_assignments_claim_requires_dispatch",
        "claimed_at IS NULL OR dispatched_at IS NOT NULL",
    ),
    (
        "ck_run_execution_assignments_request_reason_pair",
        "(cancellation_requested_at IS NULL AND cancellation_reason IS NULL) "
        "OR (cancellation_requested_at IS NOT NULL AND "
        "cancellation_reason IS NOT NULL)",
    ),
    (
        "ck_run_execution_assignments_request_requires_claim",
        "cancellation_requested_at IS NULL OR claimed_at IS NOT NULL",
    ),
}
_RECOVERY_CHECKS = _LEGACY_CHECKS | {
    (
        "ck_run_execution_assignments_requeue_requires_dispatch",
        "requeue_requested_at IS NULL OR dispatched_at IS NOT NULL",
    ),
    (
        "ck_run_execution_assignments_requeue_confirm_requires_request",
        "requeue_confirmed_at IS NULL OR requeue_requested_at IS NOT NULL",
    ),
    (
        "ck_run_execution_assignments_requeue_confirmation_order",
        "requeue_confirmed_at IS NULL OR requeue_confirmed_at >= requeue_requested_at",
    ),
}


def _normalize_sql(value: object) -> str:
    return " ".join(str(value).split())


def _table_matches(
    inspector: sa.Inspector,
    table_name: str,
    *,
    columns: tuple[tuple[str, str, bool, object, int], ...],
    checks: set[tuple[str, str]],
) -> bool:
    actual_columns = tuple(
        (
            column["name"],
            str(column["type"]).upper(),
            bool(column["nullable"]),
            column["default"],
            int(column["primary_key"]),
        )
        for column in inspector.get_columns(table_name)
    )
    actual_checks = {
        (constraint["name"], _normalize_sql(constraint["sqltext"]))
        for constraint in inspector.get_check_constraints(table_name)
    }
    primary_key = inspector.get_pk_constraint(table_name)
    unique_constraints = {
        (constraint["name"], tuple(constraint["column_names"]))
        for constraint in inspector.get_unique_constraints(table_name)
    }
    foreign_keys = {
        (
            foreign_key["name"],
            tuple(foreign_key["constrained_columns"]),
            foreign_key["referred_schema"],
            foreign_key["referred_table"],
            tuple(foreign_key["referred_columns"]),
            tuple(sorted(foreign_key["options"].items())),
        )
        for foreign_key in inspector.get_foreign_keys(table_name)
    }
    return (
        actual_columns == columns
        and actual_checks == checks
        and primary_key["name"] is None
        and tuple(primary_key["constrained_columns"]) == ("run_id",)
        and unique_constraints == {("uq_run_execution_assignments_job_id", ("job_id",))}
        and inspector.get_indexes(table_name) == []
        and foreign_keys
        == {
            (
                None,
                ("run_id",),
                None,
                "runs",
                ("run_id",),
                (("ondelete", "CASCADE"),),
            )
        }
    )


def _batch_table_has_no_external_schema_dependencies(
    connection: sa.Connection,
    inspector: sa.Inspector,
) -> bool:
    references = connection.execute(
        sa.text(
            "SELECT type, name FROM sqlite_master "
            "WHERE sql IS NOT NULL "
            "AND lower(sql) LIKE :reference "
            "ORDER BY type, name"
        ),
        {"reference": f"%{_BATCH_TEMP_TABLE.lower()}%"},
    ).all()
    if references != [("table", _BATCH_TEMP_TABLE)]:
        return False
    return not any(
        foreign_key["referred_table"] == _BATCH_TEMP_TABLE
        for table_name in inspector.get_table_names()
        if table_name != _BATCH_TEMP_TABLE
        for foreign_key in inspector.get_foreign_keys(table_name)
    )


def _batch_rows_are_redundant(
    connection: sa.Connection,
    *,
    temp_has_recovery_columns: bool,
) -> bool:
    legacy_names = tuple(column[0] for column in _LEGACY_COLUMNS)
    equality = " AND ".join(
        f'staged."{name}" IS durable."{name}"' for name in legacy_names
    )
    marker_guard = (
        'staged."requeue_requested_at" IS NOT NULL OR '
        'staged."requeue_confirmed_at" IS NOT NULL OR '
        if temp_has_recovery_columns
        else ""
    )
    unsafe_row = connection.execute(
        sa.text(
            f'SELECT 1 FROM "{_BATCH_TEMP_TABLE}" AS staged '
            f'LEFT JOIN "{_ASSIGNMENT_TABLE}" AS durable '
            f'ON durable."run_id" IS staged."run_id" '
            f'WHERE {marker_guard}durable."run_id" IS NULL OR NOT ({equality}) '
            "LIMIT 1"
        )
    ).first()
    return unsafe_row is None


def _discard_proven_batch_residue(
    *,
    expected_revision: str,
    canonical_columns: tuple[tuple[str, str, bool, object, int], ...],
    canonical_checks: set[tuple[str, str]],
    temp_columns: tuple[tuple[str, str, bool, object, int], ...],
    temp_checks: set[tuple[str, str]],
) -> None:
    """Delete only a provably redundant table from this exact batch rebuild."""
    connection = op.get_bind()
    inspector = sa.inspect(connection)
    table_names = set(inspector.get_table_names())
    if _BATCH_TEMP_TABLE not in table_names:
        return
    proven = (
        connection.dialect.name == "sqlite"
        and _ASSIGNMENT_TABLE in table_names
        and connection.scalar(sa.text("SELECT version_num FROM alembic_version"))
        == expected_revision
        and _table_matches(
            inspector,
            _ASSIGNMENT_TABLE,
            columns=canonical_columns,
            checks=canonical_checks,
        )
        and _table_matches(
            inspector,
            _BATCH_TEMP_TABLE,
            columns=temp_columns,
            checks=temp_checks,
        )
        and _batch_table_has_no_external_schema_dependencies(connection, inspector)
        and _batch_rows_are_redundant(
            connection,
            temp_has_recovery_columns=temp_columns == _RECOVERY_COLUMNS,
        )
    )
    if not proven:
        raise RuntimeError("run recovery migration found an ambiguous batch table")
    op.execute(sa.text(f'DROP TABLE "{_BATCH_TEMP_TABLE}"'))


def upgrade() -> None:
    _discard_proven_batch_residue(
        expected_revision=down_revision,
        canonical_columns=_LEGACY_COLUMNS,
        canonical_checks=_LEGACY_CHECKS,
        temp_columns=_RECOVERY_COLUMNS,
        temp_checks=_RECOVERY_CHECKS,
    )
    with op.batch_alter_table("run_execution_assignments") as batch_op:
        batch_op.add_column(
            sa.Column(
                "requeue_requested_at",
                sa.DateTime(timezone=True),
                nullable=True,
            )
        )
        batch_op.add_column(
            sa.Column(
                "requeue_confirmed_at",
                sa.DateTime(timezone=True),
                nullable=True,
            )
        )
        batch_op.create_check_constraint(
            "ck_run_execution_assignments_requeue_requires_dispatch",
            "requeue_requested_at IS NULL OR dispatched_at IS NOT NULL",
        )
        batch_op.create_check_constraint(
            "ck_run_execution_assignments_requeue_confirm_requires_request",
            "requeue_confirmed_at IS NULL OR requeue_requested_at IS NOT NULL",
        )
        batch_op.create_check_constraint(
            "ck_run_execution_assignments_requeue_confirmation_order",
            "requeue_confirmed_at IS NULL "
            "OR requeue_confirmed_at >= requeue_requested_at",
        )


def downgrade() -> None:
    _discard_proven_batch_residue(
        expected_revision=revision,
        canonical_columns=_RECOVERY_COLUMNS,
        canonical_checks=_RECOVERY_CHECKS,
        temp_columns=_LEGACY_COLUMNS,
        temp_checks=_LEGACY_CHECKS,
    )
    with op.batch_alter_table("run_execution_assignments") as batch_op:
        batch_op.drop_constraint(
            "ck_run_execution_assignments_requeue_confirmation_order",
            type_="check",
        )
        batch_op.drop_constraint(
            "ck_run_execution_assignments_requeue_confirm_requires_request",
            type_="check",
        )
        batch_op.drop_constraint(
            "ck_run_execution_assignments_requeue_requires_dispatch",
            type_="check",
        )
        batch_op.drop_column("requeue_confirmed_at")
        batch_op.drop_column("requeue_requested_at")
