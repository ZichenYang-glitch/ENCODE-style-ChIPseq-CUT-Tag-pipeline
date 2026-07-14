"""Add bounded keyset indexes for read-only run history.

Revision ID: 20260714_07
Revises: 20260714_06
"""

from __future__ import annotations

from alembic import op


revision = "20260714_07"
down_revision = "20260714_06"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.create_index(
        "ix_runs_created_run_id",
        "runs",
        ["created_at", "run_id"],
    )
    op.create_index(
        "ix_runs_workflow_created_run_id",
        "runs",
        ["workflow_id", "created_at", "run_id"],
    )
    op.create_index(
        "ix_runs_status_created_run_id",
        "runs",
        ["status", "created_at", "run_id"],
    )
    op.create_index(
        "ix_runs_workflow_status_created_run_id",
        "runs",
        ["workflow_id", "status", "created_at", "run_id"],
    )


def downgrade() -> None:
    op.drop_index("ix_runs_workflow_status_created_run_id", table_name="runs")
    op.drop_index("ix_runs_status_created_run_id", table_name="runs")
    op.drop_index("ix_runs_workflow_created_run_id", table_name="runs")
    op.drop_index("ix_runs_created_run_id", table_name="runs")
