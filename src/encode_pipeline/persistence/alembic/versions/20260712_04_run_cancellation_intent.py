"""Add durable execution cancellation intent and acknowledgement.

Revision ID: 20260712_04
Revises: 20260712_03
"""

from __future__ import annotations

from alembic import op
import sqlalchemy as sa


revision = "20260712_04"
down_revision = "20260712_03"
branch_labels = None
depends_on = None


def upgrade() -> None:
    with op.batch_alter_table("run_execution_assignments") as batch_op:
        batch_op.add_column(
            sa.Column(
                "cancellation_requested_at",
                sa.DateTime(timezone=True),
                nullable=True,
            )
        )
        batch_op.add_column(sa.Column("cancellation_reason", sa.Text(), nullable=True))
        batch_op.add_column(
            sa.Column(
                "cancellation_acknowledged_at",
                sa.DateTime(timezone=True),
                nullable=True,
            )
        )
        batch_op.create_check_constraint(
            "ck_run_execution_assignments_request_reason_pair",
            "(cancellation_requested_at IS NULL AND cancellation_reason IS NULL) "
            "OR (cancellation_requested_at IS NOT NULL AND "
            "cancellation_reason IS NOT NULL)",
        )
        batch_op.create_check_constraint(
            "ck_run_execution_assignments_request_requires_claim",
            "cancellation_requested_at IS NULL OR claimed_at IS NOT NULL",
        )
        batch_op.create_check_constraint(
            "ck_run_execution_assignments_ack_requires_request",
            "cancellation_acknowledged_at IS NULL "
            "OR cancellation_requested_at IS NOT NULL",
        )


def downgrade() -> None:
    with op.batch_alter_table("run_execution_assignments") as batch_op:
        batch_op.drop_constraint(
            "ck_run_execution_assignments_ack_requires_request",
            type_="check",
        )
        batch_op.drop_constraint(
            "ck_run_execution_assignments_request_requires_claim",
            type_="check",
        )
        batch_op.drop_constraint(
            "ck_run_execution_assignments_request_reason_pair",
            type_="check",
        )
        batch_op.drop_column("cancellation_acknowledged_at")
        batch_op.drop_column("cancellation_reason")
        batch_op.drop_column("cancellation_requested_at")
