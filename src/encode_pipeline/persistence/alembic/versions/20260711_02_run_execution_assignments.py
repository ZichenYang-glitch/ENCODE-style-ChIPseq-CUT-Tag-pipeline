"""Add durable run execution assignments.

Revision ID: 20260711_02
Revises: 20260711_01
"""

from __future__ import annotations

from alembic import op
import sqlalchemy as sa


revision = "20260711_02"
down_revision = "20260711_01"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.create_table(
        "run_execution_assignments",
        sa.Column("run_id", sa.String(length=128), nullable=False),
        sa.Column("job_id", sa.String(length=255), nullable=False),
        sa.Column("backend", sa.String(length=64), nullable=False),
        sa.Column("queue_name", sa.String(length=128), nullable=False),
        sa.Column("created_at", sa.DateTime(timezone=True), nullable=False),
        sa.Column("dispatched_at", sa.DateTime(timezone=True), nullable=True),
        sa.Column("claimed_at", sa.DateTime(timezone=True), nullable=True),
        sa.CheckConstraint(
            "claimed_at IS NULL OR dispatched_at IS NOT NULL",
            name="ck_run_execution_assignments_claim_requires_dispatch",
        ),
        sa.ForeignKeyConstraint(["run_id"], ["runs.run_id"], ondelete="CASCADE"),
        sa.PrimaryKeyConstraint("run_id"),
        sa.UniqueConstraint(
            "job_id",
            name="uq_run_execution_assignments_job_id",
        ),
    )


def downgrade() -> None:
    op.drop_table("run_execution_assignments")
