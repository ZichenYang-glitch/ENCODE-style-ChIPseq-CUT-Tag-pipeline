"""Add durable workflow build identities.

Revision ID: 20260712_03
Revises: 20260711_02
"""

from __future__ import annotations

from alembic import op
import sqlalchemy as sa


revision = "20260712_03"
down_revision = "20260711_02"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.create_table(
        "run_workflow_build_identities",
        sa.Column("run_id", sa.String(length=128), nullable=False),
        sa.Column("workflow_id", sa.String(length=255), nullable=False),
        sa.Column("adapter_version", sa.String(length=255), nullable=False),
        sa.Column("scheme", sa.String(length=64), nullable=False),
        sa.Column("logical_entrypoint", sa.String(length=512), nullable=False),
        sa.Column("digest", sa.String(length=64), nullable=False),
        sa.Column("captured_at", sa.DateTime(timezone=True), nullable=False),
        sa.ForeignKeyConstraint(["run_id"], ["runs.run_id"], ondelete="CASCADE"),
        sa.PrimaryKeyConstraint("run_id"),
    )


def downgrade() -> None:
    op.drop_table("run_workflow_build_identities")
