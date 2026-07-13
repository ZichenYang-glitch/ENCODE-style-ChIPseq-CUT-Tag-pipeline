"""Add immutable validated workflow input snapshots.

Revision ID: 20260714_06
Revises: 20260712_05
"""

from __future__ import annotations

from alembic import op
import sqlalchemy as sa


revision = "20260714_06"
down_revision = "20260712_05"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.create_table(
        "validated_input_snapshots",
        sa.Column("snapshot_id", sa.String(length=64), nullable=False),
        sa.Column("workflow_id", sa.String(length=255), nullable=False),
        sa.Column("adapter_version", sa.String(length=128), nullable=False),
        sa.Column("schema_version", sa.String(length=64), nullable=False),
        sa.Column("schema_dialect", sa.String(length=255), nullable=False),
        sa.Column("canonical_payload", sa.Text(), nullable=False),
        sa.Column("payload_digest_scheme", sa.String(length=64), nullable=False),
        sa.Column("payload_digest", sa.String(length=64), nullable=False),
        sa.Column("validation_outcome", sa.String(length=64), nullable=False),
        sa.Column("validation_issue_codes", sa.JSON(), nullable=False),
        sa.Column("validated_at", sa.DateTime(timezone=True), nullable=False),
        sa.Column("expires_at", sa.DateTime(timezone=True), nullable=False),
        sa.Column("build_adapter_version", sa.String(length=128), nullable=False),
        sa.Column("build_scheme", sa.String(length=64), nullable=False),
        sa.Column("build_logical_entrypoint", sa.String(length=512), nullable=False),
        sa.Column("build_digest", sa.String(length=64), nullable=False),
        sa.Column("build_captured_at", sa.DateTime(timezone=True), nullable=False),
        sa.Column("consumed_run_id", sa.String(length=128), nullable=True),
        sa.Column("consumed_at", sa.DateTime(timezone=True), nullable=True),
        sa.CheckConstraint(
            "(consumed_run_id IS NULL AND consumed_at IS NULL) OR "
            "(consumed_run_id IS NOT NULL AND consumed_at IS NOT NULL)",
            name="ck_validated_input_snapshots_consumption_pair",
        ),
        sa.CheckConstraint(
            "expires_at > validated_at",
            name="ck_validated_input_snapshots_expiry",
        ),
        sa.CheckConstraint(
            "validation_outcome = 'adapter_validation_succeeded'",
            name="ck_validated_input_snapshots_success",
        ),
        sa.CheckConstraint(
            "length(payload_digest) = 64",
            name="ck_validated_input_snapshots_digest_length",
        ),
        sa.ForeignKeyConstraint(
            ["consumed_run_id"],
            ["runs.run_id"],
            ondelete="RESTRICT",
        ),
        sa.PrimaryKeyConstraint("snapshot_id"),
        sa.UniqueConstraint(
            "consumed_run_id",
            name="uq_validated_input_snapshots_consumed_run_id",
        ),
    )
    op.create_index(
        "ix_validated_input_snapshots_workflow_id",
        "validated_input_snapshots",
        ["workflow_id"],
    )
    op.create_index(
        "ix_validated_input_snapshots_expires_at",
        "validated_input_snapshots",
        ["expires_at"],
    )


def downgrade() -> None:
    op.drop_index(
        "ix_validated_input_snapshots_expires_at",
        table_name="validated_input_snapshots",
    )
    op.drop_index(
        "ix_validated_input_snapshots_workflow_id",
        table_name="validated_input_snapshots",
    )
    op.drop_table("validated_input_snapshots")
