"""Bind artifact and QC results to durable monotonic generations.

Revision ID: 20260717_08
Revises: 20260714_07
"""

from __future__ import annotations

from alembic import op
import sqlalchemy as sa


revision = "20260717_08"
down_revision = "20260714_07"
branch_labels = None
depends_on = None


def upgrade() -> None:
    # Nullable by design: existing artifact rows cannot be assigned a content or
    # descriptor identity without re-opening mutable external files. They remain
    # unbound and repository reads fail closed until a complete re-extraction.
    op.add_column(
        "run_artifacts",
        sa.Column("revision", sa.String(length=76), nullable=True),
    )
    op.create_table(
        "run_result_states",
        sa.Column("run_id", sa.String(length=128), nullable=False),
        sa.Column(
            "artifact_revision", sa.Integer(), server_default="0", nullable=False
        ),
        sa.Column("artifact_generation", sa.String(length=76), nullable=True),
        sa.Column("artifact_manifest_digest", sa.String(length=64), nullable=True),
        sa.Column("artifact_attempt_id", sa.String(length=78), nullable=True),
        sa.Column("artifact_attempt_status", sa.String(length=16), nullable=True),
        sa.Column("artifact_outcome", sa.String(length=16), nullable=True),
        sa.Column("artifact_reason_code", sa.String(length=128), nullable=True),
        sa.Column("qc_revision", sa.Integer(), server_default="0", nullable=False),
        sa.Column("qc_generation", sa.String(length=70), nullable=True),
        sa.Column("qc_manifest_digest", sa.String(length=64), nullable=True),
        sa.Column("qc_attempt_id", sa.String(length=78), nullable=True),
        sa.Column("qc_attempt_status", sa.String(length=16), nullable=True),
        sa.Column(
            "qc_attempt_artifact_generation", sa.String(length=76), nullable=True
        ),
        sa.Column("qc_artifact_generation", sa.String(length=76), nullable=True),
        sa.Column("qc_outcome", sa.String(length=16), nullable=True),
        sa.Column("qc_reason_code", sa.String(length=128), nullable=True),
        sa.CheckConstraint(
            "artifact_revision >= 0 AND qc_revision >= 0",
            name="ck_run_result_states_nonnegative_revisions",
        ),
        sa.CheckConstraint(
            "(artifact_revision = 0 AND artifact_generation IS NULL AND "
            "artifact_manifest_digest IS NULL) OR "
            "(artifact_revision > 0 AND artifact_generation IS NOT NULL AND "
            "artifact_manifest_digest IS NOT NULL)",
            name="ck_run_result_states_artifact_binding",
        ),
        sa.CheckConstraint(
            "(qc_revision = 0 AND qc_generation IS NULL AND "
            "qc_manifest_digest IS NULL AND qc_artifact_generation IS NULL) OR "
            "(qc_revision > 0 AND ((qc_generation IS NULL AND "
            "qc_manifest_digest IS NULL AND qc_artifact_generation IS NULL) OR "
            "(qc_generation IS NOT NULL AND qc_manifest_digest IS NOT NULL AND "
            "qc_artifact_generation IS NOT NULL)))",
            name="ck_run_result_states_qc_binding",
        ),
        sa.CheckConstraint(
            "(artifact_attempt_id IS NULL AND artifact_attempt_status IS NULL) OR "
            "(artifact_attempt_id IS NOT NULL AND artifact_attempt_status IN "
            "('pending', 'succeeded', 'failed'))",
            name="ck_run_result_states_artifact_attempt",
        ),
        sa.CheckConstraint(
            "(qc_attempt_id IS NULL AND qc_attempt_status IS NULL AND "
            "qc_attempt_artifact_generation IS NULL) OR "
            "(qc_attempt_id IS NOT NULL AND qc_attempt_status IN "
            "('pending', 'succeeded', 'failed') AND "
            "qc_attempt_artifact_generation IS NOT NULL)",
            name="ck_run_result_states_qc_attempt",
        ),
        sa.CheckConstraint(
            "artifact_outcome IS NULL OR artifact_outcome IN ('succeeded', 'failed')",
            name="ck_run_result_states_artifact_outcome",
        ),
        sa.CheckConstraint(
            "qc_outcome IS NULL OR qc_outcome IN "
            "('succeeded', 'failed', 'invalidated')",
            name="ck_run_result_states_qc_outcome",
        ),
        sa.CheckConstraint(
            "(artifact_outcome = 'failed' AND artifact_reason_code IS NOT NULL) OR "
            "(artifact_outcome IS NULL AND artifact_reason_code IS NULL) OR "
            "(artifact_outcome = 'succeeded' AND artifact_reason_code IS NULL)",
            name="ck_run_result_states_artifact_reason",
        ),
        sa.CheckConstraint(
            "(qc_outcome = 'failed' AND qc_reason_code IS NOT NULL) OR "
            "(qc_outcome IS NULL AND qc_reason_code IS NULL) OR "
            "(qc_outcome IN ('succeeded', 'invalidated') AND qc_reason_code IS NULL)",
            name="ck_run_result_states_qc_reason",
        ),
        sa.ForeignKeyConstraint(["run_id"], ["runs.run_id"], ondelete="CASCADE"),
        sa.PrimaryKeyConstraint("run_id"),
    )
    op.execute(
        sa.text(
            "INSERT INTO run_result_states "
            "(run_id, artifact_revision, qc_revision) "
            "SELECT run_id, 0, 0 FROM runs"
        )
    )
    op.create_table(
        "run_result_attempts",
        sa.Column("attempt_id", sa.String(length=78), nullable=False),
        sa.Column("run_id", sa.String(length=128), nullable=False),
        sa.Column("result_kind", sa.String(length=16), nullable=False),
        sa.Column("artifact_generation", sa.String(length=76), nullable=True),
        sa.CheckConstraint(
            "(result_kind = 'artifact' AND artifact_generation IS NULL) OR "
            "(result_kind = 'qc' AND artifact_generation IS NOT NULL)",
            name="ck_run_result_attempts_binding",
        ),
        sa.ForeignKeyConstraint(["run_id"], ["runs.run_id"], ondelete="CASCADE"),
        sa.PrimaryKeyConstraint("attempt_id"),
    )
    op.create_index(
        "ix_run_result_attempts_run_id",
        "run_result_attempts",
        ["run_id"],
    )


def downgrade() -> None:
    op.drop_index("ix_run_result_attempts_run_id", table_name="run_result_attempts")
    op.drop_table("run_result_attempts")
    op.drop_table("run_result_states")
    with op.batch_alter_table("run_artifacts") as batch_op:
        batch_op.drop_column("revision")
