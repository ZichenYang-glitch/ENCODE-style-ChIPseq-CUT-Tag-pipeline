"""Add durable workflow-neutral QC metrics.

Revision ID: 20260712_05
Revises: 20260712_04
"""

from __future__ import annotations

from alembic import op
import sqlalchemy as sa


revision = "20260712_05"
down_revision = "20260712_04"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.create_table(
        "run_qc_metrics",
        sa.Column("id", sa.Integer(), autoincrement=True, nullable=False),
        sa.Column("metric_id", sa.String(length=128), nullable=False),
        sa.Column("run_id", sa.String(length=128), nullable=False),
        sa.Column("metric_key", sa.String(length=128), nullable=False),
        sa.Column("display_name", sa.String(length=255), nullable=False),
        sa.Column("value_text", sa.String(length=40), nullable=False),
        sa.Column("unit", sa.String(length=32), nullable=False),
        sa.Column("scope", sa.String(length=16), nullable=False),
        sa.Column("sample_id", sa.String(length=255), nullable=True),
        sa.Column("experiment_id", sa.String(length=255), nullable=True),
        sa.Column("assay", sa.String(length=255), nullable=True),
        sa.Column("qc_flag", sa.String(length=16), nullable=True),
        sa.Column("source_artifact_id", sa.String(length=128), nullable=False),
        sa.Column("produced_at", sa.DateTime(timezone=True), nullable=False),
        sa.CheckConstraint(
            "qc_flag IS NULL OR qc_flag IN ('pass', 'warning', 'fail')",
            name="ck_run_qc_metrics_flag",
        ),
        sa.CheckConstraint(
            "scope IN ('run', 'sample', 'experiment')",
            name="ck_run_qc_metrics_scope",
        ),
        sa.CheckConstraint(
            "(scope = 'run' AND sample_id IS NULL AND experiment_id IS NULL) OR "
            "(scope = 'sample' AND sample_id IS NOT NULL) OR "
            "(scope = 'experiment' AND sample_id IS NULL AND "
            "experiment_id IS NOT NULL)",
            name="ck_run_qc_metrics_scope_identifiers",
        ),
        sa.CheckConstraint(
            "length(value_text) BETWEEN 1 AND 40",
            name="ck_run_qc_metrics_value_text_length",
        ),
        sa.ForeignKeyConstraint(["run_id"], ["runs.run_id"], ondelete="CASCADE"),
        sa.PrimaryKeyConstraint("id"),
        sa.UniqueConstraint(
            "run_id",
            "metric_id",
            name="uq_run_qc_metrics_run_metric",
        ),
    )
    op.create_index("ix_run_qc_metrics_run_id", "run_qc_metrics", ["run_id"])


def downgrade() -> None:
    op.drop_index("ix_run_qc_metrics_run_id", table_name="run_qc_metrics")
    op.drop_table("run_qc_metrics")
