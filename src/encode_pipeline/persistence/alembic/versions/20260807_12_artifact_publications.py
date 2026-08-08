"""Add append-only artifact publication provenance.

Revision ID: 20260807_12
Revises: 20260803_11
"""

from __future__ import annotations

from alembic import op
import sqlalchemy as sa


revision = "20260807_12"
down_revision = "20260803_11"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.create_table(
        "artifact_publications",
        sa.Column("run_id", sa.String(length=128), nullable=False),
        sa.Column("artifact_id", sa.String(length=128), nullable=False),
        sa.Column("artifact_generation", sa.String(length=76), nullable=False),
        sa.Column("artifact_revision", sa.String(length=76), nullable=False),
        sa.Column("output_type", sa.String(length=128), nullable=False),
        sa.Column("published_at", sa.DateTime(timezone=True), nullable=False),
        sa.CheckConstraint(
            "length(artifact_id) BETWEEN 1 AND 128 AND "
            "substr(artifact_id, 1, 1) GLOB '[A-Za-z]' AND "
            "artifact_id NOT GLOB '*[^A-Za-z0-9_.-]*'",
            name="ck_artifact_publications_artifact_id",
        ),
        sa.CheckConstraint(
            "length(artifact_generation) = 76 AND "
            "substr(artifact_generation, 1, 12) = 'artifactgen-' AND "
            "substr(artifact_generation, 13) NOT GLOB '*[^0-9a-f]*'",
            name="ck_artifact_publications_generation",
        ),
        sa.CheckConstraint(
            "length(artifact_revision) = 76 AND "
            "substr(artifact_revision, 1, 12) = 'artifactrev-' AND "
            "substr(artifact_revision, 13) NOT GLOB '*[^0-9a-f]*'",
            name="ck_artifact_publications_revision",
        ),
        sa.CheckConstraint(
            "length(output_type) BETWEEN 1 AND 128 AND "
            "substr(output_type, 1, 1) GLOB '[A-Za-z]' AND "
            "output_type NOT GLOB '*[^A-Za-z0-9_.-]*'",
            name="ck_artifact_publications_output_type",
        ),
        sa.ForeignKeyConstraint(
            ["run_id"],
            ["runs.run_id"],
            name="fk_artifact_publications_run",
            ondelete="RESTRICT",
        ),
        sa.PrimaryKeyConstraint(
            "run_id",
            "artifact_id",
            "artifact_generation",
        ),
    )
    op.create_index(
        "ix_artifact_publications_published",
        "artifact_publications",
        ["published_at", "run_id", "artifact_generation", "artifact_id"],
    )
    op.create_index(
        "ix_artifact_publications_run_published",
        "artifact_publications",
        ["run_id", "published_at", "artifact_generation", "artifact_id"],
    )
    op.create_index(
        "ix_artifact_publications_output_type_published",
        "artifact_publications",
        [
            "output_type",
            "published_at",
            "run_id",
            "artifact_generation",
            "artifact_id",
        ],
    )
    op.execute(
        """
        CREATE TRIGGER trg_artifact_publications_no_update
        BEFORE UPDATE ON artifact_publications
        BEGIN
            SELECT RAISE(ABORT, 'artifact publications are append-only');
        END
        """
    )
    op.execute(
        """
        CREATE TRIGGER trg_artifact_publications_no_delete
        BEFORE DELETE ON artifact_publications
        BEGIN
            SELECT RAISE(ABORT, 'artifact publications are append-only');
        END
        """
    )


def downgrade() -> None:
    op.execute("DROP TRIGGER IF EXISTS trg_artifact_publications_no_delete")
    op.execute("DROP TRIGGER IF EXISTS trg_artifact_publications_no_update")
    op.drop_index(
        "ix_artifact_publications_output_type_published",
        table_name="artifact_publications",
    )
    op.drop_index(
        "ix_artifact_publications_run_published",
        table_name="artifact_publications",
    )
    op.drop_index(
        "ix_artifact_publications_published",
        table_name="artifact_publications",
    )
    op.drop_table("artifact_publications")
