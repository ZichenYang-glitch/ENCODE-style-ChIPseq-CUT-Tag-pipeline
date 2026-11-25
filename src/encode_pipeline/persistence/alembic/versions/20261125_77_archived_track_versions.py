"""Pin a complete, immutable archived track version to its experiment."""

from alembic import op
import sqlalchemy as sa


revision = "20261125_77"
down_revision = "20261119_39"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.create_table(
        "archived_track_versions",
        sa.Column("archive_version_id", sa.String(37), primary_key=True),
        sa.Column("project_id", sa.String(36), nullable=False),
        sa.Column("experiment_id", sa.String(36), nullable=False),
        sa.Column("experiment_revision_id", sa.String(37), nullable=False),
        sa.Column("track_key", sa.String(255), nullable=False),
        sa.Column("version_number", sa.Integer(), nullable=False),
        sa.Column("display_name", sa.String(255), nullable=False),
        sa.Column("assembly", sa.String(64), nullable=False),
        sa.Column("track_type", sa.String(16), nullable=False),
        sa.Column("source_identity", sa.String(64), nullable=False),
        sa.Column("source", sa.JSON(), nullable=False),
        sa.Column("members", sa.JSON(), nullable=False),
        sa.ForeignKeyConstraint(
            ["experiment_id", "project_id"],
            ["experiments.experiment_id", "experiments.project_id"],
            ondelete="RESTRICT",
        ),
        sa.ForeignKeyConstraint(
            ["experiment_revision_id"],
            ["experiment_revisions.experiment_revision_id"],
            ondelete="RESTRICT",
        ),
        sa.UniqueConstraint("experiment_id", "track_key", "version_number"),
        sa.UniqueConstraint("experiment_id", "track_key", "source_identity"),
        sa.CheckConstraint("version_number >= 1"),
        sa.CheckConstraint("track_type IN ('bigwig','bed','longrange')"),
    )


def downgrade() -> None:
    has_versions = (
        op.get_bind()
        .execute(sa.text("SELECT 1 FROM archived_track_versions LIMIT 1"))
        .scalar()
    )
    if has_versions:
        raise RuntimeError("Cannot downgrade a catalog with sealed track versions")
    op.drop_table("archived_track_versions")
