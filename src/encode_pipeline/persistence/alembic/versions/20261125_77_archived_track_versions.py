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

    op.execute("""
        CREATE TRIGGER archive_version_experiment_insert BEFORE INSERT ON archived_track_versions
        WHEN NOT EXISTS (
            SELECT 1 FROM experiment_revisions WHERE experiment_revision_id = NEW.experiment_revision_id
                AND experiment_id = NEW.experiment_id AND project_id = NEW.project_id
        )
        BEGIN SELECT RAISE(ABORT, 'archive experiment scope mismatch'); END
    """)
    op.execute("""
        CREATE TRIGGER archive_version_fixed_update BEFORE UPDATE OF
            archive_version_id, project_id, experiment_id, experiment_revision_id, track_key,
            version_number, display_name, assembly, track_type, source_identity, source
        ON archived_track_versions
        BEGIN SELECT RAISE(ABORT, 'sealed archive version is immutable'); END
    """)
    # The only mutable fields are the physical file identities after an operator
    # verifies a byte-identical restore. Compare each member independent of JSON key order.
    protected = [
        f"json_extract(OLD.members, '$.{role}.{field}') IS NOT json_extract(NEW.members, '$.{role}.{field}')"
        for role in ("data", "index")
        for field in ("sha256", "size_bytes", "chunks")
    ]
    op.execute(
        "CREATE TRIGGER archive_version_fixed_bytes BEFORE UPDATE OF members ON archived_track_versions "
        "WHEN "
        + " OR ".join(protected)
        + " BEGIN SELECT RAISE(ABORT, 'sealed archive members are immutable'); END"
    )
    op.execute("""
        CREATE TRIGGER archive_version_fixed_delete BEFORE DELETE ON archived_track_versions
        BEGIN SELECT RAISE(ABORT, 'sealed archive version is immutable'); END
    """)


def downgrade() -> None:
    has_versions = (
        op.get_bind()
        .execute(sa.text("SELECT 1 FROM archived_track_versions LIMIT 1"))
        .scalar()
    )
    if has_versions:
        raise RuntimeError("Cannot downgrade a catalog with sealed track versions")
    op.drop_table("archived_track_versions")
