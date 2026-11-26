"""Install the archive guards on databases already stamped by an earlier 77 build."""

from alembic import op
import sqlalchemy as sa


revision = "20261125_77b"
down_revision = "20261125_77"
branch_labels = None
depends_on = None


def upgrade() -> None:
    # Both layouts have the same tables and revision stamp. Keep all historic
    # IDs and manifests in place and reject inconsistent pre-guard bindings.
    mismatch = op.get_bind().scalar(
        sa.text("""
        SELECT 1 FROM archived_track_versions AS v
        WHERE NOT EXISTS (
            SELECT 1 FROM experiment_revisions AS r
            WHERE r.experiment_revision_id = v.experiment_revision_id
              AND r.experiment_id = v.experiment_id AND r.project_id = v.project_id
        ) LIMIT 1
    """)
    )
    if mismatch:
        raise RuntimeError(
            "Cannot protect archived track versions with inconsistent experiment scope"
        )
    op.execute("""
        CREATE TRIGGER IF NOT EXISTS archive_version_experiment_insert
        BEFORE INSERT ON archived_track_versions
        WHEN NOT EXISTS (
            SELECT 1 FROM experiment_revisions WHERE experiment_revision_id = NEW.experiment_revision_id
                AND experiment_id = NEW.experiment_id AND project_id = NEW.project_id
        )
        BEGIN SELECT RAISE(ABORT, 'archive experiment scope mismatch'); END
    """)
    op.execute("""
        CREATE TRIGGER IF NOT EXISTS archive_version_fixed_update BEFORE UPDATE OF
            archive_version_id, project_id, experiment_id, experiment_revision_id, track_key,
            version_number, display_name, assembly, track_type, source_identity, source
        ON archived_track_versions
        BEGIN SELECT RAISE(ABORT, 'sealed archive version is immutable'); END
    """)
    protected = [
        f"json_extract(OLD.members, '$.{role}.{field}') IS NOT json_extract(NEW.members, '$.{role}.{field}')"
        for role in ("data", "index")
        for field in ("sha256", "size_bytes", "chunks")
    ]
    op.execute(
        "CREATE TRIGGER IF NOT EXISTS archive_version_fixed_bytes BEFORE UPDATE OF members ON archived_track_versions "
        "WHEN "
        + " OR ".join(protected)
        + " BEGIN SELECT RAISE(ABORT, 'sealed archive members are immutable'); END"
    )
    op.execute("""
        CREATE TRIGGER IF NOT EXISTS archive_version_fixed_delete BEFORE DELETE ON archived_track_versions
        BEGIN SELECT RAISE(ABORT, 'sealed archive version is immutable'); END
    """)


def downgrade() -> None:
    # Revision 77 already ships the guards for fresh installs. Keep them on an
    # empty catalog and preserve the current stamp if a downgrade risks history.
    if op.get_bind().scalar(sa.text("SELECT 1 FROM archived_track_versions LIMIT 1")):
        raise RuntimeError("Cannot downgrade a catalog with sealed track versions")
