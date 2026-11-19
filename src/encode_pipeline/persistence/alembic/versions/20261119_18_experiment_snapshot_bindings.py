"""Add Slice B experiment bindings for snapshots and runs.

Revision ID: 20261119_18
Revises: 20261118_17

Upgrade creates empty binding tables. Existing snapshots and runs stay unbound.
"""

from __future__ import annotations

from alembic import op
import sqlalchemy as sa


revision = "20261119_18"
down_revision = "20261118_17"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.create_table(
        "snapshot_experiment_bindings",
        sa.Column("snapshot_id", sa.String(length=64), nullable=False),
        sa.Column("experiment_revision_id", sa.String(length=37), nullable=False),
        sa.Column("payload_digest", sa.String(length=64), nullable=False),
        sa.CheckConstraint(
            "length(snapshot_id) = 38 AND substr(snapshot_id, 1, 6) = 'vsnap_' "
            "AND substr(snapshot_id, 7) NOT GLOB '*[^0-9a-f]*'",
            name="ck_snapshot_experiment_bindings_snapshot_id",
        ),
        sa.CheckConstraint(
            "length(experiment_revision_id) = 37 "
            "AND substr(experiment_revision_id, 1, 5) = 'expr_' "
            "AND substr(experiment_revision_id, 6) NOT GLOB '*[^0-9a-f]*'",
            name="ck_snapshot_experiment_bindings_revision_id",
        ),
        sa.CheckConstraint(
            "length(payload_digest) = 64 AND payload_digest NOT GLOB '*[^0-9a-f]*'",
            name="ck_snapshot_experiment_bindings_digest",
        ),
        sa.ForeignKeyConstraint(
            ["snapshot_id"],
            ["validated_input_snapshots.snapshot_id"],
            name="fk_snapshot_experiment_bindings_snapshot",
            ondelete="RESTRICT",
        ),
        sa.ForeignKeyConstraint(
            ["experiment_revision_id"],
            ["experiment_revisions.experiment_revision_id"],
            name="fk_snapshot_experiment_bindings_revision",
            ondelete="RESTRICT",
        ),
        sa.PrimaryKeyConstraint("snapshot_id"),
    )
    op.create_index(
        "ix_snapshot_experiment_bindings_experiment",
        "snapshot_experiment_bindings",
        ["experiment_revision_id", "snapshot_id"],
    )
    op.create_table(
        "run_experiment_bindings",
        sa.Column("run_id", sa.String(length=128), nullable=False),
        sa.Column("experiment_revision_id", sa.String(length=37), nullable=False),
        sa.Column("payload_digest", sa.String(length=64), nullable=False),
        sa.CheckConstraint(
            "length(run_id) BETWEEN 1 AND 128",
            name="ck_run_experiment_bindings_run_id",
        ),
        sa.CheckConstraint(
            "length(experiment_revision_id) = 37 "
            "AND substr(experiment_revision_id, 1, 5) = 'expr_' "
            "AND substr(experiment_revision_id, 6) NOT GLOB '*[^0-9a-f]*'",
            name="ck_run_experiment_bindings_revision_id",
        ),
        sa.CheckConstraint(
            "length(payload_digest) = 64 AND payload_digest NOT GLOB '*[^0-9a-f]*'",
            name="ck_run_experiment_bindings_digest",
        ),
        sa.ForeignKeyConstraint(
            ["run_id"],
            ["runs.run_id"],
            name="fk_run_experiment_bindings_run",
            ondelete="RESTRICT",
        ),
        sa.ForeignKeyConstraint(
            ["experiment_revision_id"],
            ["experiment_revisions.experiment_revision_id"],
            name="fk_run_experiment_bindings_revision",
            ondelete="RESTRICT",
        ),
        sa.PrimaryKeyConstraint("run_id"),
    )
    op.create_index(
        "ix_run_experiment_bindings_experiment",
        "run_experiment_bindings",
        ["experiment_revision_id", "run_id"],
    )
    op.execute(
        """
        CREATE TRIGGER trg_snapshot_experiment_bindings_no_update
        BEFORE UPDATE ON snapshot_experiment_bindings
        BEGIN
            SELECT RAISE(ABORT, 'snapshot experiment binding is immutable');
        END
        """
    )
    op.execute(
        """
        CREATE TRIGGER trg_snapshot_experiment_bindings_no_delete
        BEFORE DELETE ON snapshot_experiment_bindings
        BEGIN
            SELECT RAISE(ABORT, 'snapshot experiment binding is immutable');
        END
        """
    )
    op.execute(
        """
        CREATE TRIGGER trg_run_experiment_bindings_no_update
        BEFORE UPDATE ON run_experiment_bindings
        BEGIN
            SELECT RAISE(ABORT, 'run experiment binding is immutable');
        END
        """
    )
    op.execute(
        """
        CREATE TRIGGER trg_run_experiment_bindings_no_delete
        BEFORE DELETE ON run_experiment_bindings
        BEGIN
            SELECT RAISE(ABORT, 'run experiment binding is immutable');
        END
        """
    )


def downgrade() -> None:
    for name in (
        "trg_run_experiment_bindings_no_delete",
        "trg_run_experiment_bindings_no_update",
        "trg_snapshot_experiment_bindings_no_delete",
        "trg_snapshot_experiment_bindings_no_update",
    ):
        op.execute(f"DROP TRIGGER IF EXISTS {name}")
    op.drop_index(
        "ix_run_experiment_bindings_experiment",
        table_name="run_experiment_bindings",
    )
    op.drop_table("run_experiment_bindings")
    op.drop_index(
        "ix_snapshot_experiment_bindings_experiment",
        table_name="snapshot_experiment_bindings",
    )
    op.drop_table("snapshot_experiment_bindings")
