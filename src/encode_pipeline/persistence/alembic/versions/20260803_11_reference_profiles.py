"""Add operator Reference Profiles and immutable snapshot/run evidence.

Revision ID: 20260803_11
Revises: 20260726_10
"""

from __future__ import annotations

from alembic import op
import sqlalchemy as sa


revision = "20260803_11"
down_revision = "20260726_10"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.create_table(
        "reference_profiles",
        sa.Column("profile_id", sa.String(length=37), nullable=False),
        sa.Column("safe_key", sa.String(length=255), nullable=False),
        sa.Column("created_at", sa.DateTime(timezone=True), nullable=False),
        sa.Column("enabled_revision_id", sa.String(length=38), nullable=True),
        sa.CheckConstraint(
            "length(profile_id) = 37 AND substr(profile_id, 1, 5) = 'refp_'",
            name="ck_reference_profiles_id",
        ),
        sa.CheckConstraint(
            "length(trim(safe_key)) BETWEEN 1 AND 255",
            name="ck_reference_profiles_safe_key",
        ),
        sa.PrimaryKeyConstraint("profile_id"),
        sa.UniqueConstraint("safe_key", name="uq_reference_profiles_safe_key"),
        sa.UniqueConstraint(
            "profile_id",
            "enabled_revision_id",
            name="uq_reference_profiles_enabled_revision",
        ),
    )
    op.create_table(
        "reference_profile_revisions",
        sa.Column("revision_id", sa.String(length=38), nullable=False),
        sa.Column("profile_id", sa.String(length=37), nullable=False),
        sa.Column("revision_number", sa.Integer(), nullable=False),
        sa.Column("display_name", sa.String(length=255), nullable=False),
        sa.Column("organism", sa.String(length=255), nullable=False),
        sa.Column("assembly", sa.String(length=255), nullable=False),
        sa.Column("config_key", sa.String(length=255), nullable=False),
        sa.Column("public_identity_scheme", sa.String(length=64), nullable=False),
        sa.Column("public_identity_sha256", sa.String(length=64), nullable=False),
        sa.Column("created_at", sa.DateTime(timezone=True), nullable=False),
        sa.CheckConstraint(
            "length(revision_id) = 38 AND substr(revision_id, 1, 6) = 'refpr_'",
            name="ck_reference_profile_revisions_id",
        ),
        sa.CheckConstraint(
            "revision_number >= 1",
            name="ck_reference_profile_revisions_positive_number",
        ),
        sa.CheckConstraint(
            "length(trim(display_name)) BETWEEN 1 AND 255",
            name="ck_reference_profile_revisions_display_name",
        ),
        sa.CheckConstraint(
            "length(trim(organism)) BETWEEN 1 AND 255",
            name="ck_reference_profile_revisions_organism",
        ),
        sa.CheckConstraint(
            "length(trim(assembly)) BETWEEN 1 AND 255",
            name="ck_reference_profile_revisions_assembly",
        ),
        sa.CheckConstraint(
            "length(trim(config_key)) BETWEEN 1 AND 255",
            name="ck_reference_profile_revisions_config_key",
        ),
        sa.CheckConstraint(
            "public_identity_scheme = 'sha256-framed-reference-profile-revision-v1'",
            name="ck_reference_profile_revisions_identity_scheme",
        ),
        sa.CheckConstraint(
            "length(public_identity_sha256) = 64",
            name="ck_reference_profile_revisions_identity_length",
        ),
        sa.ForeignKeyConstraint(
            ["profile_id"],
            ["reference_profiles.profile_id"],
            name="fk_reference_profile_revisions_profile",
            ondelete="RESTRICT",
        ),
        sa.PrimaryKeyConstraint("revision_id"),
        sa.UniqueConstraint(
            "profile_id",
            "revision_id",
            name="uq_reference_profile_revisions_profile_revision",
        ),
        sa.UniqueConstraint(
            "profile_id",
            "revision_number",
            name="uq_reference_profile_revisions_profile_number",
        ),
        sa.UniqueConstraint(
            "profile_id",
            "revision_id",
            "public_identity_sha256",
            name="uq_reference_profile_revisions_identity",
        ),
    )
    with op.batch_alter_table("reference_profiles") as batch:
        batch.create_foreign_key(
            "fk_reference_profiles_enabled_revision",
            "reference_profile_revisions",
            ["profile_id", "enabled_revision_id"],
            ["profile_id", "revision_id"],
            ondelete="RESTRICT",
        )
    op.create_index(
        "ix_reference_profiles_enabled_created",
        "reference_profiles",
        ["enabled_revision_id", "created_at", "profile_id"],
    )
    op.create_index(
        "ix_reference_profile_revisions_profile_number",
        "reference_profile_revisions",
        ["profile_id", "revision_number"],
    )
    op.create_table(
        "reference_profile_workflow_bindings",
        sa.Column("profile_id", sa.String(length=37), nullable=False),
        sa.Column("revision_id", sa.String(length=38), nullable=False),
        sa.Column("workflow_id", sa.String(length=255), nullable=False),
        sa.Column("contract_version", sa.String(length=255), nullable=False),
        sa.Column("identity_scheme", sa.String(length=64), nullable=False),
        sa.Column("identity_sha256", sa.String(length=64), nullable=False),
        sa.CheckConstraint(
            "length(trim(workflow_id)) BETWEEN 1 AND 255",
            name="ck_reference_profile_workflow_bindings_workflow",
        ),
        sa.CheckConstraint(
            "length(trim(contract_version)) BETWEEN 1 AND 255",
            name="ck_reference_profile_workflow_bindings_contract",
        ),
        sa.CheckConstraint(
            "identity_scheme = 'sha256-framed-adapter-reference-binding-v1'",
            name="ck_reference_profile_workflow_bindings_identity_scheme",
        ),
        sa.CheckConstraint(
            "length(identity_sha256) = 64",
            name="ck_reference_profile_workflow_bindings_identity_length",
        ),
        sa.ForeignKeyConstraint(
            ["profile_id", "revision_id"],
            [
                "reference_profile_revisions.profile_id",
                "reference_profile_revisions.revision_id",
            ],
            name="fk_reference_profile_workflow_bindings_revision",
            ondelete="RESTRICT",
        ),
        sa.PrimaryKeyConstraint("profile_id", "revision_id", "workflow_id"),
        sa.UniqueConstraint(
            "profile_id",
            "revision_id",
            "workflow_id",
            "contract_version",
            "identity_sha256",
            name="uq_reference_profile_workflow_bindings_identity",
        ),
    )
    op.create_index(
        "ix_reference_profile_workflow_bindings_workflow",
        "reference_profile_workflow_bindings",
        ["workflow_id", "profile_id", "revision_id"],
    )
    _create_evidence_table("snapshot")
    _create_evidence_table("run")
    _create_lifecycle_triggers()


def downgrade() -> None:
    _drop_triggers()
    op.drop_index(
        "ix_run_reference_bindings_revision",
        table_name="run_reference_bindings",
    )
    op.drop_table("run_reference_bindings")
    op.drop_index(
        "ix_snapshot_reference_bindings_revision",
        table_name="snapshot_reference_bindings",
    )
    op.drop_table("snapshot_reference_bindings")
    op.drop_index(
        "ix_reference_profile_workflow_bindings_workflow",
        table_name="reference_profile_workflow_bindings",
    )
    op.drop_table("reference_profile_workflow_bindings")
    with op.batch_alter_table("reference_profiles") as batch:
        batch.drop_constraint(
            "fk_reference_profiles_enabled_revision",
            type_="foreignkey",
        )
    op.drop_index(
        "ix_reference_profile_revisions_profile_number",
        table_name="reference_profile_revisions",
    )
    op.drop_table("reference_profile_revisions")
    op.drop_index(
        "ix_reference_profiles_enabled_created",
        table_name="reference_profiles",
    )
    op.drop_table("reference_profiles")


def _create_evidence_table(kind: str) -> None:
    is_snapshot = kind == "snapshot"
    table_name = f"{kind}_reference_bindings"
    owner_column = "snapshot_id" if is_snapshot else "run_id"
    owner_table = "validated_input_snapshots" if is_snapshot else "runs"
    owner_length = 64 if is_snapshot else 128
    op.create_table(
        table_name,
        sa.Column(owner_column, sa.String(length=owner_length), nullable=False),
        sa.Column("profile_id", sa.String(length=37), nullable=False),
        sa.Column("revision_id", sa.String(length=38), nullable=False),
        sa.Column("workflow_id", sa.String(length=255), nullable=False),
        sa.Column(
            "revision_public_identity_scheme",
            sa.String(length=64),
            nullable=False,
        ),
        sa.Column(
            "revision_public_identity_sha256",
            sa.String(length=64),
            nullable=False,
        ),
        sa.Column("adapter_contract_version", sa.String(length=255), nullable=False),
        sa.Column("adapter_identity_scheme", sa.String(length=64), nullable=False),
        sa.Column("adapter_identity_sha256", sa.String(length=64), nullable=False),
        sa.Column("binding_digest_scheme", sa.String(length=64), nullable=False),
        sa.Column("binding_digest", sa.String(length=64), nullable=False),
        sa.Column("bound_at", sa.DateTime(timezone=True), nullable=False),
        sa.CheckConstraint(
            "revision_public_identity_scheme = "
            "'sha256-framed-reference-profile-revision-v1'",
            name=f"ck_{kind}_reference_bindings_revision_scheme",
        ),
        sa.CheckConstraint(
            "length(revision_public_identity_sha256) = 64",
            name=f"ck_{kind}_reference_bindings_revision_identity",
        ),
        sa.CheckConstraint(
            "adapter_identity_scheme = 'sha256-framed-adapter-reference-binding-v1'",
            name=f"ck_{kind}_reference_bindings_adapter_scheme",
        ),
        sa.CheckConstraint(
            "length(adapter_identity_sha256) = 64",
            name=f"ck_{kind}_reference_bindings_adapter_identity",
        ),
        sa.CheckConstraint(
            "binding_digest_scheme = 'sha256-framed-reference-profile-binding-v1'",
            name=f"ck_{kind}_reference_bindings_digest_scheme",
        ),
        sa.CheckConstraint(
            "length(binding_digest) = 64",
            name=f"ck_{kind}_reference_bindings_digest_length",
        ),
        sa.ForeignKeyConstraint(
            [owner_column],
            [f"{owner_table}.{owner_column}"],
            name=f"fk_{kind}_reference_bindings_{kind}",
            ondelete="RESTRICT" if is_snapshot else "CASCADE",
        ),
        sa.ForeignKeyConstraint(
            ["profile_id", "revision_id", "revision_public_identity_sha256"],
            [
                "reference_profile_revisions.profile_id",
                "reference_profile_revisions.revision_id",
                "reference_profile_revisions.public_identity_sha256",
            ],
            name=f"fk_{kind}_reference_bindings_revision",
            ondelete="RESTRICT",
        ),
        sa.ForeignKeyConstraint(
            [
                "profile_id",
                "revision_id",
                "workflow_id",
                "adapter_contract_version",
                "adapter_identity_sha256",
            ],
            [
                "reference_profile_workflow_bindings.profile_id",
                "reference_profile_workflow_bindings.revision_id",
                "reference_profile_workflow_bindings.workflow_id",
                "reference_profile_workflow_bindings.contract_version",
                "reference_profile_workflow_bindings.identity_sha256",
            ],
            name=f"fk_{kind}_reference_bindings_workflow_binding",
            ondelete="RESTRICT",
        ),
        sa.PrimaryKeyConstraint(owner_column),
    )
    op.create_index(
        f"ix_{kind}_reference_bindings_revision",
        table_name,
        ["revision_id", owner_column],
    )


def _create_lifecycle_triggers() -> None:
    op.execute(
        """
        CREATE TRIGGER trg_reference_profiles_no_delete
        BEFORE DELETE ON reference_profiles
        BEGIN
            SELECT RAISE(ABORT, 'ReferenceProfile rows are durable');
        END
        """
    )
    op.execute(
        """
        CREATE TRIGGER trg_reference_profiles_identity_immutable
        BEFORE UPDATE ON reference_profiles
        WHEN NEW.profile_id != OLD.profile_id
          OR NEW.safe_key != OLD.safe_key
          OR NEW.created_at != OLD.created_at
        BEGIN
            SELECT RAISE(ABORT, 'ReferenceProfile identity is immutable');
        END
        """
    )
    for table, label in (
        ("reference_profile_revisions", "ReferenceProfile revision"),
        ("reference_profile_workflow_bindings", "ReferenceProfile workflow binding"),
        ("snapshot_reference_bindings", "snapshot reference binding"),
        ("run_reference_bindings", "run reference binding"),
    ):
        op.execute(
            f"""
            CREATE TRIGGER trg_{table}_no_update
            BEFORE UPDATE ON {table}
            BEGIN
                SELECT RAISE(ABORT, '{label} rows are immutable');
            END
            """
        )
        op.execute(
            f"""
            CREATE TRIGGER trg_{table}_no_delete
            BEFORE DELETE ON {table}
            BEGIN
                SELECT RAISE(ABORT, '{label} rows are immutable');
            END
            """
        )
    op.execute(
        """
        CREATE TRIGGER trg_snapshot_reference_bindings_workflow
        BEFORE INSERT ON snapshot_reference_bindings
        WHEN NOT EXISTS (
            SELECT 1 FROM validated_input_snapshots
            WHERE snapshot_id = NEW.snapshot_id
              AND workflow_id = NEW.workflow_id
        )
        BEGIN
            SELECT RAISE(ABORT, 'snapshot reference workflow mismatch');
        END
        """
    )
    op.execute(
        """
        CREATE TRIGGER trg_run_reference_bindings_workflow
        BEFORE INSERT ON run_reference_bindings
        WHEN NOT EXISTS (
            SELECT 1 FROM runs
            WHERE run_id = NEW.run_id
              AND workflow_id = NEW.workflow_id
        )
        BEGIN
            SELECT RAISE(ABORT, 'run reference workflow mismatch');
        END
        """
    )


def _drop_triggers() -> None:
    for name in (
        "trg_run_reference_bindings_workflow",
        "trg_snapshot_reference_bindings_workflow",
        "trg_run_reference_bindings_no_delete",
        "trg_run_reference_bindings_no_update",
        "trg_snapshot_reference_bindings_no_delete",
        "trg_snapshot_reference_bindings_no_update",
        "trg_reference_profile_workflow_bindings_no_delete",
        "trg_reference_profile_workflow_bindings_no_update",
        "trg_reference_profile_revisions_no_delete",
        "trg_reference_profile_revisions_no_update",
        "trg_reference_profiles_identity_immutable",
        "trg_reference_profiles_no_delete",
    ):
        op.execute(f"DROP TRIGGER IF EXISTS {name}")
