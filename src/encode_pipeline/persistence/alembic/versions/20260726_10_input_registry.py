"""Add approved storage, immutable input revisions, and input-use bindings.

Revision ID: 20260726_10
Revises: 20260726_09
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Mapping, Sequence
from hashlib import sha256
import json
from typing import Any

from alembic import op
import sqlalchemy as sa


revision = "20260726_10"
down_revision = "20260726_09"
branch_labels = None
depends_on = None


LEGACY_PROJECT_ID = "prj_00000000000000000000000000000000"
_PROJECT_SAMPLE_DIGEST_SCHEME = "sha256-framed-project-sample-binding-v1"
_INPUT_BINDING_DIGEST_SCHEME = "sha256-framed-input-use-binding-envelope-v1"
_WORKFLOW_INPUTS_DIGEST_SCHEME = "sha256-framed-workflow-inputs-v1"
_CANONICAL_JSON_SEPARATORS = (",", ":")
_PREFLIGHT_ERROR = (
    "inconsistent project/sample binding evidence blocks input registry upgrade"
)


def upgrade() -> None:
    connection = op.get_bind()
    snapshot_rows = (
        connection.execute(
            sa.text(
                "SELECT snapshot_id, workflow_id, adapter_version, payload_digest, "
                "payload_digest_scheme, canonical_payload, validated_at, "
                "consumed_run_id "
                "FROM validated_input_snapshots ORDER BY snapshot_id"
            )
        )
        .mappings()
        .all()
    )
    run_rows = (
        connection.execute(
            sa.text(
                "SELECT run_id, workflow_id, inputs, created_at "
                "FROM runs ORDER BY run_id"
            )
        )
        .mappings()
        .all()
    )
    snapshot_project_rows = (
        connection.execute(
            sa.text(
                "SELECT snapshot_id, project_id, binding_mode, provenance, "
                "workflow_inputs_digest, binding_digest_scheme, binding_digest "
                "FROM snapshot_project_bindings ORDER BY snapshot_id"
            )
        )
        .mappings()
        .all()
    )
    snapshot_sample_rows = (
        connection.execute(
            sa.text(
                "SELECT snapshot_id, project_id, sample_revision_id, "
                "payload_digest, ordinal FROM snapshot_sample_revisions "
                "ORDER BY snapshot_id, ordinal"
            )
        )
        .mappings()
        .all()
    )
    run_project_rows = (
        connection.execute(
            sa.text(
                "SELECT run_id, project_id, binding_mode, provenance, "
                "workflow_inputs_digest, binding_digest_scheme, binding_digest "
                "FROM run_project_bindings ORDER BY run_id"
            )
        )
        .mappings()
        .all()
    )
    run_sample_rows = (
        connection.execute(
            sa.text(
                "SELECT run_id, project_id, sample_revision_id, payload_digest, "
                "ordinal FROM run_samples ORDER BY run_id, ordinal"
            )
        )
        .mappings()
        .all()
    )
    snapshot_plan, run_plan = _preflight_and_plan_bindings(
        snapshot_rows=snapshot_rows,
        run_rows=run_rows,
        snapshot_project_rows=snapshot_project_rows,
        snapshot_sample_rows=snapshot_sample_rows,
        run_project_rows=run_project_rows,
        run_sample_rows=run_sample_rows,
    )

    _create_registry_tables()
    _create_input_binding_tables()
    _backfill_snapshot_input_bindings(connection, snapshot_rows, snapshot_plan)
    _backfill_run_input_bindings(connection, run_rows, run_plan)


def downgrade() -> None:
    op.execute("DROP TRIGGER IF EXISTS trg_run_input_members_managed_use")
    op.execute("DROP TRIGGER IF EXISTS trg_snapshot_input_members_managed_use")
    op.execute("DROP TRIGGER IF EXISTS trg_run_input_uses_declared_binding")
    op.execute("DROP TRIGGER IF EXISTS trg_snapshot_input_uses_declared_binding")
    _drop_immutable_triggers(
        "run_input_members",
        "run_input_uses",
        "run_input_bindings",
        "snapshot_input_members",
        "snapshot_input_uses",
        "snapshot_input_bindings",
        "input_file_revisions",
    )
    op.drop_index("ix_run_input_members_revision", table_name="run_input_members")
    op.drop_table("run_input_members")
    op.drop_index("ix_run_input_uses_key", table_name="run_input_uses")
    op.drop_table("run_input_uses")
    op.drop_index("ix_run_input_bindings_project", table_name="run_input_bindings")
    op.drop_table("run_input_bindings")
    op.drop_index(
        "ix_snapshot_input_members_revision",
        table_name="snapshot_input_members",
    )
    op.drop_table("snapshot_input_members")
    op.drop_index("ix_snapshot_input_uses_key", table_name="snapshot_input_uses")
    op.drop_table("snapshot_input_uses")
    op.drop_index(
        "ix_snapshot_input_bindings_project",
        table_name="snapshot_input_bindings",
    )
    op.drop_table("snapshot_input_bindings")
    op.drop_index(
        "uq_run_project_bindings_input_evidence",
        table_name="run_project_bindings",
    )
    op.drop_index(
        "uq_snapshot_project_bindings_input_evidence",
        table_name="snapshot_project_bindings",
    )
    op.drop_index(
        "ix_input_file_revisions_project_created",
        table_name="input_file_revisions",
    )
    op.drop_index(
        "ix_input_file_revisions_input_number",
        table_name="input_file_revisions",
    )
    op.drop_table("input_file_revisions")
    _drop_archive_only_triggers(
        "input_files",
        "storage_pools",
    )
    _drop_immutable_triggers("project_storage_pool_bindings")
    op.drop_index("ix_input_files_project_created", table_name="input_files")
    op.drop_table("input_files")
    op.drop_index(
        "ix_project_storage_pool_bindings_pool",
        table_name="project_storage_pool_bindings",
    )
    op.drop_table("project_storage_pool_bindings")
    op.drop_index(
        "ix_storage_pools_archived_created",
        table_name="storage_pools",
    )
    op.drop_table("storage_pools")


def _create_registry_tables() -> None:
    op.create_table(
        "storage_pools",
        sa.Column("storage_pool_id", sa.String(length=37), nullable=False),
        sa.Column("config_key", sa.String(length=255), nullable=False),
        sa.Column("display_name", sa.String(length=255), nullable=False),
        sa.Column("created_at", sa.DateTime(timezone=True), nullable=False),
        sa.Column("archived_at", sa.DateTime(timezone=True), nullable=True),
        sa.CheckConstraint(
            "length(storage_pool_id) = 37 AND substr(storage_pool_id, 1, 5) = 'stgp_'",
            name="ck_storage_pools_id",
        ),
        sa.CheckConstraint(
            "length(trim(config_key)) BETWEEN 1 AND 255",
            name="ck_storage_pools_config_key",
        ),
        sa.CheckConstraint(
            "length(trim(display_name)) BETWEEN 1 AND 255",
            name="ck_storage_pools_display_name",
        ),
        sa.CheckConstraint(
            "archived_at IS NULL OR archived_at >= created_at",
            name="ck_storage_pools_archive_order",
        ),
        sa.PrimaryKeyConstraint("storage_pool_id"),
        sa.UniqueConstraint("config_key", name="uq_storage_pools_config_key"),
    )
    op.create_index(
        "ix_storage_pools_archived_created",
        "storage_pools",
        ["archived_at", "created_at", "storage_pool_id"],
    )
    _create_archive_only_triggers(
        table_name="storage_pools",
        immutable_columns=(
            "storage_pool_id",
            "config_key",
            "display_name",
            "created_at",
        ),
        label="StoragePool",
    )

    op.create_table(
        "project_storage_pool_bindings",
        sa.Column("project_id", sa.String(length=36), nullable=False),
        sa.Column("storage_pool_id", sa.String(length=37), nullable=False),
        sa.Column("bound_at", sa.DateTime(timezone=True), nullable=False),
        sa.CheckConstraint(
            f"project_id != '{LEGACY_PROJECT_ID}'",
            name="ck_project_storage_pool_bindings_not_legacy",
        ),
        sa.ForeignKeyConstraint(
            ["project_id"],
            ["projects.project_id"],
            name="fk_project_storage_pool_bindings_project",
            ondelete="RESTRICT",
        ),
        sa.ForeignKeyConstraint(
            ["storage_pool_id"],
            ["storage_pools.storage_pool_id"],
            name="fk_project_storage_pool_bindings_pool",
            ondelete="RESTRICT",
        ),
        sa.PrimaryKeyConstraint("project_id"),
        sa.UniqueConstraint(
            "project_id",
            "storage_pool_id",
            name="uq_project_storage_pool_bindings_project_pool",
        ),
    )
    op.create_index(
        "ix_project_storage_pool_bindings_pool",
        "project_storage_pool_bindings",
        ["storage_pool_id", "project_id"],
    )
    _create_immutable_triggers(
        "project_storage_pool_bindings",
        "Project/StoragePool binding",
    )

    op.create_table(
        "input_files",
        sa.Column("input_file_id", sa.String(length=37), nullable=False),
        sa.Column("project_id", sa.String(length=36), nullable=False),
        sa.Column("storage_pool_id", sa.String(length=37), nullable=False),
        sa.Column("stable_key", sa.String(length=255), nullable=False),
        sa.Column("created_at", sa.DateTime(timezone=True), nullable=False),
        sa.Column("archived_at", sa.DateTime(timezone=True), nullable=True),
        sa.CheckConstraint(
            "length(input_file_id) = 37 AND substr(input_file_id, 1, 5) = 'inpf_'",
            name="ck_input_files_id",
        ),
        sa.CheckConstraint(
            "length(trim(stable_key)) BETWEEN 1 AND 255",
            name="ck_input_files_stable_key",
        ),
        sa.CheckConstraint(
            "archived_at IS NULL OR archived_at >= created_at",
            name="ck_input_files_archive_order",
        ),
        sa.ForeignKeyConstraint(
            ["project_id", "storage_pool_id"],
            [
                "project_storage_pool_bindings.project_id",
                "project_storage_pool_bindings.storage_pool_id",
            ],
            name="fk_input_files_project_pool",
            ondelete="RESTRICT",
        ),
        sa.PrimaryKeyConstraint("input_file_id"),
        sa.UniqueConstraint(
            "project_id",
            "stable_key",
            name="uq_input_files_project_stable_key",
        ),
        sa.UniqueConstraint(
            "project_id",
            "storage_pool_id",
            "input_file_id",
            name="uq_input_files_project_pool_file",
        ),
    )
    op.create_index(
        "ix_input_files_project_created",
        "input_files",
        ["project_id", "archived_at", "created_at", "input_file_id"],
    )
    _create_archive_only_triggers(
        table_name="input_files",
        immutable_columns=(
            "input_file_id",
            "project_id",
            "storage_pool_id",
            "stable_key",
            "created_at",
        ),
        label="InputFile",
    )

    op.create_table(
        "input_file_revisions",
        sa.Column("input_file_revision_id", sa.String(length=38), nullable=False),
        sa.Column("input_file_id", sa.String(length=37), nullable=False),
        sa.Column("project_id", sa.String(length=36), nullable=False),
        sa.Column("storage_pool_id", sa.String(length=37), nullable=False),
        sa.Column("revision_number", sa.Integer(), nullable=False),
        sa.Column("relative_path", sa.Text(), nullable=False),
        sa.Column("size_bytes", sa.Integer(), nullable=False),
        sa.Column("content_sha256", sa.String(length=64), nullable=False),
        sa.Column("digest_scheme", sa.String(length=64), nullable=False),
        sa.Column("digest", sa.String(length=64), nullable=False),
        sa.Column("created_at", sa.DateTime(timezone=True), nullable=False),
        sa.CheckConstraint(
            "length(input_file_revision_id) = 38 "
            "AND substr(input_file_revision_id, 1, 6) = 'inpfr_'",
            name="ck_input_file_revisions_id",
        ),
        sa.CheckConstraint(
            "revision_number >= 1",
            name="ck_input_file_revisions_positive_revision",
        ),
        sa.CheckConstraint(
            "length(relative_path) > 0 AND substr(relative_path, 1, 1) != '/' "
            "AND relative_path != '.' AND relative_path != '..' "
            "AND relative_path NOT LIKE '../%' "
            "AND relative_path NOT LIKE '%/../%' "
            "AND relative_path NOT LIKE '%/..' "
            "AND relative_path NOT LIKE './%' "
            "AND relative_path NOT LIKE '%/./%' "
            "AND relative_path NOT LIKE '%/.' "
            "AND relative_path NOT LIKE '%//%' "
            "AND instr(relative_path, char(92)) = 0",
            name="ck_input_file_revisions_safe_relative_path",
        ),
        sa.CheckConstraint(
            "size_bytes >= 0",
            name="ck_input_file_revisions_nonnegative_size",
        ),
        sa.CheckConstraint(
            "length(content_sha256) = 64",
            name="ck_input_file_revisions_content_sha256_length",
        ),
        sa.CheckConstraint(
            "digest_scheme = 'sha256-framed-input-file-revision-v1'",
            name="ck_input_file_revisions_digest_scheme",
        ),
        sa.CheckConstraint(
            "length(digest) = 64",
            name="ck_input_file_revisions_digest_length",
        ),
        sa.ForeignKeyConstraint(
            ["project_id", "storage_pool_id", "input_file_id"],
            [
                "input_files.project_id",
                "input_files.storage_pool_id",
                "input_files.input_file_id",
            ],
            name="fk_input_file_revisions_file",
            ondelete="RESTRICT",
        ),
        sa.PrimaryKeyConstraint("input_file_revision_id"),
        sa.UniqueConstraint(
            "input_file_id",
            "revision_number",
            name="uq_input_file_revisions_file_number",
        ),
        sa.UniqueConstraint(
            "project_id",
            "input_file_revision_id",
            "digest",
            name="uq_input_file_revisions_project_revision_digest",
        ),
        sa.UniqueConstraint(
            "project_id",
            "input_file_id",
            "input_file_revision_id",
            "digest",
            "size_bytes",
            "content_sha256",
            name="uq_input_file_revisions_binding_evidence",
        ),
    )
    op.create_index(
        "ix_input_file_revisions_input_number",
        "input_file_revisions",
        ["input_file_id", "revision_number"],
    )
    op.create_index(
        "ix_input_file_revisions_project_created",
        "input_file_revisions",
        ["project_id", "created_at", "input_file_revision_id"],
    )
    _create_immutable_triggers("input_file_revisions", "InputFileRevision")


def _create_input_binding_tables() -> None:
    op.create_index(
        "uq_snapshot_project_bindings_input_evidence",
        "snapshot_project_bindings",
        ["snapshot_id", "project_id", "binding_digest"],
        unique=True,
    )
    op.create_index(
        "uq_run_project_bindings_input_evidence",
        "run_project_bindings",
        ["run_id", "project_id", "binding_digest"],
        unique=True,
    )
    _create_binding_envelope_table(owner="snapshot")
    _create_input_use_table(owner="snapshot")
    _create_input_member_table(owner="snapshot")
    _create_binding_envelope_table(owner="run")
    _create_input_use_table(owner="run")
    _create_input_member_table(owner="run")


def _create_binding_envelope_table(*, owner: str) -> None:
    owner_column = "snapshot_id" if owner == "snapshot" else "run_id"
    owner_length = 64 if owner == "snapshot" else 128
    table_name = f"{owner}_input_bindings"
    parent_table = (
        "snapshot_project_bindings" if owner == "snapshot" else "run_project_bindings"
    )
    parent_delete = "RESTRICT" if owner == "snapshot" else "CASCADE"
    op.create_table(
        table_name,
        sa.Column(owner_column, sa.String(length=owner_length), nullable=False),
        sa.Column("project_id", sa.String(length=36), nullable=False),
        sa.Column("workflow_id", sa.Text(), nullable=False),
        sa.Column(
            "adapter_contract_version",
            sa.String(length=255),
            nullable=True,
        ),
        sa.Column("binding_mode", sa.String(length=40), nullable=False),
        sa.Column("workflow_inputs_digest", sa.String(length=64), nullable=False),
        sa.Column(
            "project_sample_binding_digest",
            sa.String(length=64),
            nullable=False,
        ),
        sa.Column("binding_digest_scheme", sa.String(length=64), nullable=False),
        sa.Column("binding_digest", sa.String(length=64), nullable=False),
        sa.Column("created_at", sa.DateTime(timezone=True), nullable=False),
        sa.CheckConstraint(
            "length(trim(workflow_id)) >= 1",
            name=f"ck_{table_name}_workflow_id",
        ),
        sa.CheckConstraint(
            "adapter_contract_version IS NULL "
            "OR length(trim(adapter_contract_version)) BETWEEN 1 AND 255",
            name=f"ck_{table_name}_adapter_contract_version",
        ),
        sa.CheckConstraint(
            "binding_mode IN ('compatibility_unresolved_v1', 'declared_input_uses_v1')",
            name=f"ck_{table_name}_mode",
        ),
        sa.CheckConstraint(
            "binding_mode = 'compatibility_unresolved_v1' "
            "OR adapter_contract_version IS NOT NULL",
            name=f"ck_{table_name}_declared_adapter_version",
        ),
        sa.CheckConstraint(
            "length(workflow_inputs_digest) = 64",
            name=f"ck_{table_name}_workflow_digest_length",
        ),
        sa.CheckConstraint(
            "length(project_sample_binding_digest) = 64",
            name=f"ck_{table_name}_project_sample_digest_length",
        ),
        sa.CheckConstraint(
            f"binding_digest_scheme = '{_INPUT_BINDING_DIGEST_SCHEME}'",
            name=f"ck_{table_name}_digest_scheme",
        ),
        sa.CheckConstraint(
            "length(binding_digest) = 64",
            name=f"ck_{table_name}_digest_length",
        ),
        sa.ForeignKeyConstraint(
            [
                owner_column,
                "project_id",
                "project_sample_binding_digest",
            ],
            [
                f"{parent_table}.{owner_column}",
                f"{parent_table}.project_id",
                f"{parent_table}.binding_digest",
            ],
            name=f"fk_{table_name}_project_binding",
            ondelete=parent_delete,
        ),
        sa.PrimaryKeyConstraint(owner_column),
        sa.UniqueConstraint(
            owner_column,
            "project_id",
            name=f"uq_{table_name}_project",
        ),
    )
    op.create_index(
        f"ix_{table_name}_project",
        table_name,
        ["project_id", "created_at", owner_column],
    )
    _create_immutable_triggers(
        table_name,
        f"{owner.title()} input binding",
        cascade_delete_root=("runs", "run_id") if owner == "run" else None,
    )


def _create_input_use_table(*, owner: str) -> None:
    owner_column = "snapshot_id" if owner == "snapshot" else "run_id"
    owner_length = 64 if owner == "snapshot" else 128
    table_name = f"{owner}_input_uses"
    binding_table = f"{owner}_input_bindings"
    parent_delete = "RESTRICT" if owner == "snapshot" else "CASCADE"
    op.create_table(
        table_name,
        sa.Column(owner_column, sa.String(length=owner_length), nullable=False),
        sa.Column("project_id", sa.String(length=36), nullable=False),
        sa.Column("ordinal", sa.Integer(), nullable=False),
        sa.Column("input_use_key", sa.String(length=255), nullable=False),
        sa.Column("occurrence", sa.Integer(), nullable=False),
        sa.Column("capability_version", sa.String(length=255), nullable=False),
        sa.Column(
            "closure_contract_version",
            sa.String(length=255),
            nullable=False,
        ),
        sa.Column("provenance_mode", sa.String(length=40), nullable=False),
        sa.Column("closure_digest_scheme", sa.String(length=64), nullable=True),
        sa.Column("closure_digest", sa.String(length=64), nullable=True),
        sa.CheckConstraint(
            "ordinal >= 0",
            name=f"ck_{table_name}_ordinal",
        ),
        sa.CheckConstraint(
            "length(trim(input_use_key)) BETWEEN 1 AND 255",
            name=f"ck_{table_name}_key",
        ),
        sa.CheckConstraint(
            "occurrence >= 0",
            name=f"ck_{table_name}_occurrence",
        ),
        sa.CheckConstraint(
            "length(trim(capability_version)) BETWEEN 1 AND 255",
            name=f"ck_{table_name}_capability_version",
        ),
        sa.CheckConstraint(
            "length(trim(closure_contract_version)) BETWEEN 1 AND 255",
            name=f"ck_{table_name}_closure_version",
        ),
        sa.CheckConstraint(
            "(provenance_mode = 'transitional_unmanaged_v1' "
            "AND closure_digest_scheme IS NULL AND closure_digest IS NULL) OR "
            "(provenance_mode = 'managed_revision_v1' "
            "AND closure_digest_scheme = 'sha256-framed-input-closure-v1' "
            "AND length(closure_digest) = 64)",
            name=f"ck_{table_name}_provenance_evidence",
        ),
        sa.ForeignKeyConstraint(
            [owner_column, "project_id"],
            [f"{binding_table}.{owner_column}", f"{binding_table}.project_id"],
            name=f"fk_{table_name}_binding",
            ondelete=parent_delete,
        ),
        sa.PrimaryKeyConstraint(owner_column, "ordinal"),
        sa.UniqueConstraint(
            owner_column,
            "input_use_key",
            "occurrence",
            name=f"uq_{table_name}_identity",
        ),
        sa.UniqueConstraint(
            owner_column,
            "ordinal",
            "project_id",
            name=f"uq_{table_name}_project",
        ),
    )
    op.create_index(
        f"ix_{table_name}_key",
        table_name,
        ["input_use_key", "capability_version", owner_column],
    )
    _create_immutable_triggers(
        table_name,
        f"{owner.title()} input use",
        cascade_delete_root=("runs", "run_id") if owner == "run" else None,
    )
    op.execute(
        f"""
        CREATE TRIGGER trg_{table_name}_declared_binding
        BEFORE INSERT ON {table_name}
        WHEN COALESCE(
            (
                SELECT binding_mode
                FROM {binding_table}
                WHERE {owner_column} = NEW.{owner_column}
                  AND project_id = NEW.project_id
            ),
            ''
        ) != 'declared_input_uses_v1'
        BEGIN
            SELECT RAISE(
                ABORT,
                'input uses require declared_input_uses_v1 binding'
            );
        END
        """
    )


def _create_input_member_table(*, owner: str) -> None:
    owner_column = "snapshot_id" if owner == "snapshot" else "run_id"
    owner_length = 64 if owner == "snapshot" else 128
    table_name = f"{owner}_input_members"
    use_table = f"{owner}_input_uses"
    parent_delete = "RESTRICT" if owner == "snapshot" else "CASCADE"
    op.create_table(
        table_name,
        sa.Column(owner_column, sa.String(length=owner_length), nullable=False),
        sa.Column("project_id", sa.String(length=36), nullable=False),
        sa.Column("use_ordinal", sa.Integer(), nullable=False),
        sa.Column("member_ordinal", sa.Integer(), nullable=False),
        sa.Column("logical_member_key", sa.String(length=255), nullable=False),
        sa.Column("input_file_id", sa.String(length=37), nullable=False),
        sa.Column("input_file_revision_id", sa.String(length=38), nullable=False),
        sa.Column("revision_digest", sa.String(length=64), nullable=False),
        sa.Column("size_bytes", sa.Integer(), nullable=False),
        sa.Column("content_sha256", sa.String(length=64), nullable=False),
        sa.CheckConstraint(
            "use_ordinal >= 0 AND member_ordinal >= 0",
            name=f"ck_{table_name}_ordinals",
        ),
        sa.CheckConstraint(
            "length(trim(logical_member_key)) BETWEEN 1 AND 255",
            name=f"ck_{table_name}_member_key",
        ),
        sa.CheckConstraint(
            "length(revision_digest) = 64",
            name=f"ck_{table_name}_revision_digest_length",
        ),
        sa.CheckConstraint(
            "size_bytes >= 0",
            name=f"ck_{table_name}_nonnegative_size",
        ),
        sa.CheckConstraint(
            "length(content_sha256) = 64",
            name=f"ck_{table_name}_content_sha256_length",
        ),
        sa.ForeignKeyConstraint(
            [owner_column, "use_ordinal", "project_id"],
            [
                f"{use_table}.{owner_column}",
                f"{use_table}.ordinal",
                f"{use_table}.project_id",
            ],
            name=f"fk_{table_name}_use",
            ondelete=parent_delete,
        ),
        sa.ForeignKeyConstraint(
            [
                "project_id",
                "input_file_id",
                "input_file_revision_id",
                "revision_digest",
                "size_bytes",
                "content_sha256",
            ],
            [
                "input_file_revisions.project_id",
                "input_file_revisions.input_file_id",
                "input_file_revisions.input_file_revision_id",
                "input_file_revisions.digest",
                "input_file_revisions.size_bytes",
                "input_file_revisions.content_sha256",
            ],
            name=f"fk_{table_name}_revision",
            ondelete="RESTRICT",
        ),
        sa.PrimaryKeyConstraint(owner_column, "use_ordinal", "member_ordinal"),
        sa.UniqueConstraint(
            owner_column,
            "use_ordinal",
            "logical_member_key",
            name=f"uq_{table_name}_member_key",
        ),
        sa.UniqueConstraint(
            owner_column,
            "use_ordinal",
            "input_file_revision_id",
            name=f"uq_{table_name}_revision",
        ),
    )
    op.create_index(
        f"ix_{table_name}_revision",
        table_name,
        ["input_file_revision_id", owner_column],
    )
    _create_immutable_triggers(
        table_name,
        f"{owner.title()} input member",
        cascade_delete_root=("runs", "run_id") if owner == "run" else None,
    )
    op.execute(
        f"""
        CREATE TRIGGER trg_{table_name}_managed_use
        BEFORE INSERT ON {table_name}
        WHEN COALESCE(
            (
                SELECT provenance_mode
                FROM {use_table}
                WHERE {owner_column} = NEW.{owner_column}
                  AND ordinal = NEW.use_ordinal
                  AND project_id = NEW.project_id
            ),
            ''
        ) != 'managed_revision_v1'
        BEGIN
            SELECT RAISE(
                ABORT,
                'input members require managed_revision_v1 provenance'
            );
        END
        """
    )


def _preflight_and_plan_bindings(
    *,
    snapshot_rows: Sequence[sa.RowMapping],
    run_rows: Sequence[sa.RowMapping],
    snapshot_project_rows: Sequence[sa.RowMapping],
    snapshot_sample_rows: Sequence[sa.RowMapping],
    run_project_rows: Sequence[sa.RowMapping],
    run_sample_rows: Sequence[sa.RowMapping],
) -> tuple[dict[str, dict[str, Any]], dict[str, dict[str, Any]]]:
    snapshot_projects = _one_row_per_owner(
        snapshot_project_rows,
        owner_column="snapshot_id",
    )
    run_projects = _one_row_per_owner(run_project_rows, owner_column="run_id")
    snapshot_samples = _rows_per_owner(
        snapshot_sample_rows,
        owner_column="snapshot_id",
    )
    run_samples = _rows_per_owner(run_sample_rows, owner_column="run_id")
    snapshots_by_id = {str(row["snapshot_id"]): row for row in snapshot_rows}
    runs_by_id = {str(row["run_id"]): row for row in run_rows}
    if set(snapshot_projects) != set(snapshots_by_id) or set(run_projects) != set(
        runs_by_id
    ):
        raise RuntimeError(_PREFLIGHT_ERROR)

    for snapshot_id, row in snapshots_by_id.items():
        binding = snapshot_projects[snapshot_id]
        try:
            canonical_payload = _canonical_stored_workflow_inputs(
                row["canonical_payload"],
                require_canonical_text=True,
            )
        except (TypeError, ValueError):
            raise RuntimeError(_PREFLIGHT_ERROR) from None
        if (
            row["payload_digest_scheme"] != _WORKFLOW_INPUTS_DIGEST_SCHEME
            or str(row["payload_digest"]) != _workflow_inputs_digest(canonical_payload)
            or str(binding["workflow_inputs_digest"]) != str(row["payload_digest"])
        ):
            raise RuntimeError(_PREFLIGHT_ERROR)
        _validate_project_sample_binding(binding, snapshot_samples[snapshot_id])

    consumed_snapshot_by_run: dict[str, str] = {}
    for snapshot_id, row in snapshots_by_id.items():
        consumed_run_id = row["consumed_run_id"]
        if consumed_run_id is None:
            continue
        run_id = str(consumed_run_id)
        if run_id in consumed_snapshot_by_run or run_id not in runs_by_id:
            raise RuntimeError(_PREFLIGHT_ERROR)
        consumed_snapshot_by_run[run_id] = snapshot_id

    for run_id, row in runs_by_id.items():
        binding = run_projects[run_id]
        _validate_project_sample_binding(binding, run_samples[run_id])
        consumed_snapshot_id = consumed_snapshot_by_run.get(run_id)
        if consumed_snapshot_id is not None:
            snapshot_row = snapshots_by_id[consumed_snapshot_id]
            snapshot_binding = snapshot_projects[consumed_snapshot_id]
            try:
                run_payload = _canonical_stored_workflow_inputs(row["inputs"])
            except (TypeError, ValueError):
                raise RuntimeError(_PREFLIGHT_ERROR) from None
            if (
                str(row["workflow_id"]) != str(snapshot_row["workflow_id"])
                or run_payload != str(snapshot_row["canonical_payload"])
                or _project_sample_evidence(binding, run_samples[run_id])
                != _project_sample_evidence(
                    snapshot_binding,
                    snapshot_samples[consumed_snapshot_id],
                )
            ):
                raise RuntimeError(_PREFLIGHT_ERROR)
        else:
            expected_workflow_digest = _standalone_workflow_inputs_digest(row["inputs"])
            if str(binding["workflow_inputs_digest"]) != expected_workflow_digest:
                raise RuntimeError(_PREFLIGHT_ERROR)

    snapshot_plan: dict[str, dict[str, Any]] = {}
    consumed_run_plan: dict[str, dict[str, Any]] = {}
    for snapshot_id, row in snapshots_by_id.items():
        project_binding = snapshot_projects[snapshot_id]
        plan = _compatibility_envelope_plan(
            project_id=str(project_binding["project_id"]),
            workflow_id=str(row["workflow_id"]),
            adapter_contract_version=None,
            workflow_inputs_digest=str(project_binding["workflow_inputs_digest"]),
            project_sample_binding_digest=str(project_binding["binding_digest"]),
        )
        snapshot_plan[snapshot_id] = plan
        consumed_run_id = row["consumed_run_id"]
        if consumed_run_id is not None:
            consumed_run_plan[str(consumed_run_id)] = plan

    run_plan: dict[str, dict[str, Any]] = {}
    for run_id, row in runs_by_id.items():
        plan = consumed_run_plan.get(run_id)
        if plan is None:
            project_binding = run_projects[run_id]
            plan = _compatibility_envelope_plan(
                project_id=str(project_binding["project_id"]),
                workflow_id=str(row["workflow_id"]),
                adapter_contract_version=None,
                workflow_inputs_digest=str(project_binding["workflow_inputs_digest"]),
                project_sample_binding_digest=str(project_binding["binding_digest"]),
            )
        run_plan[run_id] = plan
    return snapshot_plan, run_plan


def _one_row_per_owner(
    rows: Sequence[sa.RowMapping],
    *,
    owner_column: str,
) -> dict[str, sa.RowMapping]:
    result: dict[str, sa.RowMapping] = {}
    for row in rows:
        owner_id = str(row[owner_column])
        if owner_id in result:
            raise RuntimeError(_PREFLIGHT_ERROR)
        result[owner_id] = row
    return result


def _rows_per_owner(
    rows: Sequence[sa.RowMapping],
    *,
    owner_column: str,
) -> defaultdict[str, list[sa.RowMapping]]:
    result: defaultdict[str, list[sa.RowMapping]] = defaultdict(list)
    for row in rows:
        result[str(row[owner_column])].append(row)
    return result


def _validate_project_sample_binding(
    binding: Mapping[str, Any],
    sample_rows: Sequence[Mapping[str, Any]],
) -> None:
    if tuple(int(row["ordinal"]) for row in sample_rows) != tuple(
        range(len(sample_rows))
    ) or any(
        str(row["project_id"]) != str(binding["project_id"]) for row in sample_rows
    ):
        raise RuntimeError(_PREFLIGHT_ERROR)
    sample_revisions = [
        {
            "payload_digest": str(row["payload_digest"]),
            "sample_revision_id": str(row["sample_revision_id"]),
        }
        for row in sample_rows
    ]
    binding_mode = str(binding["binding_mode"])
    provenance = str(binding["provenance"])
    project_id = str(binding["project_id"])
    if binding["binding_digest_scheme"] != _PROJECT_SAMPLE_DIGEST_SCHEME:
        raise RuntimeError(_PREFLIGHT_ERROR)
    if (
        binding_mode == "legacy_v1"
        and (
            provenance != "unresolved"
            or project_id != LEGACY_PROJECT_ID
            or sample_revisions
        )
    ) or (
        binding_mode == "bound_v1"
        and (
            provenance != "resolved"
            or project_id == LEGACY_PROJECT_ID
            or not sample_revisions
        )
    ):
        raise RuntimeError(_PREFLIGHT_ERROR)
    if binding_mode not in {"legacy_v1", "bound_v1"}:
        raise RuntimeError(_PREFLIGHT_ERROR)
    expected_digest = _framed_digest(
        _PROJECT_SAMPLE_DIGEST_SCHEME,
        {
            "binding_mode": binding_mode,
            "project_id": project_id,
            "provenance": provenance,
            "sample_revisions": sample_revisions,
            "workflow_inputs_digest": str(binding["workflow_inputs_digest"]),
        },
    )
    if str(binding["binding_digest"]) != expected_digest:
        raise RuntimeError(_PREFLIGHT_ERROR)


def _project_sample_evidence(
    binding: Mapping[str, Any],
    sample_rows: Sequence[Mapping[str, Any]],
) -> tuple[Any, ...]:
    return (
        str(binding["project_id"]),
        str(binding["binding_mode"]),
        str(binding["provenance"]),
        str(binding["workflow_inputs_digest"]),
        str(binding["binding_digest_scheme"]),
        str(binding["binding_digest"]),
        tuple(
            (
                str(row["sample_revision_id"]),
                str(row["payload_digest"]),
                int(row["ordinal"]),
            )
            for row in sample_rows
        ),
    )


def _compatibility_envelope_plan(
    *,
    project_id: str,
    workflow_id: str,
    adapter_contract_version: str | None,
    workflow_inputs_digest: str,
    project_sample_binding_digest: str,
) -> dict[str, Any]:
    payload = {
        "adapter_contract_version": adapter_contract_version,
        "contract_mode": "compatibility_unresolved_v1",
        "input_uses": [],
        "project_id": project_id,
        "project_sample_binding_digest": project_sample_binding_digest,
        "workflow_id": workflow_id,
        "workflow_inputs_digest": workflow_inputs_digest,
    }
    return {
        "project_id": project_id,
        "workflow_id": workflow_id,
        "adapter_contract_version": adapter_contract_version,
        "binding_mode": "compatibility_unresolved_v1",
        "workflow_inputs_digest": workflow_inputs_digest,
        "project_sample_binding_digest": project_sample_binding_digest,
        "binding_digest_scheme": _INPUT_BINDING_DIGEST_SCHEME,
        "binding_digest": _framed_digest(_INPUT_BINDING_DIGEST_SCHEME, payload),
    }


def _standalone_workflow_inputs_digest(inputs: Any) -> str:
    try:
        decoded = json.loads(inputs) if isinstance(inputs, str) else inputs
        if not isinstance(decoded, dict):
            raise ValueError
        canonical = json.dumps(
            decoded,
            ensure_ascii=False,
            allow_nan=False,
            separators=_CANONICAL_JSON_SEPARATORS,
            sort_keys=True,
        )
    except (TypeError, ValueError):
        raise RuntimeError(_PREFLIGHT_ERROR) from None
    return sha256(canonical.encode()).hexdigest()


def _canonical_stored_workflow_inputs(
    value: Any,
    *,
    require_canonical_text: bool = False,
) -> str:
    decoded = json.loads(value) if isinstance(value, str) else value
    if not isinstance(decoded, dict) or set(decoded) != {
        "config",
        "options",
        "samples",
    }:
        raise ValueError("stored workflow inputs are invalid")
    canonical = json.dumps(
        decoded,
        ensure_ascii=False,
        allow_nan=False,
        separators=_CANONICAL_JSON_SEPARATORS,
        sort_keys=True,
    )
    if require_canonical_text and value != canonical:
        raise ValueError("stored snapshot payload is not canonical")
    return canonical


def _workflow_inputs_digest(canonical_payload: str) -> str:
    digest = sha256()
    for part in (_WORKFLOW_INPUTS_DIGEST_SCHEME, canonical_payload):
        encoded = part.encode()
        digest.update(len(encoded).to_bytes(8, "big"))
        digest.update(encoded)
    return digest.hexdigest()


def _framed_digest(scheme: str, payload: Mapping[str, Any]) -> str:
    canonical = json.dumps(
        payload,
        ensure_ascii=True,
        allow_nan=False,
        separators=_CANONICAL_JSON_SEPARATORS,
        sort_keys=True,
    )
    digest = sha256()
    for part in (scheme, canonical):
        encoded = part.encode("utf-8")
        digest.update(len(encoded).to_bytes(8, byteorder="big", signed=False))
        digest.update(encoded)
    return digest.hexdigest()


def _backfill_snapshot_input_bindings(
    connection: sa.Connection,
    snapshot_rows: Sequence[sa.RowMapping],
    plan: Mapping[str, Mapping[str, Any]],
) -> None:
    for row in snapshot_rows:
        snapshot_id = str(row["snapshot_id"])
        _insert_input_binding(
            connection,
            owner="snapshot",
            owner_id=snapshot_id,
            created_at=row["validated_at"],
            plan=plan[snapshot_id],
        )


def _backfill_run_input_bindings(
    connection: sa.Connection,
    run_rows: Sequence[sa.RowMapping],
    plan: Mapping[str, Mapping[str, Any]],
) -> None:
    for row in run_rows:
        run_id = str(row["run_id"])
        _insert_input_binding(
            connection,
            owner="run",
            owner_id=run_id,
            created_at=row["created_at"],
            plan=plan[run_id],
        )


def _insert_input_binding(
    connection: sa.Connection,
    *,
    owner: str,
    owner_id: str,
    created_at: Any,
    plan: Mapping[str, Any],
) -> None:
    owner_column = "snapshot_id" if owner == "snapshot" else "run_id"
    connection.execute(
        sa.text(
            f"INSERT INTO {owner}_input_bindings "
            f"({owner_column}, project_id, workflow_id, "
            "adapter_contract_version, binding_mode, workflow_inputs_digest, "
            "project_sample_binding_digest, binding_digest_scheme, "
            "binding_digest, created_at) VALUES "
            "(:owner_id, :project_id, :workflow_id, "
            ":adapter_contract_version, :binding_mode, "
            ":workflow_inputs_digest, :project_sample_binding_digest, "
            ":binding_digest_scheme, :binding_digest, :created_at)"
        ),
        {
            "owner_id": owner_id,
            "created_at": created_at,
            **plan,
        },
    )


def _create_archive_only_triggers(
    *,
    table_name: str,
    immutable_columns: tuple[str, ...],
    label: str,
) -> None:
    unchanged = " AND ".join(
        f"NEW.{column} IS OLD.{column}" for column in immutable_columns
    )
    timestamp_column = (
        "bound_at" if table_name == "project_storage_pool_bindings" else ("created_at")
    )
    op.execute(
        f"""
        CREATE TRIGGER trg_{table_name}_archive_only
        BEFORE UPDATE ON {table_name}
        WHEN NOT (
            {unchanged}
            AND OLD.archived_at IS NULL
            AND NEW.archived_at IS NOT NULL
            AND NEW.archived_at >= OLD.{timestamp_column}
        )
        BEGIN
            SELECT RAISE(ABORT, '{label} permits only one archive transition');
        END
        """
    )
    op.execute(
        f"""
        CREATE TRIGGER trg_{table_name}_no_delete
        BEFORE DELETE ON {table_name}
        BEGIN
            SELECT RAISE(ABORT, '{label} identity cannot be deleted');
        END
        """
    )


def _create_immutable_triggers(
    table_name: str,
    label: str,
    *,
    cascade_delete_root: tuple[str, str] | None = None,
) -> None:
    delete_guard = ""
    if cascade_delete_root is not None:
        root_table, root_column = cascade_delete_root
        delete_guard = (
            f"WHEN EXISTS (SELECT 1 FROM {root_table} "
            f"WHERE {root_column} = OLD.{root_column})"
        )
    op.execute(
        f"""
        CREATE TRIGGER trg_{table_name}_no_update
        BEFORE UPDATE ON {table_name}
        BEGIN
            SELECT RAISE(ABORT, '{label} is immutable');
        END
        """
    )
    op.execute(
        f"""
        CREATE TRIGGER trg_{table_name}_no_delete
        BEFORE DELETE ON {table_name}
        {delete_guard}
        BEGIN
            SELECT RAISE(ABORT, '{label} is immutable');
        END
        """
    )


def _drop_archive_only_triggers(*table_names: str) -> None:
    for table_name in table_names:
        op.execute(f"DROP TRIGGER IF EXISTS trg_{table_name}_archive_only")
        op.execute(f"DROP TRIGGER IF EXISTS trg_{table_name}_no_delete")


def _drop_immutable_triggers(*table_names: str) -> None:
    for table_name in table_names:
        op.execute(f"DROP TRIGGER IF EXISTS trg_{table_name}_no_update")
        op.execute(f"DROP TRIGGER IF EXISTS trg_{table_name}_no_delete")
