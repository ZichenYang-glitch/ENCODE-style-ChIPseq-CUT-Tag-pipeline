"""Add the Project and immutable Sample registry.

Revision ID: 20260726_09
Revises: 20260717_08
"""

from __future__ import annotations

from collections.abc import Sequence
from hashlib import sha256
import json
from typing import Any

from alembic import op
import sqlalchemy as sa


revision = "20260726_09"
down_revision = "20260717_08"
branch_labels = None
depends_on = None


LEGACY_PROJECT_ID = "prj_00000000000000000000000000000000"
_BINDING_DIGEST_SCHEME = "sha256-framed-data-binding-v1"
_WORKFLOW_INPUTS_DIGEST_SCHEME = "sha256-framed-workflow-inputs-v1"
_CANONICAL_JSON_SEPARATORS = (",", ":")


def upgrade() -> None:
    connection = op.get_bind()
    snapshot_rows = (
        connection.execute(
            sa.text(
                "SELECT snapshot_id, workflow_id, canonical_payload, "
                "payload_digest_scheme, payload_digest, validated_at, "
                "consumed_run_id, consumed_at "
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
    _validate_consumed_snapshot_links(snapshot_rows, run_rows)

    _create_registry_tables()
    _create_binding_tables()
    _backfill_legacy_project(connection)
    snapshot_digests = _backfill_snapshot_bindings(connection, snapshot_rows)
    _backfill_run_bindings(connection, run_rows, snapshot_digests)


def downgrade() -> None:
    op.drop_index("ix_run_samples_revision", table_name="run_samples")
    op.drop_table("run_samples")
    op.drop_index(
        "ix_run_project_bindings_project",
        table_name="run_project_bindings",
    )
    op.drop_table("run_project_bindings")
    op.drop_index(
        "ix_snapshot_sample_revisions_revision",
        table_name="snapshot_sample_revisions",
    )
    op.drop_table("snapshot_sample_revisions")
    op.drop_index(
        "ix_snapshot_project_bindings_project",
        table_name="snapshot_project_bindings",
    )
    op.drop_table("snapshot_project_bindings")
    op.drop_index(
        "ix_sample_revisions_project_created",
        table_name="sample_revisions",
    )
    op.drop_index(
        "ix_sample_revisions_sample_revision",
        table_name="sample_revisions",
    )
    op.execute("DROP TRIGGER IF EXISTS trg_sample_revisions_no_update")
    op.execute("DROP TRIGGER IF EXISTS trg_sample_revisions_no_delete")
    op.drop_table("sample_revisions")
    op.drop_index("ix_samples_project_created", table_name="samples")
    op.execute("DROP TRIGGER IF EXISTS trg_samples_no_update")
    op.execute("DROP TRIGGER IF EXISTS trg_samples_no_delete")
    op.drop_table("samples")
    op.execute("DROP TRIGGER IF EXISTS trg_projects_legacy_no_update")
    op.execute("DROP TRIGGER IF EXISTS trg_projects_legacy_no_delete")
    op.drop_index("ix_projects_archived_created", table_name="projects")
    op.drop_table("projects")


def _create_registry_tables() -> None:
    op.create_table(
        "projects",
        sa.Column("project_id", sa.String(length=36), nullable=False),
        sa.Column("display_name", sa.String(length=255), nullable=False),
        sa.Column("kind", sa.String(length=32), nullable=False),
        sa.Column("created_at", sa.DateTime(timezone=True), nullable=False),
        sa.Column("archived_at", sa.DateTime(timezone=True), nullable=True),
        sa.CheckConstraint(
            "length(project_id) = 36 AND substr(project_id, 1, 4) = 'prj_'",
            name="ck_projects_id",
        ),
        sa.CheckConstraint(
            "length(trim(display_name)) BETWEEN 1 AND 255",
            name="ck_projects_display_name",
        ),
        sa.CheckConstraint(
            "kind IN ('user', 'system')",
            name="ck_projects_kind",
        ),
        sa.CheckConstraint(
            "(kind = 'system' AND "
            f"project_id = '{LEGACY_PROJECT_ID}' AND "
            "display_name = 'Legacy Project' AND archived_at IS NULL) OR "
            "(kind = 'user' AND "
            f"project_id != '{LEGACY_PROJECT_ID}')",
            name="ck_projects_legacy_identity",
        ),
        sa.CheckConstraint(
            "archived_at IS NULL OR archived_at >= created_at",
            name="ck_projects_archive_order",
        ),
        sa.PrimaryKeyConstraint("project_id"),
    )
    op.create_index(
        "ix_projects_archived_created",
        "projects",
        ["archived_at", "created_at", "project_id"],
    )
    op.execute(
        f"""
        CREATE TRIGGER trg_projects_legacy_no_update
        BEFORE UPDATE ON projects
        WHEN OLD.project_id = '{LEGACY_PROJECT_ID}'
        BEGIN
            SELECT RAISE(ABORT, 'reserved Legacy Project is immutable');
        END
        """
    )
    op.execute(
        f"""
        CREATE TRIGGER trg_projects_legacy_no_delete
        BEFORE DELETE ON projects
        WHEN OLD.project_id = '{LEGACY_PROJECT_ID}'
        BEGIN
            SELECT RAISE(ABORT, 'reserved Legacy Project is immutable');
        END
        """
    )

    op.create_table(
        "samples",
        sa.Column("sample_id", sa.String(length=36), nullable=False),
        sa.Column("project_id", sa.String(length=36), nullable=False),
        sa.Column("stable_key", sa.String(length=255), nullable=False),
        sa.Column("created_at", sa.DateTime(timezone=True), nullable=False),
        sa.CheckConstraint(
            "length(sample_id) = 36 AND substr(sample_id, 1, 4) = 'smp_'",
            name="ck_samples_id",
        ),
        sa.CheckConstraint(
            "length(trim(stable_key)) BETWEEN 1 AND 255",
            name="ck_samples_stable_key",
        ),
        sa.ForeignKeyConstraint(
            ["project_id"],
            ["projects.project_id"],
            name="fk_samples_project",
            ondelete="RESTRICT",
        ),
        sa.PrimaryKeyConstraint("sample_id"),
        sa.UniqueConstraint(
            "project_id",
            "sample_id",
            name="uq_samples_project_sample",
        ),
        sa.UniqueConstraint(
            "project_id",
            "stable_key",
            name="uq_samples_project_stable_key",
        ),
    )
    op.create_index(
        "ix_samples_project_created",
        "samples",
        ["project_id", "created_at", "sample_id"],
    )
    op.execute(
        """
        CREATE TRIGGER trg_samples_no_update
        BEFORE UPDATE ON samples
        BEGIN
            SELECT RAISE(ABORT, 'Sample identity is immutable');
        END
        """
    )
    op.execute(
        """
        CREATE TRIGGER trg_samples_no_delete
        BEFORE DELETE ON samples
        BEGIN
            SELECT RAISE(ABORT, 'Sample identity is immutable');
        END
        """
    )

    op.create_table(
        "sample_revisions",
        sa.Column("sample_revision_id", sa.String(length=37), nullable=False),
        sa.Column("project_id", sa.String(length=36), nullable=False),
        sa.Column("sample_id", sa.String(length=36), nullable=False),
        sa.Column("revision_number", sa.Integer(), nullable=False),
        sa.Column("canonical_payload", sa.Text(), nullable=False),
        sa.Column(
            "payload_digest_scheme",
            sa.String(length=64),
            nullable=False,
        ),
        sa.Column("payload_digest", sa.String(length=64), nullable=False),
        sa.Column("created_at", sa.DateTime(timezone=True), nullable=False),
        sa.CheckConstraint(
            "length(sample_revision_id) = 37 "
            "AND substr(sample_revision_id, 1, 5) = 'smpr_'",
            name="ck_sample_revisions_id",
        ),
        sa.CheckConstraint(
            "revision_number >= 1",
            name="ck_sample_revisions_positive_revision",
        ),
        sa.CheckConstraint(
            "payload_digest_scheme = 'sha256-framed-sample-revision-payload-v1'",
            name="ck_sample_revisions_digest_scheme",
        ),
        sa.CheckConstraint(
            "length(payload_digest) = 64",
            name="ck_sample_revisions_digest_length",
        ),
        sa.ForeignKeyConstraint(
            ["project_id", "sample_id"],
            ["samples.project_id", "samples.sample_id"],
            name="fk_sample_revisions_sample",
            ondelete="RESTRICT",
        ),
        sa.PrimaryKeyConstraint("sample_revision_id"),
        sa.UniqueConstraint(
            "project_id",
            "sample_revision_id",
            name="uq_sample_revisions_project_revision",
        ),
        sa.UniqueConstraint(
            "project_id",
            "sample_revision_id",
            "payload_digest",
            name="uq_sample_revisions_project_revision_digest",
        ),
        sa.UniqueConstraint(
            "sample_id",
            "revision_number",
            name="uq_sample_revisions_sample_number",
        ),
    )
    op.create_index(
        "ix_sample_revisions_sample_revision",
        "sample_revisions",
        ["sample_id", "revision_number"],
    )
    op.create_index(
        "ix_sample_revisions_project_created",
        "sample_revisions",
        ["project_id", "created_at", "sample_revision_id"],
    )
    op.execute(
        """
        CREATE TRIGGER trg_sample_revisions_no_update
        BEFORE UPDATE ON sample_revisions
        BEGIN
            SELECT RAISE(ABORT, 'SampleRevision is immutable');
        END
        """
    )
    op.execute(
        """
        CREATE TRIGGER trg_sample_revisions_no_delete
        BEFORE DELETE ON sample_revisions
        BEGIN
            SELECT RAISE(ABORT, 'SampleRevision is immutable');
        END
        """
    )


def _create_binding_tables() -> None:
    op.create_table(
        "snapshot_project_bindings",
        sa.Column("snapshot_id", sa.String(length=64), nullable=False),
        sa.Column("project_id", sa.String(length=36), nullable=False),
        sa.Column("binding_mode", sa.String(length=32), nullable=False),
        sa.Column("provenance", sa.String(length=32), nullable=False),
        sa.Column("workflow_inputs_digest", sa.String(length=64), nullable=False),
        sa.Column(
            "binding_digest_scheme",
            sa.String(length=64),
            nullable=False,
        ),
        sa.Column("binding_digest", sa.String(length=64), nullable=False),
        sa.Column("created_at", sa.DateTime(timezone=True), nullable=False),
        sa.CheckConstraint(
            "binding_mode IN ('legacy_v1', 'bound_v1')",
            name="ck_snapshot_project_bindings_mode",
        ),
        sa.CheckConstraint(
            "provenance IN ('resolved', 'unresolved')",
            name="ck_snapshot_project_bindings_provenance",
        ),
        sa.CheckConstraint(
            f"(binding_mode = 'legacy_v1' AND provenance = 'unresolved' "
            f"AND project_id = '{LEGACY_PROJECT_ID}') OR "
            "(binding_mode = 'bound_v1' AND provenance = 'resolved' "
            f"AND project_id != '{LEGACY_PROJECT_ID}')",
            name="ck_snapshot_project_bindings_legacy_project",
        ),
        sa.CheckConstraint(
            f"binding_digest_scheme = '{_BINDING_DIGEST_SCHEME}'",
            name="ck_snapshot_project_bindings_digest_scheme",
        ),
        sa.CheckConstraint(
            "length(binding_digest) = 64",
            name="ck_snapshot_project_bindings_digest_length",
        ),
        sa.CheckConstraint(
            "length(workflow_inputs_digest) = 64",
            name="ck_snapshot_project_bindings_workflow_digest_length",
        ),
        sa.ForeignKeyConstraint(
            ["snapshot_id"],
            ["validated_input_snapshots.snapshot_id"],
            name="fk_snapshot_project_bindings_snapshot",
            ondelete="RESTRICT",
        ),
        sa.ForeignKeyConstraint(
            ["project_id"],
            ["projects.project_id"],
            name="fk_snapshot_project_bindings_project",
            ondelete="RESTRICT",
        ),
        sa.PrimaryKeyConstraint("snapshot_id"),
        sa.UniqueConstraint(
            "snapshot_id",
            "project_id",
            name="uq_snapshot_project_bindings_project",
        ),
    )
    op.create_index(
        "ix_snapshot_project_bindings_project",
        "snapshot_project_bindings",
        ["project_id", "created_at", "snapshot_id"],
    )

    op.create_table(
        "snapshot_sample_revisions",
        sa.Column("snapshot_id", sa.String(length=64), nullable=False),
        sa.Column("project_id", sa.String(length=36), nullable=False),
        sa.Column("sample_revision_id", sa.String(length=37), nullable=False),
        sa.Column("payload_digest", sa.String(length=64), nullable=False),
        sa.Column("ordinal", sa.Integer(), nullable=False),
        sa.CheckConstraint(
            "ordinal >= 0",
            name="ck_snapshot_sample_revisions_ordinal",
        ),
        sa.ForeignKeyConstraint(
            ["snapshot_id", "project_id"],
            [
                "snapshot_project_bindings.snapshot_id",
                "snapshot_project_bindings.project_id",
            ],
            name="fk_snapshot_sample_revisions_binding",
            ondelete="RESTRICT",
        ),
        sa.ForeignKeyConstraint(
            ["project_id", "sample_revision_id", "payload_digest"],
            [
                "sample_revisions.project_id",
                "sample_revisions.sample_revision_id",
                "sample_revisions.payload_digest",
            ],
            name="fk_snapshot_sample_revisions_revision",
            ondelete="RESTRICT",
        ),
        sa.PrimaryKeyConstraint("snapshot_id", "ordinal"),
        sa.UniqueConstraint(
            "snapshot_id",
            "sample_revision_id",
            name="uq_snapshot_sample_revisions_revision",
        ),
    )
    op.create_index(
        "ix_snapshot_sample_revisions_revision",
        "snapshot_sample_revisions",
        ["sample_revision_id", "snapshot_id"],
    )

    op.create_table(
        "run_project_bindings",
        sa.Column("run_id", sa.String(length=128), nullable=False),
        sa.Column("project_id", sa.String(length=36), nullable=False),
        sa.Column("binding_mode", sa.String(length=32), nullable=False),
        sa.Column("provenance", sa.String(length=32), nullable=False),
        sa.Column("workflow_inputs_digest", sa.String(length=64), nullable=False),
        sa.Column(
            "binding_digest_scheme",
            sa.String(length=64),
            nullable=False,
        ),
        sa.Column("binding_digest", sa.String(length=64), nullable=False),
        sa.Column("created_at", sa.DateTime(timezone=True), nullable=False),
        sa.CheckConstraint(
            "binding_mode IN ('legacy_v1', 'bound_v1')",
            name="ck_run_project_bindings_mode",
        ),
        sa.CheckConstraint(
            "provenance IN ('resolved', 'unresolved')",
            name="ck_run_project_bindings_provenance",
        ),
        sa.CheckConstraint(
            f"(binding_mode = 'legacy_v1' AND provenance = 'unresolved' "
            f"AND project_id = '{LEGACY_PROJECT_ID}') OR "
            "(binding_mode = 'bound_v1' AND provenance = 'resolved' "
            f"AND project_id != '{LEGACY_PROJECT_ID}')",
            name="ck_run_project_bindings_legacy_project",
        ),
        sa.CheckConstraint(
            f"binding_digest_scheme = '{_BINDING_DIGEST_SCHEME}'",
            name="ck_run_project_bindings_digest_scheme",
        ),
        sa.CheckConstraint(
            "length(binding_digest) = 64",
            name="ck_run_project_bindings_digest_length",
        ),
        sa.CheckConstraint(
            "length(workflow_inputs_digest) = 64",
            name="ck_run_project_bindings_workflow_digest_length",
        ),
        sa.ForeignKeyConstraint(
            ["run_id"],
            ["runs.run_id"],
            name="fk_run_project_bindings_run",
            ondelete="CASCADE",
        ),
        sa.ForeignKeyConstraint(
            ["project_id"],
            ["projects.project_id"],
            name="fk_run_project_bindings_project",
            ondelete="RESTRICT",
        ),
        sa.PrimaryKeyConstraint("run_id"),
        sa.UniqueConstraint(
            "run_id",
            "project_id",
            name="uq_run_project_bindings_project",
        ),
    )
    op.create_index(
        "ix_run_project_bindings_project",
        "run_project_bindings",
        ["project_id", "created_at", "run_id"],
    )

    op.create_table(
        "run_samples",
        sa.Column("run_id", sa.String(length=128), nullable=False),
        sa.Column("project_id", sa.String(length=36), nullable=False),
        sa.Column("sample_revision_id", sa.String(length=37), nullable=False),
        sa.Column("payload_digest", sa.String(length=64), nullable=False),
        sa.Column("ordinal", sa.Integer(), nullable=False),
        sa.CheckConstraint("ordinal >= 0", name="ck_run_samples_ordinal"),
        sa.ForeignKeyConstraint(
            ["run_id", "project_id"],
            [
                "run_project_bindings.run_id",
                "run_project_bindings.project_id",
            ],
            name="fk_run_samples_binding",
            ondelete="CASCADE",
        ),
        sa.ForeignKeyConstraint(
            ["project_id", "sample_revision_id", "payload_digest"],
            [
                "sample_revisions.project_id",
                "sample_revisions.sample_revision_id",
                "sample_revisions.payload_digest",
            ],
            name="fk_run_samples_revision",
            ondelete="RESTRICT",
        ),
        sa.PrimaryKeyConstraint("run_id", "ordinal"),
        sa.UniqueConstraint(
            "run_id",
            "sample_revision_id",
            name="uq_run_samples_revision",
        ),
    )
    op.create_index(
        "ix_run_samples_revision",
        "run_samples",
        ["sample_revision_id", "run_id"],
    )


def _validate_consumed_snapshot_links(
    snapshot_rows: Sequence[sa.RowMapping],
    run_rows: Sequence[sa.RowMapping],
) -> None:
    runs_by_id = {str(row["run_id"]): row for row in run_rows}
    consumed_run_ids: set[str] = set()
    for row in snapshot_rows:
        try:
            snapshot_payload = _canonical_stored_workflow_inputs(
                row["canonical_payload"],
                require_canonical_text=True,
            )
            snapshot_digest = _workflow_inputs_digest(snapshot_payload)
        except (TypeError, ValueError):
            raise RuntimeError(
                "inconsistent snapshot evidence blocks registry upgrade"
            ) from None
        if (
            row["payload_digest_scheme"] != _WORKFLOW_INPUTS_DIGEST_SCHEME
            or str(row["payload_digest"]) != snapshot_digest
        ):
            raise RuntimeError("inconsistent snapshot evidence blocks registry upgrade")
        consumed_run_id = row["consumed_run_id"]
        if (consumed_run_id is None) != (row["consumed_at"] is None):
            raise RuntimeError(
                "inconsistent consumed snapshot linkage blocks registry upgrade"
            )
        if consumed_run_id is None:
            continue
        normalized_run_id = str(consumed_run_id)
        run_row = runs_by_id.get(normalized_run_id)
        if run_row is None or normalized_run_id in consumed_run_ids:
            raise RuntimeError(
                "inconsistent consumed snapshot linkage blocks registry upgrade"
            )
        try:
            run_payload = _canonical_stored_workflow_inputs(run_row["inputs"])
        except (TypeError, ValueError):
            raise RuntimeError(
                "inconsistent consumed snapshot evidence blocks registry upgrade"
            ) from None
        if (
            str(row["workflow_id"]) != str(run_row["workflow_id"])
            or snapshot_payload != run_payload
        ):
            raise RuntimeError(
                "inconsistent consumed snapshot evidence blocks registry upgrade"
            )
        consumed_run_ids.add(normalized_run_id)


def _backfill_legacy_project(connection: sa.Connection) -> None:
    connection.execute(
        sa.text(
            "INSERT INTO projects "
            "(project_id, display_name, kind, created_at, archived_at) VALUES "
            "(:project_id, 'Legacy Project', 'system', "
            "'1970-01-01 00:00:00.000000', NULL)"
        ),
        {"project_id": LEGACY_PROJECT_ID},
    )


def _backfill_snapshot_bindings(
    connection: sa.Connection,
    snapshot_rows: Sequence[sa.RowMapping],
) -> dict[str, tuple[str, str]]:
    consumed_run_digests: dict[str, tuple[str, str]] = {}
    for row in snapshot_rows:
        workflow_inputs_digest = str(row["payload_digest"])
        binding_digest = _legacy_binding_digest(workflow_inputs_digest)
        connection.execute(
            sa.text(
                "INSERT INTO snapshot_project_bindings "
                "(snapshot_id, project_id, binding_mode, provenance, "
                "workflow_inputs_digest, binding_digest_scheme, binding_digest, "
                "created_at) VALUES "
                "(:snapshot_id, :project_id, 'legacy_v1', 'unresolved', "
                ":workflow_inputs_digest, :digest_scheme, :binding_digest, "
                ":created_at)"
            ),
            {
                "snapshot_id": str(row["snapshot_id"]),
                "project_id": LEGACY_PROJECT_ID,
                "workflow_inputs_digest": workflow_inputs_digest,
                "digest_scheme": _BINDING_DIGEST_SCHEME,
                "binding_digest": binding_digest,
                "created_at": row["validated_at"],
            },
        )
        consumed_run_id = row["consumed_run_id"]
        if consumed_run_id is not None:
            consumed_run_digests[str(consumed_run_id)] = (
                workflow_inputs_digest,
                binding_digest,
            )
    return consumed_run_digests


def _backfill_run_bindings(
    connection: sa.Connection,
    run_rows: Sequence[sa.RowMapping],
    snapshot_digests: dict[str, tuple[str, str]],
) -> None:
    for row in run_rows:
        run_id = str(row["run_id"])
        snapshot_binding = snapshot_digests.get(run_id)
        if snapshot_binding is None:
            workflow_inputs_digest = _standalone_workflow_inputs_digest(row["inputs"])
            binding_digest = _legacy_binding_digest(workflow_inputs_digest)
        else:
            workflow_inputs_digest, binding_digest = snapshot_binding
        connection.execute(
            sa.text(
                "INSERT INTO run_project_bindings "
                "(run_id, project_id, binding_mode, provenance, "
                "workflow_inputs_digest, binding_digest_scheme, binding_digest, "
                "created_at) VALUES "
                "(:run_id, :project_id, 'legacy_v1', 'unresolved', "
                ":workflow_inputs_digest, :digest_scheme, :binding_digest, "
                ":created_at)"
            ),
            {
                "run_id": run_id,
                "project_id": LEGACY_PROJECT_ID,
                "workflow_inputs_digest": workflow_inputs_digest,
                "digest_scheme": _BINDING_DIGEST_SCHEME,
                "binding_digest": binding_digest,
                "created_at": row["created_at"],
            },
        )


def _legacy_binding_digest(workflow_inputs_digest: str) -> str:
    return _framed_canonical_digest(
        {
            "binding_mode": "legacy_v1",
            "input_revisions": [],
            "project_id": LEGACY_PROJECT_ID,
            "provenance": "unresolved",
            "sample_revisions": [],
            "workflow_inputs_digest": workflow_inputs_digest,
        }
    )


def _standalone_workflow_inputs_digest(inputs: Any) -> str:
    try:
        decoded = json.loads(inputs) if isinstance(inputs, str) else inputs
        if not isinstance(decoded, dict):
            raise ValueError("stored run inputs are invalid")
        canonical_inputs = _canonical_stored_json(decoded)
    except (TypeError, ValueError):
        raise RuntimeError("invalid legacy run inputs block registry upgrade") from None
    return sha256(canonical_inputs.encode()).hexdigest()


def _canonical_stored_workflow_inputs(
    value: Any,
    *,
    require_canonical_text: bool = False,
) -> str:
    if isinstance(value, str):
        decoded = json.loads(value)
    else:
        decoded = value
    if not isinstance(decoded, dict) or set(decoded) != {
        "config",
        "options",
        "samples",
    }:
        raise ValueError("stored workflow inputs are invalid")
    canonical = _canonical_stored_json(decoded)
    if require_canonical_text and value != canonical:
        raise ValueError("stored snapshot payload is not canonical")
    return canonical


def _canonical_stored_json(value: Any) -> str:
    if isinstance(value, str):
        decoded = json.loads(value)
    else:
        decoded = value
    return json.dumps(
        decoded,
        ensure_ascii=False,
        allow_nan=False,
        separators=_CANONICAL_JSON_SEPARATORS,
        sort_keys=True,
    )


def _workflow_inputs_digest(canonical_payload: str) -> str:
    digest = sha256()
    for part in (_WORKFLOW_INPUTS_DIGEST_SCHEME, canonical_payload):
        encoded = part.encode()
        digest.update(len(encoded).to_bytes(8, "big"))
        digest.update(encoded)
    return digest.hexdigest()


def _framed_canonical_digest(value: dict[str, Any]) -> str:
    digest = sha256()
    for part in (_BINDING_DIGEST_SCHEME, _canonical_json(value)):
        encoded = part.encode("utf-8")
        digest.update(len(encoded).to_bytes(8, byteorder="big", signed=False))
        digest.update(encoded)
    return digest.hexdigest()


def _canonical_json(value: Any) -> str:
    return json.dumps(
        value,
        ensure_ascii=True,
        allow_nan=False,
        separators=_CANONICAL_JSON_SEPARATORS,
        sort_keys=True,
    )
