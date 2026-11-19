"""Register external results and extend download audit targets.

Revision ID: 20261119_39
Revises: 20261119_18

The merged chain is 20261118_17 then 20261119_18. Revision 20261119_18
already creates experiment snapshot bindings, so this change uses a new id.
IRENE-32 reserves 20261119_32, currently parented on 20261118_17 in its
own branch. This revision continues the merged head and does not take 32.
"""

from __future__ import annotations

from alembic import op
import sqlalchemy as sa


revision = "20261119_39"
down_revision = "20261119_18"
branch_labels = None
depends_on = None

_AUDIT_ACTIONS_PREVIOUS = (
    "action IN ('auth.login', 'auth.logout', 'run.create', 'run.start', "
    "'run.cancel', 'artifact.download', 'reference.register', "
    "'reference.enable', 'reference.disable', 'storage.register', "
    "'storage.archive', 'account.create', 'account.enable', "
    "'account.disable', 'account.password_reset', "
    "'account.sessions_revoke', 'run.fail', 'run.requeue')"
)
_AUDIT_ACTIONS = (
    "action IN ('auth.login', 'auth.logout', 'run.create', 'run.start', "
    "'run.cancel', 'artifact.download', 'reference.register', "
    "'reference.enable', 'reference.disable', 'storage.register', "
    "'storage.archive', 'account.create', 'account.enable', "
    "'account.disable', 'account.password_reset', "
    "'account.sessions_revoke', 'run.fail', 'run.requeue', "
    "'registered_result.download')"
)
_AUDIT_RESOURCES_PREVIOUS = (
    "resource_kind IS NULL OR resource_kind IN "
    "('account', 'run', 'artifact', 'reference', 'storage')"
)
_AUDIT_RESOURCES = (
    "resource_kind IS NULL OR resource_kind IN "
    "('account', 'run', 'artifact', 'reference', 'storage', "
    "'registered_result')"
)


def upgrade() -> None:
    op.create_table(
        "registered_results",
        sa.Column("result_id", sa.String(length=37), nullable=False),
        sa.Column("project_id", sa.String(length=36), nullable=False),
        sa.Column("stable_key", sa.String(length=255), nullable=False),
        sa.Column("created_at", sa.DateTime(timezone=True), nullable=False),
        sa.Column("archived_at", sa.DateTime(timezone=True), nullable=True),
        sa.CheckConstraint(
            "length(result_id) = 37 AND substr(result_id, 1, 5) = 'rreg_'",
            name="ck_registered_results_id",
        ),
        sa.CheckConstraint(
            "length(project_id) = 36 AND substr(project_id, 1, 4) = 'prj_'",
            name="ck_registered_results_project",
        ),
        sa.CheckConstraint(
            "length(stable_key) BETWEEN 1 AND 255",
            name="ck_registered_results_stable_key",
        ),
        sa.CheckConstraint(
            "archived_at IS NULL OR archived_at >= created_at",
            name="ck_registered_results_archive_order",
        ),
        sa.ForeignKeyConstraint(
            ["project_id"],
            ["projects.project_id"],
            name="fk_registered_results_project",
            ondelete="RESTRICT",
        ),
        sa.PrimaryKeyConstraint("result_id"),
        sa.UniqueConstraint(
            "project_id",
            "stable_key",
            name="uq_registered_results_project_stable_key",
        ),
    )
    op.create_table(
        "registered_result_revisions",
        sa.Column("result_revision_id", sa.String(length=37), nullable=False),
        sa.Column("result_id", sa.String(length=37), nullable=False),
        sa.Column("revision_number", sa.Integer(), nullable=False),
        sa.Column("canonical_payload", sa.Text(), nullable=False),
        sa.Column("digest_scheme", sa.String(length=80), nullable=False),
        sa.Column("revision_digest", sa.String(length=64), nullable=False),
        sa.Column("created_at", sa.DateTime(timezone=True), nullable=False),
        sa.Column("storage_pool_id", sa.String(length=37), nullable=False),
        sa.Column("relative_path", sa.Text(), nullable=False),
        sa.Column("size_bytes", sa.Integer(), nullable=False),
        sa.Column("content_sha256", sa.String(length=64), nullable=False),
        sa.CheckConstraint(
            "length(result_revision_id) = 37 AND "
            "substr(result_revision_id, 1, 5) = 'rrev_'",
            name="ck_registered_result_revisions_id",
        ),
        sa.CheckConstraint(
            "revision_number >= 1",
            name="ck_registered_result_revisions_number",
        ),
        sa.CheckConstraint(
            "digest_scheme = 'sha256-framed-registered-result-revision-v1'",
            name="ck_registered_result_revisions_scheme",
        ),
        sa.CheckConstraint(
            "length(revision_digest) = 64 AND revision_digest NOT GLOB '*[^0-9a-f]*'",
            name="ck_registered_result_revisions_digest",
        ),
        sa.CheckConstraint(
            "length(content_sha256) = 64 AND content_sha256 NOT GLOB '*[^0-9a-f]*'",
            name="ck_registered_result_revisions_content_sha256",
        ),
        sa.CheckConstraint(
            "size_bytes >= 0",
            name="ck_registered_result_revisions_size",
        ),
        sa.CheckConstraint(
            "length(storage_pool_id) = 37 AND substr(storage_pool_id, 1, 5) = 'stgp_'",
            name="ck_registered_result_revisions_pool",
        ),
        sa.CheckConstraint(
            "length(relative_path) BETWEEN 1 AND 4096",
            name="ck_registered_result_revisions_path",
        ),
        sa.ForeignKeyConstraint(
            ["result_id"],
            ["registered_results.result_id"],
            name="fk_registered_result_revisions_result",
            ondelete="RESTRICT",
        ),
        sa.ForeignKeyConstraint(
            ["storage_pool_id"],
            ["storage_pools.storage_pool_id"],
            name="fk_registered_result_revisions_pool",
            ondelete="RESTRICT",
        ),
        sa.PrimaryKeyConstraint("result_revision_id"),
        sa.UniqueConstraint(
            "result_id",
            "revision_number",
            name="uq_registered_result_revisions_number",
        ),
    )
    op.create_table(
        "registered_result_revision_heads",
        sa.Column("result_id", sa.String(length=37), nullable=False),
        sa.Column("result_revision_id", sa.String(length=37), nullable=False),
        sa.Column("revision_number", sa.Integer(), nullable=False),
        sa.CheckConstraint(
            "revision_number >= 1",
            name="ck_registered_result_revision_heads_number",
        ),
        sa.ForeignKeyConstraint(
            ["result_id"],
            ["registered_results.result_id"],
            name="fk_registered_result_revision_heads_result",
            ondelete="RESTRICT",
        ),
        sa.ForeignKeyConstraint(
            ["result_revision_id"],
            ["registered_result_revisions.result_revision_id"],
            name="fk_registered_result_revision_heads_revision",
            ondelete="RESTRICT",
        ),
        sa.PrimaryKeyConstraint("result_id"),
        sa.UniqueConstraint(
            "result_revision_id",
            name="uq_registered_result_revision_heads_revision",
        ),
    )
    op.create_table(
        "registered_result_parents",
        sa.Column("child_result_revision_id", sa.String(length=37), nullable=False),
        sa.Column("ordinal", sa.Integer(), nullable=False),
        sa.Column("parent_result_revision_id", sa.String(length=37), nullable=False),
        sa.Column("parent_revision_digest", sa.String(length=64), nullable=False),
        sa.Column("role", sa.String(length=64), nullable=False),
        sa.CheckConstraint(
            "ordinal >= 0",
            name="ck_registered_result_parents_ordinal",
        ),
        sa.CheckConstraint(
            "length(parent_revision_digest) = 64 AND "
            "parent_revision_digest NOT GLOB '*[^0-9a-f]*'",
            name="ck_registered_result_parents_digest",
        ),
        sa.CheckConstraint(
            "length(role) BETWEEN 1 AND 64",
            name="ck_registered_result_parents_role",
        ),
        sa.ForeignKeyConstraint(
            ["child_result_revision_id"],
            ["registered_result_revisions.result_revision_id"],
            name="fk_registered_result_parents_child",
            ondelete="RESTRICT",
        ),
        sa.ForeignKeyConstraint(
            ["parent_result_revision_id"],
            ["registered_result_revisions.result_revision_id"],
            name="fk_registered_result_parents_parent",
            ondelete="RESTRICT",
        ),
        sa.PrimaryKeyConstraint("child_result_revision_id", "ordinal"),
        sa.UniqueConstraint(
            "child_result_revision_id",
            "parent_result_revision_id",
            name="uq_registered_result_parents_parent",
        ),
    )
    op.create_index(
        "ix_registered_result_parents_parent",
        "registered_result_parents",
        ["parent_result_revision_id", "child_result_revision_id"],
    )
    op.create_table(
        "registered_result_write_receipts",
        sa.Column("actor_user_id", sa.String(length=36), nullable=False),
        sa.Column("idempotency_key", sa.String(length=128), nullable=False),
        sa.Column("request_sha256", sa.String(length=64), nullable=False),
        sa.Column("result_id", sa.String(length=37), nullable=False),
        sa.Column("result_revision_id", sa.String(length=37), nullable=False),
        sa.Column("created_at", sa.DateTime(timezone=True), nullable=False),
        sa.CheckConstraint(
            "length(request_sha256) = 64 AND request_sha256 NOT GLOB '*[^0-9a-f]*'",
            name="ck_registered_result_write_receipts_digest",
        ),
        sa.CheckConstraint(
            "length(idempotency_key) BETWEEN 1 AND 128",
            name="ck_registered_result_write_receipts_key",
        ),
        sa.ForeignKeyConstraint(
            ["actor_user_id"],
            ["user_accounts.user_id"],
            name="fk_registered_result_write_receipts_actor",
            ondelete="RESTRICT",
        ),
        sa.ForeignKeyConstraint(
            ["result_id"],
            ["registered_results.result_id"],
            name="fk_registered_result_write_receipts_result",
            ondelete="RESTRICT",
        ),
        sa.ForeignKeyConstraint(
            ["result_revision_id"],
            ["registered_result_revisions.result_revision_id"],
            name="fk_registered_result_write_receipts_revision",
            ondelete="RESTRICT",
        ),
        sa.PrimaryKeyConstraint("actor_user_id", "idempotency_key"),
    )
    with op.batch_alter_table("security_audit_events") as batch_op:
        batch_op.drop_constraint("ck_security_audit_events_action", type_="check")
        batch_op.create_check_constraint(
            "ck_security_audit_events_action",
            _AUDIT_ACTIONS,
        )
        batch_op.drop_constraint(
            "ck_security_audit_events_resource_kind",
            type_="check",
        )
        batch_op.create_check_constraint(
            "ck_security_audit_events_resource_kind",
            _AUDIT_RESOURCES,
        )


def downgrade() -> None:
    with op.batch_alter_table("security_audit_events") as batch_op:
        batch_op.drop_constraint("ck_security_audit_events_action", type_="check")
        batch_op.create_check_constraint(
            "ck_security_audit_events_action",
            _AUDIT_ACTIONS_PREVIOUS,
        )
        batch_op.drop_constraint(
            "ck_security_audit_events_resource_kind",
            type_="check",
        )
        batch_op.create_check_constraint(
            "ck_security_audit_events_resource_kind",
            _AUDIT_RESOURCES_PREVIOUS,
        )
    op.drop_table("registered_result_write_receipts")
    op.drop_index(
        "ix_registered_result_parents_parent",
        table_name="registered_result_parents",
    )
    op.drop_table("registered_result_parents")
    op.drop_table("registered_result_revision_heads")
    op.drop_table("registered_result_revisions")
    op.drop_table("registered_results")
