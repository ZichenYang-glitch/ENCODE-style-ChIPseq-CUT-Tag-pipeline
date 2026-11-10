"""Add restricted authoring drafts and proposal receipts.

Revision ID: 20261110_16
Revises: 20260827_15
"""

from __future__ import annotations

from alembic import op
import sqlalchemy as sa


revision = "20261110_16"
down_revision = "20260827_15"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.create_table(
        "agent_authoring_sessions",
        sa.Column("sid", sa.String(length=36), nullable=False),
        sa.Column("owner_user_id", sa.String(length=36), nullable=False),
        sa.Column("project_id", sa.String(length=36), nullable=False),
        sa.Column("workflow_id", sa.String(length=255), nullable=False),
        sa.Column("schema_digest", sa.String(length=64), nullable=False),
        sa.Column("projection_digest", sa.String(length=64), nullable=True),
        sa.Column("draft_json", sa.JSON(), nullable=False),
        sa.Column("draft_digest", sa.String(length=64), nullable=False),
        sa.Column("draft_revision", sa.Integer(), nullable=False),
        sa.Column("created_at", sa.DateTime(timezone=True), nullable=False),
        sa.Column("content_updated_at", sa.DateTime(timezone=True), nullable=False),
        sa.Column("draft_expires_at", sa.DateTime(timezone=True), nullable=False),
        sa.Column("purge_after", sa.DateTime(timezone=True), nullable=False),
        sa.Column("last_snapshot_id", sa.String(length=38), nullable=True),
        sa.Column("last_validated_revision", sa.Integer(), nullable=True),
        sa.Column("validation_issue_codes", sa.JSON(), nullable=True),
        sa.Column("validation_issues", sa.JSON(), nullable=True),
        sa.Column("connection_state", sa.String(length=16), nullable=False),
        sa.Column("connection_reason_code", sa.String(length=64), nullable=True),
        sa.Column("originating_auth_session_ref", sa.String(length=64), nullable=True),
        sa.Column("instance_id", sa.String(length=36), nullable=True),
        sa.Column("generation", sa.Integer(), nullable=False),
        sa.Column("bridge_owner_epoch", sa.String(length=64), nullable=True),
        sa.Column("connection_capability_digest", sa.String(length=64), nullable=True),
        sa.CheckConstraint(
            "length(sid) = 36 AND substr(sid, 1, 4) = 'aas_'",
            name="ck_agent_authoring_sessions_id",
        ),
        sa.CheckConstraint(
            "length(schema_digest) = 64 AND schema_digest NOT GLOB '*[^0-9a-f]*'",
            name="ck_agent_authoring_sessions_schema_digest",
        ),
        sa.CheckConstraint(
            "projection_digest IS NULL OR (length(projection_digest) = 64 "
            "AND projection_digest NOT GLOB '*[^0-9a-f]*')",
            name="ck_agent_authoring_sessions_projection_digest",
        ),
        sa.CheckConstraint(
            "length(draft_digest) = 64 AND draft_digest NOT GLOB '*[^0-9a-f]*'",
            name="ck_agent_authoring_sessions_draft_digest",
        ),
        sa.CheckConstraint(
            "draft_revision >= 0",
            name="ck_agent_authoring_sessions_revision",
        ),
        sa.CheckConstraint(
            "connection_state IN ('closed', 'blocked', 'ready')",
            name="ck_agent_authoring_sessions_connection_state",
        ),
        sa.CheckConstraint(
            "generation >= 0",
            name="ck_agent_authoring_sessions_generation",
        ),
        sa.CheckConstraint(
            "content_updated_at >= created_at",
            name="ck_agent_authoring_sessions_content_order",
        ),
        sa.CheckConstraint(
            "draft_expires_at > content_updated_at",
            name="ck_agent_authoring_sessions_expiry",
        ),
        sa.CheckConstraint(
            "purge_after > draft_expires_at",
            name="ck_agent_authoring_sessions_purge",
        ),
        sa.ForeignKeyConstraint(
            ["owner_user_id"],
            ["user_accounts.user_id"],
            name="fk_agent_authoring_sessions_owner",
            ondelete="RESTRICT",
        ),
        sa.ForeignKeyConstraint(
            ["project_id"],
            ["projects.project_id"],
            name="fk_agent_authoring_sessions_project",
            ondelete="RESTRICT",
        ),
        sa.PrimaryKeyConstraint("sid"),
    )
    op.create_index(
        "ix_agent_authoring_sessions_owner_updated",
        "agent_authoring_sessions",
        ["owner_user_id", "content_updated_at", "sid"],
    )
    op.create_index(
        "ix_agent_authoring_sessions_auth",
        "agent_authoring_sessions",
        ["originating_auth_session_ref"],
    )
    op.create_index(
        "ix_agent_authoring_sessions_epoch",
        "agent_authoring_sessions",
        ["bridge_owner_epoch"],
    )
    op.create_table(
        "agent_authoring_turns",
        sa.Column("turn_id", sa.String(length=36), nullable=False),
        sa.Column("sid", sa.String(length=36), nullable=False),
        sa.Column("client_turn_id", sa.String(length=64), nullable=False),
        sa.Column("request_digest", sa.String(length=64), nullable=False),
        sa.Column("question_json", sa.JSON(), nullable=False),
        sa.Column("originating_auth_session_ref", sa.String(length=64), nullable=False),
        sa.Column("instance_id", sa.String(length=36), nullable=False),
        sa.Column("generation", sa.Integer(), nullable=False),
        sa.Column("bridge_owner_epoch", sa.String(length=64), nullable=False),
        sa.Column("base_revision", sa.Integer(), nullable=False),
        sa.Column("schema_digest", sa.String(length=64), nullable=False),
        sa.Column("projection_digest", sa.String(length=64), nullable=False),
        sa.Column("status", sa.String(length=16), nullable=False),
        sa.Column("failure_reason", sa.String(length=64), nullable=True),
        sa.Column("completed_at", sa.DateTime(timezone=True), nullable=True),
        sa.Column("result_kind", sa.String(length=16), nullable=True),
        sa.Column("result_digest", sa.String(length=64), nullable=True),
        sa.Column("proposal_id", sa.String(length=37), nullable=True),
        sa.Column("proposal_json", sa.JSON(), nullable=True),
        sa.Column("proposal_digest", sa.String(length=64), nullable=True),
        sa.Column("proposal_expires_at", sa.DateTime(timezone=True), nullable=True),
        sa.Column("decision", sa.String(length=16), nullable=True),
        sa.Column("invalidation_reason", sa.String(length=64), nullable=True),
        sa.Column("adopt_request_digest", sa.String(length=64), nullable=True),
        sa.Column("adoption_actor_user_id", sa.String(length=36), nullable=True),
        sa.Column("adoption_auth_session_ref", sa.String(length=64), nullable=True),
        sa.Column("adopted_at", sa.DateTime(timezone=True), nullable=True),
        sa.Column("result_revision", sa.Integer(), nullable=True),
        sa.Column("adoption_changed", sa.Boolean(), nullable=True),
        sa.Column("adoption_content_digest", sa.String(length=64), nullable=True),
        sa.Column(
            "adoption_draft_expires_at", sa.DateTime(timezone=True), nullable=True
        ),
        sa.CheckConstraint(
            "length(turn_id) = 36 AND substr(turn_id, 1, 4) = 'atr_'",
            name="ck_agent_authoring_turns_id",
        ),
        sa.CheckConstraint(
            "status IN ('generating', 'completed', 'failed')",
            name="ck_agent_authoring_turns_status",
        ),
        sa.CheckConstraint(
            "decision IS NULL OR decision IN "
            "('pending', 'adopted', 'rejected', 'invalidated')",
            name="ck_agent_authoring_turns_decision",
        ),
        sa.CheckConstraint(
            "base_revision >= 0",
            name="ck_agent_authoring_turns_revision",
        ),
        sa.ForeignKeyConstraint(
            ["sid"],
            ["agent_authoring_sessions.sid"],
            name="fk_agent_authoring_turns_session",
            ondelete="CASCADE",
        ),
        sa.PrimaryKeyConstraint("turn_id"),
        sa.UniqueConstraint(
            "sid",
            "client_turn_id",
            name="uq_agent_authoring_turns_client",
        ),
        sa.UniqueConstraint("proposal_id", name="uq_agent_authoring_turns_proposal"),
    )
    op.create_index(
        "ix_agent_authoring_turns_session_status",
        "agent_authoring_turns",
        ["sid", "status"],
    )
    op.create_index(
        "ix_agent_authoring_turns_generation",
        "agent_authoring_turns",
        ["sid", "generation"],
    )


def downgrade() -> None:
    op.drop_index(
        "ix_agent_authoring_turns_generation",
        table_name="agent_authoring_turns",
    )
    op.drop_index(
        "ix_agent_authoring_turns_session_status",
        table_name="agent_authoring_turns",
    )
    op.drop_table("agent_authoring_turns")
    op.drop_index(
        "ix_agent_authoring_sessions_epoch",
        table_name="agent_authoring_sessions",
    )
    op.drop_index(
        "ix_agent_authoring_sessions_auth",
        table_name="agent_authoring_sessions",
    )
    op.drop_index(
        "ix_agent_authoring_sessions_owner_updated",
        table_name="agent_authoring_sessions",
    )
    op.drop_table("agent_authoring_sessions")
