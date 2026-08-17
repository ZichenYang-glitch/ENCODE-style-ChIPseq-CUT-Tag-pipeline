"""Add LAN authentication accounts, sessions, and security audit ledger.

Revision ID: 20260818_14
Revises: 20260809_13
"""

from __future__ import annotations

from alembic import op
import sqlalchemy as sa


revision = "20260818_14"
down_revision = "20260809_13"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.create_table(
        "user_accounts",
        sa.Column("user_id", sa.String(length=36), nullable=False),
        sa.Column("username", sa.String(length=64), nullable=False),
        sa.Column("role", sa.String(length=16), nullable=False),
        sa.Column("status", sa.String(length=16), nullable=False),
        sa.Column("password_hash", sa.String(length=256), nullable=False),
        sa.Column("created_at", sa.DateTime(timezone=True), nullable=False),
        sa.Column("updated_at", sa.DateTime(timezone=True), nullable=False),
        sa.Column("password_changed_at", sa.DateTime(timezone=True), nullable=False),
        sa.CheckConstraint(
            "length(user_id) = 36 AND substr(user_id, 1, 4) = 'usr_'",
            name="ck_user_accounts_id",
        ),
        sa.CheckConstraint(
            "length(username) BETWEEN 3 AND 64",
            name="ck_user_accounts_username",
        ),
        sa.CheckConstraint(
            "role IN ('administrator', 'member')",
            name="ck_user_accounts_role",
        ),
        sa.CheckConstraint(
            "status IN ('enabled', 'disabled')",
            name="ck_user_accounts_status",
        ),
        sa.CheckConstraint(
            "updated_at >= created_at",
            name="ck_user_accounts_update_order",
        ),
        sa.CheckConstraint(
            "password_changed_at >= created_at AND password_changed_at <= updated_at",
            name="ck_user_accounts_password_change_order",
        ),
        sa.PrimaryKeyConstraint("user_id"),
        sa.UniqueConstraint("username", name="uq_user_accounts_username"),
    )
    op.create_table(
        "auth_sessions",
        sa.Column("session_digest", sa.String(length=64), nullable=False),
        sa.Column("csrf_digest", sa.String(length=64), nullable=False),
        sa.Column("user_id", sa.String(length=36), nullable=False),
        sa.Column("created_at", sa.DateTime(timezone=True), nullable=False),
        sa.Column("expires_at", sa.DateTime(timezone=True), nullable=False),
        sa.Column("revoked_at", sa.DateTime(timezone=True), nullable=True),
        sa.Column("revocation_reason", sa.String(length=32), nullable=True),
        sa.CheckConstraint(
            "length(session_digest) = 64 AND session_digest NOT GLOB '*[^0-9a-f]*'",
            name="ck_auth_sessions_id",
        ),
        sa.CheckConstraint(
            "length(csrf_digest) = 64 AND csrf_digest NOT GLOB '*[^0-9a-f]*'",
            name="ck_auth_sessions_csrf_digest",
        ),
        sa.CheckConstraint(
            "expires_at > created_at",
            name="ck_auth_sessions_expiry_order",
        ),
        sa.CheckConstraint(
            "revocation_reason IS NULL OR revocation_reason IN "
            "('logout', 'all_sessions', 'account_disabled', 'password_reset')",
            name="ck_auth_sessions_revocation_reason",
        ),
        sa.CheckConstraint(
            "(revoked_at IS NULL AND revocation_reason IS NULL) OR "
            "(revoked_at IS NOT NULL AND revocation_reason IS NOT NULL)",
            name="ck_auth_sessions_revocation_pair",
        ),
        sa.CheckConstraint(
            "revoked_at IS NULL OR revoked_at >= created_at",
            name="ck_auth_sessions_revocation_order",
        ),
        sa.ForeignKeyConstraint(
            ["user_id"],
            ["user_accounts.user_id"],
            name="fk_auth_sessions_user",
            ondelete="CASCADE",
        ),
        sa.PrimaryKeyConstraint("session_digest"),
    )
    op.create_index(
        "ix_auth_sessions_user",
        "auth_sessions",
        ["user_id"],
    )
    op.create_table(
        "security_audit_events",
        sa.Column("event_id", sa.String(length=37), nullable=False),
        sa.Column("occurred_at", sa.DateTime(timezone=True), nullable=False),
        sa.Column("action", sa.String(length=40), nullable=False),
        sa.Column("outcome", sa.String(length=16), nullable=False),
        sa.Column("actor_kind", sa.String(length=24), nullable=False),
        sa.Column("actor_user_id", sa.String(length=36), nullable=True),
        sa.Column("resource_kind", sa.String(length=24), nullable=True),
        sa.Column("resource_id", sa.String(length=88), nullable=True),
        sa.Column("reason_code", sa.String(length=40), nullable=True),
        sa.CheckConstraint(
            "length(event_id) = 37 AND substr(event_id, 1, 5) = 'aevt_'",
            name="ck_security_audit_events_id",
        ),
        sa.CheckConstraint(
            "action IN ('auth.login', 'auth.logout', 'run.create', 'run.start', "
            "'run.cancel', 'artifact.download', 'reference.register', "
            "'reference.enable', 'reference.disable', 'storage.register', "
            "'storage.archive', 'account.create', 'account.enable', "
            "'account.disable', 'account.password_reset', "
            "'account.sessions_revoke')",
            name="ck_security_audit_events_action",
        ),
        sa.CheckConstraint(
            "outcome IN ('succeeded', 'denied', 'failed')",
            name="ck_security_audit_events_outcome",
        ),
        sa.CheckConstraint(
            "actor_kind IN ('user', 'local_operator', 'unauthenticated')",
            name="ck_security_audit_events_actor_kind",
        ),
        sa.CheckConstraint(
            "(actor_kind = 'user' AND actor_user_id IS NOT NULL) OR "
            "(actor_kind != 'user' AND actor_user_id IS NULL)",
            name="ck_security_audit_events_actor_pair",
        ),
        sa.CheckConstraint(
            "resource_kind IS NULL OR resource_kind IN "
            "('account', 'run', 'artifact', 'reference', 'storage')",
            name="ck_security_audit_events_resource_kind",
        ),
        sa.CheckConstraint(
            "(resource_kind IS NULL AND resource_id IS NULL) OR "
            "(resource_kind IS NOT NULL AND resource_id IS NOT NULL)",
            name="ck_security_audit_events_resource_pair",
        ),
        sa.CheckConstraint(
            "reason_code IS NULL OR reason_code IN "
            "('INVALID_CREDENTIALS', 'LOGIN_RATE_LIMITED', "
            "'AUTHENTICATION_REQUIRED', 'SESSION_INVALID', 'CSRF_INVALID', "
            "'ADMINISTRATOR_REQUIRED', 'ACCOUNT_DISABLED', 'SETUP_REQUIRED', "
            "'OPERATION_CONFLICT', 'RESOURCE_NOT_FOUND', 'INTERNAL_FAILURE')",
            name="ck_security_audit_events_reason_code",
        ),
        sa.CheckConstraint(
            "(outcome = 'succeeded' AND reason_code IS NULL) OR "
            "(outcome != 'succeeded' AND reason_code IS NOT NULL)",
            name="ck_security_audit_events_outcome_reason_pair",
        ),
        sa.PrimaryKeyConstraint("event_id"),
    )
    op.create_index(
        "ix_security_audit_events_occurred",
        "security_audit_events",
        ["occurred_at", "event_id"],
    )


def downgrade() -> None:
    op.drop_index(
        "ix_security_audit_events_occurred",
        table_name="security_audit_events",
    )
    op.drop_table("security_audit_events")
    op.drop_index("ix_auth_sessions_user", table_name="auth_sessions")
    op.drop_table("auth_sessions")
    op.drop_table("user_accounts")
