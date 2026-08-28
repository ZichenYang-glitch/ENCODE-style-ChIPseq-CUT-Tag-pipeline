"""Add private terminal-email preferences and immutable run requesters.

Revision ID: 20260827_15
Revises: 20260818_14
"""

from __future__ import annotations

from alembic import op
import sqlalchemy as sa


revision = "20260827_15"
down_revision = "20260818_14"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.add_column(
        "user_accounts",
        sa.Column("notification_email", sa.String(length=254), nullable=True),
    )
    op.add_column(
        "user_accounts",
        sa.Column(
            "terminal_email_enabled",
            sa.Boolean(),
            nullable=False,
            server_default=sa.true(),
        ),
    )
    connection = op.get_bind()
    if connection.dialect.name == "sqlite":
        op.execute(
            "ALTER TABLE runs ADD COLUMN requested_by_user_id VARCHAR(36) "
            "CONSTRAINT fk_runs_requester_user REFERENCES user_accounts(user_id) "
            "ON DELETE RESTRICT"
        )
    else:
        op.add_column(
            "runs",
            sa.Column("requested_by_user_id", sa.String(length=36), nullable=True),
        )
        op.create_foreign_key(
            "fk_runs_requester_user",
            "runs",
            "user_accounts",
            ["requested_by_user_id"],
            ["user_id"],
            ondelete="RESTRICT",
        )
    op.create_index(
        "ix_runs_requested_by_user",
        "runs",
        ["requested_by_user_id"],
    )


def downgrade() -> None:
    op.drop_index("ix_runs_requested_by_user", table_name="runs")
    connection = op.get_bind()
    if connection.dialect.name == "sqlite":
        op.execute("ALTER TABLE runs DROP COLUMN requested_by_user_id")
    else:
        op.drop_constraint(
            "fk_runs_requester_user",
            "runs",
            type_="foreignkey",
        )
        op.drop_column("runs", "requested_by_user_id")
    op.drop_column("user_accounts", "terminal_email_enabled")
    op.drop_column("user_accounts", "notification_email")
