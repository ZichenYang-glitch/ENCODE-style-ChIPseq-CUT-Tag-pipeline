"""Store offline catalog proposals and the one adoption receipt for each.

Revision ID: 20261119_32
Revises: 20261119_18

The existing 20261119_18 snapshot-binding revision stays in place. This
revision is the single next head.
"""

from __future__ import annotations

from alembic import op
import sqlalchemy as sa


revision = "20261119_32"
down_revision = "20261119_18"
branch_labels = None
depends_on = None


def upgrade() -> None:
    op.create_table(
        "catalog_curation_versions",
        sa.Column("version_id", sa.Integer(), nullable=False),
        sa.Column("rule_version", sa.String(length=64), nullable=False),
        sa.Column("vocabulary_version", sa.String(length=64), nullable=False),
        sa.Column("projection_version", sa.String(length=64), nullable=False),
        sa.CheckConstraint(
            "version_id = 1",
            name="ck_catalog_curation_versions_singleton",
        ),
        sa.CheckConstraint(
            "length(rule_version) BETWEEN 1 AND 64",
            name="ck_catalog_curation_versions_rule",
        ),
        sa.CheckConstraint(
            "length(vocabulary_version) BETWEEN 1 AND 64",
            name="ck_catalog_curation_versions_vocabulary",
        ),
        sa.CheckConstraint(
            "length(projection_version) BETWEEN 1 AND 64",
            name="ck_catalog_curation_versions_projection",
        ),
        sa.PrimaryKeyConstraint("version_id"),
    )
    op.execute(
        "INSERT INTO catalog_curation_versions "
        "(version_id, rule_version, vocabulary_version, projection_version) "
        "VALUES (1, 'm2-catalog-rules/1', 'terms-snapshot/2026-11-18', "
        "'hw-m2-catalog/1')"
    )
    op.create_table(
        "catalog_proposals",
        sa.Column("proposal_id", sa.String(length=38), nullable=False),
        sa.Column("project_id", sa.String(length=36), nullable=False),
        sa.Column("created_by_user_id", sa.String(length=36), nullable=False),
        sa.Column("source", sa.String(length=32), nullable=False),
        sa.Column("rule_version", sa.String(length=64), nullable=False),
        sa.Column("vocabulary_version", sa.String(length=64), nullable=False),
        sa.Column("projection_version", sa.String(length=64), nullable=False),
        sa.Column("proposal_digest", sa.String(length=64), nullable=False),
        sa.Column("document_json", sa.Text(), nullable=False),
        sa.Column("response_json", sa.Text(), nullable=False),
        sa.Column("register_actor_user_id", sa.String(length=36), nullable=False),
        sa.Column("register_idempotency_key", sa.String(length=128), nullable=False),
        sa.Column("created_at", sa.DateTime(timezone=True), nullable=False),
        sa.Column("expires_at", sa.DateTime(timezone=True), nullable=False),
        sa.Column("retention_until", sa.DateTime(timezone=True), nullable=False),
        sa.CheckConstraint(
            "length(proposal_id) = 38 AND substr(proposal_id, 1, 6) = 'cprop_' "
            "AND substr(proposal_id, 7) NOT GLOB '*[^0-9a-f]*'",
            name="ck_catalog_proposals_id",
        ),
        sa.CheckConstraint(
            "source IN ('offline_suggestion', 'agent', 'manual')",
            name="ck_catalog_proposals_source",
        ),
        sa.CheckConstraint(
            "length(proposal_digest) = 64 AND proposal_digest NOT GLOB '*[^0-9a-f]*'",
            name="ck_catalog_proposals_digest",
        ),
        sa.CheckConstraint(
            "length(register_idempotency_key) BETWEEN 1 AND 128",
            name="ck_catalog_proposals_key",
        ),
        sa.ForeignKeyConstraint(
            ["project_id"],
            ["projects.project_id"],
            name="fk_catalog_proposals_project",
            ondelete="RESTRICT",
        ),
        sa.ForeignKeyConstraint(
            ["created_by_user_id"],
            ["user_accounts.user_id"],
            name="fk_catalog_proposals_creator",
            ondelete="RESTRICT",
        ),
        sa.ForeignKeyConstraint(
            ["register_actor_user_id"],
            ["user_accounts.user_id"],
            name="fk_catalog_proposals_actor",
            ondelete="RESTRICT",
        ),
        sa.PrimaryKeyConstraint("proposal_id"),
        sa.UniqueConstraint(
            "register_actor_user_id",
            "register_idempotency_key",
            name="uq_catalog_proposals_register_key",
        ),
    )
    op.create_index(
        "ix_catalog_proposals_project_created",
        "catalog_proposals",
        ["project_id", "created_at", "proposal_id"],
    )
    op.create_table(
        "catalog_adoption_receipts",
        sa.Column("receipt_id", sa.String(length=39), nullable=False),
        sa.Column("proposal_id", sa.String(length=38), nullable=False),
        sa.Column("project_id", sa.String(length=36), nullable=False),
        sa.Column("actor_user_id", sa.String(length=36), nullable=False),
        sa.Column("idempotency_key", sa.String(length=128), nullable=False),
        sa.Column("response_json", sa.Text(), nullable=False),
        sa.Column("created_at", sa.DateTime(timezone=True), nullable=False),
        sa.CheckConstraint(
            "length(receipt_id) = 39 AND substr(receipt_id, 1, 7) = 'cadopt_' "
            "AND substr(receipt_id, 8) NOT GLOB '*[^0-9a-f]*'",
            name="ck_catalog_adoption_receipts_id",
        ),
        sa.CheckConstraint(
            "length(idempotency_key) BETWEEN 1 AND 128",
            name="ck_catalog_adoption_receipts_key",
        ),
        sa.ForeignKeyConstraint(
            ["proposal_id"],
            ["catalog_proposals.proposal_id"],
            name="fk_catalog_adoption_receipts_proposal",
            ondelete="RESTRICT",
        ),
        sa.ForeignKeyConstraint(
            ["project_id"],
            ["projects.project_id"],
            name="fk_catalog_adoption_receipts_project",
            ondelete="RESTRICT",
        ),
        sa.ForeignKeyConstraint(
            ["actor_user_id"],
            ["user_accounts.user_id"],
            name="fk_catalog_adoption_receipts_actor",
            ondelete="RESTRICT",
        ),
        sa.PrimaryKeyConstraint("receipt_id"),
        sa.UniqueConstraint(
            "proposal_id",
            name="uq_catalog_adoption_receipts_proposal",
        ),
        sa.UniqueConstraint(
            "actor_user_id",
            "idempotency_key",
            name="uq_catalog_adoption_receipts_actor_key",
        ),
    )
    op.create_index(
        "ix_catalog_adoption_receipts_project_created",
        "catalog_adoption_receipts",
        ["project_id", "created_at", "receipt_id"],
    )
    op.execute(
        """
        CREATE TRIGGER trg_catalog_proposals_no_update
        BEFORE UPDATE ON catalog_proposals
        BEGIN
            SELECT RAISE(ABORT, 'catalog proposal is immutable');
        END
        """
    )
    op.execute(
        """
        CREATE TRIGGER trg_catalog_proposals_no_delete
        BEFORE DELETE ON catalog_proposals
        BEGIN
            SELECT RAISE(ABORT, 'catalog proposal is immutable');
        END
        """
    )
    op.execute(
        """
        CREATE TRIGGER trg_catalog_adoption_receipts_no_update
        BEFORE UPDATE ON catalog_adoption_receipts
        BEGIN
            SELECT RAISE(ABORT, 'catalog adoption receipt is immutable');
        END
        """
    )
    op.execute(
        """
        CREATE TRIGGER trg_catalog_adoption_receipts_no_delete
        BEFORE DELETE ON catalog_adoption_receipts
        BEGIN
            SELECT RAISE(ABORT, 'catalog adoption receipt is immutable');
        END
        """
    )


def downgrade() -> None:
    for name in (
        "trg_catalog_adoption_receipts_no_delete",
        "trg_catalog_adoption_receipts_no_update",
        "trg_catalog_proposals_no_delete",
        "trg_catalog_proposals_no_update",
    ):
        op.execute(f"DROP TRIGGER IF EXISTS {name}")
    op.drop_index(
        "ix_catalog_adoption_receipts_project_created",
        table_name="catalog_adoption_receipts",
    )
    op.drop_table("catalog_adoption_receipts")
    op.drop_index(
        "ix_catalog_proposals_project_created",
        table_name="catalog_proposals",
    )
    op.drop_table("catalog_proposals")
    op.drop_table("catalog_curation_versions")
