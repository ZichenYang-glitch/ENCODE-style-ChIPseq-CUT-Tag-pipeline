"""Allow catalog proposals that were produced by the bound Peri model.

Revision ID: 20261126_40
Revises: 20261119_39

Offline registration still cannot submit this source. Only the classification
suggestion command inserts it.
"""

from __future__ import annotations

from alembic import op
import sqlalchemy as sa


revision = "20261126_40"
down_revision = "20261119_39"
branch_labels = None
depends_on = None


def upgrade() -> None:
    # Receipts reference these immutable rows. SQLite can append this column in
    # place without rebuilding the parent or replacing either audit trigger.
    op.execute(
        "ALTER TABLE catalog_proposals ADD COLUMN model_source VARCHAR(32) "
        "CONSTRAINT ck_catalog_proposals_model_source "
        "CHECK (model_source IS NULL OR "
        "(source = 'agent' AND model_source = 'peri_bound_model'))"
    )


def downgrade() -> None:
    incompatible = (
        op.get_bind()
        .execute(
            sa.text(
                "SELECT 1 FROM catalog_proposals WHERE model_source IS NOT NULL LIMIT 1"
            )
        )
        .first()
    )
    if incompatible:
        raise RuntimeError(
            "20261126_40 cannot downgrade while peri_bound_model proposals exist"
        )
    op.execute("ALTER TABLE catalog_proposals DROP COLUMN model_source")
