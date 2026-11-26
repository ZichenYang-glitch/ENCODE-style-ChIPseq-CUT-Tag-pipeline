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

_MODEL_SOURCE = (
    "source IN ('offline_suggestion', 'agent', 'manual', 'peri_bound_model')"
)
_OFFLINE_SOURCES = "source IN ('offline_suggestion', 'agent', 'manual')"


def upgrade() -> None:
    with op.batch_alter_table("catalog_proposals") as batch_op:
        batch_op.drop_constraint("ck_catalog_proposals_source", type_="check")
        batch_op.create_check_constraint(
            "ck_catalog_proposals_source",
            _MODEL_SOURCE,
        )


def downgrade() -> None:
    connection = op.get_bind()
    incompatible = connection.execute(
        sa.text(
            "SELECT count(*) FROM catalog_proposals WHERE source = 'peri_bound_model'"
        )
    ).scalar()
    if incompatible:
        raise RuntimeError(
            "20261126_40 cannot downgrade while peri_bound_model proposals exist"
        )
    with op.batch_alter_table("catalog_proposals") as batch_op:
        batch_op.drop_constraint("ck_catalog_proposals_source", type_="check")
        batch_op.create_check_constraint(
            "ck_catalog_proposals_source",
            _OFFLINE_SOURCES,
        )
