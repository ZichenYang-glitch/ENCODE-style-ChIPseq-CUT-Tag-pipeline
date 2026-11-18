"""Add the Slice A experiment catalog and sample revision heads.

Revision ID: 20261118_17
Revises: 20261110_16

Upgrade inserts one head pointer per existing sample and does not rewrite
sample payloads, runs, publications, or legacy bindings.
"""

from __future__ import annotations

from alembic import op
import sqlalchemy as sa


revision = "20261118_17"
down_revision = "20261110_16"
branch_labels = None
depends_on = None

_FIELD_STATUS = "'known', 'unknown', 'not_applicable', 'not_detected', 'ambiguous'"
_LEGACY_PROJECT_ID = "prj_00000000000000000000000000000000"


def upgrade() -> None:
    op.create_table(
        "catalog_terms",
        sa.Column("term_id", sa.String(length=37), nullable=False),
        sa.Column("term_type", sa.String(length=32), nullable=False),
        sa.Column("created_at", sa.DateTime(timezone=True), nullable=False),
        sa.Column("archived_at", sa.DateTime(timezone=True), nullable=True),
        sa.CheckConstraint(
            "length(term_id) = 37 AND substr(term_id, 1, 5) = 'term_' "
            "AND substr(term_id, 6) NOT GLOB '*[^0-9a-f]*'",
            name="ck_catalog_terms_id",
        ),
        sa.CheckConstraint(
            "term_type IN ('species', 'tissue', 'cell_type', 'cell_line', "
            "'treatment', 'tag', 'assay', 'target', 'batch', 'replicate_role')",
            name="ck_catalog_terms_type",
        ),
        sa.CheckConstraint(
            "archived_at IS NULL OR archived_at >= created_at",
            name="ck_catalog_terms_archive_order",
        ),
        sa.PrimaryKeyConstraint("term_id"),
    )
    op.create_index(
        "ix_catalog_terms_active_created",
        "catalog_terms",
        ["archived_at", "created_at", "term_id"],
    )
    op.create_table(
        "catalog_term_revisions",
        sa.Column("term_revision_id", sa.String(length=37), nullable=False),
        sa.Column("term_id", sa.String(length=37), nullable=False),
        sa.Column("revision_number", sa.Integer(), nullable=False),
        sa.Column("display_name", sa.String(length=255), nullable=False),
        sa.Column("namespace", sa.String(length=255), nullable=True),
        sa.Column("accession", sa.String(length=255), nullable=True),
        sa.Column("aliases_json", sa.Text(), nullable=False),
        sa.Column("change_kind", sa.String(length=16), nullable=False),
        sa.Column("created_at", sa.DateTime(timezone=True), nullable=False),
        sa.CheckConstraint(
            "length(term_revision_id) = 37 AND substr(term_revision_id, 1, 5) = 'trmv_' "
            "AND substr(term_revision_id, 6) NOT GLOB '*[^0-9a-f]*'",
            name="ck_catalog_term_revisions_id",
        ),
        sa.CheckConstraint(
            "revision_number >= 1",
            name="ck_catalog_term_revisions_number",
        ),
        sa.CheckConstraint(
            "length(trim(display_name)) BETWEEN 1 AND 255",
            name="ck_catalog_term_revisions_display_name",
        ),
        sa.CheckConstraint(
            "change_kind IN ('initial', 'editorial')",
            name="ck_catalog_term_revisions_change_kind",
        ),
        sa.CheckConstraint(
            "(revision_number = 1 AND change_kind = 'initial') OR "
            "(revision_number > 1 AND change_kind = 'editorial')",
            name="ck_catalog_term_revisions_initial_kind",
        ),
        sa.ForeignKeyConstraint(
            ["term_id"],
            ["catalog_terms.term_id"],
            name="fk_catalog_term_revisions_term",
            ondelete="RESTRICT",
        ),
        sa.PrimaryKeyConstraint("term_revision_id"),
        sa.UniqueConstraint(
            "term_id",
            "revision_number",
            name="uq_catalog_term_revisions_number",
        ),
    )
    op.create_index(
        "ix_catalog_term_revisions_term_created",
        "catalog_term_revisions",
        ["term_id", "created_at", "term_revision_id"],
    )
    op.create_table(
        "experiments",
        sa.Column("experiment_id", sa.String(length=36), nullable=False),
        sa.Column("project_id", sa.String(length=36), nullable=False),
        sa.Column("stable_key", sa.String(length=255), nullable=False),
        sa.Column("created_at", sa.DateTime(timezone=True), nullable=False),
        sa.Column("archived_at", sa.DateTime(timezone=True), nullable=True),
        sa.CheckConstraint(
            "length(experiment_id) = 36 AND substr(experiment_id, 1, 4) = 'exp_' "
            "AND substr(experiment_id, 5) NOT GLOB '*[^0-9a-f]*'",
            name="ck_experiments_id",
        ),
        sa.CheckConstraint(
            "length(trim(stable_key)) BETWEEN 1 AND 255",
            name="ck_experiments_stable_key",
        ),
        sa.CheckConstraint(
            f"project_id != '{_LEGACY_PROJECT_ID}'",
            name="ck_experiments_not_legacy",
        ),
        sa.CheckConstraint(
            "archived_at IS NULL OR archived_at >= created_at",
            name="ck_experiments_archive_order",
        ),
        sa.ForeignKeyConstraint(
            ["project_id"],
            ["projects.project_id"],
            name="fk_experiments_project",
            ondelete="RESTRICT",
        ),
        sa.PrimaryKeyConstraint("experiment_id"),
        sa.UniqueConstraint(
            "project_id",
            "stable_key",
            name="uq_experiments_project_stable_key",
        ),
        sa.UniqueConstraint(
            "experiment_id",
            "project_id",
            name="uq_experiments_id_project",
        ),
    )
    op.create_index(
        "ix_experiments_project_created",
        "experiments",
        ["project_id", "created_at", "experiment_id"],
    )
    op.create_table(
        "experiment_revisions",
        sa.Column("experiment_revision_id", sa.String(length=37), nullable=False),
        sa.Column("experiment_id", sa.String(length=36), nullable=False),
        sa.Column("project_id", sa.String(length=36), nullable=False),
        sa.Column("revision_number", sa.Integer(), nullable=False),
        sa.Column("canonical_payload", sa.Text(), nullable=False),
        sa.Column("payload_digest_scheme", sa.String(length=80), nullable=False),
        sa.Column("payload_digest", sa.String(length=64), nullable=False),
        sa.Column("created_at", sa.DateTime(timezone=True), nullable=False),
        sa.CheckConstraint(
            "length(experiment_revision_id) = 37 "
            "AND substr(experiment_revision_id, 1, 5) = 'expr_' "
            "AND substr(experiment_revision_id, 6) NOT GLOB '*[^0-9a-f]*'",
            name="ck_experiment_revisions_id",
        ),
        sa.CheckConstraint(
            "revision_number >= 1",
            name="ck_experiment_revisions_number",
        ),
        sa.CheckConstraint(
            "payload_digest_scheme = 'sha256-framed-experiment-revision-payload-v1'",
            name="ck_experiment_revisions_digest_scheme",
        ),
        sa.CheckConstraint(
            "length(payload_digest) = 64 AND payload_digest NOT GLOB '*[^0-9a-f]*'",
            name="ck_experiment_revisions_digest",
        ),
        sa.ForeignKeyConstraint(
            ["experiment_id", "project_id"],
            ["experiments.experiment_id", "experiments.project_id"],
            name="fk_experiment_revisions_experiment",
            ondelete="RESTRICT",
        ),
        sa.PrimaryKeyConstraint("experiment_revision_id"),
        sa.UniqueConstraint(
            "experiment_id",
            "revision_number",
            name="uq_experiment_revisions_number",
        ),
    )
    op.create_index(
        "ix_experiment_revisions_experiment_created",
        "experiment_revisions",
        ["experiment_id", "created_at", "experiment_revision_id"],
    )
    op.create_table(
        "experiment_revision_heads",
        sa.Column("experiment_id", sa.String(length=36), nullable=False),
        sa.Column("project_id", sa.String(length=36), nullable=False),
        sa.Column("experiment_revision_id", sa.String(length=37), nullable=False),
        sa.Column("revision_number", sa.Integer(), nullable=False),
        sa.Column("created_at", sa.DateTime(timezone=True), nullable=False),
        sa.Column("revision_created_at", sa.DateTime(timezone=True), nullable=False),
        sa.CheckConstraint(
            "revision_number >= 1",
            name="ck_experiment_revision_heads_number",
        ),
        sa.ForeignKeyConstraint(
            ["experiment_id", "project_id"],
            ["experiments.experiment_id", "experiments.project_id"],
            name="fk_experiment_revision_heads_experiment",
            ondelete="RESTRICT",
        ),
        sa.ForeignKeyConstraint(
            ["experiment_revision_id"],
            ["experiment_revisions.experiment_revision_id"],
            name="fk_experiment_revision_heads_revision",
            ondelete="RESTRICT",
        ),
        sa.PrimaryKeyConstraint("experiment_id"),
        sa.UniqueConstraint(
            "experiment_revision_id",
            name="uq_experiment_revision_heads_revision",
        ),
    )
    op.create_index(
        "ix_experiment_revision_heads_project_created",
        "experiment_revision_heads",
        ["project_id", "created_at", "experiment_id"],
    )
    op.create_table(
        "experiment_revision_samples",
        sa.Column("experiment_revision_id", sa.String(length=37), nullable=False),
        sa.Column("ordinal", sa.Integer(), nullable=False),
        sa.Column("project_id", sa.String(length=36), nullable=False),
        sa.Column("sample_id", sa.String(length=36), nullable=False),
        sa.Column("sample_revision_id", sa.String(length=37), nullable=False),
        sa.Column("payload_digest", sa.String(length=64), nullable=False),
        sa.Column("replicate_group", sa.String(length=255), nullable=True),
        sa.Column("replicate_role_status", sa.String(length=32), nullable=False),
        sa.Column("replicate_role_term_id", sa.String(length=37), nullable=True),
        sa.Column(
            "replicate_role_term_revision_id",
            sa.String(length=37),
            nullable=True,
        ),
        sa.CheckConstraint(
            "ordinal >= 0",
            name="ck_experiment_revision_samples_ordinal",
        ),
        sa.CheckConstraint(
            "replicate_role_status IN (" + _FIELD_STATUS + ")",
            name="ck_experiment_revision_samples_role_status",
        ),
        sa.CheckConstraint(
            "(replicate_role_status = 'known' AND "
            "replicate_role_term_id IS NOT NULL AND "
            "replicate_role_term_revision_id IS NOT NULL) OR "
            "(replicate_role_status != 'known' AND "
            "replicate_role_term_id IS NULL AND "
            "replicate_role_term_revision_id IS NULL)",
            name="ck_experiment_revision_samples_role_term",
        ),
        sa.ForeignKeyConstraint(
            ["experiment_revision_id"],
            ["experiment_revisions.experiment_revision_id"],
            name="fk_experiment_revision_samples_revision",
            ondelete="RESTRICT",
        ),
        sa.ForeignKeyConstraint(
            ["project_id", "sample_revision_id", "payload_digest"],
            [
                "sample_revisions.project_id",
                "sample_revisions.sample_revision_id",
                "sample_revisions.payload_digest",
            ],
            name="fk_experiment_revision_samples_sample_revision",
            ondelete="RESTRICT",
        ),
        sa.ForeignKeyConstraint(
            ["replicate_role_term_revision_id"],
            ["catalog_term_revisions.term_revision_id"],
            name="fk_experiment_revision_samples_role_revision",
            ondelete="RESTRICT",
        ),
        sa.PrimaryKeyConstraint("experiment_revision_id", "ordinal"),
        sa.UniqueConstraint(
            "experiment_revision_id",
            "sample_revision_id",
            name="uq_experiment_revision_samples_revision",
        ),
    )
    op.create_index(
        "ix_experiment_revision_samples_sample_revision",
        "experiment_revision_samples",
        ["sample_revision_id", "experiment_revision_id"],
    )
    op.create_table(
        "experiment_catalog_projection",
        sa.Column("experiment_revision_id", sa.String(length=37), nullable=False),
        sa.Column("experiment_id", sa.String(length=36), nullable=False),
        sa.Column("project_id", sa.String(length=36), nullable=False),
        sa.Column("assay_status", sa.String(length=32), nullable=False),
        sa.Column("assay_term_id", sa.String(length=37), nullable=True),
        sa.Column("assay_term_revision_id", sa.String(length=37), nullable=True),
        sa.Column("target_status", sa.String(length=32), nullable=False),
        sa.Column("target_term_id", sa.String(length=37), nullable=True),
        sa.Column("target_term_revision_id", sa.String(length=37), nullable=True),
        sa.Column("batch_status", sa.String(length=32), nullable=False),
        sa.Column("batch_term_id", sa.String(length=37), nullable=True),
        sa.Column("batch_term_revision_id", sa.String(length=37), nullable=True),
        sa.CheckConstraint(
            "assay_status IN (" + _FIELD_STATUS + ")",
            name="ck_experiment_catalog_projection_assay_status",
        ),
        sa.CheckConstraint(
            "target_status IN (" + _FIELD_STATUS + ")",
            name="ck_experiment_catalog_projection_target_status",
        ),
        sa.CheckConstraint(
            "batch_status IN (" + _FIELD_STATUS + ")",
            name="ck_experiment_catalog_projection_batch_status",
        ),
        sa.ForeignKeyConstraint(
            ["experiment_revision_id"],
            ["experiment_revisions.experiment_revision_id"],
            name="fk_experiment_catalog_projection_revision",
            ondelete="RESTRICT",
        ),
        sa.ForeignKeyConstraint(
            ["experiment_id", "project_id"],
            ["experiments.experiment_id", "experiments.project_id"],
            name="fk_experiment_catalog_projection_experiment",
            ondelete="RESTRICT",
        ),
        sa.PrimaryKeyConstraint("experiment_revision_id"),
    )
    op.create_index(
        "ix_experiment_catalog_projection_assay",
        "experiment_catalog_projection",
        ["project_id", "assay_term_id", "experiment_revision_id"],
    )
    op.create_table(
        "sample_revision_heads",
        sa.Column("sample_id", sa.String(length=36), nullable=False),
        sa.Column("project_id", sa.String(length=36), nullable=False),
        sa.Column("sample_revision_id", sa.String(length=37), nullable=False),
        sa.Column("revision_number", sa.Integer(), nullable=False),
        sa.CheckConstraint(
            "revision_number >= 1",
            name="ck_sample_revision_heads_number",
        ),
        sa.ForeignKeyConstraint(
            ["project_id", "sample_id"],
            ["samples.project_id", "samples.sample_id"],
            name="fk_sample_revision_heads_sample",
            ondelete="RESTRICT",
        ),
        sa.ForeignKeyConstraint(
            ["sample_revision_id"],
            ["sample_revisions.sample_revision_id"],
            name="fk_sample_revision_heads_revision",
            ondelete="RESTRICT",
        ),
        sa.PrimaryKeyConstraint("sample_id"),
        sa.UniqueConstraint(
            "sample_revision_id",
            name="uq_sample_revision_heads_revision",
        ),
    )
    op.create_table(
        "sample_catalog_projection",
        sa.Column("sample_revision_id", sa.String(length=37), nullable=False),
        sa.Column("project_id", sa.String(length=36), nullable=False),
        sa.Column("sample_id", sa.String(length=36), nullable=False),
        sa.Column("catalog_state", sa.String(length=16), nullable=False),
        sa.Column("species_status", sa.String(length=32), nullable=True),
        sa.Column("species_term_id", sa.String(length=37), nullable=True),
        sa.Column("species_term_revision_id", sa.String(length=37), nullable=True),
        sa.Column("tissue_status", sa.String(length=32), nullable=True),
        sa.Column("tissue_term_id", sa.String(length=37), nullable=True),
        sa.Column("tissue_term_revision_id", sa.String(length=37), nullable=True),
        sa.Column("cell_type_status", sa.String(length=32), nullable=True),
        sa.Column("cell_line_status", sa.String(length=32), nullable=True),
        sa.Column("treatment_status", sa.String(length=32), nullable=True),
        sa.Column("tag_status", sa.String(length=32), nullable=True),
        sa.CheckConstraint(
            "catalog_state IN ('absent', 'present')",
            name="ck_sample_catalog_projection_state",
        ),
        sa.ForeignKeyConstraint(
            ["sample_revision_id"],
            ["sample_revisions.sample_revision_id"],
            name="fk_sample_catalog_projection_revision",
            ondelete="RESTRICT",
        ),
        sa.PrimaryKeyConstraint("sample_revision_id"),
    )
    op.create_index(
        "ix_sample_catalog_projection_project_state",
        "sample_catalog_projection",
        ["project_id", "catalog_state", "sample_revision_id"],
    )
    op.create_table(
        "sample_catalog_values",
        sa.Column("sample_revision_id", sa.String(length=37), nullable=False),
        sa.Column("field_name", sa.String(length=32), nullable=False),
        sa.Column("term_id", sa.String(length=37), nullable=False),
        sa.Column("term_revision_id", sa.String(length=37), nullable=False),
        sa.CheckConstraint(
            "field_name IN ('cell_type', 'cell_line', 'treatment', 'tag')",
            name="ck_sample_catalog_values_field",
        ),
        sa.ForeignKeyConstraint(
            ["sample_revision_id"],
            ["sample_revisions.sample_revision_id"],
            name="fk_sample_catalog_values_revision",
            ondelete="RESTRICT",
        ),
        sa.ForeignKeyConstraint(
            ["term_revision_id"],
            ["catalog_term_revisions.term_revision_id"],
            name="fk_sample_catalog_values_term_revision",
            ondelete="RESTRICT",
        ),
        sa.PrimaryKeyConstraint("sample_revision_id", "field_name", "term_id"),
    )
    op.create_index(
        "ix_sample_catalog_values_field_term",
        "sample_catalog_values",
        ["field_name", "term_id", "sample_revision_id"],
    )
    op.create_table(
        "catalog_write_receipts",
        sa.Column("actor_user_id", sa.String(length=36), nullable=False),
        sa.Column("idempotency_key", sa.String(length=128), nullable=False),
        sa.Column("operation", sa.String(length=64), nullable=False),
        sa.Column("request_sha256", sa.String(length=64), nullable=False),
        sa.Column("subject_id", sa.String(length=37), nullable=False),
        sa.Column("result_revision_id", sa.String(length=37), nullable=True),
        sa.Column("result_archived_at", sa.DateTime(timezone=True), nullable=True),
        sa.Column("created_at", sa.DateTime(timezone=True), nullable=False),
        sa.CheckConstraint(
            "length(idempotency_key) BETWEEN 1 AND 128",
            name="ck_catalog_write_receipts_key",
        ),
        sa.CheckConstraint(
            "length(request_sha256) = 64 AND request_sha256 NOT GLOB '*[^0-9a-f]*'",
            name="ck_catalog_write_receipts_digest",
        ),
        sa.ForeignKeyConstraint(
            ["actor_user_id"],
            ["user_accounts.user_id"],
            name="fk_catalog_write_receipts_actor",
            ondelete="RESTRICT",
        ),
        sa.PrimaryKeyConstraint("actor_user_id", "idempotency_key"),
    )
    op.create_index(
        "ix_catalog_write_receipts_subject",
        "catalog_write_receipts",
        ["operation", "subject_id"],
    )
    _create_immutability_triggers()
    op.execute(
        """
        INSERT INTO sample_revision_heads (
            sample_id,
            project_id,
            sample_revision_id,
            revision_number
        )
        SELECT
            revision.sample_id,
            revision.project_id,
            revision.sample_revision_id,
            revision.revision_number
        FROM sample_revisions AS revision
        WHERE revision.revision_number = (
            SELECT MAX(latest.revision_number)
            FROM sample_revisions AS latest
            WHERE latest.sample_id = revision.sample_id
        )
        """
    )


def downgrade() -> None:
    for name in (
        "trg_catalog_write_receipts_no_delete",
        "trg_catalog_write_receipts_no_update",
        "trg_experiment_revisions_no_delete",
        "trg_experiment_revisions_no_update",
        "trg_experiments_no_delete",
        "trg_experiments_identity_no_update",
        "trg_catalog_term_revisions_no_delete",
        "trg_catalog_term_revisions_no_update",
        "trg_catalog_terms_no_delete",
        "trg_catalog_terms_identity_no_update",
    ):
        op.execute(f"DROP TRIGGER IF EXISTS {name}")
    op.drop_index(
        "ix_catalog_write_receipts_subject",
        table_name="catalog_write_receipts",
    )
    op.drop_table("catalog_write_receipts")
    op.drop_index(
        "ix_sample_catalog_values_field_term",
        table_name="sample_catalog_values",
    )
    op.drop_table("sample_catalog_values")
    op.drop_index(
        "ix_sample_catalog_projection_project_state",
        table_name="sample_catalog_projection",
    )
    op.drop_table("sample_catalog_projection")
    op.drop_table("sample_revision_heads")
    op.drop_index(
        "ix_experiment_catalog_projection_assay",
        table_name="experiment_catalog_projection",
    )
    op.drop_table("experiment_catalog_projection")
    op.drop_index(
        "ix_experiment_revision_samples_sample_revision",
        table_name="experiment_revision_samples",
    )
    op.drop_table("experiment_revision_samples")
    op.drop_index(
        "ix_experiment_revision_heads_project_created",
        table_name="experiment_revision_heads",
    )
    op.drop_table("experiment_revision_heads")
    op.drop_index(
        "ix_experiment_revisions_experiment_created",
        table_name="experiment_revisions",
    )
    op.drop_table("experiment_revisions")
    op.drop_index("ix_experiments_project_created", table_name="experiments")
    op.drop_table("experiments")
    op.drop_index(
        "ix_catalog_term_revisions_term_created",
        table_name="catalog_term_revisions",
    )
    op.drop_table("catalog_term_revisions")
    op.drop_index("ix_catalog_terms_active_created", table_name="catalog_terms")
    op.drop_table("catalog_terms")


def _create_immutability_triggers() -> None:
    op.execute(
        """
        CREATE TRIGGER trg_catalog_terms_identity_no_update
        BEFORE UPDATE ON catalog_terms
        WHEN OLD.term_id != NEW.term_id
          OR OLD.term_type != NEW.term_type
          OR OLD.created_at != NEW.created_at
        BEGIN
            SELECT RAISE(ABORT, 'catalog term identity is immutable');
        END
        """
    )
    op.execute(
        """
        CREATE TRIGGER trg_catalog_terms_no_delete
        BEFORE DELETE ON catalog_terms
        BEGIN
            SELECT RAISE(ABORT, 'catalog term cannot be deleted');
        END
        """
    )
    op.execute(
        """
        CREATE TRIGGER trg_catalog_term_revisions_no_update
        BEFORE UPDATE ON catalog_term_revisions
        BEGIN
            SELECT RAISE(ABORT, 'catalog term revision is immutable');
        END
        """
    )
    op.execute(
        """
        CREATE TRIGGER trg_catalog_term_revisions_no_delete
        BEFORE DELETE ON catalog_term_revisions
        BEGIN
            SELECT RAISE(ABORT, 'catalog term revision is immutable');
        END
        """
    )
    op.execute(
        """
        CREATE TRIGGER trg_experiments_identity_no_update
        BEFORE UPDATE ON experiments
        WHEN OLD.experiment_id != NEW.experiment_id
          OR OLD.project_id != NEW.project_id
          OR OLD.stable_key != NEW.stable_key
          OR OLD.created_at != NEW.created_at
        BEGIN
            SELECT RAISE(ABORT, 'experiment identity is immutable');
        END
        """
    )
    op.execute(
        """
        CREATE TRIGGER trg_experiments_no_delete
        BEFORE DELETE ON experiments
        BEGIN
            SELECT RAISE(ABORT, 'experiment cannot be deleted');
        END
        """
    )
    op.execute(
        """
        CREATE TRIGGER trg_experiment_revisions_no_update
        BEFORE UPDATE ON experiment_revisions
        BEGIN
            SELECT RAISE(ABORT, 'experiment revision is immutable');
        END
        """
    )
    op.execute(
        """
        CREATE TRIGGER trg_experiment_revisions_no_delete
        BEFORE DELETE ON experiment_revisions
        BEGIN
            SELECT RAISE(ABORT, 'experiment revision is immutable');
        END
        """
    )
    op.execute(
        """
        CREATE TRIGGER trg_catalog_write_receipts_no_update
        BEFORE UPDATE ON catalog_write_receipts
        BEGIN
            SELECT RAISE(ABORT, 'catalog write receipt is immutable');
        END
        """
    )
    op.execute(
        """
        CREATE TRIGGER trg_catalog_write_receipts_no_delete
        BEFORE DELETE ON catalog_write_receipts
        BEGIN
            SELECT RAISE(ABORT, 'catalog write receipt is immutable');
        END
        """
    )
