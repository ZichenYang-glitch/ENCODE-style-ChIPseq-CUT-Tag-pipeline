"""SQLAlchemy rows kept behind the persistence repository boundary."""

from __future__ import annotations

from datetime import datetime
from typing import Any

from sqlalchemy import (
    CheckConstraint,
    DateTime,
    ForeignKey,
    ForeignKeyConstraint,
    Index,
    Integer,
    JSON,
    String,
    Text,
    UniqueConstraint,
)
from sqlalchemy.orm import DeclarativeBase, Mapped, mapped_column


class Base(DeclarativeBase):
    """Declarative base for platform persistence rows."""


class ProjectRow(Base):
    __tablename__ = "projects"
    __table_args__ = (
        CheckConstraint(
            "length(project_id) = 36 AND substr(project_id, 1, 4) = 'prj_'",
            name="ck_projects_id",
        ),
        CheckConstraint(
            "length(trim(display_name)) BETWEEN 1 AND 255",
            name="ck_projects_display_name",
        ),
        CheckConstraint(
            "kind IN ('user', 'system')",
            name="ck_projects_kind",
        ),
        CheckConstraint(
            "(kind = 'system' AND "
            "project_id = 'prj_00000000000000000000000000000000' AND "
            "display_name = 'Legacy Project' AND archived_at IS NULL) OR "
            "(kind = 'user' AND "
            "project_id != 'prj_00000000000000000000000000000000')",
            name="ck_projects_legacy_identity",
        ),
        CheckConstraint(
            "archived_at IS NULL OR archived_at >= created_at",
            name="ck_projects_archive_order",
        ),
        Index(
            "ix_projects_archived_created",
            "archived_at",
            "created_at",
            "project_id",
        ),
    )

    project_id: Mapped[str] = mapped_column(String(36), primary_key=True)
    display_name: Mapped[str] = mapped_column(String(255), nullable=False)
    kind: Mapped[str] = mapped_column(String(32), nullable=False)
    created_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), nullable=False
    )
    archived_at: Mapped[datetime | None] = mapped_column(DateTime(timezone=True))


class SampleRow(Base):
    __tablename__ = "samples"
    __table_args__ = (
        CheckConstraint(
            "length(sample_id) = 36 AND substr(sample_id, 1, 4) = 'smp_'",
            name="ck_samples_id",
        ),
        CheckConstraint(
            "length(trim(stable_key)) BETWEEN 1 AND 255",
            name="ck_samples_stable_key",
        ),
        CheckConstraint(
            "project_id != 'prj_00000000000000000000000000000000'",
            name="ck_samples_not_legacy",
        ),
        ForeignKeyConstraint(
            ["project_id"],
            ["projects.project_id"],
            name="fk_samples_project",
            ondelete="RESTRICT",
        ),
        UniqueConstraint(
            "project_id",
            "sample_id",
            name="uq_samples_project_sample",
        ),
        UniqueConstraint(
            "project_id",
            "stable_key",
            name="uq_samples_project_stable_key",
        ),
        Index(
            "ix_samples_project_created",
            "project_id",
            "created_at",
            "sample_id",
        ),
    )

    sample_id: Mapped[str] = mapped_column(String(36), primary_key=True)
    project_id: Mapped[str] = mapped_column(String(36), nullable=False)
    stable_key: Mapped[str] = mapped_column(String(255), nullable=False)
    created_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), nullable=False
    )


class SampleRevisionRow(Base):
    __tablename__ = "sample_revisions"
    __table_args__ = (
        CheckConstraint(
            "length(sample_revision_id) = 37 "
            "AND substr(sample_revision_id, 1, 5) = 'smpr_'",
            name="ck_sample_revisions_id",
        ),
        CheckConstraint(
            "revision_number >= 1",
            name="ck_sample_revisions_positive_revision",
        ),
        CheckConstraint(
            "payload_digest_scheme = 'sha256-framed-sample-revision-payload-v1'",
            name="ck_sample_revisions_digest_scheme",
        ),
        CheckConstraint(
            "length(payload_digest) = 64",
            name="ck_sample_revisions_digest_length",
        ),
        ForeignKeyConstraint(
            ["project_id", "sample_id"],
            ["samples.project_id", "samples.sample_id"],
            name="fk_sample_revisions_sample",
            ondelete="RESTRICT",
        ),
        UniqueConstraint(
            "project_id",
            "sample_revision_id",
            name="uq_sample_revisions_project_revision",
        ),
        UniqueConstraint(
            "project_id",
            "sample_revision_id",
            "payload_digest",
            name="uq_sample_revisions_project_revision_digest",
        ),
        UniqueConstraint(
            "sample_id",
            "revision_number",
            name="uq_sample_revisions_sample_number",
        ),
        Index(
            "ix_sample_revisions_sample_revision",
            "sample_id",
            "revision_number",
        ),
        Index(
            "ix_sample_revisions_project_created",
            "project_id",
            "created_at",
            "sample_revision_id",
        ),
    )

    sample_revision_id: Mapped[str] = mapped_column(String(37), primary_key=True)
    project_id: Mapped[str] = mapped_column(String(36), nullable=False)
    sample_id: Mapped[str] = mapped_column(String(36), nullable=False)
    revision_number: Mapped[int] = mapped_column(Integer, nullable=False)
    canonical_payload: Mapped[str] = mapped_column(Text, nullable=False)
    payload_digest_scheme: Mapped[str] = mapped_column(String(64), nullable=False)
    payload_digest: Mapped[str] = mapped_column(String(64), nullable=False)
    created_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), nullable=False
    )


class SnapshotProjectBindingRow(Base):
    __tablename__ = "snapshot_project_bindings"
    __table_args__ = (
        CheckConstraint(
            "binding_mode IN ('legacy_v1', 'bound_v1')",
            name="ck_snapshot_project_bindings_mode",
        ),
        CheckConstraint(
            "provenance IN ('resolved', 'unresolved')",
            name="ck_snapshot_project_bindings_provenance",
        ),
        CheckConstraint(
            "(binding_mode = 'legacy_v1' AND provenance = 'unresolved' AND "
            "project_id = 'prj_00000000000000000000000000000000') OR "
            "(binding_mode = 'bound_v1' AND provenance = 'resolved' AND "
            "project_id != 'prj_00000000000000000000000000000000')",
            name="ck_snapshot_project_bindings_legacy_project",
        ),
        CheckConstraint(
            "binding_digest_scheme = 'sha256-framed-project-sample-binding-v1'",
            name="ck_snapshot_project_bindings_digest_scheme",
        ),
        CheckConstraint(
            "length(binding_digest) = 64",
            name="ck_snapshot_project_bindings_digest_length",
        ),
        CheckConstraint(
            "length(workflow_inputs_digest) = 64",
            name="ck_snapshot_project_bindings_workflow_digest_length",
        ),
        ForeignKeyConstraint(
            ["snapshot_id"],
            ["validated_input_snapshots.snapshot_id"],
            name="fk_snapshot_project_bindings_snapshot",
            ondelete="RESTRICT",
        ),
        ForeignKeyConstraint(
            ["project_id"],
            ["projects.project_id"],
            name="fk_snapshot_project_bindings_project",
            ondelete="RESTRICT",
        ),
        UniqueConstraint(
            "snapshot_id",
            "project_id",
            name="uq_snapshot_project_bindings_project",
        ),
        UniqueConstraint(
            "snapshot_id",
            "project_id",
            "binding_digest",
            name="uq_snapshot_project_bindings_input_evidence",
        ),
        Index(
            "ix_snapshot_project_bindings_project",
            "project_id",
            "created_at",
            "snapshot_id",
        ),
    )

    snapshot_id: Mapped[str] = mapped_column(String(64), primary_key=True)
    project_id: Mapped[str] = mapped_column(String(36), nullable=False)
    binding_mode: Mapped[str] = mapped_column(String(32), nullable=False)
    provenance: Mapped[str] = mapped_column(String(32), nullable=False)
    workflow_inputs_digest: Mapped[str] = mapped_column(String(64), nullable=False)
    binding_digest_scheme: Mapped[str] = mapped_column(String(64), nullable=False)
    binding_digest: Mapped[str] = mapped_column(String(64), nullable=False)
    created_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), nullable=False
    )


class SnapshotSampleRevisionRow(Base):
    __tablename__ = "snapshot_sample_revisions"
    __table_args__ = (
        CheckConstraint(
            "ordinal >= 0",
            name="ck_snapshot_sample_revisions_ordinal",
        ),
        ForeignKeyConstraint(
            ["snapshot_id", "project_id"],
            [
                "snapshot_project_bindings.snapshot_id",
                "snapshot_project_bindings.project_id",
            ],
            name="fk_snapshot_sample_revisions_binding",
            ondelete="RESTRICT",
        ),
        ForeignKeyConstraint(
            ["project_id", "sample_revision_id", "payload_digest"],
            [
                "sample_revisions.project_id",
                "sample_revisions.sample_revision_id",
                "sample_revisions.payload_digest",
            ],
            name="fk_snapshot_sample_revisions_revision",
            ondelete="RESTRICT",
        ),
        UniqueConstraint(
            "snapshot_id",
            "sample_revision_id",
            name="uq_snapshot_sample_revisions_revision",
        ),
        Index(
            "ix_snapshot_sample_revisions_revision",
            "sample_revision_id",
            "snapshot_id",
        ),
    )

    snapshot_id: Mapped[str] = mapped_column(String(64), primary_key=True)
    project_id: Mapped[str] = mapped_column(String(36), nullable=False)
    sample_revision_id: Mapped[str] = mapped_column(String(37), nullable=False)
    payload_digest: Mapped[str] = mapped_column(String(64), nullable=False)
    ordinal: Mapped[int] = mapped_column(Integer, primary_key=True)


class RunProjectBindingRow(Base):
    __tablename__ = "run_project_bindings"
    __table_args__ = (
        CheckConstraint(
            "binding_mode IN ('legacy_v1', 'bound_v1')",
            name="ck_run_project_bindings_mode",
        ),
        CheckConstraint(
            "provenance IN ('resolved', 'unresolved')",
            name="ck_run_project_bindings_provenance",
        ),
        CheckConstraint(
            "(binding_mode = 'legacy_v1' AND provenance = 'unresolved' AND "
            "project_id = 'prj_00000000000000000000000000000000') OR "
            "(binding_mode = 'bound_v1' AND provenance = 'resolved' AND "
            "project_id != 'prj_00000000000000000000000000000000')",
            name="ck_run_project_bindings_legacy_project",
        ),
        CheckConstraint(
            "binding_digest_scheme = 'sha256-framed-project-sample-binding-v1'",
            name="ck_run_project_bindings_digest_scheme",
        ),
        CheckConstraint(
            "length(binding_digest) = 64",
            name="ck_run_project_bindings_digest_length",
        ),
        CheckConstraint(
            "length(workflow_inputs_digest) = 64",
            name="ck_run_project_bindings_workflow_digest_length",
        ),
        ForeignKeyConstraint(
            ["run_id"],
            ["runs.run_id"],
            name="fk_run_project_bindings_run",
            ondelete="CASCADE",
        ),
        ForeignKeyConstraint(
            ["project_id"],
            ["projects.project_id"],
            name="fk_run_project_bindings_project",
            ondelete="RESTRICT",
        ),
        UniqueConstraint(
            "run_id",
            "project_id",
            name="uq_run_project_bindings_project",
        ),
        UniqueConstraint(
            "run_id",
            "project_id",
            "binding_digest",
            name="uq_run_project_bindings_input_evidence",
        ),
        Index(
            "ix_run_project_bindings_project",
            "project_id",
            "created_at",
            "run_id",
        ),
    )

    run_id: Mapped[str] = mapped_column(String(128), primary_key=True)
    project_id: Mapped[str] = mapped_column(String(36), nullable=False)
    binding_mode: Mapped[str] = mapped_column(String(32), nullable=False)
    provenance: Mapped[str] = mapped_column(String(32), nullable=False)
    workflow_inputs_digest: Mapped[str] = mapped_column(String(64), nullable=False)
    binding_digest_scheme: Mapped[str] = mapped_column(String(64), nullable=False)
    binding_digest: Mapped[str] = mapped_column(String(64), nullable=False)
    created_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), nullable=False
    )


class RunSampleRow(Base):
    __tablename__ = "run_samples"
    __table_args__ = (
        CheckConstraint("ordinal >= 0", name="ck_run_samples_ordinal"),
        ForeignKeyConstraint(
            ["run_id", "project_id"],
            [
                "run_project_bindings.run_id",
                "run_project_bindings.project_id",
            ],
            name="fk_run_samples_binding",
            ondelete="CASCADE",
        ),
        ForeignKeyConstraint(
            ["project_id", "sample_revision_id", "payload_digest"],
            [
                "sample_revisions.project_id",
                "sample_revisions.sample_revision_id",
                "sample_revisions.payload_digest",
            ],
            name="fk_run_samples_revision",
            ondelete="RESTRICT",
        ),
        UniqueConstraint(
            "run_id",
            "sample_revision_id",
            name="uq_run_samples_revision",
        ),
        Index(
            "ix_run_samples_revision",
            "sample_revision_id",
            "run_id",
        ),
    )

    run_id: Mapped[str] = mapped_column(String(128), primary_key=True)
    project_id: Mapped[str] = mapped_column(String(36), nullable=False)
    sample_revision_id: Mapped[str] = mapped_column(String(37), nullable=False)
    payload_digest: Mapped[str] = mapped_column(String(64), nullable=False)
    ordinal: Mapped[int] = mapped_column(Integer, primary_key=True)


class StoragePoolRow(Base):
    __tablename__ = "storage_pools"
    __table_args__ = (
        CheckConstraint(
            "length(storage_pool_id) = 37 AND substr(storage_pool_id, 1, 5) = 'stgp_'",
            name="ck_storage_pools_id",
        ),
        CheckConstraint(
            "length(trim(config_key)) BETWEEN 1 AND 255",
            name="ck_storage_pools_config_key",
        ),
        CheckConstraint(
            "length(trim(display_name)) BETWEEN 1 AND 255",
            name="ck_storage_pools_display_name",
        ),
        CheckConstraint(
            "archived_at IS NULL OR archived_at >= created_at",
            name="ck_storage_pools_archive_order",
        ),
        UniqueConstraint("config_key", name="uq_storage_pools_config_key"),
        Index(
            "ix_storage_pools_archived_created",
            "archived_at",
            "created_at",
            "storage_pool_id",
        ),
    )

    storage_pool_id: Mapped[str] = mapped_column(String(37), primary_key=True)
    config_key: Mapped[str] = mapped_column(String(255), nullable=False)
    display_name: Mapped[str] = mapped_column(String(255), nullable=False)
    created_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), nullable=False
    )
    archived_at: Mapped[datetime | None] = mapped_column(DateTime(timezone=True))


class ProjectStoragePoolBindingRow(Base):
    __tablename__ = "project_storage_pool_bindings"
    __table_args__ = (
        CheckConstraint(
            "project_id != 'prj_00000000000000000000000000000000'",
            name="ck_project_storage_pool_bindings_not_legacy",
        ),
        ForeignKeyConstraint(
            ["project_id"],
            ["projects.project_id"],
            name="fk_project_storage_pool_bindings_project",
            ondelete="RESTRICT",
        ),
        ForeignKeyConstraint(
            ["storage_pool_id"],
            ["storage_pools.storage_pool_id"],
            name="fk_project_storage_pool_bindings_pool",
            ondelete="RESTRICT",
        ),
        UniqueConstraint(
            "project_id",
            "storage_pool_id",
            name="uq_project_storage_pool_bindings_project_pool",
        ),
        Index(
            "ix_project_storage_pool_bindings_pool",
            "storage_pool_id",
            "project_id",
        ),
    )

    project_id: Mapped[str] = mapped_column(String(36), primary_key=True)
    storage_pool_id: Mapped[str] = mapped_column(String(37), nullable=False)
    bound_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), nullable=False)


class InputFileRow(Base):
    __tablename__ = "input_files"
    __table_args__ = (
        CheckConstraint(
            "length(input_file_id) = 37 AND substr(input_file_id, 1, 5) = 'inpf_'",
            name="ck_input_files_id",
        ),
        CheckConstraint(
            "length(trim(stable_key)) BETWEEN 1 AND 255",
            name="ck_input_files_stable_key",
        ),
        CheckConstraint(
            "archived_at IS NULL OR archived_at >= created_at",
            name="ck_input_files_archive_order",
        ),
        ForeignKeyConstraint(
            ["project_id", "storage_pool_id"],
            [
                "project_storage_pool_bindings.project_id",
                "project_storage_pool_bindings.storage_pool_id",
            ],
            name="fk_input_files_project_pool",
            ondelete="RESTRICT",
        ),
        UniqueConstraint(
            "project_id",
            "stable_key",
            name="uq_input_files_project_stable_key",
        ),
        UniqueConstraint(
            "project_id",
            "storage_pool_id",
            "input_file_id",
            name="uq_input_files_project_pool_file",
        ),
        Index(
            "ix_input_files_project_created",
            "project_id",
            "archived_at",
            "created_at",
            "input_file_id",
        ),
    )

    input_file_id: Mapped[str] = mapped_column(String(37), primary_key=True)
    project_id: Mapped[str] = mapped_column(String(36), nullable=False)
    storage_pool_id: Mapped[str] = mapped_column(String(37), nullable=False)
    stable_key: Mapped[str] = mapped_column(String(255), nullable=False)
    created_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), nullable=False
    )
    archived_at: Mapped[datetime | None] = mapped_column(DateTime(timezone=True))


class InputFileRevisionRow(Base):
    __tablename__ = "input_file_revisions"
    __table_args__ = (
        CheckConstraint(
            "length(input_file_revision_id) = 38 "
            "AND substr(input_file_revision_id, 1, 6) = 'inpfr_'",
            name="ck_input_file_revisions_id",
        ),
        CheckConstraint(
            "revision_number >= 1",
            name="ck_input_file_revisions_positive_revision",
        ),
        CheckConstraint(
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
        CheckConstraint(
            "size_bytes >= 0",
            name="ck_input_file_revisions_nonnegative_size",
        ),
        CheckConstraint(
            "length(content_sha256) = 64",
            name="ck_input_file_revisions_content_sha256_length",
        ),
        CheckConstraint(
            "digest_scheme = 'sha256-framed-input-file-revision-v1'",
            name="ck_input_file_revisions_digest_scheme",
        ),
        CheckConstraint(
            "length(digest) = 64",
            name="ck_input_file_revisions_digest_length",
        ),
        ForeignKeyConstraint(
            ["project_id", "storage_pool_id", "input_file_id"],
            [
                "input_files.project_id",
                "input_files.storage_pool_id",
                "input_files.input_file_id",
            ],
            name="fk_input_file_revisions_file",
            ondelete="RESTRICT",
        ),
        UniqueConstraint(
            "input_file_id",
            "revision_number",
            name="uq_input_file_revisions_file_number",
        ),
        UniqueConstraint(
            "project_id",
            "input_file_revision_id",
            "digest",
            name="uq_input_file_revisions_project_revision_digest",
        ),
        UniqueConstraint(
            "project_id",
            "input_file_id",
            "input_file_revision_id",
            "digest",
            "size_bytes",
            "content_sha256",
            name="uq_input_file_revisions_binding_evidence",
        ),
        Index(
            "ix_input_file_revisions_input_number",
            "input_file_id",
            "revision_number",
        ),
        Index(
            "ix_input_file_revisions_project_created",
            "project_id",
            "created_at",
            "input_file_revision_id",
        ),
    )

    input_file_revision_id: Mapped[str] = mapped_column(String(38), primary_key=True)
    input_file_id: Mapped[str] = mapped_column(String(37), nullable=False)
    project_id: Mapped[str] = mapped_column(String(36), nullable=False)
    storage_pool_id: Mapped[str] = mapped_column(String(37), nullable=False)
    revision_number: Mapped[int] = mapped_column(Integer, nullable=False)
    relative_path: Mapped[str] = mapped_column(Text, nullable=False)
    size_bytes: Mapped[int] = mapped_column(Integer, nullable=False)
    content_sha256: Mapped[str] = mapped_column(String(64), nullable=False)
    digest_scheme: Mapped[str] = mapped_column(String(64), nullable=False)
    digest: Mapped[str] = mapped_column(String(64), nullable=False)
    created_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), nullable=False
    )


class SnapshotInputBindingRow(Base):
    __tablename__ = "snapshot_input_bindings"
    __table_args__ = (
        CheckConstraint(
            "length(trim(workflow_id)) >= 1",
            name="ck_snapshot_input_bindings_workflow_id",
        ),
        CheckConstraint(
            "adapter_contract_version IS NULL "
            "OR length(trim(adapter_contract_version)) BETWEEN 1 AND 255",
            name="ck_snapshot_input_bindings_adapter_contract_version",
        ),
        CheckConstraint(
            "binding_mode IN ('compatibility_unresolved_v1', 'declared_input_uses_v1')",
            name="ck_snapshot_input_bindings_mode",
        ),
        CheckConstraint(
            "binding_mode = 'compatibility_unresolved_v1' "
            "OR adapter_contract_version IS NOT NULL",
            name="ck_snapshot_input_bindings_declared_adapter_version",
        ),
        CheckConstraint(
            "length(workflow_inputs_digest) = 64",
            name="ck_snapshot_input_bindings_workflow_digest_length",
        ),
        CheckConstraint(
            "length(project_sample_binding_digest) = 64",
            name="ck_snapshot_input_bindings_project_sample_digest_length",
        ),
        CheckConstraint(
            "binding_digest_scheme = 'sha256-framed-input-use-binding-envelope-v1'",
            name="ck_snapshot_input_bindings_digest_scheme",
        ),
        CheckConstraint(
            "length(binding_digest) = 64",
            name="ck_snapshot_input_bindings_digest_length",
        ),
        ForeignKeyConstraint(
            [
                "snapshot_id",
                "project_id",
                "project_sample_binding_digest",
            ],
            [
                "snapshot_project_bindings.snapshot_id",
                "snapshot_project_bindings.project_id",
                "snapshot_project_bindings.binding_digest",
            ],
            name="fk_snapshot_input_bindings_project_binding",
            ondelete="RESTRICT",
        ),
        UniqueConstraint(
            "snapshot_id",
            "project_id",
            name="uq_snapshot_input_bindings_project",
        ),
        Index(
            "ix_snapshot_input_bindings_project",
            "project_id",
            "created_at",
            "snapshot_id",
        ),
    )

    snapshot_id: Mapped[str] = mapped_column(String(64), primary_key=True)
    project_id: Mapped[str] = mapped_column(String(36), nullable=False)
    workflow_id: Mapped[str] = mapped_column(Text, nullable=False)
    adapter_contract_version: Mapped[str | None] = mapped_column(String(255))
    binding_mode: Mapped[str] = mapped_column(String(40), nullable=False)
    workflow_inputs_digest: Mapped[str] = mapped_column(String(64), nullable=False)
    project_sample_binding_digest: Mapped[str] = mapped_column(
        String(64), nullable=False
    )
    binding_digest_scheme: Mapped[str] = mapped_column(String(64), nullable=False)
    binding_digest: Mapped[str] = mapped_column(String(64), nullable=False)
    created_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), nullable=False
    )


class SnapshotInputUseRow(Base):
    __tablename__ = "snapshot_input_uses"
    __table_args__ = (
        CheckConstraint("ordinal >= 0", name="ck_snapshot_input_uses_ordinal"),
        CheckConstraint(
            "length(trim(input_use_key)) BETWEEN 1 AND 255",
            name="ck_snapshot_input_uses_key",
        ),
        CheckConstraint(
            "occurrence >= 0",
            name="ck_snapshot_input_uses_occurrence",
        ),
        CheckConstraint(
            "length(trim(capability_version)) BETWEEN 1 AND 255",
            name="ck_snapshot_input_uses_capability_version",
        ),
        CheckConstraint(
            "length(trim(closure_contract_version)) BETWEEN 1 AND 255",
            name="ck_snapshot_input_uses_closure_version",
        ),
        CheckConstraint(
            "(provenance_mode = 'transitional_unmanaged_v1' "
            "AND closure_digest_scheme IS NULL AND closure_digest IS NULL) OR "
            "(provenance_mode = 'managed_revision_v1' "
            "AND closure_digest_scheme = 'sha256-framed-input-closure-v1' "
            "AND length(closure_digest) = 64)",
            name="ck_snapshot_input_uses_provenance_evidence",
        ),
        ForeignKeyConstraint(
            ["snapshot_id", "project_id"],
            [
                "snapshot_input_bindings.snapshot_id",
                "snapshot_input_bindings.project_id",
            ],
            name="fk_snapshot_input_uses_binding",
            ondelete="RESTRICT",
        ),
        UniqueConstraint(
            "snapshot_id",
            "input_use_key",
            "occurrence",
            name="uq_snapshot_input_uses_identity",
        ),
        UniqueConstraint(
            "snapshot_id",
            "ordinal",
            "project_id",
            name="uq_snapshot_input_uses_project",
        ),
        Index(
            "ix_snapshot_input_uses_key",
            "input_use_key",
            "capability_version",
            "snapshot_id",
        ),
    )

    snapshot_id: Mapped[str] = mapped_column(String(64), primary_key=True)
    project_id: Mapped[str] = mapped_column(String(36), nullable=False)
    ordinal: Mapped[int] = mapped_column(Integer, primary_key=True)
    input_use_key: Mapped[str] = mapped_column(String(255), nullable=False)
    occurrence: Mapped[int] = mapped_column(Integer, nullable=False)
    capability_version: Mapped[str] = mapped_column(String(255), nullable=False)
    closure_contract_version: Mapped[str] = mapped_column(String(255), nullable=False)
    provenance_mode: Mapped[str] = mapped_column(String(40), nullable=False)
    closure_digest_scheme: Mapped[str | None] = mapped_column(String(64))
    closure_digest: Mapped[str | None] = mapped_column(String(64))


class SnapshotInputMemberRow(Base):
    __tablename__ = "snapshot_input_members"
    __table_args__ = (
        CheckConstraint(
            "use_ordinal >= 0 AND member_ordinal >= 0",
            name="ck_snapshot_input_members_ordinals",
        ),
        CheckConstraint(
            "length(trim(logical_member_key)) BETWEEN 1 AND 255",
            name="ck_snapshot_input_members_member_key",
        ),
        CheckConstraint(
            "length(revision_digest) = 64",
            name="ck_snapshot_input_members_revision_digest_length",
        ),
        CheckConstraint(
            "size_bytes >= 0",
            name="ck_snapshot_input_members_nonnegative_size",
        ),
        CheckConstraint(
            "length(content_sha256) = 64",
            name="ck_snapshot_input_members_content_sha256_length",
        ),
        ForeignKeyConstraint(
            ["snapshot_id", "use_ordinal", "project_id"],
            [
                "snapshot_input_uses.snapshot_id",
                "snapshot_input_uses.ordinal",
                "snapshot_input_uses.project_id",
            ],
            name="fk_snapshot_input_members_use",
            ondelete="RESTRICT",
        ),
        ForeignKeyConstraint(
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
            name="fk_snapshot_input_members_revision",
            ondelete="RESTRICT",
        ),
        UniqueConstraint(
            "snapshot_id",
            "use_ordinal",
            "logical_member_key",
            name="uq_snapshot_input_members_member_key",
        ),
        UniqueConstraint(
            "snapshot_id",
            "use_ordinal",
            "input_file_revision_id",
            name="uq_snapshot_input_members_revision",
        ),
        Index(
            "ix_snapshot_input_members_revision",
            "input_file_revision_id",
            "snapshot_id",
        ),
    )

    snapshot_id: Mapped[str] = mapped_column(String(64), primary_key=True)
    project_id: Mapped[str] = mapped_column(String(36), nullable=False)
    use_ordinal: Mapped[int] = mapped_column(Integer, primary_key=True)
    member_ordinal: Mapped[int] = mapped_column(Integer, primary_key=True)
    logical_member_key: Mapped[str] = mapped_column(String(255), nullable=False)
    input_file_id: Mapped[str] = mapped_column(String(37), nullable=False)
    input_file_revision_id: Mapped[str] = mapped_column(String(38), nullable=False)
    revision_digest: Mapped[str] = mapped_column(String(64), nullable=False)
    size_bytes: Mapped[int] = mapped_column(Integer, nullable=False)
    content_sha256: Mapped[str] = mapped_column(String(64), nullable=False)


class RunInputBindingRow(Base):
    __tablename__ = "run_input_bindings"
    __table_args__ = (
        CheckConstraint(
            "length(trim(workflow_id)) >= 1",
            name="ck_run_input_bindings_workflow_id",
        ),
        CheckConstraint(
            "adapter_contract_version IS NULL "
            "OR length(trim(adapter_contract_version)) BETWEEN 1 AND 255",
            name="ck_run_input_bindings_adapter_contract_version",
        ),
        CheckConstraint(
            "binding_mode IN ('compatibility_unresolved_v1', 'declared_input_uses_v1')",
            name="ck_run_input_bindings_mode",
        ),
        CheckConstraint(
            "binding_mode = 'compatibility_unresolved_v1' "
            "OR adapter_contract_version IS NOT NULL",
            name="ck_run_input_bindings_declared_adapter_version",
        ),
        CheckConstraint(
            "length(workflow_inputs_digest) = 64",
            name="ck_run_input_bindings_workflow_digest_length",
        ),
        CheckConstraint(
            "length(project_sample_binding_digest) = 64",
            name="ck_run_input_bindings_project_sample_digest_length",
        ),
        CheckConstraint(
            "binding_digest_scheme = 'sha256-framed-input-use-binding-envelope-v1'",
            name="ck_run_input_bindings_digest_scheme",
        ),
        CheckConstraint(
            "length(binding_digest) = 64",
            name="ck_run_input_bindings_digest_length",
        ),
        ForeignKeyConstraint(
            [
                "run_id",
                "project_id",
                "project_sample_binding_digest",
            ],
            [
                "run_project_bindings.run_id",
                "run_project_bindings.project_id",
                "run_project_bindings.binding_digest",
            ],
            name="fk_run_input_bindings_project_binding",
            ondelete="CASCADE",
        ),
        UniqueConstraint(
            "run_id",
            "project_id",
            name="uq_run_input_bindings_project",
        ),
        Index(
            "ix_run_input_bindings_project",
            "project_id",
            "created_at",
            "run_id",
        ),
    )

    run_id: Mapped[str] = mapped_column(String(128), primary_key=True)
    project_id: Mapped[str] = mapped_column(String(36), nullable=False)
    workflow_id: Mapped[str] = mapped_column(Text, nullable=False)
    adapter_contract_version: Mapped[str | None] = mapped_column(String(255))
    binding_mode: Mapped[str] = mapped_column(String(40), nullable=False)
    workflow_inputs_digest: Mapped[str] = mapped_column(String(64), nullable=False)
    project_sample_binding_digest: Mapped[str] = mapped_column(
        String(64), nullable=False
    )
    binding_digest_scheme: Mapped[str] = mapped_column(String(64), nullable=False)
    binding_digest: Mapped[str] = mapped_column(String(64), nullable=False)
    created_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), nullable=False
    )


class RunInputUseRow(Base):
    __tablename__ = "run_input_uses"
    __table_args__ = (
        CheckConstraint("ordinal >= 0", name="ck_run_input_uses_ordinal"),
        CheckConstraint(
            "length(trim(input_use_key)) BETWEEN 1 AND 255",
            name="ck_run_input_uses_key",
        ),
        CheckConstraint(
            "occurrence >= 0",
            name="ck_run_input_uses_occurrence",
        ),
        CheckConstraint(
            "length(trim(capability_version)) BETWEEN 1 AND 255",
            name="ck_run_input_uses_capability_version",
        ),
        CheckConstraint(
            "length(trim(closure_contract_version)) BETWEEN 1 AND 255",
            name="ck_run_input_uses_closure_version",
        ),
        CheckConstraint(
            "(provenance_mode = 'transitional_unmanaged_v1' "
            "AND closure_digest_scheme IS NULL AND closure_digest IS NULL) OR "
            "(provenance_mode = 'managed_revision_v1' "
            "AND closure_digest_scheme = 'sha256-framed-input-closure-v1' "
            "AND length(closure_digest) = 64)",
            name="ck_run_input_uses_provenance_evidence",
        ),
        ForeignKeyConstraint(
            ["run_id", "project_id"],
            [
                "run_input_bindings.run_id",
                "run_input_bindings.project_id",
            ],
            name="fk_run_input_uses_binding",
            ondelete="CASCADE",
        ),
        UniqueConstraint(
            "run_id",
            "input_use_key",
            "occurrence",
            name="uq_run_input_uses_identity",
        ),
        UniqueConstraint(
            "run_id",
            "ordinal",
            "project_id",
            name="uq_run_input_uses_project",
        ),
        Index(
            "ix_run_input_uses_key",
            "input_use_key",
            "capability_version",
            "run_id",
        ),
    )

    run_id: Mapped[str] = mapped_column(String(128), primary_key=True)
    project_id: Mapped[str] = mapped_column(String(36), nullable=False)
    ordinal: Mapped[int] = mapped_column(Integer, primary_key=True)
    input_use_key: Mapped[str] = mapped_column(String(255), nullable=False)
    occurrence: Mapped[int] = mapped_column(Integer, nullable=False)
    capability_version: Mapped[str] = mapped_column(String(255), nullable=False)
    closure_contract_version: Mapped[str] = mapped_column(String(255), nullable=False)
    provenance_mode: Mapped[str] = mapped_column(String(40), nullable=False)
    closure_digest_scheme: Mapped[str | None] = mapped_column(String(64))
    closure_digest: Mapped[str | None] = mapped_column(String(64))


class RunInputMemberRow(Base):
    __tablename__ = "run_input_members"
    __table_args__ = (
        CheckConstraint(
            "use_ordinal >= 0 AND member_ordinal >= 0",
            name="ck_run_input_members_ordinals",
        ),
        CheckConstraint(
            "length(trim(logical_member_key)) BETWEEN 1 AND 255",
            name="ck_run_input_members_member_key",
        ),
        CheckConstraint(
            "length(revision_digest) = 64",
            name="ck_run_input_members_revision_digest_length",
        ),
        CheckConstraint(
            "size_bytes >= 0",
            name="ck_run_input_members_nonnegative_size",
        ),
        CheckConstraint(
            "length(content_sha256) = 64",
            name="ck_run_input_members_content_sha256_length",
        ),
        ForeignKeyConstraint(
            ["run_id", "use_ordinal", "project_id"],
            [
                "run_input_uses.run_id",
                "run_input_uses.ordinal",
                "run_input_uses.project_id",
            ],
            name="fk_run_input_members_use",
            ondelete="CASCADE",
        ),
        ForeignKeyConstraint(
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
            name="fk_run_input_members_revision",
            ondelete="RESTRICT",
        ),
        UniqueConstraint(
            "run_id",
            "use_ordinal",
            "logical_member_key",
            name="uq_run_input_members_member_key",
        ),
        UniqueConstraint(
            "run_id",
            "use_ordinal",
            "input_file_revision_id",
            name="uq_run_input_members_revision",
        ),
        Index(
            "ix_run_input_members_revision",
            "input_file_revision_id",
            "run_id",
        ),
    )

    run_id: Mapped[str] = mapped_column(String(128), primary_key=True)
    project_id: Mapped[str] = mapped_column(String(36), nullable=False)
    use_ordinal: Mapped[int] = mapped_column(Integer, primary_key=True)
    member_ordinal: Mapped[int] = mapped_column(Integer, primary_key=True)
    logical_member_key: Mapped[str] = mapped_column(String(255), nullable=False)
    input_file_id: Mapped[str] = mapped_column(String(37), nullable=False)
    input_file_revision_id: Mapped[str] = mapped_column(String(38), nullable=False)
    revision_digest: Mapped[str] = mapped_column(String(64), nullable=False)
    size_bytes: Mapped[int] = mapped_column(Integer, nullable=False)
    content_sha256: Mapped[str] = mapped_column(String(64), nullable=False)


class RunRow(Base):
    __tablename__ = "runs"
    __table_args__ = (
        Index("ix_runs_created_run_id", "created_at", "run_id"),
        Index(
            "ix_runs_workflow_created_run_id",
            "workflow_id",
            "created_at",
            "run_id",
        ),
        Index(
            "ix_runs_status_created_run_id",
            "status",
            "created_at",
            "run_id",
        ),
        Index(
            "ix_runs_workflow_status_created_run_id",
            "workflow_id",
            "status",
            "created_at",
            "run_id",
        ),
    )

    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    run_id: Mapped[str] = mapped_column(String(128), nullable=False, unique=True)
    workflow_id: Mapped[str] = mapped_column(String(255), nullable=False, index=True)
    inputs: Mapped[dict[str, Any]] = mapped_column(JSON, nullable=False)
    status: Mapped[str] = mapped_column(String(32), nullable=False, index=True)
    created_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), nullable=False
    )
    updated_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), nullable=False
    )
    started_at: Mapped[datetime | None] = mapped_column(DateTime(timezone=True))
    ended_at: Mapped[datetime | None] = mapped_column(DateTime(timezone=True))
    current_stage: Mapped[str | None] = mapped_column(String(255))
    cancellation_reason: Mapped[str | None] = mapped_column(Text)
    error: Mapped[dict[str, Any] | None] = mapped_column(JSON)
    tags: Mapped[dict[str, str]] = mapped_column(JSON, nullable=False)


class ValidatedInputSnapshotRow(Base):
    __tablename__ = "validated_input_snapshots"
    __table_args__ = (
        CheckConstraint(
            "(consumed_run_id IS NULL AND consumed_at IS NULL) OR "
            "(consumed_run_id IS NOT NULL AND consumed_at IS NOT NULL)",
            name="ck_validated_input_snapshots_consumption_pair",
        ),
        CheckConstraint(
            "expires_at > validated_at",
            name="ck_validated_input_snapshots_expiry",
        ),
        CheckConstraint(
            "validation_outcome = 'adapter_validation_succeeded'",
            name="ck_validated_input_snapshots_success",
        ),
        CheckConstraint(
            "length(payload_digest) = 64",
            name="ck_validated_input_snapshots_digest_length",
        ),
        Index("ix_validated_input_snapshots_workflow_id", "workflow_id"),
        Index("ix_validated_input_snapshots_expires_at", "expires_at"),
    )

    snapshot_id: Mapped[str] = mapped_column(String(64), primary_key=True)
    workflow_id: Mapped[str] = mapped_column(String(255), nullable=False)
    adapter_version: Mapped[str] = mapped_column(String(128), nullable=False)
    schema_version: Mapped[str] = mapped_column(String(64), nullable=False)
    schema_dialect: Mapped[str] = mapped_column(String(255), nullable=False)
    canonical_payload: Mapped[str] = mapped_column(Text, nullable=False)
    payload_digest_scheme: Mapped[str] = mapped_column(String(64), nullable=False)
    payload_digest: Mapped[str] = mapped_column(String(64), nullable=False)
    validation_outcome: Mapped[str] = mapped_column(String(64), nullable=False)
    validation_issue_codes: Mapped[list[str]] = mapped_column(JSON, nullable=False)
    validated_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), nullable=False
    )
    expires_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), nullable=False
    )
    build_adapter_version: Mapped[str] = mapped_column(String(128), nullable=False)
    build_scheme: Mapped[str] = mapped_column(String(64), nullable=False)
    build_logical_entrypoint: Mapped[str] = mapped_column(String(512), nullable=False)
    build_digest: Mapped[str] = mapped_column(String(64), nullable=False)
    build_captured_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), nullable=False
    )
    consumed_run_id: Mapped[str | None] = mapped_column(
        String(128),
        ForeignKey("runs.run_id", ondelete="RESTRICT"),
        unique=True,
    )
    consumed_at: Mapped[datetime | None] = mapped_column(DateTime(timezone=True))


class RunExecutionAssignmentRow(Base):
    __tablename__ = "run_execution_assignments"
    __table_args__ = (
        UniqueConstraint(
            "job_id",
            name="uq_run_execution_assignments_job_id",
        ),
        CheckConstraint(
            "claimed_at IS NULL OR dispatched_at IS NOT NULL",
            name="ck_run_execution_assignments_claim_requires_dispatch",
        ),
        CheckConstraint(
            "(cancellation_requested_at IS NULL AND cancellation_reason IS NULL) "
            "OR (cancellation_requested_at IS NOT NULL AND "
            "cancellation_reason IS NOT NULL)",
            name="ck_run_execution_assignments_request_reason_pair",
        ),
        CheckConstraint(
            "cancellation_requested_at IS NULL OR claimed_at IS NOT NULL",
            name="ck_run_execution_assignments_request_requires_claim",
        ),
        CheckConstraint(
            "cancellation_acknowledged_at IS NULL "
            "OR cancellation_requested_at IS NOT NULL",
            name="ck_run_execution_assignments_ack_requires_request",
        ),
        CheckConstraint(
            "requeue_requested_at IS NULL OR dispatched_at IS NOT NULL",
            name="ck_run_execution_assignments_requeue_requires_dispatch",
        ),
        CheckConstraint(
            "requeue_confirmed_at IS NULL OR requeue_requested_at IS NOT NULL",
            name="ck_run_execution_assignments_requeue_confirm_requires_request",
        ),
        CheckConstraint(
            "requeue_confirmed_at IS NULL "
            "OR requeue_confirmed_at >= requeue_requested_at",
            name="ck_run_execution_assignments_requeue_confirmation_order",
        ),
    )

    run_id: Mapped[str] = mapped_column(
        String(128),
        ForeignKey("runs.run_id", ondelete="CASCADE"),
        primary_key=True,
    )
    job_id: Mapped[str] = mapped_column(String(255), nullable=False)
    backend: Mapped[str] = mapped_column(String(64), nullable=False)
    queue_name: Mapped[str] = mapped_column(String(128), nullable=False)
    created_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), nullable=False
    )
    dispatched_at: Mapped[datetime | None] = mapped_column(DateTime(timezone=True))
    claimed_at: Mapped[datetime | None] = mapped_column(DateTime(timezone=True))
    cancellation_requested_at: Mapped[datetime | None] = mapped_column(
        DateTime(timezone=True)
    )
    cancellation_reason: Mapped[str | None] = mapped_column(Text)
    cancellation_acknowledged_at: Mapped[datetime | None] = mapped_column(
        DateTime(timezone=True)
    )
    requeue_requested_at: Mapped[datetime | None] = mapped_column(
        DateTime(timezone=True)
    )
    requeue_confirmed_at: Mapped[datetime | None] = mapped_column(
        DateTime(timezone=True)
    )


class RunWorkflowBuildIdentityRow(Base):
    __tablename__ = "run_workflow_build_identities"

    run_id: Mapped[str] = mapped_column(
        String(128),
        ForeignKey("runs.run_id", ondelete="CASCADE"),
        primary_key=True,
    )
    workflow_id: Mapped[str] = mapped_column(String(255), nullable=False)
    adapter_version: Mapped[str] = mapped_column(String(255), nullable=False)
    scheme: Mapped[str] = mapped_column(String(64), nullable=False)
    logical_entrypoint: Mapped[str] = mapped_column(String(512), nullable=False)
    digest: Mapped[str] = mapped_column(String(64), nullable=False)
    captured_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), nullable=False
    )


class RunEventRow(Base):
    __tablename__ = "run_events"
    __table_args__ = (
        UniqueConstraint("run_id", "sequence", name="uq_run_events_run_sequence"),
        UniqueConstraint("run_id", "event_id", name="uq_run_events_run_event_id"),
        Index("ix_run_events_run_sequence", "run_id", "sequence"),
    )

    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    event_id: Mapped[str] = mapped_column(String(128), nullable=False)
    run_id: Mapped[str] = mapped_column(
        String(128),
        ForeignKey("runs.run_id", ondelete="CASCADE"),
        nullable=False,
    )
    sequence: Mapped[int] = mapped_column(Integer, nullable=False)
    event_type: Mapped[str] = mapped_column(String(128), nullable=False)
    timestamp: Mapped[datetime] = mapped_column(DateTime(timezone=True), nullable=False)
    status: Mapped[str | None] = mapped_column(String(32))
    stage: Mapped[str | None] = mapped_column(String(255))
    message: Mapped[str] = mapped_column(Text, nullable=False)
    context: Mapped[dict[str, Any]] = mapped_column(JSON, nullable=False)
    issue: Mapped[dict[str, Any] | None] = mapped_column(JSON)


class RunLogRow(Base):
    __tablename__ = "run_logs"
    __table_args__ = (
        UniqueConstraint(
            "run_id",
            "stream_name",
            "sequence",
            name="uq_run_logs_run_stream_sequence",
        ),
        UniqueConstraint(
            "run_id",
            "stream_name",
            "chunk_id",
            name="uq_run_logs_run_stream_chunk_id",
        ),
        Index("ix_run_logs_run_stream_sequence", "run_id", "stream_name", "sequence"),
    )

    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    chunk_id: Mapped[str] = mapped_column(String(128), nullable=False)
    run_id: Mapped[str] = mapped_column(
        String(128),
        ForeignKey("runs.run_id", ondelete="CASCADE"),
        nullable=False,
    )
    stream_name: Mapped[str] = mapped_column(String(32), nullable=False)
    sequence: Mapped[int] = mapped_column(Integer, nullable=False)
    timestamp: Mapped[datetime] = mapped_column(DateTime(timezone=True), nullable=False)
    lines: Mapped[list[str]] = mapped_column(JSON, nullable=False)


class RunArtifactRow(Base):
    __tablename__ = "run_artifacts"
    __table_args__ = (
        UniqueConstraint("run_id", "artifact_id", name="uq_run_artifacts_run_artifact"),
        Index("ix_run_artifacts_run_id", "run_id"),
    )

    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    artifact_id: Mapped[str] = mapped_column(String(128), nullable=False)
    run_id: Mapped[str] = mapped_column(
        String(128),
        ForeignKey("runs.run_id", ondelete="CASCADE"),
        nullable=False,
    )
    artifact_type: Mapped[str] = mapped_column(String(128), nullable=False)
    name: Mapped[str] = mapped_column(String(255), nullable=False)
    uri: Mapped[str] = mapped_column(Text, nullable=False)
    mime_type: Mapped[str | None] = mapped_column(String(255))
    produced_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), nullable=False
    )
    revision: Mapped[str | None] = mapped_column(String(76))
    artifact_metadata: Mapped[dict[str, Any]] = mapped_column(JSON, nullable=False)


class ArtifactPublicationRow(Base):
    __tablename__ = "artifact_publications"
    __table_args__ = (
        CheckConstraint(
            "length(artifact_id) BETWEEN 1 AND 128 AND "
            "substr(artifact_id, 1, 1) GLOB '[A-Za-z]' AND "
            "artifact_id NOT GLOB '*[^A-Za-z0-9_.-]*'",
            name="ck_artifact_publications_artifact_id",
        ),
        CheckConstraint(
            "length(artifact_generation) = 76 AND "
            "substr(artifact_generation, 1, 12) = 'artifactgen-' AND "
            "substr(artifact_generation, 13) NOT GLOB '*[^0-9a-f]*'",
            name="ck_artifact_publications_generation",
        ),
        CheckConstraint(
            "length(artifact_revision) = 76 AND "
            "substr(artifact_revision, 1, 12) = 'artifactrev-' AND "
            "substr(artifact_revision, 13) NOT GLOB '*[^0-9a-f]*'",
            name="ck_artifact_publications_revision",
        ),
        CheckConstraint(
            "length(output_type) BETWEEN 1 AND 128 AND "
            "substr(output_type, 1, 1) GLOB '[A-Za-z]' AND "
            "output_type NOT GLOB '*[^A-Za-z0-9_.-]*'",
            name="ck_artifact_publications_output_type",
        ),
        ForeignKeyConstraint(
            ["run_id"],
            ["runs.run_id"],
            name="fk_artifact_publications_run",
            ondelete="RESTRICT",
        ),
        Index(
            "ix_artifact_publications_published",
            "published_at",
            "run_id",
            "artifact_generation",
            "artifact_id",
        ),
        Index(
            "ix_artifact_publications_run_published",
            "run_id",
            "published_at",
            "artifact_generation",
            "artifact_id",
        ),
        Index(
            "ix_artifact_publications_output_type_published",
            "output_type",
            "published_at",
            "run_id",
            "artifact_generation",
            "artifact_id",
        ),
    )

    run_id: Mapped[str] = mapped_column(String(128), primary_key=True)
    artifact_id: Mapped[str] = mapped_column(String(128), primary_key=True)
    artifact_generation: Mapped[str] = mapped_column(String(76), primary_key=True)
    artifact_revision: Mapped[str] = mapped_column(String(76), nullable=False)
    output_type: Mapped[str] = mapped_column(String(128), nullable=False)
    published_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), nullable=False
    )


class RunResultStateRow(Base):
    __tablename__ = "run_result_states"
    __table_args__ = (
        CheckConstraint(
            "artifact_revision >= 0 AND qc_revision >= 0",
            name="ck_run_result_states_nonnegative_revisions",
        ),
        CheckConstraint(
            "(artifact_revision = 0 AND artifact_generation IS NULL AND "
            "artifact_manifest_digest IS NULL) OR "
            "(artifact_revision > 0 AND artifact_generation IS NOT NULL AND "
            "artifact_manifest_digest IS NOT NULL)",
            name="ck_run_result_states_artifact_binding",
        ),
        CheckConstraint(
            "(qc_revision = 0 AND qc_generation IS NULL AND "
            "qc_manifest_digest IS NULL AND qc_artifact_generation IS NULL) OR "
            "(qc_revision > 0 AND ((qc_generation IS NULL AND "
            "qc_manifest_digest IS NULL AND qc_artifact_generation IS NULL) OR "
            "(qc_generation IS NOT NULL AND qc_manifest_digest IS NOT NULL AND "
            "qc_artifact_generation IS NOT NULL)))",
            name="ck_run_result_states_qc_binding",
        ),
        CheckConstraint(
            "(artifact_attempt_id IS NULL AND artifact_attempt_status IS NULL) OR "
            "(artifact_attempt_id IS NOT NULL AND artifact_attempt_status IN "
            "('pending', 'succeeded', 'failed'))",
            name="ck_run_result_states_artifact_attempt",
        ),
        CheckConstraint(
            "(qc_attempt_id IS NULL AND qc_attempt_status IS NULL AND "
            "qc_attempt_artifact_generation IS NULL) OR "
            "(qc_attempt_id IS NOT NULL AND qc_attempt_status IN "
            "('pending', 'succeeded', 'failed') AND "
            "qc_attempt_artifact_generation IS NOT NULL)",
            name="ck_run_result_states_qc_attempt",
        ),
        CheckConstraint(
            "artifact_outcome IS NULL OR artifact_outcome IN ('succeeded', 'failed')",
            name="ck_run_result_states_artifact_outcome",
        ),
        CheckConstraint(
            "qc_outcome IS NULL OR qc_outcome IN "
            "('succeeded', 'failed', 'invalidated')",
            name="ck_run_result_states_qc_outcome",
        ),
        CheckConstraint(
            "(artifact_outcome = 'failed' AND artifact_reason_code IS NOT NULL) OR "
            "(artifact_outcome IS NULL AND artifact_reason_code IS NULL) OR "
            "(artifact_outcome = 'succeeded' AND artifact_reason_code IS NULL)",
            name="ck_run_result_states_artifact_reason",
        ),
        CheckConstraint(
            "(qc_outcome = 'failed' AND qc_reason_code IS NOT NULL) OR "
            "(qc_outcome IS NULL AND qc_reason_code IS NULL) OR "
            "(qc_outcome IN ('succeeded', 'invalidated') AND qc_reason_code IS NULL)",
            name="ck_run_result_states_qc_reason",
        ),
    )

    run_id: Mapped[str] = mapped_column(
        String(128),
        ForeignKey("runs.run_id", ondelete="CASCADE"),
        primary_key=True,
    )
    artifact_revision: Mapped[int] = mapped_column(Integer, nullable=False, default=0)
    artifact_generation: Mapped[str | None] = mapped_column(String(76))
    artifact_manifest_digest: Mapped[str | None] = mapped_column(String(64))
    artifact_attempt_id: Mapped[str | None] = mapped_column(String(78))
    artifact_attempt_status: Mapped[str | None] = mapped_column(String(16))
    artifact_outcome: Mapped[str | None] = mapped_column(String(16))
    artifact_reason_code: Mapped[str | None] = mapped_column(String(128))
    qc_revision: Mapped[int] = mapped_column(Integer, nullable=False, default=0)
    qc_generation: Mapped[str | None] = mapped_column(String(70))
    qc_manifest_digest: Mapped[str | None] = mapped_column(String(64))
    qc_attempt_id: Mapped[str | None] = mapped_column(String(78))
    qc_attempt_status: Mapped[str | None] = mapped_column(String(16))
    qc_attempt_artifact_generation: Mapped[str | None] = mapped_column(String(76))
    qc_artifact_generation: Mapped[str | None] = mapped_column(String(76))
    qc_outcome: Mapped[str | None] = mapped_column(String(16))
    qc_reason_code: Mapped[str | None] = mapped_column(String(128))


class RunResultAttemptRow(Base):
    __tablename__ = "run_result_attempts"
    __table_args__ = (
        CheckConstraint(
            "(result_kind = 'artifact' AND artifact_generation IS NULL) OR "
            "(result_kind = 'qc' AND artifact_generation IS NOT NULL)",
            name="ck_run_result_attempts_binding",
        ),
        Index("ix_run_result_attempts_run_id", "run_id"),
    )

    attempt_id: Mapped[str] = mapped_column(String(78), primary_key=True)
    run_id: Mapped[str] = mapped_column(
        String(128),
        ForeignKey("runs.run_id", ondelete="CASCADE"),
        nullable=False,
    )
    result_kind: Mapped[str] = mapped_column(String(16), nullable=False)
    artifact_generation: Mapped[str | None] = mapped_column(String(76))


class RunQcMetricRow(Base):
    __tablename__ = "run_qc_metrics"
    __table_args__ = (
        UniqueConstraint(
            "run_id",
            "metric_id",
            name="uq_run_qc_metrics_run_metric",
        ),
        CheckConstraint(
            "scope IN ('run', 'sample', 'experiment')",
            name="ck_run_qc_metrics_scope",
        ),
        CheckConstraint(
            "qc_flag IS NULL OR qc_flag IN ('pass', 'warning', 'fail')",
            name="ck_run_qc_metrics_flag",
        ),
        CheckConstraint(
            "(scope = 'run' AND sample_id IS NULL AND experiment_id IS NULL) OR "
            "(scope = 'sample' AND sample_id IS NOT NULL) OR "
            "(scope = 'experiment' AND sample_id IS NULL AND "
            "experiment_id IS NOT NULL)",
            name="ck_run_qc_metrics_scope_identifiers",
        ),
        CheckConstraint(
            "length(value_text) BETWEEN 1 AND 40",
            name="ck_run_qc_metrics_value_text_length",
        ),
        Index("ix_run_qc_metrics_run_id", "run_id"),
    )

    id: Mapped[int] = mapped_column(Integer, primary_key=True, autoincrement=True)
    metric_id: Mapped[str] = mapped_column(String(128), nullable=False)
    run_id: Mapped[str] = mapped_column(
        String(128),
        ForeignKey("runs.run_id", ondelete="CASCADE"),
        nullable=False,
    )
    metric_key: Mapped[str] = mapped_column(String(128), nullable=False)
    display_name: Mapped[str] = mapped_column(String(255), nullable=False)
    value_text: Mapped[str] = mapped_column(String(40), nullable=False)
    unit: Mapped[str] = mapped_column(String(32), nullable=False)
    scope: Mapped[str] = mapped_column(String(16), nullable=False)
    sample_id: Mapped[str | None] = mapped_column(String(255))
    experiment_id: Mapped[str | None] = mapped_column(String(255))
    assay: Mapped[str | None] = mapped_column(String(255))
    qc_flag: Mapped[str | None] = mapped_column(String(16))
    source_artifact_id: Mapped[str] = mapped_column(String(128), nullable=False)
    produced_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), nullable=False
    )


class ReferenceProfileRow(Base):
    """Stable logical profile; private configuration lives on revisions."""

    __tablename__ = "reference_profiles"
    __table_args__ = (
        CheckConstraint(
            "length(profile_id) = 37 AND substr(profile_id, 1, 5) = 'refp_'",
            name="ck_reference_profiles_id",
        ),
        CheckConstraint(
            "length(trim(safe_key)) BETWEEN 1 AND 255",
            name="ck_reference_profiles_safe_key",
        ),
        ForeignKeyConstraint(
            ["profile_id", "enabled_revision_id"],
            [
                "reference_profile_revisions.profile_id",
                "reference_profile_revisions.revision_id",
            ],
            name="fk_reference_profiles_enabled_revision",
            ondelete="RESTRICT",
            use_alter=True,
        ),
        UniqueConstraint("safe_key", name="uq_reference_profiles_safe_key"),
        UniqueConstraint(
            "profile_id",
            "enabled_revision_id",
            name="uq_reference_profiles_enabled_revision",
        ),
        Index(
            "ix_reference_profiles_enabled_created",
            "enabled_revision_id",
            "created_at",
            "profile_id",
        ),
    )

    profile_id: Mapped[str] = mapped_column(String(37), primary_key=True)
    safe_key: Mapped[str] = mapped_column(String(255), nullable=False)
    created_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), nullable=False
    )
    enabled_revision_id: Mapped[str | None] = mapped_column(String(38))


class ReferenceProfileRevisionRow(Base):
    """Append-only profile metadata plus aggregate public identity."""

    __tablename__ = "reference_profile_revisions"
    __table_args__ = (
        CheckConstraint(
            "length(revision_id) = 38 AND substr(revision_id, 1, 6) = 'refpr_'",
            name="ck_reference_profile_revisions_id",
        ),
        CheckConstraint(
            "revision_number >= 1",
            name="ck_reference_profile_revisions_positive_number",
        ),
        CheckConstraint(
            "length(trim(display_name)) BETWEEN 1 AND 255",
            name="ck_reference_profile_revisions_display_name",
        ),
        CheckConstraint(
            "length(trim(organism)) BETWEEN 1 AND 255",
            name="ck_reference_profile_revisions_organism",
        ),
        CheckConstraint(
            "length(trim(assembly)) BETWEEN 1 AND 255",
            name="ck_reference_profile_revisions_assembly",
        ),
        CheckConstraint(
            "length(trim(config_key)) BETWEEN 1 AND 255",
            name="ck_reference_profile_revisions_config_key",
        ),
        CheckConstraint(
            "public_identity_scheme = 'sha256-framed-reference-profile-revision-v1'",
            name="ck_reference_profile_revisions_identity_scheme",
        ),
        CheckConstraint(
            "length(public_identity_sha256) = 64",
            name="ck_reference_profile_revisions_identity_length",
        ),
        ForeignKeyConstraint(
            ["profile_id"],
            ["reference_profiles.profile_id"],
            name="fk_reference_profile_revisions_profile",
            ondelete="RESTRICT",
        ),
        UniqueConstraint(
            "profile_id",
            "revision_id",
            name="uq_reference_profile_revisions_profile_revision",
        ),
        UniqueConstraint(
            "profile_id",
            "revision_number",
            name="uq_reference_profile_revisions_profile_number",
        ),
        UniqueConstraint(
            "profile_id",
            "revision_id",
            "public_identity_sha256",
            name="uq_reference_profile_revisions_identity",
        ),
        Index(
            "ix_reference_profile_revisions_profile_number",
            "profile_id",
            "revision_number",
        ),
    )

    revision_id: Mapped[str] = mapped_column(String(38), primary_key=True)
    profile_id: Mapped[str] = mapped_column(String(37), nullable=False)
    revision_number: Mapped[int] = mapped_column(Integer, nullable=False)
    display_name: Mapped[str] = mapped_column(String(255), nullable=False)
    organism: Mapped[str] = mapped_column(String(255), nullable=False)
    assembly: Mapped[str] = mapped_column(String(255), nullable=False)
    config_key: Mapped[str] = mapped_column(String(255), nullable=False)
    public_identity_scheme: Mapped[str] = mapped_column(String(64), nullable=False)
    public_identity_sha256: Mapped[str] = mapped_column(String(64), nullable=False)
    created_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), nullable=False
    )


class ReferenceProfileWorkflowBindingRow(Base):
    """Adapter-owned identity for one revision/workflow compatibility edge."""

    __tablename__ = "reference_profile_workflow_bindings"
    __table_args__ = (
        CheckConstraint(
            "length(trim(workflow_id)) BETWEEN 1 AND 255",
            name="ck_reference_profile_workflow_bindings_workflow",
        ),
        CheckConstraint(
            "length(trim(contract_version)) BETWEEN 1 AND 255",
            name="ck_reference_profile_workflow_bindings_contract",
        ),
        CheckConstraint(
            "identity_scheme = 'sha256-framed-adapter-reference-binding-v1'",
            name="ck_reference_profile_workflow_bindings_identity_scheme",
        ),
        CheckConstraint(
            "length(identity_sha256) = 64",
            name="ck_reference_profile_workflow_bindings_identity_length",
        ),
        ForeignKeyConstraint(
            ["profile_id", "revision_id"],
            [
                "reference_profile_revisions.profile_id",
                "reference_profile_revisions.revision_id",
            ],
            name="fk_reference_profile_workflow_bindings_revision",
            ondelete="RESTRICT",
        ),
        UniqueConstraint(
            "profile_id",
            "revision_id",
            "workflow_id",
            "contract_version",
            "identity_sha256",
            name="uq_reference_profile_workflow_bindings_identity",
        ),
        Index(
            "ix_reference_profile_workflow_bindings_workflow",
            "workflow_id",
            "profile_id",
            "revision_id",
        ),
    )

    profile_id: Mapped[str] = mapped_column(String(37), primary_key=True)
    revision_id: Mapped[str] = mapped_column(String(38), primary_key=True)
    workflow_id: Mapped[str] = mapped_column(String(255), primary_key=True)
    contract_version: Mapped[str] = mapped_column(String(255), nullable=False)
    identity_scheme: Mapped[str] = mapped_column(String(64), nullable=False)
    identity_sha256: Mapped[str] = mapped_column(String(64), nullable=False)


class SnapshotReferenceBindingRow(Base):
    """Exact path-free Reference Profile evidence frozen into a snapshot."""

    __tablename__ = "snapshot_reference_bindings"
    __table_args__ = (
        CheckConstraint(
            "revision_public_identity_scheme = "
            "'sha256-framed-reference-profile-revision-v1'",
            name="ck_snapshot_reference_bindings_revision_scheme",
        ),
        CheckConstraint(
            "length(revision_public_identity_sha256) = 64",
            name="ck_snapshot_reference_bindings_revision_identity",
        ),
        CheckConstraint(
            "adapter_identity_scheme = 'sha256-framed-adapter-reference-binding-v1'",
            name="ck_snapshot_reference_bindings_adapter_scheme",
        ),
        CheckConstraint(
            "length(adapter_identity_sha256) = 64",
            name="ck_snapshot_reference_bindings_adapter_identity",
        ),
        CheckConstraint(
            "binding_digest_scheme = 'sha256-framed-reference-profile-binding-v1'",
            name="ck_snapshot_reference_bindings_digest_scheme",
        ),
        CheckConstraint(
            "length(binding_digest) = 64",
            name="ck_snapshot_reference_bindings_digest_length",
        ),
        ForeignKeyConstraint(
            ["snapshot_id"],
            ["validated_input_snapshots.snapshot_id"],
            name="fk_snapshot_reference_bindings_snapshot",
            ondelete="RESTRICT",
        ),
        ForeignKeyConstraint(
            ["profile_id", "revision_id", "revision_public_identity_sha256"],
            [
                "reference_profile_revisions.profile_id",
                "reference_profile_revisions.revision_id",
                "reference_profile_revisions.public_identity_sha256",
            ],
            name="fk_snapshot_reference_bindings_revision",
            ondelete="RESTRICT",
        ),
        ForeignKeyConstraint(
            [
                "profile_id",
                "revision_id",
                "workflow_id",
                "adapter_contract_version",
                "adapter_identity_sha256",
            ],
            [
                "reference_profile_workflow_bindings.profile_id",
                "reference_profile_workflow_bindings.revision_id",
                "reference_profile_workflow_bindings.workflow_id",
                "reference_profile_workflow_bindings.contract_version",
                "reference_profile_workflow_bindings.identity_sha256",
            ],
            name="fk_snapshot_reference_bindings_workflow_binding",
            ondelete="RESTRICT",
        ),
        Index(
            "ix_snapshot_reference_bindings_revision",
            "revision_id",
            "snapshot_id",
        ),
    )

    snapshot_id: Mapped[str] = mapped_column(String(64), primary_key=True)
    profile_id: Mapped[str] = mapped_column(String(37), nullable=False)
    revision_id: Mapped[str] = mapped_column(String(38), nullable=False)
    workflow_id: Mapped[str] = mapped_column(String(255), nullable=False)
    revision_public_identity_scheme: Mapped[str] = mapped_column(
        String(64), nullable=False
    )
    revision_public_identity_sha256: Mapped[str] = mapped_column(
        String(64), nullable=False
    )
    adapter_contract_version: Mapped[str] = mapped_column(String(255), nullable=False)
    adapter_identity_scheme: Mapped[str] = mapped_column(String(64), nullable=False)
    adapter_identity_sha256: Mapped[str] = mapped_column(String(64), nullable=False)
    binding_digest_scheme: Mapped[str] = mapped_column(String(64), nullable=False)
    binding_digest: Mapped[str] = mapped_column(String(64), nullable=False)
    bound_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), nullable=False)


class RunReferenceBindingRow(Base):
    """Exact snapshot evidence copied atomically to one durable run."""

    __tablename__ = "run_reference_bindings"
    __table_args__ = (
        CheckConstraint(
            "revision_public_identity_scheme = "
            "'sha256-framed-reference-profile-revision-v1'",
            name="ck_run_reference_bindings_revision_scheme",
        ),
        CheckConstraint(
            "length(revision_public_identity_sha256) = 64",
            name="ck_run_reference_bindings_revision_identity",
        ),
        CheckConstraint(
            "adapter_identity_scheme = 'sha256-framed-adapter-reference-binding-v1'",
            name="ck_run_reference_bindings_adapter_scheme",
        ),
        CheckConstraint(
            "length(adapter_identity_sha256) = 64",
            name="ck_run_reference_bindings_adapter_identity",
        ),
        CheckConstraint(
            "binding_digest_scheme = 'sha256-framed-reference-profile-binding-v1'",
            name="ck_run_reference_bindings_digest_scheme",
        ),
        CheckConstraint(
            "length(binding_digest) = 64",
            name="ck_run_reference_bindings_digest_length",
        ),
        ForeignKeyConstraint(
            ["run_id"],
            ["runs.run_id"],
            name="fk_run_reference_bindings_run",
            ondelete="CASCADE",
        ),
        ForeignKeyConstraint(
            ["profile_id", "revision_id", "revision_public_identity_sha256"],
            [
                "reference_profile_revisions.profile_id",
                "reference_profile_revisions.revision_id",
                "reference_profile_revisions.public_identity_sha256",
            ],
            name="fk_run_reference_bindings_revision",
            ondelete="RESTRICT",
        ),
        ForeignKeyConstraint(
            [
                "profile_id",
                "revision_id",
                "workflow_id",
                "adapter_contract_version",
                "adapter_identity_sha256",
            ],
            [
                "reference_profile_workflow_bindings.profile_id",
                "reference_profile_workflow_bindings.revision_id",
                "reference_profile_workflow_bindings.workflow_id",
                "reference_profile_workflow_bindings.contract_version",
                "reference_profile_workflow_bindings.identity_sha256",
            ],
            name="fk_run_reference_bindings_workflow_binding",
            ondelete="RESTRICT",
        ),
        Index("ix_run_reference_bindings_revision", "revision_id", "run_id"),
    )

    run_id: Mapped[str] = mapped_column(String(128), primary_key=True)
    profile_id: Mapped[str] = mapped_column(String(37), nullable=False)
    revision_id: Mapped[str] = mapped_column(String(38), nullable=False)
    workflow_id: Mapped[str] = mapped_column(String(255), nullable=False)
    revision_public_identity_scheme: Mapped[str] = mapped_column(
        String(64), nullable=False
    )
    revision_public_identity_sha256: Mapped[str] = mapped_column(
        String(64), nullable=False
    )
    adapter_contract_version: Mapped[str] = mapped_column(String(255), nullable=False)
    adapter_identity_scheme: Mapped[str] = mapped_column(String(64), nullable=False)
    adapter_identity_sha256: Mapped[str] = mapped_column(String(64), nullable=False)
    binding_digest_scheme: Mapped[str] = mapped_column(String(64), nullable=False)
    binding_digest: Mapped[str] = mapped_column(String(64), nullable=False)
    bound_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), nullable=False)
