"""SQLAlchemy rows kept behind the persistence repository boundary."""

from __future__ import annotations

from datetime import datetime
from typing import Any

from sqlalchemy import (
    CheckConstraint,
    DateTime,
    ForeignKey,
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
    artifact_metadata: Mapped[dict[str, Any]] = mapped_column(JSON, nullable=False)


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
