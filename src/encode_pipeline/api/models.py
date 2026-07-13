"""Pydantic v2 request/response models for the workflow platform API."""

from __future__ import annotations

from collections.abc import Mapping
from datetime import datetime, timezone
from pathlib import PurePosixPath
import re
from typing import Any, Literal
from urllib.parse import quote

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

from encode_pipeline.platform.runs import (
    build_qc_metric_id,
    validate_qc_identifier_token,
)
from encode_pipeline.services.run_repositories import canonical_decimal_text


_LOGICAL_ID_PATTERN = r"^[A-Za-z][A-Za-z0-9_.-]{0,127}$"
_MIME_TYPE_PATTERN = re.compile(
    r"^[A-Za-z0-9][A-Za-z0-9!#$&^_.+-]{0,126}/"
    r"[A-Za-z0-9][A-Za-z0-9!#$&^_.+-]{0,126}$"
)
_WINDOWS_ABSOLUTE_PATTERN = re.compile(r"^[A-Za-z]:[\\/]")
_ENVIRONMENT_ASSIGNMENT_PATTERN = re.compile(r"(?:^|\s)[A-Za-z_][A-Za-z0-9_]*=")
_EMBEDDED_POSIX_PATH_PATTERN = re.compile(r"(?<![A-Za-z0-9._-])/[^/\s]")
_ARTIFACT_DOWNLOAD_REASON_CODE_PATTERN = re.compile(r"^[A-Z][A-Z0-9_]{0,127}$")


class IssueResponse(BaseModel):
    """JSON-ready Issue shape from the PR84 contract."""

    code: str
    message: str
    severity: str = "error"
    path: str | None = None
    source: str | None = None
    technical_message: str | None = None
    hint: str | None = None
    context: dict[str, Any] = Field(default_factory=dict)


class WorkflowMetadataResponse(BaseModel):
    """WorkflowMetadata JSON shape."""

    workflow_id: str
    name: str
    version: str
    description: str | None = None
    engines: list[str] = Field(default_factory=list)
    tags: list[str] = Field(default_factory=list)


class WorkflowCapabilityResponse(BaseModel):
    """WorkflowCapabilities JSON shape."""

    supports: list[str] = Field(default_factory=list)


class WorkflowListItem(BaseModel):
    """One entry in the workflow list response."""

    metadata: WorkflowMetadataResponse
    capabilities: WorkflowCapabilityResponse


class WorkflowListResponse(BaseModel):
    """Envelope for GET /api/v1/workflows."""

    ok: bool
    workflows: list[WorkflowListItem]
    issues: list[IssueResponse] = Field(default_factory=list)


class SchemaResponse(BaseModel):
    """Envelope for GET /api/v1/workflows/{workflow_id}/schema."""

    ok: bool
    workflow_id: str
    schema_hints: dict[str, Any] = Field(default_factory=dict)
    issues: list[IssueResponse] = Field(default_factory=list)


class ValidationRequest(BaseModel):
    """Request body for POST /api/v1/workflows/{workflow_id}/validate."""

    config: dict[str, Any]
    samples: str | list[dict[str, str]] | None = None
    options: dict[str, Any] = Field(default_factory=dict)


class ValidationResponse(BaseModel):
    """Envelope for POST /api/v1/workflows/{workflow_id}/validate."""

    ok: bool
    workflow_id: str | None
    value: Any | None
    issues: list[IssueResponse] = Field(default_factory=list)


class AgentToolCall(BaseModel):
    """Transparent record of a read-only tool invocation by the agent."""

    tool_name: str
    input_summary: dict[str, Any] = Field(default_factory=dict)
    output_summary: str = ""
    read_only: Literal[True] = True


class AgentSuggestion(BaseModel):
    """A model-generated suggestion returned by the agent."""

    type: Literal["config_edit", "schema_hint", "assay_selection", "general"] = (
        "general"
    )
    description: str = Field(..., min_length=1)
    target_path: str | None = None
    current_value: Any | None = None
    proposed_value: Any | None = None
    rationale: str | None = None
    disclaimer: str = Field(
        default="This is a model-generated suggestion. Please verify against your experimental design and ENCODE guidelines.",
        min_length=1,
    )


class AgentContext(BaseModel):
    """Snapshot of platform state sent to the agent with a user message."""

    current_issues: list[IssueResponse] = Field(default_factory=list)
    current_config: dict[str, Any] = Field(default_factory=dict)
    current_schema: dict[str, Any] = Field(default_factory=dict)


class AgentRequest(BaseModel):
    """Request body for POST /api/v1/workflows/{workflow_id}/agent/chat."""

    session_id: str | None = None
    message: str = Field(..., min_length=1)
    context: AgentContext | None = None


class AgentResponse(BaseModel):
    """Envelope for POST /api/v1/workflows/{workflow_id}/agent/chat."""

    ok: bool
    session_id: str | None = None
    message: str
    suggestions: list[AgentSuggestion] = Field(default_factory=list)
    tool_calls: list[AgentToolCall] = Field(default_factory=list)
    issues: list[IssueResponse] = Field(default_factory=list)


class RunRecordResponse(BaseModel):
    """JSON-ready RunRecord shape."""

    run_id: str
    workflow_id: str
    inputs: dict[str, Any]
    status: str
    created_at: datetime
    updated_at: datetime
    started_at: datetime | None
    ended_at: datetime | None
    current_stage: str | None
    cancellation_reason: str | None
    error: IssueResponse | None = None
    tags: dict[str, str] = Field(default_factory=dict)


class RunEventResponse(BaseModel):
    """JSON-ready RunEvent shape."""

    event_id: str
    run_id: str
    sequence: int
    event_type: str
    timestamp: datetime
    status: str | None = None
    stage: str | None = None
    message: str
    context: dict[str, Any] = Field(default_factory=dict)
    issue: IssueResponse | None = None


class RunLogChunkResponse(BaseModel):
    """JSON-ready RunLogChunk shape."""

    chunk_id: str
    run_id: str
    stream_name: str
    sequence: int
    timestamp: datetime
    lines: list[str]


class ArtifactMetadataResponse(BaseModel):
    """Whitelisted adapter metadata safe for the public read API."""

    model_config = ConfigDict(extra="ignore")

    catalog_id: str | None = Field(default=None, max_length=512)
    scope: str | None = Field(default=None, max_length=512)
    level: str | None = Field(default=None, max_length=512)
    sample_id: str | None = Field(default=None, max_length=512)
    experiment_id: str | None = Field(default=None, max_length=512)
    assay: str | None = Field(default=None, max_length=512)
    target: str | None = Field(default=None, max_length=512)
    genome: str | None = Field(default=None, max_length=512)
    method: str | None = Field(default=None, max_length=512)
    qc_flag: str | None = Field(default=None, max_length=512)

    @field_validator(
        "catalog_id",
        "scope",
        "level",
        "sample_id",
        "experiment_id",
        "assay",
        "target",
        "genome",
        "method",
        "qc_flag",
    )
    @classmethod
    def validate_public_metadata_text(cls, value: str | None) -> str | None:
        if value is not None and _is_unsafe_public_text(value):
            raise ValueError("artifact metadata value is not public-safe")
        return value


class _PersistedArtifactMetadata(ArtifactMetadataResponse):
    """Validated PR127 persistence shape before public projection."""

    relative_path: str = Field(min_length=1, max_length=2048)
    output_type: str = Field(pattern=_LOGICAL_ID_PATTERN, max_length=128)
    size_bytes: int = Field(ge=0, strict=True)

    @field_validator("relative_path")
    @classmethod
    def validate_relative_path(cls, value: str) -> str:
        path = PurePosixPath(value)
        if (
            "\x00" in value
            or "\\" in value
            or value.startswith(("/", "~"))
            or _WINDOWS_ABSOLUTE_PATTERN.match(value) is not None
            or path.as_posix() != value
            or len(path.parts) < 2
            or path.parts[0] != "results"
            or any(part in {"", ".", ".."} for part in value.split("/"))
        ):
            raise ValueError("artifact relative path is not public-safe")
        return value


class ArtifactReferenceResponse(BaseModel):
    """Disclosure-safe API projection of a persisted RunArtifactRef."""

    artifact_id: str = Field(pattern=_LOGICAL_ID_PATTERN, max_length=128)
    run_id: str = Field(min_length=1, max_length=128)
    artifact_type: Literal["file"]
    name: str = Field(min_length=1, max_length=255)
    uri: str = Field(min_length=1, max_length=2048)
    mime_type: str | None = Field(default=None, max_length=255)
    produced_at: datetime
    relative_path: str
    output_type: str
    size_bytes: int
    metadata: ArtifactMetadataResponse

    @classmethod
    def from_artifact(
        cls,
        artifact: Any,
        *,
        expected_run_id: str,
    ) -> "ArtifactReferenceResponse":
        """Build a strict public projection without domain ``to_dict`` output."""
        if artifact.run_id != expected_run_id:
            raise ValueError("artifact does not belong to the requested run")
        persisted = _PersistedArtifactMetadata.model_validate(artifact.metadata)
        controlled = ArtifactMetadataResponse.model_validate(artifact.metadata)
        return cls(
            artifact_id=artifact.artifact_id,
            run_id=artifact.run_id,
            artifact_type=artifact.artifact_type,
            name=artifact.name,
            uri=artifact.uri,
            mime_type=artifact.mime_type,
            produced_at=artifact.produced_at,
            relative_path=persisted.relative_path,
            output_type=persisted.output_type,
            size_bytes=persisted.size_bytes,
            metadata=controlled,
        )

    @field_validator("name")
    @classmethod
    def validate_name(cls, value: str) -> str:
        if "/" in value or "\\" in value or _has_control_character(value):
            raise ValueError("artifact name is not public-safe")
        return value

    @field_validator("mime_type")
    @classmethod
    def validate_mime_type(cls, value: str | None) -> str | None:
        if value is not None and _MIME_TYPE_PATTERN.fullmatch(value) is None:
            raise ValueError("artifact MIME type is not public-safe")
        return value

    @model_validator(mode="after")
    def validate_derived_fields(self) -> "ArtifactReferenceResponse":
        expected_uri = (
            f"run://runs/{quote(self.run_id, safe='')}/artifacts/{self.artifact_id}"
        )
        if self.uri != expected_uri:
            raise ValueError("artifact URI is not the expected opaque reference")
        if self.name != PurePosixPath(self.relative_path).name:
            raise ValueError("artifact name does not match its relative path")
        return self


class RunArtifactsResponse(BaseModel):
    """Envelope for GET /api/v1/runs/{run_id}/artifacts."""

    ok: bool
    run_id: str
    artifacts: list[ArtifactReferenceResponse] = Field(default_factory=list)
    next_cursor: str | None = None
    issues: list[IssueResponse] = Field(default_factory=list)


class RunArtifactDetailResponse(BaseModel):
    """Envelope for one run-scoped artifact reference."""

    ok: bool
    run_id: str
    artifact: ArtifactReferenceResponse | None = None
    issues: list[IssueResponse] = Field(default_factory=list)


class ArtifactDownloadIssueContextResponse(BaseModel):
    """Strict public context projection for artifact download failures."""

    model_config = ConfigDict(extra="forbid")

    reason_code: str | None = Field(
        default=None,
        pattern=_ARTIFACT_DOWNLOAD_REASON_CODE_PATTERN.pattern,
        max_length=128,
    )

    @classmethod
    def from_issue_context(
        cls,
        context: Mapping[str, object],
    ) -> ArtifactDownloadIssueContextResponse:
        """Retain only a valid stable reason code from an internal Issue."""
        reason_code = context.get("reason_code")
        if (
            not isinstance(reason_code, str)
            or _ARTIFACT_DOWNLOAD_REASON_CODE_PATTERN.fullmatch(reason_code) is None
        ):
            reason_code = None
        return cls(reason_code=reason_code)


class ArtifactDownloadIssueResponse(BaseModel):
    """Allowlisted public Issue fields for artifact download failures."""

    code: str
    message: str
    severity: str = "error"
    path: str | None = None
    source: str | None = None
    hint: str | None = None
    context: ArtifactDownloadIssueContextResponse = Field(
        default_factory=ArtifactDownloadIssueContextResponse
    )


class RunArtifactDownloadErrorResponse(BaseModel):
    """Disclosure-safe envelope for an artifact download failure."""

    ok: Literal[False] = False
    run_id: str
    artifact_id: str
    issues: list[ArtifactDownloadIssueResponse] = Field(default_factory=list)


class QcMetricResponse(BaseModel):
    """Lossless disclosure-safe projection of one persisted QC metric."""

    metric_id: str = Field(pattern=r"^qcmetric-[0-9a-f]{64}$", max_length=73)
    metric_key: str = Field(
        pattern=r"^[a-z][a-z0-9_]*(?:\.[a-z][a-z0-9_]*)*$",
        max_length=128,
    )
    display_name: str = Field(min_length=1, max_length=255)
    value: str = Field(
        pattern=r"^-?(?:0|[1-9]\d{0,25})(?:\.\d{1,12})?$",
        max_length=40,
    )
    unit: Literal["count", "fraction", "ratio"]
    scope: Literal["run", "sample", "experiment"]
    sample_id: str | None = Field(max_length=255)
    experiment_id: str | None = Field(max_length=255)
    assay: str | None = Field(max_length=255)
    qc_flag: Literal["pass", "warning", "fail"] | None
    source_artifact_id: str = Field(pattern=_LOGICAL_ID_PATTERN, max_length=128)
    produced_at: datetime

    @classmethod
    def from_metric(
        cls,
        metric: Any,
        *,
        expected_run_id: str,
    ) -> "QcMetricResponse":
        """Build a strict public projection from a run-scoped domain value."""
        if metric.run_id != expected_run_id:
            raise ValueError("QC metric does not belong to the requested run")
        projected = cls(
            metric_id=metric.metric_id,
            metric_key=metric.metric_key,
            display_name=metric.display_name,
            value=canonical_decimal_text(metric.value),
            unit=metric.unit,
            scope=metric.scope,
            sample_id=metric.sample_id,
            experiment_id=metric.experiment_id,
            assay=metric.assay,
            qc_flag=metric.qc_flag,
            source_artifact_id=metric.source_artifact_id,
            produced_at=metric.produced_at,
        )
        expected_metric_id = build_qc_metric_id(
            projected.metric_key,
            projected.scope,
            projected.sample_id,
            projected.experiment_id,
        )
        if projected.metric_id != expected_metric_id:
            raise ValueError("QC metric ID does not match its semantic coordinates")
        return projected

    @field_validator("display_name")
    @classmethod
    def validate_display_name(cls, value: str) -> str:
        if _is_unsafe_public_text(value):
            raise ValueError("QC metric display name is not public-safe")
        return value

    @field_validator("sample_id", "experiment_id", "assay")
    @classmethod
    def validate_identifier(cls, value: str | None) -> str | None:
        if value is not None:
            validate_qc_identifier_token(value)
        return value

    @field_validator("produced_at")
    @classmethod
    def normalize_produced_at(cls, value: datetime) -> datetime:
        if value.tzinfo is None or value.utcoffset() is None:
            raise ValueError("QC metric produced_at must be timezone-aware")
        return value.astimezone(timezone.utc)

    @model_validator(mode="after")
    def validate_scope_identifiers(self) -> "QcMetricResponse":
        if self.scope == "run" and (
            self.sample_id is not None or self.experiment_id is not None
        ):
            raise ValueError("run QC scope cannot include sample or experiment IDs")
        if self.scope == "sample" and self.sample_id is None:
            raise ValueError("sample QC scope requires sample_id")
        if self.scope == "experiment" and (
            self.sample_id is not None or self.experiment_id is None
        ):
            raise ValueError("experiment QC scope requires only experiment_id")
        return self


class RunQcMetricsResponse(BaseModel):
    """Envelope for GET /api/v1/runs/{run_id}/qc-metrics."""

    ok: bool
    run_id: str
    qc_metrics: list[QcMetricResponse] = Field(default_factory=list)
    next_cursor: str | None = None
    issues: list[IssueResponse] = Field(default_factory=list)


def _has_control_character(value: str) -> bool:
    return any(ord(character) < 32 or ord(character) == 127 for character in value)


def _is_unsafe_public_text(value: str) -> bool:
    lowered = value.lower()
    return (
        _has_control_character(value)
        or value.startswith(("/", "~", "\\"))
        or "\\" in value
        or _WINDOWS_ABSOLUTE_PATTERN.match(value) is not None
        or lowered.startswith("file:")
        or "://" in lowered
        or _ENVIRONMENT_ASSIGNMENT_PATTERN.search(value) is not None
        or _EMBEDDED_POSIX_PATH_PATTERN.search(value) is not None
    )


class RunCreateRequest(BaseModel):
    """Request body for POST /api/v1/workflows/{workflow_id}/runs."""

    config: dict[str, Any]
    samples: str | list[dict[str, str]] | None = None
    options: dict[str, Any] = Field(default_factory=dict)
    tags: dict[str, str] = Field(default_factory=dict)


class RunResponse(BaseModel):
    """Envelope for single-run endpoints."""

    ok: bool
    run: RunRecordResponse | None = None
    issues: list[IssueResponse] = Field(default_factory=list)


class RunEventsResponse(BaseModel):
    """Envelope for GET /api/v1/runs/{run_id}/events."""

    ok: bool
    run_id: str
    events: list[RunEventResponse] = Field(default_factory=list)
    next_cursor: str | None = None
    issues: list[IssueResponse] = Field(default_factory=list)


class RunLogsResponse(BaseModel):
    """Envelope for GET /api/v1/runs/{run_id}/logs."""

    ok: bool
    run_id: str
    stream_name: str
    chunks: list[RunLogChunkResponse] = Field(default_factory=list)
    next_cursor: str | None = None
    issues: list[IssueResponse] = Field(default_factory=list)
