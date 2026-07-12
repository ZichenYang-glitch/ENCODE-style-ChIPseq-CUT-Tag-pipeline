"""Pydantic v2 request/response models for the workflow platform API."""

from __future__ import annotations

from datetime import datetime
from pathlib import PurePosixPath
import re
from typing import Any, Literal
from urllib.parse import quote

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator


_LOGICAL_ID_PATTERN = r"^[A-Za-z][A-Za-z0-9_.-]{0,127}$"
_MIME_TYPE_PATTERN = re.compile(
    r"^[A-Za-z0-9][A-Za-z0-9!#$&^_.+-]{0,126}/"
    r"[A-Za-z0-9][A-Za-z0-9!#$&^_.+-]{0,126}$"
)
_WINDOWS_ABSOLUTE_PATTERN = re.compile(r"^[A-Za-z]:[\\/]")
_PUBLIC_METADATA_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.:+-]{0,511}$")


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
        if value is not None and (
            _is_unsafe_public_text(value)
            or _PUBLIC_METADATA_PATTERN.fullmatch(value) is None
        ):
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


def _has_control_character(value: str) -> bool:
    return any(ord(character) < 32 or ord(character) == 127 for character in value)


def _is_unsafe_public_text(value: str) -> bool:
    lowered = value.lower()
    return (
        _has_control_character(value)
        or value.startswith(("/", "~", "\\\\"))
        or _WINDOWS_ABSOLUTE_PATTERN.match(value) is not None
        or lowered.startswith("file:")
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
