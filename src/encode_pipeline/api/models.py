"""Pydantic v2 request/response models for the workflow platform API."""

from __future__ import annotations

from datetime import datetime
from typing import Any, Literal

from pydantic import BaseModel, Field


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

    type: Literal["config_edit", "schema_hint", "assay_selection", "general"] = "general"
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
