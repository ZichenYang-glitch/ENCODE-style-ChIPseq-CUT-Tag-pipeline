"""Pydantic v2 request/response models for the workflow platform API."""

from __future__ import annotations

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
    current_schema: dict[str, Any] | None = None


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
