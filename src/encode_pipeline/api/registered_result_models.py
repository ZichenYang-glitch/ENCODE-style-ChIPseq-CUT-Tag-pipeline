"""Public HTTP models for registered results. Private coordinates stay out."""

from __future__ import annotations

from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field

from encode_pipeline.api.models import IssueResponse


class RegisteredResultByteResponse(BaseModel):
    """Size and digest the server computed from the registered file."""

    model_config = ConfigDict(extra="forbid")

    size_bytes: int
    content_sha256: str


class RegisteredResultRevisionResponse(BaseModel):
    """One exact revision. Storage coordinates are not part of this object."""

    model_config = ConfigDict(extra="forbid")

    project_id: str
    result_id: str
    result_revision_id: str
    revision_number: int
    revision_status: Literal["current", "superseded"]
    display_name: str
    format: str
    schema_version: Literal[1]
    origin: dict[str, Any]
    method: dict[str, Any]
    reference: dict[str, Any]
    bytes: RegisteredResultByteResponse
    created_at: str


class RegisteredResultListResponse(BaseModel):
    """One page of current revisions, or one result's history."""

    model_config = ConfigDict(extra="forbid")

    ok: bool
    results: list[RegisteredResultRevisionResponse]
    next_cursor: str | None
    issues: list[IssueResponse] = Field(default_factory=list)


class RegisteredResultDetailResponse(BaseModel):
    """One revision read, registration, or rename."""

    model_config = ConfigDict(extra="forbid")

    ok: bool
    revision: RegisteredResultRevisionResponse
    issues: list[IssueResponse] = Field(default_factory=list)


class RegisteredResultArchiveResponse(BaseModel):
    """Archive acknowledgement. Archiving does not create a revision."""

    model_config = ConfigDict(extra="forbid")

    ok: bool
    result_id: str
    archived_at: str
    issues: list[IssueResponse] = Field(default_factory=list)


class RegisterRegisteredResultRequest(BaseModel):
    """Metadata plus the pool-relative path of a file that is already present."""

    model_config = ConfigDict(extra="forbid")

    stable_key: str
    display_name: str
    format: str
    origin: dict[str, Any]
    method: dict[str, Any]
    reference: dict[str, Any]
    pool_relative_path: str
    schema_version: Literal[1] | None = None


class ReviseRegisteredResultRequest(BaseModel):
    """A new revision. Omitted description fields stay as they are on the head."""

    model_config = ConfigDict(extra="forbid")

    expected_revision_number: int
    display_name: str | None = None
    format: str | None = None
    origin: dict[str, Any] | None = None
    method: dict[str, Any] | None = None
    reference: dict[str, Any] | None = None
    pool_relative_path: str | None = None


class ArchiveRegisteredResultRequest(BaseModel):
    """Retire one result at its current head without appending a revision."""

    model_config = ConfigDict(extra="forbid")

    expected_revision_number: int
