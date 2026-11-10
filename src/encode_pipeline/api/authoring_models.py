"""Request and response models for restricted authoring sessions."""

from __future__ import annotations

from datetime import datetime
from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

from encode_pipeline.platform.authoring import normalize_draft_content
from encode_pipeline.api.models import (
    InputFileRevisionSelectionRequest,
    IssueResponse,
    JsonValue,
    ProjectId,
    ReferenceProfileRevisionId,
    SampleId,
    SampleRequestPayload,
    SampleRevisionId,
    ValidatedInputSnapshotResponse,
)


WorkflowId = Field(pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,254}$", max_length=255)


AuthoringSid = Field(pattern=r"^aas_[0-9a-f]{32}$", max_length=36)
ProposalId = Field(pattern=r"^prop_[0-9a-f]{32}$", max_length=37)
TurnId = Field(pattern=r"^atr_[0-9a-f]{32}$", max_length=36)


class AuthoringErrorResponse(BaseModel):
    """Stable authoring failure. Issue text is fixed and path-free."""

    model_config = ConfigDict(extra="forbid")

    ok: Literal[False] = False
    reason_code: str
    issues: list[IssueResponse]


class AuthoringFieldConstraintResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    type: str
    minimum: int | None = None
    maximum: int | None = None


class AuthoringFieldResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    token: str
    label: str
    type: str
    constraints: AuthoringFieldConstraintResponse


class AuthoringIssueLabelResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    code: str
    label: str


class AuthoringCapabilityResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    available: bool
    reason_code: str | None = None
    projection_version: str
    projection_digest: str | None = None
    allowed_intents: list[Literal["explain_issues", "suggest"]]
    allowed_issue_codes: list[AuthoringIssueLabelResponse]
    fields: list[AuthoringFieldResponse]


class AuthoringConnectionResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    state: Literal["closed", "blocked", "ready"]
    reason_code: str | None = None


class AuthoringDraftDocument(BaseModel):
    model_config = ConfigDict(extra="forbid")

    config: dict[str, JsonValue]
    samples: SampleRequestPayload = None
    options: dict[str, JsonValue]
    sample_revision_ids: list[SampleRevisionId]
    input_selections: list[InputFileRevisionSelectionRequest]
    reference_profile_revision_id: ReferenceProfileRevisionId | None = None


class AuthoringSessionResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    ok: Literal[True] = True
    sid: str = AuthoringSid
    project_id: ProjectId
    workflow_id: str
    draft_state: Literal["active", "expired"]
    draft_revision: int
    expires_at: datetime
    purge_after: datetime
    schema_digest: str
    projection_digest: str | None = None
    draft: AuthoringDraftDocument | None
    authoring_capability: AuthoringCapabilityResponse
    connection: AuthoringConnectionResponse
    requires_new_connection: bool
    changed: bool | None = None


class AuthoringSessionSummary(BaseModel):
    model_config = ConfigDict(extra="forbid")

    sid: str = AuthoringSid
    project_id: ProjectId
    workflow_id: str
    draft_state: Literal["active", "expired"]
    draft_revision: int
    expires_at: datetime
    purge_after: datetime
    connection: AuthoringConnectionResponse
    requires_new_connection: bool


class AuthoringSessionListResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    ok: Literal[True] = True
    sessions: list[AuthoringSessionSummary]
    next_cursor: str | None = None


class AuthoringSessionCreateRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    project_id: ProjectId
    workflow_id: str = WorkflowId


class AuthoringDraftWriteRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    expected_revision: int = Field(ge=0, le=1_000_000)
    config: dict[str, JsonValue]
    samples: SampleRequestPayload = None
    options: dict[str, JsonValue] = Field(default_factory=dict)
    sample_revision_ids: list[SampleRevisionId] = Field(default_factory=list)
    input_selections: list[InputFileRevisionSelectionRequest] = Field(
        default_factory=list,
        max_length=256,
    )
    reference_profile_revision_id: ReferenceProfileRevisionId | None = None

    @model_validator(mode="after")
    def draft_is_canonical_json(self) -> "AuthoringDraftWriteRequest":
        normalize_draft_content(
            {
                "config": self.config,
                "input_selections": [
                    selection.model_dump() for selection in self.input_selections
                ],
                "options": self.options,
                "reference_profile_revision_id": self.reference_profile_revision_id,
                "sample_revision_ids": list(self.sample_revision_ids),
                "samples": self.samples,
            }
        )
        return self


class AuthoringPreferenceRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    field_token: str = Field(min_length=1, max_length=64)
    value: int

    @field_validator("value")
    @classmethod
    def value_is_integer(cls, value: int) -> int:
        if isinstance(value, bool):
            raise ValueError("preference value must be an integer")
        return value


class AuthoringQuestionRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    intent: Literal["suggest", "explain_issues"]
    field_tokens: list[str] = Field(max_length=16)
    preferences: list[AuthoringPreferenceRequest] = Field(
        default_factory=list,
        max_length=16,
    )


class AuthoringTurnCreateRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    client_turn_id: str = Field(pattern=r"^[A-Za-z0-9_-]{1,64}$", max_length=64)
    expected_revision: int = Field(ge=0, le=1_000_000)
    question: AuthoringQuestionRequest


class AuthoringProposalChangeResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    field_token: str
    value: int
    before: dict[str, Any]
    after: dict[str, Any]


class AuthoringProposalResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    changes: list[AuthoringProposalChangeResponse]
    reason_codes: list[str]


class AuthoringExplanationResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    reason_codes: list[str] = Field(max_length=16)


class AuthoringReceiptResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    proposal_id: str = ProposalId
    base_revision: int
    result_revision: int
    changed: bool
    adopted_at: datetime
    draft_expires_at: datetime


class AuthoringTurnResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    ok: Literal[True] = True
    turn_id: str = TurnId
    client_turn_id: str
    status: Literal["generating", "completed", "failed"]
    failure_reason: str | None = None
    result_kind: Literal["proposal", "explanation"] | None = None
    proposal_id: str | None = None
    proposal: AuthoringProposalResponse | None = None
    proposal_digest: str | None = Field(
        default=None,
        pattern=r"^[0-9a-f]{64}$",
        max_length=64,
    )
    explanation: AuthoringExplanationResponse | None = None
    proposal_state: Literal["pending", "adopted", "rejected", "invalidated"] | None = (
        None
    )
    adoptable: bool
    reason_code: str | None = None
    proposal_expires_at: datetime | None = None
    adoption_receipt: AuthoringReceiptResponse | None = None


class AuthoringAdoptRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    expected_revision: int = Field(ge=0, le=1_000_000)
    proposal_digest: str = Field(pattern=r"^[0-9a-f]{64}$", max_length=64)


class AuthoringRejectRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    proposal_digest: str = Field(pattern=r"^[0-9a-f]{64}$", max_length=64)


class AuthoringAdoptionResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    ok: Literal[True] = True
    replayed: bool
    adoption_receipt: AuthoringReceiptResponse
    draft_state: Literal["active", "expired"]
    draft_revision: int
    expires_at: datetime
    purge_after: datetime


class AuthoringValidateRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    expected_revision: int = Field(ge=0, le=1_000_000)


class AuthoringValidationResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    draft_revision: int
    ok: bool
    workflow_id: str | None = None
    value: None = None
    snapshot: ValidatedInputSnapshotResponse | None = None
    issues: list[IssueResponse] = Field(default_factory=list)


class ProjectSummaryResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    project_id: ProjectId
    display_name: str


class ProjectListResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    ok: Literal[True] = True
    projects: list[ProjectSummaryResponse]
    next_cursor: str | None = None


class ProjectSampleSummaryResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    sample_id: SampleId
    display_name: str
    sample_revision_id: SampleRevisionId


class ProjectSampleListResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    ok: Literal[True] = True
    project_id: ProjectId
    samples: list[ProjectSampleSummaryResponse]
    next_cursor: str | None = None
