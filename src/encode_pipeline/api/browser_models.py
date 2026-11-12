"""Public response models for the WashU browser read API."""

from __future__ import annotations

from typing import Any

from pydantic import BaseModel, ConfigDict, Field


class BrowserIssueContextResponse(BaseModel):
    """Allowlisted failure context. Host paths are not a field."""

    model_config = ConfigDict(extra="ignore")

    reason_code: str | None = None


class BrowserIssueResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    code: str
    message: str
    severity: str = "error"
    path: str
    source: str = "api"
    context: BrowserIssueContextResponse = Field(
        default_factory=BrowserIssueContextResponse
    )


class BrowserErrorResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    ok: bool = False
    run_id: str
    issues: list[BrowserIssueResponse]


class BrowserProvenanceResponse(BaseModel):
    model_config = ConfigDict(extra="ignore")

    converter_id: str | None = None
    kind: str | None = None
    source_artifact_id: str | None = None
    source_revision: str | None = None
    source_generation: str | None = None
    source_output_type: str | None = None
    coordinate_convention: str | None = None
    score_meaning: str | None = None
    score_policy: str | None = None
    pair_count: int | None = None
    anchor_row_count: int | None = None
    data_sha256: str | None = None
    index_sha256: str | None = None
    data_size: int | None = None
    index_size: int | None = None
    suggested_region: str | None = None


class BrowserTrackResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    artifact_id: str
    name: str
    output_type: str
    revision: str
    role: str
    kind: str
    wash_u_type: str
    label: str
    detail: str
    score_meaning: str
    source_form: str
    preparation: str
    content_path: str | None = None
    index_path: str | None = None
    provenance: BrowserProvenanceResponse | None = None


class BrowserUnavailableTrackResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    artifact_id: str
    name: str
    output_type: str
    reason: str | None = None
    detail: str


class BrowserReferenceResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    status: str
    profile_id: str | None = None
    revision_id: str | None = None
    identity_sha256: str | None = None
    binding_digest: str | None = None
    assembly_label: str | None = None
    viewer_scaffold: str | None = None


class BrowserSessionResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    ok: bool
    run_id: str
    workflow_id: str
    artifact_generation: str | None
    reference: BrowserReferenceResponse
    annotation_status: str
    interaction_status: str
    tracks: list[BrowserTrackResponse]
    unavailable: list[BrowserUnavailableTrackResponse]
    draw_record_limit: int
    issues: list[BrowserIssueResponse] = Field(default_factory=list)


class BrowserPrepareRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    generation: str = Field(pattern=r"^artifactgen-[0-9a-f]{64}$", max_length=76)
    revision: str = Field(pattern=r"^artifactrev-[0-9a-f]{64}$", max_length=76)


class BrowserPrepareResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    ok: bool
    run_id: str
    track: BrowserTrackResponse
    issues: list[BrowserIssueResponse] = Field(default_factory=list)


class BrowserRegionSummaryResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    ok: bool
    kind: str
    records: int | None
    truncated: bool
    complete: bool
    detail: str
    issues: list[BrowserIssueResponse] = Field(default_factory=list)


def browser_error(
    *,
    run_id: str,
    code: str,
    message: str,
    path: str,
    reason: str | None,
) -> dict[str, Any]:
    context = {"reason_code": reason} if reason else {}
    return BrowserErrorResponse(
        run_id=run_id,
        issues=[
            BrowserIssueResponse(
                code=code,
                message=message,
                path=path,
                context=BrowserIssueContextResponse(**context),
            )
        ],
    ).model_dump(mode="json")
