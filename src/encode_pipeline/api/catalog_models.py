"""Request and response models for the Slice A experiment catalog."""

from __future__ import annotations

from typing import Any

from pydantic import BaseModel, ConfigDict, Field


class CatalogIssueResponse(BaseModel):
    """Stable catalog failure fields. Codes stay lowercase except auth codes."""

    model_config = ConfigDict(extra="forbid")

    code: str
    message: str
    severity: str = "error"
    path: str | None = None
    source: str | None = None
    hint: str | None = None


class CatalogErrorResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    issues: list[CatalogIssueResponse]


class CatalogCreateTermRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    term_type: str
    display_name: str
    namespace: str | None = None
    accession: str | None = None
    aliases: list[str] = Field(default_factory=list)


class CatalogReviseTermRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    expected_revision_number: Any
    change_kind: str
    display_name: str
    namespace: str | None = None
    accession: str | None = None
    aliases: list[str] = Field(default_factory=list)


class CatalogArchiveRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    expected_revision_number: Any


class CatalogCreateExperimentRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    stable_key: str
    display_name: str
    assay: Any
    target: Any
    batch: Any
    samples: Any


class CatalogReviseExperimentRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    expected_revision_number: Any
    display_name: str
    assay: Any
    target: Any
    batch: Any
    samples: Any


class CatalogReviseSampleRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    expected_revision_number: Any
    catalog: Any
    display_name: str | None = None


class CatalogSlotResponse(BaseModel):
    """Known, negative, or ambiguous slot. Extra keys carry pins and candidates."""

    model_config = ConfigDict(extra="allow")

    status: str


class CatalogSampleMaterialResponse(BaseModel):
    model_config = ConfigDict(extra="allow")

    sample_id: str
    sample_revision_id: str
    catalog_state: str | None = None


class CatalogTermSummaryResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    term_id: str
    term_type: str
    term_revision_id: str
    revision_number: int
    display_name: str
    archived_at: str | None = None


class CatalogTermRevisionResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    term_revision_id: str
    revision_number: int
    display_name: str
    namespace: str | None = None
    accession: str | None = None
    aliases: list[str]
    change_kind: str
    created_at: str


class CatalogTermResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    term_id: str
    term_type: str
    archived_at: str | None = None
    revision: CatalogTermRevisionResponse


class CatalogTermListResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    terms: list[CatalogTermSummaryResponse]
    next_cursor: str | None = None


class CatalogExperimentResponse(BaseModel):
    model_config = ConfigDict(extra="allow")

    experiment_id: str
    experiment_revision_id: str
    revision_number: int
    display_name: str
    assay: CatalogSlotResponse
    target: CatalogSlotResponse
    batch: CatalogSlotResponse
    samples: list[dict[str, Any]]


class CatalogExperimentListResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    project_id: str
    revision_view: str
    experiments: list[CatalogExperimentResponse]
    next_cursor: str | None = None


class CatalogExperimentRevisionListResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    experiment_revisions: list[CatalogExperimentResponse]
    next_cursor: str | None = None


class CatalogArchiveExperimentResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    experiment_id: str
    archived_at: str
    revision_number: int


class CatalogSampleRevisionResponse(BaseModel):
    model_config = ConfigDict(extra="allow")

    sample_id: str
    sample_revision_id: str
    revision_number: int
    display_name: str
    created_at: str


class CatalogSampleListResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    project_id: str
    samples: list[CatalogSampleRevisionResponse]
    next_cursor: str | None = None


class CatalogSampleRevisionListResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    sample_revisions: list[CatalogSampleRevisionResponse]
    next_cursor: str | None = None


class CatalogSampleRunListResponse(BaseModel):
    model_config = ConfigDict(extra="allow")

    project_id: str
    sample_revision_id: str
    results: list[dict[str, Any]]
    next_cursor: str | None = None


class CatalogExperimentRunListResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    project_id: str
    experiment_revision_id: str
    trace_status: str
    reason: str
    results: list[dict[str, Any]]
    next_cursor: str | None = None
