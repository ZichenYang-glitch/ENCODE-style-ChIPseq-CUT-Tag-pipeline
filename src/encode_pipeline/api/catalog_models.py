"""Request and response models for the Slice A experiment catalog."""

from __future__ import annotations

from typing import Annotated, Any, Literal

from pydantic import BaseModel, ConfigDict, Field, WithJsonSchema, model_validator


class CatalogIssueResponse(BaseModel):
    """Stable catalog failure fields. Codes stay lowercase except auth codes."""

    model_config = ConfigDict(extra="forbid")

    code: str
    message: str
    severity: str = "error"
    path: str | None = None
    source: str | None = None
    hint: str | None = None


class CatalogRevisionConflictDetail(BaseModel):
    """One selected target whose head revision differs from the proposal baseline."""

    model_config = ConfigDict(extra="forbid")

    target_kind: str
    target_id: str
    expected_revision_number: int
    current_revision_number: int
    current_revision_id: str
    expected_revision_id: str | None = None


class CatalogErrorResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    issues: list[CatalogIssueResponse]
    conflicts: list[CatalogRevisionConflictDetail] | None = None


_REQUEST_PIN_SCHEMA: dict[str, Any] = {
    "type": "object",
    "additionalProperties": False,
    "required": ["term_id", "term_revision_id"],
    "properties": {
        "term_id": {"type": "string"},
        "term_revision_id": {"type": "string"},
    },
}
_NEGATIVE_STATUS = ["unknown", "not_applicable", "not_detected"]
_REQUEST_SINGLE_SCHEMA: dict[str, Any] = {
    "oneOf": [
        {
            "type": "object",
            "additionalProperties": False,
            "required": ["status", "term_id", "term_revision_id"],
            "properties": {
                "status": {"const": "known"},
                "term_id": {"type": "string"},
                "term_revision_id": {"type": "string"},
            },
        },
        {
            "type": "object",
            "additionalProperties": False,
            "required": ["status"],
            "properties": {"status": {"enum": _NEGATIVE_STATUS}},
        },
        {
            "type": "object",
            "additionalProperties": False,
            "required": ["status"],
            "properties": {
                "status": {"const": "ambiguous"},
                "candidates": {
                    "type": "array",
                    "items": _REQUEST_PIN_SCHEMA,
                },
            },
        },
    ]
}
_REQUEST_MULTI_SCHEMA: dict[str, Any] = {
    "oneOf": [
        {
            "type": "object",
            "additionalProperties": False,
            "required": ["status", "values"],
            "properties": {
                "status": {"const": "known"},
                "values": {
                    "type": "array",
                    "minItems": 1,
                    "items": _REQUEST_PIN_SCHEMA,
                },
            },
        },
        {
            "type": "object",
            "additionalProperties": False,
            "required": ["status", "values"],
            "properties": {
                "status": {"enum": _NEGATIVE_STATUS},
                "values": {
                    "type": "array",
                    "maxItems": 0,
                    "items": _REQUEST_PIN_SCHEMA,
                },
            },
        },
        {
            "type": "object",
            "additionalProperties": False,
            "required": ["status", "values"],
            "properties": {
                "status": {"const": "ambiguous"},
                "values": {
                    "type": "array",
                    "maxItems": 0,
                    "items": _REQUEST_PIN_SCHEMA,
                },
                "candidates": {
                    "type": "array",
                    "items": _REQUEST_PIN_SCHEMA,
                },
            },
        },
    ]
}
_REQUEST_SAMPLE_PIN_SCHEMA: dict[str, Any] = {
    "type": "object",
    "additionalProperties": False,
    "required": ["sample_revision_id", "replicate_role"],
    "properties": {
        "sample_revision_id": {"type": "string"},
        "replicate_group": {"type": "string"},
        "replicate_role": _REQUEST_SINGLE_SCHEMA,
    },
}
_REQUEST_CATALOG_SCHEMA: dict[str, Any] = {
    "anyOf": [
        {"type": "null"},
        {
            "type": "object",
            "additionalProperties": False,
            "required": [
                "schema_version",
                "species",
                "tissue",
                "cell_type",
                "cell_line",
                "treatment",
                "tag",
            ],
            "properties": {
                "schema_version": {"type": "integer", "enum": [1]},
                "species": _REQUEST_SINGLE_SCHEMA,
                "tissue": _REQUEST_SINGLE_SCHEMA,
                "cell_type": _REQUEST_MULTI_SCHEMA,
                "cell_line": _REQUEST_MULTI_SCHEMA,
                "treatment": _REQUEST_MULTI_SCHEMA,
                "tag": _REQUEST_MULTI_SCHEMA,
            },
        },
    ]
}
_REQUEST_SAMPLES_SCHEMA: dict[str, Any] = {
    "type": "array",
    "minItems": 1,
    "items": _REQUEST_SAMPLE_PIN_SCHEMA,
}

CatalogRequestSingle = Annotated[Any, WithJsonSchema(_REQUEST_SINGLE_SCHEMA)]
CatalogRequestMulti = Annotated[Any, WithJsonSchema(_REQUEST_MULTI_SCHEMA)]
CatalogRequestSamples = Annotated[Any, WithJsonSchema(_REQUEST_SAMPLES_SCHEMA)]
CatalogRequestCatalog = Annotated[Any, WithJsonSchema(_REQUEST_CATALOG_SCHEMA)]
RevisionNumber = Annotated[int, Field(ge=1)]


class CatalogCreateTermRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    term_type: str
    display_name: str
    namespace: str | None = None
    accession: str | None = None
    aliases: list[str] = Field(default_factory=list)


class CatalogReviseTermRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    expected_revision_number: RevisionNumber
    change_kind: str
    display_name: str
    namespace: str | None = None
    accession: str | None = None
    aliases: list[str] = Field(default_factory=list)


class CatalogArchiveRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    expected_revision_number: RevisionNumber


class CatalogCreateExperimentRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    stable_key: str
    display_name: str
    assay: CatalogRequestSingle
    target: CatalogRequestSingle
    batch: CatalogRequestSingle
    samples: CatalogRequestSamples


class CatalogReviseExperimentRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    expected_revision_number: RevisionNumber
    display_name: str
    assay: CatalogRequestSingle
    target: CatalogRequestSingle
    batch: CatalogRequestSingle
    samples: CatalogRequestSamples


class CatalogReviseSampleRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    expected_revision_number: RevisionNumber
    catalog: CatalogRequestCatalog
    display_name: str | None = None


class CatalogTermPinResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    term_id: str
    term_revision_id: str
    display_name: str


class CatalogKnownSingleSlotResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    status: Literal["known"]
    term_id: str
    term_revision_id: str
    display_name: str


class CatalogNegativeSingleSlotResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    status: Literal["unknown", "not_applicable", "not_detected"]


class CatalogAmbiguousSingleSlotResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    status: Literal["ambiguous"]
    candidates: list[CatalogTermPinResponse] | None = None


CatalogSingleSlotResponse = Annotated[
    CatalogKnownSingleSlotResponse
    | CatalogNegativeSingleSlotResponse
    | CatalogAmbiguousSingleSlotResponse,
    Field(union_mode="left_to_right"),
]


class CatalogKnownMultiSlotResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    status: Literal["known"]
    values: list[CatalogTermPinResponse] = Field(min_length=1)


class CatalogNegativeMultiSlotResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    status: Literal["unknown", "not_applicable", "not_detected"]
    values: list[CatalogTermPinResponse] = Field(max_length=0)


class CatalogAmbiguousMultiSlotResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    status: Literal["ambiguous"]
    values: list[CatalogTermPinResponse] = Field(max_length=0)
    candidates: list[CatalogTermPinResponse] | None = None


CatalogMultiSlotResponse = Annotated[
    CatalogKnownMultiSlotResponse
    | CatalogNegativeMultiSlotResponse
    | CatalogAmbiguousMultiSlotResponse,
    Field(union_mode="left_to_right"),
]


class _ExperimentSampleCore(BaseModel):
    model_config = ConfigDict(extra="forbid")

    sample_id: str
    sample_revision_id: str
    ordinal: int
    replicate_group: str | None = None
    replicate_role: CatalogSingleSlotResponse
    linked_sample_revision_id: str | None = None
    current_sample_revision_id: str | None = None


class CatalogClassifiedExperimentSampleResponse(_ExperimentSampleCore):
    species: CatalogSingleSlotResponse
    tissue: CatalogSingleSlotResponse
    cell_type: CatalogMultiSlotResponse
    cell_line: CatalogMultiSlotResponse
    treatment: CatalogMultiSlotResponse
    tag: CatalogMultiSlotResponse


class CatalogUnclassifiedExperimentSampleResponse(_ExperimentSampleCore):
    catalog_state: Literal["absent"]


CatalogExperimentSampleResponse = Annotated[
    CatalogClassifiedExperimentSampleResponse
    | CatalogUnclassifiedExperimentSampleResponse,
    Field(union_mode="left_to_right"),
]


class CatalogTermSummaryResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    term_id: str
    term_type: str
    term_revision_id: str
    revision_number: RevisionNumber
    display_name: str
    archived_at: str | None = None


class CatalogTermRevisionResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    term_revision_id: str
    revision_number: RevisionNumber
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
    model_config = ConfigDict(extra="forbid")

    experiment_id: str
    experiment_revision_id: str
    revision_number: RevisionNumber
    display_name: str
    assay: CatalogSingleSlotResponse
    target: CatalogSingleSlotResponse
    batch: CatalogSingleSlotResponse
    samples: list[CatalogExperimentSampleResponse]


class CatalogExperimentListResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    ok: Literal[True]
    project_id: str
    revision_view: str
    experiments: list[CatalogExperimentResponse]
    next_cursor: str | None


class CatalogExperimentRevisionListResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    experiment_revisions: list[CatalogExperimentResponse]
    next_cursor: str | None = None


class CatalogArchiveExperimentResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    experiment_id: str
    archived_at: str
    revision_number: RevisionNumber


class _SampleRevisionCore(BaseModel):
    model_config = ConfigDict(extra="forbid")

    sample_id: str
    sample_revision_id: str
    revision_number: RevisionNumber
    display_name: str
    created_at: str


class CatalogUnclassifiedSampleRevisionResponse(_SampleRevisionCore):
    catalog_state: Literal["absent"]


class CatalogClassifiedSampleRevisionResponse(_SampleRevisionCore):
    catalog_state: Literal["present"]
    species: CatalogSingleSlotResponse
    tissue: CatalogSingleSlotResponse
    cell_type: CatalogMultiSlotResponse
    cell_line: CatalogMultiSlotResponse
    treatment: CatalogMultiSlotResponse
    tag: CatalogMultiSlotResponse


CatalogSampleRevisionResponse = Annotated[
    CatalogClassifiedSampleRevisionResponse | CatalogUnclassifiedSampleRevisionResponse,
    Field(union_mode="left_to_right"),
]


class CatalogSampleListResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    project_id: str
    samples: list[CatalogSampleRevisionResponse]
    next_cursor: str | None = None


class CatalogSampleRevisionListResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    sample_revisions: list[CatalogSampleRevisionResponse]
    next_cursor: str | None = None


class CatalogAssociatedRunSampleResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    sample_id: str
    sample_revision_id: str
    revision_number: RevisionNumber
    ordinal: int


class CatalogRunSampleBindingResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    binding_mode: str
    provenance: str
    associated_run_samples: list[CatalogAssociatedRunSampleResponse]


class CatalogPublishedRunResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    run_id: str
    provenance: str
    publication_evidence: Literal["present"]
    artifact_association_kind: Literal["run_input_scope"]
    workflow_id: str
    artifact_id: str
    output_type: str
    artifact_generation: str
    artifact_revision: str
    published_at: str
    current_artifact_generation: str
    generation_status: Literal["current", "superseded"]
    run_sample_binding: CatalogRunSampleBindingResponse


class CatalogUnpublishedRunResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    run_id: str
    provenance: str
    publication_evidence: Literal["absent"]


CatalogRunResultResponse = Annotated[
    CatalogPublishedRunResponse | CatalogUnpublishedRunResponse,
    Field(union_mode="left_to_right"),
]


class CatalogSampleRunListResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    ok: Literal[True]
    project_id: str
    sample_revision_id: str
    results: list[CatalogRunResultResponse]
    next_cursor: str | None


class CatalogExperimentRunListResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    ok: Literal[True]
    project_id: str
    experiment_revision_id: str
    trace_status: Literal["unresolved", "experiment_snapshot"]
    reason: Literal["experiment_snapshot_binding_absent"] | None = None
    results: list[CatalogRunResultResponse]
    next_cursor: str | None

    @model_validator(mode="after")
    def trace_shape(self) -> "CatalogExperimentRunListResponse":
        if self.trace_status == "unresolved":
            if self.reason != "experiment_snapshot_binding_absent" or self.results:
                raise ValueError(
                    "unresolved experiment trace is an empty unbound envelope"
                )
        elif self.reason is not None or any(
            item.publication_evidence != "present" for item in self.results
        ):
            raise ValueError("bound experiment trace lists published rows only")
        return self


class CatalogRegisterProposalRequest(BaseModel):
    """Offline suggestion package. The server assigns the proposal id and digest."""

    model_config = ConfigDict(extra="forbid")

    source: Literal["offline_suggestion", "agent", "manual"]
    rule_version: str | None = None
    vocabulary_version: str | None = None
    projection_version: str | None = None
    items: list[dict[str, Any]]


class CatalogPreviewProposalRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    selected_item_ids: list[str]


class CatalogAdoptProposalRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")

    preview_digest: str
    selected_item_ids: list[str]
    expected_revisions: dict[str, Any]


class CatalogProposalItemView(BaseModel):
    model_config = ConfigDict(extra="forbid")

    item_id: str
    target_kind: str
    target_id: str
    base_revision_id: str
    base_revision_number: int
    current_revision_id: str
    current_revision_number: int
    allowed_fields: list[str]
    before: dict[str, Any]
    after: dict[str, Any]
    citations: dict[str, Any]


class CatalogProposalResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    ok: Literal[True]
    contract_version: str
    proposal_id: str
    project_id: str
    source: str
    proposal_state: str
    adoptable_by_caller: bool
    rule_version: str
    vocabulary_version: str
    projection_version: str
    proposal_digest: str
    created_at: str
    expires_at: str
    retention_until: str
    adoption_receipt: dict[str, Any] | None = None
    items: list[CatalogProposalItemView]


class CatalogPreviewItemResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    item_id: str
    target_kind: str
    target_id: str
    expected_revision_number: int
    base_revision_id: str
    current_revision_number: int
    current_revision_id: str
    conflict_detected: bool
    allowed_fields: list[str]
    before: dict[str, Any]
    after: dict[str, Any]
    citations: dict[str, Any]


class CatalogProposalPreviewResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    ok: Literal[True]
    contract_version: str
    proposal_id: str
    project_id: str
    proposal_digest: str
    preview_digest: str
    selected_item_ids: list[str]
    excluded_item_ids: list[str]
    expected_revisions: dict[str, int]
    has_conflicts: bool
    conflicts: list[CatalogRevisionConflictDetail]
    expires_at: str
    items: list[CatalogPreviewItemResponse]


class CatalogAppliedTargetResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    item_id: str
    target_kind: str
    target_id: str
    base_revision_id: str
    base_revision_number: int
    result_revision_id: str
    result_revision_number: int
    changed_fields: list[str]


class CatalogAdoptionReceiptResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    ok: Literal[True]
    replayed: bool
    contract_version: str
    receipt_id: str
    proposal_id: str
    project_id: str
    proposal_digest: str
    preview_digest: str
    actor_user_id: str
    idempotency_key: str
    selected_item_ids: list[str]
    excluded_item_ids: list[str]
    applied_targets: list[CatalogAppliedTargetResponse]
    adopted_at: str


class CatalogAdoptionReceiptListResponse(BaseModel):
    model_config = ConfigDict(extra="forbid")

    ok: Literal[True]
    project_id: str
    receipts: list[CatalogAdoptionReceiptResponse]


class CatalogPublicationExperimentTraceResponse(BaseModel):
    """One published result traced back to the revision copied onto its run."""

    model_config = ConfigDict(extra="forbid")

    ok: Literal[True]
    project_id: str
    run_id: str
    artifact_id: str
    artifact_generation: str
    trace_status: Literal["unresolved", "experiment_snapshot"]
    reason: (
        Literal[
            "experiment_snapshot_binding_absent",
            "experiment_revision_evidence_inconsistent",
        ]
        | None
    ) = None
    experiment_revision_id: str | None = None
    payload_digest: str | None = None
    result: CatalogPublishedRunResponse

    @model_validator(mode="after")
    def trace_shape(self) -> "CatalogPublicationExperimentTraceResponse":
        if (
            self.result.run_id != self.run_id
            or self.result.artifact_id != self.artifact_id
            or self.result.artifact_generation != self.artifact_generation
            or self.result.artifact_association_kind != "run_input_scope"
        ):
            raise ValueError("publication trace result does not match its envelope")
        if self.trace_status == "unresolved":
            if (
                self.reason
                not in {
                    "experiment_snapshot_binding_absent",
                    "experiment_revision_evidence_inconsistent",
                }
                or self.experiment_revision_id is not None
                or self.payload_digest is not None
            ):
                raise ValueError("unresolved publication trace has no revision")
        elif (
            self.reason is not None
            or self.experiment_revision_id is None
            or self.payload_digest is None
        ):
            raise ValueError("frozen publication trace requires the copied revision")
        return self
