"""Restricted authoring session and project catalog routes."""

from __future__ import annotations

from fastapi import APIRouter, Depends, Query, Request
from fastapi.responses import JSONResponse

from encode_pipeline.api.authoring_models import (
    AuthoringAdoptRequest,
    AuthoringAdoptionResponse,
    AuthoringCapabilityResponse,
    AuthoringConnectionResponse,
    AuthoringDraftDocument,
    AuthoringDraftWriteRequest,
    AuthoringErrorResponse,
    AuthoringExplanationResponse,
    AuthoringProposalResponse,
    AuthoringReceiptListResponse,
    AuthoringReceiptResponse,
    AuthoringRejectRequest,
    AuthoringSessionCreateRequest,
    AuthoringSessionListResponse,
    AuthoringSessionResponse,
    AuthoringSessionSummary,
    AuthoringTurnCreateRequest,
    AuthoringTurnResponse,
    AuthoringValidateRequest,
    AuthoringValidationResponse,
    ProjectListResponse,
    ProjectSampleListResponse,
    ProjectSampleSummaryResponse,
    ProjectSummaryResponse,
)
from encode_pipeline.api.dependencies import require_authoring_actor
from encode_pipeline.api.models import (
    IssueResponse,
    ValidatedInputSnapshotResponse,
    ValidationResponse,
)
from encode_pipeline.platform.authoring import AuthoringActor, AuthoringRejected
from encode_pipeline.services.authoring import (
    AdoptionView,
    AuthoringService,
    PreparedValidation,
    ReceiptPage,
    SessionView,
    TurnView,
)
from encode_pipeline.services.validated_inputs import ValidatedInputService


_INTERNAL_ERROR = {500: {"model": ValidationResponse}}

projects_router = APIRouter(prefix="/projects", tags=["authoring"])
router = APIRouter(prefix="/agent/authoring-sessions", tags=["authoring"])


def get_authoring_service(request: Request) -> AuthoringService:
    return request.app.state.authoring_service


def get_validated_inputs(request: Request) -> ValidatedInputService:
    return request.app.state.validated_input_service


@projects_router.get(
    "",
    response_model=ProjectListResponse,
    operation_id="authoringListProjects",
    responses=_INTERNAL_ERROR,
)
def list_projects(
    actor: AuthoringActor = Depends(require_authoring_actor),
    service: AuthoringService = Depends(get_authoring_service),
    cursor: str | None = Query(default=None, max_length=256),
    limit: int = Query(default=20, ge=1, le=50),
) -> ProjectListResponse:
    page = service.list_projects(actor, cursor=cursor, limit=limit)
    return ProjectListResponse(
        projects=[
            ProjectSummaryResponse(project_id=project_id, display_name=name)
            for project_id, name in page.projects
        ],
        next_cursor=page.next_cursor,
    )


@projects_router.get(
    "/{project_id}/samples",
    response_model=ProjectSampleListResponse,
    operation_id="authoringListProjectSamples",
    responses=_INTERNAL_ERROR,
)
def list_project_samples(
    project_id: str,
    actor: AuthoringActor = Depends(require_authoring_actor),
    service: AuthoringService = Depends(get_authoring_service),
    cursor: str | None = Query(default=None, max_length=256),
    limit: int = Query(default=20, ge=1, le=50),
) -> ProjectSampleListResponse:
    page = service.list_samples(actor, project_id, cursor=cursor, limit=limit)
    return ProjectSampleListResponse(
        project_id=page.project_id,
        samples=[
            ProjectSampleSummaryResponse(
                sample_id=sample_id,
                display_name=name,
                sample_revision_id=revision_id,
            )
            for sample_id, name, revision_id in page.samples
        ],
        next_cursor=page.next_cursor,
    )


@router.get(
    "",
    response_model=AuthoringSessionListResponse,
    operation_id="authoringListSessions",
    responses=_INTERNAL_ERROR,
)
def list_sessions(
    actor: AuthoringActor = Depends(require_authoring_actor),
    service: AuthoringService = Depends(get_authoring_service),
    project_id: str | None = Query(default=None, max_length=36),
    workflow_id: str | None = Query(default=None, max_length=255),
    cursor: str | None = Query(default=None, max_length=256),
    limit: int = Query(default=20, ge=1, le=50),
) -> AuthoringSessionListResponse:
    sessions, next_cursor = service.list_sessions(
        actor,
        project_id=project_id,
        workflow_id=workflow_id,
        cursor=cursor,
        limit=limit,
    )
    return AuthoringSessionListResponse(
        sessions=[_summary(session) for session in sessions],
        next_cursor=next_cursor,
    )


@router.post(
    "",
    response_model=AuthoringSessionResponse,
    operation_id="authoringCreateSession",
    responses=_INTERNAL_ERROR,
)
def create_session(
    body: AuthoringSessionCreateRequest,
    actor: AuthoringActor = Depends(require_authoring_actor),
    service: AuthoringService = Depends(get_authoring_service),
) -> AuthoringSessionResponse:
    return _session(
        service.create_session(
            actor,
            project_id=body.project_id,
            workflow_id=body.workflow_id,
        )
    )


@router.get(
    "/{sid}",
    response_model=AuthoringSessionResponse,
    operation_id="authoringGetSession",
    responses=_INTERNAL_ERROR,
)
def get_session(
    sid: str,
    actor: AuthoringActor = Depends(require_authoring_actor),
    service: AuthoringService = Depends(get_authoring_service),
) -> AuthoringSessionResponse:
    return _session(service.get_session(actor, sid))


@router.get(
    "/{sid}/receipts",
    response_model=AuthoringReceiptListResponse,
    operation_id="authoringListReceipts",
    responses=_INTERNAL_ERROR,
)
def list_receipts(
    sid: str,
    actor: AuthoringActor = Depends(require_authoring_actor),
    service: AuthoringService = Depends(get_authoring_service),
    cursor: str | None = Query(default=None, max_length=256),
    limit: int = Query(default=20, ge=1, le=50),
) -> AuthoringReceiptListResponse:
    return _receipts(
        service.list_adoption_receipts(actor, sid, cursor=cursor, limit=limit)
    )


@router.put(
    "/{sid}/draft",
    response_model=AuthoringSessionResponse,
    operation_id="authoringSaveDraft",
    responses=_INTERNAL_ERROR,
)
def save_draft(
    sid: str,
    body: AuthoringDraftWriteRequest,
    actor: AuthoringActor = Depends(require_authoring_actor),
    service: AuthoringService = Depends(get_authoring_service),
) -> AuthoringSessionResponse:
    return _session(
        service.save_draft(
            actor,
            sid,
            expected_revision=body.expected_revision,
            content=body.model_dump(),
        )
    )


@router.post(
    "/{sid}/connection",
    response_model=AuthoringSessionResponse,
    operation_id="authoringOpenConnection",
    responses=_INTERNAL_ERROR,
)
def open_connection(
    sid: str,
    actor: AuthoringActor = Depends(require_authoring_actor),
    service: AuthoringService = Depends(get_authoring_service),
) -> AuthoringSessionResponse:
    return _session(service.open_connection(actor, sid))


@router.delete(
    "/{sid}/connection",
    response_model=AuthoringSessionResponse,
    operation_id="authoringCloseConnection",
    responses=_INTERNAL_ERROR,
)
def close_connection(
    sid: str,
    actor: AuthoringActor = Depends(require_authoring_actor),
    service: AuthoringService = Depends(get_authoring_service),
) -> AuthoringSessionResponse:
    return _session(service.close_connection(actor, sid))


@router.post(
    "/{sid}/turns",
    response_model=AuthoringTurnResponse,
    operation_id="authoringCreateTurn",
    responses=_INTERNAL_ERROR,
)
def create_turn(
    sid: str,
    body: AuthoringTurnCreateRequest,
    actor: AuthoringActor = Depends(require_authoring_actor),
    service: AuthoringService = Depends(get_authoring_service),
) -> AuthoringTurnResponse:
    return _turn(
        service.start_turn(
            actor,
            sid,
            client_turn_id=body.client_turn_id,
            expected_revision=body.expected_revision,
            question=body.question.model_dump(),
        )
    )


@router.get(
    "/{sid}/turns/{turn_id}",
    response_model=AuthoringTurnResponse,
    operation_id="authoringGetTurn",
    responses=_INTERNAL_ERROR,
)
def get_turn(
    sid: str,
    turn_id: str,
    actor: AuthoringActor = Depends(require_authoring_actor),
    service: AuthoringService = Depends(get_authoring_service),
) -> AuthoringTurnResponse:
    return _turn(service.get_turn(actor, sid, turn_id))


@router.post(
    "/{sid}/proposals/{proposal_id}/adopt",
    response_model=AuthoringAdoptionResponse,
    operation_id="authoringAdoptProposal",
    responses=_INTERNAL_ERROR,
)
def adopt_proposal(
    sid: str,
    proposal_id: str,
    body: AuthoringAdoptRequest,
    actor: AuthoringActor = Depends(require_authoring_actor),
    service: AuthoringService = Depends(get_authoring_service),
) -> AuthoringAdoptionResponse:
    return _adoption(
        service.adopt(
            actor,
            sid,
            proposal_id,
            expected_revision=body.expected_revision,
            proposal_digest_value=body.proposal_digest,
        )
    )


@router.post(
    "/{sid}/proposals/{proposal_id}/reject",
    response_model=AuthoringTurnResponse,
    operation_id="authoringRejectProposal",
    responses=_INTERNAL_ERROR,
)
def reject_proposal(
    sid: str,
    proposal_id: str,
    body: AuthoringRejectRequest,
    actor: AuthoringActor = Depends(require_authoring_actor),
    service: AuthoringService = Depends(get_authoring_service),
) -> AuthoringTurnResponse:
    return _turn(
        service.reject(
            actor,
            sid,
            proposal_id,
            proposal_digest_value=body.proposal_digest,
        )
    )


@router.post(
    "/{sid}/validate",
    response_model=AuthoringValidationResponse,
    operation_id="authoringValidateDraft",
    responses={
        404: {"model": AuthoringValidationResponse},
        409: {"model": AuthoringValidationResponse},
        410: {"model": AuthoringErrorResponse},
        413: {"model": AuthoringErrorResponse},
        422: {"model": AuthoringErrorResponse},
        500: {"model": AuthoringValidationResponse},
        503: {"model": AuthoringValidationResponse},
    },
)
def validate_draft(
    sid: str,
    body: AuthoringValidateRequest,
    actor: AuthoringActor = Depends(require_authoring_actor),
    service: AuthoringService = Depends(get_authoring_service),
    validated: ValidatedInputService = Depends(get_validated_inputs),
) -> AuthoringValidationResponse | JSONResponse:
    prepared = service.validate(
        actor,
        sid,
        expected_revision=body.expected_revision,
    )
    status, response = _validation_response(prepared, validated)
    if status == 200:
        return response
    return JSONResponse(status_code=status, content=response.model_dump(mode="json"))


def _session(view: SessionView) -> AuthoringSessionResponse:
    draft = None
    if view.draft is not None:
        draft = AuthoringDraftDocument.model_validate(view.draft)
    return AuthoringSessionResponse(
        sid=view.sid,
        project_id=view.project_id,
        workflow_id=view.workflow_id,
        draft_state=view.draft_state,  # type: ignore[arg-type]
        draft_revision=view.draft_revision,
        expires_at=view.expires_at,
        purge_after=view.purge_after,
        schema_digest=view.schema_digest,
        projection_digest=view.projection_digest,
        draft=draft,
        authoring_capability=AuthoringCapabilityResponse.model_validate(
            view.capability
        ),
        connection=AuthoringConnectionResponse(
            state=view.connection_state,  # type: ignore[arg-type]
            reason_code=view.connection_reason_code,
        ),
        requires_new_connection=view.requires_new_connection,
        changed=view.changed,
    )


def _summary(view: SessionView) -> AuthoringSessionSummary:
    return AuthoringSessionSummary(
        sid=view.sid,
        project_id=view.project_id,
        workflow_id=view.workflow_id,
        draft_state=view.draft_state,  # type: ignore[arg-type]
        draft_revision=view.draft_revision,
        expires_at=view.expires_at,
        purge_after=view.purge_after,
        connection=AuthoringConnectionResponse(
            state=view.connection_state,  # type: ignore[arg-type]
            reason_code=view.connection_reason_code,
        ),
        requires_new_connection=view.requires_new_connection,
    )


def _receipts(page: ReceiptPage) -> AuthoringReceiptListResponse:
    return AuthoringReceiptListResponse(
        sid=page.sid,
        draft_state=page.draft_state,  # type: ignore[arg-type]
        receipts=[
            {
                "turn_id": entry.turn_id,
                "adoption_receipt": entry.receipt,
            }
            for entry in page.entries
        ],
        next_cursor=page.next_cursor,
    )


def _turn(view: TurnView) -> AuthoringTurnResponse:
    proposal = None
    if view.proposal is not None:
        proposal = AuthoringProposalResponse.model_validate(view.proposal)
    explanation = None
    if view.explanation is not None:
        explanation = AuthoringExplanationResponse.model_validate(view.explanation)
    receipt = None
    if view.receipt is not None:
        receipt = AuthoringReceiptResponse.model_validate(view.receipt)
    return AuthoringTurnResponse(
        turn_id=view.turn_id,
        client_turn_id=view.client_turn_id,
        status=view.status,  # type: ignore[arg-type]
        failure_reason=view.failure_reason,
        result_kind=view.result_kind,  # type: ignore[arg-type]
        proposal_id=view.proposal_id,
        proposal=proposal,
        proposal_digest=view.proposal_digest,
        explanation=explanation,
        proposal_state=view.proposal_state,  # type: ignore[arg-type]
        adoptable=view.adoptable,
        reason_code=view.reason_code,
        proposal_expires_at=view.proposal_expires_at,
        adoption_receipt=receipt,
    )


def _adoption(view: AdoptionView) -> AuthoringAdoptionResponse:
    return AuthoringAdoptionResponse(
        replayed=view.replayed,
        adoption_receipt=AuthoringReceiptResponse.model_validate(view.receipt),
        draft_state=view.draft_state,  # type: ignore[arg-type]
        draft_revision=view.draft_revision,
        expires_at=view.expires_at,
        purge_after=view.purge_after,
    )


def _validation_response(
    prepared: PreparedValidation,
    validated: ValidatedInputService,
) -> tuple[int, AuthoringValidationResponse]:
    result = prepared.result
    snapshot = None
    if result.value is not None:
        snapshot = ValidatedInputSnapshotResponse.from_snapshot(
            result.value,
            validated.get_validated_input_binding(result.value.snapshot_id),
            validated.get_validated_input_use_binding(result.value.snapshot_id),
            validated.get_validated_reference_summary(result.value.snapshot_id),
        )
    body = AuthoringValidationResponse(
        draft_revision=prepared.draft_revision,
        ok=result.is_success,
        workflow_id=None if result.value is None else result.value.workflow_id,
        value=None,
        snapshot=snapshot,
        issues=[
            IssueResponse.model_validate(issue.to_dict()) for issue in result.issues
        ],
    )
    if result.is_failure and result.issues:
        status = {
            "WORKFLOW_NOT_FOUND": 404,
            "WORKFLOW_CAPABILITY_UNSUPPORTED": 409,
            "DATA_BINDING_SELECTION_INVALID": 409,
            "INPUT_BINDING_SELECTION_INVALID": 409,
            "INPUT_USE_CAPABILITY_UNAVAILABLE": 409,
            "VALIDATION_WORKFLOW_BUILD_CHANGED": 409,
            "VALIDATION_WORKFLOW_BUILD_UNAVAILABLE": 503,
            "VALIDATION_WORKFLOW_SCHEMA_UNAVAILABLE": 503,
            "VALIDATED_SNAPSHOT_PERSISTENCE_FAILED": 500,
            "REFERENCE_PROFILE_REQUIRED": 409,
            "REFERENCE_PROFILE_DISABLED": 409,
            "REFERENCE_PROFILE_STALE": 409,
            "REFERENCE_PROFILE_INCOMPATIBLE": 409,
            "REFERENCE_PROFILE_IDENTITY_MISMATCH": 409,
            "REFERENCE_PROFILE_UNAVAILABLE": 503,
            "REFERENCE_PROFILE_CONFIG_INVALID": 503,
            "REFERENCE_PROFILE_BINDING_INVALID": 503,
        }.get(result.issues[0].code, 200)
        return status, body
    return 200, body


def authoring_error_response(exc: AuthoringRejected) -> JSONResponse:
    issue = IssueResponse(
        code=exc.reason_code,
        message=exc.message,
        severity="error",
        path=None,
        source="authoring",
        technical_message=None,
        hint=None,
        context=dict(exc.detail),
    )
    body = AuthoringErrorResponse(reason_code=exc.reason_code, issues=[issue])
    return JSONResponse(
        status_code=exc.status_code,
        content=body.model_dump(mode="json"),
    )
