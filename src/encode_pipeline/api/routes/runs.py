"""Workflow run lifecycle API routes."""

from __future__ import annotations

from typing import Any

from fastapi import APIRouter, Depends, Query
from fastapi.responses import JSONResponse

from encode_pipeline.api.dependencies import (
    get_run_cancellation_service,
    get_run_service,
    get_run_submission_service,
    get_validated_run_creation_service,
)
from encode_pipeline.api.models import (
    IssueResponse,
    RunCreateRequest,
    RunEventResponse,
    RunEventsResponse,
    RunLogChunkResponse,
    RunLogsResponse,
    RunRecordResponse,
    RunResponse,
)
from encode_pipeline.services.run_cancellation import (
    RunCancellationConflictError,
    RunCancellationService,
    RunCancellationUnavailableError,
)
from encode_pipeline.services.runs import RunService
from encode_pipeline.services.run_submission import (
    RunBuildIdentityMissingError,
    RunNotReadyError,
    RunStartConflictError,
    RunSubmissionService,
    RunSubmissionUnavailableError,
)
from encode_pipeline.services.validated_inputs import (
    ValidatedRunCreationService,
    ValidatedSnapshotBuildUnavailableError,
    ValidatedSnapshotDataInvalidError,
    ValidatedSnapshotExpiredError,
    ValidatedSnapshotNotFoundError,
    ValidatedSnapshotReplayConflictError,
    ValidatedSnapshotStaleError,
)


router = APIRouter(tags=["runs"])


def _issue(
    code: str,
    message: str,
    *,
    path: str | None = None,
    context: dict[str, Any] | None = None,
) -> IssueResponse:
    return IssueResponse(
        code=code,
        message=message,
        severity="error",
        path=path,
        source="api",
        context=context or {},
    )


def _workflow_not_found_issue(workflow_id: str) -> IssueResponse:
    return _issue(
        code="WORKFLOW_NOT_FOUND",
        message="Workflow was not found.",
        path="workflow_id",
        context={"workflow_id": workflow_id},
    )


def _run_not_found_issue(run_id: str) -> IssueResponse:
    return _issue(
        code="RUN_NOT_FOUND",
        message="Run was not found.",
        path="run_id",
        context={"run_id": run_id},
    )


def _run_cursor_not_found_issue() -> IssueResponse:
    return _issue(
        code="RUN_CURSOR_NOT_FOUND",
        message="Cursor does not exist for this run/stream.",
        path="after",
    )


def _run_not_ready_issue(current_status: str) -> IssueResponse:
    return _issue(
        code="RUN_NOT_READY",
        message="Run must complete preflight before it can start.",
        path="run_id",
        context={"current_status": current_status},
    )


def _run_workflow_build_identity_missing_issue(
    current_status: str,
) -> IssueResponse:
    return _issue(
        code="RUN_WORKFLOW_BUILD_IDENTITY_MISSING",
        message="Run has no durable workflow build identity.",
        path="run_id",
        context={"current_status": current_status},
    )


def _run_start_conflict_issue(current_status: str) -> IssueResponse:
    return _issue(
        code="RUN_START_CONFLICT",
        message="Run execution submission conflicts with durable state.",
        path="run_id",
        context={"current_status": current_status},
    )


def _run_queue_unavailable_issue(current_status: str) -> IssueResponse:
    return _issue(
        code="RUN_QUEUE_UNAVAILABLE",
        message="The execution queue could not confirm submission.",
        path="run_id",
        context={"current_status": current_status, "retryable": True},
    )


def _run_cancellation_conflict_issue(current_status: str) -> IssueResponse:
    return _issue(
        code="RUN_CANCELLATION_CONFLICT",
        message="Run cancellation conflicts with durable execution state.",
        path="run_id",
        context={"current_status": current_status},
    )


def _run_cancellation_unavailable_issue(current_status: str) -> IssueResponse:
    return _issue(
        code="RUN_CANCELLATION_UNAVAILABLE",
        message="The execution queue could not confirm cancellation.",
        path="run_id",
        context={"current_status": current_status, "retryable": True},
    )


def _api_request_invalid_issue(
    message: str, context: dict[str, Any] | None = None
) -> IssueResponse:
    return _issue(
        code="API_REQUEST_INVALID",
        message=message,
        path="workflow_id",
        context=context or {},
    )


def _validated_snapshot_issue(code: str, message: str) -> IssueResponse:
    return _issue(
        code=code,
        message=message,
        path="snapshot_id",
    )


def _run_record_response(record: Any) -> RunRecordResponse:
    return RunRecordResponse(**record.to_dict())


def _run_event_response(event: Any) -> RunEventResponse:
    return RunEventResponse(**event.to_dict())


def _run_log_chunk_response(chunk: Any) -> RunLogChunkResponse:
    return RunLogChunkResponse(**chunk.to_dict())


@router.post(
    "/workflows/{workflow_id}/runs",
    response_model=RunResponse,
    status_code=201,
    operation_id="createRun",
    responses={
        200: {
            "model": RunResponse,
            "description": "Idempotent replay of the snapshot's canonical run.",
        },
        400: {"model": RunResponse},
        404: {"model": RunResponse},
        409: {"model": RunResponse},
        413: {
            "model": RunResponse,
            "description": "Request body too large.",
        },
        503: {"model": RunResponse},
    },
)
def create_run(
    workflow_id: str,
    request_body: RunCreateRequest,
    creation_service: ValidatedRunCreationService = Depends(
        get_validated_run_creation_service
    ),
) -> RunResponse | JSONResponse:
    """Create or replay one run using a server-owned validated snapshot."""
    try:
        creation = creation_service.create_run(
            workflow_id,
            request_body.snapshot_id,
            tags=request_body.tags,
        )
    except ValidatedSnapshotNotFoundError:
        return JSONResponse(
            status_code=404,
            content=RunResponse(
                ok=False,
                run=None,
                issues=[
                    _validated_snapshot_issue(
                        "VALIDATED_SNAPSHOT_NOT_FOUND",
                        "Validated input snapshot was not found.",
                    )
                ],
            ).model_dump(mode="json"),
        )
    except ValidatedSnapshotExpiredError:
        return JSONResponse(
            status_code=409,
            content=RunResponse(
                ok=False,
                run=None,
                issues=[
                    _validated_snapshot_issue(
                        "VALIDATED_SNAPSHOT_EXPIRED",
                        "Validated input snapshot expired. Validate again.",
                    )
                ],
            ).model_dump(mode="json"),
        )
    except ValidatedSnapshotStaleError:
        return JSONResponse(
            status_code=409,
            content=RunResponse(
                ok=False,
                run=None,
                issues=[
                    _validated_snapshot_issue(
                        "VALIDATED_SNAPSHOT_STALE",
                        "Workflow source changed after validation. Validate again.",
                    )
                ],
            ).model_dump(mode="json"),
        )
    except ValidatedSnapshotReplayConflictError:
        return JSONResponse(
            status_code=409,
            content=RunResponse(
                ok=False,
                run=None,
                issues=[
                    _validated_snapshot_issue(
                        "VALIDATED_SNAPSHOT_REPLAY_CONFLICT",
                        "Validated input snapshot is already linked to different run metadata.",
                    )
                ],
            ).model_dump(mode="json"),
        )
    except ValidatedSnapshotDataInvalidError:
        return JSONResponse(
            status_code=409,
            content=RunResponse(
                ok=False,
                run=None,
                issues=[
                    _validated_snapshot_issue(
                        "VALIDATED_SNAPSHOT_DATA_INVALID",
                        "Validated input snapshot cannot be used safely.",
                    )
                ],
            ).model_dump(mode="json"),
        )
    except ValidatedSnapshotBuildUnavailableError:
        return JSONResponse(
            status_code=503,
            content=RunResponse(
                ok=False,
                run=None,
                issues=[
                    _validated_snapshot_issue(
                        "VALIDATED_SNAPSHOT_BUILD_UNAVAILABLE",
                        "Workflow source identity could not be confirmed. Retry later.",
                    )
                ],
            ).model_dump(mode="json"),
        )

    response = RunResponse(
        ok=True,
        run=_run_record_response(creation.record),
        issues=[],
    )
    if not creation.created:
        return JSONResponse(status_code=200, content=response.model_dump(mode="json"))
    return response


@router.get("/runs/{run_id}", response_model=RunResponse, operation_id="getRun")
async def get_run(
    run_id: str,
    run_service: RunService = Depends(get_run_service),
) -> RunResponse | JSONResponse:
    """Return the run record."""
    try:
        record = run_service.get_run(run_id)
    except KeyError:
        return JSONResponse(
            status_code=404,
            content=RunResponse(
                ok=False,
                run=None,
                issues=[_run_not_found_issue(run_id)],
            ).model_dump(),
        )

    return RunResponse(ok=True, run=_run_record_response(record), issues=[])


@router.post(
    "/runs/{run_id}/start",
    response_model=RunResponse,
    status_code=202,
    operation_id="startRun",
    responses={
        404: {"model": RunResponse},
        409: {"model": RunResponse},
        503: {"model": RunResponse},
    },
)
def start_run(
    run_id: str,
    submission_service: RunSubmissionService = Depends(get_run_submission_service),
) -> RunResponse | JSONResponse:
    """Explicitly submit a planned run for durable worker execution."""
    try:
        record = submission_service.start_run(run_id)
    except KeyError:
        return JSONResponse(
            status_code=404,
            content=RunResponse(
                ok=False,
                run=None,
                issues=[_run_not_found_issue(run_id)],
            ).model_dump(mode="json"),
        )
    except RunNotReadyError as exc:
        return JSONResponse(
            status_code=409,
            content=RunResponse(
                ok=False,
                run=_run_record_response(exc.record),
                issues=[_run_not_ready_issue(exc.record.status.value)],
            ).model_dump(mode="json"),
        )
    except RunBuildIdentityMissingError as exc:
        return JSONResponse(
            status_code=409,
            content=RunResponse(
                ok=False,
                run=_run_record_response(exc.record),
                issues=[
                    _run_workflow_build_identity_missing_issue(exc.record.status.value)
                ],
            ).model_dump(mode="json"),
        )
    except RunStartConflictError as exc:
        return JSONResponse(
            status_code=409,
            content=RunResponse(
                ok=False,
                run=_run_record_response(exc.record),
                issues=[_run_start_conflict_issue(exc.record.status.value)],
            ).model_dump(mode="json"),
        )
    except RunSubmissionUnavailableError as exc:
        return JSONResponse(
            status_code=503,
            content=RunResponse(
                ok=False,
                run=_run_record_response(exc.record),
                issues=[_run_queue_unavailable_issue(exc.record.status.value)],
            ).model_dump(mode="json"),
        )

    return RunResponse(ok=True, run=_run_record_response(record), issues=[])


@router.post(
    "/runs/{run_id}/cancel",
    response_model=RunResponse,
    operation_id="cancelRun",
    responses={
        202: {"model": RunResponse},
        404: {"model": RunResponse},
        409: {"model": RunResponse},
        503: {"model": RunResponse},
    },
)
def cancel_run(
    run_id: str,
    cancellation_service: RunCancellationService = Depends(
        get_run_cancellation_service
    ),
) -> RunResponse | JSONResponse:
    """Cancel before execution or request an acknowledged RQ process stop."""
    try:
        result = cancellation_service.cancel_run(
            run_id,
            reason="User requested cancellation.",
        )
    except KeyError:
        return JSONResponse(
            status_code=404,
            content=RunResponse(
                ok=False,
                run=None,
                issues=[_run_not_found_issue(run_id)],
            ).model_dump(mode="json"),
        )
    except RunCancellationConflictError as exc:
        return JSONResponse(
            status_code=409,
            content=RunResponse(
                ok=False,
                run=_run_record_response(exc.record),
                issues=[_run_cancellation_conflict_issue(exc.record.status.value)],
            ).model_dump(mode="json"),
        )
    except RunCancellationUnavailableError as exc:
        return JSONResponse(
            status_code=503,
            content=RunResponse(
                ok=False,
                run=_run_record_response(exc.record),
                issues=[_run_cancellation_unavailable_issue(exc.record.status.value)],
            ).model_dump(mode="json"),
        )

    response = RunResponse(
        ok=True,
        run=_run_record_response(result.record),
        issues=[],
    )
    if result.stop_requested:
        return JSONResponse(status_code=202, content=response.model_dump(mode="json"))
    return response


@router.get(
    "/runs/{run_id}/events",
    response_model=RunEventsResponse,
    operation_id="listRunEvents",
)
async def list_run_events(
    run_id: str,
    after: str | None = None,
    limit: int = Query(default=50, ge=1),
    run_service: RunService = Depends(get_run_service),
) -> RunEventsResponse | JSONResponse:
    """List run events with cursor pagination."""
    try:
        run_service.get_run(run_id)
    except KeyError:
        return JSONResponse(
            status_code=404,
            content=RunEventsResponse(
                ok=False,
                run_id=run_id,
                events=[],
                next_cursor=None,
                issues=[_run_not_found_issue(run_id)],
            ).model_dump(),
        )

    try:
        events = run_service.list_events(run_id, after=after, limit=limit + 1)
    except KeyError:
        return JSONResponse(
            status_code=400,
            content=RunEventsResponse(
                ok=False,
                run_id=run_id,
                events=[],
                next_cursor=None,
                issues=[_run_cursor_not_found_issue()],
            ).model_dump(),
        )

    next_cursor: str | None = None
    if len(events) > limit:
        events = events[:limit]
        next_cursor = events[-1].event_id

    return RunEventsResponse(
        ok=True,
        run_id=run_id,
        events=[_run_event_response(event) for event in events],
        next_cursor=next_cursor,
        issues=[],
    )


@router.get(
    "/runs/{run_id}/logs", response_model=RunLogsResponse, operation_id="listRunLogs"
)
async def list_run_logs(
    run_id: str,
    stream_name: str = "stdout",
    after: str | None = None,
    limit: int = Query(default=50, ge=1),
    run_service: RunService = Depends(get_run_service),
) -> RunLogsResponse | JSONResponse:
    """List log chunks for a run stream with cursor pagination."""
    try:
        run_service.get_run(run_id)
    except KeyError:
        return JSONResponse(
            status_code=404,
            content=RunLogsResponse(
                ok=False,
                run_id=run_id,
                stream_name=stream_name,
                chunks=[],
                next_cursor=None,
                issues=[_run_not_found_issue(run_id)],
            ).model_dump(),
        )

    try:
        chunks = run_service.list_logs(
            run_id, stream_name, after=after, limit=limit + 1
        )
    except KeyError:
        return JSONResponse(
            status_code=400,
            content=RunLogsResponse(
                ok=False,
                run_id=run_id,
                stream_name=stream_name,
                chunks=[],
                next_cursor=None,
                issues=[_run_cursor_not_found_issue()],
            ).model_dump(),
        )

    next_cursor: str | None = None
    if len(chunks) > limit:
        chunks = chunks[:limit]
        next_cursor = chunks[-1].chunk_id

    return RunLogsResponse(
        ok=True,
        run_id=run_id,
        stream_name=stream_name,
        chunks=[_run_log_chunk_response(chunk) for chunk in chunks],
        next_cursor=next_cursor,
        issues=[],
    )
