"""Workflow run lifecycle API routes."""

from __future__ import annotations

from typing import Any

from fastapi import APIRouter, Depends, Query
from fastapi.responses import JSONResponse

from encode_pipeline.api.dependencies import get_run_service
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
from encode_pipeline.platform.adapters import WorkflowInputs
from encode_pipeline.platform.results import Issue
from encode_pipeline.services.runs import RunService


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


def _api_request_invalid_issue(message: str, context: dict[str, Any] | None = None) -> IssueResponse:
    return _issue(
        code="API_REQUEST_INVALID",
        message=message,
        path="workflow_id",
        context=context or {},
    )


def _run_record_response(record: Any) -> RunRecordResponse:
    return RunRecordResponse(**record.to_dict())


def _run_event_response(event: Any) -> RunEventResponse:
    return RunEventResponse(**event.to_dict())


def _run_log_chunk_response(chunk: Any) -> RunLogChunkResponse:
    return RunLogChunkResponse(**chunk.to_dict())


@router.post("/workflows/{workflow_id}/runs", response_model=RunResponse, status_code=201)
def create_run(
    workflow_id: str,
    request_body: RunCreateRequest,
    run_service: RunService = Depends(get_run_service),
) -> RunResponse | JSONResponse:
    """Create a new run for the given workflow."""
    inputs = WorkflowInputs(
        config=request_body.config,
        samples=request_body.samples,
        options=request_body.options,
    )
    try:
        record = run_service.create_run(workflow_id, inputs, tags=request_body.tags)
    except KeyError:
        return JSONResponse(
            status_code=404,
            content=RunResponse(
                ok=False,
                run=None,
                issues=[_workflow_not_found_issue(workflow_id)],
            ).model_dump(),
        )
    except ValueError as exc:
        return JSONResponse(
            status_code=400,
            content=RunResponse(
                ok=False,
                run=None,
                issues=[_api_request_invalid_issue(str(exc), context={"workflow_id": workflow_id})],
            ).model_dump(),
        )

    return RunResponse(ok=True, run=_run_record_response(record), issues=[])


@router.get("/runs/{run_id}", response_model=RunResponse)
def get_run(
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


@router.post("/runs/{run_id}/cancel", response_model=RunResponse)
def cancel_run(
    run_id: str,
    run_service: RunService = Depends(get_run_service),
) -> RunResponse | JSONResponse:
    """Cancel an active run, or return an already-terminal run unchanged."""
    try:
        record = run_service.cancel_run(run_id, reason="User requested cancellation.")
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


@router.get("/runs/{run_id}/events", response_model=RunEventsResponse)
def list_run_events(
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


@router.get("/runs/{run_id}/logs", response_model=RunLogsResponse)
def list_run_logs(
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
        chunks = run_service.list_logs(run_id, stream_name, after=after, limit=limit + 1)
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
