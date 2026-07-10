"""Preflight trigger route for workflow runs."""

from __future__ import annotations

import asyncio

from fastapi import APIRouter, BackgroundTasks, Depends
from fastapi.responses import JSONResponse

from encode_pipeline.api.dependencies import (
    get_preflight_service,
    get_run_service,
)
from encode_pipeline.api.models import RunResponse
from encode_pipeline.api.routes.runs import _run_not_found_issue, _run_record_response
from encode_pipeline.platform.runs import RunStatus
from encode_pipeline.services.preflight import LocalPreflightService
from encode_pipeline.services.runs import RunService


router = APIRouter(tags=["preflight"])


async def _run_preflight_in_background(
    preflight_service: LocalPreflightService,
    run_id: str,
) -> None:
    """Async wrapper for the synchronous preflight worker.

    Avoids passing a synchronous callable directly to ``BackgroundTasks``,
    which can hang under Python 3.13 with the current ASGI test client.
    """
    await asyncio.to_thread(preflight_service.run_preflight, run_id)


def _preflight_already_triggered_issue(current_status: str):
    from encode_pipeline.api.models import IssueResponse

    return IssueResponse(
        code="PREFLIGHT_ALREADY_TRIGGERED",
        message="Preflight has already been triggered for this run.",
        severity="error",
        path="run_id",
        source="api",
        context={"current_status": current_status},
    )


@router.post("/runs/{run_id}/preflight", response_model=RunResponse, status_code=202)
async def trigger_preflight(
    run_id: str,
    background_tasks: BackgroundTasks,
    preflight_service: LocalPreflightService = Depends(get_preflight_service),
    run_service: RunService = Depends(get_run_service),
) -> RunResponse | JSONResponse:
    """Trigger local preflight for an existing run.

    The actual planning/materialization/dry-run work runs in a background task
    so the HTTP response returns immediately with the run in ``VALIDATING``.
    """
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

    if record.status is not RunStatus.CREATED:
        return JSONResponse(
            status_code=409,
            content=RunResponse(
                ok=False,
                run=None,
                issues=[_preflight_already_triggered_issue(record.status.value)],
            ).model_dump(),
        )

    run_service.transition_run(
        run_id,
        RunStatus.VALIDATING,
        stage="preflight",
        message="Local preflight accepted.",
    )
    background_tasks.add_task(
        _run_preflight_in_background,
        preflight_service,
        run_id,
    )
    updated = run_service.get_run(run_id)
    return RunResponse(ok=True, run=_run_record_response(updated), issues=[])
