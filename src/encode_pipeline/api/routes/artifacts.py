"""Read-only API routes for persisted run artifact references."""

from __future__ import annotations

from typing import Any

from fastapi import APIRouter, Depends, Query
from fastapi.responses import JSONResponse

from encode_pipeline.api.dependencies import get_run_service
from encode_pipeline.api.models import (
    ArtifactReferenceResponse,
    IssueResponse,
    RunArtifactDetailResponse,
    RunArtifactsResponse,
)
from encode_pipeline.services.runs import RunService


router = APIRouter(tags=["artifacts"])


def _issue(
    code: str,
    message: str,
    *,
    path: str,
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


def _run_not_found_issue(run_id: str) -> IssueResponse:
    return _issue(
        "RUN_NOT_FOUND",
        "Run was not found.",
        path="run_id",
        context={"run_id": run_id},
    )


def _artifact_not_found_issue(run_id: str, artifact_id: str) -> IssueResponse:
    return _issue(
        "RUN_ARTIFACT_NOT_FOUND",
        "Artifact was not found for this run.",
        path="artifact_id",
        context={"run_id": run_id, "artifact_id": artifact_id},
    )


def _artifact_cursor_not_found_issue() -> IssueResponse:
    return _issue(
        "RUN_ARTIFACT_CURSOR_NOT_FOUND",
        "Artifact cursor does not exist for this run.",
        path="after",
    )


def _artifact_data_invalid_issue() -> IssueResponse:
    return _issue(
        "RUN_ARTIFACT_DATA_INVALID",
        "Persisted artifact metadata is unavailable.",
        path="artifacts",
    )


@router.get(
    "/runs/{run_id}/artifacts",
    response_model=RunArtifactsResponse,
    response_model_exclude_none=True,
    operation_id="listRunArtifacts",
    responses={
        400: {"model": RunArtifactsResponse},
        "4XX": {"model": RunArtifactsResponse},
        404: {"model": RunArtifactsResponse},
        500: {"model": RunArtifactsResponse},
    },
)
async def list_run_artifacts(
    run_id: str,
    after: str | None = None,
    limit: int = Query(default=50, ge=1, le=100),
    run_service: RunService = Depends(get_run_service),
) -> RunArtifactsResponse | JSONResponse:
    """List persisted artifact references with a run-scoped keyset cursor."""
    try:
        run_service.get_run(run_id)
    except KeyError:
        return JSONResponse(
            status_code=404,
            content=RunArtifactsResponse(
                ok=False,
                run_id=run_id,
                artifacts=[],
                next_cursor=None,
                issues=[_run_not_found_issue(run_id)],
            ).model_dump(mode="json", exclude_none=True),
        )

    try:
        artifacts = run_service.list_artifacts(
            run_id,
            after=after,
            limit=limit + 1,
        )
    except KeyError:
        return JSONResponse(
            status_code=400,
            content=RunArtifactsResponse(
                ok=False,
                run_id=run_id,
                artifacts=[],
                next_cursor=None,
                issues=[_artifact_cursor_not_found_issue()],
            ).model_dump(mode="json", exclude_none=True),
        )
    except (AttributeError, TypeError, ValueError):
        return _invalid_list_response(run_id)

    next_cursor: str | None = None
    if len(artifacts) > limit:
        artifacts = artifacts[:limit]
        next_cursor = artifacts[-1].artifact_id
    try:
        projected = [
            ArtifactReferenceResponse.from_artifact(
                artifact,
                expected_run_id=run_id,
            )
            for artifact in artifacts
        ]
    except (AttributeError, TypeError, ValueError):
        return _invalid_list_response(run_id)

    return RunArtifactsResponse(
        ok=True,
        run_id=run_id,
        artifacts=projected,
        next_cursor=next_cursor,
        issues=[],
    )


@router.get(
    "/runs/{run_id}/artifacts/{artifact_id}",
    response_model=RunArtifactDetailResponse,
    response_model_exclude_none=True,
    operation_id="getRunArtifact",
    responses={
        "4XX": {"model": RunArtifactDetailResponse},
        404: {"model": RunArtifactDetailResponse},
        500: {"model": RunArtifactDetailResponse},
    },
)
async def get_run_artifact(
    run_id: str,
    artifact_id: str,
    run_service: RunService = Depends(get_run_service),
) -> RunArtifactDetailResponse | JSONResponse:
    """Return one persisted artifact reference scoped to ``run_id``."""
    try:
        run_service.get_run(run_id)
    except KeyError:
        return JSONResponse(
            status_code=404,
            content=RunArtifactDetailResponse(
                ok=False,
                run_id=run_id,
                artifact=None,
                issues=[_run_not_found_issue(run_id)],
            ).model_dump(mode="json"),
        )

    try:
        artifact = run_service.get_artifact(run_id, artifact_id)
    except KeyError:
        return JSONResponse(
            status_code=404,
            content=RunArtifactDetailResponse(
                ok=False,
                run_id=run_id,
                artifact=None,
                issues=[_artifact_not_found_issue(run_id, artifact_id)],
            ).model_dump(mode="json"),
        )
    except (AttributeError, TypeError, ValueError):
        return _invalid_detail_response(run_id)

    try:
        projected = ArtifactReferenceResponse.from_artifact(
            artifact,
            expected_run_id=run_id,
        )
    except (AttributeError, TypeError, ValueError):
        return _invalid_detail_response(run_id)

    return RunArtifactDetailResponse(
        ok=True,
        run_id=run_id,
        artifact=projected,
        issues=[],
    )


def _invalid_list_response(run_id: str) -> JSONResponse:
    return JSONResponse(
        status_code=500,
        content=RunArtifactsResponse(
            ok=False,
            run_id=run_id,
            artifacts=[],
            next_cursor=None,
            issues=[_artifact_data_invalid_issue()],
        ).model_dump(mode="json", exclude_none=True),
    )


def _invalid_detail_response(run_id: str) -> JSONResponse:
    return JSONResponse(
        status_code=500,
        content=RunArtifactDetailResponse(
            ok=False,
            run_id=run_id,
            artifact=None,
            issues=[_artifact_data_invalid_issue()],
        ).model_dump(mode="json"),
    )
