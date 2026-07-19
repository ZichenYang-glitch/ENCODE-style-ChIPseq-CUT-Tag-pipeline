"""Read-only API routes for persisted run artifact references."""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

from fastapi import APIRouter, Depends, Query
from fastapi.responses import JSONResponse, StreamingResponse
from starlette.background import BackgroundTask
from starlette.types import Receive, Scope, Send

from encode_pipeline.api.dependencies import (
    get_artifact_download_service,
    get_run_service,
)
from encode_pipeline.api.models import (
    ArtifactDownloadIssueContextResponse,
    ArtifactDownloadIssueResponse,
    ArtifactReferenceResponse,
    IssueResponse,
    RunArtifactDetailResponse,
    RunArtifactDownloadErrorResponse,
    RunArtifactsResponse,
)
from encode_pipeline.platform.result_generations import (
    ArtifactCursor,
    decode_artifact_cursor,
    encode_artifact_cursor,
)
from encode_pipeline.services.artifact_downloads import ArtifactDownloadService
from encode_pipeline.services.run_repositories import ResultGenerationChangedError
from encode_pipeline.services.runs import RunService


router = APIRouter(tags=["artifacts"])


_DOWNLOAD_STATUS_BY_CODE = {
    "RUN_NOT_FOUND": 404,
    "RUN_ARTIFACT_NOT_FOUND": 404,
    "RUN_ARTIFACT_DOWNLOAD_CONFLICT": 409,
    "RUN_ARTIFACT_DOWNLOAD_DATA_INVALID": 500,
    "RUN_ARTIFACT_DOWNLOAD_UNAVAILABLE": 500,
}


class _ArtifactStreamingResponse(StreamingResponse):
    """Close owned descriptors even when response transmission raises."""

    def __init__(
        self,
        *args: Any,
        close_callback: Callable[[], None],
        **kwargs: Any,
    ) -> None:
        super().__init__(*args, **kwargs)
        self._close_callback = close_callback

    async def __call__(self, scope: Scope, receive: Receive, send: Send) -> None:
        try:
            await super().__call__(scope, receive, send)
        finally:
            self._close_callback()


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


def _artifact_generation_changed_issue() -> IssueResponse:
    return _issue(
        "RUN_ARTIFACT_GENERATION_CHANGED",
        "Artifact generation changed. Reload from the first page.",
        path="generation",
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
        409: {"model": RunArtifactsResponse},
        500: {"model": RunArtifactsResponse},
    },
)
async def list_run_artifacts(
    run_id: str,
    after: str | None = Query(
        default=None,
        pattern=r"^artifactcur_[A-Za-z0-9_-]{1,1012}$",
        max_length=1024,
    ),
    generation: str | None = Query(
        default=None,
        pattern=r"^artifactgen-[0-9a-f]{64}$",
        max_length=76,
    ),
    limit: int = Query(default=50, ge=1, le=100),
    run_service: RunService = Depends(get_run_service),
) -> RunArtifactsResponse | JSONResponse:
    """List persisted artifact references with a run-scoped keyset cursor."""
    try:
        after_artifact_id = None
        if after is not None:
            if generation is None:
                raise KeyError(after)
            try:
                after_artifact_id = decode_artifact_cursor(
                    after,
                    run_id=run_id,
                    artifact_generation=generation,
                ).after_artifact_id
            except ValueError as exc:
                raise KeyError(after) from exc
        run_service.get_run(run_id)
        artifact_generation, artifacts = run_service.list_artifacts_page(
            run_id,
            expected_generation=generation,
            after=after_artifact_id,
            limit=limit + 1,
        )
    except ResultGenerationChangedError:
        return _artifact_generation_changed_list_response(run_id)
    except KeyError:
        try:
            run_service.get_run(run_id)
        except KeyError:
            return JSONResponse(
                status_code=404,
                content=_artifact_list_error_content(
                    RunArtifactsResponse(
                        ok=False,
                        run_id=run_id,
                        artifact_generation=None,
                        artifacts=[],
                        next_cursor=None,
                        issues=[_run_not_found_issue(run_id)],
                    )
                ),
            )
        return JSONResponse(
            status_code=400,
            content=_artifact_list_error_content(
                RunArtifactsResponse(
                    ok=False,
                    run_id=run_id,
                    artifact_generation=None,
                    artifacts=[],
                    next_cursor=None,
                    issues=[_artifact_cursor_not_found_issue()],
                )
            ),
        )
    except (AttributeError, TypeError, ValueError):
        return _invalid_list_response(run_id)

    next_cursor: str | None = None
    if len(artifacts) > limit:
        artifacts = artifacts[:limit]
        if artifact_generation is None:
            return _invalid_list_response(run_id)
        next_cursor = encode_artifact_cursor(
            ArtifactCursor(
                run_id=run_id,
                artifact_generation=artifact_generation,
                after_artifact_id=artifacts[-1].artifact_id,
            )
        )
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

    return JSONResponse(
        content=_artifact_list_content(
            RunArtifactsResponse(
                ok=True,
                run_id=run_id,
                artifact_generation=artifact_generation,
                artifacts=projected,
                next_cursor=next_cursor,
                issues=[],
            )
        )
    )


@router.get(
    "/runs/{run_id}/artifacts/{artifact_id}",
    response_model=RunArtifactDetailResponse,
    response_model_exclude_none=True,
    operation_id="getRunArtifact",
    responses={
        "4XX": {"model": RunArtifactDetailResponse},
        404: {"model": RunArtifactDetailResponse},
        409: {"model": RunArtifactDetailResponse},
        500: {"model": RunArtifactDetailResponse},
    },
)
async def get_run_artifact(
    run_id: str,
    artifact_id: str,
    generation: str = Query(
        pattern=r"^artifactgen-[0-9a-f]{64}$",
        max_length=76,
    ),
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
                artifact_generation=None,
                artifact=None,
                issues=[_run_not_found_issue(run_id)],
            ).model_dump(mode="json"),
        )

    try:
        artifact_generation, artifact = run_service.get_artifact_at_generation(
            run_id,
            artifact_id,
            expected_generation=generation,
        )
    except ResultGenerationChangedError:
        return JSONResponse(
            status_code=409,
            content=_artifact_detail_error_content(
                RunArtifactDetailResponse(
                    ok=False,
                    run_id=run_id,
                    artifact_generation=None,
                    artifact=None,
                    issues=[_artifact_generation_changed_issue()],
                )
            ),
        )
    except KeyError:
        return JSONResponse(
            status_code=404,
            content=RunArtifactDetailResponse(
                ok=False,
                run_id=run_id,
                artifact_generation=None,
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

    return JSONResponse(
        content=_artifact_detail_content(
            RunArtifactDetailResponse(
                ok=True,
                run_id=run_id,
                artifact_generation=artifact_generation,
                artifact=projected,
                issues=[],
            )
        )
    )


@router.get(
    "/runs/{run_id}/artifacts/{artifact_id}/download",
    response_class=StreamingResponse,
    response_model=None,
    operation_id="downloadRunArtifact",
    responses={
        200: {
            "description": "The persisted artifact byte stream.",
            "content": {
                "application/octet-stream": {
                    "schema": {"type": "string", "format": "binary"}
                }
            },
        },
        "4XX": {"model": RunArtifactDownloadErrorResponse},
        404: {"model": RunArtifactDownloadErrorResponse},
        409: {"model": RunArtifactDownloadErrorResponse},
        500: {"model": RunArtifactDownloadErrorResponse},
    },
)
def download_run_artifact(
    run_id: str,
    artifact_id: str,
    generation: str = Query(
        pattern=r"^artifactgen-[0-9a-f]{64}$",
        max_length=76,
    ),
    revision: str = Query(
        pattern=r"^artifactrev-[0-9a-f]{64}$",
        max_length=76,
    ),
    download_service: ArtifactDownloadService = Depends(get_artifact_download_service),
) -> StreamingResponse | JSONResponse:
    """Stream one persisted run-scoped artifact through a safe descriptor."""
    result = download_service.prepare(
        run_id,
        artifact_id,
        expected_generation=generation,
        expected_revision=revision,
    )
    if result.is_failure or result.value is None:
        issues = [
            ArtifactDownloadIssueResponse(
                code=issue.code,
                message=issue.message,
                severity=issue.severity.value,
                path=issue.path,
                source=issue.source,
                hint=issue.hint,
                context=ArtifactDownloadIssueContextResponse.from_issue_context(
                    issue.context
                ),
            )
            for issue in result.issues
        ]
        code = issues[0].code if issues else "RUN_ARTIFACT_DOWNLOAD_UNAVAILABLE"
        body = RunArtifactDownloadErrorResponse(
            run_id=run_id,
            artifact_id=artifact_id,
            issues=issues,
        )
        return JSONResponse(
            status_code=_DOWNLOAD_STATUS_BY_CODE.get(code, 500),
            content=body.model_dump(mode="json", exclude_none=True),
        )

    plan = result.value
    try:
        return _ArtifactStreamingResponse(
            plan.iter_bytes(),
            close_callback=plan.close,
            media_type=plan.media_type,
            headers={
                "Content-Disposition": plan.content_disposition,
                "Content-Length": str(plan.size_bytes),
                "Cache-Control": "private, no-store",
                "X-Content-Type-Options": "nosniff",
                "Accept-Ranges": "none",
            },
            background=BackgroundTask(plan.close),
        )
    except BaseException:
        plan.close()
        raise


def _invalid_list_response(run_id: str) -> JSONResponse:
    return JSONResponse(
        status_code=500,
        content=_artifact_list_error_content(
            RunArtifactsResponse(
                ok=False,
                run_id=run_id,
                artifact_generation=None,
                artifacts=[],
                next_cursor=None,
                issues=[_artifact_data_invalid_issue()],
            )
        ),
    )


def _invalid_detail_response(run_id: str) -> JSONResponse:
    return JSONResponse(
        status_code=500,
        content=RunArtifactDetailResponse(
            ok=False,
            run_id=run_id,
            artifact_generation=None,
            artifact=None,
            issues=[_artifact_data_invalid_issue()],
        ).model_dump(mode="json"),
    )


def _artifact_generation_changed_list_response(run_id: str) -> JSONResponse:
    return JSONResponse(
        status_code=409,
        content=_artifact_list_error_content(
            RunArtifactsResponse(
                ok=False,
                run_id=run_id,
                artifact_generation=None,
                artifacts=[],
                next_cursor=None,
                issues=[_artifact_generation_changed_issue()],
            )
        ),
    )


def _artifact_list_error_content(response: RunArtifactsResponse) -> dict[str, Any]:
    return _artifact_list_content(response)


def _artifact_list_content(response: RunArtifactsResponse) -> dict[str, Any]:
    content = response.model_dump(mode="json", exclude_none=True)
    content["artifact_generation"] = response.artifact_generation
    return content


def _artifact_detail_error_content(
    response: RunArtifactDetailResponse,
) -> dict[str, Any]:
    return _artifact_detail_content(response)


def _artifact_detail_content(response: RunArtifactDetailResponse) -> dict[str, Any]:
    content = response.model_dump(mode="json", exclude_none=True)
    content["artifact_generation"] = response.artifact_generation
    return content
