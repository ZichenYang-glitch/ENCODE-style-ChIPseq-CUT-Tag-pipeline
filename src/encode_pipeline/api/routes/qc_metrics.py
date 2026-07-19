"""Read-only API route for persisted run QC metrics."""

from __future__ import annotations

from typing import Any

from fastapi import APIRouter, Depends, Query
from fastapi.responses import JSONResponse

from encode_pipeline.api.dependencies import get_run_service
from encode_pipeline.api.models import (
    IssueResponse,
    QcMetricResponse,
    RunQcMetricsResponse,
)
from encode_pipeline.platform.result_generations import (
    QcMetricCursor,
    decode_qc_metric_cursor,
    encode_qc_metric_cursor,
)
from encode_pipeline.services.run_repositories import ResultGenerationChangedError
from encode_pipeline.services.runs import RunService


router = APIRouter(tags=["qc-metrics"])


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


def _cursor_not_found_issue() -> IssueResponse:
    return _issue(
        "RUN_QC_METRIC_CURSOR_NOT_FOUND",
        "QC metric cursor does not exist for this run.",
        path="after",
    )


def _data_invalid_issue() -> IssueResponse:
    return _issue(
        "RUN_QC_METRIC_DATA_INVALID",
        "Persisted QC metric data is unavailable.",
        path="qc_metrics",
    )


def _generation_changed_issue() -> IssueResponse:
    return _issue(
        "RUN_QC_METRIC_GENERATION_CHANGED",
        "QC metric generation changed. Reload from the first page.",
        path="generation",
    )


@router.get(
    "/runs/{run_id}/qc-metrics",
    response_model=RunQcMetricsResponse,
    operation_id="listRunQcMetrics",
    responses={
        400: {"model": RunQcMetricsResponse},
        "4XX": {"model": RunQcMetricsResponse},
        404: {"model": RunQcMetricsResponse},
        409: {"model": RunQcMetricsResponse},
        500: {"model": RunQcMetricsResponse},
    },
)
def list_run_qc_metrics(
    run_id: str,
    after: str | None = Query(
        default=None,
        pattern=r"^qccur_[A-Za-z0-9_-]{1,1018}$",
        max_length=1024,
    ),
    generation: str | None = Query(
        default=None,
        pattern=r"^qcgen-[0-9a-f]{64}$",
        max_length=70,
    ),
    limit: int = Query(default=50, ge=1, le=100),
    run_service: RunService = Depends(get_run_service),
) -> RunQcMetricsResponse | JSONResponse:
    """List persisted QC metrics with a run-scoped keyset cursor."""
    try:
        after_metric_id = None
        if after is not None:
            if generation is None:
                raise KeyError(after)
            try:
                after_metric_id = decode_qc_metric_cursor(
                    after,
                    run_id=run_id,
                    qc_generation=generation,
                ).after_metric_id
            except ValueError as exc:
                raise KeyError(after) from exc
        qc_generation, metrics = _load_metric_page(
            run_service,
            run_id,
            generation,
            after_metric_id,
            limit + 1,
        )
    except ResultGenerationChangedError:
        return JSONResponse(
            status_code=409,
            content=_error_content(
                RunQcMetricsResponse(
                    ok=False,
                    run_id=run_id,
                    qc_generation=None,
                    qc_metrics=[],
                    next_cursor=None,
                    issues=[_generation_changed_issue()],
                )
            ),
        )
    except KeyError:
        try:
            run_service.get_run(run_id)
        except KeyError:
            return JSONResponse(
                status_code=404,
                content=_error_content(
                    RunQcMetricsResponse(
                        ok=False,
                        run_id=run_id,
                        qc_generation=None,
                        qc_metrics=[],
                        next_cursor=None,
                        issues=[_run_not_found_issue(run_id)],
                    )
                ),
            )
        return JSONResponse(
            status_code=400,
            content=_error_content(
                RunQcMetricsResponse(
                    ok=False,
                    run_id=run_id,
                    qc_generation=None,
                    qc_metrics=[],
                    next_cursor=None,
                    issues=[_cursor_not_found_issue()],
                )
            ),
        )
    except (AttributeError, TypeError, ValueError):
        return _invalid_response(run_id)

    next_cursor: str | None = None
    if len(metrics) > limit:
        metrics = metrics[:limit]
        if qc_generation is None:
            return _invalid_response(run_id)
        next_cursor = encode_qc_metric_cursor(
            QcMetricCursor(
                run_id=run_id,
                qc_generation=qc_generation,
                after_metric_id=metrics[-1].metric_id,
            )
        )
    try:
        projected = [
            QcMetricResponse.from_metric(metric, expected_run_id=run_id)
            for metric in metrics
        ]
    except (AttributeError, TypeError, ValueError):
        return _invalid_response(run_id)

    return RunQcMetricsResponse(
        ok=True,
        run_id=run_id,
        qc_generation=qc_generation,
        qc_metrics=projected,
        next_cursor=next_cursor,
        issues=[],
    )


def _load_metric_page(
    run_service: RunService,
    run_id: str,
    expected_generation: str | None,
    after_metric_id: str | None,
    limit: int,
) -> tuple[str | None, tuple[Any, ...]]:
    run_service.get_run(run_id)
    return run_service.list_qc_metrics_page(
        run_id,
        expected_generation=expected_generation,
        after=after_metric_id,
        limit=limit,
    )


def _invalid_response(run_id: str) -> JSONResponse:
    return JSONResponse(
        status_code=500,
        content=_error_content(
            RunQcMetricsResponse(
                ok=False,
                run_id=run_id,
                qc_generation=None,
                qc_metrics=[],
                next_cursor=None,
                issues=[_data_invalid_issue()],
            )
        ),
    )


def _error_content(response: RunQcMetricsResponse) -> dict[str, Any]:
    content = response.model_dump(mode="json", exclude_none=True)
    # qc_generation is required-but-nullable in the versioned response model.
    # Preserve that wire field while continuing to omit unrelated optional nulls.
    content["qc_generation"] = None
    return content
