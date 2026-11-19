"""Read-only published-result context for a member and the result page."""

from __future__ import annotations

from fastapi import APIRouter, Depends, Query
from fastapi.responses import JSONResponse

from encode_pipeline.api.catalog_models import CatalogErrorResponse
from encode_pipeline.api.dependencies import (
    get_reference_profile_service,
    get_run_service,
)
from encode_pipeline.api.models import (
    PublishedResultContextErrorResponse,
    PublishedResultContextResponse,
    QcMetricResponse,
    ReferenceProfileRevisionResponse,
)
from encode_pipeline.api.routes.catalog import get_catalog_service
from encode_pipeline.services.experiment_catalog import ExperimentCatalogService
from encode_pipeline.services.published_result_context import (
    PublishedResultContextReader,
    PublishedResultEvidenceError,
    RecordedPreflightFailure,
)
from encode_pipeline.services.reference_profiles import ReferenceProfileService
from encode_pipeline.services.runs import RunService


router = APIRouter(prefix="/agent", tags=["result-context"])


@router.get(
    "/projects/{project_id}/artifact-publications/{run_id}/{artifact_id}/context",
    response_model=PublishedResultContextResponse,
    operation_id="readPublishedResultContext",
    responses={
        404: {"model": CatalogErrorResponse},
        500: {"model": PublishedResultContextErrorResponse},
    },
)
def read_published_result_context(
    project_id: str,
    run_id: str,
    artifact_id: str,
    generation: str = Query(
        min_length=76,
        max_length=76,
        pattern=r"^artifactgen-[0-9a-f]{64}$",
    ),
    catalog: ExperimentCatalogService = Depends(get_catalog_service),
    runs: RunService = Depends(get_run_service),
    reference_profiles: ReferenceProfileService = Depends(
        get_reference_profile_service
    ),
) -> JSONResponse:
    """Read one published result with the caller's existing project visibility.

    This is a GET. It does not trigger preflight, create a run, or start one.
    Missing reference, QC, and preflight evidence stay explicit.
    """

    reader = PublishedResultContextReader(catalog, runs, reference_profiles)
    try:
        context = reader.read(
            project_id=project_id,
            run_id=run_id,
            artifact_id=artifact_id,
            artifact_generation=generation,
        )
        reference = (
            None
            if context.reference is None
            else ReferenceProfileRevisionResponse.from_summary(context.reference)
        )
        metrics = [
            QcMetricResponse.from_metric(metric, expected_run_id=run_id)
            for metric in context.qc_metrics
        ]
        body = PublishedResultContextResponse(
            ok=True,
            project_id=project_id,
            run_id=run_id,
            artifact_id=artifact_id,
            artifact_generation=generation,
            experiment_trace=context.trace,
            run_status=context.run_status,
            current_stage=context.current_stage,
            reference_evidence=context.reference_evidence,
            reference_profile=reference,
            qc_evidence=context.qc_evidence,
            qc_generation=context.qc_generation,
            qc_metrics=metrics,
            qc_next_cursor=context.qc_next_cursor,
            preflight_failure_evidence=context.preflight_failure_evidence,
            preflight_failures=[
                _preflight_payload(failure) for failure in context.preflight_failures
            ],
        )
    except PublishedResultEvidenceError as exc:
        return _error(exc.code, exc.message)
    except (AttributeError, TypeError, ValueError):
        return _error(
            "RUN_RESULT_CONTEXT_INCONSISTENT",
            "Published result context could not be read.",
        )
    content = body.model_dump(mode="json", exclude_none=True)
    if body.qc_evidence == "absent":
        content.pop("qc_metrics", None)
    else:
        # QC metrics keep explicit null coordinates; exclude_none would drop them.
        content["qc_metrics"] = [
            metric.model_dump(mode="json") for metric in body.qc_metrics
        ]
    if body.preflight_failure_evidence == "absent":
        content.pop("preflight_failures", None)
    return JSONResponse(content=content)


def _preflight_payload(failure: RecordedPreflightFailure) -> dict[str, str]:
    payload = {
        "event_id": failure.event_id,
        "event_type": failure.event_type,
        "message": failure.message,
    }
    if failure.stage is not None:
        payload["stage"] = failure.stage
    if failure.code is not None:
        payload["code"] = failure.code
    return payload


def _error(code: str, message: str) -> JSONResponse:
    payload = PublishedResultContextErrorResponse(
        ok=False,
        code=code,
        message=message,
    )
    return JSONResponse(status_code=500, content=payload.model_dump(mode="json"))
