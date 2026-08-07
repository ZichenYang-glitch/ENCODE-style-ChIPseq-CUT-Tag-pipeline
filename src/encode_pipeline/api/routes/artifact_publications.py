"""Read-only artifact-publication API routes."""

from __future__ import annotations

from datetime import datetime
from typing import Annotated, Any

from fastapi import APIRouter, Depends, Path, Query
from fastapi.responses import JSONResponse

from encode_pipeline.api.dependencies import get_artifact_publication_service
from encode_pipeline.api.models import (
    ArtifactPublicationDetailResponse,
    ArtifactPublicationIssueResponse,
    ArtifactPublicationListResponse,
    ArtifactPublicationResponse,
)
from encode_pipeline.services.artifact_publications import (
    ArtifactPublicationCursorInvalidError,
    ArtifactPublicationCursorNotFoundError,
    ArtifactPublicationDataInvalidError,
    ArtifactPublicationFilterInvalidError,
    ArtifactPublicationNotFoundError,
    ArtifactPublicationQueryService,
)


router = APIRouter(tags=["artifact-publications"])

RunIdPath = Annotated[
    str,
    Path(min_length=1, max_length=128, pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$"),
]
ArtifactIdPath = Annotated[
    str,
    Path(min_length=1, max_length=128, pattern=r"^[A-Za-z][A-Za-z0-9_.-]*$"),
]


def _issue(code: str, message: str, *, path: str) -> ArtifactPublicationIssueResponse:
    return ArtifactPublicationIssueResponse(
        code=code,
        message=message,
        severity="error",
        path=path,
        source="api",
        hint=None,
    )


def _list_content(response: ArtifactPublicationListResponse) -> dict[str, Any]:
    content = response.model_dump(mode="json", exclude_none=True)
    content["next_cursor"] = response.next_cursor
    return content


def _detail_content(response: ArtifactPublicationDetailResponse) -> dict[str, Any]:
    content = response.model_dump(mode="json", exclude_none=True)
    content["publication"] = (
        None
        if response.publication is None
        else response.publication.model_dump(mode="json")
    )
    return content


@router.get(
    "/artifact-publications",
    response_model=ArtifactPublicationListResponse,
    operation_id="listArtifactPublications",
    responses={
        400: {"model": ArtifactPublicationListResponse},
        500: {"model": ArtifactPublicationListResponse},
        "4XX": {"model": ArtifactPublicationListResponse},
    },
)
def list_artifact_publications(
    project_id: str | None = Query(
        default=None,
        min_length=36,
        max_length=36,
        pattern=r"^prj_[0-9a-f]{32}$",
    ),
    run_id: str | None = Query(
        default=None,
        min_length=1,
        max_length=128,
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    ),
    workflow_id: str | None = Query(
        default=None,
        min_length=1,
        max_length=255,
        pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    ),
    output_type: str | None = Query(
        default=None,
        min_length=1,
        max_length=128,
        pattern=r"^[A-Za-z][A-Za-z0-9_.-]*$",
    ),
    associated_run_sample_revision_id: str | None = Query(
        default=None,
        min_length=37,
        max_length=37,
        pattern=r"^smpr_[0-9a-f]{32}$",
    ),
    published_from: datetime | None = Query(
        default=None,
        description="Inclusive UTC publication timestamp lower bound.",
    ),
    published_before: datetime | None = Query(
        default=None,
        description="Exclusive UTC publication timestamp upper bound.",
    ),
    current_only: bool = Query(default=True),
    after: str | None = Query(
        default=None,
        min_length=16,
        max_length=1024,
        pattern=r"^artifactpubcur_[A-Za-z0-9_-]+$",
    ),
    limit: int = Query(default=50, ge=1, le=100),
    service: ArtifactPublicationQueryService = Depends(
        get_artifact_publication_service
    ),
) -> ArtifactPublicationListResponse | JSONResponse:
    """List indexed publications with relational filters and stable pagination."""
    try:
        page = service.list_artifact_publications(
            project_id=project_id,
            run_id=run_id,
            workflow_id=workflow_id,
            output_type=output_type,
            associated_run_sample_revision_id=(associated_run_sample_revision_id),
            published_from=published_from,
            published_before=published_before,
            current_only=current_only,
            after=after,
            limit=limit,
        )
        publications = [
            ArtifactPublicationResponse.from_summary(item)
            for item in page.artifact_publications
        ]
    except ArtifactPublicationFilterInvalidError:
        response = ArtifactPublicationListResponse(
            ok=False,
            publications=[],
            next_cursor=None,
            issues=[
                _issue(
                    "API_REQUEST_INVALID",
                    "Artifact publication query parameters are invalid.",
                    path="query",
                )
            ],
        )
        return JSONResponse(status_code=400, content=_list_content(response))
    except ArtifactPublicationCursorInvalidError:
        response = ArtifactPublicationListResponse(
            ok=False,
            publications=[],
            next_cursor=None,
            issues=[
                _issue(
                    "ARTIFACT_PUBLICATION_CURSOR_INVALID",
                    "Artifact publication cursor is invalid for this query.",
                    path="after",
                )
            ],
        )
        return JSONResponse(status_code=400, content=_list_content(response))
    except ArtifactPublicationCursorNotFoundError:
        response = ArtifactPublicationListResponse(
            ok=False,
            publications=[],
            next_cursor=None,
            issues=[
                _issue(
                    "ARTIFACT_PUBLICATION_CURSOR_NOT_FOUND",
                    "Artifact publication cursor is no longer available.",
                    path="after",
                )
            ],
        )
        return JSONResponse(status_code=400, content=_list_content(response))
    except (ArtifactPublicationDataInvalidError, ValueError):
        response = ArtifactPublicationListResponse(
            ok=False,
            publications=[],
            next_cursor=None,
            issues=[
                _issue(
                    "ARTIFACT_PUBLICATION_DATA_INVALID",
                    "Artifact publication data could not be read safely.",
                    path="publications",
                )
            ],
        )
        return JSONResponse(status_code=500, content=_list_content(response))

    return ArtifactPublicationListResponse(
        ok=True,
        publications=publications,
        next_cursor=page.next_cursor,
        issues=[],
    )


@router.get(
    "/artifact-publications/{run_id}/{artifact_id}",
    response_model=ArtifactPublicationDetailResponse,
    operation_id="getArtifactPublication",
    responses={
        400: {"model": ArtifactPublicationDetailResponse},
        404: {"model": ArtifactPublicationDetailResponse},
        500: {"model": ArtifactPublicationDetailResponse},
        "4XX": {"model": ArtifactPublicationDetailResponse},
    },
)
def get_artifact_publication(
    run_id: RunIdPath,
    artifact_id: ArtifactIdPath,
    generation: str = Query(
        min_length=76,
        max_length=76,
        pattern=r"^artifactgen-[0-9a-f]{64}$",
    ),
    service: ArtifactPublicationQueryService = Depends(
        get_artifact_publication_service
    ),
) -> ArtifactPublicationDetailResponse | JSONResponse:
    """Return one exact run/artifact/generation publication identity."""
    try:
        summary = service.get_artifact_publication(
            run_id=run_id,
            artifact_id=artifact_id,
            artifact_generation=generation,
        )
        publication = ArtifactPublicationResponse.from_summary(summary)
    except ArtifactPublicationFilterInvalidError:
        response = ArtifactPublicationDetailResponse(
            ok=False,
            publication=None,
            issues=[
                _issue(
                    "API_REQUEST_INVALID",
                    "Artifact publication identity is invalid.",
                    path="query",
                )
            ],
        )
        return JSONResponse(status_code=400, content=_detail_content(response))
    except ArtifactPublicationNotFoundError:
        response = ArtifactPublicationDetailResponse(
            ok=False,
            publication=None,
            issues=[
                _issue(
                    "ARTIFACT_PUBLICATION_NOT_FOUND",
                    "Artifact publication was not found.",
                    path="generation",
                )
            ],
        )
        return JSONResponse(status_code=404, content=_detail_content(response))
    except (ArtifactPublicationDataInvalidError, ValueError):
        response = ArtifactPublicationDetailResponse(
            ok=False,
            publication=None,
            issues=[
                _issue(
                    "ARTIFACT_PUBLICATION_DATA_INVALID",
                    "Artifact publication data could not be read safely.",
                    path="publication",
                )
            ],
        )
        return JSONResponse(status_code=500, content=_detail_content(response))

    return ArtifactPublicationDetailResponse(
        ok=True,
        publication=publication,
        issues=[],
    )
