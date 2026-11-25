"""Read sealed track versions for one experiment."""

from __future__ import annotations

from collections.abc import Iterator

from fastapi import APIRouter, Depends, Path, Query, Request
from pydantic import ValidationError
from fastapi.responses import JSONResponse, Response, StreamingResponse
from starlette.types import Receive, Scope, Send

from encode_pipeline.api.archived_tracks import (
    PRIVATE_NO_STORE,
    ArchivedTrackListResponse,
    ArchivedTrackResponse,
    archived_track_json,
)
from encode_pipeline.api.catalog_models import CatalogErrorResponse
from encode_pipeline.api.dependencies import AuthApiError, require_principal
from encode_pipeline.platform.authentication import AuthenticatedPrincipal
from encode_pipeline.services.archived_tracks import (
    ArchivedByteBody,
    ArchivedTrackError,
    ArchivedTrackService,
)


router = APIRouter(prefix="/catalog", tags=["catalog"])

_PROJECT_ID = Path(pattern=r"^prj_[0-9a-f]{32}$")
_EXPERIMENT_ID = Path(pattern=r"^exp_[0-9a-f]{32}$")
_ARCHIVE_VERSION_ID = Path(pattern=r"^aver_[0-9a-f]{32}$")
_JSON_ERRORS = {
    400: {"model": CatalogErrorResponse},
    401: {"model": CatalogErrorResponse},
    404: {"model": CatalogErrorResponse},
    409: {"model": CatalogErrorResponse},
    500: {"model": CatalogErrorResponse},
}
_BYTE_ERRORS = {
    **_JSON_ERRORS,
    416: {
        "description": "The requested byte range cannot be satisfied.",
        "headers": {
            "Content-Range": {
                "schema": {"type": "string"},
                "description": "Unsatisfied range, bytes */size.",
            },
            "Cache-Control": {
                "schema": {"type": "string"},
                "description": "private, no-store",
            },
        },
    },
}
_BYTE_SUCCESS = {
    "description": (
        "Bytes of the requested sealed member. One Range returns 206 with "
        "Content-Range. An unsatisfiable range returns 416. Responses use "
        "private, no-store caching."
    ),
    "content": {
        "application/octet-stream": {
            "schema": {"type": "string", "format": "binary"}
        }
    },
    "headers": {
        "Accept-Ranges": {"schema": {"type": "string"}},
        "Cache-Control": {"schema": {"type": "string"}},
        "Content-Range": {"schema": {"type": "string"}},
    },
}


class _ClosingStream(StreamingResponse):
    """Close the sealed reader when transmission ends or fails."""

    def __init__(self, content: Iterator[bytes], close_callback, **kwargs) -> None:
        super().__init__(content, **kwargs)
        self._close_callback = close_callback

    async def __call__(self, scope: Scope, receive: Receive, send: Send) -> None:
        try:
            await super().__call__(scope, receive, send)
        finally:
            self._close_callback()


def get_archived_track_service(request: Request) -> ArchivedTrackService:
    service = getattr(request.app.state, "archived_track_service", None)
    if isinstance(service, ArchivedTrackService):
        return service
    raise ArchivedTrackError(
        409,
        "archive_offline",
        "The archived track is unavailable.",
    )


def _dump_track(payload: dict) -> dict:
    return ArchivedTrackResponse.model_validate(payload).model_dump(
        mode="json",
        exclude_none=True,
    )


@router.get(
    "/projects/{project_id}/experiments/{experiment_id}/archived-tracks",
    response_model=ArchivedTrackListResponse,
    response_model_exclude_none=True,
    operation_id="catalogListArchivedTracks",
    responses=_JSON_ERRORS,
)
async def catalog_list_archived_tracks(
    project_id: str = _PROJECT_ID,
    experiment_id: str = _EXPERIMENT_ID,
    limit: int = Query(default=50, ge=1, le=100),
    offset: int = Query(default=0, ge=0),
    principal: AuthenticatedPrincipal = Depends(require_principal),
    service: ArchivedTrackService = Depends(get_archived_track_service),
) -> JSONResponse:
    """List sealed track versions for one experiment."""

    try:
        page = service.list_versions(
            principal,
            project_id,
            experiment_id,
            limit=limit,
            offset=offset,
        )
        content = {
            "items": [_dump_track(item) for item in page["items"]],
            "next_offset": page["next_offset"],
        }
        ArchivedTrackListResponse.model_validate(content)
    except ArchivedTrackError as exc:
        return _failure(exc)
    except ValidationError:
        return archived_track_json(
            400,
            "archive_invalid",
            "The archived track request is invalid.",
        )
    return JSONResponse(
        content=content,
        headers={"Cache-Control": PRIVATE_NO_STORE},
    )


@router.get(
    "/projects/{project_id}/experiments/{experiment_id}/archived-tracks/{archive_version_id}",
    response_model=ArchivedTrackResponse,
    response_model_exclude_none=True,
    operation_id="catalogGetArchivedTrack",
    responses=_JSON_ERRORS,
)
async def catalog_get_archived_track(
    project_id: str = _PROJECT_ID,
    experiment_id: str = _EXPERIMENT_ID,
    archive_version_id: str = _ARCHIVE_VERSION_ID,
    principal: AuthenticatedPrincipal = Depends(require_principal),
    service: ArchivedTrackService = Depends(get_archived_track_service),
) -> JSONResponse:
    """Return one sealed version. The shape matches a list item."""

    try:
        payload = service.get_version(
            principal,
            project_id,
            experiment_id,
            archive_version_id,
        )
    except ArchivedTrackError as exc:
        return _failure(exc)
    try:
        content = _dump_track(payload)
    except ValidationError:
        return archived_track_json(
            400,
            "archive_invalid",
            "The archived track request is invalid.",
        )
    return JSONResponse(
        content=content,
        headers={"Cache-Control": PRIVATE_NO_STORE},
    )


def _read_bytes(
    *,
    project_id: str,
    experiment_id: str,
    archive_version_id: str,
    request: Request,
    principal: AuthenticatedPrincipal,
    service: ArchivedTrackService,
    index: bool,
    head: bool,
) -> Response:
    try:
        body = service.open_bytes(
            principal,
            project_id,
            experiment_id,
            archive_version_id,
            index=index,
            range_header=request.headers.get("range"),
            head=head,
        )
    except ArchivedTrackError as exc:
        return _failure(exc)
    return _byte_response(body)


@router.get(
    "/projects/{project_id}/experiments/{experiment_id}/archived-tracks/{archive_version_id}/content",
    response_class=StreamingResponse,
    response_model=None,
    operation_id="catalogReadArchivedTrackContent",
    responses={200: _BYTE_SUCCESS, 206: _BYTE_SUCCESS, **_BYTE_ERRORS},
)
async def catalog_read_archived_track_content(
    request: Request,
    project_id: str = _PROJECT_ID,
    experiment_id: str = _EXPERIMENT_ID,
    archive_version_id: str = _ARCHIVE_VERSION_ID,
    principal: AuthenticatedPrincipal = Depends(require_principal),
    service: ArchivedTrackService = Depends(get_archived_track_service),
) -> Response:
    """Send one byte window of the sealed data member."""

    return _read_bytes(
        project_id=project_id,
        experiment_id=experiment_id,
        archive_version_id=archive_version_id,
        request=request,
        principal=principal,
        service=service,
        index=False,
        head=False,
    )


@router.head(
    "/projects/{project_id}/experiments/{experiment_id}/archived-tracks/{archive_version_id}/content",
    response_class=Response,
    response_model=None,
    operation_id="catalogHeadArchivedTrackContent",
    responses={200: _BYTE_SUCCESS, 206: _BYTE_SUCCESS, **_BYTE_ERRORS},
)
async def catalog_head_archived_track_content(
    request: Request,
    project_id: str = _PROJECT_ID,
    experiment_id: str = _EXPERIMENT_ID,
    archive_version_id: str = _ARCHIVE_VERSION_ID,
    principal: AuthenticatedPrincipal = Depends(require_principal),
    service: ArchivedTrackService = Depends(get_archived_track_service),
) -> Response:
    """Report the sealed data member's length and range without a body."""

    return _read_bytes(
        project_id=project_id,
        experiment_id=experiment_id,
        archive_version_id=archive_version_id,
        request=request,
        principal=principal,
        service=service,
        index=False,
        head=True,
    )


@router.get(
    "/projects/{project_id}/experiments/{experiment_id}/archived-tracks/{archive_version_id}/content.tbi",
    response_class=StreamingResponse,
    response_model=None,
    operation_id="catalogReadArchivedTrackIndex",
    responses={200: _BYTE_SUCCESS, 206: _BYTE_SUCCESS, **_BYTE_ERRORS},
)
async def catalog_read_archived_track_index(
    request: Request,
    project_id: str = _PROJECT_ID,
    experiment_id: str = _EXPERIMENT_ID,
    archive_version_id: str = _ARCHIVE_VERSION_ID,
    principal: AuthenticatedPrincipal = Depends(require_principal),
    service: ArchivedTrackService = Depends(get_archived_track_service),
) -> Response:
    """Send one byte window of the sealed index member."""

    return _read_bytes(
        project_id=project_id,
        experiment_id=experiment_id,
        archive_version_id=archive_version_id,
        request=request,
        principal=principal,
        service=service,
        index=True,
        head=False,
    )


@router.head(
    "/projects/{project_id}/experiments/{experiment_id}/archived-tracks/{archive_version_id}/content.tbi",
    response_class=Response,
    response_model=None,
    operation_id="catalogHeadArchivedTrackIndex",
    responses={200: _BYTE_SUCCESS, 206: _BYTE_SUCCESS, **_BYTE_ERRORS},
)
async def catalog_head_archived_track_index(
    request: Request,
    project_id: str = _PROJECT_ID,
    experiment_id: str = _EXPERIMENT_ID,
    archive_version_id: str = _ARCHIVE_VERSION_ID,
    principal: AuthenticatedPrincipal = Depends(require_principal),
    service: ArchivedTrackService = Depends(get_archived_track_service),
) -> Response:
    """Report the sealed index member's length and range without a body."""

    return _read_bytes(
        project_id=project_id,
        experiment_id=experiment_id,
        archive_version_id=archive_version_id,
        request=request,
        principal=principal,
        service=service,
        index=True,
        head=True,
    )


def _failure(exc: ArchivedTrackError) -> JSONResponse:
    if exc.status_code == 401:
        raise AuthApiError(exc.status_code, exc.code, exc.message)
    return archived_track_json(exc.status_code, exc.code, exc.message)


def _byte_response(body: ArchivedByteBody) -> Response:
    if body.chunks is None:
        return Response(
            status_code=body.status_code,
            headers=body.headers,
            media_type="application/octet-stream",
        )
    return _ClosingStream(
        body.chunks,
        body.close,
        status_code=body.status_code,
        headers=body.headers,
        media_type="application/octet-stream",
    )
