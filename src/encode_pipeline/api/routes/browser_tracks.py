"""Read-only WashU track routes. They do not change artifact download."""

from __future__ import annotations

from collections.abc import Iterator

from fastapi import APIRouter, Body, Depends, Query, Request
from fastapi.responses import JSONResponse, Response, StreamingResponse
from starlette.types import Receive, Scope, Send

from encode_pipeline.api.browser_models import (
    BrowserErrorResponse,
    BrowserPrepareRequest,
    BrowserPrepareResponse,
    BrowserRegionSummaryResponse,
    BrowserSessionResponse,
    BrowserTrackResponse,
    browser_error,
)
from encode_pipeline.services.browser_tracks import (
    BrowserByteBody,
    BrowserTrackError,
    BrowserTrackService,
)


router = APIRouter(tags=["browser"])


class _ClosingStream(StreamingResponse):
    """Close owned descriptors even when transmission fails."""

    def __init__(self, content: Iterator[bytes], close_callback, **kwargs) -> None:
        super().__init__(content, **kwargs)
        self._close_callback = close_callback

    async def __call__(self, scope: Scope, receive: Receive, send: Send) -> None:
        try:
            await super().__call__(scope, receive, send)
        finally:
            self._close_callback()


def get_browser_track_service(request: Request) -> BrowserTrackService:
    service = getattr(request.app.state, "browser_track_service", None)
    if isinstance(service, BrowserTrackService):
        return service
    service = BrowserTrackService(
        run_service=request.app.state.run_service,
        download_service=request.app.state.artifact_download_service,
        workspace_root=request.app.state.workspace_root,
        reference_profile_service=request.app.state.reference_profile_service,
    )
    request.app.state.browser_track_service = service
    return service


@router.get(
    "/runs/{run_id}/browser-session",
    response_model=BrowserSessionResponse,
    response_model_exclude_none=True,
    operation_id="getRunBrowserSession",
    responses={
        404: {"model": BrowserErrorResponse},
        500: {"model": BrowserErrorResponse},
    },
)
def get_run_browser_session(
    run_id: str,
    service: BrowserTrackService = Depends(get_browser_track_service),
) -> BrowserSessionResponse | JSONResponse:
    try:
        payload = service.session(run_id)
    except BrowserTrackError as exc:
        return _error(run_id, exc)
    except Exception:
        return _unexpected(run_id)
    return BrowserSessionResponse.model_validate(payload)


@router.post(
    "/runs/{run_id}/browser-tracks/{artifact_id}/derivative",
    response_model=BrowserPrepareResponse,
    response_model_exclude_none=True,
    operation_id="prepareRunBrowserTrack",
    responses={
        404: {"model": BrowserErrorResponse},
        409: {"model": BrowserErrorResponse},
        500: {"model": BrowserErrorResponse},
    },
)
def prepare_run_browser_track(
    run_id: str,
    artifact_id: str,
    body: BrowserPrepareRequest = Body(),
    service: BrowserTrackService = Depends(get_browser_track_service),
) -> BrowserPrepareResponse | JSONResponse:
    try:
        track = service.prepare(
            run_id,
            artifact_id,
            generation=body.generation,
            revision=body.revision,
        )
    except BrowserTrackError as exc:
        return _error(run_id, exc)
    except Exception:
        return _unexpected(run_id)
    return BrowserPrepareResponse(
        ok=True,
        run_id=run_id,
        track=BrowserTrackResponse.model_validate(track),
    )


@router.get(
    "/runs/{run_id}/browser-tracks/{artifact_id}/region-summary/{generation}/{revision}",
    response_model=BrowserRegionSummaryResponse,
    operation_id="summarizeRunBrowserTrackRegion",
    responses={
        400: {"model": BrowserErrorResponse},
        404: {"model": BrowserErrorResponse},
        409: {"model": BrowserErrorResponse},
        500: {"model": BrowserErrorResponse},
    },
)
def summarize_run_browser_track_region(
    run_id: str,
    artifact_id: str,
    generation: str,
    revision: str,
    chrom: str = Query(min_length=1, max_length=128, pattern=r"^[A-Za-z0-9_.:-]+$"),
    start: int = Query(ge=0),
    end: int = Query(ge=1),
    service: BrowserTrackService = Depends(get_browser_track_service),
) -> BrowserRegionSummaryResponse | JSONResponse:
    try:
        payload = service.region_summary(
            run_id,
            artifact_id,
            generation=generation,
            revision=revision,
            chrom=chrom,
            start=start,
            end=end,
        )
    except BrowserTrackError as exc:
        return _error(run_id, exc)
    except Exception:
        return _unexpected(run_id)
    return BrowserRegionSummaryResponse.model_validate(payload)


_BYTE_RESPONSES = {
    200: {
        "description": (
            "Artifact or display-derivative bytes. A Range header returns 206 "
            "with Content-Range; an invalid range returns 416."
        ),
        "content": {
            "application/octet-stream": {
                "schema": {"type": "string", "format": "binary"}
            }
        },
    },
    404: {"model": BrowserErrorResponse},
    409: {"model": BrowserErrorResponse},
    500: {"model": BrowserErrorResponse},
}


def _read_track_bytes(
    *,
    run_id: str,
    artifact_id: str,
    generation: str,
    file_token: str,
    request: Request,
    service: BrowserTrackService,
    head: bool,
) -> Response:
    index = file_token.endswith(".tbi")
    revision = file_token[: -len(".tbi")] if index else file_token
    try:
        body = service.open_bytes(
            run_id,
            artifact_id,
            generation=generation,
            revision=revision,
            index=index,
            range_header=request.headers.get("range"),
            head=head,
        )
    except BrowserTrackError as exc:
        return _error(run_id, exc)
    except Exception:
        return _unexpected(run_id)
    return _byte_response(body)


@router.get(
    "/runs/{run_id}/browser-tracks/{artifact_id}/bytes/{generation}/{file_token}",
    response_class=StreamingResponse,
    response_model=None,
    operation_id="readRunBrowserTrackBytes",
    responses=_BYTE_RESPONSES,
)
def read_run_browser_track_bytes(
    run_id: str,
    artifact_id: str,
    generation: str,
    file_token: str,
    request: Request,
    service: BrowserTrackService = Depends(get_browser_track_service),
) -> Response:
    return _read_track_bytes(
        run_id=run_id,
        artifact_id=artifact_id,
        generation=generation,
        file_token=file_token,
        request=request,
        service=service,
        head=False,
    )


@router.head(
    "/runs/{run_id}/browser-tracks/{artifact_id}/bytes/{generation}/{file_token}",
    response_class=Response,
    response_model=None,
    operation_id="headRunBrowserTrackBytes",
    responses=_BYTE_RESPONSES,
)
def head_run_browser_track_bytes(
    run_id: str,
    artifact_id: str,
    generation: str,
    file_token: str,
    request: Request,
    service: BrowserTrackService = Depends(get_browser_track_service),
) -> Response:
    return _read_track_bytes(
        run_id=run_id,
        artifact_id=artifact_id,
        generation=generation,
        file_token=file_token,
        request=request,
        service=service,
        head=True,
    )


def _unexpected(run_id: str) -> JSONResponse:
    """Redacted 500 for an unexpected failure. The exception text stays off the wire."""

    return JSONResponse(
        status_code=500,
        content=browser_error(
            run_id=run_id,
            code="INTERNAL_SERVER_ERROR",
            message="The browser view is temporarily unavailable.",
            path="run_id",
            reason=None,
            source="runtime",
        ),
    )


def _error(run_id: str, exc: BrowserTrackError) -> JSONResponse:
    return JSONResponse(
        status_code=exc.status,
        content=browser_error(
            run_id=run_id,
            code=exc.code,
            message=exc.message,
            path=exc.path,
            reason=exc.reason,
        ),
    )


def _byte_response(body: BrowserByteBody) -> Response:
    if body.chunks is None:
        return Response(
            status_code=body.status,
            headers=body.headers,
            media_type=body.media_type,
        )
    return _ClosingStream(
        body.chunks,
        body.close,
        status_code=body.status,
        headers=body.headers,
        media_type=body.media_type,
    )
