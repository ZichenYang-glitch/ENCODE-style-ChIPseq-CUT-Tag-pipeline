"""HTTP routes for registered results that are not platform runs."""

from __future__ import annotations

import os
from contextlib import suppress

from fastapi import APIRouter, Depends, Header, Path, Query, Request
from fastapi.responses import JSONResponse, StreamingResponse

from encode_pipeline.api.dependencies import (
    AuthApiError,
    require_administrator,
    require_principal,
)
from encode_pipeline.api.models import IssueResponse, ValidationResponse
from encode_pipeline.api.registered_result_models import (
    ArchiveRegisteredResultRequest,
    RegisterRegisteredResultRequest,
    RegisteredResultArchiveResponse,
    RegisteredResultDetailResponse,
    RegisteredResultListResponse,
    RegisteredResultRevisionResponse,
    ReviseRegisteredResultRequest,
)
from encode_pipeline.platform.authentication import AuthenticatedPrincipal
from encode_pipeline.services.registered_result_bytes import discard_private_copy
from encode_pipeline.services.registered_results import (
    RegisteredResultAPIError,
    RegisteredResultDownload,
    RegisteredResultService,
)


router = APIRouter(tags=["registered-results"])

_PROJECT_ID = Path(pattern=r"^prj_[0-9a-f]{32}$")
_RESULT_ID = Path(pattern=r"^rreg_[0-9a-f]{32}$")
_REVISION_ID = Path(pattern=r"^rrev_[0-9a-f]{32}$")
_IDEMPOTENCY_KEY = Header(
    alias="Idempotency-Key",
    min_length=1,
    max_length=128,
    pattern=r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,127}$",
)
_JSON_ERRORS = {
    400: {"model": ValidationResponse},
    404: {"model": ValidationResponse},
    409: {"model": ValidationResponse},
    500: {"model": ValidationResponse},
}


def get_registered_result_service(request: Request) -> RegisteredResultService:
    """Return the catalog service bound to this process."""

    return request.app.state.registered_result_service


def _failure(exc: RegisteredResultAPIError) -> JSONResponse:
    if exc.status_code in {401, 403}:
        raise AuthApiError(exc.status_code, exc.code, exc.message)
    body = ValidationResponse(
        ok=False,
        workflow_id=None,
        value=None,
        snapshot=None,
        issues=[
            IssueResponse(code=exc.code, message=exc.message, severity="error"),
        ],
    )
    return JSONResponse(
        status_code=exc.status_code,
        content=body.model_dump(mode="json"),
    )


@router.get(
    "/catalog/projects/{project_id}/registered-results",
    operation_id="catalogListRegisteredResults",
    response_model=RegisteredResultListResponse,
    responses=_JSON_ERRORS,
)
async def catalog_list_registered_results(
    project_id: str = _PROJECT_ID,
    cursor: str | None = Query(default=None, max_length=4096),
    limit: int = Query(default=20, ge=1, le=50),
    include_archived: bool = False,
    principal: AuthenticatedPrincipal = Depends(require_principal),
    service: RegisteredResultService = Depends(get_registered_result_service),
) -> RegisteredResultListResponse | JSONResponse:
    """List current revisions for one visible project."""

    try:
        results, next_cursor = service.list_current(
            principal,
            project_id,
            cursor=cursor,
            limit=limit,
            include_archived=include_archived,
        )
    except RegisteredResultAPIError as exc:
        return _failure(exc)
    return RegisteredResultListResponse(
        ok=True,
        results=[
            RegisteredResultRevisionResponse.model_validate(item) for item in results
        ],
        next_cursor=next_cursor,
        issues=[],
    )


@router.post(
    "/catalog/projects/{project_id}/registered-results",
    operation_id="catalogRegisterRegisteredResult",
    response_model=RegisteredResultDetailResponse,
    responses=_JSON_ERRORS,
)
async def catalog_register_registered_result(
    body: RegisterRegisteredResultRequest,
    project_id: str = _PROJECT_ID,
    idempotency_key: str = _IDEMPOTENCY_KEY,
    principal: AuthenticatedPrincipal = Depends(require_administrator),
    service: RegisteredResultService = Depends(get_registered_result_service),
) -> RegisteredResultDetailResponse | JSONResponse:
    """Record one external or downstream file without copying or running it."""

    payload = body.model_dump(exclude_none=True)
    relative_path = payload.pop("pool_relative_path")
    stable_key = payload.pop("stable_key")
    try:
        revision = service.register(
            principal,
            project_id,
            idempotency_key=idempotency_key,
            stable_key=stable_key,
            metadata_payload=payload,
            pool_relative_path=relative_path,
        )
    except RegisteredResultAPIError as exc:
        return _failure(exc)
    return _detail(revision)


@router.get(
    "/catalog/projects/{project_id}/registered-results/{result_id}/revisions",
    operation_id="catalogListRegisteredResultRevisions",
    response_model=RegisteredResultListResponse,
    responses=_JSON_ERRORS,
)
async def catalog_list_registered_result_revisions(
    project_id: str = _PROJECT_ID,
    result_id: str = _RESULT_ID,
    cursor: str | None = Query(default=None, max_length=4096),
    limit: int = Query(default=20, ge=1, le=50),
    principal: AuthenticatedPrincipal = Depends(require_principal),
    service: RegisteredResultService = Depends(get_registered_result_service),
) -> RegisteredResultListResponse | JSONResponse:
    """List every revision of one result, oldest first."""

    try:
        results, next_cursor = service.list_revisions(
            principal,
            project_id,
            result_id,
            cursor=cursor,
            limit=limit,
        )
    except RegisteredResultAPIError as exc:
        return _failure(exc)
    return RegisteredResultListResponse(
        ok=True,
        results=[
            RegisteredResultRevisionResponse.model_validate(item) for item in results
        ],
        next_cursor=next_cursor,
        issues=[],
    )


@router.post(
    "/catalog/projects/{project_id}/registered-results/{result_id}/revisions",
    operation_id="catalogReviseRegisteredResult",
    response_model=RegisteredResultDetailResponse,
    responses=_JSON_ERRORS,
)
async def catalog_revise_registered_result(
    body: ReviseRegisteredResultRequest,
    project_id: str = _PROJECT_ID,
    result_id: str = _RESULT_ID,
    idempotency_key: str = _IDEMPOTENCY_KEY,
    principal: AuthenticatedPrincipal = Depends(require_administrator),
    service: RegisteredResultService = Depends(get_registered_result_service),
) -> RegisteredResultDetailResponse | JSONResponse:
    """Append one revision. A rename re-reads the same bytes before it commits."""

    try:
        revision = service.revise(
            principal,
            project_id,
            result_id,
            idempotency_key=idempotency_key,
            expected_revision_number=body.expected_revision_number,
            display_name=body.display_name,
            result_format=body.format,
            origin=body.origin,
            method=body.method,
            reference=body.reference,
            pool_relative_path=body.pool_relative_path,
        )
    except RegisteredResultAPIError as exc:
        return _failure(exc)
    return _detail(revision)


@router.post(
    "/catalog/projects/{project_id}/registered-results/{result_id}/archive",
    operation_id="catalogArchiveRegisteredResult",
    response_model=RegisteredResultArchiveResponse,
    responses=_JSON_ERRORS,
)
async def catalog_archive_registered_result(
    body: ArchiveRegisteredResultRequest,
    project_id: str = _PROJECT_ID,
    result_id: str = _RESULT_ID,
    idempotency_key: str = _IDEMPOTENCY_KEY,
    principal: AuthenticatedPrincipal = Depends(require_administrator),
    service: RegisteredResultService = Depends(get_registered_result_service),
) -> RegisteredResultArchiveResponse | JSONResponse:
    """Stop new revisions. Historical reads stay available."""

    try:
        archived = service.archive(
            principal,
            project_id,
            result_id,
            idempotency_key=idempotency_key,
            expected_revision_number=body.expected_revision_number,
        )
    except RegisteredResultAPIError as exc:
        return _failure(exc)
    return RegisteredResultArchiveResponse(
        ok=True,
        result_id=archived["result_id"],
        archived_at=archived["archived_at"],
        issues=[],
    )


@router.get(
    "/catalog/projects/{project_id}/registered-result-revisions/{result_revision_id}",
    operation_id="catalogGetRegisteredResultRevision",
    response_model=RegisteredResultDetailResponse,
    responses=_JSON_ERRORS,
)
async def catalog_get_registered_result_revision(
    project_id: str = _PROJECT_ID,
    result_revision_id: str = _REVISION_ID,
    principal: AuthenticatedPrincipal = Depends(require_principal),
    service: RegisteredResultService = Depends(get_registered_result_service),
) -> RegisteredResultDetailResponse | JSONResponse:
    """Read one exact revision, including a superseded name."""

    try:
        revision = service.get_revision(principal, project_id, result_revision_id)
    except RegisteredResultAPIError as exc:
        return _failure(exc)
    return _detail(revision)


@router.get(
    "/catalog/projects/{project_id}/registered-result-revisions/{result_revision_id}/content",
    response_class=StreamingResponse,
    response_model=None,
    operation_id="catalogDownloadRegisteredResultRevision",
    responses={
        200: {
            "description": "Verified private copy of the requested revision.",
            "content": {
                "application/octet-stream": {
                    "schema": {"type": "string", "format": "binary"}
                }
            },
        },
        400: {"model": ValidationResponse},
        404: {"model": ValidationResponse},
        409: {"model": ValidationResponse},
        500: {"model": ValidationResponse},
    },
)
async def catalog_download_registered_result_revision(
    project_id: str = _PROJECT_ID,
    result_revision_id: str = _REVISION_ID,
    principal: AuthenticatedPrincipal = Depends(require_principal),
    service: RegisteredResultService = Depends(get_registered_result_service),
) -> StreamingResponse | JSONResponse:
    """Send the request's verified copy. The source file is not the body."""

    try:
        plan = service.prepare_download(principal, project_id, result_revision_id)
    except RegisteredResultAPIError as exc:
        return _failure(exc)
    try:
        return StreamingResponse(
            _iter_copy(plan),
            media_type="application/octet-stream",
            headers={
                "Accept-Ranges": "none",
                "Cache-Control": "no-store",
                "Content-Disposition": f'attachment; filename="{plan.filename}"',
                "Content-Length": str(plan.size_bytes),
                "X-Content-Type-Options": "nosniff",
            },
        )
    except BaseException:
        _release_download(plan)
        raise


def _detail(revision: dict) -> RegisteredResultDetailResponse:
    return RegisteredResultDetailResponse(
        ok=True,
        revision=RegisteredResultRevisionResponse.model_validate(revision),
        issues=[],
    )


def _iter_copy(plan: RegisteredResultDownload):
    try:
        while True:
            chunk = os.read(plan.descriptor, 1024 * 1024)
            if not chunk:
                break
            yield chunk
    finally:
        _release_download(plan)


def _release_download(plan: RegisteredResultDownload) -> None:
    with suppress(OSError):
        os.close(plan.descriptor)
    discard_private_copy(plan.copy_path, plan.directory)
