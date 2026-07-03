"""FastAPI app factory and exception handlers for the workflow platform API."""

from __future__ import annotations

from fastapi import FastAPI, Request
from fastapi.exceptions import RequestValidationError
from fastapi.responses import JSONResponse

from encode_pipeline.api.models import ValidationResponse
from encode_pipeline.api.routes import api_v1_router
from encode_pipeline.platform.results import Issue
from encode_pipeline.services.defaults import (
    create_default_agent_service,
    create_default_run_service,
    create_default_validation_service,
    create_default_workflow_registry,
)


def create_app() -> FastAPI:
    """Return a configured FastAPI app with default platform services."""
    app = FastAPI(
        title="Workflow Platform API",
        version="0.3.0",
        description="Validation-first workflow platform API.",
    )

    registry = create_default_workflow_registry()
    app.state.registry = registry
    app.state.validation_service = create_default_validation_service(registry=registry)
    app.state.agent_service = create_default_agent_service(registry=registry)
    app.state.run_service = create_default_run_service(registry=registry)

    app.include_router(api_v1_router, prefix="/api/v1")
    app.add_exception_handler(RequestValidationError, _handle_request_validation_error)
    app.add_exception_handler(Exception, _handle_internal_server_error)

    return app


def _issue_from_request_validation_error(
    exc: RequestValidationError,
    workflow_id: str | None = None,
) -> Issue:
    """Build a stable API_REQUEST_INVALID issue from a Pydantic validation error."""
    errors = exc.errors()
    sanitized = "; ".join(
        f"{'.'.join(str(loc) for loc in error.get('loc', []))}: {error.get('msg', '')}"
        for error in errors
    ) if errors else None

    return Issue(
        code="API_REQUEST_INVALID",
        message="Request body is invalid.",
        severity="error",
        path="body",
        source="api",
        technical_message=sanitized,
        hint="Submit an object with config, samples, and options fields.",
        context={"workflow_id": workflow_id} if workflow_id is not None else {},
    )


async def _handle_request_validation_error(
    request: Request,
    exc: RequestValidationError,
) -> JSONResponse:
    """Return 400 with the PR84 API_REQUEST_INVALID envelope."""
    workflow_id = request.path_params.get("workflow_id")
    issue = _issue_from_request_validation_error(exc, workflow_id=workflow_id)
    body = ValidationResponse(
        ok=False,
        workflow_id=workflow_id,
        value=None,
        issues=[issue.to_dict()],
    )
    return JSONResponse(status_code=400, content=body.model_dump())


async def _handle_internal_server_error(
    request: Request,
    exc: Exception,
) -> JSONResponse:
    """Return 500 with the PR84 INTERNAL_SERVER_ERROR envelope."""
    workflow_id = request.path_params.get("workflow_id")
    issue = Issue(
        code="INTERNAL_SERVER_ERROR",
        message="An internal server error occurred.",
        severity="error",
        path=None,
        source="runtime",
        technical_message=None,
        hint=None,
        context={},
    )
    body = ValidationResponse(
        ok=False,
        workflow_id=workflow_id,
        value=None,
        issues=[issue.to_dict()],
    )
    return JSONResponse(status_code=500, content=body.model_dump())
