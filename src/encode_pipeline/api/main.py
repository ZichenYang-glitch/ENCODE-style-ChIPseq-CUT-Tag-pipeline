"""FastAPI app factory and exception handlers for the workflow platform API."""

from __future__ import annotations

import os
from contextlib import asynccontextmanager
from pathlib import Path
from typing import AsyncIterator

from fastapi import FastAPI, Request
from fastapi.exceptions import RequestValidationError
from fastapi.responses import JSONResponse

from encode_pipeline.api.models import (
    RunArtifactDownloadErrorResponse,
    RunArtifactsResponse,
    RunQcMetricsResponse,
    RunResponse,
    ValidationResponse,
)
from encode_pipeline.api.request_limits import AuthoringRequestLimitMiddleware
from encode_pipeline.api.routes import api_v1_router
from encode_pipeline.persistence import DATABASE_URL_ENV, open_run_persistence
from encode_pipeline.platform.results import Issue
from encode_pipeline.platform.adapters import MAX_AUTHORING_REQUEST_BYTES
from encode_pipeline.services.defaults import (
    create_default_agent_service,
    create_default_command_builder,
    create_default_local_run_driver,
    create_default_run_service,
    create_default_validation_service,
    create_default_workflow_build_identity_provider,
    create_default_workspace_planner,
    create_default_workflow_registry,
)
from encode_pipeline.services.artifact_downloads import ArtifactDownloadService
from encode_pipeline.services.planning import ExecutionPlanner
from encode_pipeline.services.preflight import LocalPreflightService
from encode_pipeline.services.run_cancellation import RunCancellationService
from encode_pipeline.services.run_submission import RunSubmissionService
from encode_pipeline.services.validated_inputs import (
    ValidatedInputService,
    ValidatedRunCreationService,
)
from encode_pipeline.workers.rq_queue import RqRunQueue
from encode_pipeline.workers.settings import (
    WORKSPACE_ROOT_ENV,
    load_worker_settings,
)


def create_app(
    *,
    database_url: str | None = None,
    workspace_root: Path | None = None,
    project_root: Path | None = None,
) -> FastAPI:
    """Return a configured FastAPI app with default platform services."""
    settings_environment = dict(os.environ)
    if database_url is not None:
        settings_environment[DATABASE_URL_ENV] = database_url
    if workspace_root is not None:
        settings_environment[WORKSPACE_ROOT_ENV] = str(workspace_root)
    worker_settings = load_worker_settings(settings_environment)
    persistence = open_run_persistence(worker_settings.database_url)
    run_queue = RqRunQueue(worker_settings)

    @asynccontextmanager
    async def lifespan(_app: FastAPI) -> AsyncIterator[None]:
        try:
            yield
        finally:
            run_queue.close()
            persistence.close()

    app = FastAPI(
        title="Workflow Platform API",
        version="0.3.0",
        description="Validation-first workflow platform API.",
        lifespan=lifespan,
    )
    app.add_middleware(
        AuthoringRequestLimitMiddleware,
        max_request_bytes=MAX_AUTHORING_REQUEST_BYTES,
    )

    registry = create_default_workflow_registry()
    build_identity_provider = create_default_workflow_build_identity_provider(
        registry=registry,
        project_root=project_root,
    )
    run_service = create_default_run_service(
        registry=registry,
        repository=persistence.repository,
    )
    validation_service = create_default_validation_service(registry=registry)
    validated_input_service = ValidatedInputService(
        registry=registry,
        validation_service=validation_service,
        build_identity_provider=build_identity_provider,
        repository=persistence.repository,
    )
    validated_run_creation_service = ValidatedRunCreationService(
        run_service=run_service,
        build_identity_provider=build_identity_provider,
    )
    run_submission_service = RunSubmissionService(
        run_service=run_service,
        run_queue=run_queue,
    )
    run_cancellation_service = RunCancellationService(
        run_service=run_service,
        run_queue=run_queue,
    )
    recovered_runs = run_service.recover_interrupted_runs()
    command_builder = create_default_command_builder(
        registry=registry,
        project_root=project_root,
    )
    local_run_driver = create_default_local_run_driver(
        run_service=run_service,
        workspace_root=worker_settings.workspace_root,
        command_builder=command_builder,
    )
    preflight_service = LocalPreflightService(
        run_service=run_service,
        execution_planner=ExecutionPlanner(run_service=run_service),
        workspace_planner=create_default_workspace_planner(registry=registry),
        local_run_driver=local_run_driver,
        build_identity_provider=build_identity_provider,
    )

    app.state.registry = registry
    app.state.persistence = persistence
    app.state.database_url = persistence.database_url
    app.state.workspace_root = worker_settings.workspace_root
    app.state.worker_settings = worker_settings
    app.state.run_queue = run_queue
    app.state.recovered_run_ids = tuple(run.run_id for run in recovered_runs)
    app.state.validation_service = validation_service
    app.state.validated_input_service = validated_input_service
    app.state.validated_run_creation_service = validated_run_creation_service
    app.state.agent_service = create_default_agent_service(registry=registry)
    app.state.run_service = run_service
    app.state.artifact_download_service = ArtifactDownloadService(
        run_service=run_service,
        workspace_root=worker_settings.workspace_root,
    )
    app.state.build_identity_provider = build_identity_provider
    app.state.run_submission_service = run_submission_service
    app.state.run_cancellation_service = run_cancellation_service
    app.state.local_run_driver = local_run_driver
    app.state.preflight_service = preflight_service
    # Stub driver is intentionally not attached in production.

    app.include_router(api_v1_router, prefix="/api/v1")
    app.add_exception_handler(RequestValidationError, _handle_request_validation_error)
    app.add_exception_handler(Exception, _handle_internal_server_error)

    return app


def _issue_from_request_validation_error(
    exc: RequestValidationError,
    workflow_id: str | None = None,
) -> Issue:
    """Build a stable API_REQUEST_INVALID issue from a Pydantic validation error."""
    return Issue(
        code="API_REQUEST_INVALID",
        message="Request body is invalid.",
        severity="error",
        path="body",
        source="api",
        technical_message=None,
        hint="Submit the fields required by this operation's OpenAPI schema.",
        context={"workflow_id": workflow_id} if workflow_id is not None else {},
    )


async def _handle_request_validation_error(
    request: Request,
    exc: RequestValidationError,
) -> JSONResponse:
    """Return 400 with the PR84 API_REQUEST_INVALID envelope."""
    route = request.scope.get("route")
    if getattr(route, "operation_id", None) == "listRunArtifacts":
        run_id = request.path_params.get("run_id", "")
        issue = Issue(
            code="API_REQUEST_INVALID",
            message="Artifact query parameters are invalid.",
            severity="error",
            path="limit",
            source="api",
            technical_message=None,
            hint="Use an integer limit between 1 and 100.",
            context={},
        )
        body = RunArtifactsResponse(
            ok=False,
            run_id=run_id,
            artifacts=[],
            next_cursor=None,
            issues=[issue.to_dict()],
        )
        return JSONResponse(
            status_code=400,
            content=body.model_dump(mode="json", exclude_none=True),
        )
    if getattr(route, "operation_id", None) == "listRunQcMetrics":
        run_id = request.path_params.get("run_id", "")
        issue = Issue(
            code="API_REQUEST_INVALID",
            message="QC metric query parameters are invalid.",
            severity="error",
            path="query",
            source="api",
            technical_message=None,
            hint="Use a limit between 1 and 100 and a cursor returned by this endpoint.",
            context={},
        )
        body = RunQcMetricsResponse(
            ok=False,
            run_id=run_id,
            qc_metrics=[],
            next_cursor=None,
            issues=[issue.to_dict()],
        )
        return JSONResponse(
            status_code=400,
            content=body.model_dump(mode="json", exclude_none=True),
        )
    if getattr(route, "operation_id", None) == "createRun":
        issue = _issue_from_request_validation_error(
            exc,
            workflow_id=request.path_params.get("workflow_id"),
        )
        body = RunResponse(ok=False, run=None, issues=[issue.to_dict()])
        return JSONResponse(
            status_code=400,
            content=body.model_dump(mode="json"),
        )
    workflow_id = request.path_params.get("workflow_id")
    issue = _issue_from_request_validation_error(exc, workflow_id=workflow_id)
    body = ValidationResponse(
        ok=False,
        workflow_id=workflow_id,
        value=None,
        snapshot=None,
        issues=[issue.to_dict()],
    )
    return JSONResponse(status_code=400, content=body.model_dump())


async def _handle_internal_server_error(
    request: Request,
    exc: Exception,
) -> JSONResponse:
    """Return 500 with the PR84 INTERNAL_SERVER_ERROR envelope."""
    route = request.scope.get("route")
    if getattr(route, "operation_id", None) == "createRun":
        issue = Issue(
            code="INTERNAL_SERVER_ERROR",
            message="Run creation is temporarily unavailable.",
            severity="error",
            path="snapshot_id",
            source="runtime",
            technical_message=None,
            hint=None,
            context={},
        )
        body = RunResponse(ok=False, run=None, issues=[issue.to_dict()])
        return JSONResponse(
            status_code=500,
            content=body.model_dump(mode="json"),
        )
    if getattr(route, "operation_id", None) == "downloadRunArtifact":
        issue = Issue(
            code="INTERNAL_SERVER_ERROR",
            message="Artifact download is temporarily unavailable.",
            severity="error",
            path="artifact_id",
            source="runtime",
            technical_message=None,
            hint=None,
            context={},
        )
        body = RunArtifactDownloadErrorResponse(
            run_id=request.path_params.get("run_id", ""),
            artifact_id=request.path_params.get("artifact_id", ""),
            issues=[issue.to_dict()],
        )
        return JSONResponse(
            status_code=500,
            content=body.model_dump(mode="json", exclude_none=True),
        )
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
        snapshot=None,
        issues=[issue.to_dict()],
    )
    return JSONResponse(status_code=500, content=body.model_dump())
