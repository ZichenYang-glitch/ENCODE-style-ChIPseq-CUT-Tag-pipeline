"""FastAPI app factory and exception handlers for the workflow platform API."""

from __future__ import annotations

import os
from contextlib import asynccontextmanager
from pathlib import Path
from typing import AsyncIterator, Callable

from fastapi import FastAPI, Request
from fastapi.exceptions import RequestValidationError
from fastapi.responses import JSONResponse

from encode_pipeline.api.dependencies import AuthApiError
from encode_pipeline.api.models import (
    ArtifactPublicationDetailResponse,
    ArtifactPublicationIssueResponse,
    ArtifactPublicationListResponse,
    AuthErrorResponse,
    IssueResponse,
    RunArtifactDetailResponse,
    RunArtifactDownloadErrorResponse,
    RunArtifactsResponse,
    RunHistoryResponse,
    RunQcMetricsResponse,
    RunResponse,
    ValidationResponse,
)
from encode_pipeline.api.request_limits import AuthoringRequestLimitMiddleware
from encode_pipeline.api.routes import api_v1_router
from encode_pipeline.persistence import (
    DATABASE_URL_ENV,
    RunPersistence,
    create_session_factory,
    open_run_persistence,
)
from encode_pipeline.persistence.authentication import (
    SqlAlchemyAuthenticationRepository,
)
from encode_pipeline.platform.results import Issue
from encode_pipeline.platform.adapters import MAX_AUTHORING_REQUEST_BYTES
from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.services.authentication import BrowserSessionCookiePolicy
from encode_pipeline.services.authentication_service import (
    AccountAdministrationService,
    AuthenticationService,
)
from encode_pipeline.services.defaults import (
    create_default_agent_service,
    create_default_command_builder,
    create_default_local_run_driver,
    create_default_process_runner,
    create_default_run_service,
    create_default_validation_service,
    create_default_workflow_build_identity_provider,
    create_default_workspace_planner,
    create_default_workflow_registry,
)
from encode_pipeline.services.artifact_downloads import ArtifactDownloadService
from encode_pipeline.services.artifact_publications import (
    ArtifactPublicationQueryService,
)
from encode_pipeline.services.planning import ExecutionPlanner
from encode_pipeline.services.preflight import LocalPreflightService
from encode_pipeline.services.private_reference_profiles import (
    PrivateReferenceProfileConfigError,
    load_private_reference_profile_config,
)
from encode_pipeline.services.reference_profile_runtime import (
    ReferenceProfileBindingService,
    ReferenceProfileRuntimeResolver,
)
from encode_pipeline.services.reference_profiles import ReferenceProfileService
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


AUTH_SECURE_COOKIES_ENV = "ENCODE_PIPELINE_AUTH_SECURE_COOKIES"


def create_app(
    *,
    database_url: str | None = None,
    workspace_root: Path | None = None,
    project_root: Path | None = None,
    persistence_opener: Callable[[str | None], RunPersistence] | None = None,
    registry: WorkflowRegistry | None = None,
) -> FastAPI:
    """Return a configured FastAPI app with default platform services."""
    settings_environment = dict(os.environ)
    if database_url is not None:
        settings_environment[DATABASE_URL_ENV] = database_url
    if workspace_root is not None:
        settings_environment[WORKSPACE_ROOT_ENV] = str(workspace_root)
    worker_settings = load_worker_settings(settings_environment)
    opener = open_run_persistence if persistence_opener is None else persistence_opener
    if not callable(opener):
        raise ValueError("persistence_opener must be callable")
    persistence = opener(worker_settings.database_url)
    run_queue = RqRunQueue(worker_settings)
    resolved_project_root = (
        worker_settings.encode_runtime_root if project_root is None else project_root
    )

    @asynccontextmanager
    async def lifespan(_app: FastAPI) -> AsyncIterator[None]:
        try:
            yield
        finally:
            run_queue.close()
            persistence.close()

    app = FastAPI(
        title="HelixWeave API",
        version="0.3.0",
        description="Reproducible omics workflows, from inputs to evidence.",
        lifespan=lifespan,
    )
    app.add_middleware(
        AuthoringRequestLimitMiddleware,
        max_request_bytes=MAX_AUTHORING_REQUEST_BYTES,
    )

    registry = (
        create_default_workflow_registry(environ=settings_environment)
        if registry is None
        else registry
    )
    if not isinstance(registry, WorkflowRegistry):
        persistence.close()
        run_queue.close()
        raise ValueError("registry must be a WorkflowRegistry")
    build_identity_provider = create_default_workflow_build_identity_provider(
        registry=registry,
        project_root=resolved_project_root,
    )
    run_service = create_default_run_service(
        registry=registry,
        repository=persistence.repository,
    )
    artifact_publication_service = ArtifactPublicationQueryService(
        repository=persistence.repository,
    )
    validation_service = create_default_validation_service(registry=registry)

    def _private_reference_config():
        config_path = worker_settings.reference_profile_config
        if config_path is None:
            raise PrivateReferenceProfileConfigError
        return load_private_reference_profile_config(config_path)

    reference_profile_service = ReferenceProfileService(
        repository=persistence.reference_profile_repository,
        private_config_provider=_private_reference_config,
        adapter_provider=registry.get,
    )
    reference_profile_binding_service = ReferenceProfileBindingService(
        repository=persistence.reference_profile_repository,
        private_config_provider=_private_reference_config,
        adapter_provider=registry.get,
    )
    reference_profile_resolver = ReferenceProfileRuntimeResolver(
        persistence.repository,
        registry,
        reference_profile_binding_service,
    )
    validated_input_service = ValidatedInputService(
        registry=registry,
        validation_service=validation_service,
        build_identity_provider=build_identity_provider,
        repository=persistence.repository,
        reference_profile_binding_service=reference_profile_binding_service,
        reference_profile_catalog=reference_profile_service,
    )
    validated_run_creation_service = ValidatedRunCreationService(
        run_service=run_service,
        build_identity_provider=build_identity_provider,
        reference_profile_binding_service=reference_profile_binding_service,
    )
    run_submission_service = RunSubmissionService(
        run_service=run_service,
        run_queue=run_queue,
        build_identity_provider=build_identity_provider,
        reference_profile_resolver=reference_profile_resolver,
    )
    run_cancellation_service = RunCancellationService(
        run_service=run_service,
        run_queue=run_queue,
    )
    recovered_runs = run_service.recover_interrupted_runs()
    command_builder = create_default_command_builder(
        registry=registry,
        project_root=resolved_project_root,
        reference_profile_resolver=reference_profile_resolver,
        snakemake_executable=(
            None
            if worker_settings.encode_runner_root is None
            else worker_settings.encode_runner_root / "bin" / "snakemake"
        ),
        conda_prefix=worker_settings.encode_conda_prefix,
    )
    local_run_driver = create_default_local_run_driver(
        run_service=run_service,
        workspace_root=worker_settings.workspace_root,
        command_builder=command_builder,
        process_runner=create_default_process_runner(
            registry=registry,
            settings=worker_settings,
        ),
    )
    preflight_service = LocalPreflightService(
        run_service=run_service,
        execution_planner=ExecutionPlanner(run_service=run_service),
        workspace_planner=create_default_workspace_planner(
            registry=registry,
            reference_profile_resolver=reference_profile_resolver,
        ),
        local_run_driver=local_run_driver,
        build_identity_provider=build_identity_provider,
        reference_profile_resolver=reference_profile_resolver,
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
    app.state.reference_profile_service = reference_profile_service
    app.state.reference_profile_binding_service = reference_profile_binding_service
    app.state.reference_profile_resolver = reference_profile_resolver
    app.state.validated_run_creation_service = validated_run_creation_service
    app.state.agent_service = create_default_agent_service(registry=registry)
    app.state.run_service = run_service
    app.state.artifact_publication_service = artifact_publication_service
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

    authentication_repository = SqlAlchemyAuthenticationRepository(
        create_session_factory(persistence.engine)
    )
    app.state.authentication_repository = authentication_repository
    app.state.authentication_service = AuthenticationService(
        repository=authentication_repository
    )
    app.state.account_administration_service = AccountAdministrationService(
        repository=authentication_repository
    )
    app.state.auth_cookie_policy = BrowserSessionCookiePolicy(
        secure=settings_environment.get(AUTH_SECURE_COOKIES_ENV) == "1"
    )

    app.include_router(api_v1_router, prefix="/api/v1")
    app.add_exception_handler(RequestValidationError, _handle_request_validation_error)
    app.add_exception_handler(AuthApiError, _handle_auth_api_error)
    app.add_exception_handler(Exception, _handle_internal_server_error)

    return app


async def _handle_auth_api_error(
    request: Request,
    exc: AuthApiError,
) -> JSONResponse:
    """Return the stable authentication/authorization failure envelope."""
    del request
    body = AuthErrorResponse(
        issues=[
            IssueResponse(
                code=exc.code,
                message=exc.message,
                source="api",
            )
        ]
    )
    return JSONResponse(
        status_code=exc.status_code,
        content=body.model_dump(mode="json"),
    )


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
    if getattr(route, "operation_id", None) == "listArtifactPublications":
        issue = ArtifactPublicationIssueResponse(
            code="API_REQUEST_INVALID",
            message="Artifact publication query parameters are invalid.",
            severity="error",
            path="query",
            source="api",
            hint="Use the filters and pagination bounds documented by this endpoint.",
        )
        body = ArtifactPublicationListResponse(
            ok=False,
            publications=[],
            next_cursor=None,
            issues=[issue],
        )
        content = body.model_dump(mode="json", exclude_none=True)
        content["next_cursor"] = None
        return JSONResponse(status_code=400, content=content)
    if getattr(route, "operation_id", None) == "getArtifactPublication":
        issue = ArtifactPublicationIssueResponse(
            code="API_REQUEST_INVALID",
            message="Artifact publication identity is invalid.",
            severity="error",
            path="query",
            source="api",
            hint="Use an exact generation returned by the publication list endpoint.",
        )
        body = ArtifactPublicationDetailResponse(
            ok=False,
            publication=None,
            issues=[issue],
        )
        content = body.model_dump(mode="json", exclude_none=True)
        content["publication"] = None
        return JSONResponse(status_code=400, content=content)
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
            artifact_generation=None,
            artifacts=[],
            next_cursor=None,
            issues=[issue.to_dict()],
        )
        content = body.model_dump(mode="json", exclude_none=True)
        content["artifact_generation"] = None
        return JSONResponse(
            status_code=400,
            content=content,
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
            qc_generation=None,
            qc_metrics=[],
            next_cursor=None,
            issues=[issue.to_dict()],
        )
        content = body.model_dump(mode="json", exclude_none=True)
        content["qc_generation"] = None
        return JSONResponse(
            status_code=400,
            content=content,
        )
    if getattr(route, "operation_id", None) == "getRunArtifact":
        run_id = request.path_params.get("run_id", "")
        issue = Issue(
            code="API_REQUEST_INVALID",
            message="Artifact detail query parameters are invalid.",
            severity="error",
            path="query",
            source="api",
            technical_message=None,
            hint="Use an artifact generation returned by the artifact list endpoint.",
            context={},
        )
        body = RunArtifactDetailResponse(
            ok=False,
            run_id=run_id,
            artifact_generation=None,
            artifact=None,
            issues=[issue.to_dict()],
        )
        return JSONResponse(status_code=400, content=body.model_dump(mode="json"))
    if getattr(route, "operation_id", None) == "downloadRunArtifact":
        issue = Issue(
            code="API_REQUEST_INVALID",
            message="Artifact download query parameters are invalid.",
            severity="error",
            path="query",
            source="api",
            technical_message=None,
            hint="Use the generation and revision shown by the artifact detail endpoint.",
            context={},
        )
        body = RunArtifactDownloadErrorResponse(
            run_id=request.path_params.get("run_id", ""),
            artifact_id=request.path_params.get("artifact_id", ""),
            issues=[issue.to_dict()],
        )
        return JSONResponse(status_code=400, content=body.model_dump(mode="json"))
    if getattr(route, "operation_id", None) == "listRuns":
        issue = Issue(
            code="API_REQUEST_INVALID",
            message="Run history query parameters are invalid.",
            severity="error",
            path="query",
            source="api",
            technical_message=None,
            hint="Use supported workflow/status filters, a limit from 1 to 100, and a returned cursor.",
            context={},
        )
        body = RunHistoryResponse(
            ok=False,
            runs=[],
            next_cursor=None,
            issues=[issue.to_dict()],
        )
        content = body.model_dump(mode="json", exclude_none=True)
        content["next_cursor"] = None
        return JSONResponse(
            status_code=400,
            content=content,
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
    if getattr(route, "operation_id", None) == "listArtifactPublications":
        issue = ArtifactPublicationIssueResponse(
            code="INTERNAL_SERVER_ERROR",
            message="Artifact publications are temporarily unavailable.",
            severity="error",
            path="publications",
            source="runtime",
            hint=None,
        )
        body = ArtifactPublicationListResponse(
            ok=False,
            publications=[],
            next_cursor=None,
            issues=[issue],
        )
        content = body.model_dump(mode="json", exclude_none=True)
        content["next_cursor"] = None
        return JSONResponse(status_code=500, content=content)
    if getattr(route, "operation_id", None) == "getArtifactPublication":
        issue = ArtifactPublicationIssueResponse(
            code="INTERNAL_SERVER_ERROR",
            message="Artifact publication is temporarily unavailable.",
            severity="error",
            path="publication",
            source="runtime",
            hint=None,
        )
        body = ArtifactPublicationDetailResponse(
            ok=False,
            publication=None,
            issues=[issue],
        )
        content = body.model_dump(mode="json", exclude_none=True)
        content["publication"] = None
        return JSONResponse(status_code=500, content=content)
    if getattr(route, "operation_id", None) == "listRuns":
        issue = Issue(
            code="INTERNAL_SERVER_ERROR",
            message="Run history is temporarily unavailable.",
            severity="error",
            path="runs",
            source="runtime",
            technical_message=None,
            hint=None,
            context={},
        )
        body = RunHistoryResponse(
            ok=False,
            runs=[],
            next_cursor=None,
            issues=[issue.to_dict()],
        )
        content = body.model_dump(mode="json", exclude_none=True)
        content["next_cursor"] = None
        return JSONResponse(
            status_code=500,
            content=content,
        )
    if getattr(route, "operation_id", None) == "listRunQcMetrics":
        issue = Issue(
            code="INTERNAL_SERVER_ERROR",
            message="Run QC metrics are temporarily unavailable.",
            severity="error",
            path="qc_metrics",
            source="runtime",
            technical_message=None,
            hint=None,
            context={},
        )
        body = RunQcMetricsResponse(
            ok=False,
            run_id=request.path_params.get("run_id", ""),
            qc_generation=None,
            qc_metrics=[],
            next_cursor=None,
            issues=[issue.to_dict()],
        )
        content = body.model_dump(mode="json", exclude_none=True)
        content["qc_generation"] = None
        return JSONResponse(status_code=500, content=content)
    if getattr(route, "operation_id", None) in {
        "listRunArtifacts",
        "getRunArtifact",
    }:
        operation_id = getattr(route, "operation_id", None)
        issue = Issue(
            code="INTERNAL_SERVER_ERROR",
            message="Run artifacts are temporarily unavailable.",
            severity="error",
            path="artifacts" if operation_id == "listRunArtifacts" else "artifact_id",
            source="runtime",
            technical_message=None,
            hint=None,
            context={},
        )
        if operation_id == "listRunArtifacts":
            body = RunArtifactsResponse(
                ok=False,
                run_id=request.path_params.get("run_id", ""),
                artifact_generation=None,
                artifacts=[],
                next_cursor=None,
                issues=[issue.to_dict()],
            )
        else:
            body = RunArtifactDetailResponse(
                ok=False,
                run_id=request.path_params.get("run_id", ""),
                artifact_generation=None,
                artifact=None,
                issues=[issue.to_dict()],
            )
        content = body.model_dump(mode="json", exclude_none=True)
        content["artifact_generation"] = None
        if operation_id == "listRunArtifacts":
            content["artifacts"] = []
        else:
            content["artifact"] = None
        return JSONResponse(status_code=500, content=content)
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
