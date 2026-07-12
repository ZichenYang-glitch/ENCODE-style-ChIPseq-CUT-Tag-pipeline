"""Fresh worker-side dependency composition for one durable execution job."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from types import TracebackType

from encode_pipeline.persistence.runtime import RunPersistence, open_run_persistence
from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.services.command_builder import CommandBuilder
from encode_pipeline.services.artifact_extraction import ArtifactExtractionService
from encode_pipeline.services.local_execution import LocalExecutionService
from encode_pipeline.services.local_run_driver import LocalRunDriver
from encode_pipeline.services.materialization import WorkspaceMaterializer
from encode_pipeline.services.planning import ExecutionPlanner, WorkspacePlanner
from encode_pipeline.services.preflight import LocalPreflightService
from encode_pipeline.services.process_runner import ProcessRunner
from encode_pipeline.services.runs import RunService
from encode_pipeline.services.workflow_builds import WorkflowBuildIdentityProvider
from encode_pipeline.workers.settings import WorkerSettings, load_worker_settings
from encode_pipeline.workers.timeouts import WorkerHardTimeout


@dataclass(frozen=True)
class WorkerRuntime:
    """Job-local dependencies rebuilt from environment and durable state."""

    settings: WorkerSettings
    persistence: RunPersistence
    registry: WorkflowRegistry
    run_service: RunService
    build_identity_provider: WorkflowBuildIdentityProvider
    execution_planner: ExecutionPlanner
    workspace_planner: WorkspacePlanner
    materializer: WorkspaceMaterializer
    command_builder: CommandBuilder
    local_run_driver: LocalRunDriver
    local_execution_service: LocalExecutionService
    artifact_extraction_service: ArtifactExtractionService
    preflight_service: LocalPreflightService

    def close(self) -> None:
        """Release database resources owned by this job-local runtime."""
        self.persistence.close()

    def __enter__(self) -> WorkerRuntime:
        return self

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc_value: BaseException | None,
        traceback: TracebackType | None,
    ) -> None:
        self.close()


def open_worker_runtime(
    settings: WorkerSettings | None = None,
    *,
    project_root: Path | None = None,
    build_identity_provider: WorkflowBuildIdentityProvider | None = None,
) -> WorkerRuntime:
    """Reopen SQLite and reconstruct all adapter/execution dependencies."""
    from encode_pipeline.services.defaults import (
        create_default_command_builder,
        create_default_artifact_extraction_service,
        create_default_execution_planner,
        create_default_local_execution_service,
        create_default_local_run_driver,
        create_default_run_service,
        create_default_workflow_registry,
        create_default_workspace_materializer,
        create_default_workspace_planner,
        create_default_workflow_build_identity_provider,
    )

    resolved_settings = settings if settings is not None else load_worker_settings()
    if not isinstance(resolved_settings, WorkerSettings):
        raise ValueError("settings must be a WorkerSettings instance or None")

    persistence = open_run_persistence(resolved_settings.database_url)
    try:
        source_hint = (
            build_identity_provider.project_root
            if build_identity_provider is not None
            else project_root
        )
        registry = create_default_workflow_registry(project_root=source_hint)
        if build_identity_provider is None:
            build_identity_provider = create_default_workflow_build_identity_provider(
                registry=registry,
                project_root=project_root,
            )
        elif not isinstance(
            build_identity_provider,
            WorkflowBuildIdentityProvider,
        ):
            raise ValueError(
                "build_identity_provider must be a WorkflowBuildIdentityProvider"
            )
        elif (
            project_root is not None
            and build_identity_provider.project_root != project_root
        ):
            raise ValueError(
                "project_root must match build_identity_provider.project_root"
            )
        source_project_root = build_identity_provider.project_root
        run_service = create_default_run_service(
            registry=registry,
            repository=persistence.repository,
        )
        execution_planner = create_default_execution_planner(run_service)
        workspace_planner = create_default_workspace_planner(registry=registry)
        materializer = create_default_workspace_materializer()
        command_builder = create_default_command_builder(
            registry=registry,
            project_root=source_project_root,
        )
        local_run_driver = create_default_local_run_driver(
            run_service,
            workspace_root=resolved_settings.workspace_root,
            materializer=materializer,
            command_builder=command_builder,
            process_runner=ProcessRunner(
                timeout_seconds=resolved_settings.job_timeout_seconds,
                passthrough_exceptions=(WorkerHardTimeout,),
            ),
        )
        local_execution_service = create_default_local_execution_service(
            run_service,
            execution_planner=execution_planner,
            workspace_planner=workspace_planner,
            command_builder=command_builder,
            local_run_driver=local_run_driver,
        )
        artifact_extraction_service = create_default_artifact_extraction_service(
            run_service=run_service,
            registry=registry,
            build_identity_provider=build_identity_provider,
            workspace_root=resolved_settings.workspace_root,
        )
        preflight_service = LocalPreflightService(
            run_service=run_service,
            execution_planner=execution_planner,
            workspace_planner=workspace_planner,
            local_run_driver=local_run_driver,
            build_identity_provider=build_identity_provider,
        )
        return WorkerRuntime(
            settings=resolved_settings,
            persistence=persistence,
            registry=registry,
            run_service=run_service,
            build_identity_provider=build_identity_provider,
            execution_planner=execution_planner,
            workspace_planner=workspace_planner,
            materializer=materializer,
            command_builder=command_builder,
            local_run_driver=local_run_driver,
            local_execution_service=local_execution_service,
            artifact_extraction_service=artifact_extraction_service,
            preflight_service=preflight_service,
        )
    except BaseException:
        persistence.close()
        raise
