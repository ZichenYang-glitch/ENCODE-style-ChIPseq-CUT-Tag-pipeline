"""Fresh worker-side dependency composition for one durable execution job."""

from __future__ import annotations

from dataclasses import dataclass
from types import TracebackType

from encode_pipeline.persistence.runtime import RunPersistence, open_run_persistence
from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.services.command_builder import CommandBuilder
from encode_pipeline.services.local_run_driver import LocalRunDriver
from encode_pipeline.services.materialization import WorkspaceMaterializer
from encode_pipeline.services.planning import ExecutionPlanner, WorkspacePlanner
from encode_pipeline.services.preflight import LocalPreflightService
from encode_pipeline.services.runs import RunService
from encode_pipeline.workers.settings import WorkerSettings, load_worker_settings


@dataclass(frozen=True)
class WorkerRuntime:
    """Job-local dependencies rebuilt from environment and durable state."""

    settings: WorkerSettings
    persistence: RunPersistence
    registry: WorkflowRegistry
    run_service: RunService
    execution_planner: ExecutionPlanner
    workspace_planner: WorkspacePlanner
    materializer: WorkspaceMaterializer
    command_builder: CommandBuilder
    local_run_driver: LocalRunDriver
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


def open_worker_runtime(settings: WorkerSettings | None = None) -> WorkerRuntime:
    """Reopen SQLite and reconstruct all adapter/execution dependencies."""
    from encode_pipeline.services.defaults import (
        create_default_command_builder,
        create_default_execution_planner,
        create_default_local_run_driver,
        create_default_run_service,
        create_default_workflow_registry,
        create_default_workspace_materializer,
        create_default_workspace_planner,
    )

    resolved_settings = settings if settings is not None else load_worker_settings()
    if not isinstance(resolved_settings, WorkerSettings):
        raise ValueError("settings must be a WorkerSettings instance or None")

    persistence = open_run_persistence(resolved_settings.database_url)
    try:
        registry = create_default_workflow_registry()
        run_service = create_default_run_service(
            registry=registry,
            repository=persistence.repository,
        )
        execution_planner = create_default_execution_planner(run_service)
        workspace_planner = create_default_workspace_planner(registry=registry)
        materializer = create_default_workspace_materializer()
        command_builder = create_default_command_builder(registry=registry)
        local_run_driver = create_default_local_run_driver(
            run_service,
            workspace_root=resolved_settings.workspace_root,
            materializer=materializer,
            command_builder=command_builder,
        )
        preflight_service = LocalPreflightService(
            run_service=run_service,
            execution_planner=execution_planner,
            workspace_planner=workspace_planner,
            local_run_driver=local_run_driver,
        )
        return WorkerRuntime(
            settings=resolved_settings,
            persistence=persistence,
            registry=registry,
            run_service=run_service,
            execution_planner=execution_planner,
            workspace_planner=workspace_planner,
            materializer=materializer,
            command_builder=command_builder,
            local_run_driver=local_run_driver,
            preflight_service=preflight_service,
        )
    except BaseException:
        persistence.close()
        raise
