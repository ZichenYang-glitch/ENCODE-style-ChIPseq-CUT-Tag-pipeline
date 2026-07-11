"""Default service composition for the bundled workflow adapter."""

from __future__ import annotations

from typing import TYPE_CHECKING

from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.services.validation import ValidationService


if TYPE_CHECKING:
    from pathlib import Path

    from encode_pipeline.services.agent import AgentService
    from encode_pipeline.services.command_builder import CommandBuilder
    from encode_pipeline.services.local_run_driver import LocalRunDriver
    from encode_pipeline.services.local_execution import LocalExecutionService
    from encode_pipeline.services.materialization import WorkspaceMaterializer
    from encode_pipeline.services.planning import ExecutionPlanner, WorkspacePlanner
    from encode_pipeline.services.runs import RunService
    from encode_pipeline.services.run_repositories import RunRepository
    from encode_pipeline.services.process_runner import ProcessRunner
    from encode_pipeline.services.stub_execution_driver import StubExecutionDriver
    from encode_pipeline.services.workflow_builds import WorkflowBuildIdentityProvider


def create_default_workflow_registry() -> WorkflowRegistry:
    """Return a fresh registry containing the bundled ENCODE-style adapter."""
    from encode_pipeline.adapters import EncodeStyleWorkflowAdapter

    return WorkflowRegistry(adapters=[EncodeStyleWorkflowAdapter()])


def create_default_validation_service(
    registry: WorkflowRegistry | None = None,
) -> ValidationService:
    """Return a fresh validation service wired to the default registry.

    Args:
        registry: Optional existing registry. When omitted, a fresh default
            registry is created. Passing a shared registry lets ``create_app``
            compose all services from one adapter instance.
    """
    if registry is None:
        registry = create_default_workflow_registry()
    return ValidationService(registry=registry)


def create_default_workflow_build_identity_provider(
    registry: WorkflowRegistry | None = None,
    *,
    project_root: Path | None = None,
) -> WorkflowBuildIdentityProvider:
    """Return a build identity provider for the bundled source tree."""
    from encode_pipeline.services.workflow_builds import WorkflowBuildIdentityProvider

    if registry is None:
        registry = create_default_workflow_registry()
    return WorkflowBuildIdentityProvider(registry, project_root=project_root)


def create_default_agent_service(
    registry: WorkflowRegistry | None = None,
) -> "AgentService":
    """Return a fresh agent service wired to the default registry and LLM client.

    Composes the default registry, workflow info service, validation service,
    an LLM client selected from environment variables by ``create_llm_client()``
    (defaulting to a deterministic mock client), and the PR96 safety components
    (redaction policy, output filter, and bounded in-memory audit sink). This
    function uses lazy imports to preserve the lazy-import guarantees of the
    services package.

    Args:
        registry: Optional existing registry. When omitted, a fresh default
            registry is created. Passing a shared registry lets ``create_app``
            compose all services from one adapter instance.
    """
    from encode_pipeline.services.agent import AgentService
    from encode_pipeline.services.agent_audit import InMemoryAuditSink
    from encode_pipeline.services.agent_output_filter import OutputFilter
    from encode_pipeline.services.agent_redaction import RedactionPolicy
    from encode_pipeline.services.llm_factory import create_llm_client
    from encode_pipeline.services.validation import ValidationService
    from encode_pipeline.services.workflow_info import WorkflowInfoService

    if registry is None:
        registry = create_default_workflow_registry()
    workflow_info = WorkflowInfoService(registry=registry)
    validation_service = ValidationService(registry=registry)
    llm_client = create_llm_client()
    return AgentService(
        workflow_info,
        validation_service,
        llm_client,
        redaction_policy=RedactionPolicy(),
        output_filter=OutputFilter(),
        audit_sink=InMemoryAuditSink(),
    )


def create_default_run_service(
    registry: WorkflowRegistry | None = None,
    *,
    repository: "RunRepository | None" = None,
) -> "RunService":
    """Return a fresh run service wired to the default registry.

    Args:
        registry: Optional existing registry. When omitted, a fresh default
            registry is created. Passing a shared registry lets ``create_app``
            compose all services from one adapter instance.
        repository: Optional run storage backend. Omitted callers retain the
            lightweight in-memory service used by isolated domain tests.
    """
    from encode_pipeline.services.runs import RunService

    if registry is None:
        registry = create_default_workflow_registry()
    return RunService(registry=registry, repository=repository)


def create_default_local_run_driver(
    run_service: RunService,
    *,
    workspace_root: Path | None = None,
    materializer: WorkspaceMaterializer | None = None,
    command_builder: CommandBuilder | None = None,
    process_runner: ProcessRunner | None = None,
) -> LocalRunDriver:
    """Return a local run driver wired to the given run service.

    Args:
        run_service: Required run service for lifeccle management.
        workspace_root: Parent directory for per-run workspace directories.
            Defaults to ``Path.home() / ".encode-pipeline" / "workspaces"``.
        materializer: Optional workspace materializer. Defaults to
            ``WorkspaceMaterializer()``.
        command_builder: Optional command builder. Defaults to
            ``create_default_command_builder()``.
        process_runner: Optional controlled subprocess runner.
    """
    from pathlib import Path

    from encode_pipeline.services.local_run_driver import LocalRunDriver
    from encode_pipeline.services.materialization import WorkspaceMaterializer

    if workspace_root is None:
        workspace_root = Path.home() / ".encode-pipeline" / "workspaces"
    if materializer is None:
        materializer = WorkspaceMaterializer()
    if command_builder is None:
        command_builder = create_default_command_builder()
    return LocalRunDriver(
        run_service=run_service,
        materializer=materializer,
        command_builder=command_builder,
        workspace_root=workspace_root,
        process_runner=process_runner,
    )


def create_default_local_execution_service(
    run_service: RunService,
    *,
    execution_planner: ExecutionPlanner,
    workspace_planner: WorkspacePlanner,
    command_builder: CommandBuilder,
    local_run_driver: LocalRunDriver,
) -> LocalExecutionService:
    """Return durable execution orchestration from worker-local dependencies."""
    from encode_pipeline.services.local_execution import LocalExecutionService

    return LocalExecutionService(
        run_service=run_service,
        execution_planner=execution_planner,
        workspace_planner=workspace_planner,
        command_builder=command_builder,
        local_run_driver=local_run_driver,
    )


def create_default_stub_execution_driver(
    run_service: RunService,
) -> StubExecutionDriver:
    """Return a stub execution driver wired to the given run service."""
    from encode_pipeline.services.stub_execution_driver import StubExecutionDriver

    return StubExecutionDriver(run_service=run_service)


def create_default_execution_planner(
    run_service: RunService,
) -> ExecutionPlanner:
    """Return an execution planner wired to the given run service."""
    from encode_pipeline.services.planning import ExecutionPlanner

    return ExecutionPlanner(run_service=run_service)


def create_default_workspace_planner(
    registry: WorkflowRegistry | None = None,
) -> WorkspacePlanner:
    """Return a fresh workspace planner wired to the given registry."""
    from encode_pipeline.services.planning import WorkspacePlanner

    if registry is None:
        registry = create_default_workflow_registry()
    return WorkspacePlanner(registry=registry)


def create_default_workspace_materializer() -> "WorkspaceMaterializer":
    """Return a fresh workspace materializer."""
    from encode_pipeline.services.materialization import WorkspaceMaterializer

    return WorkspaceMaterializer()


def create_default_command_builder(
    registry: WorkflowRegistry | None = None,
    *,
    project_root: Path | None = None,
) -> "CommandBuilder":
    """Return a fresh command builder wired to the default registry.

    Args:
        registry: Optional existing registry. When omitted, a fresh default
            registry is created. Passing a shared registry lets ``create_app``
            compose all services from one adapter instance.
        project_root: Optional absolute root containing ``workflow/Snakefile``.
            Omitted callers use the bundled project root.
    """
    from encode_pipeline.services.command_builder import CommandBuilder

    if registry is None:
        registry = create_default_workflow_registry()
    return CommandBuilder(registry=registry, project_root=project_root)
