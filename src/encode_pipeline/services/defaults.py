"""Default service composition for the bundled workflow adapter."""

from __future__ import annotations

from typing import TYPE_CHECKING

from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.services.validation import ValidationService


if TYPE_CHECKING:
    from encode_pipeline.services.agent import AgentService


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
) -> "RunService":
    """Return a fresh run service wired to the default registry.

    Args:
        registry: Optional existing registry. When omitted, a fresh default
            registry is created. Passing a shared registry lets ``create_app``
            compose all services from one adapter instance.
    """
    from encode_pipeline.services.runs import RunService

    if registry is None:
        registry = create_default_workflow_registry()
    return RunService(registry=registry)
