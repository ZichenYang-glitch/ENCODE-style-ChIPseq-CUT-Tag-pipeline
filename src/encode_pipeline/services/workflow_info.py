"""Workflow info service for listing and describing registered workflows."""

from __future__ import annotations

from encode_pipeline.platform.adapters import (
    WorkflowAdapter,
    WorkflowCapabilities,
    WorkflowMetadata,
    WorkflowSchema,
)
from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.platform.results import Issue, Result


class WorkflowInfoService:
    """Resolve registered workflow adapters and expose their metadata."""

    def __init__(self, registry: WorkflowRegistry) -> None:
        if not isinstance(registry, WorkflowRegistry):
            raise ValueError("WorkflowInfoService registry must be WorkflowRegistry")
        self._registry = registry

    def list_workflows(self) -> list[WorkflowMetadata]:
        """Return registered workflow metadata in registration order."""
        return list(self._registry.list_metadata())

    def get_schema(self, workflow_id: str) -> Result[WorkflowSchema]:
        """Return the schema for a registered workflow."""
        result = self._resolve_adapter(workflow_id)
        if result.is_failure:
            return result
        return Result.success(result.value.schema())

    def get_capabilities(self, workflow_id: str) -> Result[WorkflowCapabilities]:
        """Return the capabilities for a registered workflow."""
        result = self._resolve_adapter(workflow_id)
        if result.is_failure:
            return result
        return Result.success(result.value.capabilities)

    def _resolve_adapter(
        self,
        workflow_id: str,
    ) -> Result[WorkflowAdapter]:
        """Fetch an adapter from the registry or return a failure result."""
        try:
            return Result.success(self._registry.get(workflow_id))
        except KeyError:
            return Result.failure([
                Issue(
                    code="WORKFLOW_NOT_FOUND",
                    message="Workflow was not found.",
                    source="workflow_info",
                    path="workflow_id",
                    context={"workflow_id": workflow_id},
                )
            ])
