"""Workflow info service for listing and describing registered workflows."""

from __future__ import annotations

from encode_pipeline.platform.adapters import (
    EXECUTION_CAPABILITY_NAMES,
    WorkflowAdapter,
    WorkflowAvailability,
    WorkflowAvailabilityProvidingAdapter,
    WorkflowCapabilities,
    WorkflowDescriptor,
    WorkflowMetadata,
    WorkflowSchema,
    WorkflowUpstreamIdentityProvidingAdapter,
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

    def list_descriptors(self) -> list[WorkflowDescriptor]:
        """Return safe product descriptors in registration order."""
        descriptors: list[WorkflowDescriptor] = []
        for metadata in self._registry.list_metadata():
            result = self.get_descriptor(metadata.workflow_id)
            if result.is_failure or result.value is None:
                continue
            descriptors.append(result.value)
        return descriptors

    def get_descriptor(self, workflow_id: str) -> Result[WorkflowDescriptor]:
        """Return one safe product descriptor for a registered workflow."""
        result = self._resolve_adapter(workflow_id)
        if result.is_failure:
            return result
        adapter = result.value
        try:
            schema_version = adapter.schema().schema_version
            availability = resolve_workflow_availability(adapter)
            capabilities = effective_workflow_capabilities(adapter, availability)
            upstream_identity = (
                adapter.upstream_identity
                if isinstance(adapter, WorkflowUpstreamIdentityProvidingAdapter)
                else None
            )
            descriptor = WorkflowDescriptor(
                metadata=adapter.metadata,
                schema_version=schema_version,
                capabilities=capabilities,
                upstream_identity=upstream_identity,
                availability=availability,
            )
        except Exception:
            return Result.failure(
                [
                    Issue(
                        code="WORKFLOW_DESCRIPTOR_UNAVAILABLE",
                        message="Workflow product information is unavailable.",
                        source="workflow_info",
                        path="workflow",
                        context={"workflow_id": workflow_id},
                    )
                ]
            )
        return Result.success(descriptor)

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
            return Result.failure(
                [
                    Issue(
                        code="WORKFLOW_NOT_FOUND",
                        message="Workflow was not found.",
                        source="workflow_info",
                        path="workflow_id",
                        context={"workflow_id": workflow_id},
                    )
                ]
            )


def resolve_workflow_availability(adapter: WorkflowAdapter) -> WorkflowAvailability:
    """Return current safe availability, failing closed for provider faults."""
    if not isinstance(adapter, WorkflowAvailabilityProvidingAdapter):
        return WorkflowAvailability()
    try:
        availability = adapter.execution_availability()
    except Exception:
        availability = None
    if not isinstance(availability, WorkflowAvailability):
        return WorkflowAvailability(
            execution="unavailable",
            reason_code="WORKFLOW_EXECUTION_UNAVAILABLE",
        )
    return availability


def effective_workflow_capabilities(
    adapter: WorkflowAdapter,
    availability: WorkflowAvailability | None = None,
) -> WorkflowCapabilities:
    """Hide execution capabilities unless current admission is available."""
    current = availability or resolve_workflow_availability(adapter)
    if current.execution == "available":
        return adapter.capabilities
    return WorkflowCapabilities(
        supports=tuple(
            capability
            for capability in adapter.capabilities.supports
            if capability not in EXECUTION_CAPABILITY_NAMES
        )
    )
