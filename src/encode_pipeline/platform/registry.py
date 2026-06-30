"""Workflow-platform adapter registry primitives."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from types import MappingProxyType

from encode_pipeline.platform.adapters import WorkflowAdapter, WorkflowMetadata


class WorkflowRegistry:
    """Immutable adapter lookup by workflow ID."""

    __slots__ = ("__adapters_by_workflow_id", "__metadata")

    def __init__(self, adapters: Iterable[WorkflowAdapter] = ()) -> None:
        adapters_by_workflow_id: dict[str, WorkflowAdapter] = {}
        metadata: list[WorkflowMetadata] = []

        for adapter in adapters:
            if not isinstance(adapter, WorkflowAdapter):
                raise ValueError(
                    "WorkflowRegistry adapters must satisfy WorkflowAdapter"
                )
            if not isinstance(adapter.metadata, WorkflowMetadata):
                raise ValueError(
                    "WorkflowRegistry adapter metadata must be WorkflowMetadata"
                )

            workflow_id = adapter.metadata.workflow_id
            if workflow_id in adapters_by_workflow_id:
                raise ValueError(f"Duplicate workflow_id: {workflow_id}")

            adapters_by_workflow_id[workflow_id] = adapter
            metadata.append(adapter.metadata)

        self.__adapters_by_workflow_id: Mapping[str, WorkflowAdapter] = (
            MappingProxyType(adapters_by_workflow_id)
        )
        self.__metadata = tuple(metadata)

    def get(self, workflow_id: str) -> WorkflowAdapter:
        """Return the adapter for a workflow ID."""
        normalized = _normalize_workflow_id(workflow_id)
        try:
            return self.__adapters_by_workflow_id[normalized]
        except KeyError as exc:
            raise KeyError(normalized) from exc

    def has(self, workflow_id: str) -> bool:
        """Return whether a workflow ID is registered."""
        normalized = _normalize_workflow_id(workflow_id)
        return normalized in self.__adapters_by_workflow_id

    def list_metadata(self) -> tuple[WorkflowMetadata, ...]:
        """Return registered workflow metadata in registration order."""
        return self.__metadata


def _normalize_workflow_id(workflow_id: str) -> str:
    if not isinstance(workflow_id, str):
        raise ValueError("workflow_id must be a string")
    normalized = workflow_id.strip()
    if not normalized:
        raise ValueError("workflow_id must be non-empty")
    return normalized
