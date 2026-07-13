"""Workflow validation service."""

from __future__ import annotations

from encode_pipeline.platform.adapters import VALIDATION_CAPABILITY, WorkflowInputs
from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.platform.results import Issue, Result


class ValidationService:
    """Resolve workflow adapters and delegate input validation."""

    def __init__(self, registry: WorkflowRegistry) -> None:
        if not isinstance(registry, WorkflowRegistry):
            raise ValueError("ValidationService registry must be WorkflowRegistry")
        self._registry = registry

    def validate(self, workflow_id: str, inputs: WorkflowInputs) -> Result[object]:
        """Validate workflow inputs through the registered adapter."""
        try:
            adapter = self._registry.get(workflow_id)
        except KeyError:
            return Result.failure(
                [
                    Issue(
                        code="WORKFLOW_NOT_FOUND",
                        message="Workflow was not found.",
                        source="registry",
                        path="workflow_id",
                        context={"workflow_id": workflow_id},
                    )
                ]
            )

        if VALIDATION_CAPABILITY not in adapter.capabilities.supports:
            return Result.failure(
                [
                    Issue(
                        code="WORKFLOW_CAPABILITY_UNSUPPORTED",
                        message="Workflow does not support validation.",
                        source="registry",
                        path="workflow.capabilities",
                        context={
                            "workflow_id": workflow_id,
                            "capability": VALIDATION_CAPABILITY,
                        },
                    )
                ]
            )

        return adapter.validate(inputs)
