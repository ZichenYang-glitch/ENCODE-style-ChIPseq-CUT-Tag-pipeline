"""Default service composition for the bundled workflow adapter."""

from __future__ import annotations

from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.services.validation import ValidationService


def create_default_workflow_registry() -> WorkflowRegistry:
    """Return a fresh registry containing the bundled ENCODE-style adapter."""
    from encode_pipeline.adapters import EncodeStyleWorkflowAdapter

    return WorkflowRegistry(adapters=[EncodeStyleWorkflowAdapter()])


def create_default_validation_service() -> ValidationService:
    """Return a fresh validation service wired to the default registry."""
    return ValidationService(registry=create_default_workflow_registry())
