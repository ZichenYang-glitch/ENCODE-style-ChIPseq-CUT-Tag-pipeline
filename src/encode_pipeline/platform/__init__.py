"""Workflow-platform shared primitives."""

from encode_pipeline.platform.adapters import (
    CommandSpec,
    DagEdge,
    DagNode,
    DagPreview,
    WorkflowAdapter,
    WorkflowCapabilities,
    WorkflowInputs,
    WorkflowMetadata,
    WorkflowSchema,
    WorkspacePlan,
)
from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.platform.results import Issue, IssueSeverity, Result

__all__ = [
    "CommandSpec",
    "DagEdge",
    "DagNode",
    "DagPreview",
    "Issue",
    "IssueSeverity",
    "Result",
    "WorkflowAdapter",
    "WorkflowCapabilities",
    "WorkflowInputs",
    "WorkflowMetadata",
    "WorkflowRegistry",
    "WorkflowSchema",
    "WorkspacePlan",
]
