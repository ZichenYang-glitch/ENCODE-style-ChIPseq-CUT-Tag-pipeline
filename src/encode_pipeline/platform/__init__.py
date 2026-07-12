"""Workflow-platform shared primitives."""

from encode_pipeline.platform.adapters import (
    CommandSpec,
    DagEdge,
    DagNode,
    DagPreview,
    ExtractedArtifactCandidate,
    WorkflowAdapter,
    WorkflowCapabilities,
    WorkflowInputs,
    WorkflowMetadata,
    WorkflowSchema,
    WorkspacePlan,
)
from encode_pipeline.platform.builds import WorkflowBuildIdentity
from encode_pipeline.platform.execution import (
    RunExecutionClaim,
    RunExecutionAssignment,
    RunExecutionCancellationRequest,
    RunExecutionOwnership,
    RunExecutionStopAcknowledgement,
    build_execution_job_id,
)
from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.platform.results import Issue, IssueSeverity, Result
from encode_pipeline.platform.runs import (
    RunArtifactRef,
    RunEvent,
    RunLogChunk,
    RunRecord,
    RunStatus,
    can_transition,
    require_transition,
)

__all__ = [
    "CommandSpec",
    "DagEdge",
    "DagNode",
    "DagPreview",
    "ExtractedArtifactCandidate",
    "Issue",
    "IssueSeverity",
    "Result",
    "RunArtifactRef",
    "RunEvent",
    "RunExecutionClaim",
    "RunExecutionAssignment",
    "RunExecutionCancellationRequest",
    "RunExecutionOwnership",
    "RunExecutionStopAcknowledgement",
    "RunLogChunk",
    "RunRecord",
    "RunStatus",
    "WorkflowAdapter",
    "WorkflowBuildIdentity",
    "WorkflowCapabilities",
    "WorkflowInputs",
    "WorkflowMetadata",
    "WorkflowRegistry",
    "WorkflowSchema",
    "WorkspacePlan",
    "build_execution_job_id",
    "can_transition",
    "require_transition",
]
