"""Contract-only adapter for pinned nf-core/rnaseq input semantics."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from encode_pipeline import __version__
from encode_pipeline.adapters.bulk_rnaseq.authoring import (
    build_bulk_rnaseq_authoring_schema,
)
from encode_pipeline.adapters.bulk_rnaseq.validation import (
    validate_bulk_rnaseq_inputs,
)
from encode_pipeline.adapters.bulk_rnaseq.execution import (
    BulkRnaSeqExecutionBinding,
    build_bulk_rnaseq_command,
    capture_bulk_rnaseq_build_identity,
    doctor_bulk_rnaseq_runtime,
    plan_bulk_rnaseq_workspace,
)
from encode_pipeline.platform.adapters import (
    COMMAND_CAPABILITY,
    INPUT_AUTHORING_CAPABILITY,
    VALIDATION_CAPABILITY,
    WORKSPACE_PLAN_CAPABILITY,
    CommandSpec,
    DagPreview,
    ExtractedArtifactCandidate,
    WorkflowCapabilities,
    WorkflowInputs,
    WorkflowMetadata,
    WorkflowSchema,
    WorkspacePlan,
)
from encode_pipeline.platform.results import Issue, Result


WORKFLOW_ID = "bulk-rnaseq"


class BulkRnaSeqWorkflowAdapter:
    """Expose contract-only or explicitly composed offline runtime behavior."""

    metadata = WorkflowMetadata(
        workflow_id=WORKFLOW_ID,
        name="Bulk RNA-seq",
        version=__version__ or "0.0.0",
        description=(
            "Offline-composable HelixWeave adapter for pinned nf-core/rnaseq 3.26.0."
        ),
        engines=("nextflow",),
        tags=("bulk-rnaseq", "illumina", "nf-core"),
    )
    capabilities = WorkflowCapabilities(
        supports=(VALIDATION_CAPABILITY, INPUT_AUTHORING_CAPABILITY)
    )

    def __init__(
        self,
        *,
        execution: BulkRnaSeqExecutionBinding | None = None,
    ) -> None:
        if execution is not None and not isinstance(
            execution, BulkRnaSeqExecutionBinding
        ):
            raise ValueError("execution must be a BulkRnaSeqExecutionBinding or None")
        self._execution = execution
        supports = [VALIDATION_CAPABILITY, INPUT_AUTHORING_CAPABILITY]
        if execution is not None:
            supports.extend((WORKSPACE_PLAN_CAPABILITY, COMMAND_CAPABILITY))
        self.capabilities = WorkflowCapabilities(supports=tuple(supports))

    def schema(self) -> WorkflowSchema:
        """Return the fresh versioned authoring contract."""
        return build_bulk_rnaseq_authoring_schema()

    def validate(self, inputs: WorkflowInputs) -> Result[object]:
        """Validate structure and science without probing submitted paths."""
        return validate_bulk_rnaseq_inputs(inputs)

    def preview_dag(self, inputs: WorkflowInputs) -> Result[DagPreview]:
        """Remain unsupported until the pinned runtime is composed."""
        return _unsupported("preview_dag")

    def plan_workspace(
        self,
        inputs: WorkflowInputs,
        workspace: str | Path,
    ) -> Result[WorkspacePlan]:
        """Build deterministic workspace bytes only for a composed runtime."""
        if self._execution is None:
            return _unsupported("plan_workspace")
        return plan_bulk_rnaseq_workspace(
            inputs,
            workspace,
            binding=self._execution,
            adapter_version=self.metadata.version,
        )

    def build_command(
        self,
        plan: WorkspacePlan,
        workspace: str | Path,
    ) -> Result[CommandSpec]:
        """Build the fixed offline Nextflow argv for a composed runtime."""
        if self._execution is None:
            return _unsupported("build_command")
        return build_bulk_rnaseq_command(
            plan,
            workspace,
            binding=self._execution,
            adapter_version=self.metadata.version,
        )

    def capture_build_identity(self):
        """Capture the composed source/engine/plugin/container build identity."""
        if self._execution is None:
            return _unsupported("capture_build_identity")
        return capture_bulk_rnaseq_build_identity(
            binding=self._execution,
            adapter_version=self.metadata.version,
        )

    def doctor(self):
        """Return redacted runtime availability diagnostics."""
        if self._execution is None:
            return _unsupported("doctor")
        return doctor_bulk_rnaseq_runtime(self._execution)

    def extract_artifacts(
        self,
        inputs: WorkflowInputs,
        workspace: str | Path,
    ) -> Result[tuple[ExtractedArtifactCandidate, ...]]:
        """Remain unsupported until deterministic output contracts are fixed."""
        return _unsupported("extract_artifacts")


def _unsupported(method: str) -> Result[Any]:
    return Result.failure(
        [
            Issue(
                code="BULK_RNASEQ_CAPABILITY_UNSUPPORTED",
                message="The bulk RNA-seq capability is not implemented.",
                severity="error",
                path=method,
                source="adapter",
                context={"method": method},
            )
        ]
    )
