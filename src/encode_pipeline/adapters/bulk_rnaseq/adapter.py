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
from encode_pipeline.platform.adapters import (
    INPUT_AUTHORING_CAPABILITY,
    VALIDATION_CAPABILITY,
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
    """Expose authoring and validation without claiming runtime readiness."""

    metadata = WorkflowMetadata(
        workflow_id=WORKFLOW_ID,
        name="Bulk RNA-seq",
        version=__version__ or "0.0.0",
        description=(
            "Contract-only HelixWeave adapter for pinned nf-core/rnaseq 3.26.0."
        ),
        engines=("nextflow",),
        tags=("bulk-rnaseq", "illumina", "nf-core"),
    )
    capabilities = WorkflowCapabilities(
        supports=(VALIDATION_CAPABILITY, INPUT_AUTHORING_CAPABILITY)
    )

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
        """Remain unsupported until runtime assets and layout are fixed."""
        return _unsupported("plan_workspace")

    def build_command(self, plan: WorkspacePlan) -> Result[CommandSpec]:
        """Remain unsupported until the Nextflow command contract is fixed."""
        return _unsupported("build_command")

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
