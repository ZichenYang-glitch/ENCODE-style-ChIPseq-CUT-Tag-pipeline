"""Composable adapters for pinned nf-core/rnaseq input and result semantics."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from encode_pipeline import __version__
from encode_pipeline.adapters.bulk_rnaseq.authoring import (
    build_bulk_rnaseq_authoring_schema,
)
from encode_pipeline.adapters.bulk_rnaseq.artifacts import (
    discover_bulk_rnaseq_artifacts,
)
from encode_pipeline.adapters.bulk_rnaseq.qc import (
    BULK_RNASEQ_QC_SOURCE_TYPES,
    extract_bulk_rnaseq_qc_metrics,
)
from encode_pipeline.adapters.bulk_rnaseq.qualification import (
    BulkRnaSeqExecutionMode,
)
from encode_pipeline.adapters.bulk_rnaseq.results_contract import (
    load_bulk_rnaseq_results_contract,
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
from encode_pipeline.adapters.bulk_rnaseq.upstream import (
    NFCORE_RNASEQ_COMMIT,
    NFCORE_RNASEQ_RELEASE,
)
from encode_pipeline.platform.adapters import (
    ARTIFACT_EXTRACT_CAPABILITY,
    COMMAND_CAPABILITY,
    INPUT_AUTHORING_CAPABILITY,
    QC_SUMMARY_EXTRACT_CAPABILITY,
    VALIDATION_CAPABILITY,
    WORKSPACE_PLAN_CAPABILITY,
    CommandSpec,
    DagPreview,
    ExtractedArtifactCandidate,
    ExtractedQcMetricCandidate,
    QcSourceDocument,
    WorkflowCapabilities,
    WorkflowAvailability,
    WorkflowInputs,
    WorkflowMetadata,
    WorkflowSchema,
    WorkflowUpstreamIdentity,
    WorkspacePlan,
)
from encode_pipeline.platform.results import Issue, Result


WORKFLOW_ID = "bulk-rnaseq"


class BulkRnaSeqWorkflowAdapter:
    """Expose contract-only or explicitly composed offline runtime behavior."""

    _execution_variant = "runtime-v1"
    _execution_mode = BulkRnaSeqExecutionMode.STANDARD

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
    upstream_identity = WorkflowUpstreamIdentity(
        name="nf-core/rnaseq",
        version=NFCORE_RNASEQ_RELEASE,
        revision=NFCORE_RNASEQ_COMMIT,
    )

    def __init__(
        self,
        *,
        execution: BulkRnaSeqExecutionBinding | None = None,
        configured_availability: WorkflowAvailability | None = None,
    ) -> None:
        if execution is not None and not isinstance(
            execution, BulkRnaSeqExecutionBinding
        ):
            raise ValueError("execution must be a BulkRnaSeqExecutionBinding or None")
        if configured_availability is None:
            configured_availability = WorkflowAvailability(
                execution=("available" if execution is not None else "not_configured"),
                reason_code=(
                    "WORKFLOW_EXECUTION_READY"
                    if execution is not None
                    else "WORKFLOW_EXECUTION_NOT_CONFIGURED"
                ),
            )
        if not isinstance(configured_availability, WorkflowAvailability):
            raise ValueError("configured_availability must be WorkflowAvailability")
        if execution is None and configured_availability.execution == "available":
            raise ValueError("execution availability requires an execution binding")
        self._execution = execution
        self._configured_availability = configured_availability
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
            adapter_variant=self._execution_variant,
            execution_mode=self._execution_mode,
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
            adapter_variant=self._execution_variant,
            execution_mode=self._execution_mode,
        )

    def capture_build_identity(self):
        """Capture the composed source/engine/plugin/container build identity."""
        if self._execution is None:
            return _unsupported("capture_build_identity")
        return capture_bulk_rnaseq_build_identity(
            binding=self._execution,
            adapter_version=self.metadata.version,
            adapter_variant=self._execution_variant,
            execution_mode=self._execution_mode,
        )

    def doctor(self):
        """Return redacted runtime availability diagnostics."""
        if self._execution is None:
            return _unsupported("doctor")
        return doctor_bulk_rnaseq_runtime(self._execution)

    def execution_availability(self) -> WorkflowAvailability:
        """Return path-free readiness while reusing the live admission authority."""
        if self._execution is None:
            return self._configured_availability
        try:
            report = doctor_bulk_rnaseq_runtime(self._execution)
        except Exception:
            return WorkflowAvailability(
                execution="unavailable",
                reason_code="WORKFLOW_EXECUTION_UNAVAILABLE",
            )
        if report.ready:
            return WorkflowAvailability()
        return WorkflowAvailability(
            execution="unavailable",
            reason_code="WORKFLOW_EXECUTION_UNAVAILABLE",
        )

    @property
    def execution_binding(self) -> BulkRnaSeqExecutionBinding | None:
        """Return the deployment-owned binding for local service composition."""
        return self._execution

    def disable_execution(self) -> None:
        """Permanently fail closed when local runner composition is unavailable."""
        self._execution = None
        self._configured_availability = WorkflowAvailability(
            execution="unavailable",
            reason_code="WORKFLOW_EXECUTION_UNAVAILABLE",
        )

    def extract_artifacts(
        self,
        inputs: WorkflowInputs,
        workspace: str | Path,
    ) -> Result[tuple[ExtractedArtifactCandidate, ...]]:
        """Remain unsupported until deterministic output contracts are fixed."""
        return _unsupported("extract_artifacts")


class BulkRnaSeqRapidQuantQualificationAdapter(BulkRnaSeqWorkflowAdapter):
    """Run only the pinned upstream rapid_quant qualification route."""

    _execution_variant = "rapid-quant-qualification-v1"
    _execution_mode = BulkRnaSeqExecutionMode.RAPID_QUANT

    def __init__(self, *, execution: BulkRnaSeqExecutionBinding) -> None:
        if not isinstance(execution, BulkRnaSeqExecutionBinding):
            raise ValueError("execution must be a BulkRnaSeqExecutionBinding")
        super().__init__(execution=execution)


class BulkRnaSeqResultsWorkflowAdapter(BulkRnaSeqWorkflowAdapter):
    """Add the pinned STAR+Salmon artifact and machine-QC result boundary."""

    _execution_variant = "results-v1"

    def __init__(self, *, execution: BulkRnaSeqExecutionBinding) -> None:
        if not isinstance(execution, BulkRnaSeqExecutionBinding):
            raise ValueError("execution must be a BulkRnaSeqExecutionBinding")
        load_bulk_rnaseq_results_contract()
        super().__init__(execution=execution)
        self.capabilities = WorkflowCapabilities(
            supports=(
                *self.capabilities.supports,
                ARTIFACT_EXTRACT_CAPABILITY,
                QC_SUMMARY_EXTRACT_CAPABILITY,
            )
        )

    def extract_artifacts(
        self,
        inputs: WorkflowInputs,
        workspace: str | Path,
    ) -> Result[tuple[ExtractedArtifactCandidate, ...]]:
        """Discover only fixed, allowlisted nf-core/rnaseq 3.26.0 outputs."""
        return discover_bulk_rnaseq_artifacts(inputs, workspace)

    def qc_source_output_types(self) -> tuple[str, ...]:
        """Return the closed machine-readable QC source catalog."""
        return BULK_RNASEQ_QC_SOURCE_TYPES

    def extract_qc_metrics(
        self,
        inputs: WorkflowInputs,
        sources: tuple[QcSourceDocument, ...],
    ) -> Result[tuple[ExtractedQcMetricCandidate, ...]]:
        """Extract deterministic Decimal metrics from vetted source bytes."""
        return extract_bulk_rnaseq_qc_metrics(inputs, sources)


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
