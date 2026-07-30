"""Workflow-platform adapter registry primitives."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from types import MappingProxyType

from encode_pipeline.platform.adapters import (
    ARTIFACT_EXTRACT_CAPABILITY,
    COMMAND_CAPABILITY,
    EXECUTION_CAPABILITY_NAMES,
    INPUT_BUNDLE_IMPORT_CAPABILITY,
    QC_SUMMARY_EXTRACT_CAPABILITY,
    InputBundleImportingAdapter,
    QcSummaryExtractingAdapter,
    WorkflowAdapter,
    WorkflowAvailabilityProvidingAdapter,
    WorkflowBuildIdentityProvidingAdapter,
    WorkflowCapabilities,
    WorkflowMetadata,
    WORKSPACE_PLAN_CAPABILITY,
)

ENCODE_WORKFLOW_ID = "encode-style-chipseq-cuttag-atac-mnase"
_ENCODE_EXECUTION_CAPABILITIES = frozenset(
    {
        "validation",
        "workspace_plan",
        "input_authoring",
        "input_bundle_import",
        "artifact_extract",
        "qc_summary_extract",
    }
)


class WorkflowRegistry:
    """Immutable adapter lookup by workflow ID."""

    __slots__ = (
        "__adapters_by_workflow_id",
        "__legacy_execution_fallbacks",
        "__metadata",
    )

    def __init__(
        self,
        adapters: Iterable[WorkflowAdapter] = (),
        *,
        legacy_execution_fallbacks: Iterable[WorkflowAdapter] = (),
    ) -> None:
        adapter_sequence = tuple(adapters)
        fallback_sequence = tuple(legacy_execution_fallbacks)
        if len({id(adapter) for adapter in fallback_sequence}) != len(
            fallback_sequence
        ):
            raise ValueError(
                "WorkflowRegistry legacy execution fallbacks must be unique"
            )
        if any(
            not any(fallback is adapter for adapter in adapter_sequence)
            for fallback in fallback_sequence
        ):
            raise ValueError(
                "WorkflowRegistry legacy execution fallbacks must be registered "
                "adapter instances"
            )
        adapters_by_workflow_id: dict[str, WorkflowAdapter] = {}
        metadata: list[WorkflowMetadata] = []

        for adapter in adapter_sequence:
            if not isinstance(adapter, WorkflowAdapter):
                capabilities = getattr(adapter, "capabilities", None)
                if (
                    isinstance(capabilities, WorkflowCapabilities)
                    and ARTIFACT_EXTRACT_CAPABILITY in capabilities.supports
                    and not callable(getattr(adapter, "extract_artifacts", None))
                ):
                    raise ValueError(
                        "WorkflowRegistry artifact_extract capability requires "
                        "an extractor"
                    )
                raise ValueError(
                    "WorkflowRegistry adapters must satisfy WorkflowAdapter"
                )
            if not isinstance(adapter.metadata, WorkflowMetadata):
                raise ValueError(
                    "WorkflowRegistry adapter metadata must be WorkflowMetadata"
                )
            if not isinstance(adapter.capabilities, WorkflowCapabilities):
                raise ValueError(
                    "WorkflowRegistry adapter capabilities must be WorkflowCapabilities"
                )
            declares_qc = QC_SUMMARY_EXTRACT_CAPABILITY in adapter.capabilities.supports
            implements_qc = isinstance(adapter, QcSummaryExtractingAdapter)
            if declares_qc != implements_qc:
                raise ValueError(
                    "WorkflowRegistry qc_summary_extract capability and "
                    "protocol must agree"
                )
            declares_bundle_import = (
                INPUT_BUNDLE_IMPORT_CAPABILITY in adapter.capabilities.supports
            )
            implements_bundle_import = isinstance(adapter, InputBundleImportingAdapter)
            if declares_bundle_import != implements_bundle_import:
                raise ValueError(
                    "WorkflowRegistry input_bundle_import capability and "
                    "protocol must agree"
                )
            _validate_execution_declaration(
                adapter,
                uses_encode_execution_fallback=any(
                    adapter is fallback for fallback in fallback_sequence
                ),
            )

            workflow_id = adapter.metadata.workflow_id
            if workflow_id in adapters_by_workflow_id:
                raise ValueError(f"Duplicate workflow_id: {workflow_id}")

            adapters_by_workflow_id[workflow_id] = adapter
            metadata.append(adapter.metadata)

        self.__adapters_by_workflow_id: Mapping[str, WorkflowAdapter] = (
            MappingProxyType(adapters_by_workflow_id)
        )
        self.__legacy_execution_fallbacks = fallback_sequence
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

    def uses_encode_execution_fallback(self, adapter: WorkflowAdapter) -> bool:
        """Return whether composition trusted this exact legacy ENCODE instance."""
        return any(
            adapter is fallback for fallback in self.__legacy_execution_fallbacks
        ) and _matches_encode_execution_fallback(adapter)


def _normalize_workflow_id(workflow_id: str) -> str:
    if not isinstance(workflow_id, str):
        raise ValueError("workflow_id must be a string")
    normalized = workflow_id.strip()
    if not normalized:
        raise ValueError("workflow_id must be non-empty")
    return normalized


def declares_workflow_execution(adapter: WorkflowAdapter) -> bool:
    """Return whether an adapter declares any execution-owned capability."""
    return bool(EXECUTION_CAPABILITY_NAMES.intersection(adapter.capabilities.supports))


def _matches_encode_execution_fallback(adapter: WorkflowAdapter) -> bool:
    """Return whether an explicitly trusted adapter has legacy ENCODE identity."""
    try:
        return (
            adapter.metadata.workflow_id == ENCODE_WORKFLOW_ID
            and adapter.metadata.engines == ("snakemake",)
            and frozenset(adapter.capabilities.supports)
            == _ENCODE_EXECUTION_CAPABILITIES
            and not isinstance(adapter, WorkflowBuildIdentityProvidingAdapter)
            and not isinstance(adapter, WorkflowAvailabilityProvidingAdapter)
        )
    except Exception:
        return False


def _validate_execution_declaration(
    adapter: WorkflowAdapter,
    *,
    uses_encode_execution_fallback: bool,
) -> None:
    supports = frozenset(adapter.capabilities.supports)
    if ARTIFACT_EXTRACT_CAPABILITY in supports and not callable(
        getattr(adapter, "extract_artifacts", None)
    ):
        raise ValueError(
            "WorkflowRegistry artifact_extract capability requires an extractor"
        )
    if not declares_workflow_execution(adapter):
        return
    if uses_encode_execution_fallback and _matches_encode_execution_fallback(adapter):
        return
    declares_workspace = WORKSPACE_PLAN_CAPABILITY in supports
    declares_command = COMMAND_CAPABILITY in supports
    if not declares_workspace or not declares_command:
        raise ValueError(
            "WorkflowRegistry execution requires both workspace_plan and command"
        )
    if not isinstance(adapter, WorkflowAvailabilityProvidingAdapter):
        raise ValueError("WorkflowRegistry execution requires an availability provider")
    if not isinstance(adapter, WorkflowBuildIdentityProvidingAdapter):
        raise ValueError(
            "WorkflowRegistry execution requires an adapter build identity provider"
        )
