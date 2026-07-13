"""Pytest-independent conformance assertions for workflow adapters."""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
from typing import Callable, Never, TypeVar, cast

from encode_pipeline.platform.adapters import (
    ARTIFACT_EXTRACT_CAPABILITY,
    COMMAND_CAPABILITY,
    DAG_PREVIEW_CAPABILITY,
    INPUT_AUTHORING_CAPABILITY,
    QC_SUMMARY_EXTRACT_CAPABILITY,
    VALIDATION_CAPABILITY,
    WORKSPACE_PLAN_CAPABILITY,
    CommandSpec,
    DagPreview,
    ExtractedArtifactCandidate,
    ExtractedQcMetricCandidate,
    QcSourceDocument,
    QcSummaryExtractingAdapter,
    WorkflowAdapter,
    WorkflowCapabilities,
    WorkflowInputs,
    WorkflowMetadata,
    WorkflowSchema,
    WorkspacePlan,
)
from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.platform.results import Result


_T = TypeVar("_T")


class AdapterConformanceError(AssertionError):
    """Raised when an adapter violates a public contract coordinate."""


@dataclass(frozen=True)
class AdapterConformanceCase:
    """Adapter-owned fixtures needed to exercise the neutral contract."""

    adapter: WorkflowAdapter
    valid_inputs: WorkflowInputs
    invalid_inputs: WorkflowInputs
    planning_workspace: Path
    artifact_workspace: Path
    qc_sources: tuple[QcSourceDocument, ...] = ()

    def __post_init__(self) -> None:
        if not isinstance(self.valid_inputs, WorkflowInputs):
            raise ValueError("valid_inputs must be WorkflowInputs")
        if not isinstance(self.invalid_inputs, WorkflowInputs):
            raise ValueError("invalid_inputs must be WorkflowInputs")
        for name in ("planning_workspace", "artifact_workspace"):
            value = getattr(self, name)
            if not isinstance(value, Path) or not value.is_absolute():
                raise ValueError(f"{name} must be an absolute Path")
        if not isinstance(self.qc_sources, tuple) or not all(
            isinstance(source, QcSourceDocument) for source in self.qc_sources
        ):
            raise ValueError("qc_sources must be a tuple of QcSourceDocument")


def verify_adapter_conformance(case: AdapterConformanceCase) -> None:
    """Raise a bounded error when one adapter violates the current contract."""
    if not isinstance(case, AdapterConformanceCase):
        _fail("case", "must use AdapterConformanceCase")
    adapter = case.adapter
    if not _invoke(
        "adapter.protocol",
        lambda: isinstance(adapter, WorkflowAdapter),
    ):
        _fail("adapter.protocol", "does not satisfy WorkflowAdapter")
    metadata = _invoke("adapter.metadata", lambda: adapter.metadata)
    capabilities = _invoke("adapter.capabilities", lambda: adapter.capabilities)
    if not isinstance(metadata, WorkflowMetadata):
        _fail("adapter.metadata", "must use WorkflowMetadata")
    if not isinstance(capabilities, WorkflowCapabilities):
        _fail("adapter.capabilities", "must use WorkflowCapabilities")
    workflow_id = _invoke(
        "adapter.metadata.workflow_id",
        lambda: metadata.workflow_id,
    )
    engines = _invoke("adapter.metadata.engines", lambda: metadata.engines)
    tags = _invoke("adapter.metadata.tags", lambda: metadata.tags)
    supports = _invoke(
        "adapter.capabilities.supports",
        lambda: capabilities.supports,
    )

    registry = _invoke(
        "adapter.registry",
        lambda: WorkflowRegistry((adapter,)),
        failure_message="rejected the adapter declaration",
    )
    registered = _invoke(
        "adapter.registry",
        lambda: registry.get(workflow_id),
        failure_message="rejected the adapter declaration",
    )
    if registered is not adapter:
        _fail("adapter.registry", "did not preserve adapter identity")

    if not _invoke(
        "adapter.metadata.engines",
        lambda: len(engines) == len(set(engines)),
    ):
        _fail("adapter.metadata.engines", "must be unique")
    if not _invoke(
        "adapter.metadata.tags",
        lambda: len(tags) == len(set(tags)),
    ):
        _fail("adapter.metadata.tags", "must be unique")

    first_schema = _invoke("adapter.schema", lambda: adapter.schema())
    second_schema = _invoke("adapter.schema", lambda: adapter.schema())
    if not isinstance(first_schema, WorkflowSchema) or not isinstance(
        second_schema, WorkflowSchema
    ):
        _fail("adapter.schema", "must return WorkflowSchema")
    if first_schema is second_schema:
        _fail("adapter.schema", "must return a fresh contract instance")
    first_document = _invoke("adapter.schema", lambda: first_schema.to_dict())
    second_document = _invoke("adapter.schema", lambda: second_schema.to_dict())
    if not isinstance(first_document, dict) or not isinstance(second_document, dict):
        _fail("adapter.schema", "must return JSON-safe data")
    first_serialized = _invoke(
        "adapter.schema",
        lambda: json.dumps(first_document, allow_nan=False, sort_keys=True),
        failure_message="must return JSON-safe data",
    )
    second_serialized = _invoke(
        "adapter.schema",
        lambda: json.dumps(second_document, allow_nan=False, sort_keys=True),
        failure_message="must return JSON-safe data",
    )
    if first_serialized != second_serialized:
        _fail("adapter.schema", "must be stable across calls")
    _verify_json_schema(first_document)
    if INPUT_AUTHORING_CAPABILITY not in supports:
        _fail(
            "capability.input_authoring",
            "must be declared for the renderable schema contract",
        )

    _verify_validation(adapter, supports, case)
    workspace_plan = _verify_workspace(adapter, supports, case)
    _verify_dag(adapter, supports, case.valid_inputs)
    _verify_command(adapter, supports, workspace_plan)
    _verify_artifacts(adapter, supports, case)
    _verify_qc(adapter, supports, case)


def _verify_validation(
    adapter: WorkflowAdapter,
    supports: tuple[str, ...],
    case: AdapterConformanceCase,
) -> None:
    supported = VALIDATION_CAPABILITY in supports
    valid_result = _result(
        "capability.validation",
        lambda: adapter.validate(case.valid_inputs),
    )
    invalid_result = _result(
        "capability.validation",
        lambda: adapter.validate(case.invalid_inputs),
    )
    if supported:
        _require_success(valid_result, "capability.validation")
        _require_failure(invalid_result, "capability.validation.invalid_inputs")
    else:
        _require_failure(valid_result, "capability.validation")
        _require_failure(invalid_result, "capability.validation")


def _verify_workspace(
    adapter: WorkflowAdapter,
    supports: tuple[str, ...],
    case: AdapterConformanceCase,
) -> WorkspacePlan:
    supported = WORKSPACE_PLAN_CAPABILITY in supports
    if _invoke(
        "capability.workspace_plan.fixture",
        lambda: case.planning_workspace.exists(),
    ):
        _fail("capability.workspace_plan.fixture", "must start absent")
    result = _result(
        "capability.workspace_plan",
        lambda: adapter.plan_workspace(
            case.valid_inputs,
            case.planning_workspace,
        ),
    )
    if _invoke(
        "capability.workspace_plan",
        lambda: case.planning_workspace.exists(),
    ):
        _fail("capability.workspace_plan", "must not materialize the workspace")
    if not supported:
        _require_failure(result, "capability.workspace_plan")
        return WorkspacePlan()
    _require_success(result, "capability.workspace_plan")
    value = _result_value(result, "capability.workspace_plan")
    if not isinstance(value, WorkspacePlan):
        _fail("capability.workspace_plan", "must return WorkspacePlan")
    return value


def _verify_dag(
    adapter: WorkflowAdapter,
    supports: tuple[str, ...],
    inputs: WorkflowInputs,
) -> None:
    supported = DAG_PREVIEW_CAPABILITY in supports
    result = _result("capability.dag_preview", lambda: adapter.preview_dag(inputs))
    if not supported:
        _require_failure(result, "capability.dag_preview")
        return
    _require_success(result, "capability.dag_preview")
    if not isinstance(
        _result_value(result, "capability.dag_preview"),
        DagPreview,
    ):
        _fail("capability.dag_preview", "must return DagPreview")


def _verify_command(
    adapter: WorkflowAdapter,
    supports: tuple[str, ...],
    plan: WorkspacePlan,
) -> None:
    supported = COMMAND_CAPABILITY in supports
    result = _result("capability.command", lambda: adapter.build_command(plan))
    if not supported:
        _require_failure(result, "capability.command")
        return
    _require_success(result, "capability.command")
    if not isinstance(_result_value(result, "capability.command"), CommandSpec):
        _fail("capability.command", "must return CommandSpec")


def _verify_artifacts(
    adapter: WorkflowAdapter,
    supports: tuple[str, ...],
    case: AdapterConformanceCase,
) -> None:
    supported = ARTIFACT_EXTRACT_CAPABILITY in supports
    result = _result(
        "capability.artifact_extract",
        lambda: adapter.extract_artifacts(
            case.valid_inputs,
            case.artifact_workspace,
        ),
    )
    if not supported:
        _require_failure(result, "capability.artifact_extract")
        return
    _require_success(result, "capability.artifact_extract")
    value = _result_value(result, "capability.artifact_extract")
    candidates_valid = _invoke(
        "capability.artifact_extract",
        lambda: (
            isinstance(value, tuple)
            and all(
                isinstance(candidate, ExtractedArtifactCandidate) for candidate in value
            )
        ),
    )
    if not candidates_valid:
        _fail(
            "capability.artifact_extract",
            "must return a tuple of ExtractedArtifactCandidate",
        )


def _verify_qc(
    adapter: WorkflowAdapter,
    supports: tuple[str, ...],
    case: AdapterConformanceCase,
) -> None:
    supported = QC_SUMMARY_EXTRACT_CAPABILITY in supports
    implements = _invoke(
        "capability.qc_summary_extract",
        lambda: isinstance(adapter, QcSummaryExtractingAdapter),
    )
    if supported != implements:
        _fail(
            "capability.qc_summary_extract",
            "declaration and optional protocol must agree",
        )
    if not supported:
        return
    qc_adapter = cast(QcSummaryExtractingAdapter, adapter)
    source_types = _invoke(
        "capability.qc_summary_extract.source_types",
        lambda: qc_adapter.qc_source_output_types(),
    )
    source_types_valid = _invoke(
        "capability.qc_summary_extract.source_types",
        lambda: (
            isinstance(source_types, tuple)
            and bool(source_types)
            and all(_safe_output_type(value) for value in source_types)
            and len(source_types) == len(set(source_types))
        ),
    )
    if not source_types_valid:
        _fail(
            "capability.qc_summary_extract.source_types",
            "must return unique non-empty string names",
        )
    result = _result(
        "capability.qc_summary_extract",
        lambda: qc_adapter.extract_qc_metrics(case.valid_inputs, case.qc_sources),
    )
    _require_success(result, "capability.qc_summary_extract")
    value = _result_value(result, "capability.qc_summary_extract")
    candidates_valid = _invoke(
        "capability.qc_summary_extract",
        lambda: (
            isinstance(value, tuple)
            and all(
                isinstance(candidate, ExtractedQcMetricCandidate) for candidate in value
            )
        ),
    )
    if not candidates_valid:
        _fail(
            "capability.qc_summary_extract",
            "must return a tuple of ExtractedQcMetricCandidate",
        )


def _result(coordinate: str, callback: Callable[[], object]) -> Result:
    value = _invoke(coordinate, callback)
    if not isinstance(value, Result):
        _fail(coordinate, "must return Result")
    return value


def _verify_json_schema(document: dict[str, object]) -> None:
    validator = _invoke(
        "adapter.schema",
        _load_json_schema_validator,
        failure_message="must satisfy JSON Schema 2020-12",
    )
    for name in ("config_schema", "sample_schema", "option_schema"):
        value = _invoke(
            f"adapter.schema.{name}",
            lambda name=name: document.get(name),
        )
        if not isinstance(value, dict):
            _fail(f"adapter.schema.{name}", "must be a JSON Schema object")
        _invoke(
            "adapter.schema",
            lambda value=value: validator.check_schema(value),
            failure_message="must satisfy JSON Schema 2020-12",
        )


def _load_json_schema_validator():
    from jsonschema import Draft202012Validator

    return Draft202012Validator


def _safe_output_type(value: object) -> bool:
    if not isinstance(value, str) or not 1 <= len(value) <= 128:
        return False
    if not value[0].isalpha() or not value[0].isascii():
        return False
    return all(
        character.isascii() and (character.isalnum() or character in {"_", ".", "-"})
        for character in value
    )


def _require_success(result: Result, coordinate: str) -> None:
    if _invoke(coordinate, lambda: result.is_failure):
        _fail(coordinate, "declared support returned failure")


def _require_failure(result: Result, coordinate: str) -> None:
    if _invoke(coordinate, lambda: result.is_success):
        _fail(coordinate, "unsupported or invalid input returned success")


def _result_value(result: Result, coordinate: str) -> object:
    return _invoke(coordinate, lambda: result.value)


def _capture(callback: Callable[[], _T]) -> tuple[bool, _T | None]:
    try:
        value = callback()
    except Exception:
        return False, None
    return True, value


def _invoke(
    coordinate: str,
    callback: Callable[[], _T],
    *,
    failure_message: str = "raised an exception",
) -> _T:
    succeeded, value = _capture(callback)
    if not succeeded:
        _fail(coordinate, failure_message)
    return cast(_T, value)


def _fail(coordinate: str, message: str) -> Never:
    raise AdapterConformanceError(f"{coordinate}: {message}") from None
