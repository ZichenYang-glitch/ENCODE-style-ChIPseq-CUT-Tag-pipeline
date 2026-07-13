"""Workflow-platform adapter contract primitives."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from copy import deepcopy
from dataclasses import dataclass, field
from decimal import Decimal
import json
import math
from pathlib import Path
import re
from typing import TYPE_CHECKING, Any, Protocol, runtime_checkable

from encode_pipeline.platform.results import Result

if TYPE_CHECKING:
    from encode_pipeline.platform.planning import ExecutionPlan
    from encode_pipeline.platform.runs import RunRecord


SamplePayload = str | Path | list[dict[str, str]] | None

JSON_SCHEMA_DIALECT = "https://json-schema.org/draft/2020-12/schema"
MAX_AUTHORING_REQUEST_BYTES = 2 * 1024 * 1024
MAX_SAMPLE_ROWS = 1_000
MAX_SAMPLE_COLUMNS = 64
MAX_SAMPLE_COLUMN_NAME_LENGTH = 128
MAX_SAMPLE_CELL_LENGTH = 4_096

_SCHEMA_VERSION_PATTERN = re.compile(r"^[1-9]\d*\.(?:0|[1-9]\d*)\.(?:0|[1-9]\d*)$")
_MODE_PATTERN = re.compile(r"^[a-z][a-z0-9_]{0,63}$")
_SCHEMA_COVERAGE_VALUES = frozenset({"partial", "complete"})

VALIDATION_CAPABILITY = "validation"
DAG_PREVIEW_CAPABILITY = "dag_preview"
WORKSPACE_PLAN_CAPABILITY = "workspace_plan"
COMMAND_CAPABILITY = "command"
INPUT_AUTHORING_CAPABILITY = "input_authoring"
ARTIFACT_EXTRACT_CAPABILITY = "artifact_extract"
QC_SUMMARY_EXTRACT_CAPABILITY = "qc_summary_extract"
WORKFLOW_CAPABILITY_NAMES = frozenset(
    {
        VALIDATION_CAPABILITY,
        DAG_PREVIEW_CAPABILITY,
        WORKSPACE_PLAN_CAPABILITY,
        COMMAND_CAPABILITY,
        INPUT_AUTHORING_CAPABILITY,
        ARTIFACT_EXTRACT_CAPABILITY,
        QC_SUMMARY_EXTRACT_CAPABILITY,
    }
)


@dataclass(frozen=True)
class WorkflowSchemaCoverage:
    """Coverage of each canonical adapter-owned authoring surface.

    Complete means every canonical field is representable. It does not claim
    to enumerate every legacy spelling that a scientific validator accepts.
    """

    config: str = "partial"
    samples: str = "partial"
    options: str = "partial"

    def __post_init__(self) -> None:
        for name in ("config", "samples", "options"):
            value = getattr(self, name)
            if value not in _SCHEMA_COVERAGE_VALUES:
                raise ValueError(f"Workflow schema {name} coverage is invalid")

    def to_dict(self) -> dict[str, str]:
        """Return a JSON-ready coverage mapping."""
        return {
            "config": self.config,
            "samples": self.samples,
            "options": self.options,
        }


@dataclass(frozen=True)
class WorkflowAuthoringModes:
    """Adapter-declared editor modes for config, samples, and options."""

    config: tuple[str, ...] = ()
    samples: tuple[str, ...] = ()
    options: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        for name in ("config", "samples", "options"):
            object.__setattr__(
                self,
                name,
                _normalize_mode_tuple(getattr(self, name), f"authoring_modes.{name}"),
            )

    def to_dict(self) -> dict[str, list[str]]:
        """Return fresh JSON-ready mode lists."""
        return {
            "config": list(self.config),
            "samples": list(self.samples),
            "options": list(self.options),
        }


@dataclass(frozen=True)
class WorkflowInputModes:
    """Wire-level input representations accepted by an adapter."""

    config: tuple[str, ...] = ("object",)
    samples: tuple[str, ...] = ("server_path",)
    options: tuple[str, ...] = ("object",)

    def __post_init__(self) -> None:
        for name in ("config", "samples", "options"):
            object.__setattr__(
                self,
                name,
                _normalize_mode_tuple(getattr(self, name), f"input_modes.{name}"),
            )

    def to_dict(self) -> dict[str, list[str]]:
        """Return fresh JSON-ready mode lists."""
        return {
            "config": list(self.config),
            "samples": list(self.samples),
            "options": list(self.options),
        }


@dataclass(frozen=True)
class WorkflowInputLimits:
    """Public projection of the platform-wide authoring ceilings.

    Contract version 1.0.0 only permits the exact platform values because the
    transport and domain boundaries enforce one uniform set of ceilings.
    """

    max_request_bytes: int = MAX_AUTHORING_REQUEST_BYTES
    max_sample_rows: int = MAX_SAMPLE_ROWS
    max_sample_columns: int = MAX_SAMPLE_COLUMNS
    max_sample_column_name_length: int = MAX_SAMPLE_COLUMN_NAME_LENGTH
    max_sample_cell_length: int = MAX_SAMPLE_CELL_LENGTH

    def __post_init__(self) -> None:
        ceilings = {
            "max_request_bytes": MAX_AUTHORING_REQUEST_BYTES,
            "max_sample_rows": MAX_SAMPLE_ROWS,
            "max_sample_columns": MAX_SAMPLE_COLUMNS,
            "max_sample_column_name_length": MAX_SAMPLE_COLUMN_NAME_LENGTH,
            "max_sample_cell_length": MAX_SAMPLE_CELL_LENGTH,
        }
        for name, ceiling in ceilings.items():
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int):
                raise ValueError(f"Workflow input limit {name} must be an integer")
            if value != ceiling:
                raise ValueError(
                    f"Workflow input limit {name} must equal the platform "
                    f"ceiling {ceiling}"
                )

    def to_dict(self) -> dict[str, int]:
        """Return a JSON-ready limit mapping."""
        return {
            "max_request_bytes": self.max_request_bytes,
            "max_sample_rows": self.max_sample_rows,
            "max_sample_columns": self.max_sample_columns,
            "max_sample_column_name_length": self.max_sample_column_name_length,
            "max_sample_cell_length": self.max_sample_cell_length,
        }


@dataclass(frozen=True)
class WorkflowMetadata:
    """Stable identity and description for one workflow adapter."""

    workflow_id: str
    name: str
    version: str
    description: str | None = None
    engines: tuple[str, ...] = ("snakemake",)
    tags: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "workflow_id",
            _normalize_required_string(self.workflow_id, "workflow_id"),
        )
        object.__setattr__(
            self,
            "name",
            _normalize_required_string(self.name, "name"),
        )
        object.__setattr__(
            self,
            "version",
            _normalize_required_string(self.version, "version"),
        )
        if self.description is not None and not isinstance(self.description, str):
            raise ValueError("WorkflowMetadata description must be a string or None")
        object.__setattr__(
            self,
            "engines",
            _normalize_string_tuple(self.engines, "engines"),
        )
        object.__setattr__(
            self,
            "tags",
            _normalize_string_tuple(self.tags, "tags"),
        )

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-ready representation."""
        return {
            "workflow_id": self.workflow_id,
            "name": self.name,
            "version": self.version,
            "description": self.description,
            "engines": list(self.engines),
            "tags": list(self.tags),
        }


@dataclass(frozen=True)
class WorkflowCapabilities:
    """Canonical capability labels supported by a workflow adapter."""

    supports: tuple[str, ...]

    def __post_init__(self) -> None:
        supports = _normalize_mode_tuple(self.supports, "supports")
        unknown = set(supports).difference(WORKFLOW_CAPABILITY_NAMES)
        if unknown:
            raise ValueError("Workflow capabilities contain an unknown name")
        object.__setattr__(
            self,
            "supports",
            supports,
        )

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-ready representation."""
        return {"supports": list(self.supports)}


@dataclass(frozen=True)
class WorkflowSchema:
    """Top-level frozen authoring contract with copied JSON mappings.

    Constructor mappings are defensively deep-copied and ``to_dict`` returns
    fresh deep copies. The stored JSON Schema mappings remain ordinary mutable
    dictionaries; top-level dataclass fields alone are frozen.
    """

    schema_version: str = "1.0.0"
    schema_dialect: str = JSON_SCHEMA_DIALECT
    coverage: WorkflowSchemaCoverage = field(default_factory=WorkflowSchemaCoverage)
    authoring_modes: WorkflowAuthoringModes = field(
        default_factory=WorkflowAuthoringModes
    )
    input_modes: WorkflowInputModes = field(default_factory=WorkflowInputModes)
    limits: WorkflowInputLimits = field(default_factory=WorkflowInputLimits)
    config_schema: Mapping[str, object] = field(default_factory=dict)
    sample_schema: Mapping[str, object] = field(default_factory=dict)
    option_schema: Mapping[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if (
            not isinstance(self.schema_version, str)
            or _SCHEMA_VERSION_PATTERN.fullmatch(self.schema_version) is None
        ):
            raise ValueError("WorkflowSchema schema_version must be semantic version")
        if self.schema_dialect != JSON_SCHEMA_DIALECT:
            raise ValueError("WorkflowSchema schema_dialect is not supported")
        for name, expected_type in (
            ("coverage", WorkflowSchemaCoverage),
            ("authoring_modes", WorkflowAuthoringModes),
            ("input_modes", WorkflowInputModes),
            ("limits", WorkflowInputLimits),
        ):
            if not isinstance(getattr(self, name), expected_type):
                raise ValueError(
                    f"WorkflowSchema {name} must be {expected_type.__name__}"
                )
        object.__setattr__(
            self,
            "config_schema",
            _copy_json_schema(
                self.config_schema,
                "config_schema",
                dialect=self.schema_dialect,
                root_type="object",
            ),
        )
        object.__setattr__(
            self,
            "sample_schema",
            _copy_json_schema(
                self.sample_schema,
                "sample_schema",
                dialect=self.schema_dialect,
                root_type="array",
            ),
        )
        object.__setattr__(
            self,
            "option_schema",
            _copy_json_schema(
                self.option_schema,
                "option_schema",
                dialect=self.schema_dialect,
                root_type="object",
            ),
        )

    def to_dict(self) -> dict[str, Any]:
        """Return a fresh JSON-ready authoring contract."""
        return {
            "schema_version": self.schema_version,
            "schema_dialect": self.schema_dialect,
            "coverage": self.coverage.to_dict(),
            "authoring_modes": self.authoring_modes.to_dict(),
            "input_modes": self.input_modes.to_dict(),
            "limits": self.limits.to_dict(),
            "config_schema": deepcopy(self.config_schema),
            "sample_schema": deepcopy(self.sample_schema),
            "option_schema": deepcopy(self.option_schema),
        }


@dataclass(frozen=True)
class WorkflowInputs:
    """Submitted workflow inputs before adapter-specific normalization."""

    config: Mapping[str, object]
    samples: SamplePayload = None
    options: Mapping[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(self, "config", _copy_mapping(self.config, "config"))
        object.__setattr__(self, "samples", _copy_samples(self.samples))
        object.__setattr__(self, "options", _copy_mapping(self.options, "options"))

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-ready representation with fresh copies."""
        samples: str | list[dict[str, str]] | None
        if isinstance(self.samples, Path):
            samples = str(self.samples)
        elif isinstance(self.samples, list):
            samples = [dict(row) for row in self.samples]
        else:
            samples = self.samples
        return {
            "config": deepcopy(self.config),
            "samples": samples,
            "options": deepcopy(self.options),
        }


@dataclass(frozen=True)
class ExtractedArtifactCandidate:
    """Adapter-neutral logical description of one discovered output file."""

    output_type: str
    relative_path: str
    mime_type: str | None = None
    metadata: Mapping[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not isinstance(self.output_type, str):
            raise ValueError("output_type must be a string")
        if not isinstance(self.relative_path, str):
            raise ValueError("relative_path must be a string")
        if self.mime_type is not None and not isinstance(self.mime_type, str):
            raise ValueError("mime_type must be a string or None")
        object.__setattr__(
            self,
            "metadata",
            _copy_mapping(self.metadata, "metadata"),
        )


@dataclass(frozen=True)
class QcSourceArtifact:
    """Persisted artifact identity offered to an adapter as a QC source."""

    artifact_id: str
    output_type: str
    relative_path: str
    metadata: Mapping[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        for name in ("artifact_id", "output_type", "relative_path"):
            if not isinstance(getattr(self, name), str):
                raise ValueError(f"{name} must be a string")
        object.__setattr__(self, "metadata", _copy_mapping(self.metadata, "metadata"))


@dataclass(frozen=True)
class QcSourceDocument:
    """Platform-vetted bounded bytes for one persisted QC source artifact."""

    source: QcSourceArtifact
    content: bytes

    def __post_init__(self) -> None:
        if not isinstance(self.source, QcSourceArtifact):
            raise ValueError("source must be a QcSourceArtifact")
        if not isinstance(self.content, bytes):
            raise ValueError("content must be bytes")


@dataclass(frozen=True)
class ExtractedQcMetricCandidate:
    """Adapter-neutral numeric QC metric before durable identity assignment."""

    metric_key: str
    display_name: str
    value: Decimal
    unit: str
    scope: str
    source_artifact_id: str
    sample_id: str | None = None
    experiment_id: str | None = None
    assay: str | None = None
    qc_flag: str | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.value, Decimal):
            raise ValueError("value must be a Decimal")


@dataclass(frozen=True)
class DagNode:
    """Engine-neutral DAG preview node."""

    id: str
    label: str
    inputs: tuple[str, ...] = ()
    outputs: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        object.__setattr__(self, "id", _normalize_required_string(self.id, "id"))
        object.__setattr__(
            self,
            "label",
            _normalize_required_string(self.label, "label"),
        )
        object.__setattr__(
            self,
            "inputs",
            _normalize_string_tuple(self.inputs, "inputs"),
        )
        object.__setattr__(
            self,
            "outputs",
            _normalize_string_tuple(self.outputs, "outputs"),
        )

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-ready representation."""
        return {
            "id": self.id,
            "label": self.label,
            "inputs": list(self.inputs),
            "outputs": list(self.outputs),
        }


@dataclass(frozen=True)
class DagEdge:
    """Engine-neutral DAG preview edge."""

    source: str
    target: str

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "source",
            _normalize_required_string(self.source, "source"),
        )
        object.__setattr__(
            self,
            "target",
            _normalize_required_string(self.target, "target"),
        )

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-ready representation."""
        return {"source": self.source, "target": self.target}


@dataclass(frozen=True)
class DagPreview:
    """Engine-neutral DAG preview graph."""

    nodes: tuple[DagNode, ...] = ()
    edges: tuple[DagEdge, ...] = ()

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "nodes",
            _normalize_instance_tuple(self.nodes, DagNode, "nodes"),
        )
        object.__setattr__(
            self,
            "edges",
            _normalize_instance_tuple(self.edges, DagEdge, "edges"),
        )

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-ready representation."""
        return {
            "nodes": [node.to_dict() for node in self.nodes],
            "edges": [edge.to_dict() for edge in self.edges],
        }


@dataclass(frozen=True)
class WorkspacePlan:
    """A pure data description of workspace directories and file contents."""

    directories: tuple[str, ...] = ()
    files: tuple[tuple[str, bytes], ...] = ()

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "directories",
            _normalize_string_tuple(self.directories, "directories"),
        )
        object.__setattr__(self, "files", _normalize_files(self.files))


@dataclass(frozen=True)
class CommandSpec:
    """A pure data command description. It does not execute anything."""

    argv: tuple[str, ...]
    cwd: str | None = None
    env: Mapping[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        argv = _normalize_string_tuple(self.argv, "argv")
        if not argv:
            raise ValueError("CommandSpec argv must be non-empty")
        if self.cwd is not None and not isinstance(self.cwd, str):
            raise ValueError("CommandSpec cwd must be a string or None")
        object.__setattr__(self, "argv", argv)
        object.__setattr__(self, "env", _copy_string_mapping(self.env, "env"))

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-ready representation."""
        return {
            "argv": list(self.argv),
            "cwd": self.cwd,
            "env": dict(self.env),
        }


@runtime_checkable
class WorkflowAdapter(Protocol):
    """Structural protocol implemented by workflow adapter packages."""

    metadata: WorkflowMetadata
    capabilities: WorkflowCapabilities

    def schema(self) -> WorkflowSchema:
        """Return the versioned adapter-owned authoring contract."""

    def validate(self, inputs: WorkflowInputs) -> Result[object]:
        """Validate submitted inputs and return adapter-owned validated data."""

    def preview_dag(self, inputs: WorkflowInputs) -> Result[DagPreview]:
        """Return a workflow-neutral DAG preview."""

    def plan_workspace(
        self,
        inputs: WorkflowInputs,
        workspace: str | Path,
    ) -> Result[WorkspacePlan]:
        """Plan workspace directories and files without writing them."""

    def build_command(self, plan: WorkspacePlan) -> Result[CommandSpec]:
        """Build an engine-neutral command description without executing it."""

    def extract_artifacts(
        self,
        inputs: WorkflowInputs,
        workspace: str | Path,
    ) -> Result[tuple[ExtractedArtifactCandidate, ...]]:
        """Return logical artifact candidates without persistence identities."""


@runtime_checkable
class QcSummaryExtractingAdapter(Protocol):
    """Optional adapter contract for trusted machine-readable QC summaries."""

    def qc_source_output_types(self) -> tuple[str, ...]:
        """Return exact artifact output types accepted as QC source documents."""

    def extract_qc_metrics(
        self,
        inputs: WorkflowInputs,
        sources: tuple[QcSourceDocument, ...],
    ) -> Result[tuple[ExtractedQcMetricCandidate, ...]]:
        """Map platform-vetted source bytes to neutral QC candidates."""


@runtime_checkable
class LocalRunDriver(Protocol):
    """Structural protocol for local workflow execution drivers."""

    def run(self, run_id: str, plan: "ExecutionPlan") -> "Result[RunRecord]":
        """Attempt to execute the planned run locally.

        Returns a Result wrapping the run record. Concrete drivers may return
        failures for unsupported plans, missing command specs, or runtime errors.
        """


def _normalize_required_string(value: str, name: str) -> str:
    if not isinstance(value, str):
        raise ValueError(f"{name} must be a string")
    normalized = value.strip()
    if not normalized:
        raise ValueError(f"{name} must be non-empty")
    return normalized


def _normalize_string_tuple(values: Iterable[str], name: str) -> tuple[str, ...]:
    if isinstance(values, str):
        raw_values = (values,)
    else:
        try:
            raw_values = tuple(values)
        except TypeError as exc:
            raise ValueError(f"{name} must be an iterable of strings") from exc

    normalized = []
    for value in raw_values:
        if not isinstance(value, str):
            raise ValueError(f"{name} entries must be strings")
        if not value.strip():
            raise ValueError(f"{name} entries must be non-empty")
        normalized.append(value)
    return tuple(normalized)


def _normalize_mode_tuple(values: Iterable[str], name: str) -> tuple[str, ...]:
    normalized = _normalize_string_tuple(values, name)
    if len(normalized) != len(set(normalized)):
        raise ValueError(f"{name} entries must be unique")
    for value in normalized:
        if _MODE_PATTERN.fullmatch(value) is None:
            raise ValueError(f"{name} entries must be stable lowercase tokens")
    return normalized


def _copy_mapping(mapping: Mapping[str, object], name: str) -> dict[str, object]:
    if not isinstance(mapping, Mapping):
        raise ValueError(f"{name} must be a mapping")
    return deepcopy(dict(mapping))


def _copy_json_schema(
    mapping: Mapping[str, object],
    name: str,
    *,
    dialect: str,
    root_type: str,
) -> dict[str, object]:
    if not isinstance(mapping, Mapping):
        raise ValueError(f"{name} must be a mapping")
    copied = deepcopy(dict(mapping))
    copied.setdefault("$schema", dialect)
    copied.setdefault("type", root_type)
    if copied.get("$schema") != dialect:
        raise ValueError(f"{name} must use the declared schema dialect")
    if copied.get("type") != root_type:
        raise ValueError(f"{name} must have a {root_type} root")
    _validate_json_value(copied, name)
    try:
        json.dumps(copied, allow_nan=False, ensure_ascii=False)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must contain only JSON-safe values") from exc
    return copied


def _validate_json_value(value: object, name: str) -> None:
    if value is None or isinstance(value, (str, bool, int)):
        return
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError(f"{name} must not contain NaN or Infinity")
        return
    if isinstance(value, list):
        for item in value:
            _validate_json_value(item, name)
        return
    if isinstance(value, Mapping):
        for key, item in value.items():
            if not isinstance(key, str):
                raise ValueError(f"{name} keys must be strings")
            _validate_json_value(item, name)
        return
    raise ValueError(f"{name} must contain only JSON-safe values")


def _copy_string_mapping(mapping: Mapping[str, str], name: str) -> dict[str, str]:
    if not isinstance(mapping, Mapping):
        raise ValueError(f"{name} must be a mapping")
    copied = dict(mapping)
    for key, value in copied.items():
        if not isinstance(key, str) or not isinstance(value, str):
            raise ValueError(f"{name} keys and values must be strings")
    return copied


def _copy_samples(samples: SamplePayload) -> SamplePayload:
    if samples is None:
        return samples
    if isinstance(samples, (str, Path)):
        sample_path = str(samples)
        if len(sample_path) > MAX_SAMPLE_CELL_LENGTH:
            raise ValueError("WorkflowInputs sample path exceeds the platform limit")
        if any(character in sample_path for character in ("\x00", "\t", "\n", "\r")):
            raise ValueError(
                "WorkflowInputs sample path contains a forbidden control character"
            )
        return samples
    if not isinstance(samples, list):
        raise ValueError(
            "WorkflowInputs samples must be str, Path, list[dict[str, str]], or None"
        )
    if not samples or len(samples) > MAX_SAMPLE_ROWS:
        raise ValueError(
            f"WorkflowInputs inline samples must contain 1 to {MAX_SAMPLE_ROWS} rows"
        )

    copied = []
    for row in samples:
        if not isinstance(row, dict):
            raise ValueError("WorkflowInputs sample rows must be dicts")
        copied_row = dict(row)
        if not copied_row or len(copied_row) > MAX_SAMPLE_COLUMNS:
            raise ValueError(
                "WorkflowInputs sample rows exceed the platform column limit"
            )
        for key, value in copied_row.items():
            if not isinstance(key, str) or not isinstance(value, str):
                raise ValueError(
                    "WorkflowInputs sample row keys and values must be strings"
                )
            if not key or len(key) > MAX_SAMPLE_COLUMN_NAME_LENGTH:
                raise ValueError(
                    "WorkflowInputs sample column name exceeds the platform limit"
                )
            if any(character in key for character in ("\x00", "\t", "\n", "\r")):
                raise ValueError(
                    "WorkflowInputs sample column name contains a forbidden control character"
                )
            if len(value) > MAX_SAMPLE_CELL_LENGTH:
                raise ValueError(
                    "WorkflowInputs sample cell exceeds the platform limit"
                )
            if any(character in value for character in ("\x00", "\t", "\n", "\r")):
                raise ValueError(
                    "WorkflowInputs sample cell contains a forbidden control character"
                )
        copied.append(copied_row)
    return copied


def _normalize_instance_tuple(
    values: Iterable[object],
    expected_type: type,
    name: str,
) -> tuple[Any, ...]:
    try:
        value_tuple = tuple(values)
    except TypeError as exc:
        raise ValueError(f"{name} must be an iterable") from exc
    for value in value_tuple:
        if not isinstance(value, expected_type):
            raise ValueError(
                f"{name} entries must be {expected_type.__name__} instances"
            )
    return value_tuple


def _normalize_files(
    files: Iterable[tuple[str, bytes]],
) -> tuple[tuple[str, bytes], ...]:
    try:
        file_tuple = tuple(files)
    except TypeError as exc:
        raise ValueError("WorkspacePlan files must be an iterable") from exc

    normalized = []
    for entry in file_tuple:
        if not isinstance(entry, tuple) or len(entry) != 2:
            raise ValueError("WorkspacePlan files entries must be (path, bytes)")
        path, contents = entry
        path = _normalize_required_string(path, "WorkspacePlan file path")
        if not isinstance(contents, bytes):
            raise ValueError("WorkspacePlan file contents must be bytes")
        normalized.append((path, contents))
    return tuple(normalized)
