"""Workflow-platform adapter contract primitives."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from copy import deepcopy
from dataclasses import dataclass, field
from decimal import Decimal
from pathlib import Path
from typing import TYPE_CHECKING, Any, Protocol, runtime_checkable

from encode_pipeline.platform.results import Result

if TYPE_CHECKING:
    from encode_pipeline.platform.planning import ExecutionPlan
    from encode_pipeline.platform.runs import RunRecord


SamplePayload = str | Path | list[dict[str, str]] | None


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
    """String capability labels supported by a workflow adapter."""

    supports: tuple[str, ...]

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "supports",
            _normalize_string_tuple(self.supports, "supports"),
        )

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-ready representation."""
        return {"supports": list(self.supports)}


@dataclass(frozen=True)
class WorkflowSchema:
    """Adapter-owned schemas for config, samples, and options."""

    config_schema: Mapping[str, object] = field(default_factory=dict)
    sample_schema: Mapping[str, object] = field(default_factory=dict)
    option_schema: Mapping[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "config_schema",
            _copy_mapping(self.config_schema, "config_schema"),
        )
        object.__setattr__(
            self,
            "sample_schema",
            _copy_mapping(self.sample_schema, "sample_schema"),
        )
        object.__setattr__(
            self,
            "option_schema",
            _copy_mapping(self.option_schema, "option_schema"),
        )

    def to_dict(self) -> dict[str, Any]:
        """Return fresh schema dict copies."""
        return {
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
        """Return adapter-owned config, sample, and option schemas."""

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


def _copy_mapping(mapping: Mapping[str, object], name: str) -> dict[str, object]:
    if not isinstance(mapping, Mapping):
        raise ValueError(f"{name} must be a mapping")
    return deepcopy(dict(mapping))


def _copy_string_mapping(mapping: Mapping[str, str], name: str) -> dict[str, str]:
    if not isinstance(mapping, Mapping):
        raise ValueError(f"{name} must be a mapping")
    copied = dict(mapping)
    for key, value in copied.items():
        if not isinstance(key, str) or not isinstance(value, str):
            raise ValueError(f"{name} keys and values must be strings")
    return copied


def _copy_samples(samples: SamplePayload) -> SamplePayload:
    if samples is None or isinstance(samples, (str, Path)):
        return samples
    if not isinstance(samples, list):
        raise ValueError(
            "WorkflowInputs samples must be str, Path, list[dict[str, str]], or None"
        )

    copied = []
    for row in samples:
        if not isinstance(row, dict):
            raise ValueError("WorkflowInputs sample rows must be dicts")
        copied_row = dict(row)
        for key, value in copied_row.items():
            if not isinstance(key, str) or not isinstance(value, str):
                raise ValueError(
                    "WorkflowInputs sample row keys and values must be strings"
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
