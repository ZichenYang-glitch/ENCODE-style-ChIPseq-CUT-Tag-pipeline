"""ENCODE-style workflow adapter wrapper."""

from __future__ import annotations

import csv
import io
import re
import tempfile
from collections.abc import Mapping
from copy import deepcopy
from pathlib import Path
from typing import Any

import yaml

from encode_pipeline import __version__
from encode_pipeline.artifacts import (
    Artifact,
    artifacts_by_manifest_output_type,
    load_catalog,
)
from encode_pipeline.platform.adapters import (
    ARTIFACT_EXTRACT_CAPABILITY,
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
    WorkflowInputs,
    WorkflowMetadata,
    WorkflowSchema,
    WorkspacePlan,
)
from encode_pipeline.platform.results import Issue, Result


_WORKFLOW_ID = "encode-style-chipseq-cuttag-atac-mnase"

_WORKSPACE_SAMPLE_COLUMNS = (
    "sample",
    "fastq_1",
    "fastq_2",
    "layout",
    "assay",
    "target",
    "peak_mode",
    "genome",
    "bowtie2_index",
    "control_bam",
    "role",
    "control_sample",
    "experiment",
    "condition",
    "replicate",
    "biological_replicate",
    "technical_replicate",
)

_SAMPLE_DICT_TO_TSV = {
    "id": "sample",
    "fq1": "fastq_1",
    "fq2": "fastq_2",
    "layout": "layout",
    "assay": "assay",
    "target": "target",
    "peak_mode": "peak_mode",
    "genome": "genome",
    "bt2_idx": "bowtie2_index",
    "control_bam": "control_bam",
    "role": "role",
    "control_sample": "control_sample",
    "experiment": "experiment",
    "condition": "condition",
    "replicate": "replicate",
    "biological_replicate": "biological_replicate",
    "technical_replicate": "technical_replicate",
}


def _render_config_yaml(config: dict[str, Any]) -> bytes:
    """Serialize a validated config dict to canonical UTF-8 YAML bytes.

    Validation builds some nested mappings from sets of known keys. Sorting
    every mapping keeps preflight workspace bytes stable across API and worker
    processes with different Python hash seeds.
    """
    try:
        return yaml.safe_dump(
            config,
            default_flow_style=False,
            sort_keys=True,
            allow_unicode=True,
        ).encode("utf-8")
    except yaml.YAMLError as exc:
        raise ValueError("Config could not be serialized to YAML.") from exc


def _prepare_workspace_config(validated_config: dict[str, Any]) -> dict[str, Any]:
    """Return a workspace-local copy of the config with samples/outdir overrides."""
    workspace_config = dict(validated_config)
    workspace_config["samples"] = "config/samples.tsv"
    workspace_config["outdir"] = "results"
    return workspace_config


def _render_samples_tsv(sample_rows: list[dict[str, Any]]) -> bytes:
    """Serialize validated sample rows to canonical TSV bytes."""
    output = io.StringIO()
    writer = csv.DictWriter(
        output,
        fieldnames=_WORKSPACE_SAMPLE_COLUMNS,
        extrasaction="ignore",
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    for row in sample_rows:
        mapped: dict[str, str] = {}
        for src_key, tsv_col in _SAMPLE_DICT_TO_TSV.items():
            value = row.get(src_key)
            mapped[tsv_col] = "" if value is None else str(value)
        for value in mapped.values():
            if "\t" in value or "\n" in value or "\r" in value:
                raise ValueError(
                    "Sample value contains forbidden tab or newline character."
                )
        writer.writerow(mapped)
    return output.getvalue().encode("utf-8")


def _enforce_external_input_policy(
    validated_config: dict[str, Any],
    sample_rows: list[dict[str, Any]],
) -> Issue | None:
    """Fail fast if any external reference is a relative path."""
    sample_fields = ("fq1", "fq2", "bt2_idx", "control_bam")
    for index, row in enumerate(sample_rows):
        for field in sample_fields:
            value = row.get(field)
            if isinstance(value, str) and value and not Path(value).is_absolute():
                return Issue(
                    code="ENCODE_WORKSPACE_RELATIVE_EXTERNAL_PATH",
                    message="External input path must be absolute.",
                    severity="error",
                    path=f"samples[{index}].{_SAMPLE_DICT_TO_TSV[field]}",
                    source="adapter",
                    context={"field": _SAMPLE_DICT_TO_TSV[field]},
                )

    genome_resources = validated_config.get("genome_resources")
    if isinstance(genome_resources, Mapping):
        for genome in sorted(genome_resources.keys()):
            resource = genome_resources[genome]
            if not isinstance(resource, Mapping):
                continue
            for field in ("chrom_sizes", "blacklist", "gtf", "reference_fasta"):
                value = resource.get(field)
                if isinstance(value, str) and value and not Path(value).is_absolute():
                    return Issue(
                        code="ENCODE_WORKSPACE_RELATIVE_EXTERNAL_PATH",
                        message="External input path must be absolute.",
                        severity="error",
                        path=f"config.genome_resources.{genome}.{field}",
                        source="adapter",
                        context={"field": field},
                    )

    return None


_GENERIC_ADAPTER_MESSAGES = {
    "ENCODE_ADAPTER_UNSUPPORTED": "The requested feature is not supported by the current adapter.",
    "ENCODE_CONFIG_INVALID": "Workflow configuration is invalid.",
    "ENCODE_SAMPLES_INVALID": "Sample sheet is invalid.",
    "ENCODE_RESOURCES_INVALID": "Genome resources are invalid.",
    "ENCODE_OPTIONS_INVALID": "Adapter options are invalid.",
}

_DEFAULT_ADAPTER_VALIDATION_MESSAGE = "Workflow input validation failed."


def _sanitize_issue(issue: Issue) -> Issue:
    """Return a path-free, strictly bounded copy of a validate() issue."""
    safe_message = _GENERIC_ADAPTER_MESSAGES.get(
        issue.code, _DEFAULT_ADAPTER_VALIDATION_MESSAGE
    )
    return Issue(
        code=issue.code,
        message=safe_message,
        severity=issue.severity,
        path=issue.path,
        source=issue.source,
        technical_message=safe_message,
        hint=None,
        context={},
    )


class EncodeStyleWorkflowAdapter:
    """Adapter wrapper for the current ENCODE-style Snakemake workflow."""

    metadata = WorkflowMetadata(
        workflow_id=_WORKFLOW_ID,
        name="ENCODE-style ChIP-seq/CUT&Tag/ATAC/MNase",
        version=__version__ or "0.0.0",
        description=(
            "Adapter wrapper for the current ENCODE-style ChIP-seq, "
            "CUT&Tag, ATAC-seq, and MNase-seq workflow."
        ),
        engines=("snakemake",),
        tags=("chipseq", "cuttag", "atac", "mnase", "encode-style"),
    )
    capabilities = WorkflowCapabilities(
        supports=(
            VALIDATION_CAPABILITY,
            WORKSPACE_PLAN_CAPABILITY,
            INPUT_AUTHORING_CAPABILITY,
            ARTIFACT_EXTRACT_CAPABILITY,
            QC_SUMMARY_EXTRACT_CAPABILITY,
        )
    )

    def __init__(
        self,
        *,
        catalog: tuple[Artifact, ...] | list[Artifact] | None = None,
        catalog_path: str | None = None,
    ) -> None:
        if catalog is not None and catalog_path is not None:
            raise ValueError("catalog and catalog_path are mutually exclusive")
        loaded = load_catalog(catalog_path) if catalog is None else list(catalog)
        self._artifact_catalog = tuple(loaded)
        self._artifacts_by_manifest_type = artifacts_by_manifest_output_type(loaded)
        self._artifacts_by_id = {artifact.id: artifact for artifact in loaded}
        if len(self._artifacts_by_id) != len(loaded):
            raise ValueError("artifact catalog contains duplicate ids")

    def schema(self) -> WorkflowSchema:
        """Return the versioned renderable authoring contract."""
        from encode_pipeline.adapters.encode_authoring import (
            build_encode_authoring_schema,
        )

        return build_encode_authoring_schema()

    def validate(self, inputs: WorkflowInputs) -> Result[object]:
        """Validate inputs using current package behavior.

        The success value shape is adapter-private. Platform code must not
        depend on ``{"config": ..., "samples": ...}``.
        """
        option_issue = _validate_options(inputs.options)
        if option_issue is not None:
            return Result.failure([option_issue])
        strict_inputs = bool(inputs.options.get("strict_inputs", False))
        config = deepcopy(dict(inputs.config))
        if isinstance(inputs.samples, list):
            result = _validate_inline_inputs(
                config,
                inputs.samples,
                strict_inputs=strict_inputs,
            )
        else:
            if isinstance(inputs.samples, (str, Path)):
                config["samples"] = str(inputs.samples)
            result = _validate_scientific_inputs(
                config,
                strict_inputs=strict_inputs,
                sanitize_errors=False,
            )
        if result.is_failure:
            return result

        validated_config = result.value["config"]
        samples = result.value["samples"]
        policy_issue = _enforce_external_input_policy(validated_config, samples)
        if policy_issue is not None:
            return Result.failure([policy_issue])
        return result

    def preview_dag(self, inputs: WorkflowInputs) -> Result[DagPreview]:
        """Return unsupported until DAG preview is wired through the adapter."""
        return _unsupported_method("preview_dag")

    def plan_workspace(
        self,
        inputs: WorkflowInputs,
        workspace: str | Path,
    ) -> Result[WorkspacePlan]:
        """Plan workspace directories and files for the ENCODE workflow."""
        validated = self.validate(inputs)
        if validated.is_failure:
            return Result.failure(
                [
                    issue
                    if issue.code == "ENCODE_WORKSPACE_RELATIVE_EXTERNAL_PATH"
                    else _sanitize_issue(issue)
                    for issue in validated.issues
                ]
            )

        validated_value = validated.value
        validated_config = validated_value["config"]
        sample_rows = validated_value["samples"]

        workspace_config = _prepare_workspace_config(validated_config)

        try:
            config_yaml = _render_config_yaml(workspace_config)
        except ValueError:
            return Result.failure(
                [
                    Issue(
                        code="ENCODE_WORKSPACE_RENDER_FAILED",
                        message="Config could not be rendered to YAML.",
                        severity="error",
                        path="config",
                        source="adapter",
                    )
                ]
            )

        try:
            samples_tsv = _render_samples_tsv(sample_rows)
        except ValueError:
            return Result.failure(
                [
                    Issue(
                        code="ENCODE_WORKSPACE_RENDER_FAILED",
                        message="Samples could not be rendered to TSV.",
                        severity="error",
                        path="samples",
                        source="adapter",
                    )
                ]
            )

        workspace_plan = WorkspacePlan(
            directories=("logs", "results"),
            files=(
                ("config/config.yaml", config_yaml),
                ("config/samples.tsv", samples_tsv),
            ),
        )

        return Result.success(
            workspace_plan,
            issues=[
                Issue(
                    code="ENCODE_WORKSPACE_PLANNING_COMPLETE",
                    message="Workspace plan created.",
                    severity="info",
                    path="workspace_plan",
                    source="adapter",
                    context={"file_count": 2, "directory_count": 2},
                )
            ],
        )

    def build_command(self, plan: WorkspacePlan) -> Result[CommandSpec]:
        """Return unsupported until command construction is designed."""
        return _unsupported_method("build_command")

    def extract_artifacts(
        self,
        inputs: WorkflowInputs,
        workspace: str | Path,
    ) -> Result[tuple[ExtractedArtifactCandidate, ...]]:
        """Map the existing manifest vocabulary to regular-file candidates."""
        if not isinstance(inputs, WorkflowInputs):
            return _artifact_extraction_failure("ENCODE_ARTIFACT_INPUTS_INVALID")
        try:
            from encode_pipeline.manifest.make import build_manifest_rows

            workspace_path = Path(workspace).absolute()
            config_value = yaml.safe_load(
                (workspace_path / "config" / "config.yaml").read_text(encoding="utf-8")
            )
            if not isinstance(config_value, dict):
                raise ValueError("workspace config is not a mapping")
            manifest_config = deepcopy(config_value)
            manifest_config["outdir"] = str(workspace_path / "results")
            samples_path = workspace_path / "config" / "samples.tsv"
            manifest_config["samples"] = str(samples_path)
            rows, _missing, _not_applicable = build_manifest_rows(
                manifest_config,
                samples_path=str(samples_path),
            )
        except (OSError, ValueError, yaml.YAMLError):
            return _artifact_extraction_failure("ENCODE_ARTIFACT_MANIFEST_FAILED")

        candidates: list[ExtractedArtifactCandidate] = []
        try:
            for row in rows:
                output_type = row.get("output_type")
                if not isinstance(output_type, str) or not output_type:
                    raise ValueError("manifest output type is invalid")
                artifact = self._catalog_artifact_for_output_type(output_type)
                status = row.get("status")
                if status not in {"present", "missing", "not_applicable"}:
                    raise ValueError("manifest status is invalid")
                if status != "present" or artifact.path_template.endswith("/"):
                    continue
                relative_path = (
                    Path(str(row.get("path", "")))
                    .relative_to(workspace_path)
                    .as_posix()
                )
                candidates.append(
                    ExtractedArtifactCandidate(
                        output_type=output_type,
                        relative_path=relative_path,
                        mime_type=_artifact_mime_type(relative_path),
                        metadata=_artifact_candidate_metadata(artifact, row),
                    )
                )

            manifest_path = workspace_path / "results/multiqc/result_manifest.tsv"
            if manifest_path.exists():
                artifact = self._artifacts_by_id["result_manifest"]
                candidates.append(
                    ExtractedArtifactCandidate(
                        output_type="result_manifest",
                        relative_path="results/multiqc/result_manifest.tsv",
                        mime_type="text/tab-separated-values",
                        metadata=_artifact_candidate_metadata(artifact, {}),
                    )
                )
        except (KeyError, TypeError, ValueError):
            return _artifact_extraction_failure("ENCODE_ARTIFACT_CATALOG_MISMATCH")

        return Result.success(
            tuple(
                sorted(
                    candidates, key=lambda item: (item.output_type, item.relative_path)
                )
            )
        )

    def _catalog_artifact_for_output_type(self, output_type: str) -> Artifact:
        exact = self._artifacts_by_manifest_type.get(output_type)
        if exact is not None:
            return exact
        matches = []
        for template, artifact in self._artifacts_by_manifest_type.items():
            if "<N>" not in template:
                continue
            pattern = re.escape(template).replace(re.escape("<N>"), r"[1-9]\d*")
            if re.fullmatch(pattern, output_type):
                matches.append(artifact)
        if len(matches) != 1:
            raise ValueError("manifest output type has no unique catalog mapping")
        return matches[0]

    def qc_source_output_types(self) -> tuple[str, ...]:
        """Return exact persisted artifact types accepted by the QC parser."""
        from encode_pipeline.adapters.encode_qc import QC_SOURCE_OUTPUT_TYPES

        return QC_SOURCE_OUTPUT_TYPES

    def extract_qc_metrics(
        self,
        inputs: WorkflowInputs,
        sources: tuple[QcSourceDocument, ...],
    ) -> Result[tuple[ExtractedQcMetricCandidate, ...]]:
        """Map platform-vetted summary bytes to neutral numeric metrics."""
        if not isinstance(inputs, WorkflowInputs):
            return _qc_summary_failure()
        try:
            from encode_pipeline.adapters.encode_qc import parse_encode_qc_sources

            return Result.success(parse_encode_qc_sources(sources))
        except (TypeError, UnicodeError, ValueError):
            return _qc_summary_failure()


def _write_inline_samples_tsv(
    path: Path,
    sample_rows: list[dict[str, str]],
) -> None:
    """Write raw external-column rows for the canonical scientific loader."""
    if not sample_rows:
        raise ValueError("inline sample rows must be non-empty")
    allowed_columns = frozenset(_WORKSPACE_SAMPLE_COLUMNS)
    for row in sample_rows:
        if not set(row).issubset(allowed_columns):
            raise ValueError("inline sample row contains an unknown column")
        for value in row.values():
            if any(character in value for character in ("\x00", "\t", "\n", "\r")):
                raise ValueError("inline sample cell contains a control character")

    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=_WORKSPACE_SAMPLE_COLUMNS,
            extrasaction="raise",
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in sample_rows:
            writer.writerow(
                {column: row.get(column, "") for column in writer.fieldnames}
            )


def _validate_inline_inputs(
    config: dict[str, Any],
    sample_rows: list[dict[str, str]],
    *,
    strict_inputs: bool,
) -> Result[dict[str, Any]]:
    """Validate inline rows through one private temporary TSV bridge."""
    try:
        with tempfile.TemporaryDirectory(
            prefix="encode-platform-inline-samples-"
        ) as temporary_root:
            samples_path = Path(temporary_root) / "samples.tsv"
            _write_inline_samples_tsv(samples_path, sample_rows)
            config["samples"] = str(samples_path)
            result = _validate_scientific_inputs(
                config,
                strict_inputs=strict_inputs,
                sanitize_errors=True,
            )
    except (OSError, ValueError, csv.Error):
        return _inline_samples_failure()

    if result.is_failure:
        return result
    validated_config = dict(result.value["config"])
    validated_config.pop("samples", None)
    return Result.success(
        {
            "config": validated_config,
            "samples": result.value["samples"],
        }
    )


def _validate_scientific_inputs(
    config: dict[str, Any],
    *,
    strict_inputs: bool,
    sanitize_errors: bool,
) -> Result[dict[str, Any]]:
    """Run the existing config, sample, and resource validators unchanged."""
    from encode_pipeline.config.validate import (
        ValidationError,
        validate_config,
        validate_picard_reference_resources,
        validate_tss_annotation_resources,
    )
    from encode_pipeline.samples.load import load_and_validate_samples

    try:
        validated_config = validate_config(config)
    except ValidationError as exc:
        return _scientific_validation_failure(
            exc,
            code="ENCODE_CONFIG_INVALID",
            source="config",
            path="config",
            sanitize=sanitize_errors,
        )

    flags = _sample_validation_flags(validated_config)
    try:
        samples = load_and_validate_samples(
            validated_config["samples"],
            use_control=validated_config["use_control"],
            stage5_enabled=validated_config.get("stage5", False),
            strict_inputs=strict_inputs,
            reproducibility_idr_atac_narrow=flags["atac_narrow"],
            reproducibility_idr_cuttag_narrow=flags["cuttag_narrow"],
            reproducibility_idr_chipseq_broad=flags["chipseq_broad"],
            reproducibility_idr_cuttag_broad=flags["cuttag_broad"],
        )
    except ValidationError as exc:
        return _scientific_validation_failure(
            exc,
            code="ENCODE_SAMPLES_INVALID",
            source="samples",
            path="samples",
            sanitize=sanitize_errors,
        )

    try:
        validate_picard_reference_resources(validated_config, samples)
        validate_tss_annotation_resources(validated_config, samples)
    except ValidationError as exc:
        return _scientific_validation_failure(
            exc,
            code="ENCODE_RESOURCES_INVALID",
            source="config",
            path="config.genome_resources",
            sanitize=sanitize_errors,
        )

    return Result.success({"config": validated_config, "samples": samples})


def _scientific_validation_failure(
    exc: Exception,
    *,
    code: str,
    source: str,
    path: str,
    sanitize: bool,
) -> Result[dict[str, Any]]:
    issue = _issue_from_exception(exc, code=code, source=source, path=path)
    return Result.failure([_sanitize_issue(issue) if sanitize else issue])


def _inline_samples_failure() -> Result[dict[str, Any]]:
    return Result.failure(
        [
            _sanitize_issue(
                Issue(
                    code="ENCODE_SAMPLES_INVALID",
                    message="Sample sheet is invalid.",
                    source="samples",
                    path="samples",
                )
            )
        ]
    )


def _validate_options(options: dict[str, object]) -> Issue | None:
    unsupported = sorted(str(key) for key in options if key != "strict_inputs")
    if unsupported:
        return Issue(
            code="ENCODE_OPTIONS_INVALID",
            message="Unsupported adapter options were provided.",
            source="adapter",
            path="options",
            context={"unsupported_options": unsupported},
        )

    if "strict_inputs" in options and not isinstance(options["strict_inputs"], bool):
        return Issue(
            code="ENCODE_OPTIONS_INVALID",
            message="strict_inputs must be a boolean",
            source="adapter",
            path="options.strict_inputs",
        )
    return None


def _artifact_candidate_metadata(
    artifact: Artifact,
    row: Mapping[str, object],
) -> dict[str, object]:
    metadata: dict[str, object] = {
        "catalog_id": artifact.id,
        "scope": artifact.scope,
        "level": artifact.level,
    }
    for key in (
        "sample_id",
        "experiment_id",
        "assay",
        "target",
        "genome",
        "method",
        "qc_flag",
    ):
        value = row.get(key)
        if isinstance(value, str) and value:
            metadata[key] = value
    return metadata


def _artifact_mime_type(relative_path: str) -> str:
    lowered = relative_path.lower()
    if lowered.endswith(".tsv"):
        return "text/tab-separated-values"
    if lowered.endswith(".html"):
        return "text/html"
    if lowered.endswith((".bed", ".bedgraph", ".narrowpeak", ".broadpeak")):
        return "text/plain"
    return "application/octet-stream"


def _artifact_extraction_failure(code: str) -> Result[Any]:
    return Result.failure(
        [
            Issue(
                code=code,
                message="Workflow artifacts could not be discovered.",
                severity="error",
                path="artifacts",
                source="adapter",
            )
        ]
    )


def _qc_summary_failure() -> Result[Any]:
    return Result.failure(
        [
            Issue(
                code="ENCODE_QC_SUMMARY_INVALID",
                message="Workflow QC summaries could not be interpreted.",
                severity="error",
                path="qc_metrics",
                source="adapter",
            )
        ]
    )


def _sample_validation_flags(validated_config: dict[str, Any]) -> dict[str, bool]:
    repro = validated_config.get("reproducibility", {})
    idr = repro.get("idr", {})
    enabled = repro.get("enabled", False)
    return {
        "atac_narrow": bool(enabled and idr.get("atac_narrow", False)),
        "cuttag_narrow": bool(enabled and idr.get("cuttag_narrow", False)),
        "chipseq_broad": bool(enabled and idr.get("chipseq_broad_experimental", False)),
        "cuttag_broad": bool(enabled and idr.get("cuttag_broad_experimental", False)),
    }


def _issue_from_exception(
    exc: Exception,
    *,
    code: str,
    source: str,
    path: str,
) -> Issue:
    return Issue.from_exception(
        exc,
        code=code,
        source=source,
        path=path,
    )


def _unsupported_method(method: str) -> Result[Any]:
    return Result.failure(
        [
            Issue(
                code="ENCODE_ADAPTER_UNSUPPORTED",
                message=f"{method} is not supported by the current adapter.",
                source="adapter",
                path=method,
                context={"method": method},
            )
        ]
    )
