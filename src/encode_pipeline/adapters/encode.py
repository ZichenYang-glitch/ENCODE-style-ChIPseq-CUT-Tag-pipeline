"""ENCODE-style workflow adapter wrapper."""

from __future__ import annotations

import csv
import io
from collections.abc import Mapping
from copy import deepcopy
from pathlib import Path
from typing import Any

import yaml

from encode_pipeline import __version__
from encode_pipeline.platform.adapters import (
    CommandSpec,
    DagPreview,
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
    capabilities = WorkflowCapabilities(supports=("validation", "workspace_plan"))

    def schema(self) -> WorkflowSchema:
        """Return minimal UI/agent discovery hints, not complete JSON Schema."""
        from encode_pipeline.config import defaults

        return WorkflowSchema(
            config_schema={
                "schema_kind": "hints",
                "description": (
                    "Minimal configuration hints for the ENCODE-style "
                    "workflow adapter. This is not complete JSON Schema."
                ),
                "required": ["samples"],
                "properties": {
                    "samples": {
                        "type": "string",
                        "description": "Path to the sample sheet TSV.",
                    },
                    "use_control": {
                        "type": "boolean",
                        "default": False,
                    },
                    "outdir": {
                        "type": "string",
                        "default": "results",
                    },
                },
            },
            sample_schema={
                "schema_kind": "hints",
                "description": (
                    "Minimal sample-sheet hints for the ENCODE-style "
                    "workflow adapter. This is not complete JSON Schema."
                ),
                "format": "tsv",
                "required_columns": list(defaults.SAMPLE_REQUIRED_COLUMNS),
                "allowed_values": {
                    "assay": list(defaults.ASSAYS),
                    "layout": list(defaults.LAYOUTS),
                    "peak_mode": list(defaults.PEAK_MODES),
                    "role": list(defaults.ROLES),
                },
            },
            option_schema={
                "schema_kind": "hints",
                "description": (
                    "Adapter options supported by validate(). This is not "
                    "complete JSON Schema."
                ),
                "properties": {
                    "strict_inputs": {
                        "type": "boolean",
                        "default": False,
                        "description": (
                            "Validate FASTQ and Bowtie2 index file existence."
                        ),
                    },
                },
            },
        )

    def validate(self, inputs: WorkflowInputs) -> Result[object]:
        """Validate inputs using current package behavior.

        The success value shape is adapter-private. Platform code must not
        depend on ``{"config": ..., "samples": ...}``.
        """
        config = deepcopy(dict(inputs.config))

        option_issue = _validate_options(inputs.options)
        if option_issue is not None:
            return Result.failure([option_issue])
        strict_inputs = bool(inputs.options.get("strict_inputs", False))

        if isinstance(inputs.samples, list):
            return Result.failure(
                [
                    Issue(
                        code="ENCODE_ADAPTER_UNSUPPORTED",
                        message=(
                            "The current ENCODE-style adapter supports TSV/CSV "
                            "sample paths only; inline sample rows are not "
                            "supported yet."
                        ),
                        source="adapter",
                        path="samples",
                        context={"feature": "inline_samples"},
                    )
                ]
            )
        if isinstance(inputs.samples, (str, Path)):
            config["samples"] = str(inputs.samples)

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
            return Result.failure(
                [
                    _issue_from_exception(
                        exc,
                        code="ENCODE_CONFIG_INVALID",
                        source="config",
                        path="config",
                    )
                ]
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
            return Result.failure(
                [
                    _issue_from_exception(
                        exc,
                        code="ENCODE_SAMPLES_INVALID",
                        source="samples",
                        path="samples",
                    )
                ]
            )

        try:
            validate_picard_reference_resources(validated_config, samples)
            validate_tss_annotation_resources(validated_config, samples)
        except ValidationError as exc:
            return Result.failure(
                [
                    _issue_from_exception(
                        exc,
                        code="ENCODE_RESOURCES_INVALID",
                        source="config",
                        path="config.genome_resources",
                    )
                ]
            )

        return Result.success({"config": validated_config, "samples": samples})

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
                [_sanitize_issue(issue) for issue in validated.issues]
            )

        validated_value = validated.value
        validated_config = validated_value["config"]
        sample_rows = validated_value["samples"]

        policy_result = _enforce_external_input_policy(validated_config, sample_rows)
        if policy_result is not None:
            return Result.failure([policy_result])

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
