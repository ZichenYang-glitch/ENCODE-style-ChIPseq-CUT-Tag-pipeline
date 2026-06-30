"""ENCODE-style workflow adapter wrapper."""

from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from typing import Any

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
    capabilities = WorkflowCapabilities(supports=("validation",))

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
            return Result.failure([
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
            ])
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
            return Result.failure([
                _issue_from_exception(
                    exc,
                    code="ENCODE_CONFIG_INVALID",
                    source="config",
                    path="config",
                )
            ])

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
            return Result.failure([
                _issue_from_exception(
                    exc,
                    code="ENCODE_SAMPLES_INVALID",
                    source="samples",
                    path="samples",
                )
            ])

        try:
            validate_picard_reference_resources(validated_config, samples)
            validate_tss_annotation_resources(validated_config, samples)
        except ValidationError as exc:
            return Result.failure([
                _issue_from_exception(
                    exc,
                    code="ENCODE_RESOURCES_INVALID",
                    source="config",
                    path="config.genome_resources",
                )
            ])

        return Result.success({"config": validated_config, "samples": samples})

    def preview_dag(self, inputs: WorkflowInputs) -> Result[DagPreview]:
        """Return unsupported until DAG preview is wired through the adapter."""
        return _unsupported_method("preview_dag")

    def plan_workspace(
        self,
        inputs: WorkflowInputs,
        workspace: str | Path,
    ) -> Result[WorkspacePlan]:
        """Return unsupported until run workspace planning is designed."""
        return _unsupported_method("plan_workspace")

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
        "chipseq_broad": bool(
            enabled and idr.get("chipseq_broad_experimental", False)
        ),
        "cuttag_broad": bool(
            enabled and idr.get("cuttag_broad_experimental", False)
        ),
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
    return Result.failure([
        Issue(
            code="ENCODE_ADAPTER_UNSUPPORTED",
            message=f"{method} is not supported by the current adapter.",
            source="adapter",
            path=method,
            context={"method": method},
        )
    ])
