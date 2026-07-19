"""Strict bulk RNA-seq input validation without submitted-path I/O."""

from __future__ import annotations

from collections.abc import Mapping
from copy import deepcopy
import math
from pathlib import PurePosixPath
import re
from typing import Any

import jsonschema

from encode_pipeline.adapters.bulk_rnaseq.authoring import (
    ANALYSIS_DEFAULTS,
    MULTIQC_SAMPLE_CLEAN_TOKENS,
    OUTPUT_DEFAULTS,
    QC_DEFAULTS,
    RIBOSOMAL_RNA_REMOVAL_DEFAULTS,
    SCHEMA_VERSION,
    TRIMMING_DEFAULTS,
    build_bulk_rnaseq_authoring_schema,
)
from encode_pipeline.adapters.bulk_rnaseq.upstream import (
    ADVANCED_NATIVE_PARAMETERS,
    DANGEROUS_RUNTIME_PARAMETERS,
    KNOWN_UPSTREAM_PARAMETERS,
    NFCORE_RNASEQ_COMMIT,
    NFCORE_RNASEQ_RELEASE,
    PARAMETER_IMPACT_BY_NATIVE_NAME,
    PLATFORM_OWNED_NATIVE_PARAMETERS,
    RAW_ARGUMENT_NATIVE_PARAMETERS,
    STANDARD_NATIVE_PARAMETERS,
    UNSUPPORTED_NATIVE_PARAMETERS,
    UPSTREAM_PARAMETER_SCHEMA_SHA256,
    UPSTREAM_SAMPLESHEET_SCHEMA_SHA256,
    upstream_parameter_properties,
)
from encode_pipeline.adapters.bulk_rnaseq.results_contract import (
    effective_downstream_layout,
    effective_rseqc_modules,
)
from encode_pipeline.platform.adapters import WorkflowInputs
from encode_pipeline.platform.results import Issue, Result


_WORKFLOW_ID = "bulk-rnaseq"
_RSEQC_MODULES = frozenset(
    {
        "bam_stat",
        "inner_distance",
        "infer_experiment",
        "junction_annotation",
        "junction_saturation",
        "read_distribution",
        "read_duplication",
        "tin",
    }
)
_SAFE_UMI_REGEX = re.compile(r"^[A-Za-z0-9_?P<>(){}\[\].,+*^:\\-]+$")
_MULTIQC_EXTRA_TRUNCATE_TOKENS = (
    ".sortmerna",
    ".umi_dedup",
    "_val",
    ".markdup",
    "_primary",
    ".namesorted",
    ".umi_extract",
)
_MULTIQC_EXTRA_REMOVE_TOKENS = ("_raw", "_trimmed", "_filtered")
_MULTIQC_EXTRA_TOKENS = frozenset(
    (*_MULTIQC_EXTRA_TRUNCATE_TOKENS, *_MULTIQC_EXTRA_REMOVE_TOKENS)
)
_MULTIQC_FN_CLEAN_TRIM = (
    ".",
    ":",
    "_",
    "-",
    ".r",
    "_val",
    ".idxstats",
    "_trimmed",
    ".trimmed",
    ".csv",
    ".yaml",
    ".yml",
    ".json",
    "_mqc",
    "short_summary_",
    "_summary",
    ".summary",
    ".align",
    ".h5",
    "_matrix",
    ".stats",
    ".hist",
    ".phased",
    ".tar",
    "runs_",
    ".qc",
)


def validate_bulk_rnaseq_inputs(inputs: WorkflowInputs) -> Result[object]:
    """Validate and deterministically normalize one contract-only submission."""
    if not isinstance(inputs, WorkflowInputs):
        return _failure("BULK_RNASEQ_INPUTS_INVALID", "inputs")

    config = inputs.config
    if not _is_json_value(config):
        return _failure("BULK_RNASEQ_CONFIG_INVALID", "config")
    if not _is_json_value(inputs.options):
        return _failure("BULK_RNASEQ_OPTIONS_INVALID", "options")
    unknown_config = set(config).difference({"standard", "advanced"})
    if unknown_config:
        return _failure("BULK_RNASEQ_CONFIG_UNKNOWN", "config")
    standard = config.get("standard")
    if not isinstance(standard, Mapping):
        return _failure("BULK_RNASEQ_STANDARD_INVALID", "config.standard")
    advanced = config.get("advanced", {})
    if not isinstance(advanced, Mapping):
        return _failure("BULK_RNASEQ_ADVANCED_INVALID", "config.advanced")

    schema = build_bulk_rnaseq_authoring_schema().to_dict()
    standard_schema = schema["config_schema"]["properties"]["standard"]
    standard_error = _first_schema_error(standard_schema, standard)
    if standard_error is not None:
        return _failure(
            "BULK_RNASEQ_STANDARD_INVALID",
            _schema_path("config.standard", standard_error),
        )

    if inputs.options:
        return _failure("BULK_RNASEQ_OPTIONS_INVALID", "options")
    if not isinstance(inputs.samples, list):
        return _failure("BULK_RNASEQ_SAMPLES_INVALID", "samples")
    sample_error = _first_schema_error(schema["sample_schema"], inputs.samples)
    if sample_error is not None:
        return _failure(
            "BULK_RNASEQ_SAMPLES_INVALID",
            _schema_path("samples", sample_error),
        )

    semantic_issue = _validate_standard_semantics(standard, inputs.samples)
    if semantic_issue is not None:
        return Result.failure([semantic_issue])
    sample_issue = _validate_sample_semantics(inputs.samples)
    if sample_issue is not None:
        return Result.failure([sample_issue])

    advanced_result = _validate_advanced(advanced, standard, inputs.samples)
    if advanced_result.is_failure:
        return Result.failure(advanced_result.issues)

    normalized = _normalize_inputs(standard, advanced_result.value, inputs.samples)
    result_identity_issue = _validate_normalized_result_identities(normalized)
    if result_identity_issue is not None:
        return Result.failure([result_identity_issue])
    return Result.success(normalized)


def _validate_standard_semantics(
    standard: Mapping[str, object],
    samples: list[dict[str, str]],
) -> Issue | None:
    reference = standard["reference"]
    if reference["fasta"] == reference["gtf"]:
        return _issue(
            "BULK_RNASEQ_REFERENCE_CONFLICT",
            "config.standard.reference",
        )
    if reference["fasta_sha256"] == reference["gtf_sha256"]:
        return _issue(
            "BULK_RNASEQ_REFERENCE_CONFLICT",
            "config.standard.reference",
        )
    if "salmon_index" in reference and not any(
        row["strandedness"] == "auto" for row in samples
    ):
        return _issue(
            "BULK_RNASEQ_REFERENCE_CONTEXT_CONFLICT",
            "config.standard.reference.salmon_index",
        )

    umi = standard.get("umi", {"enabled": False})
    enabled = umi["enabled"]
    keys = set(umi)
    if not enabled and keys != {"enabled"}:
        return _issue("BULK_RNASEQ_UMI_CONFLICT", "config.standard.umi")
    if enabled and umi["mode"] == "read_name":
        incompatible = {
            "extraction_method",
            "barcode_pattern",
            "barcode_pattern_2",
            "discard_read",
        }
        if keys.intersection(incompatible):
            return _issue("BULK_RNASEQ_UMI_CONFLICT", "config.standard.umi")
    if enabled and umi["mode"] == "read_sequence" and "read_name_separator" in keys:
        return _issue("BULK_RNASEQ_UMI_CONFLICT", "config.standard.umi")
    if enabled and umi["mode"] == "read_sequence":
        patterns = [umi["barcode_pattern"]]
        if "barcode_pattern_2" in umi:
            patterns.append(umi["barcode_pattern_2"])
        if umi["extraction_method"] == "string":
            if any(re.fullmatch(r"[NCX]+", pattern) is None for pattern in patterns):
                return _issue(
                    "BULK_RNASEQ_UMI_PATTERN_INVALID",
                    "config.standard.umi",
                )
        elif any(not _valid_umi_regex(pattern) for pattern in patterns):
            return _issue(
                "BULK_RNASEQ_UMI_PATTERN_INVALID",
                "config.standard.umi",
            )
        if ("barcode_pattern_2" in umi or "discard_read" in umi) and any(
            row["layout"] != "PE" for row in samples
        ):
            return _issue(
                "BULK_RNASEQ_UMI_LAYOUT_CONFLICT",
                "config.standard.umi",
            )

    ribosomal_rna_removal = standard.get(
        "ribosomal_rna_removal",
        RIBOSOMAL_RNA_REMOVAL_DEFAULTS,
    )
    ribosomal_rna_enabled = ribosomal_rna_removal["enabled"]
    if not ribosomal_rna_enabled and set(ribosomal_rna_removal) != {"enabled"}:
        return _issue(
            "BULK_RNASEQ_RRNA_CONFLICT",
            "config.standard.ribosomal_rna_removal",
        )
    if (
        ribosomal_rna_enabled
        and ribosomal_rna_removal["tool"] != "sortmerna"
        and "sortmerna_index" in ribosomal_rna_removal
    ):
        return _issue(
            "BULK_RNASEQ_RRNA_CONFLICT",
            "config.standard.ribosomal_rna_removal.sortmerna_index",
        )
    if (
        ribosomal_rna_enabled
        and ribosomal_rna_removal["tool"] == "bowtie2"
        and ribosomal_rna_removal.get("save_filtered_reads", False)
        and any(row["layout"] != "SE" for row in samples)
    ):
        return _issue(
            "BULK_RNASEQ_RRNA_CONFLICT",
            "config.standard.ribosomal_rna_removal.save_filtered_reads",
        )

    qc = standard.get("qc", {})
    if qc.get("enabled") is False and any(
        value is True for name, value in qc.items() if name != "enabled"
    ):
        return _issue("BULK_RNASEQ_QC_CONFLICT", "config.standard.qc")

    trimming = {**TRIMMING_DEFAULTS, **standard.get("trimming", {})}
    outputs = standard.get("outputs", {})
    if outputs.get("umi_intermediates") is True and not enabled:
        return _issue("BULK_RNASEQ_OUTPUT_CONFLICT", "config.standard.outputs")
    if outputs.get("trimmed_reads") is True and not trimming["enabled"]:
        return _issue("BULK_RNASEQ_OUTPUT_CONFLICT", "config.standard.outputs")
    return None


def _validate_sample_semantics(samples: list[dict[str, str]]) -> Issue | None:
    sample_semantics: dict[str, tuple[str, str, str]] = {}
    coordinates: set[tuple[str, str, str]] = set()
    fastqs: set[str] = set()
    for index, row in enumerate(samples):
        sample = row["sample"]
        semantics = (row["layout"], row["strandedness"], row["platform"])
        previous = sample_semantics.setdefault(sample, semantics)
        if previous != semantics:
            return _issue(
                "BULK_RNASEQ_LANE_SEMANTICS_CONFLICT",
                f"samples[{index}]",
            )
        coordinate = (sample, row["library"], row["lane"])
        if coordinate in coordinates:
            return _issue("BULK_RNASEQ_LANE_DUPLICATE", f"samples[{index}]")
        coordinates.add(coordinate)
        row_fastqs = [row["fastq_1"]]
        if row.get("fastq_2"):
            row_fastqs.append(row["fastq_2"])
        for fastq in row_fastqs:
            if fastq in fastqs:
                return _issue("BULK_RNASEQ_FASTQ_DUPLICATE", f"samples[{index}]")
            fastqs.add(fastq)
    return None


def _validate_advanced(
    advanced: Mapping[str, object],
    standard: Mapping[str, object],
    samples: list[dict[str, str]],
) -> Result[dict[str, object]]:
    upstream = upstream_parameter_properties()
    for name in sorted(advanced):
        if name in DANGEROUS_RUNTIME_PARAMETERS:
            return _failure(
                "BULK_RNASEQ_PLATFORM_PARAMETER_FORBIDDEN",
                f"config.advanced.{name}",
            )
        if name not in KNOWN_UPSTREAM_PARAMETERS:
            return _failure(
                "BULK_RNASEQ_ADVANCED_UNKNOWN",
                "config.advanced",
            )
        if name in STANDARD_NATIVE_PARAMETERS:
            return _failure(
                "BULK_RNASEQ_PARAMETER_SURFACE_CONFLICT",
                f"config.advanced.{name}",
            )
        if name in PLATFORM_OWNED_NATIVE_PARAMETERS:
            return _failure(
                "BULK_RNASEQ_PLATFORM_PARAMETER_FORBIDDEN",
                f"config.advanced.{name}",
            )
        if name in RAW_ARGUMENT_NATIVE_PARAMETERS:
            return _failure(
                "BULK_RNASEQ_RAW_ARGUMENT_FORBIDDEN",
                f"config.advanced.{name}",
            )
        if name in UNSUPPORTED_NATIVE_PARAMETERS:
            return _failure(
                "BULK_RNASEQ_ADVANCED_UNSUPPORTED",
                f"config.advanced.{name}",
            )
        if name not in ADVANCED_NATIVE_PARAMETERS:
            return _failure(
                "BULK_RNASEQ_ADVANCED_UNSUPPORTED",
                f"config.advanced.{name}",
            )
        exact_schema = {
            "$schema": "https://json-schema.org/draft/2020-12/schema",
            **upstream[name],
        }
        if _first_schema_error(exact_schema, advanced[name]) is not None:
            return _failure(
                "BULK_RNASEQ_ADVANCED_INVALID",
                f"config.advanced.{name}",
            )

    public_advanced_schema = build_bulk_rnaseq_authoring_schema().to_dict()[
        "config_schema"
    ]["properties"]["advanced"]
    projected_error = _first_schema_error(public_advanced_schema, advanced)
    if projected_error is not None:
        return _failure(
            "BULK_RNASEQ_ADVANCED_INVALID",
            _schema_path("config.advanced", projected_error),
        )

    semantic_issue = _validate_advanced_semantics(advanced, standard, samples)
    if semantic_issue is not None:
        return Result.failure([semantic_issue])
    return Result.success({name: deepcopy(advanced[name]) for name in sorted(advanced)})


def _validate_advanced_semantics(
    advanced: Mapping[str, object],
    standard: Mapping[str, object],
    samples: list[dict[str, str]],
) -> Issue | None:
    trimming = {**TRIMMING_DEFAULTS, **standard.get("trimming", {})}
    qc = {**QC_DEFAULTS, **standard.get("qc", {})}
    outputs = {**OUTPUT_DEFAULTS, **standard.get("outputs", {})}
    reference = standard["reference"]

    if "min_trimmed_reads" in advanced and not trimming["enabled"]:
        return _issue(
            "BULK_RNASEQ_ADVANCED_CONTEXT_CONFLICT",
            "config.advanced.min_trimmed_reads",
        )
    if "star_ignore_sjdbgtf" in advanced and "star_index" not in reference:
        return _issue(
            "BULK_RNASEQ_ADVANCED_CONTEXT_CONFLICT",
            "config.advanced.star_ignore_sjdbgtf",
        )
    if "stringtie_ignore_gtf" in advanced and not outputs["stringtie"]:
        return _issue(
            "BULK_RNASEQ_ADVANCED_CONTEXT_CONFLICT",
            "config.advanced.stringtie_ignore_gtf",
        )
    if "rseqc_modules" in advanced:
        if not qc["enabled"] or not qc["rseqc"]:
            return _issue(
                "BULK_RNASEQ_ADVANCED_CONTEXT_CONFLICT",
                "config.advanced.rseqc_modules",
            )
        modules = str(advanced["rseqc_modules"]).split(",")
        if len(modules) != len(set(modules)) or not set(modules).issubset(
            _RSEQC_MODULES
        ):
            return _issue(
                "BULK_RNASEQ_ADVANCED_INVALID",
                "config.advanced.rseqc_modules",
            )
    if "deseq2_vst" in advanced and (not qc["enabled"] or not qc["deseq2_pca"]):
        return _issue(
            "BULK_RNASEQ_ADVANCED_CONTEXT_CONFLICT",
            "config.advanced.deseq2_vst",
        )
    biotype_parameters = {
        "featurecounts_group_type",
        "featurecounts_feature_type",
    }
    if biotype_parameters.intersection(advanced) and (
        not qc["enabled"] or not qc["biotype"]
    ):
        return _issue(
            "BULK_RNASEQ_ADVANCED_CONTEXT_CONFLICT",
            "config.advanced",
        )
    if (
        "featurecounts_group_type" in advanced
        and reference.get("annotation_style", "ensembl") == "gencode"
    ):
        return _issue(
            "BULK_RNASEQ_ADVANCED_CONTEXT_CONFLICT",
            "config.advanced.featurecounts_group_type",
        )
    threshold_names = {"stranded_threshold", "unstranded_threshold"}
    if threshold_names.intersection(advanced) and not any(
        row["strandedness"] == "auto" for row in samples
    ):
        return _issue(
            "BULK_RNASEQ_ADVANCED_CONTEXT_CONFLICT",
            "config.advanced",
        )
    stranded = advanced.get("stranded_threshold", 0.8)
    unstranded = advanced.get("unstranded_threshold", 0.1)
    if unstranded >= stranded:
        return _issue(
            "BULK_RNASEQ_ADVANCED_INVALID",
            "config.advanced",
        )
    return None


def _normalize_inputs(
    standard: Mapping[str, object],
    advanced: Mapping[str, object],
    samples: list[dict[str, str]],
) -> dict[str, object]:
    analysis = {**ANALYSIS_DEFAULTS, **standard.get("analysis", {})}
    trimming = {**TRIMMING_DEFAULTS, **standard.get("trimming", {})}
    qc = {**QC_DEFAULTS, **standard.get("qc", {})}
    outputs = {**OUTPUT_DEFAULTS, **standard.get("outputs", {})}
    reference = standard["reference"]
    umi = standard.get("umi", {"enabled": False})
    ribosomal_rna_removal = standard.get(
        "ribosomal_rna_removal",
        RIBOSOMAL_RNA_REMOVAL_DEFAULTS,
    )

    params: dict[str, object] = {
        "aligner": f"{analysis['alignment']}_{analysis['quantification']}",
        "fasta": reference["fasta"],
        "gencode": reference.get("annotation_style", "ensembl") == "gencode",
        "gtf": reference["gtf"],
        "igenomes_ignore": True,
        "remove_ribo_rna": ribosomal_rna_removal["enabled"],
        "save_align_intermeds": outputs["alignment_intermediates"],
        "save_merged_fastq": outputs["merged_fastq"],
        "save_non_ribo_reads": ribosomal_rna_removal.get("save_filtered_reads", False),
        "save_reference": outputs["reference_files"],
        "save_trimmed": outputs["trimmed_reads"],
        "save_umi_intermeds": outputs["umi_intermediates"],
        "save_unaligned": outputs["unaligned_reads"],
        "skip_alignment": False,
        "skip_bbsplit": True,
        "skip_bigwig": not outputs["bigwig"],
        "skip_biotype_qc": not (qc["enabled"] and qc["biotype"]),
        "skip_deseq2_qc": not (qc["enabled"] and qc["deseq2_pca"]),
        "skip_dupradar": not (qc["enabled"] and qc["dupradar"]),
        "skip_fastqc": not (qc["enabled"] and qc["fastqc"]),
        "skip_linting": False,
        "skip_markduplicates": not (
            qc["enabled"] and qc["mark_duplicates"] and not umi["enabled"]
        ),
        "skip_multiqc": not (qc["enabled"] and qc["multiqc"]),
        "skip_preseq": not (qc["enabled"] and qc["preseq"]),
        "skip_pseudo_alignment": True,
        "skip_qc": not qc["enabled"],
        "skip_qualimap": not (qc["enabled"] and qc["qualimap"]),
        "skip_quantification_merge": False,
        "skip_rseqc": not (qc["enabled"] and qc["rseqc"]),
        "skip_stringtie": not outputs["stringtie"],
        "skip_trimming": not trimming["enabled"],
        "trimmer": trimming["tool"],
        "with_umi": umi["enabled"],
    }
    if "star_index" in reference:
        params["star_index"] = reference["star_index"]["path"]
    if "salmon_index" in reference:
        params["salmon_index"] = reference["salmon_index"]["path"]
    if ribosomal_rna_removal["enabled"]:
        params["ribo_removal_tool"] = ribosomal_rna_removal["tool"]
        params["ribo_database_manifest"] = ribosomal_rna_removal["database_manifest"][
            "path"
        ]
        if "sortmerna_index" in ribosomal_rna_removal:
            params["sortmerna_index"] = ribosomal_rna_removal["sortmerna_index"]["path"]
    if umi["enabled"]:
        params.update(_normalize_umi(umi))
    params.update(advanced)
    ordered_params = {name: deepcopy(params[name]) for name in sorted(params)}
    impacts = [
        {"native_parameter": name, "impact": PARAMETER_IMPACT_BY_NATIVE_NAME[name]}
        for name in sorted(params)
        if name in PARAMETER_IMPACT_BY_NATIVE_NAME
    ]

    normalized_samples = []
    for row in samples:
        normalized_samples.append(
            {
                "sample": row["sample"],
                "library": row["library"],
                "lane": row["lane"],
                "layout": row["layout"],
                "fastq_1": row["fastq_1"],
                "fastq_2": row.get("fastq_2", ""),
                "strandedness": row["strandedness"],
                "seq_platform": "ILLUMINA",
            }
        )

    return {
        "contract": {
            "workflow_id": _WORKFLOW_ID,
            "schema_version": SCHEMA_VERSION,
            "nfcore_rnaseq_release": NFCORE_RNASEQ_RELEASE,
            "nfcore_rnaseq_commit": NFCORE_RNASEQ_COMMIT,
            "upstream_parameter_schema_sha256": UPSTREAM_PARAMETER_SCHEMA_SHA256,
            "upstream_samplesheet_schema_sha256": (UPSTREAM_SAMPLESHEET_SCHEMA_SHA256),
        },
        "reference_identity": {
            "reference_id": reference["reference_id"],
            "fasta_sha256": reference["fasta_sha256"],
            "gtf_sha256": reference["gtf_sha256"],
            **(
                {"star_index_sha256": reference["star_index"]["identity_sha256"]}
                if "star_index" in reference
                else {}
            ),
            **(
                {"salmon_index_sha256": reference["salmon_index"]["identity_sha256"]}
                if "salmon_index" in reference
                else {}
            ),
            **(
                {
                    "ribo_database_manifest_sha256": ribosomal_rna_removal[
                        "database_manifest"
                    ]["identity_sha256"]
                }
                if ribosomal_rna_removal["enabled"]
                else {}
            ),
            **(
                {
                    "sortmerna_index_sha256": ribosomal_rna_removal["sortmerna_index"][
                        "identity_sha256"
                    ]
                }
                if "sortmerna_index" in ribosomal_rna_removal
                else {}
            ),
        },
        "nfcore_params": ordered_params,
        "parameter_impacts": impacts,
        "samples": normalized_samples,
    }


def _normalize_umi(umi: Mapping[str, object]) -> dict[str, object]:
    params: dict[str, object] = {
        "umi_dedup_tool": umi["deduplication_tool"],
        "umitools_dedup_primary_only": umi.get("primary_alignments_only", False),
        "umitools_dedup_stats": umi.get("emit_dedup_stats", False),
        "umitools_grouping_method": umi.get("grouping_method", "directional"),
    }
    if umi["mode"] == "read_name":
        params["skip_umi_extract"] = True
        params["umitools_umi_separator"] = umi["read_name_separator"]
    else:
        params["skip_umi_extract"] = False
        params["umitools_extract_method"] = umi["extraction_method"]
        params["umitools_bc_pattern"] = umi["barcode_pattern"]
        if "barcode_pattern_2" in umi:
            params["umitools_bc_pattern2"] = umi["barcode_pattern_2"]
        if "discard_read" in umi:
            params["umi_discard_read"] = umi["discard_read"]
    return params


def _validate_normalized_result_identities(
    normalized: Mapping[str, object],
) -> Issue | None:
    params = normalized["nfcore_params"]
    samples = normalized["samples"]
    if not isinstance(params, Mapping) or not isinstance(samples, list):
        return _issue("BULK_RNASEQ_INPUTS_INVALID", "inputs")

    if params.get("skip_multiqc") is False:
        multiqc_issue = _validate_multiqc_sample_identities(samples, params)
        if multiqc_issue is not None:
            return multiqc_issue

    if (
        params.get("skip_rseqc") is False
        and "tin" in effective_rseqc_modules(params)
        and any("bam" in sample for sample in _canonical_sample_ids(samples))
    ):
        return _issue(
            "BULK_RNASEQ_RSEQC_TIN_SAMPLE_ID_UNSUPPORTED",
            "samples",
            context={
                "module": "tin",
                "reason_code": "upstream_filename_rewrite",
                "tool": "rseqc",
                "tool_version": "5.0.4",
            },
        )
    return None


def _validate_multiqc_sample_identities(
    samples: list[object],
    params: Mapping[str, object],
) -> Issue | None:
    rows_by_sample: dict[str, list[Mapping[str, object]]] = {}
    layouts: dict[str, str] = {}
    for row in samples:
        if not isinstance(row, Mapping):
            return _issue("BULK_RNASEQ_INPUTS_INVALID", "inputs")
        sample = row.get("sample")
        layout = row.get("layout")
        if not isinstance(sample, str) or layout not in {"SE", "PE"}:
            return _issue("BULK_RNASEQ_INPUTS_INVALID", "inputs")
        layouts[sample] = layout
        rows_by_sample.setdefault(sample, []).append(row)

    identity_owners: dict[str, tuple[str, str]] = {}
    for sample in sorted(rows_by_sample):
        if _multiqc_clean_identity(sample) != sample:
            return _multiqc_identity_issue()
        layout = layouts[sample]
        downstream_layout = effective_downstream_layout(layout, params)
        claims: list[tuple[str, str]] = [(sample, "sample")]

        # The pinned nf-core MultiQC subworkflow builds run-global table merge
        # aliases for every authored PE sample, independently of UMI discard.
        if layout == "PE":
            claims.extend(((f"{sample}_1", "read1"), (f"{sample}_2", "read2")))

        # Raw FastQC follows the authored layout. Trim Galore, trimmed FastQC,
        # and filtered FastQC follow the effective post-UMI layout.
        if params.get("skip_fastqc") is False and layout == "PE":
            claims.extend(((f"{sample}_1", "read1"), (f"{sample}_2", "read2")))
        downstream_mate_identities = downstream_layout == "PE" and (
            (
                params.get("skip_trimming") is False
                and params.get("trimmer") == "trimgalore"
            )
            or (
                params.get("skip_fastqc") is False
                and (
                    params.get("skip_trimming") is False
                    or params.get("remove_ribo_rna") is True
                )
            )
        )
        if downstream_mate_identities:
            claims.extend(
                (
                    (f"{sample}_1", "read1"),
                    (f"{sample}_2", "read2"),
                )
            )

        for identity, role in claims:
            owner = (sample, role)
            previous = identity_owners.setdefault(identity, owner)
            if previous != owner:
                return _multiqc_identity_issue()

    # The merged report also supplies exact FASTQ simpleName replacements.
    # Replicate Nextflow Path.simpleName (the basename before its first dot)
    # and reject source aliases that can rename another sample or map to two
    # biological owners/read roles. Considering every lane is conservative
    # and removes input-order dependence from nf-core's first-lane selection.
    replacements: dict[str, tuple[str, tuple[str, str]]] = {}
    observed_sources: list[tuple[str, tuple[str, str], str | None, bool]] = []
    for sample in sorted(rows_by_sample):
        layout = layouts[sample]
        for row in rows_by_sample[sample]:
            read1 = _nextflow_simple_name(row.get("fastq_1"))
            emits_replacements = read1 != sample
            sources = [(read1, "single" if layout == "SE" else "read1")]
            if layout == "PE":
                sources.append((_nextflow_simple_name(row.get("fastq_2")), "read2"))
            for index, (source_identity, role) in enumerate(sources):
                if not source_identity:
                    return _multiqc_identity_issue()
                owner = (sample, role)
                target_identity = (
                    sample
                    if role == "single"
                    else f"{sample}_{1 if role == 'read1' else 2}"
                )
                cleaned_identity = _multiqc_clean_identity(source_identity)
                if not emits_replacements and cleaned_identity != source_identity:
                    return _multiqc_identity_issue()
                expected_target: str | None = None
                if emits_replacements:
                    expected_target = target_identity
                    replacement = (target_identity, owner)
                    previous_replacement = replacements.setdefault(
                        source_identity, replacement
                    )
                    if previous_replacement != replacement:
                        return _multiqc_identity_issue()
                observed_sources.append(
                    (
                        cleaned_identity,
                        owner,
                        expected_target,
                        not emits_replacements and index == 0,
                    )
                )

    # A hard replacement is global. Its raw key must not already denote a
    # canonical/derived identity owned by another sample or mate role.
    for source_identity, (_target_identity, target_owner) in replacements.items():
        source_owner = identity_owners.get(source_identity)
        if source_owner is not None and source_owner != target_owner:
            return _multiqc_identity_issue()

    # MultiQC cleans every observed identity before looking it up in the
    # global exact-replacement map. Resolve that final graph before accepting
    # any owner; this also crosses upstream no-op rows with mappings emitted by
    # other samples.
    final_owners = dict(identity_owners)
    for cleaned_identity, owner, expected_target, canonical_noop in observed_sources:
        resolved_replacement = replacements.get(cleaned_identity)
        if resolved_replacement is None:
            final_identity = cleaned_identity
        else:
            final_identity, replacement_owner = resolved_replacement
            if replacement_owner != owner:
                return _multiqc_identity_issue()
        if expected_target is not None and final_identity != expected_target:
            return _multiqc_identity_issue()
        previous_owner = final_owners.setdefault(final_identity, owner)
        if previous_owner == owner:
            continue
        canonical_owner = (owner[0], "sample")
        if previous_owner == canonical_owner and (
            owner[1] == "single" or canonical_noop
        ):
            continue
        return _multiqc_identity_issue()
    return None


def _canonical_sample_ids(samples: list[object]) -> tuple[str, ...]:
    return tuple(
        sorted(
            {
                sample
                for row in samples
                if isinstance(row, Mapping)
                and isinstance((sample := row.get("sample")), str)
            }
        )
    )


def _nextflow_simple_name(value: object) -> str:
    if not isinstance(value, str):
        return ""
    name = PurePosixPath(value).name
    return name.split(".", 1)[0]


def _multiqc_clean_identity(value: str) -> str:
    original = value
    cleaned = value
    for token in _MULTIQC_EXTRA_TRUNCATE_TOKENS:
        cleaned = cleaned.split(token, 1)[0]
    for token in _MULTIQC_EXTRA_REMOVE_TOKENS:
        cleaned = cleaned.replace(token, "")
    for token in MULTIQC_SAMPLE_CLEAN_TOKENS:
        if token not in _MULTIQC_EXTRA_TOKENS:
            cleaned = cleaned.split(token, 1)[0]
    for token in _MULTIQC_FN_CLEAN_TRIM:
        if cleaned.endswith(token):
            cleaned = cleaned[: -len(token)]
        if cleaned.startswith(token):
            cleaned = cleaned[len(token) :]
    cleaned = cleaned.strip()
    return original if not cleaned else cleaned


def _multiqc_identity_issue() -> Issue:
    return _issue(
        "BULK_RNASEQ_MULTIQC_SAMPLE_IDENTITY_CONFLICT",
        "samples",
        context={
            "identity_space": "multiqc-1.33",
            "reason_code": "derived_sample_identity_not_unique",
        },
    )


def _first_schema_error(
    schema: Mapping[str, object],
    value: object,
) -> jsonschema.ValidationError | None:
    validator = jsonschema.Draft202012Validator(schema)
    errors = sorted(
        validator.iter_errors(value),
        key=lambda error: (
            tuple(str(part) for part in error.absolute_path),
            str(error.validator),
        ),
    )
    return errors[0] if errors else None


def _schema_path(prefix: str, error: jsonschema.ValidationError) -> str:
    suffix = "".join(
        f"[{part}]" if isinstance(part, int) else f".{part}"
        for part in error.absolute_path
    )
    return f"{prefix}{suffix}"


def _issue(
    code: str,
    path: str,
    *,
    context: Mapping[str, object] | None = None,
) -> Issue:
    return Issue(
        code=code,
        message="Bulk RNA-seq input validation failed.",
        severity="error",
        path=path,
        source="adapter",
        context={} if context is None else dict(context),
    )


def _failure(code: str, path: str) -> Result[Any]:
    return Result.failure([_issue(code, path)])


def _is_json_value(value: object) -> bool:
    if value is None or isinstance(value, (str, bool, int)):
        return True
    if isinstance(value, float):
        return math.isfinite(value)
    if isinstance(value, list):
        return all(_is_json_value(item) for item in value)
    if isinstance(value, Mapping):
        return all(
            isinstance(key, str) and _is_json_value(item) for key, item in value.items()
        )
    return False


def _valid_umi_regex(pattern: object) -> bool:
    if (
        not isinstance(pattern, str)
        or "(?P<umi" not in pattern
        or _SAFE_UMI_REGEX.fullmatch(pattern) is None
    ):
        return False
    try:
        re.compile(pattern)
    except re.error:
        return False
    return True
