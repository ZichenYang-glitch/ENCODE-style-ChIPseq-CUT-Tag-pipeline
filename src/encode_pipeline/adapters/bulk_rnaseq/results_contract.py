"""Immutable identity loader for the pinned bulk RNA-seq results contract."""

from __future__ import annotations

from copy import deepcopy
from hashlib import sha256
from importlib import resources
import json
from typing import Any, Mapping


RESULTS_CONTRACT_FILE = "results-contract-3.26.0.json"
RESULTS_CONTRACT_SIZE = 26_889
RESULTS_CONTRACT_SHA256 = (
    "a099f495de027a385580cf14e2316b7e96ea6d67b971af02f3180cbf751e859f"
)
_CONTRACT_PACKAGE = "encode_pipeline.contracts.nfcore_rnaseq"
_DEFAULT_RSEQC_MODULES = (
    "bam_stat",
    "inner_distance",
    "infer_experiment",
    "junction_annotation",
    "junction_saturation",
    "read_distribution",
    "read_duplication",
)
_CSI_INCOMPATIBLE_RSEQC_MODULES = frozenset(
    {"inner_distance", "read_distribution", "tin"}
)

ARTIFACT_FAMILY_OUTPUT_TYPES: Mapping[str, tuple[str, ...]] = {
    "pipeline_software_versions": ("bulk_rnaseq.pipeline.*",),
    "fastqc_html_zip": ("bulk_rnaseq.fastqc.*",),
    "trimming_machine_reports": ("bulk_rnaseq.trim.fastp.json",),
    "saved_trimmed_reads": (
        "bulk_rnaseq.trim.fastp.filtered.*",
        "bulk_rnaseq.trim.trimgalore.filtered.*",
    ),
    "saved_merged_reads": ("bulk_rnaseq.fastq.merged.*",),
    "saved_rrna_filtered_reads": ("bulk_rnaseq.rrna.*",),
    "star_final_summary_and_junctions": (
        "bulk_rnaseq.star.log_final",
        "bulk_rnaseq.star.splice_junctions",
    ),
    "final_and_enabled_intermediate_bam_indexes": (
        "bulk_rnaseq.star.bam*",
        "bulk_rnaseq.star.sorted_intermediate_bam*",
        "bulk_rnaseq.star.transcriptome_sorted_intermediate_bam*",
        "bulk_rnaseq.star.aligned_genome_bam",
        "bulk_rnaseq.star.aligned_transcriptome_bam",
    ),
    "star_unaligned_reads": ("bulk_rnaseq.star.unaligned.*",),
    "samtools_machine_summaries": ("bulk_rnaseq.samtools.*",),
    "umi_intermediates_and_statistics": ("bulk_rnaseq.umi.*",),
    "combined_and_strand_specific_bigwig": ("bulk_rnaseq.star.bigwig.*",),
    "salmon_sample_quantification": (
        "bulk_rnaseq.salmon.quant_transcript",
        "bulk_rnaseq.salmon.quant_gene",
        "bulk_rnaseq.salmon.meta_info",
    ),
    "salmon_merged_expression_matrices": ("bulk_rnaseq.salmon.merged.*",),
    "featurecounts_biotype_qc": ("bulk_rnaseq.featurecounts.*",),
    "enabled_rseqc_machine_summaries": ("bulk_rnaseq.rseqc.*",),
    "multiqc_allowlisted_machine_tables": ("bulk_rnaseq.multiqc.*",),
}


def effective_downstream_layout(
    original_layout: str,
    params: Mapping[str, object],
) -> str:
    """Return the fixed post-UMI-extraction layout used by downstream tools."""
    if original_layout not in {"SE", "PE"}:
        raise ValueError("invalid bulk RNA-seq sample layout")
    if (
        original_layout == "PE"
        and params.get("with_umi") is True
        and params.get("skip_umi_extract") is False
        and type(params.get("umi_discard_read")) is int
        and params.get("umi_discard_read") in {1, 2}
    ):
        return "SE"
    return original_layout


def trimmed_fastqc_enabled(params: Mapping[str, object]) -> bool:
    """Return whether the pinned route emits post-trim FastQC outputs."""
    if params.get("skip_trimming") is True:
        return False
    if params.get("trimmer") == "trimgalore":
        # nf-core/rnaseq 3.26.0 always supplies Trim Galore ``--fastqc_args``;
        # ``skip_fastqc`` only controls the separate raw/filtered FASTQC route.
        return True
    return params.get("trimmer") == "fastp" and params.get("skip_fastqc") is False


def effective_rseqc_modules(params: Mapping[str, object]) -> tuple[str, ...]:
    """Mirror the pinned nf-core RSeQC module removal for CSI-indexed BAMs."""
    configured = params.get("rseqc_modules")
    modules = (
        tuple(str(configured).split(","))
        if configured is not None
        else _DEFAULT_RSEQC_MODULES
    )
    if params.get("bam_csi_index") is True:
        modules = tuple(
            module
            for module in modules
            if module not in _CSI_INCOMPATIBLE_RSEQC_MODULES
        )
    return modules


def artifact_family_for_output_type(output_type: str) -> str:
    """Return the sole declared public family for one closed output type."""
    if not isinstance(output_type, str):
        raise ValueError("invalid bulk RNA-seq artifact output type")
    matches = [
        family
        for family, patterns in ARTIFACT_FAMILY_OUTPUT_TYPES.items()
        if any(
            output_type.startswith(pattern[:-1])
            if pattern.endswith("*")
            else output_type == pattern
            for pattern in patterns
        )
    ]
    if len(matches) != 1:
        raise ValueError("invalid bulk RNA-seq artifact family membership")
    return matches[0]


def load_bulk_rnaseq_results_contract() -> dict[str, Any]:
    """Return a fresh verified copy of the exact results contract document."""
    content = (
        resources.files(_CONTRACT_PACKAGE).joinpath(RESULTS_CONTRACT_FILE).read_bytes()
    )
    if (
        len(content) != RESULTS_CONTRACT_SIZE
        or sha256(content).hexdigest() != RESULTS_CONTRACT_SHA256
    ):
        raise ValueError("bulk RNA-seq results contract identity mismatch")
    try:
        value = json.loads(content)
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError("bulk RNA-seq results contract is invalid") from exc
    if (
        not isinstance(value, dict)
        or value.get("schema_version") != "1.0.0"
        or value.get("workflow_id") != "bulk-rnaseq"
        or not isinstance(value.get("upstream"), dict)
        or value["upstream"].get("release") != "3.26.0"
        or value["upstream"].get("commit") != "e7ca46272c8f9d5ceee3f71759f4ba551d3217a4"
        or value["upstream"].get("route") != "star_salmon"
    ):
        raise ValueError("bulk RNA-seq results contract is invalid")
    raw_membership = value.get("artifact_family_output_types")
    expected_membership = {
        family: list(patterns)
        for family, patterns in ARTIFACT_FAMILY_OUTPUT_TYPES.items()
    }
    if (
        value.get("artifact_families") != list(ARTIFACT_FAMILY_OUTPUT_TYPES)
        or raw_membership != expected_membership
    ):
        raise ValueError("bulk RNA-seq results contract family catalog is invalid")
    return deepcopy(value)
