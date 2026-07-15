"""Pinned nf-core/rnaseq identity and native-parameter policy."""

from __future__ import annotations

from copy import deepcopy
import hashlib
from importlib import resources
import json
from typing import Any


NFCORE_RNASEQ_RELEASE = "3.26.0"
NFCORE_RNASEQ_COMMIT = "e7ca46272c8f9d5ceee3f71759f4ba551d3217a4"
NFCORE_RNASEQ_LICENSE = "MIT"

UPSTREAM_PARAMETER_SCHEMA_FILE = "nextflow-schema-3.26.0.json"
UPSTREAM_PARAMETER_SCHEMA_SIZE = 60_620
UPSTREAM_PARAMETER_SCHEMA_SHA256 = (
    "8f2f84a25c0aec65a18234cf01acdd74f2385e8dfac8417e4bad23a70bfb4388"
)
UPSTREAM_SAMPLESHEET_SCHEMA_FILE = "samplesheet-schema-3.26.0.json"
UPSTREAM_SAMPLESHEET_SCHEMA_SIZE = 3_218
UPSTREAM_SAMPLESHEET_SCHEMA_SHA256 = (
    "013b669b1a3d38709f548a2548350ecdabb3f2ba578f2b7de843105e2ce87a7d"
)
UPSTREAM_LICENSE_FILE = "LICENSE-3.26.0.txt"
UPSTREAM_LICENSE_SIZE = 1_075
UPSTREAM_LICENSE_SHA256 = (
    "49147db9f580fc819494d436920fd8db97589e67989d8a82d5d18a4deb8d49f9"
)

_CONTRACT_PACKAGE = "encode_pipeline.contracts.nfcore_rnaseq"

# Stable product fields own these native targets. Supplying them through the
# Advanced surface is always a conflict, even when the values happen to match.
STANDARD_NATIVE_PARAMETERS = frozenset(
    {
        "aligner",
        "fasta",
        "gencode",
        "gtf",
        "salmon_index",
        "salmon_quant_libtype",
        "save_align_intermeds",
        "save_merged_fastq",
        "save_reference",
        "save_trimmed",
        "save_umi_intermeds",
        "save_unaligned",
        "skip_bigwig",
        "skip_biotype_qc",
        "skip_deseq2_qc",
        "skip_dupradar",
        "skip_fastqc",
        "skip_markduplicates",
        "skip_multiqc",
        "skip_preseq",
        "skip_qc",
        "skip_qualimap",
        "skip_rseqc",
        "skip_stringtie",
        "skip_trimming",
        "skip_umi_extract",
        "star_index",
        "trimmer",
        "umi_dedup_tool",
        "umi_discard_read",
        "umitools_bc_pattern",
        "umitools_bc_pattern2",
        "umitools_dedup_primary_only",
        "umitools_dedup_stats",
        "umitools_extract_method",
        "umitools_grouping_method",
        "umitools_umi_separator",
        "with_umi",
    }
)

# This is deliberately small. Every value is validated against the exact
# upstream property schema and then against adapter-owned semantic constraints.
ADVANCED_NATIVE_PARAMETERS = frozenset(
    {
        "bam_csi_index",
        "deseq2_vst",
        "featurecounts_feature_type",
        "featurecounts_group_type",
        "gffread_transcript_fasta",
        "gtf_extra_attributes",
        "gtf_group_features",
        "min_mapped_reads",
        "min_trimmed_reads",
        "rseqc_modules",
        "skip_gtf_filter",
        "skip_gtf_transcript_filter",
        "star_ignore_sjdbgtf",
        "stranded_threshold",
        "stringtie_ignore_gtf",
        "unstranded_threshold",
    }
)

# Paths, runtime policy, remote configuration, notification, hardware, and
# publication controls are never user-owned native parameters.
PLATFORM_OWNED_NATIVE_PARAMETERS = frozenset(
    {
        "additional_fasta",
        "arm",
        "bbsplit_fasta_list",
        "bbsplit_index",
        "bowtie2_index",
        "config_profile_contact",
        "config_profile_description",
        "config_profile_name",
        "config_profile_url",
        "custom_config_base",
        "custom_config_version",
        "email",
        "email_on_fail",
        "gene_bed",
        "genome",
        "gff",
        "gpu_container_options",
        "help",
        "help_full",
        "hisat2_build_memory",
        "hisat2_index",
        "igenomes_base",
        "igenomes_ignore",
        "input",
        "kallisto_index",
        "kraken_db",
        "max_multiqc_email_size",
        "monochrome_logs",
        "multiqc_config",
        "multiqc_logo",
        "multiqc_methods_description",
        "multiqc_title",
        "outdir",
        "pipelines_testdata_base_path",
        "plaintext_email",
        "prokaryotic",
        "pseudo_aligner",
        "publish_dir_mode",
        "ribo_database_manifest",
        "rsem_index",
        "seq_center",
        "seq_platform",
        "show_hidden",
        "skip_alignment",
        "skip_bbsplit",
        "skip_linting",
        "skip_pseudo_alignment",
        "skip_quantification_merge",
        "sortmerna_index",
        "splicesites",
        "sylph_db",
        "sylph_taxonomy",
        "trace_report_suffix",
        "transcript_fasta",
        "use_gpu_ribodetector",
        "use_parabricks_star",
        "use_sentieon_star",
        "validate_params",
        "version",
    }
)

# These values are strings forwarded to underlying command-line tools. The
# contract never accepts arbitrary token strings, even through a params file.
RAW_ARGUMENT_NATIVE_PARAMETERS = frozenset(
    {
        "extra_bowtie2_align_args",
        "extra_fastp_args",
        "extra_fqlint_args",
        "extra_kallisto_quant_args",
        "extra_salmon_quant_args",
        "extra_star_align_args",
        "extra_trimgalore_args",
    }
)

# Nextflow controls are not nf-core params, so the upstream schema cannot name
# them. Keeping an explicit deny set gives callers a stable failure class.
DANGEROUS_RUNTIME_PARAMETERS = frozenset(
    {
        "-c",
        "-resume",
        "cache",
        "cli_args",
        "command",
        "config",
        "config_file",
        "dag",
        "executor",
        "groovy_config",
        "launchDir",
        "launch_dir",
        "nextflow_args",
        "params_file",
        "plugin",
        "plugins",
        "profile",
        "raw_config",
        "report",
        "resume",
        "shell",
        "timeline",
        "tower",
        "trace",
        "wave",
        "with_tower",
        "workDir",
        "work_dir",
        "workdir",
    }
)

# These native values are emitted by normalization but cannot be authored on
# either user surface. They fix the Standard product route for later runtime
# and artifact contracts.
FIXED_PRODUCT_NATIVE_PARAMETERS = frozenset(
    {
        "igenomes_ignore",
        "skip_alignment",
        "skip_bbsplit",
        "skip_linting",
        "skip_pseudo_alignment",
        "skip_quantification_merge",
    }
)

PARAMETER_IMPACT_CATEGORIES = frozenset(
    {
        "additive_artifacts",
        "artifact_set",
        "content_only",
        "filename_or_format",
        "route_namespace",
        "sample_set",
        "subtractive_artifacts",
    }
)

# Artifact/QC extraction is deferred, but accepted parameter effects are not.
# This closed map is the evidence PR #152 will consume when fixing artifacts.
PARAMETER_IMPACT_BY_NATIVE_NAME = {
    "aligner": "route_namespace",
    "bam_csi_index": "filename_or_format",
    "deseq2_vst": "content_only",
    "fasta": "content_only",
    "featurecounts_feature_type": "content_only",
    "featurecounts_group_type": "content_only",
    "gencode": "content_only",
    "gffread_transcript_fasta": "content_only",
    "gtf": "content_only",
    "gtf_extra_attributes": "content_only",
    "gtf_group_features": "content_only",
    "igenomes_ignore": "content_only",
    "min_mapped_reads": "sample_set",
    "min_trimmed_reads": "sample_set",
    "rseqc_modules": "artifact_set",
    "salmon_index": "content_only",
    "salmon_quant_libtype": "content_only",
    "save_align_intermeds": "additive_artifacts",
    "save_merged_fastq": "additive_artifacts",
    "save_reference": "additive_artifacts",
    "save_trimmed": "additive_artifacts",
    "save_umi_intermeds": "additive_artifacts",
    "save_unaligned": "additive_artifacts",
    "skip_bigwig": "subtractive_artifacts",
    "skip_biotype_qc": "subtractive_artifacts",
    "skip_deseq2_qc": "subtractive_artifacts",
    "skip_dupradar": "subtractive_artifacts",
    "skip_fastqc": "subtractive_artifacts",
    "skip_gtf_filter": "content_only",
    "skip_gtf_transcript_filter": "content_only",
    "skip_alignment": "route_namespace",
    "skip_bbsplit": "subtractive_artifacts",
    "skip_linting": "content_only",
    "skip_markduplicates": "route_namespace",
    "skip_multiqc": "subtractive_artifacts",
    "skip_preseq": "subtractive_artifacts",
    "skip_pseudo_alignment": "route_namespace",
    "skip_qc": "subtractive_artifacts",
    "skip_qualimap": "subtractive_artifacts",
    "skip_quantification_merge": "route_namespace",
    "skip_rseqc": "subtractive_artifacts",
    "skip_stringtie": "subtractive_artifacts",
    "skip_trimming": "subtractive_artifacts",
    "skip_umi_extract": "route_namespace",
    "star_ignore_sjdbgtf": "content_only",
    "star_index": "content_only",
    "stranded_threshold": "route_namespace",
    "stringtie_ignore_gtf": "content_only",
    "trimmer": "route_namespace",
    "umi_dedup_tool": "route_namespace",
    "umi_discard_read": "route_namespace",
    "umitools_bc_pattern": "route_namespace",
    "umitools_bc_pattern2": "route_namespace",
    "umitools_dedup_primary_only": "content_only",
    "umitools_dedup_stats": "additive_artifacts",
    "umitools_extract_method": "route_namespace",
    "umitools_grouping_method": "content_only",
    "umitools_umi_separator": "route_namespace",
    "unstranded_threshold": "route_namespace",
    "with_umi": "route_namespace",
}


def load_upstream_parameter_schema() -> dict[str, Any]:
    """Load and verify the exact immutable upstream parameter schema."""
    return _load_json_contract(
        UPSTREAM_PARAMETER_SCHEMA_FILE,
        expected_size=UPSTREAM_PARAMETER_SCHEMA_SIZE,
        expected_sha256=UPSTREAM_PARAMETER_SCHEMA_SHA256,
    )


def load_upstream_samplesheet_schema() -> dict[str, Any]:
    """Load and verify the exact immutable upstream samplesheet schema."""
    return _load_json_contract(
        UPSTREAM_SAMPLESHEET_SCHEMA_FILE,
        expected_size=UPSTREAM_SAMPLESHEET_SCHEMA_SIZE,
        expected_sha256=UPSTREAM_SAMPLESHEET_SCHEMA_SHA256,
    )


def upstream_parameter_properties() -> dict[str, dict[str, Any]]:
    """Return fresh exact property schemas indexed by native parameter name."""
    document = load_upstream_parameter_schema()
    properties: dict[str, dict[str, Any]] = {}
    for group in document["$defs"].values():
        for name, schema in group.get("properties", {}).items():
            if name in properties:
                raise RuntimeError(
                    "Pinned upstream schema contains duplicate parameters."
                )
            properties[name] = deepcopy(schema)
    if len(properties) != 133:
        raise RuntimeError("Pinned upstream schema parameter count changed.")
    return properties


def projected_advanced_properties() -> dict[str, dict[str, Any]]:
    """Return UI-safe projections of the exact allowlisted property schemas."""
    upstream = upstream_parameter_properties()
    return {
        name: _strip_upstream_annotations(upstream[name])
        for name in sorted(ADVANCED_NATIVE_PARAMETERS)
    }


def _load_json_contract(
    filename: str,
    *,
    expected_size: int,
    expected_sha256: str,
) -> dict[str, Any]:
    data = resources.files(_CONTRACT_PACKAGE).joinpath(filename).read_bytes()
    if (
        len(data) != expected_size
        or hashlib.sha256(data).hexdigest() != expected_sha256
    ):
        raise RuntimeError("Pinned nf-core/rnaseq contract identity mismatch.")
    value = json.loads(data)
    if not isinstance(value, dict):
        raise RuntimeError("Pinned nf-core/rnaseq contract is not an object.")
    return value


def _strip_upstream_annotations(value: Any) -> Any:
    if isinstance(value, dict):
        omitted = {
            "default",
            "errorMessage",
            "exists",
            "fa_icon",
            "help_text",
            "hidden",
            "mimetype",
            "schema",
        }
        return {
            key: _strip_upstream_annotations(item)
            for key, item in value.items()
            if key not in omitted
        }
    if isinstance(value, list):
        return [_strip_upstream_annotations(item) for item in value]
    return deepcopy(value)


KNOWN_UPSTREAM_PARAMETERS = frozenset(upstream_parameter_properties())
UNSUPPORTED_NATIVE_PARAMETERS = KNOWN_UPSTREAM_PARAMETERS.difference(
    STANDARD_NATIVE_PARAMETERS,
    ADVANCED_NATIVE_PARAMETERS,
    PLATFORM_OWNED_NATIVE_PARAMETERS,
    RAW_ARGUMENT_NATIVE_PARAMETERS,
)

_POLICY_CLASSES = (
    STANDARD_NATIVE_PARAMETERS,
    ADVANCED_NATIVE_PARAMETERS,
    PLATFORM_OWNED_NATIVE_PARAMETERS,
    RAW_ARGUMENT_NATIVE_PARAMETERS,
    UNSUPPORTED_NATIVE_PARAMETERS,
)
if set().union(*_POLICY_CLASSES) != set(KNOWN_UPSTREAM_PARAMETERS):
    raise RuntimeError("nf-core/rnaseq parameter policy is incomplete.")
for index, group in enumerate(_POLICY_CLASSES):
    if any(group.intersection(other) for other in _POLICY_CLASSES[index + 1 :]):
        raise RuntimeError("nf-core/rnaseq parameter policy overlaps.")
if set(PARAMETER_IMPACT_BY_NATIVE_NAME) != set(
    STANDARD_NATIVE_PARAMETERS
    | ADVANCED_NATIVE_PARAMETERS
    | FIXED_PRODUCT_NATIVE_PARAMETERS
):
    raise RuntimeError("Accepted native parameter impact policy is incomplete.")
if not set(PARAMETER_IMPACT_BY_NATIVE_NAME.values()).issubset(
    PARAMETER_IMPACT_CATEGORIES
):
    raise RuntimeError("Native parameter impact category is invalid.")
