"""Versioned authoring contract for the bulk RNA-seq adapter."""

from __future__ import annotations

from copy import deepcopy

from encode_pipeline.adapters.bulk_rnaseq.upstream import (
    NFCORE_RNASEQ_COMMIT,
    NFCORE_RNASEQ_RELEASE,
    projected_advanced_properties,
)
from encode_pipeline.platform.adapters import (
    JSON_SCHEMA_DIALECT,
    MAX_SAMPLE_CELL_LENGTH,
    MAX_SAMPLE_COLUMNS,
    MAX_SAMPLE_ROWS,
    WorkflowAuthoringModes,
    WorkflowInputLimits,
    WorkflowInputModes,
    WorkflowSchema,
    WorkflowSchemaCoverage,
)


SCHEMA_VERSION = "1.0.0"
_SCHEMA_ID_ROOT = "https://helixweave.org/schemas/bulk-rnaseq"
_IDENTIFIER_PATTERN = r"^[A-Za-z_][A-Za-z0-9_.-]*(?:,[A-Za-z_][A-Za-z0-9_.-]*)*$"
_PATH_SEGMENT = r"(?!\.{1,2}(?:/|$))[A-Za-z0-9._+@%=-]+"
_ABSOLUTE_PATH_BODY = rf"/{_PATH_SEGMENT}(?:/{_PATH_SEGMENT})*"
_ABSOLUTE_PATH_PATTERN = rf"^{_ABSOLUTE_PATH_BODY}$"
_FASTA_PATH_PATTERN = rf"^{_ABSOLUTE_PATH_BODY}\.(?:fa|fasta|fna)(?:\.gz)?$"
_GTF_PATH_PATTERN = rf"^{_ABSOLUTE_PATH_BODY}\.gtf(?:\.gz)?$"
_FASTQ_PATH_PATTERN = rf"^{_ABSOLUTE_PATH_BODY}\.(?:fastq|fq)\.gz$"
_OPTIONAL_FASTQ_PATH_PATTERN = rf"^(?:|{_ABSOLUTE_PATH_BODY}\.(?:fastq|fq)\.gz)$"
_SHA256_PATTERN = "^[0-9a-f]{64}$"

ANALYSIS_DEFAULTS = {"alignment": "star", "quantification": "salmon"}
TRIMMING_DEFAULTS = {"enabled": True, "tool": "trimgalore"}
RIBOSOMAL_RNA_REMOVAL_DEFAULTS = {"enabled": False}
QC_DEFAULTS = {
    "enabled": True,
    "fastqc": True,
    "multiqc": True,
    "rseqc": True,
    "qualimap": True,
    "dupradar": True,
    "biotype": True,
    "deseq2_pca": True,
    "preseq": False,
    "mark_duplicates": True,
}
OUTPUT_DEFAULTS = {
    "bigwig": True,
    "stringtie": False,
    "trimmed_reads": False,
    "unaligned_reads": False,
    "reference_files": False,
    "alignment_intermediates": False,
    "merged_fastq": False,
    "umi_intermediates": False,
}


def build_bulk_rnaseq_authoring_schema() -> WorkflowSchema:
    """Return a fresh, stable adapter-owned authoring contract."""
    return WorkflowSchema(
        schema_version=SCHEMA_VERSION,
        schema_dialect=JSON_SCHEMA_DIALECT,
        coverage=WorkflowSchemaCoverage(
            config="complete",
            samples="complete",
            options="complete",
        ),
        authoring_modes=WorkflowAuthoringModes(
            config=("schema_form", "yaml"),
            samples=("tsv_upload", "inline_table"),
            options=("schema_form",),
        ),
        input_modes=WorkflowInputModes(
            config=("object",),
            samples=("inline_rows",),
            options=("object",),
        ),
        limits=WorkflowInputLimits(),
        config_schema=_config_schema(),
        sample_schema=_sample_schema(),
        option_schema=_option_schema(),
    )


def _document_id(surface: str) -> str:
    return f"{_SCHEMA_ID_ROOT}/{surface}/{SCHEMA_VERSION}"


def _path_schema(
    *, title: str, pattern: str = _ABSOLUTE_PATH_PATTERN
) -> dict[str, object]:
    return {
        "type": "string",
        "title": title,
        "minLength": 2,
        "maxLength": MAX_SAMPLE_CELL_LENGTH,
        "pattern": pattern,
    }


def _index_schema(*, title: str) -> dict[str, object]:
    return {
        "type": "object",
        "title": title,
        "properties": {
            "path": _path_schema(title=f"{title} path"),
            "identity_sha256": {
                "type": "string",
                "pattern": _SHA256_PATTERN,
                "description": "SHA-256 identity of the immutable index manifest.",
            },
        },
        "required": ["path", "identity_sha256"],
        "additionalProperties": False,
    }


def _ribosomal_rna_removal_schema() -> dict[str, object]:
    return {
        "type": "object",
        "default": RIBOSOMAL_RNA_REMOVAL_DEFAULTS,
        "description": (
            "Optional reference-based ribosomal RNA removal. Database and index "
            "resources carry immutable identities; submitted paths are not probed "
            "during contract validation."
        ),
        "properties": {
            "enabled": {"type": "boolean", "default": False},
            "tool": {
                "type": "string",
                "enum": ["sortmerna", "bowtie2"],
            },
            "save_filtered_reads": {"type": "boolean", "default": False},
            "database_manifest": {
                "type": "object",
                "title": "rRNA database manifest",
                "description": (
                    "Adapter-owned reference to the manifest consumed by the "
                    "selected removal tool."
                ),
                "properties": {
                    "path": _path_schema(title="rRNA database manifest path"),
                    "identity_sha256": {
                        "type": "string",
                        "pattern": _SHA256_PATTERN,
                        "description": (
                            "SHA-256 identity of the exact database manifest bytes."
                        ),
                    },
                },
                "required": ["path", "identity_sha256"],
                "additionalProperties": False,
            },
            "sortmerna_index": _index_schema(title="SortMeRNA index"),
        },
        "required": ["enabled"],
        "allOf": [
            {
                "if": {
                    "required": ["enabled"],
                    "properties": {"enabled": {"const": True}},
                },
                "then": {"required": ["tool", "database_manifest"]},
            }
        ],
        "additionalProperties": False,
    }


def _config_schema() -> dict[str, object]:
    standard_schema: dict[str, object] = {
        "type": "object",
        "title": "Standard bulk RNA-seq parameters",
        "description": (
            "Stable HelixWeave scientific semantics. Native nf-core names belong "
            "only in the sibling Advanced object when explicitly allowlisted."
        ),
        "properties": {
            "analysis": {
                "type": "object",
                "default": ANALYSIS_DEFAULTS,
                "properties": {
                    "alignment": {
                        "type": "string",
                        "const": "star",
                        "default": "star",
                    },
                    "quantification": {
                        "type": "string",
                        "const": "salmon",
                        "default": "salmon",
                    },
                },
                "required": ["alignment", "quantification"],
                "additionalProperties": False,
            },
            "reference": {
                "type": "object",
                "description": (
                    "Explicit reference files and immutable identities; iGenomes "
                    "and implicit remote references are not supported."
                ),
                "properties": {
                    "reference_id": {
                        "type": "string",
                        "minLength": 1,
                        "maxLength": 128,
                        "pattern": "^[A-Za-z0-9][A-Za-z0-9_.-]*$",
                    },
                    "fasta": _path_schema(
                        title="Genome FASTA",
                        pattern=_FASTA_PATH_PATTERN,
                    ),
                    "fasta_sha256": {
                        "type": "string",
                        "pattern": _SHA256_PATTERN,
                    },
                    "gtf": _path_schema(
                        title="Gene annotation GTF",
                        pattern=_GTF_PATH_PATTERN,
                    ),
                    "gtf_sha256": {
                        "type": "string",
                        "pattern": _SHA256_PATTERN,
                    },
                    "annotation_style": {
                        "type": "string",
                        "enum": ["ensembl", "gencode"],
                        "default": "ensembl",
                    },
                    "star_index": _index_schema(title="STAR index"),
                    "salmon_index": {
                        **_index_schema(title="Salmon index"),
                        "description": (
                            "Optional prebuilt index for samples with auto "
                            "strandedness. STAR+Salmon quantification itself uses "
                            "the STAR transcriptome alignment route."
                        ),
                    },
                },
                "required": [
                    "reference_id",
                    "fasta",
                    "fasta_sha256",
                    "gtf",
                    "gtf_sha256",
                ],
                "additionalProperties": False,
            },
            "trimming": {
                "type": "object",
                "default": TRIMMING_DEFAULTS,
                "properties": {
                    "enabled": {"type": "boolean", "default": True},
                    "tool": {
                        "type": "string",
                        "enum": ["trimgalore", "fastp"],
                        "default": "trimgalore",
                    },
                },
                "required": ["enabled", "tool"],
                "additionalProperties": False,
            },
            "ribosomal_rna_removal": _ribosomal_rna_removal_schema(),
            "umi": _umi_schema(),
            "qc": _qc_schema(),
            "outputs": _outputs_schema(),
        },
        "required": ["reference"],
        "additionalProperties": False,
    }

    return {
        "$schema": JSON_SCHEMA_DIALECT,
        "$id": _document_id("config"),
        "title": "HelixWeave bulk RNA-seq configuration",
        "description": (
            f"Contract for nf-core/rnaseq {NFCORE_RNASEQ_RELEASE} at immutable "
            f"commit {NFCORE_RNASEQ_COMMIT}."
        ),
        "type": "object",
        "properties": {
            "standard": standard_schema,
            "advanced": _advanced_schema(),
        },
        "required": ["standard"],
        "additionalProperties": False,
    }


def _umi_schema() -> dict[str, object]:
    return {
        "type": "object",
        "default": {"enabled": False},
        "properties": {
            "enabled": {"type": "boolean", "default": False},
            "mode": {
                "type": "string",
                "enum": ["read_name", "read_sequence"],
            },
            "deduplication_tool": {
                "type": "string",
                "const": "umitools",
                "default": "umitools",
            },
            "extraction_method": {
                "type": "string",
                "enum": ["string", "regex"],
            },
            "barcode_pattern": {
                "type": "string",
                "minLength": 1,
                "maxLength": 256,
                "pattern": "^[^\\s'\"`$;|&]+$",
            },
            "barcode_pattern_2": {
                "type": "string",
                "minLength": 1,
                "maxLength": 256,
                "pattern": "^[^\\s'\"`$;|&]+$",
            },
            "read_name_separator": {
                "type": "string",
                "enum": [":", "_", "-", "+", "."],
            },
            "discard_read": {"type": "integer", "enum": [1, 2]},
            "grouping_method": {
                "type": "string",
                "enum": [
                    "unique",
                    "percentile",
                    "cluster",
                    "adjacency",
                    "directional",
                ],
                "default": "directional",
            },
            "emit_dedup_stats": {"type": "boolean", "default": False},
            "primary_alignments_only": {"type": "boolean", "default": False},
        },
        "required": ["enabled"],
        "allOf": [
            {
                "if": {"properties": {"enabled": {"const": True}}},
                "then": {
                    "required": ["mode", "deduplication_tool"],
                },
            },
            {
                "if": {
                    "required": ["mode"],
                    "properties": {"mode": {"const": "read_name"}},
                },
                "then": {"required": ["read_name_separator"]},
            },
            {
                "if": {
                    "required": ["mode"],
                    "properties": {"mode": {"const": "read_sequence"}},
                },
                "then": {"required": ["extraction_method", "barcode_pattern"]},
            },
        ],
        "additionalProperties": False,
    }


def _qc_schema() -> dict[str, object]:
    properties = {
        name: {"type": "boolean", "default": default}
        for name, default in QC_DEFAULTS.items()
    }
    properties["mark_duplicates"]["description"] = (
        "Run Picard MarkDuplicates for non-UMI data. UMI-enabled data always "
        "uses the upstream UMI-aware deduplication branch instead."
    )
    return {
        "type": "object",
        "default": QC_DEFAULTS,
        "description": (
            "DESeq2 here is limited to exploratory QC transformation/PCA; this "
            "workflow does not perform differential-expression significance tests."
        ),
        "properties": properties,
        "additionalProperties": False,
    }


def _outputs_schema() -> dict[str, object]:
    return {
        "type": "object",
        "default": OUTPUT_DEFAULTS,
        "properties": {
            name: {"type": "boolean", "default": default}
            for name, default in OUTPUT_DEFAULTS.items()
        },
        "additionalProperties": False,
    }


def _advanced_schema() -> dict[str, object]:
    properties = projected_advanced_properties()
    properties["min_trimmed_reads"]["minimum"] = 1
    properties["min_mapped_reads"].update({"minimum": 0, "maximum": 100})
    for name in (
        "gtf_extra_attributes",
        "gtf_group_features",
        "featurecounts_group_type",
        "featurecounts_feature_type",
    ):
        properties[name].update({"pattern": _IDENTIFIER_PATTERN, "maxLength": 256})
    properties["rseqc_modules"].update(
        {
            "pattern": (
                "^(?:bam_stat|inner_distance|infer_experiment|junction_annotation|"
                "junction_saturation|read_distribution|read_duplication|tin)"
                "(?:,(?:bam_stat|inner_distance|infer_experiment|"
                "junction_annotation|junction_saturation|read_distribution|"
                "read_duplication|tin))*$"
            ),
            "maxLength": 256,
        }
    )
    return {
        "type": "object",
        "title": "Advanced nf-core/rnaseq parameters",
        "description": (
            "Closed allowlist validated against the exact pinned upstream schema. "
            "Runtime, path ownership, raw CLI, and Standard aliases are forbidden."
        ),
        "properties": properties,
        "additionalProperties": False,
    }


def _sample_schema() -> dict[str, object]:
    identifier = {
        "type": "string",
        "minLength": 1,
        "maxLength": 128,
        "pattern": "^[A-Za-z0-9][A-Za-z0-9_.-]*$",
    }
    fastq = {
        "type": "string",
        "maxLength": MAX_SAMPLE_CELL_LENGTH,
        "pattern": _FASTQ_PATH_PATTERN,
    }
    fastq_2 = {
        "type": "string",
        "maxLength": MAX_SAMPLE_CELL_LENGTH,
        "pattern": _OPTIONAL_FASTQ_PATH_PATTERN,
    }
    return {
        "$schema": JSON_SCHEMA_DIALECT,
        "$id": _document_id("samples"),
        "title": "HelixWeave bulk RNA-seq sample rows",
        "description": (
            "Repeated sample IDs represent technical libraries or lanes that "
            "nf-core concatenates before analysis. Biological replicates use "
            "distinct sample IDs."
        ),
        "type": "array",
        "minItems": 1,
        "maxItems": MAX_SAMPLE_ROWS,
        "items": {
            "type": "object",
            "minProperties": 7,
            "maxProperties": MAX_SAMPLE_COLUMNS,
            "properties": {
                "sample": deepcopy(identifier),
                "library": deepcopy(identifier),
                "lane": deepcopy(identifier),
                "layout": {"type": "string", "enum": ["SE", "PE"]},
                "fastq_1": fastq,
                "fastq_2": fastq_2,
                "strandedness": {
                    "type": "string",
                    "enum": ["auto", "forward", "reverse", "unstranded"],
                },
                "platform": {"type": "string", "const": "ILLUMINA"},
            },
            "required": [
                "sample",
                "library",
                "lane",
                "layout",
                "fastq_1",
                "strandedness",
                "platform",
            ],
            "allOf": [
                {
                    "if": {"properties": {"layout": {"const": "PE"}}},
                    "then": {
                        "required": ["fastq_2"],
                        "properties": {"fastq_2": {"minLength": 1}},
                    },
                    "else": {
                        "properties": {"fastq_2": {"maxLength": 0}},
                    },
                }
            ],
            "additionalProperties": False,
        },
    }


def _option_schema() -> dict[str, object]:
    return {
        "$schema": JSON_SCHEMA_DIALECT,
        "$id": _document_id("options"),
        "title": "HelixWeave bulk RNA-seq adapter options",
        "description": "No caller-owned platform options are defined in schema 1.0.0.",
        "type": "object",
        "properties": {},
        "additionalProperties": False,
    }
