"""Versioned authoring contract for the ENCODE-style workflow adapter."""

from __future__ import annotations

from copy import deepcopy

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
_SCHEMA_ID_ROOT = (
    "https://encode-pipeline.org/schemas/encode-style-chipseq-cuttag-atac-mnase"
)


def build_encode_authoring_schema() -> WorkflowSchema:
    """Return a fresh versioned authoring contract for the ENCODE adapter."""
    return WorkflowSchema(
        schema_version=SCHEMA_VERSION,
        schema_dialect=JSON_SCHEMA_DIALECT,
        coverage=WorkflowSchemaCoverage(
            config="partial",
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
            samples=("inline_rows", "server_path"),
            options=("object",),
        ),
        limits=WorkflowInputLimits(),
        config_schema=_config_schema(),
        sample_schema=_sample_schema(),
        option_schema=_option_schema(),
    )


def _document_id(surface: str) -> str:
    return f"{_SCHEMA_ID_ROOT}/{surface}/{SCHEMA_VERSION}"


def _boolean_or_string_schema(*, default: bool) -> dict[str, object]:
    return {
        "oneOf": [
            {"type": "boolean"},
            {"type": "string", "enum": ["true", "false"]},
        ],
        "default": default,
    }


def _config_schema() -> dict[str, object]:
    qc_properties = {
        key: _boolean_or_string_schema(default=default)
        for key, default in (
            ("blacklist_filter", True),
            ("frip", True),
            ("library_complexity", True),
            ("nrf_pbc", True),
            ("signal_tracks", True),
            ("summary", True),
            ("cuttag_fragment_size", True),
            ("cross_correlation", False),
            ("preseq_complexity", False),
            ("picard_metrics", False),
            ("tss_enrichment", False),
        )
    }
    schema: dict[str, object] = {
        "$schema": JSON_SCHEMA_DIALECT,
        "$id": _document_id("config"),
        "title": "ENCODE-style workflow configuration",
        "description": (
            "A partial form contract covering the deterministic local profile. "
            "Advanced YAML keys remain backend-authoritative. Samples are "
            "submitted through the sibling samples input."
        ),
        "type": "object",
        "properties": {
            "outdir": {
                "type": "string",
                "default": "results",
                "maxLength": MAX_SAMPLE_CELL_LENGTH,
            },
            "threads": {"type": "integer", "minimum": 1, "default": 8},
            "mapq": {"type": "integer", "minimum": 0, "default": 30},
            "binsize": {"type": "integer", "minimum": 1, "default": 10},
            "remove_dup": {
                "type": "string",
                "enum": ["auto", "yes", "no"],
                "default": "auto",
            },
            "trim": _boolean_or_string_schema(default=True),
            "extend_reads": {
                "oneOf": [
                    {"type": "string", "enum": ["auto", "yes", "no"]},
                    {"type": "string", "pattern": "^[1-9][0-9]*$"},
                ],
                "default": "auto",
            },
            "use_control": {
                "oneOf": [
                    {"type": "boolean"},
                    {
                        "type": "string",
                        "enum": ["true", "false", "yes", "no", "1", "0"],
                    },
                ],
                "default": False,
            },
            "multiqc": _boolean_or_string_schema(default=True),
            "stage4b": _boolean_or_string_schema(default=True),
            "stage5": _boolean_or_string_schema(default=False),
            "genome_resources": {
                "type": "object",
                "default": {},
                "additionalProperties": {
                    "type": "object",
                    "properties": {
                        "effective_genome_size": {
                            "oneOf": [
                                {"type": "string", "enum": ["hs", "mm"]},
                                {"type": "integer", "minimum": 1},
                            ]
                        },
                        "chrom_sizes": {
                            "type": "string",
                            "maxLength": MAX_SAMPLE_CELL_LENGTH,
                        },
                        "blacklist": {
                            "type": "string",
                            "maxLength": MAX_SAMPLE_CELL_LENGTH,
                        },
                        "gtf": {
                            "type": "string",
                            "maxLength": MAX_SAMPLE_CELL_LENGTH,
                        },
                        "reference_fasta": {
                            "type": "string",
                            "maxLength": MAX_SAMPLE_CELL_LENGTH,
                        },
                    },
                    "additionalProperties": False,
                },
            },
            "qc": {
                "type": "object",
                "default": {},
                "properties": qc_properties,
                "additionalProperties": True,
            },
        },
        "additionalProperties": True,
    }
    return deepcopy(schema)


def _sample_schema() -> dict[str, object]:
    text = {"type": "string", "maxLength": MAX_SAMPLE_CELL_LENGTH}
    required_text = {
        "type": "string",
        "minLength": 1,
        "maxLength": MAX_SAMPLE_CELL_LENGTH,
    }
    optional_positive_integer = {
        "type": "string",
        "pattern": "^(?:|[1-9][0-9]*)$",
        "maxLength": MAX_SAMPLE_CELL_LENGTH,
    }
    properties: dict[str, object] = {
        "sample": {
            **required_text,
            "pattern": "^[A-Za-z0-9_.-]+$",
            "description": "Unique sample identifier.",
        },
        "fastq_1": {**required_text, "description": "Absolute R1 FASTQ path."},
        "fastq_2": {**text, "description": "Absolute R2 FASTQ path for PE data."},
        "layout": {"type": "string", "enum": ["PE", "SE"]},
        "assay": {
            "type": "string",
            "enum": ["chipseq", "cuttag", "atac", "mnase"],
        },
        "target": deepcopy(required_text),
        "peak_mode": {
            "type": "string",
            "enum": ["narrow", "broad", "nucleosome"],
        },
        "genome": deepcopy(required_text),
        "bowtie2_index": {
            **required_text,
            "description": "Absolute Bowtie2 index prefix.",
        },
        "control_bam": deepcopy(text),
        "role": {
            "type": "string",
            "enum": ["", "treatment", "control"],
            "default": "treatment",
        },
        "control_sample": deepcopy(text),
        "experiment": deepcopy(text),
        "condition": deepcopy(text),
        "replicate": deepcopy(optional_positive_integer),
        "biological_replicate": deepcopy(optional_positive_integer),
        "technical_replicate": deepcopy(optional_positive_integer),
    }
    return {
        "$schema": JSON_SCHEMA_DIALECT,
        "$id": _document_id("samples"),
        "title": "ENCODE-style sample rows",
        "type": "array",
        "minItems": 1,
        "maxItems": MAX_SAMPLE_ROWS,
        "items": {
            "type": "object",
            "minProperties": 1,
            "maxProperties": MAX_SAMPLE_COLUMNS,
            "required": [
                "sample",
                "fastq_1",
                "layout",
                "assay",
                "target",
                "peak_mode",
                "genome",
                "bowtie2_index",
            ],
            "properties": properties,
            "additionalProperties": False,
        },
    }


def _option_schema() -> dict[str, object]:
    return {
        "$schema": JSON_SCHEMA_DIALECT,
        "$id": _document_id("options"),
        "title": "ENCODE-style adapter options",
        "type": "object",
        "properties": {
            "strict_inputs": {
                "type": "boolean",
                "default": False,
                "description": "Validate FASTQ and Bowtie2 index file existence.",
            }
        },
        "additionalProperties": False,
    }
