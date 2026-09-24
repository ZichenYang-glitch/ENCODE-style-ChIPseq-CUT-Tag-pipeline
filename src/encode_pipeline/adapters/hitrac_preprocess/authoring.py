"""Adapter-owned authoring for the qualified Hi-TrAC preprocessing route."""

from __future__ import annotations

from encode_pipeline.platform.adapters import (
    JSON_SCHEMA_DIALECT,
    MAX_SAMPLE_CELL_LENGTH,
    MAX_SAMPLE_ROWS,
    WorkflowAuthoringModes,
    WorkflowInputModes,
    WorkflowSchema,
    WorkflowSchemaCoverage,
)

SCHEMA_VERSION = "1.0.0"
DEFAULT_THREADS = 2
DEFAULT_MAPQ = 10
SAMPLE_FIELDS = ("sample_id", "fastq_1", "fastq_2")
SAMPLE_ID_PATTERN = r"^[A-Za-z0-9][A-Za-z0-9_. -]{0,127}(?![\s\S])"
_PATH_SEGMENT = r"(?!\.{1,2}(?:/|$))[^/\x00-\x1f\x7f\\]+"
FASTQ_PATH_PATTERN = rf"^/{_PATH_SEGMENT}(?:/{_PATH_SEGMENT})*(?![\s\S])"
_SCHEMA_ROOT = "https://helixweave.org/schemas/hitrac-preprocess"


def build_hitrac_authoring_schema() -> WorkflowSchema:
    """Describe author input only; runtime and reference bindings are server-owned."""
    path_schema = {
        "type": "string",
        "minLength": 2,
        "maxLength": MAX_SAMPLE_CELL_LENGTH,
        "pattern": FASTQ_PATH_PATTERN,
        "description": (
            "Absolute server path to one gzip FASTQ. Content and paired reads "
            "are verified at execution admission; validation does not read this path."
        ),
    }
    return WorkflowSchema(
        schema_version=SCHEMA_VERSION,
        coverage=WorkflowSchemaCoverage(
            config="complete", samples="complete", options="complete"
        ),
        authoring_modes=WorkflowAuthoringModes(
            config=("schema_form", "yaml"),
            samples=("inline_table", "tsv_upload"),
            options=("schema_form",),
        ),
        input_modes=WorkflowInputModes(
            config=("object",), samples=("inline_rows",), options=("object",)
        ),
        config_schema={
            **_document("config"),
            "type": "object",
            "properties": {},
            "additionalProperties": False,
            "description": (
                "Select a server-owned Reference Profile revision separately. "
                "No runtime paths, reference digests or scientific overrides are accepted."
            ),
        },
        sample_schema={
            **_document("samples"),
            "type": "array",
            "minItems": 1,
            "maxItems": MAX_SAMPLE_ROWS,
            "uniqueItems": True,
            "description": (
                "Ordered samples with unique sample_id values. Each sample has one "
                "paired gzip FASTQ; merge lanes before submission. Sample IDs are "
                "preserved while private execution tokens follow this row order."
            ),
            "items": {
                "type": "object",
                "properties": {
                    "sample_id": {
                        "type": "string",
                        "title": "Sample ID",
                        "pattern": SAMPLE_ID_PATTERN,
                        "minLength": 1,
                        "maxLength": 128,
                    },
                    "fastq_1": {**path_schema, "title": "Read 1 gzip FASTQ"},
                    "fastq_2": {**path_schema, "title": "Read 2 gzip FASTQ"},
                },
                "required": list(SAMPLE_FIELDS),
                "additionalProperties": False,
            },
        },
        option_schema={
            **_document("options"),
            "type": "object",
            "properties": {
                "threads": {
                    "type": "integer",
                    "minimum": 2,
                    "maximum": 8,
                    "default": DEFAULT_THREADS,
                    "title": "Tool threads",
                    "description": (
                        "Threads passed to the qualified original script; other tools "
                        "also create threads. This is not a total CPU resource limit."
                    ),
                },
                "mapq": {
                    "type": "integer",
                    "minimum": 0,
                    "maximum": 255,
                    "default": DEFAULT_MAPQ,
                    "title": "MAPQ threshold",
                },
            },
            "additionalProperties": False,
        },
    )


def _document(surface: str) -> dict[str, str]:
    return {
        "$schema": JSON_SCHEMA_DIALECT,
        "$id": f"{_SCHEMA_ROOT}/{surface}/{SCHEMA_VERSION}",
    }
