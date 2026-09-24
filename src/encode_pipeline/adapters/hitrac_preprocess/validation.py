"""Validate Hi-TrAC author input without reading submitted filesystem paths."""

from __future__ import annotations

from collections.abc import Mapping
import math
import re

from encode_pipeline.adapters.hitrac_preprocess.authoring import (
    DEFAULT_MAPQ,
    DEFAULT_THREADS,
    FASTQ_PATH_PATTERN,
    SAMPLE_FIELDS,
    SAMPLE_ID_PATTERN,
)
from encode_pipeline.platform.adapters import (
    MAX_SAMPLE_CELL_LENGTH,
    MAX_SAMPLE_ROWS,
    WorkflowInputs,
)
from encode_pipeline.platform.results import Issue, Result

_SAMPLE_ID = re.compile(SAMPLE_ID_PATTERN)
_FASTQ_PATH = re.compile(FASTQ_PATH_PATTERN)


def validate_hitrac_inputs(inputs: WorkflowInputs) -> Result[WorkflowInputs]:
    """Return copied inputs with defaults, retaining row order and authored values."""
    if not isinstance(inputs, WorkflowInputs):
        return _failure("HITRAC_INPUTS_INVALID", "inputs", "Submit workflow inputs.")
    if inputs.config:
        return _failure(
            "HITRAC_CONFIG_INVALID",
            "config",
            "Keep configuration empty and select a server-owned Reference Profile.",
        )
    if set(inputs.options).difference({"threads", "mapq"}):
        return _failure(
            "HITRAC_OPTIONS_INVALID", "options", "Use only threads and mapq options."
        )
    options = dict(inputs.options)
    options.setdefault("threads", DEFAULT_THREADS)
    options.setdefault("mapq", DEFAULT_MAPQ)
    for field, lower, upper in (("threads", 2, 8), ("mapq", 0, 255)):
        value = options[field]
        integer = type(value) is int or (
            type(value) is float and math.isfinite(value) and value.is_integer()
        )
        if not integer or not lower <= value <= upper:
            return _failure(
                "HITRAC_OPTION_INVALID",
                f"options.{field}",
                (
                    "Choose an integer thread count from 2 to 8."
                    if field == "threads"
                    else "Choose an integer MAPQ threshold from 0 to 255."
                ),
            )
        options[field] = int(value)
    samples = inputs.samples
    if not isinstance(samples, list) or not 1 <= len(samples) <= MAX_SAMPLE_ROWS:
        return _failure(
            "HITRAC_SAMPLES_INVALID",
            "samples",
            "Provide one or more inline sample rows with both gzip FASTQ mates.",
        )
    ids: set[str] = set()
    normalized = []
    for index, row in enumerate(samples):
        if not isinstance(row, Mapping) or set(row).difference(SAMPLE_FIELDS):
            return _failure(
                "HITRAC_SAMPLE_FIELDS_INVALID",
                f"samples[{index}]",
                "Use only sample_id, fastq_1 and fastq_2 columns.",
            )
        for field in SAMPLE_FIELDS:
            if field not in row:
                return _failure(
                    "HITRAC_SAMPLE_FIELD_REQUIRED",
                    f"samples[{index}].{field}",
                    "Provide a sample ID and both gzip FASTQ paths for every sample.",
                )
        sample_id = row["sample_id"]
        if not isinstance(sample_id, str) or _SAMPLE_ID.fullmatch(sample_id) is None:
            return _failure(
                "HITRAC_SAMPLE_ID_INVALID",
                f"samples[{index}].sample_id",
                "Use 1 to 128 characters: start with a letter or digit, then use "
                "letters, digits, spaces, dots, underscores or hyphens.",
            )
        if sample_id in ids:
            return _failure(
                "HITRAC_SAMPLE_ID_DUPLICATE",
                f"samples[{index}].sample_id",
                "Assign a distinct sample ID to each row; premerge lanes per sample.",
            )
        ids.add(sample_id)
        for field in ("fastq_1", "fastq_2"):
            value = row[field]
            if (
                not isinstance(value, str)
                or len(value) > MAX_SAMPLE_CELL_LENGTH
                or _FASTQ_PATH.fullmatch(value) is None
            ):
                return _failure(
                    "HITRAC_FASTQ_PATH_INVALID",
                    f"samples[{index}].{field}",
                    "Provide an absolute POSIX path to a gzip FASTQ without empty, "
                    "dot or parent segments, backslashes or control characters.",
                )
        normalized.append({field: row[field] for field in SAMPLE_FIELDS})
    return Result.success(
        WorkflowInputs(config={}, samples=normalized, options=options)
    )


def _failure(code: str, path: str, hint: str) -> Result[WorkflowInputs]:
    return Result.failure(
        [
            Issue(
                code=code,
                message="Hi-TrAC preprocessing inputs do not match the authoring contract.",
                path=path,
                source="hitrac_preprocess",
                hint=hint,
            )
        ]
    )
