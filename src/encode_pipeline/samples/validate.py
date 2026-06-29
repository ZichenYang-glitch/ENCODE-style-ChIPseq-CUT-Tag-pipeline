"""Per-sample validation helpers."""

import csv
import os
import tempfile

from encode_pipeline.errors import ValidationError
from encode_pipeline.samples.load import load_and_validate_samples
from encode_pipeline.samples.models import SampleRecord

__all__ = [
    "ValidationError",
    "validate_sample_row",
    "SampleRecord",
]


def validate_sample_row(row: dict, *, use_control: bool = False) -> SampleRecord:
    """Validate a single sample row dict and return a typed SampleRecord.

    This is a convenience wrapper that builds a one-row sample sheet in
    memory and runs the sample loader. It is intended for unit tests and
    future CLI usage, not for hot-path Snakefile parsing.
    """
    required = [
        "sample", "fastq_1", "layout", "assay",
        "target", "peak_mode", "genome", "bowtie2_index",
    ]
    header = list(row.keys())
    for col in required:
        if col not in header:
            raise ValidationError(f"Sample row missing required key: {col!r}")

    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".tsv", delete=False, newline=""
    ) as fh:
        path = fh.name
        writer = csv.DictWriter(fh, fieldnames=header, delimiter="\t")
        writer.writeheader()
        writer.writerow(row)

    try:
        samples = load_and_validate_samples(
            path,
            use_control=use_control,
            strict_inputs=False,
        )
    finally:
        os.unlink(path)

    if len(samples) != 1:
        raise ValidationError("Expected exactly one sample from row validation")
    return SampleRecord.from_dict(samples[0])
