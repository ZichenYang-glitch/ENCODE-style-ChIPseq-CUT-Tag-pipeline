"""Strict input (file existence) validation for sample sheets."""

import os

from encode_pipeline.config import defaults
from encode_pipeline.errors import ValidationError

__all__ = [
    "validate_strict_inputs",
]


# Backward-compatible module aliases.
_BT2_STANDARD = defaults.BT2_STANDARD
_BT2_LARGE = defaults.BT2_LARGE


def _check_fastq_exists(
    path: str, sample_id: str, label: str, error_cls=ValidationError
) -> None:
    """Raise *error_cls* if *path* does not exist as a regular file."""
    if not os.path.isfile(path):
        raise error_cls(f"Sample {sample_id!r}: {label} file not found: {path}")


def _check_bowtie2_index(
    prefix: str, sample_id: str, error_cls=ValidationError
) -> None:
    """Raise *error_cls* if neither complete .bt2 nor .bt2l index set exists."""
    standard_set = [f.format(prefix=prefix) for f in _BT2_STANDARD]
    large_set = [f.format(prefix=prefix) for f in _BT2_LARGE]

    standard_ok = all(os.path.isfile(f) for f in standard_set)
    large_ok = all(os.path.isfile(f) for f in large_set)

    if standard_ok or large_ok:
        return

    missing = [f for f in standard_set if not os.path.isfile(f)]
    raise error_cls(
        f"Sample {sample_id!r}: Bowtie2 index not found at {prefix!r}. "
        f"Missing {len(missing)} of {len(standard_set)} expected .bt2 files. "
        f"Neither complete .bt2 nor .bt2l set exists. "
        f"First missing: {missing[0] if missing else 'n/a'}"
    )


def validate_strict_inputs(
    samples: list, strict_inputs: bool, error_cls=ValidationError
) -> None:
    """If *strict_inputs*, validate FASTQ and Bowtie2 index file existence."""
    if not strict_inputs:
        return
    for s in samples:
        sid = s.get("id", s.get("sample", "?"))
        fq1 = s.get("fq1", "")
        fq2 = s.get("fq2", "")
        bt2 = s.get("bt2_idx", "")
        layout = s.get("layout", "SE")

        if fq1:
            _check_fastq_exists(fq1, sid, "fastq_1", error_cls=error_cls)
        else:
            raise error_cls(f"Sample {sid!r}: fastq_1 is empty")
        if layout == "PE":
            if fq2:
                _check_fastq_exists(fq2, sid, "fastq_2", error_cls=error_cls)
            else:
                raise error_cls(f"Sample {sid!r}: PE layout requires fastq_2")
        if bt2:
            _check_bowtie2_index(bt2, sid, error_cls=error_cls)
        else:
            raise error_cls(f"Sample {sid!r}: bowtie2_index is empty")
