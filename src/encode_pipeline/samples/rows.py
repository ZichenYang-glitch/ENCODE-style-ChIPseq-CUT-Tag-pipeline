"""Per-row sample sheet normalization and validation."""

import os

from encode_pipeline.config import defaults
from encode_pipeline.errors import ValidationError

__all__ = [
    "validate_and_build_sample",
]


# Backward-compatible module aliases.
_SAMPLE_ID_RE = defaults.SAMPLE_ID_RE
_SANITIZE_RE = defaults.SANITIZE_RE


def _sanitize_identifier(value: str) -> str:
    """Replace any character not in the safe identifier set with ``_``."""
    return _SANITIZE_RE.sub("_", value)


def _parse_positive_int(value, *, default, name, sample, row, error_cls=ValueError):
    """Parse an optional positive integer with a default.

    Blanks and missing values use *default*. Strictly rejects
    zero, negative, and non-integer values.
    """
    raw = (value or "").strip()
    if not raw:
        return default
    try:
        parsed = int(raw)
    except (ValueError, TypeError):
        raise error_cls(
            f"Row {row}: sample {sample!r}: {name} must be a positive "
            f"integer, got {raw!r}"
        )
    if parsed <= 0:
        raise error_cls(
            f"Row {row}: sample {sample!r}: {name} must be a positive "
            f"integer, got {parsed}"
        )
    return parsed


def validate_and_build_sample(
    row: dict,
    row_index: int,
    *,
    use_control: bool = False,
    error_cls=ValidationError,
) -> dict:
    """Normalize and validate one TSV row, returning the canonical sample dict.

    The returned dict keys and insertion order match the legacy
    ``load_and_validate_samples`` implementation.
    """
    sid = (row.get("sample") or "").strip()
    fq1 = (row.get("fastq_1") or "").strip()
    fq2 = (row.get("fastq_2") or "").strip()
    lo = (row.get("layout") or "").strip().upper()
    assay = (row.get("assay") or "").strip().lower()
    tgt = (row.get("target") or "").strip()
    pmode = (row.get("peak_mode") or "").strip().lower()
    gnm = (row.get("genome") or "").strip()
    bt2 = (row.get("bowtie2_index") or "").strip()

    ctrl_bam = (row.get("control_bam") or "").strip()
    role = (row.get("role") or "treatment").strip().lower()
    ctrl_sample = (row.get("control_sample") or "").strip()

    experiment = _sanitize_identifier((row.get("experiment") or "").strip() or sid)
    condition = _sanitize_identifier((row.get("condition") or "").strip() or tgt)
    replicate = _parse_positive_int(
        row.get("replicate"),
        default=1,
        name="replicate",
        sample=sid,
        row=row_index,
        error_cls=error_cls,
    )
    bio_rep = _parse_positive_int(
        row.get("biological_replicate"),
        default=replicate,
        name="biological_replicate",
        sample=sid,
        row=row_index,
        error_cls=error_cls,
    )
    tech_rep = _parse_positive_int(
        row.get("technical_replicate"),
        default=1,
        name="technical_replicate",
        sample=sid,
        row=row_index,
        error_cls=error_cls,
    )

    if not sid:
        raise error_cls(f"Row {row_index} in sample sheet has empty 'sample'")
    if not _SAMPLE_ID_RE.match(sid):
        raise error_cls(
            f"Row {row_index}: sample ID {sid!r} contains invalid "
            f"characters. Allowed: A-Z, a-z, 0-9, underscore, "
            f"period, hyphen."
        )
    if not fq1:
        raise error_cls(f"Sample {sid!r} has empty 'fastq_1'")
    if lo not in defaults.LAYOUTS:
        raise error_cls(f"Sample {sid!r}: layout must be PE or SE, got {lo!r}")
    if lo == "PE" and not fq2:
        raise error_cls(f"Sample {sid!r}: PE layout requires 'fastq_2'")
    if assay == "mnase" and lo != "PE":
        raise error_cls(
            f"Sample {sid!r}: assay=mnase requires paired-end layout (PE), got {lo!r}"
        )
    if assay not in defaults.ASSAYS:
        raise error_cls(
            f"Sample {sid!r}: assay must be chipseq, cuttag, "
            f"atac, or mnase, got {assay!r}"
        )
    if not tgt:
        raise error_cls(f"Sample {sid!r} has empty 'target'")
    if pmode not in defaults.PEAK_MODES:
        raise error_cls(
            f"Sample {sid!r}: peak_mode must be narrow, broad, "
            f"or nucleosome, got {pmode!r}"
        )
    if assay == "atac" and pmode != "narrow":
        raise error_cls(
            f"Sample {sid!r}: assay=atac currently supports "
            f"peak_mode=narrow only, got {pmode!r}"
        )
    if assay == "mnase" and pmode != "nucleosome":
        raise error_cls(
            f"Sample {sid!r}: assay=mnase requires peak_mode=nucleosome, got {pmode!r}"
        )
    if assay != "mnase" and pmode == "nucleosome":
        raise error_cls(
            f"Sample {sid!r}: peak_mode=nucleosome is only "
            f"allowed for assay=mnase, got assay={assay!r}"
        )
    if not gnm:
        raise error_cls(f"Sample {sid!r} has empty 'genome'")
    if not bt2:
        raise error_cls(f"Sample {sid!r} has empty 'bowtie2_index'")
    if not _SAMPLE_ID_RE.match(experiment):
        raise error_cls(
            f"Row {row_index}: sample {sid!r}: experiment {experiment!r} "
            f"contains invalid characters. "
            f"Allowed: A-Z, a-z, 0-9, underscore, period, hyphen."
        )
    if not _SAMPLE_ID_RE.match(condition):
        raise error_cls(
            f"Row {row_index}: sample {sid!r}: condition {condition!r} "
            f"contains invalid characters. "
            f"Allowed: A-Z, a-z, 0-9, underscore, period, hyphen."
        )
    if not experiment:
        raise error_cls(f"Sample {sid!r} has empty 'experiment' after defaulting")
    if not condition:
        raise error_cls(f"Sample {sid!r} has empty 'condition' after defaulting")
    if role not in defaults.ROLES:
        raise error_cls(
            f"Sample {sid!r}: role must be treatment or control, got {role!r}"
        )

    if use_control:
        if ctrl_sample == sid:
            raise error_cls(f"Sample {sid!r}: control_sample cannot reference itself")
        if ctrl_sample and ctrl_bam:
            raise error_cls(
                f"Sample {sid!r}: cannot set both control_sample and control_bam"
            )
        if ctrl_bam and not os.path.isfile(ctrl_bam):
            raise error_cls(f"Sample {sid!r}: control_bam file not found: {ctrl_bam}")

    return {
        "id": sid,
        "fq1": fq1,
        "fq2": fq2,
        "layout": lo,
        "assay": assay,
        "target": tgt,
        "peak_mode": pmode,
        "genome": gnm,
        "bt2_idx": bt2,
        "control_bam": ctrl_bam,
        "role": role,
        "control_sample": ctrl_sample,
        "experiment": experiment,
        "condition": condition,
        "replicate": replicate,
        "biological_replicate": bio_rep,
        "technical_replicate": tech_rep,
    }
