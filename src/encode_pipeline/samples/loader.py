"""Sample sheet loading orchestration."""

import csv
import os

from encode_pipeline.config import defaults
from encode_pipeline.errors import ValidationError
from encode_pipeline.samples import replicates as replicates_validation
from encode_pipeline.samples import rows as rows_validation
from encode_pipeline.samples import strict as strict_validation

__all__ = [
    "load_and_validate_samples",
]


def load_and_validate_samples(
    sample_tsv: str,
    *,
    use_control: bool = False,
    stage5_enabled: bool = False,
    strict_inputs: bool = False,
    reproducibility_idr_atac_narrow: bool = False,
    reproducibility_idr_cuttag_narrow: bool = False,
    reproducibility_idr_chipseq_broad: bool = False,
    reproducibility_idr_cuttag_broad: bool = False,
) -> list[dict]:
    """Load and validate a sample sheet TSV.

    Returns a list of sample dicts. Raises ValidationError on invalid input.

    Sample dict keys:
        id, fq1, fq2, layout, assay, target, peak_mode, genome,
        bt2_idx, control_bam, role, control_sample,
        experiment, condition, replicate,
        biological_replicate, technical_replicate
    """
    if not os.path.isfile(sample_tsv):
        raise ValidationError(f"Sample sheet not found: {sample_tsv}")

    samples: list[dict] = []
    seen_ids: set[str] = set()

    with open(sample_tsv) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fieldnames = reader.fieldnames or []

        # Required columns
        required = list(defaults.SAMPLE_REQUIRED_COLUMNS)
        for col in required:
            if col not in fieldnames:
                raise ValidationError(f"Sample sheet missing required column: {col!r}")

        for i, row in enumerate(reader, start=2):
            sid = (row.get("sample") or "").strip()
            sample = rows_validation.validate_and_build_sample(
                row, i, use_control=use_control, error_cls=ValidationError
            )

            if sid in seen_ids:
                raise ValidationError(f"Duplicate sample ID in sample sheet: {sid!r}")
            seen_ids.add(sid)
            samples.append(sample)

    # --- Pass 2: cross-reference validation ---
    # Control references are intentionally ignored when use_control is false.
    if use_control:
        all_ids = {s["id"] for s in samples}
        for s in samples:
            cs = s.get("control_sample", "")
            if cs and cs not in all_ids:
                raise ValidationError(
                    f"Sample {s['id']!r}: control_sample {cs!r} "
                    f"not found in sample sheet"
                )
            if cs:
                cs_row = next(r for r in samples if r["id"] == cs)
                if cs_row["role"] != "control":
                    raise ValidationError(
                        f"Sample {s['id']!r}: control_sample {cs!r} "
                        f"has role={cs_row['role']!r}, expected role=control"
                    )

    # --- Pass 3: replicate group validation (Stage 4b) ---
    replicates_validation.validate_replicate_groups(
        samples,
        use_control,
        stage5_enabled,
        reproducibility_idr_atac_narrow=reproducibility_idr_atac_narrow,
        reproducibility_idr_cuttag_narrow=reproducibility_idr_cuttag_narrow,
        reproducibility_idr_chipseq_broad=reproducibility_idr_chipseq_broad,
        reproducibility_idr_cuttag_broad=reproducibility_idr_cuttag_broad,
        error_cls=ValidationError,
    )

    # --- Pass 4: strict input validation (Stage 30, optional) ---
    strict_validation.validate_strict_inputs(
        samples, strict_inputs, error_cls=ValidationError
    )

    return samples
