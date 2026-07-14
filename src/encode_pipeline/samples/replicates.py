"""Cross-replicate validation for Stage 4b replicate-aware outputs."""

from encode_pipeline.errors import ValidationError

__all__ = [
    "validate_replicate_groups",
]


def _group_samples_by_experiment(samples):
    """Group samples into treatment and control dicts keyed by experiment."""
    exp_treatments: dict[str, list[dict]] = {}
    exp_controls: dict[str, list[dict]] = {}
    for s in samples:
        if s["role"] == "treatment":
            exp_treatments.setdefault(s["experiment"], []).append(s)
        else:
            exp_controls.setdefault(s["experiment"], []).append(s)
    return exp_treatments, exp_controls


def _biological_replicates(rows):
    """Return sorted biological replicate ids for *rows*."""
    return sorted({r["biological_replicate"] for r in rows})


def _validate_treatment_consistency(exp, rows, error_cls):
    """Check that treatment rows in one experiment agree on key fields."""
    for field in ("assay", "target", "condition", "genome", "peak_mode", "layout"):
        values = {r[field] for r in rows}
        if len(values) > 1:
            raise error_cls(
                f"Experiment {exp!r}: treatment samples disagree on "
                f"{field}: {sorted(values)}"
            )


def _validate_duplicate_replicate_combos(exp, rows, error_cls):
    """Check for duplicate (biological_replicate, technical_replicate) combos."""
    seen_bt: set[tuple[int, int]] = set()
    for r in rows:
        bt = (r["biological_replicate"], r["technical_replicate"])
        if bt in seen_bt:
            raise error_cls(
                f"Experiment {exp!r}: duplicate "
                f"(biological_replicate={bt[0]}, "
                f"technical_replicate={bt[1]}) combination"
            )
        seen_bt.add(bt)


def _validate_control_consistency(exp, rows, all_samples, error_cls):
    """Check control usage is consistent within one experiment's treatment rows."""
    has_bam = [r for r in rows if r.get("control_bam")]
    has_sample = [r for r in rows if r.get("control_sample")]

    # Partial controls: not all have controls and not all lack them
    n_with_ctrl = len(has_bam) + len(has_sample)
    if n_with_ctrl not in (0, len(rows)):
        raise error_cls(
            f"Experiment {exp!r}: partial controls detected — "
            f"{n_with_ctrl} of {len(rows)} treatment rows have "
            f"controls. All treatment rows in an experiment must "
            f"either all have controls or all lack controls."
        )

    # Mixed control types
    if has_bam and has_sample:
        raise error_cls(
            f"Experiment {exp!r}: mixed control types — some rows "
            f"use control_bam and others use control_sample. "
            f"All treatment rows in an experiment must use the "
            f"same control type."
        )

    # Control sample references must belong to the same experiment
    if has_sample:
        for r in has_sample:
            cs = r["control_sample"]
            cs_row = next(s for s in all_samples if s["id"] == cs)
            if cs_row["experiment"] != exp:
                raise error_cls(
                    f"Experiment {exp!r}: sample {r['id']!r} references "
                    f"control_sample {cs!r} which belongs to experiment "
                    f"{cs_row['experiment']!r}, not {exp!r}"
                )


def validate_replicate_groups(samples, use_control, stage5_enabled=False,
                              reproducibility_idr_atac_narrow=False,
                              reproducibility_idr_cuttag_narrow=False,
                              reproducibility_idr_chipseq_broad=False,
                              reproducibility_idr_cuttag_broad=False,
                              error_cls=ValidationError):
    """Pass 3: cross-replicate validation for Stage 4b replicate-aware outputs.

    Raises *error_cls* on:
    - Inconsistent assay/target/genome/peak_mode/layout within an experiment
      (treatment samples only)
    - Duplicate (biological_replicate, technical_replicate) combos within
      an experiment
    - Mixed control types (control_bam vs control_sample) within an experiment
    - Partial controls: some treatment rows have controls and others do not
    """
    exp_treatments, exp_controls = _group_samples_by_experiment(samples)
    # exp_controls is reserved for future experiment-scoped control checks.
    _ = exp_controls

    for exp, rows in exp_treatments.items():
        # --- Consistency checks ---
        if len(rows) < 2:
            continue

        _validate_treatment_consistency(exp, rows, error_cls)
        _validate_duplicate_replicate_combos(exp, rows, error_cls)

        # --- Control consistency (only when use_control is true) ---
        if not use_control:
            continue

        _validate_control_consistency(exp, rows, samples, error_cls)

    # --- Stage 5 IDR eligibility (only when stage5_enabled) ---
    if stage5_enabled:
        chipseq_narrow_exps = []
        for exp, rows in exp_treatments.items():
            if len(rows) == 0:
                continue
            first = rows[0]

            if first["assay"] != "chipseq":
                continue       # skip non-chipseq silently

            if first["peak_mode"] != "narrow":
                continue       # skip chipseq broad

            bio_reps = _biological_replicates(rows)
            if len(bio_reps) != 2:
                raise error_cls(
                    f"Stage 5 IDR: experiment {exp!r} has "
                    f"{len(bio_reps)} biological replicate(s) ({bio_reps}). "
                    f"Stage 5 requires exactly 2."
                )
            chipseq_narrow_exps.append(exp)

        if not chipseq_narrow_exps:
            raise error_cls(
                "stage5=true but no eligible ChIP-seq narrow experiments "
                "were found."
            )

    # --- reproducibility.idr.atac_narrow eligibility (Stage 55) ---
    if reproducibility_idr_atac_narrow:
        atac_narrow_exps = []
        for exp, rows in exp_treatments.items():
            if len(rows) == 0:
                continue
            first = rows[0]
            if first["assay"] != "atac":
                continue       # skip non-ATAC silently
            if first["peak_mode"] != "narrow":
                continue       # ATAC broad — no IDR
            bio_reps = _biological_replicates(rows)
            if len(bio_reps) != 2:
                raise error_cls(
                    f"reproducibility.idr.atac_narrow: ATAC narrow "
                    f"experiment {exp!r} has {len(bio_reps)} biological "
                    f"replicate(s). IDR requires exactly 2."
                )
            atac_narrow_exps.append(exp)

        if not atac_narrow_exps:
            raise error_cls(
                "reproducibility.idr.atac_narrow is true but no eligible "
                "ATAC narrow experiments with exactly 2 biological "
                "replicates were found."
            )

    # --- reproducibility.idr.cuttag_narrow eligibility (Stage 64) ---
    if reproducibility_idr_cuttag_narrow:
        cuttag_narrow_exps = []
        for exp, rows in exp_treatments.items():
            if len(rows) == 0:
                continue
            first = rows[0]

            if first["assay"] != "cuttag":
                continue       # skip non-CUT&Tag silently

            if first["peak_mode"] != "narrow":
                continue       # skip CUT&Tag broad

            bio_reps = _biological_replicates(rows)
            if len(bio_reps) != 2:
                raise error_cls(
                    f"reproducibility.idr.cuttag_narrow: CUT&Tag narrow "
                    f"experiment {exp!r} has {len(bio_reps)} biological "
                    f"replicate(s). IDR requires exactly 2."
                )
            cuttag_narrow_exps.append(exp)

        if not cuttag_narrow_exps:
            raise error_cls(
                "reproducibility.idr.cuttag_narrow is true but no eligible "
                "CUT&Tag narrow experiments with exactly 2 biological "
                "replicates were found."
            )

    # --- reproducibility.idr.chipseq_broad_experimental eligibility (Stage 65) ---
    if reproducibility_idr_chipseq_broad:
        chipseq_broad_exps = []
        for exp, rows in exp_treatments.items():
            if len(rows) == 0:
                continue
            first = rows[0]

            if first["assay"] != "chipseq":
                continue       # skip non-chipseq silently

            if first["peak_mode"] != "broad":
                continue       # skip chipseq narrow

            bio_reps = _biological_replicates(rows)
            if len(bio_reps) != 2:
                raise error_cls(
                    f"reproducibility.idr.chipseq_broad_experimental: "
                    f"ChIP-seq broad experiment {exp!r} has "
                    f"{len(bio_reps)} biological replicate(s). "
                    f"IDR requires exactly 2."
                )
            chipseq_broad_exps.append(exp)

        if not chipseq_broad_exps:
            raise error_cls(
                "reproducibility.idr.chipseq_broad_experimental is true "
                "but no eligible ChIP-seq broad experiments with exactly "
                "2 biological replicates were found."
            )

    # --- reproducibility.idr.cuttag_broad_experimental eligibility (Stage 65) ---
    if reproducibility_idr_cuttag_broad:
        cuttag_broad_exps = []
        for exp, rows in exp_treatments.items():
            if len(rows) == 0:
                continue
            first = rows[0]

            if first["assay"] != "cuttag":
                continue       # skip non-CUT&Tag silently

            if first["peak_mode"] != "broad":
                continue       # skip CUT&Tag narrow

            bio_reps = _biological_replicates(rows)
            if len(bio_reps) != 2:
                raise error_cls(
                    f"reproducibility.idr.cuttag_broad_experimental: "
                    f"CUT&Tag broad experiment {exp!r} has "
                    f"{len(bio_reps)} biological replicate(s). "
                    f"IDR requires exactly 2."
                )
            cuttag_broad_exps.append(exp)

        if not cuttag_broad_exps:
            raise error_cls(
                "reproducibility.idr.cuttag_broad_experimental is true "
                "but no eligible CUT&Tag broad experiments with exactly "
                "2 biological replicates were found."
            )
