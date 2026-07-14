"""Direct-API characterization tests for replicate-group validation.

These tests pin the behavior of ``encode_pipeline.samples.replicates.validate_replicate_groups``
so future refactors can safely move or split the function without changing behavior.
"""

import pytest

from encode_pipeline.errors import ValidationError
from encode_pipeline.samples.replicates import validate_replicate_groups


def _sample(
    *,
    sid="S1",
    experiment="EXP1",
    role="treatment",
    assay="chipseq",
    target="T",
    condition="T",
    genome="hs",
    peak_mode="narrow",
    layout="PE",
    biological_replicate=1,
    technical_replicate=1,
    control_sample="",
    control_bam="",
):
    return {
        "id": sid,
        "experiment": experiment,
        "role": role,
        "assay": assay,
        "target": target,
        "condition": condition,
        "genome": genome,
        "peak_mode": peak_mode,
        "layout": layout,
        "biological_replicate": biological_replicate,
        "technical_replicate": technical_replicate,
        "control_sample": control_sample,
        "control_bam": control_bam,
    }


def _treatment_group(
    *,
    experiment,
    assay,
    peak_mode,
    biological_replicates=(1, 2),
    layout="PE",
):
    return [
        _sample(
            sid=f"{experiment}_{biological_replicate}",
            experiment=experiment,
            assay=assay,
            peak_mode=peak_mode,
            layout=layout,
            biological_replicate=biological_replicate,
        )
        for biological_replicate in biological_replicates
    ]


# ---------------------------------------------------------------------------
# Baseline / grouping
# ---------------------------------------------------------------------------


def test_two_biological_replicates_same_experiment_passes():
    samples = [
        _sample(sid="S1", biological_replicate=1),
        _sample(sid="S2", biological_replicate=2),
    ]
    assert validate_replicate_groups(samples, use_control=False) is None


def test_controls_are_excluded_from_treatment_consistency_checks():
    samples = [
        _sample(sid="S1", biological_replicate=1),
        _sample(
            sid="C1",
            role="control",
            assay="atac",
            target="IgG",
            condition="IgG",
            genome="mm",
            peak_mode="broad",
            layout="SE",
        ),
    ]
    assert validate_replicate_groups(samples, use_control=False) is None


def test_multiple_experiments_validated_independently():
    samples = [
        _sample(sid="S1", experiment="EXP1", biological_replicate=1),
        _sample(sid="S2", experiment="EXP1", biological_replicate=2),
        _sample(
            sid="S3",
            experiment="EXP2",
            assay="atac",
            biological_replicate=1,
        ),
        _sample(
            sid="S4",
            experiment="EXP2",
            assay="atac",
            biological_replicate=2,
        ),
    ]
    assert validate_replicate_groups(samples, use_control=False) is None


def test_single_treatment_experiment_skips_consistency_checks():
    samples = [
        _sample(sid="S1", experiment="EXP1"),
    ]
    assert validate_replicate_groups(samples, use_control=False) is None


# ---------------------------------------------------------------------------
# Consistency fields
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "field,value1,value2",
    [
        ("assay", "chipseq", "atac"),
        ("target", "T1", "T2"),
        ("condition", "C1", "C2"),
        ("genome", "hs", "mm"),
        ("peak_mode", "narrow", "broad"),
        ("layout", "PE", "SE"),
    ],
)
def test_inconsistent_treatment_field_rejected(field, value1, value2):
    s1 = _sample(sid="S1", biological_replicate=1, **{field: value1})
    s2 = _sample(sid="S2", biological_replicate=2, **{field: value2})
    with pytest.raises(ValidationError, match=f"disagree on {field}"):
        validate_replicate_groups([s1, s2], use_control=False)


# ---------------------------------------------------------------------------
# Duplicate replicate combinations
# ---------------------------------------------------------------------------


def test_duplicate_biological_technical_replicate_combo_rejected():
    samples = [
        _sample(sid="S1", biological_replicate=1, technical_replicate=1),
        _sample(sid="S2", biological_replicate=1, technical_replicate=1),
    ]
    with pytest.raises(
        ValidationError,
        match="duplicate .*biological_replicate=1.*technical_replicate=1",
    ):
        validate_replicate_groups(samples, use_control=False)


# ---------------------------------------------------------------------------
# use_control=False ignores control consistency
# ---------------------------------------------------------------------------


def test_use_control_false_ignores_partial_controls():
    samples = [
        _sample(sid="S1", biological_replicate=1, control_sample="C1"),
        _sample(sid="S2", biological_replicate=2),
    ]
    assert validate_replicate_groups(samples, use_control=False) is None


def test_use_control_false_ignores_mixed_control_types():
    samples = [
        _sample(sid="S1", biological_replicate=1, control_sample="C1"),
        _sample(sid="S2", biological_replicate=2, control_bam="/tmp/a.bam"),
    ]
    assert validate_replicate_groups(samples, use_control=False) is None


# ---------------------------------------------------------------------------
# use_control=True control consistency
# ---------------------------------------------------------------------------


def test_use_control_true_valid_control_sample_passes():
    samples = [
        _sample(
            sid="S1",
            experiment="EXP1",
            biological_replicate=1,
            control_sample="C1",
        ),
        _sample(
            sid="S2",
            experiment="EXP1",
            biological_replicate=2,
            control_sample="C1",
        ),
        _sample(
            sid="C1",
            experiment="EXP1",
            role="control",
        ),
    ]
    assert validate_replicate_groups(samples, use_control=True) is None


def test_use_control_true_all_control_bam_passes():
    samples = [
        _sample(
            sid="S1",
            biological_replicate=1,
            control_bam="/path/to/control.bam",
        ),
        _sample(
            sid="S2",
            biological_replicate=2,
            control_bam="/path/to/control.bam",
        ),
    ]
    assert validate_replicate_groups(samples, use_control=True) is None


def test_use_control_true_partial_controls_rejected():
    samples = [
        _sample(sid="S1", biological_replicate=1, control_sample="C1"),
        _sample(sid="S2", biological_replicate=2),
    ]
    with pytest.raises(ValidationError, match="partial controls detected"):
        validate_replicate_groups(samples, use_control=True)


def test_use_control_true_mixed_control_types_rejected():
    samples = [
        _sample(sid="S1", biological_replicate=1, control_sample="C1"),
        _sample(sid="S2", biological_replicate=2, control_bam="/tmp/a.bam"),
    ]
    with pytest.raises(ValidationError, match="mixed control types"):
        validate_replicate_groups(samples, use_control=True)


def test_control_sample_must_belong_to_same_experiment():
    samples = [
        _sample(
            sid="S1",
            experiment="EXP1",
            biological_replicate=1,
            control_sample="C1",
        ),
        _sample(
            sid="S2",
            experiment="EXP1",
            biological_replicate=2,
            control_sample="C1",
        ),
        _sample(
            sid="C1",
            experiment="EXP2",
            role="control",
        ),
    ]
    with pytest.raises(
        ValidationError,
        match="belongs to experiment 'EXP2', not 'EXP1'",
    ):
        validate_replicate_groups(samples, use_control=True)


# ---------------------------------------------------------------------------
# Stage 5 IDR eligibility
# ---------------------------------------------------------------------------


def test_chipseq_narrow_idr_two_bioreps_passes():
    samples = [
        _sample(sid="S1", biological_replicate=1),
        _sample(sid="S2", biological_replicate=2),
    ]
    assert (
        validate_replicate_groups(samples, use_control=False, stage5_enabled=True)
        is None
    )


@pytest.mark.parametrize(
    "assay,peak_mode",
    [
        ("chipseq", "broad"),
        ("atac", "narrow"),
        ("cuttag", "narrow"),
    ],
)
def test_chipseq_narrow_idr_skips_other_assay_modes(assay, peak_mode):
    samples = [
        _sample(sid="S1", assay=assay, peak_mode=peak_mode, biological_replicate=1),
        _sample(sid="S2", assay=assay, peak_mode=peak_mode, biological_replicate=2),
    ]
    with pytest.raises(
        ValidationError,
        match="stage5=true but no eligible ChIP-seq narrow experiments",
    ):
        validate_replicate_groups(samples, use_control=False, stage5_enabled=True)


def test_chipseq_narrow_idr_wrong_biorep_count_rejected():
    samples = [
        _sample(sid="S1", biological_replicate=1),
        _sample(sid="S2", biological_replicate=1, technical_replicate=2),
    ]
    with pytest.raises(ValidationError, match="Stage 5 IDR"):
        validate_replicate_groups(samples, use_control=False, stage5_enabled=True)


def test_chipseq_narrow_idr_without_eligible_experiment_rejected():
    samples = [
        _sample(sid="S1", assay="atac", peak_mode="narrow", biological_replicate=1),
        _sample(sid="S2", assay="atac", peak_mode="narrow", biological_replicate=2),
    ]
    with pytest.raises(
        ValidationError,
        match="stage5=true but no eligible ChIP-seq narrow experiments",
    ):
        validate_replicate_groups(samples, use_control=False, stage5_enabled=True)


# ---------------------------------------------------------------------------
# Reproducibility IDR eligibility gates
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "flag,assay,peak_mode,label",
    [
        ("reproducibility_idr_atac_narrow", "atac", "narrow", "ATAC narrow"),
        ("reproducibility_idr_cuttag_narrow", "cuttag", "narrow", "CUT&Tag narrow"),
        ("reproducibility_idr_chipseq_broad", "chipseq", "broad", "ChIP-seq broad"),
        ("reproducibility_idr_cuttag_broad", "cuttag", "broad", "CUT&Tag broad"),
    ],
)
def test_reproducibility_idr_gate_valid_eligible_experiment_passes(
    flag, assay, peak_mode, label
):
    samples = [
        _sample(sid="S1", assay=assay, peak_mode=peak_mode, biological_replicate=1),
        _sample(sid="S2", assay=assay, peak_mode=peak_mode, biological_replicate=2),
    ]
    kwargs = {flag: True}
    assert validate_replicate_groups(samples, use_control=False, **kwargs) is None


@pytest.mark.parametrize(
    "flag,assay,peak_mode,label",
    [
        ("reproducibility_idr_atac_narrow", "atac", "narrow", "ATAC narrow"),
        ("reproducibility_idr_cuttag_narrow", "cuttag", "narrow", "CUT&Tag narrow"),
        ("reproducibility_idr_chipseq_broad", "chipseq", "broad", "ChIP-seq broad"),
        ("reproducibility_idr_cuttag_broad", "cuttag", "broad", "CUT&Tag broad"),
    ],
)
def test_reproducibility_idr_gate_wrong_biorep_count_rejected(
    flag, assay, peak_mode, label
):
    samples = [
        _sample(sid="S1", assay=assay, peak_mode=peak_mode, biological_replicate=1),
        _sample(
            sid="S2",
            assay=assay,
            peak_mode=peak_mode,
            biological_replicate=1,
            technical_replicate=2,
        ),
    ]
    kwargs = {flag: True}
    with pytest.raises(ValidationError, match=label):
        validate_replicate_groups(samples, use_control=False, **kwargs)


@pytest.mark.parametrize(
    "flag,assay,peak_mode,expected_message",
    [
        (
            "reproducibility_idr_atac_narrow",
            "chipseq",
            "narrow",
            "no eligible ATAC narrow experiments",
        ),
        (
            "reproducibility_idr_cuttag_narrow",
            "chipseq",
            "narrow",
            "no eligible CUT&Tag narrow experiments",
        ),
        (
            "reproducibility_idr_chipseq_broad",
            "atac",
            "narrow",
            "no eligible ChIP-seq broad experiments",
        ),
        (
            "reproducibility_idr_cuttag_broad",
            "atac",
            "narrow",
            "no eligible CUT&Tag broad experiments",
        ),
    ],
)
def test_reproducibility_idr_gate_no_eligible_experiment_rejected(
    flag, assay, peak_mode, expected_message
):
    samples = [
        _sample(sid="S1", assay=assay, peak_mode=peak_mode, biological_replicate=1),
        _sample(sid="S2", assay=assay, peak_mode=peak_mode, biological_replicate=2),
    ]
    kwargs = {flag: True}
    with pytest.raises(ValidationError, match=expected_message):
        validate_replicate_groups(samples, use_control=False, **kwargs)


@pytest.mark.parametrize(
    "flag,wrong_assay,wrong_peak_mode",
    [
        ("reproducibility_idr_atac_narrow", "chipseq", "narrow"),
        ("reproducibility_idr_cuttag_narrow", "chipseq", "narrow"),
        ("reproducibility_idr_chipseq_broad", "chipseq", "narrow"),
        ("reproducibility_idr_cuttag_broad", "cuttag", "narrow"),
    ],
)
def test_reproducibility_idr_gate_wrong_assay_or_mode_skipped(
    flag, wrong_assay, wrong_peak_mode
):
    samples = [
        _sample(
            sid="S1",
            assay=wrong_assay,
            peak_mode=wrong_peak_mode,
            biological_replicate=1,
        ),
        _sample(
            sid="S2",
            assay=wrong_assay,
            peak_mode=wrong_peak_mode,
            biological_replicate=2,
        ),
    ]
    kwargs = {flag: True}
    with pytest.raises(ValidationError, match="no eligible"):
        validate_replicate_groups(samples, use_control=False, **kwargs)


def test_all_idr_modes_validate_independently_in_one_mixed_sample_sheet():
    samples = [
        *_treatment_group(
            experiment="CHIP_NARROW",
            assay="chipseq",
            peak_mode="narrow",
        ),
        *_treatment_group(
            experiment="ATAC_NARROW",
            assay="atac",
            peak_mode="narrow",
        ),
        *_treatment_group(
            experiment="CUTTAG_NARROW",
            assay="cuttag",
            peak_mode="narrow",
        ),
        *_treatment_group(
            experiment="CHIP_BROAD",
            assay="chipseq",
            peak_mode="broad",
        ),
        *_treatment_group(
            experiment="CUTTAG_BROAD",
            assay="cuttag",
            peak_mode="broad",
        ),
        *_treatment_group(
            experiment="MNASE",
            assay="mnase",
            peak_mode="nucleosome",
        ),
    ]

    assert (
        validate_replicate_groups(
            samples,
            use_control=False,
            stage5_enabled=True,
            reproducibility_idr_atac_narrow=True,
            reproducibility_idr_cuttag_narrow=True,
            reproducibility_idr_chipseq_broad=True,
            reproducibility_idr_cuttag_broad=True,
        )
        is None
    )


def test_cuttag_narrow_idr_accepts_single_end_replicate_groups():
    samples = _treatment_group(
        experiment="CUTTAG_SE",
        assay="cuttag",
        peak_mode="narrow",
        layout="SE",
    )

    assert (
        validate_replicate_groups(
            samples,
            use_control=False,
            reproducibility_idr_cuttag_narrow=True,
        )
        is None
    )


@pytest.mark.parametrize(
    "flag,assay,eligible_mode,irrelevant_mode",
    [
        ("reproducibility_idr_atac_narrow", "atac", "narrow", "broad"),
        ("reproducibility_idr_cuttag_narrow", "cuttag", "narrow", "broad"),
        ("reproducibility_idr_chipseq_broad", "chipseq", "broad", "narrow"),
        ("reproducibility_idr_cuttag_broad", "cuttag", "broad", "narrow"),
    ],
)
def test_idr_modes_ignore_wrong_peak_mode_in_mixed_sample_sheets(
    flag, assay, eligible_mode, irrelevant_mode
):
    samples = [
        *_treatment_group(
            experiment="ELIGIBLE",
            assay=assay,
            peak_mode=eligible_mode,
        ),
        *_treatment_group(
            experiment="IRRELEVANT",
            assay=assay,
            peak_mode=irrelevant_mode,
            biological_replicates=(1, 2, 3),
        ),
    ]

    assert validate_replicate_groups(samples, use_control=False, **{flag: True}) is None


@pytest.mark.parametrize(
    "flag,assay,peak_mode,error_pattern",
    [
        (
            "reproducibility_idr_atac_narrow",
            "atac",
            "narrow",
            "ATAC narrow experiment 'INVALID'.*exactly 2",
        ),
        (
            "reproducibility_idr_cuttag_narrow",
            "cuttag",
            "narrow",
            "CUT&Tag narrow experiment 'INVALID'.*exactly 2",
        ),
        (
            "reproducibility_idr_chipseq_broad",
            "chipseq",
            "broad",
            "ChIP-seq broad experiment 'INVALID'.*exactly 2",
        ),
        (
            "reproducibility_idr_cuttag_broad",
            "cuttag",
            "broad",
            "CUT&Tag broad experiment 'INVALID'.*exactly 2",
        ),
    ],
)
def test_idr_modes_reject_any_invalid_eligible_experiment_in_mixed_sample_sheets(
    flag, assay, peak_mode, error_pattern
):
    samples = [
        *_treatment_group(
            experiment="VALID",
            assay=assay,
            peak_mode=peak_mode,
        ),
        *_treatment_group(
            experiment="INVALID",
            assay=assay,
            peak_mode=peak_mode,
            biological_replicates=(1, 2, 3),
        ),
    ]

    with pytest.raises(ValidationError, match=error_pattern):
        validate_replicate_groups(samples, use_control=False, **{flag: True})


# ---------------------------------------------------------------------------
# Exception type
# ---------------------------------------------------------------------------


def test_default_exception_type_is_validation_error():
    samples = [
        _sample(sid="S1", biological_replicate=1),
        _sample(sid="S2", biological_replicate=1),
    ]
    with pytest.raises(ValidationError):
        validate_replicate_groups(samples, use_control=False)


def test_custom_error_cls_is_respected():
    class CustomError(ValueError):
        pass

    samples = [
        _sample(sid="S1", biological_replicate=1),
        _sample(sid="S2", biological_replicate=1),
    ]
    with pytest.raises(CustomError):
        validate_replicate_groups(samples, use_control=False, error_cls=CustomError)
