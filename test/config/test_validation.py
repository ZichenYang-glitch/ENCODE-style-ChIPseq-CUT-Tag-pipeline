"""Pytest-based validation tests for config/sample sheets.

This replaces test/test_validation_stress.py with native pytest assertions
and shared conftest.py fixtures.
"""

import pytest


VALID_SAMPLES_HEADER = (
    "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\tpeak_mode\tgenome\tbowtie2_index\n"
)


def _valid_samples(sid="S1"):
    return VALID_SAMPLES_HEADER + f"{sid}\tR1.fq\tR2.fq\tPE\tchipseq\tT\tnarrow\ths\tidx\n"


def _base_config():
    return {"use_control": False}


@pytest.fixture
def make_validation_case(tmp_config):
    """Return a helper that builds a validation case and returns the result."""

    def _make(config=None, samples=None, placeholders=None):
        cfg = _base_config()
        if config:
            cfg.update(config)
        return tmp_config(
            config=cfg,
            samples=samples or _valid_samples(),
            placeholders=placeholders,
        )

    return _make


# ---------------------------------------------------------------------------
# Baseline
# ---------------------------------------------------------------------------


def test_baseline_valid(make_validation_case, run_validator):
    """A minimal valid config and sample sheet passes validation."""
    _, config_path, _ = make_validation_case()
    result = run_validator(config_path)
    assert result.returncode == 0, f"Expected pass, got:\n{result.stderr}"
    assert "OK:" in result.stdout


# ---------------------------------------------------------------------------
# Configuration validation
# ---------------------------------------------------------------------------


def test_threads_not_int_rejected(make_validation_case, run_validator):
    _, config_path, _ = make_validation_case(config={"threads": "abc"})
    result = run_validator(config_path)
    assert result.returncode != 0
    assert "must be an integer" in (result.stdout + result.stderr)


def test_threads_zero_rejected(make_validation_case, run_validator):
    _, config_path, _ = make_validation_case(config={"threads": 0})
    result = run_validator(config_path)
    assert result.returncode != 0
    assert "must be positive" in (result.stdout + result.stderr)


def test_mapq_negative_rejected(make_validation_case, run_validator):
    _, config_path, _ = make_validation_case(config={"mapq": -5})
    result = run_validator(config_path)
    assert result.returncode != 0
    assert "must be non-negative" in (result.stdout + result.stderr)


def test_invalid_remove_dup_rejected(make_validation_case, run_validator):
    _, config_path, _ = make_validation_case(config={"remove_dup": "whatever"})
    result = run_validator(config_path)
    assert result.returncode != 0
    assert "auto, yes, or no" in (result.stdout + result.stderr)


# ---------------------------------------------------------------------------
# Genome resources
# ---------------------------------------------------------------------------


def test_missing_genome_resources_egs_rejected(make_validation_case, run_validator):
    _, config_path, _ = make_validation_case(
        config={"genome_resources": {"hs": {"chrom_sizes": ""}}},
    )
    result = run_validator(config_path)
    assert result.returncode != 0
    assert "effective_genome_size is required" in (result.stdout + result.stderr)


def test_invalid_genome_resources_egs_rejected(make_validation_case, run_validator):
    _, config_path, _ = make_validation_case(
        config={"genome_resources": {"hs": {"effective_genome_size": "hg38"}}},
    )
    result = run_validator(config_path)
    assert result.returncode != 0
    assert "'hs', 'mm', or a positive integer" in (result.stdout + result.stderr)


# ---------------------------------------------------------------------------
# Sample sheet structure
# ---------------------------------------------------------------------------


def test_missing_fastq_1_column_rejected(make_validation_case, run_validator):
    samples = "sample\tfastq_2\tlayout\nS1\tR2\tSE\n"
    _, config_path, _ = make_validation_case(samples=samples)
    result = run_validator(config_path)
    assert result.returncode != 0
    assert "missing required column: 'fastq_1'" in (result.stdout + result.stderr)


def test_invalid_sample_id_rejected(make_validation_case, run_validator):
    samples = VALID_SAMPLES_HEADER + "S1 Bad\tR1\tR2\tPE\tchipseq\tT\tnarrow\ths\tidx\n"
    _, config_path, _ = make_validation_case(samples=samples)
    result = run_validator(config_path)
    assert result.returncode != 0
    assert "invalid characters" in (result.stdout + result.stderr)


def test_pe_missing_fastq_2_rejected(make_validation_case, run_validator):
    samples = VALID_SAMPLES_HEADER + "S1\tR1\t\tPE\tchipseq\tT\tnarrow\ths\tidx\n"
    _, config_path, _ = make_validation_case(samples=samples)
    result = run_validator(config_path)
    assert result.returncode != 0
    assert "PE layout requires 'fastq_2'" in (result.stdout + result.stderr)


def test_invalid_assay_rejected(make_validation_case, run_validator):
    samples = VALID_SAMPLES_HEADER + "S1\tR1\t\tSE\trnaseq\tT\tnarrow\ths\tidx\n"
    _, config_path, _ = make_validation_case(samples=samples)
    result = run_validator(config_path)
    assert result.returncode != 0
    assert "assay must be chipseq, cuttag, atac, or mnase" in (result.stdout + result.stderr)


# ---------------------------------------------------------------------------
# Control logic
# ---------------------------------------------------------------------------


def test_control_referencing_missing_sample_rejected(make_validation_case, run_validator):
    header = (
        "sample\tfastq_1\tlayout\tassay\ttarget\tpeak_mode\tgenome\tbowtie2_index\t"
        "role\tcontrol_sample\n"
    )
    samples = header + "S1\tR1\tSE\tchipseq\tT\tnarrow\ths\tidx\ttreatment\tCtrl1\n"
    _, config_path, _ = make_validation_case(
        config={"use_control": True},
        samples=samples,
    )
    result = run_validator(config_path)
    assert result.returncode != 0
    assert "not found in sample sheet" in (result.stdout + result.stderr)


def test_control_bad_role_rejected(make_validation_case, run_validator):
    header = (
        "sample\tfastq_1\tlayout\tassay\ttarget\tpeak_mode\tgenome\tbowtie2_index\t"
        "role\tcontrol_sample\n"
    )
    samples = (
        header
        + "S1\tR1\tSE\tchipseq\tT\tnarrow\ths\tidx\ttreatment\tCtrl1\n"
        + "Ctrl1\tR1\tSE\tchipseq\tT\tnarrow\ths\tidx\ttreatment\t\n"
    )
    _, config_path, _ = make_validation_case(
        config={"use_control": True},
        samples=samples,
    )
    result = run_validator(config_path)
    assert result.returncode != 0
    assert "expected role=control" in (result.stdout + result.stderr)


def test_mutually_exclusive_controls_rejected(make_validation_case, run_validator):
    header = (
        "sample\tfastq_1\tlayout\tassay\ttarget\tpeak_mode\tgenome\tbowtie2_index\t"
        "role\tcontrol_sample\tcontrol_bam\n"
    )
    samples = (
        header
        + "S1\tR1\tSE\tchipseq\tT\tnarrow\ths\tidx\ttreatment\tCtrl1\t/tmp/a.bam\n"
        + "Ctrl1\tR1\tSE\tchipseq\tT\tnarrow\ths\tidx\tcontrol\t\t\n"
    )
    _, config_path, _ = make_validation_case(
        config={"use_control": True},
        samples=samples,
    )
    result = run_validator(config_path)
    assert result.returncode != 0
    assert "cannot set both control_sample and control_bam" in (result.stdout + result.stderr)


def test_missing_control_bam_rejected(make_validation_case, run_validator):
    header = (
        "sample\tfastq_1\tlayout\tassay\ttarget\tpeak_mode\tgenome\tbowtie2_index\t"
        "role\tcontrol_bam\n"
    )
    samples = header + "S1\tR1\tSE\tchipseq\tT\tnarrow\ths\tidx\ttreatment\t/does/not/exist.bam\n"
    _, config_path, _ = make_validation_case(
        config={"use_control": True},
        samples=samples,
    )
    result = run_validator(config_path)
    assert result.returncode != 0
    assert "control_bam file not found" in (result.stdout + result.stderr)


# ---------------------------------------------------------------------------
# MNase validation
# ---------------------------------------------------------------------------


def _mnase_samples(layout="PE", peak_mode="nucleosome", sid="M1"):
    return VALID_SAMPLES_HEADER + f"{sid}\tR1.fq\tR2.fq\t{layout}\tmnase\tH3\t{peak_mode}\ths\tidx\n"


def test_mnase_pe_nucleosome_valid(make_validation_case, run_validator):
    _, config_path, _ = make_validation_case(samples=_mnase_samples())
    result = run_validator(config_path)
    assert result.returncode == 0, f"Expected pass, got:\n{result.stderr}"


def test_mnase_se_rejected(make_validation_case, run_validator):
    _, config_path, _ = make_validation_case(samples=_mnase_samples(layout="SE"))
    result = run_validator(config_path)
    assert result.returncode != 0
    assert "requires paired-end layout" in (result.stdout + result.stderr)


def test_mnase_narrow_peak_mode_rejected(make_validation_case, run_validator):
    _, config_path, _ = make_validation_case(samples=_mnase_samples(peak_mode="narrow"))
    result = run_validator(config_path)
    assert result.returncode != 0
    assert "requires peak_mode=nucleosome" in (result.stdout + result.stderr)


def test_chipseq_nucleosome_rejected(make_validation_case, run_validator):
    samples = VALID_SAMPLES_HEADER + "S1\tR1.fq\tR2.fq\tPE\tchipseq\tT\tnucleosome\ths\tidx\n"
    _, config_path, _ = make_validation_case(samples=samples)
    result = run_validator(config_path)
    assert result.returncode != 0
    assert "peak_mode=nucleosome is only allowed for assay=mnase" in (result.stdout + result.stderr)


@pytest.mark.parametrize(
    "mono_range,expected_error",
    [
        ([200, 100], "min must be < max"),
        ([140], "exactly 2 elements"),
    ],
)
def test_mnase_mono_range_invalid(make_validation_case, run_validator, mono_range, expected_error):
    _, config_path, _ = make_validation_case(
        config={"mnase": {"mono_range": mono_range}},
        samples=_mnase_samples(),
    )
    result = run_validator(config_path)
    assert result.returncode != 0
    assert expected_error in (result.stdout + result.stderr)


def test_mnase_mono_range_valid(make_validation_case, run_validator):
    _, config_path, _ = make_validation_case(
        config={"mnase": {"mono_range": [100, 180]}},
        samples=_mnase_samples(),
    )
    result = run_validator(config_path)
    assert result.returncode == 0, f"Expected pass, got:\n{result.stderr}"


# ---------------------------------------------------------------------------
# MNase Stage 40 validation
# ---------------------------------------------------------------------------


def test_mnase_fragments_valid(make_validation_case, run_validator):
    _, config_path, _ = make_validation_case(
        config={
            "mnase": {
                "fragments": {
                    "sub": [1, 139],
                    "mono": [140, 200],
                    "di": [300, 400],
                },
            },
        },
        samples=_mnase_samples(),
    )
    result = run_validator(config_path)
    assert result.returncode == 0, f"Expected pass, got:\n{result.stderr}"


def test_mnase_fragments_mono_overrides_mono_range(make_validation_case, run_validator):
    _, config_path, _ = make_validation_case(
        config={
            "mnase": {
                "mono_range": [100, 180],
                "fragments": {"mono": [140, 200]},
            },
        },
        samples=_mnase_samples(),
    )
    result = run_validator(config_path)
    assert result.returncode == 0, f"Expected pass, got:\n{result.stderr}"


@pytest.mark.parametrize(
    "fragments,expected_error",
    [
        ({"sub": [200, 100]}, "min must be < max"),
        ({"mono": [0, 200]}, "values must be positive"),
    ],
)
def test_mnase_fragments_invalid(make_validation_case, run_validator, fragments, expected_error):
    _, config_path, _ = make_validation_case(
        config={"mnase": {"fragments": fragments}},
        samples=_mnase_samples(),
    )
    result = run_validator(config_path)
    assert result.returncode != 0
    assert expected_error in (result.stdout + result.stderr)


@pytest.mark.parametrize(
    "dyad_range,expected_error",
    [
        ([100], "exactly 2 elements"),
        (["a", 200], "must be integers"),
        ([130.5, 200], "must be integers"),
    ],
)
def test_mnase_dyad_range_invalid(make_validation_case, run_validator, dyad_range, expected_error):
    _, config_path, _ = make_validation_case(
        config={"mnase": {"dyad_range": dyad_range}},
        samples=_mnase_samples(),
    )
    result = run_validator(config_path)
    assert result.returncode != 0
    assert expected_error in (result.stdout + result.stderr)


def test_mnase_dyad_range_valid(make_validation_case, run_validator):
    _, config_path, _ = make_validation_case(
        config={"mnase": {"dyad_range": [140, 200]}},
        samples=_mnase_samples(),
    )
    result = run_validator(config_path)
    assert result.returncode == 0, f"Expected pass, got:\n{result.stderr}"


def test_mnase_dyad_range_defaults_ok(make_validation_case, run_validator):
    _, config_path, _ = make_validation_case(
        config={"mnase": {"mono_range": [140, 200]}},
        samples=_mnase_samples(),
    )
    result = run_validator(config_path)
    assert result.returncode == 0, f"Expected pass, got:\n{result.stderr}"


@pytest.mark.parametrize(
    "sub_value",
    [True, "true"],
)
def test_mnase_fragments_sub_bool_rejected(make_validation_case, run_validator, sub_value):
    _, config_path, _ = make_validation_case(
        config={"mnase": {"fragments": {"sub": [sub_value, 139]}}},
        samples=_mnase_samples(),
    )
    result = run_validator(config_path)
    assert result.returncode != 0
    assert "must be integers" in (result.stdout + result.stderr)


@pytest.mark.parametrize(
    "callers,expected_error",
    [
        ({"danpos3": True}, "not implemented"),
        ({"inps": True}, "not implemented"),
        ({"sem": True}, "not implemented"),
        ({"danpos3": "false"}, "must be boolean"),
        ({"unknown": False}, "unknown key"),
    ],
)
def test_mnase_callers_invalid(make_validation_case, run_validator, callers, expected_error):
    _, config_path, _ = make_validation_case(
        config={"mnase": {"callers": callers}},
        samples=_mnase_samples(),
    )
    result = run_validator(config_path)
    assert result.returncode != 0
    assert expected_error in (result.stdout + result.stderr)


def test_mnase_callers_all_false_valid(make_validation_case, run_validator):
    _, config_path, _ = make_validation_case(
        config={
            "mnase": {
                "callers": {"danpos3": False, "inps": False, "sem": False},
            },
        },
        samples=_mnase_samples(),
    )
    result = run_validator(config_path)
    assert result.returncode == 0, f"Expected pass, got:\n{result.stderr}"


def test_mnase_se_still_rejected(make_validation_case, run_validator):
    _, config_path, _ = make_validation_case(samples=_mnase_samples(layout="SE"))
    result = run_validator(config_path)
    assert result.returncode != 0
    assert "requires paired-end layout" in (result.stdout + result.stderr)


def test_mnase_narrow_still_rejected(make_validation_case, run_validator):
    _, config_path, _ = make_validation_case(samples=_mnase_samples(peak_mode="narrow"))
    result = run_validator(config_path)
    assert result.returncode != 0
    assert "requires peak_mode=nucleosome" in (result.stdout + result.stderr)
