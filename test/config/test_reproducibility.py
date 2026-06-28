"""Direct-API characterization tests for reproducibility/IDR validation."""

import re
from contextlib import nullcontext

import pytest

from encode_pipeline.config.validate import ValidationError, validate_config


SAMPLES_HEADER = (
    "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\tpeak_mode\tgenome\tbowtie2_index\n"
)
SAMPLES_ROW = "S1\tR1.fq\tR2.fq\tPE\tchipseq\tT\tnarrow\ths\tidx\n"


def _make_config(tmp_path, **overrides):
    samples_path = tmp_path / "samples.tsv"
    samples_path.write_text(SAMPLES_HEADER + SAMPLES_ROW, encoding="utf-8")
    config = {"samples": str(samples_path), "use_control": False}
    config.update(overrides)
    return config


def _validate(tmp_path, **overrides):
    return validate_config(_make_config(tmp_path, **overrides))


# ---------------------------------------------------------------------------
# Reproducibility structure and defaults
# ---------------------------------------------------------------------------


def test_reproducibility_omitted_defaults_disabled(tmp_path):
    validated = _validate(tmp_path)
    assert validated["reproducibility"] == {"enabled": False}
    assert validated["idr"] == {"threshold": 0.05, "rank": "p.value", "seed": 42}


def test_disabled_reproducibility_ignores_nested_keys(tmp_path):
    validated = _validate(
        tmp_path,
        reproducibility={
            "enabled": False,
            "unknown": "ignored",
            "consensus": "not-validated-while-disabled",
            "idr": "not-validated-while-disabled",
        },
    )
    assert validated["reproducibility"] == {"enabled": False}


def test_enabled_reproducibility_expands_defaults(tmp_path):
    validated = _validate(tmp_path, reproducibility={"enabled": True})
    assert validated["reproducibility"] == {
        "enabled": True,
        "consensus": {"enabled": True, "min_replicates": 2, "reciprocal_overlap": 0.5},
        "idr": {
            "chipseq_narrow": None,
            "atac_narrow": False,
            "cuttag_narrow": False,
            "chipseq_broad_experimental": False,
            "cuttag_broad_experimental": False,
        },
    }
    assert validated["idr"] == {"threshold": 0.05, "rank": "p.value", "seed": 42}


@pytest.mark.parametrize(
    "reproducibility,pattern",
    [
        (["bad"], "reproducibility must be a mapping"),
        ({"enabled": True, "bad": True}, "reproducibility: unknown key 'bad'"),
        (
            {"enabled": True, "consensus": "bad"},
            "reproducibility.consensus must be a mapping",
        ),
        (
            {"enabled": True, "consensus": {"bad": True}},
            "reproducibility.consensus: unknown key 'bad'",
        ),
        (
            {"enabled": True, "idr": "bad"},
            "reproducibility.idr must be a mapping",
        ),
        (
            {"enabled": True, "idr": {"bad": True}},
            "reproducibility.idr: unknown key 'bad'",
        ),
    ],
)
def test_reproducibility_structure_and_unknown_keys_rejected(
    tmp_path, reproducibility, pattern
):
    with pytest.raises(ValidationError, match=re.escape(pattern)):
        _validate(tmp_path, reproducibility=reproducibility)


# ---------------------------------------------------------------------------
# Consensus settings
# ---------------------------------------------------------------------------


def test_consensus_numeric_fields_normalize_strings(tmp_path):
    validated = _validate(
        tmp_path,
        reproducibility={
            "enabled": True,
            "consensus": {"min_replicates": "3", "reciprocal_overlap": "0.75"},
        },
    )
    assert validated["reproducibility"]["consensus"] == {
        "enabled": True,
        "min_replicates": 3,
        "reciprocal_overlap": 0.75,
    }


@pytest.mark.parametrize("raw", [True, 1, "1", "abc"])
def test_consensus_min_replicates_rejects_invalid_values(tmp_path, raw):
    with pytest.raises(
        ValidationError,
        match=(
            "reproducibility.consensus.min_replicates must be an integer|"
            "reproducibility.consensus.min_replicates must be >= 2"
        ),
    ):
        _validate(
            tmp_path,
            reproducibility={
                "enabled": True,
                "consensus": {"min_replicates": raw},
            },
        )


@pytest.mark.parametrize("raw", [True, 0, -0.1, 1.1, "abc"])
def test_consensus_reciprocal_overlap_rejects_invalid_values(tmp_path, raw):
    with pytest.raises(
        ValidationError,
        match=(
            "reproducibility.consensus.reciprocal_overlap must be a float|"
            "reproducibility.consensus.reciprocal_overlap must be in"
        ),
    ):
        _validate(
            tmp_path,
            reproducibility={
                "enabled": True,
                "consensus": {"reciprocal_overlap": raw},
            },
        )


# ---------------------------------------------------------------------------
# IDR settings
# ---------------------------------------------------------------------------


def test_idr_settings_normalize_when_stage5_enabled(tmp_path):
    validated = _validate(
        tmp_path,
        stage5=True,
        idr={"threshold": "0.01", "rank": "signal.value", "seed": "123"},
    )
    assert validated["idr"] == {
        "threshold": 0.01,
        "rank": "signal.value",
        "seed": 123,
    }


def test_idr_settings_normalize_when_reproducibility_idr_enabled(tmp_path):
    validated = _validate(
        tmp_path,
        reproducibility={"enabled": True, "idr": {"atac_narrow": True}},
        idr={"threshold": "0.02"},
    )
    assert validated["idr"] == {"threshold": 0.02, "rank": "p.value", "seed": 42}


def test_idr_settings_ignored_when_no_idr_mode_enabled(tmp_path):
    validated = _validate(tmp_path, idr={"bad": "ignored"})
    assert validated["idr"] == {"threshold": 0.05, "rank": "p.value", "seed": 42}


def test_idr_settings_must_be_mapping_when_idr_mode_enabled(tmp_path):
    with pytest.raises(ValidationError, match="idr must be a mapping"):
        _validate(tmp_path, stage5=True, idr="bad")


@pytest.mark.parametrize("raw", [True, 0, 1, 1.5, "bad"])
def test_idr_threshold_rejects_invalid_values(tmp_path, raw):
    with pytest.raises(ValidationError, match="idr.threshold"):
        _validate(tmp_path, stage5=True, idr={"threshold": raw})


def test_idr_rank_rejects_invalid_value(tmp_path):
    with pytest.raises(ValidationError, match="idr.rank must be"):
        _validate(tmp_path, stage5=True, idr={"rank": "score"})


@pytest.mark.parametrize("raw", [True, 0, -1, "0", "-1", "abc", 1.5])
def test_idr_seed_rejects_invalid_values(tmp_path, raw):
    with pytest.raises(ValidationError, match="idr.seed must be"):
        _validate(tmp_path, stage5=True, idr={"seed": raw})


def test_idr_unknown_key_rejected_when_idr_mode_enabled(tmp_path):
    with pytest.raises(ValidationError, match="idr: unknown key 'bad'"):
        _validate(tmp_path, stage5=True, idr={"bad": True})


# ---------------------------------------------------------------------------
# Stage gates and warnings
# ---------------------------------------------------------------------------


def test_stage5_requires_stage4b(tmp_path):
    with pytest.raises(ValidationError, match="stage5=true requires stage4b=true"):
        _validate(tmp_path, stage4b=False, stage5=True)


@pytest.mark.parametrize(
    "flag,pattern",
    [
        ("atac_narrow", "reproducibility.idr.atac_narrow=true requires stage4b=true"),
        ("cuttag_narrow", "reproducibility.idr.cuttag_narrow=true requires stage4b=true"),
        (
            "chipseq_broad_experimental",
            "chipseq_broad_experimental=true or cuttag_broad_experimental=true",
        ),
        (
            "cuttag_broad_experimental",
            "chipseq_broad_experimental=true or cuttag_broad_experimental=true",
        ),
    ],
)
def test_reproducibility_idr_modes_require_stage4b(tmp_path, flag, pattern):
    broad_flags = {"chipseq_broad_experimental", "cuttag_broad_experimental"}
    warning_context = (
        pytest.warns(UserWarning, match="Experimental IDR flag enabled")
        if flag in broad_flags
        else nullcontext()
    )
    with warning_context:
        with pytest.raises(ValidationError, match=re.escape(pattern)):
            _validate(
                tmp_path,
                stage4b=False,
                reproducibility={"enabled": True, "idr": {flag: True}},
            )


@pytest.mark.parametrize(
    "flag", ["chipseq_broad_experimental", "cuttag_broad_experimental"]
)
def test_experimental_broad_idr_flags_emit_warning(tmp_path, flag):
    with pytest.warns(UserWarning, match="Experimental IDR flag enabled"):
        validated = _validate(
            tmp_path,
            reproducibility={"enabled": True, "idr": {flag: True}},
        )
    assert validated["reproducibility"]["idr"][flag] is True


def test_stage5_with_explicit_chipseq_narrow_false_warns(tmp_path):
    with pytest.warns(UserWarning, match="Config contradiction"):
        validated = _validate(
            tmp_path,
            stage5=True,
            reproducibility={
                "enabled": True,
                "idr": {"chipseq_narrow": False},
            },
        )
    assert validated["reproducibility"]["idr"]["chipseq_narrow"] is False
