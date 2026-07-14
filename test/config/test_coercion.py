"""Direct-API characterization tests for config coercion behavior.

These tests pin the current behavior of
encode_pipeline.config.validate.validate_config so that future extractions
(e.g., PR54 coercion helpers) can be verified behavior-preserving.
"""

import pytest

from encode_pipeline.config.validate import ValidationError, validate_config


SAMPLES_HEADER = (
    "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\tpeak_mode\tgenome\tbowtie2_index\n"
)
SAMPLES_ROW = "S1\tR1.fq\tR2.fq\tPE\tchipseq\tT\tnarrow\ths\tidx\n"


def _make_config(tmp_path, **overrides):
    """Return a minimal valid config dict pointing at a temporary sample TSV."""
    samples_path = tmp_path / "samples.tsv"
    samples_path.write_text(SAMPLES_HEADER + SAMPLES_ROW, encoding="utf-8")
    config = {"samples": str(samples_path)}
    config.update(overrides)
    return config


# ---------------------------------------------------------------------------
# threads
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("raw,expected", [(4, 4), ("8", 8)])
def test_threads_accepts_positive_int_and_numeric_string(tmp_path, raw, expected):
    config = _make_config(tmp_path, threads=raw)
    validated = validate_config(config)
    assert validated["threads"] == expected


@pytest.mark.parametrize("raw", [True, False, 0, -1, "-5", "3.5", "abc"])
def test_threads_rejects_invalid_values(tmp_path, raw):
    config = _make_config(tmp_path, threads=raw)
    with pytest.raises(ValidationError, match="must be an integer|must be positive"):
        validate_config(config)


# ---------------------------------------------------------------------------
# mapq
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("raw,expected", [(0, 0), (30, 30), ("20", 20)])
def test_mapq_accepts_non_negative_int_and_numeric_string(tmp_path, raw, expected):
    config = _make_config(tmp_path, mapq=raw)
    validated = validate_config(config)
    assert validated["mapq"] == expected


@pytest.mark.parametrize("raw", [True, False, -1, "-5", "3.5"])
def test_mapq_rejects_invalid_values(tmp_path, raw):
    config = _make_config(tmp_path, mapq=raw)
    with pytest.raises(
        ValidationError, match="must be an integer|must be non-negative"
    ):
        validate_config(config)


# ---------------------------------------------------------------------------
# binsize
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("raw,expected", [(10, 10), ("50", 50)])
def test_binsize_accepts_positive_int_and_numeric_string(tmp_path, raw, expected):
    config = _make_config(tmp_path, binsize=raw)
    validated = validate_config(config)
    assert validated["binsize"] == expected


@pytest.mark.parametrize("raw", [True, False, 0, -1, "-5", "3.5", "abc"])
def test_binsize_rejects_invalid_values(tmp_path, raw):
    config = _make_config(tmp_path, binsize=raw)
    with pytest.raises(ValidationError, match="must be an integer|must be positive"):
        validate_config(config)


# ---------------------------------------------------------------------------
# trim
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "raw,expected",
    [
        (True, "true"),
        (False, "false"),
        ("true", "true"),
        ("false", "false"),
        ("TRUE", "true"),
    ],
)
def test_trim_accepts_bool_and_string_booleans(tmp_path, raw, expected):
    config = _make_config(tmp_path, trim=raw)
    validated = validate_config(config)
    assert validated["trim"] == expected


def test_trim_rejects_invalid_string(tmp_path):
    config = _make_config(tmp_path, trim="yes")
    with pytest.raises(ValidationError, match="must be true or false"):
        validate_config(config)


# ---------------------------------------------------------------------------
# use_control
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "raw,expected",
    [
        (True, True),
        (False, False),
        ("true", True),
        ("yes", True),
        ("1", True),
        ("false", False),
        ("no", False),
        ("0", False),
    ],
)
def test_use_control_accepts_bool_and_synonyms(tmp_path, raw, expected):
    config = _make_config(tmp_path, use_control=raw)
    validated = validate_config(config)
    assert validated["use_control"] == expected


def test_use_control_rejects_invalid_string(tmp_path):
    config = _make_config(tmp_path, use_control="maybe")
    with pytest.raises(ValidationError, match="must be true or false"):
        validate_config(config)


# ---------------------------------------------------------------------------
# multiqc / stage4b / stage5
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("field", ["multiqc", "stage4b", "stage5"])
@pytest.mark.parametrize(
    "raw,expected",
    [
        (True, True),
        (False, False),
        ("true", True),
        ("false", False),
        ("TRUE", True),
    ],
)
def test_bool_flag_accepts_bool_and_string_booleans(tmp_path, field, raw, expected):
    config = _make_config(tmp_path, **{field: raw})
    validated = validate_config(config)
    assert validated[field] == expected


@pytest.mark.parametrize("field", ["multiqc", "stage4b", "stage5"])
def test_bool_flag_rejects_invalid_string(tmp_path, field):
    config = _make_config(tmp_path, **{field: "yes"})
    with pytest.raises(ValidationError, match="must be true or false"):
        validate_config(config)


# ---------------------------------------------------------------------------
# reproducibility boolean coercion (_coerce_bool call sites)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "raw,expected",
    [
        (True, True),
        (False, False),
        ("true", True),
        ("false", False),
        ("TRUE", True),
    ],
)
def test_reproducibility_enabled_accepts_bool_and_string_booleans(
    tmp_path, raw, expected
):
    config = _make_config(tmp_path, reproducibility={"enabled": raw})
    validated = validate_config(config)
    assert validated["reproducibility"]["enabled"] == expected


def test_reproducibility_enabled_rejects_invalid_string(tmp_path):
    config = _make_config(tmp_path, reproducibility={"enabled": "yes"})
    with pytest.raises(
        ValidationError, match="reproducibility.enabled must be true or false"
    ):
        validate_config(config)


def test_reproducibility_consensus_enabled_accepts_string_false(tmp_path):
    config = _make_config(
        tmp_path,
        reproducibility={"enabled": True, "consensus": {"enabled": "false"}},
    )
    validated = validate_config(config)
    assert validated["reproducibility"]["consensus"]["enabled"] is False


@pytest.mark.parametrize(
    "flag,raw,expected",
    [
        ("atac_narrow", "true", True),
        ("cuttag_narrow", "false", False),
        ("chipseq_broad_experimental", "true", True),
        ("cuttag_broad_experimental", "false", False),
    ],
)
def test_reproducibility_idr_flag_accepts_string_booleans(
    tmp_path, flag, raw, expected
):
    config = _make_config(
        tmp_path,
        reproducibility={"enabled": True, "idr": {flag: raw}},
    )
    validated = validate_config(config)
    assert validated["reproducibility"]["idr"][flag] == expected


@pytest.mark.parametrize(
    "flag",
    [
        "atac_narrow",
        "cuttag_narrow",
        "chipseq_broad_experimental",
        "cuttag_broad_experimental",
    ],
)
def test_reproducibility_idr_flag_rejects_invalid_string(tmp_path, flag):
    config = _make_config(
        tmp_path,
        reproducibility={"enabled": True, "idr": {flag: "yes"}},
    )
    with pytest.raises(
        ValidationError,
        match=f"reproducibility.idr.{flag} must be true or false",
    ):
        validate_config(config)
