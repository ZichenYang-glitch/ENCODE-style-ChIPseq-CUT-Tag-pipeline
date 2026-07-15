"""Direct-API characterization tests for CUT&Tag config validation behavior.

These tests pin the current behavior of `_validate_cuttag_config` so that
PR64 can safely extract CUT&Tag helpers without changing behavior.
"""

import pytest

from encode_pipeline.config.validate import ValidationError, validate_config


SAMPLES_HEADER = "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\tpeak_mode\tgenome\tbowtie2_index\n"
SAMPLES_ROW = "S1\tR1.fq\tR2.fq\tPE\tcuttag\tT\tnarrow\ths\tidx\n"


CUTTAG_DEFAULTS = {
    "peak_caller": "macs3",
    "seacr": {
        "enabled": False,
        "mode": "stringent",
        "normalization": "non",
        "threshold": 0.01,
    },
}


SEACR_DEFAULTS = CUTTAG_DEFAULTS["seacr"]


def _make_config(tmp_path, **overrides):
    """Return a minimal valid config dict pointing at a temporary sample TSV."""
    samples_path = tmp_path / "samples.tsv"
    samples_path.write_text(SAMPLES_HEADER + SAMPLES_ROW, encoding="utf-8")
    config = {"samples": str(samples_path), "use_control": False}
    config.update(overrides)
    return config


# ---------------------------------------------------------------------------
# Defaults and structure
# ---------------------------------------------------------------------------


def test_cuttag_omitted_expands_to_defaults(tmp_path):
    validated = validate_config(_make_config(tmp_path))
    assert validated["cuttag"] == CUTTAG_DEFAULTS


def test_cuttag_empty_mapping_expands_to_defaults(tmp_path):
    validated = validate_config(_make_config(tmp_path, cuttag={}))
    assert validated["cuttag"] == CUTTAG_DEFAULTS


@pytest.mark.parametrize("raw", ["yes", ["x"], True])
def test_cuttag_must_be_mapping(tmp_path, raw):
    with pytest.raises(ValidationError, match="cuttag must be a mapping"):
        validate_config(_make_config(tmp_path, cuttag=raw))


# ---------------------------------------------------------------------------
# Unknown key rejection
# ---------------------------------------------------------------------------


def test_cuttag_unknown_top_level_key_rejected(tmp_path):
    with pytest.raises(ValidationError, match="cuttag: unknown key"):
        validate_config(_make_config(tmp_path, cuttag={"unknown_key": True}))


def test_cuttag_unknown_seacr_key_rejected(tmp_path):
    with pytest.raises(ValidationError, match="cuttag.seacr: unknown key"):
        validate_config(_make_config(tmp_path, cuttag={"seacr": {"unknown_key": True}}))


# ---------------------------------------------------------------------------
# peak_caller
# ---------------------------------------------------------------------------


def test_cuttag_peak_caller_defaults_to_macs3(tmp_path):
    validated = validate_config(_make_config(tmp_path))
    assert validated["cuttag"]["peak_caller"] == "macs3"


def test_cuttag_peak_caller_accepts_macs3(tmp_path):
    validated = validate_config(_make_config(tmp_path, cuttag={"peak_caller": "macs3"}))
    assert validated["cuttag"]["peak_caller"] == "macs3"


@pytest.mark.parametrize("raw", ["seacr", "macs2", ""])
def test_cuttag_peak_caller_rejects_non_macs3(tmp_path, raw):
    with pytest.raises(ValidationError, match="cuttag.peak_caller must be 'macs3'"):
        validate_config(_make_config(tmp_path, cuttag={"peak_caller": raw}))


# ---------------------------------------------------------------------------
# seacr structure
# ---------------------------------------------------------------------------


def test_cuttag_seacr_empty_mapping_expands_to_defaults(tmp_path):
    validated = validate_config(_make_config(tmp_path, cuttag={"seacr": {}}))
    assert validated["cuttag"]["seacr"] == SEACR_DEFAULTS


@pytest.mark.parametrize("raw", ["yes", ["x"], True, False])
def test_cuttag_seacr_must_be_mapping(tmp_path, raw):
    with pytest.raises(ValidationError, match="cuttag.seacr must be a mapping"):
        validate_config(_make_config(tmp_path, cuttag={"seacr": raw}))


# ---------------------------------------------------------------------------
# seacr.enabled
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "raw,expected",
    [
        (True, True),
        (False, False),
        ("true", True),
        ("false", False),
        ("TRUE", True),
        ("FALSE", False),
    ],
)
def test_cuttag_seacr_enabled_accepts_bool_and_string_booleans(tmp_path, raw, expected):
    validated = validate_config(
        _make_config(tmp_path, cuttag={"seacr": {"enabled": raw}})
    )
    assert validated["cuttag"]["seacr"]["enabled"] is expected


def test_cuttag_seacr_enabled_rejects_invalid_strings(tmp_path):
    with pytest.raises(
        ValidationError, match="cuttag.seacr.enabled must be true or false"
    ):
        validate_config(_make_config(tmp_path, cuttag={"seacr": {"enabled": "maybe"}}))


# ---------------------------------------------------------------------------
# seacr.mode
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("raw", ["stringent", "relaxed"])
def test_cuttag_seacr_mode_accepts_allowed_values(tmp_path, raw):
    validated = validate_config(_make_config(tmp_path, cuttag={"seacr": {"mode": raw}}))
    assert validated["cuttag"]["seacr"]["mode"] == raw


@pytest.mark.parametrize("raw", ["non", "strict", "", "stringentish"])
def test_cuttag_seacr_mode_rejects_invalid_values(tmp_path, raw):
    with pytest.raises(
        ValidationError,
        match="cuttag.seacr.mode must be 'stringent' or 'relaxed'",
    ):
        validate_config(_make_config(tmp_path, cuttag={"seacr": {"mode": raw}}))


# ---------------------------------------------------------------------------
# seacr.normalization
# ---------------------------------------------------------------------------


def test_cuttag_seacr_normalization_defaults_to_non(tmp_path):
    validated = validate_config(_make_config(tmp_path))
    assert validated["cuttag"]["seacr"]["normalization"] == "non"


def test_cuttag_seacr_normalization_accepts_non(tmp_path):
    validated = validate_config(
        _make_config(tmp_path, cuttag={"seacr": {"normalization": "non"}})
    )
    assert validated["cuttag"]["seacr"]["normalization"] == "non"


@pytest.mark.parametrize("raw", ["log", "tmm", "cpm", ""])
def test_cuttag_seacr_normalization_rejects_non_non_values(tmp_path, raw):
    with pytest.raises(
        ValidationError,
        match="cuttag.seacr.normalization must be 'non'",
    ):
        validate_config(
            _make_config(tmp_path, cuttag={"seacr": {"normalization": raw}})
        )


# ---------------------------------------------------------------------------
# seacr.threshold
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("raw", [0.01, "0.05", 0.99])
def test_cuttag_seacr_threshold_accepts_floats_in_open_interval(tmp_path, raw):
    validated = validate_config(
        _make_config(tmp_path, cuttag={"seacr": {"threshold": raw}})
    )
    assert validated["cuttag"]["seacr"]["threshold"] == float(raw)


@pytest.mark.parametrize("raw", [0.0, 1.0, -0.01, 1.01])
def test_cuttag_seacr_threshold_rejects_out_of_range_values(tmp_path, raw):
    with pytest.raises(
        ValidationError, match="cuttag.seacr.threshold must be in \\(0, 1\\)"
    ):
        validate_config(_make_config(tmp_path, cuttag={"seacr": {"threshold": raw}}))


@pytest.mark.parametrize("raw", [True, False, "abc", None])
def test_cuttag_seacr_threshold_rejects_non_numeric_values(tmp_path, raw):
    with pytest.raises(
        ValidationError,
        match="cuttag.seacr.threshold must be a float in \\(0, 1\\)",
    ):
        validate_config(_make_config(tmp_path, cuttag={"seacr": {"threshold": raw}}))


# ---------------------------------------------------------------------------
# Explicit overrides
# ---------------------------------------------------------------------------


def test_cuttag_explicit_values_override_defaults(tmp_path):
    validated = validate_config(
        _make_config(
            tmp_path,
            cuttag={
                "peak_caller": "macs3",
                "seacr": {
                    "enabled": "true",
                    "mode": "relaxed",
                    "normalization": "non",
                    "threshold": "0.05",
                },
            },
        )
    )
    assert validated["cuttag"]["peak_caller"] == "macs3"
    assert validated["cuttag"]["seacr"]["enabled"] is True
    assert validated["cuttag"]["seacr"]["mode"] == "relaxed"
    assert validated["cuttag"]["seacr"]["normalization"] == "non"
    assert validated["cuttag"]["seacr"]["threshold"] == 0.05
