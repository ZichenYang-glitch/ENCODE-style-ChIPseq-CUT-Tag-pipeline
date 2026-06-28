"""Direct-API characterization tests for QC config validation behavior.

These tests pin the current behavior of `_validate_qc_config` so that
PR62 can safely extract QC helpers without changing behavior.
"""

import pytest

from encode_pipeline.config.validate import ValidationError, validate_config


SAMPLES_HEADER = (
    "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\tpeak_mode\tgenome\tbowtie2_index\n"
)
SAMPLES_ROW = "S1\tR1.fq\tR2.fq\tPE\tchipseq\tT\tnarrow\ths\tidx\n"


QC_DEFAULTS = {
    "blacklist_filter": True,
    "frip": True,
    "library_complexity": True,
    "nrf_pbc": True,
    "signal_tracks": True,
    "summary": True,
    "cuttag_fragment_size": True,
    "cross_correlation": False,
    "preseq_complexity": False,
    "picard_metrics": False,
    "tss_enrichment": False,
}


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


def test_qc_omitted_expands_to_defaults(tmp_path):
    validated = validate_config(_make_config(tmp_path))
    assert validated["qc"] == QC_DEFAULTS


def test_qc_empty_mapping_expands_to_defaults(tmp_path):
    validated = validate_config(_make_config(tmp_path, qc={}))
    assert validated["qc"] == QC_DEFAULTS


@pytest.mark.parametrize("raw", ["yes", ["x"], True])
def test_qc_must_be_mapping(tmp_path, raw):
    with pytest.raises(ValidationError, match="qc must be a mapping"):
        validate_config(_make_config(tmp_path, qc=raw))


def test_qc_unknown_keys_are_silently_ignored(tmp_path):
    validated = validate_config(
        _make_config(tmp_path, qc={"unknown_qc_key": True, "frip": False})
    )
    assert "unknown_qc_key" not in validated["qc"]
    assert validated["qc"]["frip"] is False


# ---------------------------------------------------------------------------
# Known key normalization
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("key", list(QC_DEFAULTS))
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
def test_qc_known_keys_accept_bool_and_string_booleans(
    tmp_path, key, raw, expected
):
    validated = validate_config(_make_config(tmp_path, qc={key: raw}))
    assert validated["qc"][key] == expected


@pytest.mark.parametrize("key", list(QC_DEFAULTS))
def test_qc_known_keys_reject_invalid_strings(tmp_path, key):
    with pytest.raises(
        ValidationError, match=f"qc.{key} must be true or false"
    ):
        validate_config(_make_config(tmp_path, qc={key: "maybe"}))


# ---------------------------------------------------------------------------
# Explicit overrides
# ---------------------------------------------------------------------------


def test_qc_explicit_values_override_defaults(tmp_path):
    validated = validate_config(
        _make_config(
            tmp_path,
            qc={
                "blacklist_filter": "false",
                "cross_correlation": "true",
                "preseq_complexity": False,
            },
        )
    )
    assert validated["qc"]["blacklist_filter"] is False
    assert validated["qc"]["cross_correlation"] is True
    assert validated["qc"]["preseq_complexity"] is False
    assert validated["qc"]["frip"] is True
