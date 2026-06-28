"""Unit tests for QC validation helpers."""

import pytest

from encode_pipeline.config.qc import validate_qc_config


class CustomError(Exception):
    """Sentinel exception for testing error injection."""


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


def test_validate_qc_config_expands_defaults():
    validated = validate_qc_config({})
    assert validated == QC_DEFAULTS


def test_validate_qc_config_normalizes_known_values():
    validated = validate_qc_config(
        {
            "frip": "false",
            "cross_correlation": "true",
            "preseq_complexity": False,
        }
    )
    assert validated["frip"] is False
    assert validated["cross_correlation"] is True
    assert validated["preseq_complexity"] is False
    assert validated["blacklist_filter"] is True


def test_validate_qc_config_silently_ignores_unknown_keys():
    validated = validate_qc_config({"unknown_qc_key": True, "frip": "false"})
    assert "unknown_qc_key" not in validated
    assert validated["frip"] is False


def test_validate_qc_config_rejects_invalid_values():
    with pytest.raises(ValueError, match="qc.frip must be true or false"):
        validate_qc_config({"frip": "maybe"})


def test_validate_qc_config_rejects_non_mapping():
    with pytest.raises(ValueError, match="qc must be a mapping"):
        validate_qc_config(["qc"])


def test_validate_qc_config_uses_custom_error_class():
    with pytest.raises(CustomError, match="qc.frip must be true or false"):
        validate_qc_config({"frip": "maybe"}, error_cls=CustomError)
