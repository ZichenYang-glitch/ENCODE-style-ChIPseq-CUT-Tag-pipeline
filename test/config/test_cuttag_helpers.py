"""Unit tests for CUT&Tag validation helpers."""

import pytest

from encode_pipeline.config.cuttag import validate_cuttag_config


class CustomError(Exception):
    """Sentinel exception for testing error injection."""


CUTTAG_DEFAULTS = {
    "peak_caller": "macs3",
    "seacr": {
        "enabled": False,
        "mode": "stringent",
        "normalization": "non",
        "threshold": 0.01,
    },
}


def test_validate_cuttag_config_expands_defaults():
    validated = validate_cuttag_config({})
    assert validated == CUTTAG_DEFAULTS


def test_validate_cuttag_config_explicit_seacr_values_normalize():
    validated = validate_cuttag_config(
        {
            "peak_caller": "macs3",
            "seacr": {
                "enabled": "true",
                "mode": "relaxed",
                "normalization": "non",
                "threshold": "0.05",
            },
        }
    )
    assert validated["peak_caller"] == "macs3"
    assert validated["seacr"]["enabled"] is True
    assert validated["seacr"]["mode"] == "relaxed"
    assert validated["seacr"]["normalization"] == "non"
    assert validated["seacr"]["threshold"] == 0.05


def test_validate_cuttag_config_rejects_unknown_top_level_keys():
    with pytest.raises(ValueError, match="cuttag: unknown key"):
        validate_cuttag_config({"unknown_key": True})


def test_validate_cuttag_config_rejects_unknown_seacr_keys():
    with pytest.raises(ValueError, match="cuttag.seacr: unknown key"):
        validate_cuttag_config({"seacr": {"unknown_key": True}})


def test_validate_cuttag_config_rejects_bool_threshold():
    with pytest.raises(
        ValueError, match="cuttag.seacr.threshold must be a float in \\(0, 1\\)"
    ):
        validate_cuttag_config({"seacr": {"threshold": True}})


def test_validate_cuttag_config_uses_custom_error_class():
    with pytest.raises(CustomError, match="cuttag: unknown key"):
        validate_cuttag_config({"unknown_key": True}, error_cls=CustomError)
