"""Unit tests for reproducibility/IDR validation helpers."""

import pytest

from encode_pipeline.config.reproducibility import (
    validate_idr_settings,
    validate_reproducibility,
)


class CustomError(Exception):
    """Sentinel exception for testing error injection."""


def test_validate_reproducibility_expands_defaults():
    validated = validate_reproducibility({"enabled": True}, {"stage5": False})
    assert validated == {
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


def test_validate_reproducibility_uses_custom_error_class():
    with pytest.raises(CustomError, match="unknown key 'bad'"):
        validate_reproducibility(
            {"enabled": True, "idr": {"bad": True}},
            {"stage5": False},
            error_cls=CustomError,
        )


def test_validate_idr_settings_normalizes_values():
    validated = validate_idr_settings(
        {"threshold": "0.01", "rank": "signal.value", "seed": "123"}
    )
    assert validated == {
        "threshold": 0.01,
        "rank": "signal.value",
        "seed": 123,
    }


def test_validate_idr_settings_uses_custom_error_class():
    with pytest.raises(CustomError, match="idr.rank must be"):
        validate_idr_settings({"rank": "score"}, error_cls=CustomError)
