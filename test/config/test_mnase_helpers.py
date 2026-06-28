"""Unit tests for MNase validation helpers."""

import pytest

from encode_pipeline.config.mnase import validate_mnase_config


class CustomError(Exception):
    """Sentinel exception for testing error injection."""


MNASE_DEFAULTS = {
    "mono_range": [140, 200],
    "fragments": {
        "sub": [1, 139],
        "mono": [140, 200],
        "di": [300, 400],
    },
    "dyad_range": [130, 200],
    "callers": {"danpos3": False, "inps": False, "sem": False},
}


def test_validate_mnase_config_expands_defaults():
    validated = validate_mnase_config({})
    assert validated == MNASE_DEFAULTS


def test_validate_mnase_config_tuple_and_string_endpoints_normalize():
    validated = validate_mnase_config(
        {
            "mono_range": ("140", "200"),
            "fragments": {"sub": (10, 80)},
            "dyad_range": ["130", "200"],
        }
    )
    assert validated["mono_range"] == [140, 200]
    assert validated["fragments"]["sub"] == [10, 80]
    assert validated["dyad_range"] == [130, 200]
    assert isinstance(validated["mono_range"], list)
    assert isinstance(validated["fragments"]["sub"], list)
    assert isinstance(validated["dyad_range"], list)


def test_validate_mnase_config_fragments_mono_precedence():
    validated = validate_mnase_config(
        {
            "mono_range": [100, 180],
            "fragments": {"mono": [140, 200]},
        }
    )
    assert validated["mono_range"] == [100, 180]
    assert validated["fragments"]["mono"] == [140, 200]


def test_validate_mnase_config_invalid_range_raises():
    with pytest.raises(ValueError, match="mnase.mono_range.*min must be < max"):
        validate_mnase_config({"mono_range": [200, 100]})


def test_validate_mnase_config_unknown_top_level_key_rejected():
    with pytest.raises(ValueError, match="mnase: unknown key"):
        validate_mnase_config({"unknown": True})


def test_validate_mnase_config_callers_true_rejected_as_not_implemented():
    with pytest.raises(ValueError, match="caller execution is not implemented"):
        validate_mnase_config({"callers": {"danpos3": True}})


def test_validate_mnase_config_uses_custom_error_class():
    with pytest.raises(CustomError, match="mnase: unknown key"):
        validate_mnase_config({"unknown": True}, error_cls=CustomError)
