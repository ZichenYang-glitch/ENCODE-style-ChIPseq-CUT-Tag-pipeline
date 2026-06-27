"""Unit tests for primitive config coercion helpers."""

import pytest

from encode_pipeline.config.coercion import coerce_bool, coerce_int


class CustomError(Exception):
    """Sentinel exception for testing error injection."""


@pytest.mark.parametrize("raw,expected", [(4, 4), ("8", 8)])
def test_coerce_int_accepts_int_and_numeric_string(raw, expected):
    assert coerce_int(raw, name="threads", minimum=1) == expected


@pytest.mark.parametrize("raw", [True, False, "3.5", "abc"])
def test_coerce_int_rejects_non_integer_values(raw):
    with pytest.raises(ValueError, match="config threads must be an integer"):
        coerce_int(raw, name="threads", minimum=1)


def test_coerce_int_uses_custom_error_class():
    with pytest.raises(CustomError, match="config threads must be positive"):
        coerce_int(0, name="threads", minimum=1, error_cls=CustomError)


@pytest.mark.parametrize(
    "raw,expected",
    [(True, True), (False, False), ("true", True), ("false", False)],
)
def test_coerce_bool_accepts_bool_and_string_bool(raw, expected):
    assert coerce_bool(raw, "feature.enabled") is expected


def test_coerce_bool_rejects_invalid_values():
    with pytest.raises(ValueError, match="feature.enabled must be true or false"):
        coerce_bool("yes", "feature.enabled")


def test_coerce_bool_uses_custom_error_class():
    with pytest.raises(CustomError, match="feature.enabled must be true or false"):
        coerce_bool("yes", "feature.enabled", error_cls=CustomError)
