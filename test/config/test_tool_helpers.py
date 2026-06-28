"""Unit tests for tool parameter validation helpers."""

import pytest

from encode_pipeline.config.tools import validate_tool_params


class CustomError(Exception):
    """Sentinel exception for testing error injection."""


def test_validate_tool_params_expands_defaults():
    validated = validate_tool_params({})
    assert validated["fastqc"]["extra_args"] == ""
    assert validated["bamcoverage"]["normalize_using"] == "CPM"
    assert validated["macs3"]["qvalue"] == 0.01


def test_validate_tool_params_normalizes_values():
    validated = validate_tool_params(
        {
            "bowtie2": {"dovetail": "true"},
            "samtools_filter": {"filter_flags": "0x904"},
            "macs3": {"qvalue": "0.05"},
        }
    )
    assert validated["bowtie2"]["dovetail"] is True
    assert validated["samtools_filter"]["filter_flags"] == 2308
    assert validated["macs3"]["qvalue"] == 0.05


def test_validate_tool_params_uses_custom_error_class():
    with pytest.raises(CustomError, match="unknown tool block"):
        validate_tool_params({"bad_tool": {}}, error_cls=CustomError)


def test_validate_tool_params_rejects_invalid_nested_values():
    with pytest.raises(ValueError, match="tool_parameters.bowtie2.mode"):
        validate_tool_params({"bowtie2": {"mode": "ultra"}})
