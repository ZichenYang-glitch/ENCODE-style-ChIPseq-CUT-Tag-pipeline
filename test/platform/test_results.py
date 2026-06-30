"""Tests for workflow-platform Result / Issue primitives."""

import os
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

from encode_pipeline.errors import ValidationError
from encode_pipeline.platform.results import Issue, IssueSeverity, Result


SRC_ROOT = Path(__file__).resolve().parents[2] / "src"


def test_issue_defaults_and_to_dict():
    issue = Issue(code="INVALID_CONFIG", message="Config is invalid")

    assert issue.code == "INVALID_CONFIG"
    assert issue.message == "Config is invalid"
    assert issue.severity is IssueSeverity.ERROR
    assert issue.to_dict() == {
        "code": "INVALID_CONFIG",
        "message": "Config is invalid",
        "severity": "error",
        "path": None,
        "source": None,
        "technical_message": None,
        "hint": None,
        "context": {},
    }


@pytest.mark.parametrize(
    "kwargs",
    [
        {"code": "", "message": "Config is invalid"},
        {"code": "   ", "message": "Config is invalid"},
        {"code": "INVALID_CONFIG", "message": ""},
        {"code": "INVALID_CONFIG", "message": "   "},
    ],
)
def test_issue_rejects_empty_code_and_empty_message(kwargs):
    with pytest.raises(ValueError):
        Issue(**kwargs)


def test_issue_severity_serializes_as_lowercase_string():
    issue = Issue(
        code="CHECK_INPUT",
        message="Check input",
        severity=IssueSeverity.WARNING,
    )

    assert issue.to_dict()["severity"] == "warning"


def test_valid_string_severity_converts_to_issue_severity():
    issue = Issue(
        code="CHECK_INPUT",
        message="Check input",
        severity="warning",
    )

    assert issue.severity is IssueSeverity.WARNING


def test_invalid_severity_raises_value_error():
    with pytest.raises(ValueError):
        Issue(code="CHECK_INPUT", message="Check input", severity="fatal")


def test_issue_from_exception_uses_generic_exception_without_platform_importing_validator():
    issue = Issue.from_exception(
        ValidationError("config threads must be positive"),
        source="config",
        path="config.threads",
        hint="Use an integer greater than zero.",
        context={"field": "threads"},
    )

    assert issue == Issue(
        code="VALIDATION_ERROR",
        message="config threads must be positive",
        source="config",
        path="config.threads",
        technical_message="config threads must be positive",
        hint="Use an integer greater than zero.",
        context={"field": "threads"},
    )


def test_issue_context_is_copied_and_to_dict_returns_copy():
    context = {"sample": "S1"}
    issue = Issue(
        code="INVALID_SAMPLE",
        message="Sample is invalid",
        context=context,
    )

    context["sample"] = "S2"
    as_dict = issue.to_dict()
    as_dict["context"]["sample"] = "S3"

    assert issue.context == {"sample": "S1"}
    assert issue.to_dict()["context"] == {"sample": "S1"}


def test_source_is_preserved_as_component_level_string():
    issue = Issue(
        code="INVALID_SAMPLE",
        message="Sample is invalid",
        source="samples",
    )

    assert issue.source == "samples"
    assert issue.to_dict()["source"] == "samples"


def test_result_success_with_no_issues_is_success():
    result = Result.success({"ok": True})

    assert result.is_success
    assert not result.is_failure
    assert result.value == {"ok": True}
    assert result.issues == ()


def test_result_success_with_warning_issue_is_still_success_and_exposes_warnings():
    warning = Issue(
        code="CHECK_INPUT",
        message="Input may need review",
        severity=IssueSeverity.WARNING,
    )
    result = Result.success("value", issues=[warning])

    assert result.is_success
    assert not result.is_failure
    assert result.warnings == (warning,)
    assert result.errors == ()


def test_result_success_with_error_issue_raises_value_error():
    error = Issue(code="INVALID_CONFIG", message="Config is invalid")

    with pytest.raises(ValueError):
        Result.success("value", issues=[error])


def test_result_failure_with_error_issue_is_failure_and_exposes_errors():
    error = Issue(code="INVALID_CONFIG", message="Config is invalid")
    result = Result.failure([error])

    assert result.is_failure
    assert not result.is_success
    assert result.value is None
    assert result.errors == (error,)
    assert result.warnings == ()


def test_result_failure_with_warning_only_issues_raises_value_error():
    warning = Issue(
        code="CHECK_INPUT",
        message="Input may need review",
        severity=IssueSeverity.WARNING,
    )

    with pytest.raises(ValueError):
        Result.failure([warning])


def test_result_rejects_non_issue_entries():
    with pytest.raises(ValueError):
        Result.success("value", issues=["not an issue"])


def test_result_failure_rejects_non_issue_entries():
    with pytest.raises(ValueError):
        Result.failure(["not an issue"])


def test_result_to_dict_serializes_value_and_issues_with_serializer():
    warning = Issue(
        code="CHECK_INPUT",
        message="Input may need review",
        severity=IssueSeverity.WARNING,
    )
    result = Result.success(5, issues=[warning])

    assert result.to_dict(value_serializer=lambda value: {"number": value}) == {
        "ok": True,
        "value": {"number": 5},
        "issues": [warning.to_dict()],
    }


def test_result_generic_value_type_is_preserved():
    result = Result[int].success(5)

    assert isinstance(result.value, int)
    assert result.value == 5


def test_importing_platform_results_does_not_import_config_validator():
    code = """
        import sys
        import encode_pipeline.platform.results
        print('encode_pipeline.config.validator' in sys.modules)
    """
    env = dict(os.environ)
    env["PYTHONPATH"] = str(SRC_ROOT)
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    proc = subprocess.run(
        [sys.executable, "-c", textwrap.dedent(code)],
        capture_output=True,
        text=True,
        env=env,
    )

    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.strip() == "False"
