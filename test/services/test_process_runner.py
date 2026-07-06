"""Tests for the ProcessRunner subprocess boundary."""

import ast
import os
import sys
from pathlib import Path

import pytest

from encode_pipeline.platform.adapters import CommandSpec
from encode_pipeline.services.process_runner import ProcessResult, ProcessRunner


def _make_runner(**kwargs):
    """Return a ProcessRunner that allows the current Python executable."""
    exe = sys.executable
    name = Path(exe).name
    defaults = {"allowed_executables": (exe, name)}
    defaults.update(kwargs)
    return ProcessRunner(**defaults)


# ---------------------------------------------------------------------------
# Constructor validation
# ---------------------------------------------------------------------------


def test_process_runner_default_construction():
    runner = ProcessRunner()
    assert runner._allowed_executables == ("snakemake",)
    assert runner._timeout_seconds == 300.0
    assert runner._max_output_bytes == 10_000_000


def test_process_runner_rejects_empty_allowed_executables():
    with pytest.raises(ValueError, match="non-empty"):
        ProcessRunner(allowed_executables=())


def test_process_runner_rejects_empty_string_in_allowed_executables():
    with pytest.raises(ValueError, match="non-empty strings"):
        ProcessRunner(allowed_executables=("snakemake", ""))


def test_process_runner_rejects_non_tuple_allowed_executables():
    with pytest.raises(ValueError, match="must be a tuple"):
        ProcessRunner(allowed_executables=["snakemake"])


def test_process_runner_rejects_zero_timeout():
    with pytest.raises(ValueError, match="positive"):
        ProcessRunner(timeout_seconds=0)


def test_process_runner_rejects_negative_timeout():
    with pytest.raises(ValueError, match="positive"):
        ProcessRunner(timeout_seconds=-1.0)


def test_process_runner_rejects_bool_timeout():
    with pytest.raises(ValueError, match="must be a number"):
        ProcessRunner(timeout_seconds=True)


def test_process_runner_rejects_zero_max_output_bytes():
    with pytest.raises(ValueError, match="positive"):
        ProcessRunner(max_output_bytes=0)


def test_process_runner_rejects_negative_max_output_bytes():
    with pytest.raises(ValueError, match="positive"):
        ProcessRunner(max_output_bytes=-1)


def test_process_runner_rejects_bool_max_output_bytes():
    with pytest.raises(ValueError, match="must be an int"):
        ProcessRunner(max_output_bytes=True)


# ---------------------------------------------------------------------------
# Allowlist enforcement
# ---------------------------------------------------------------------------


def test_process_runner_rejects_executable_not_in_allowlist():
    runner = ProcessRunner(allowed_executables=("snakemake",))
    spec = CommandSpec(argv=("python", "-c", "print(1)"))
    result = runner.run(spec)
    assert result.is_failure is True
    issue = result.issues[0]
    assert issue.code == "PROCESS_RUNNER_EXECUTABLE_NOT_ALLOWED"
    assert issue.path == "command_spec.argv[0]"
    assert issue.source == "process_runner"


def test_process_runner_allows_executable_in_allowlist():
    runner = _make_runner()
    spec = CommandSpec(argv=(sys.executable, "-c", "print(42)"))
    result = runner.run(spec)
    assert result.is_success is True
    assert result.value.exit_code == 0


# ---------------------------------------------------------------------------
# Invalid spec
# ---------------------------------------------------------------------------


def test_process_runner_rejects_non_command_spec():
    runner = ProcessRunner()
    result = runner.run("not-a-command-spec")
    assert result.is_failure is True
    issue = result.issues[0]
    assert issue.code == "PROCESS_RUNNER_INVALID_COMMAND_SPEC"
    assert issue.path == "command_spec"
    assert issue.source == "process_runner"


# ---------------------------------------------------------------------------
# Successful execution
# ---------------------------------------------------------------------------


def test_process_runner_captures_stdout():
    runner = _make_runner()
    spec = CommandSpec(argv=(sys.executable, "-c", "print('hello')"))
    result = runner.run(spec)
    assert result.is_success is True
    assert result.value.exit_code == 0
    assert result.value.stdout.strip() == "hello"
    assert result.value.stderr == ""


def test_process_runner_captures_stderr():
    runner = _make_runner()
    spec = CommandSpec(
        argv=(sys.executable, "-c", "import sys; sys.stderr.write('err msg\\n')")
    )
    result = runner.run(spec)
    assert result.is_success is True
    assert result.value.exit_code == 0
    assert "err msg" in result.value.stderr


def test_process_runner_no_issues_on_clean_exit():
    runner = _make_runner()
    spec = CommandSpec(argv=(sys.executable, "-c", "pass"))
    result = runner.run(spec)
    assert result.is_success is True
    assert result.value.exit_code == 0
    assert len(result.issues) == 0


def test_process_runner_nonzero_exit_is_success_with_warning():
    runner = _make_runner()
    spec = CommandSpec(argv=(sys.executable, "-c", "import sys; sys.exit(3)"))
    result = runner.run(spec)
    assert result.is_success is True
    assert result.value.exit_code == 3
    nonzero_issues = [
        i for i in result.issues if i.code == "PROCESS_RUNNER_NONZERO_EXIT"
    ]
    assert len(nonzero_issues) == 1
    issue = nonzero_issues[0]
    assert issue.severity.value == "warning"
    assert issue.context == {"exit_code": 3}
    assert issue.source == "process_runner"


def test_process_runner_respects_cwd(tmp_path):
    runner = _make_runner()
    spec = CommandSpec(
        argv=(sys.executable, "-c", "import os; print(os.getcwd())"),
        cwd=str(tmp_path),
    )
    result = runner.run(spec)
    assert result.is_success is True
    assert result.value.stdout.strip() == str(tmp_path)


def test_process_runner_passes_env_overrides(monkeypatch):
    monkeypatch.delenv("PR112_TEST_VAR", raising=False)
    runner = _make_runner()
    spec = CommandSpec(
        argv=(
            sys.executable,
            "-c",
            "import os; print(os.environ.get('PR112_TEST_VAR', 'NOT_SET'))",
        ),
        env={"PR112_TEST_VAR": "hello_from_spec"},
    )
    result = runner.run(spec)
    assert result.is_success is True
    assert result.value.stdout.strip() == "hello_from_spec"
