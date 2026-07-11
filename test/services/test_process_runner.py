"""Tests for the ProcessRunner subprocess boundary."""

import ast
import json
import sys
import time
from pathlib import Path

import pytest

import encode_pipeline.services.process_runner as process_runner_module
from encode_pipeline.platform.adapters import CommandSpec
from encode_pipeline.services.process_runner import ProcessRunner
from encode_pipeline.workers.timeouts import WorkerHardTimeout


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
    assert runner._passthrough_exceptions == ()


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


def test_process_runner_rejects_invalid_passthrough_exceptions():
    with pytest.raises(ValueError, match="must be a tuple"):
        ProcessRunner(passthrough_exceptions=[WorkerHardTimeout])
    with pytest.raises(ValueError, match="exception types"):
        ProcessRunner(passthrough_exceptions=("JobTimeoutException",))


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


def test_process_runner_scrubs_platform_agent_and_provider_secrets(monkeypatch):
    secret_names = (
        "ENCODE_PIPELINE_DATABASE_URL",
        "ENCODE_PIPELINE_REDIS_URL",
        "ENCODE_AGENT_LLM_API_KEY",
        "OPENAI_API_KEY",
        "DEEPSEEK_API_KEY",
        "GEMINI_CLI_IDE_AUTH_TOKEN",
        "DATABASE_URL",
        "LC_PRIVATE_VALUE",
    )
    for name in secret_names:
        monkeypatch.setenv(name, f"secret-for-{name}")
    monkeypatch.setenv("CONDA_PREFIX", "/safe/conda-prefix")
    runner = _make_runner()
    script = (
        "import json, os; "
        f"print(json.dumps({{name: os.environ.get(name) for name in {secret_names!r}}})); "
        "print(os.environ.get('CONDA_PREFIX'))"
    )

    result = runner.run(CommandSpec(argv=(sys.executable, "-c", script)))

    assert result.is_success
    lines = result.value.stdout.splitlines()
    assert json.loads(lines[0]) == {name: None for name in secret_names}
    assert lines[1] == "/safe/conda-prefix"


def test_process_runner_inherits_only_standard_locale_variables(monkeypatch):
    monkeypatch.setenv("LC_ALL", "C")
    monkeypatch.setenv("LC_PRIVATE_VALUE", "private-locale-value")
    runner = _make_runner()
    script = (
        "import json, os; "
        "print(json.dumps({name: os.environ.get(name) "
        "for name in ('LC_ALL', 'LC_PRIVATE_VALUE')}))"
    )

    result = runner.run(CommandSpec(argv=(sys.executable, "-c", script)))

    assert result.is_success
    assert json.loads(result.value.stdout) == {
        "LC_ALL": "C",
        "LC_PRIVATE_VALUE": None,
    }


def test_process_runner_rejects_protected_command_environment():
    runner = _make_runner()
    result = runner.run(
        CommandSpec(
            argv=(sys.executable, "-c", "pass"),
            env={"ENCODE_PIPELINE_REDIS_URL": "redis://secret@localhost/0"},
        )
    )

    assert result.is_failure
    assert result.issues[0].code == "PROCESS_RUNNER_PROTECTED_ENVIRONMENT"
    assert "secret" not in str(result.issues[0].to_dict())


@pytest.mark.parametrize(
    "protected_name",
    (
        "DATABASE_URL",
        "MY_REDIS_URL",
        "MY_SECRET_FILE",
        "SERVICE_PRIVATE_KEY_PATH",
    ),
)
def test_process_runner_rejects_broader_sensitive_environment_names(
    protected_name,
):
    result = _make_runner().run(
        CommandSpec(
            argv=(sys.executable, "-c", "pass"),
            env={protected_name: "private-value"},
        )
    )

    assert result.is_failure
    assert result.issues[0].code == "PROCESS_RUNNER_PROTECTED_ENVIRONMENT"
    assert "private-value" not in str(result.issues[0].to_dict())


# ---------------------------------------------------------------------------
# Failure modes
# ---------------------------------------------------------------------------


def test_process_runner_timeout_returns_failure():
    runner = ProcessRunner(
        allowed_executables=(sys.executable, Path(sys.executable).name),
        timeout_seconds=0.1,
    )
    spec = CommandSpec(argv=(sys.executable, "-c", "import time; time.sleep(10)"))
    result = runner.run(spec)
    assert result.is_failure is True
    issue = result.issues[0]
    assert issue.code == "PROCESS_RUNNER_TIMEOUT"
    assert issue.severity.value == "error"
    assert issue.context == {"timeout_seconds": 0.1}


def test_process_runner_timeout_terminates_direct_child(monkeypatch):
    terminated = []
    original_terminate = ProcessRunner._terminate

    def track_termination(process):
        terminated.append(process)
        original_terminate(process)

    monkeypatch.setattr(
        ProcessRunner,
        "_terminate",
        staticmethod(track_termination),
    )
    runner = _make_runner(timeout_seconds=0.1)

    result = runner.run(
        CommandSpec(argv=(sys.executable, "-c", "import time; time.sleep(10)"))
    )

    assert result.is_failure
    assert result.issues[0].code == "PROCESS_RUNNER_TIMEOUT"
    assert len(terminated) == 1
    assert terminated[0].poll() is not None


def test_process_runner_timeout_reports_output_truncation():
    chunks: list[str] = []
    result = _make_runner(timeout_seconds=0.1, max_output_bytes=8).run(
        CommandSpec(
            argv=(
                sys.executable,
                "-c",
                "import time; print('x' * 1000, flush=True); time.sleep(10)",
            )
        ),
        output_callback=lambda _stream, lines: chunks.extend(lines),
    )

    assert result.is_failure
    assert [issue.code for issue in result.issues] == [
        "PROCESS_RUNNER_TIMEOUT",
        "PROCESS_RUNNER_OUTPUT_LIMIT_EXCEEDED",
    ]
    assert chunks == ["xxxxxxxx"]


def test_process_runner_timeout_drains_bytes_already_buffered_at_kill(monkeypatch):
    real_selector_factory = process_runner_module.selectors.DefaultSelector

    class HideReadsUntilNonBlockingDrain:
        def __init__(self):
            self._selector = real_selector_factory()

        def register(self, *args, **kwargs):
            return self._selector.register(*args, **kwargs)

        def unregister(self, *args, **kwargs):
            return self._selector.unregister(*args, **kwargs)

        def get_map(self):
            return self._selector.get_map()

        def select(self, timeout=None):
            if timeout != 0:
                time.sleep(timeout or 0)
                return []
            return self._selector.select(timeout=0)

        def close(self):
            self._selector.close()

    monkeypatch.setattr(
        process_runner_module.selectors,
        "DefaultSelector",
        HideReadsUntilNonBlockingDrain,
    )
    chunks: list[str] = []

    result = _make_runner(timeout_seconds=0.1).run(
        CommandSpec(
            argv=(
                sys.executable,
                "-c",
                "import time; print('last-before-timeout', flush=True); time.sleep(10)",
            )
        ),
        output_callback=lambda _stream, lines: chunks.extend(lines),
    )

    assert result.is_failure
    assert result.issues[0].code == "PROCESS_RUNNER_TIMEOUT"
    assert chunks == ["last-before-timeout"]


def test_post_timeout_drain_is_bounded_for_continuously_ready_descendant(
    monkeypatch,
):
    real_selector_factory = process_runner_module.selectors.DefaultSelector
    chunks: list[str] = []

    class ContinuouslyReadyAfterTimeout:
        def __init__(self):
            self._selector = real_selector_factory()

        def register(self, *args, **kwargs):
            return self._selector.register(*args, **kwargs)

        def unregister(self, *args, **kwargs):
            return self._selector.unregister(*args, **kwargs)

        def get_map(self):
            return self._selector.get_map()

        def select(self, timeout=None):
            if timeout != 0:
                time.sleep(timeout or 0)
                return []
            return self._selector.select(timeout=0)

        def close(self):
            self._selector.close()

    monkeypatch.setattr(
        process_runner_module.selectors,
        "DefaultSelector",
        ContinuouslyReadyAfterTimeout,
    )
    grandchild = "import os; data=b'x'*65536\nwhile True: os.write(1, data)"
    parent = (
        "import subprocess, sys, time; "
        f"subprocess.Popen([sys.executable, '-c', {grandchild!r}]); "
        "time.sleep(10)"
    )
    started = time.monotonic()

    result = _make_runner(timeout_seconds=0.1).run(
        CommandSpec(argv=(sys.executable, "-c", parent)),
        output_callback=lambda _stream, lines: chunks.extend(lines),
    )
    elapsed = time.monotonic() - started

    assert result.is_failure
    assert result.issues[0].code == "PROCESS_RUNNER_TIMEOUT"
    assert (
        0
        < len("".join(chunks).encode("utf-8"))
        <= (process_runner_module._POST_TERMINATION_DRAIN_MAX_BYTES)
    )
    assert elapsed < 1.0


def test_process_runner_executable_not_found():
    runner = ProcessRunner(allowed_executables=("nonexistent_binary_xyz_123",))
    spec = CommandSpec(argv=("nonexistent_binary_xyz_123", "--help"))
    result = runner.run(spec)
    assert result.is_failure is True
    issue = result.issues[0]
    assert issue.code == "PROCESS_RUNNER_EXECUTABLE_NOT_FOUND"
    assert issue.path == "command_spec.argv[0]"
    assert issue.source == "process_runner"


def test_process_runner_output_truncation():
    runner = ProcessRunner(
        allowed_executables=(sys.executable, Path(sys.executable).name),
        max_output_bytes=100,
    )
    spec = CommandSpec(argv=(sys.executable, "-c", "print('x' * 200_000)"))
    result = runner.run(spec)
    assert result.is_success is True
    assert result.value.exit_code == 0
    assert len(result.value.stdout) <= 100
    trunc_issues = [
        i for i in result.issues if i.code == "PROCESS_RUNNER_OUTPUT_LIMIT_EXCEEDED"
    ]
    assert len(trunc_issues) == 1
    assert trunc_issues[0].severity.value == "warning"


def test_process_runner_stderr_truncation():
    runner = ProcessRunner(
        allowed_executables=(sys.executable, Path(sys.executable).name),
        max_output_bytes=50,
    )
    spec = CommandSpec(
        argv=(
            sys.executable,
            "-c",
            "import sys; sys.stderr.write('e' * 100_000)",
        )
    )
    result = runner.run(spec)
    assert result.is_success is True
    assert len(result.value.stderr) <= 50
    trunc_issues = [
        i for i in result.issues if i.code == "PROCESS_RUNNER_OUTPUT_LIMIT_EXCEEDED"
    ]
    assert len(trunc_issues) == 1


def test_streamed_output_is_bounded_by_capture_limit():
    runner = _make_runner(max_output_bytes=64)
    chunks: list[tuple[str, tuple[str, ...]]] = []
    spec = CommandSpec(argv=(sys.executable, "-c", "print('x' * 10000)"))

    result = runner.run(
        spec, output_callback=lambda stream, lines: chunks.append((stream, lines))
    )

    assert result.is_success
    persisted = "".join(
        line for stream, lines in chunks if stream == "stdout" for line in lines
    )
    assert len(persisted.encode("utf-8")) <= 64
    assert len(result.value.stdout.encode("utf-8")) <= 64
    assert any(
        issue.code == "PROCESS_RUNNER_OUTPUT_LIMIT_EXCEEDED" for issue in result.issues
    )


def test_streaming_preserves_utf8_split_across_pipe_reads():
    runner = _make_runner()
    chunks: list[str] = []
    script = (
        "import sys, time; data='你好'.encode(); "
        "sys.stdout.buffer.write(data[:2]); sys.stdout.flush(); time.sleep(.05); "
        "sys.stdout.buffer.write(data[2:] + b'\\n'); sys.stdout.flush()"
    )

    result = runner.run(
        CommandSpec(argv=(sys.executable, "-c", script)),
        output_callback=lambda _stream, lines: chunks.extend(lines),
    )

    assert result.is_success
    assert chunks == ["你好"]


def test_streaming_preserves_long_utf8_line_across_read_boundaries():
    runner = _make_runner(max_output_bytes=200_000)
    chunks: list[str] = []
    expected = "你" * 30_000

    result = runner.run(
        CommandSpec(
            argv=(
                sys.executable,
                "-c",
                "import sys; sys.stdout.write('你' * 30000); sys.stdout.flush()",
            )
        ),
        output_callback=lambda _stream, lines: chunks.extend(lines),
    )

    assert result.is_success
    assert "".join(chunks) == expected
    assert result.value.stdout == expected
    assert "\ufffd" not in "".join(chunks)


def test_streaming_callback_runs_before_process_exit():
    runner = _make_runner()
    callback_times: list[float] = []
    started = time.monotonic()
    script = "import time; print('ready', flush=True); time.sleep(.2)"

    result = runner.run(
        CommandSpec(argv=(sys.executable, "-c", script)),
        output_callback=lambda _stream, _lines: callback_times.append(time.monotonic()),
    )
    ended = time.monotonic()

    assert result.is_success
    assert callback_times
    assert callback_times[0] - started < ended - callback_times[0]


def test_streaming_callback_failure_is_public_safe():
    runner = _make_runner()

    def fail_callback(_stream, _lines):
        raise RuntimeError("private database detail")

    result = runner.run(
        CommandSpec(argv=(sys.executable, "-c", "print('output', flush=True)")),
        output_callback=fail_callback,
    )

    assert result.is_failure
    assert result.issues[0].code == "PROCESS_RUNNER_OUTPUT_CALLBACK_FAILED"
    assert "private database detail" not in str(result.issues[0].to_dict())


def test_streaming_reraises_injected_worker_hard_timeout():
    runner = _make_runner(passthrough_exceptions=(WorkerHardTimeout,))

    def timeout_callback(_stream, _lines):
        raise WorkerHardTimeout("RQ deadline reached")

    with pytest.raises(WorkerHardTimeout, match="RQ deadline reached"):
        runner.run(
            CommandSpec(
                argv=(
                    sys.executable,
                    "-c",
                    "import time; print('ready', flush=True); time.sleep(10)",
                )
            ),
            output_callback=timeout_callback,
        )


def test_base_exception_unwind_terminates_direct_child(monkeypatch):
    class AbortExecution(BaseException):
        pass

    terminated = []
    original_terminate = ProcessRunner._terminate

    def track_termination(process):
        terminated.append(process)
        original_terminate(process)

    monkeypatch.setattr(
        ProcessRunner,
        "_terminate",
        staticmethod(track_termination),
    )

    def abort_callback(_stream, _lines):
        raise AbortExecution

    with pytest.raises(AbortExecution):
        _make_runner().run(
            CommandSpec(
                argv=(
                    sys.executable,
                    "-c",
                    "import time; print('ready', flush=True); time.sleep(10)",
                )
            ),
            output_callback=abort_callback,
        )

    assert len(terminated) == 1
    assert terminated[0].poll() is not None


# ---------------------------------------------------------------------------
# Safety: no env or path leakage in issues
# ---------------------------------------------------------------------------


def test_process_runner_issues_contain_no_env_values(monkeypatch):
    monkeypatch.setenv("PR112_SECRET", "do-not-leak")
    runner = _make_runner()
    spec = CommandSpec(
        argv=(sys.executable, "-c", "import sys; sys.exit(1)"),
    )
    result = runner.run(spec)
    for issue in result.issues:
        assert "do-not-leak" not in issue.message
        assert "do-not-leak" not in str(issue.context)
        assert "PR112_SECRET" not in issue.message


def test_process_runner_issues_contain_no_filesystem_paths(tmp_path):
    runner = _make_runner()
    spec = CommandSpec(
        argv=(sys.executable, "-c", "import sys; sys.exit(1)"),
        cwd=str(tmp_path),
    )
    result = runner.run(spec)
    cwd_str = str(tmp_path)
    for issue in result.issues:
        assert cwd_str not in issue.message
        assert cwd_str not in str(issue.context or {})


# ---------------------------------------------------------------------------
# Import boundary
# ---------------------------------------------------------------------------


def test_process_runner_import_boundary():
    source_path = (
        Path(__file__).resolve().parents[2]
        / "src/encode_pipeline/services/process_runner.py"
    )
    source = source_path.read_text(encoding="utf-8")
    tree = ast.parse(source)

    forbidden_modules = {
        "encode_pipeline.api",
        "encode_pipeline.frontend",
        "fastapi",
        "pydantic",
        "snakemake",
        "encode_pipeline.services.runs",
        "encode_pipeline.services.local_run_driver",
    }
    imported_modules: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported_modules.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module is not None:
            imported_modules.add(node.module)

    violations = [
        m
        for m in imported_modules
        for f in forbidden_modules
        if m == f or m.startswith(f"{f}.")
    ]
    assert not violations, f"Forbidden imports found: {violations}"


def test_process_runner_allows_required_stdlib_imports():
    source_path = (
        Path(__file__).resolve().parents[2]
        / "src/encode_pipeline/services/process_runner.py"
    )
    source = source_path.read_text(encoding="utf-8")
    tree = ast.parse(source)

    imported_modules: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported_modules.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module is not None:
            imported_modules.add(node.module)

    assert "subprocess" in imported_modules
    assert "os" in imported_modules
