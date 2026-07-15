"""Tests for the ProcessRunner subprocess boundary."""

import ast
import ctypes
import json
import os
import sys
import time
import tracemalloc
from pathlib import Path

import pytest

import encode_pipeline.services.process_runner as process_runner_module
from encode_pipeline.platform.adapters import CommandSpec
from encode_pipeline.platform.managed_containers import (
    managed_container_endpoint_identity,
)
from encode_pipeline.platform.results import Issue, Result
from encode_pipeline.services.managed_containers import ManagedContainerCleaner
from encode_pipeline.services.process_runner import (
    ProcessRunner,
    ProcessRunnerCleanupError,
)
from encode_pipeline.workers.timeouts import WorkerHardTimeout


_MANAGED_ENDPOINT_IDENTITY = managed_container_endpoint_identity(
    Path("/server/bin/docker"),
    Path("/server/run/docker.sock"),
)


def _make_runner(**kwargs):
    """Return a ProcessRunner that allows the current Python executable."""
    exe = sys.executable
    name = Path(exe).name
    defaults = {"allowed_executables": (exe, name)}
    defaults.update(kwargs)
    return ProcessRunner(**defaults)


class _RecordingContainerCleaner(ManagedContainerCleaner):
    def __init__(self, *, fail: bool = False, endpoint_changed: bool = False):
        object.__setattr__(self, "calls", [])
        object.__setattr__(self, "fail", fail)
        object.__setattr__(self, "endpoint_changed", endpoint_changed)
        object.__setattr__(self, "executable", Path("/server/bin/docker"))
        object.__setattr__(self, "unix_socket", Path("/server/run/docker.sock"))

    def verify_endpoint(self):
        if self.endpoint_changed:
            return Result.failure(
                [
                    Issue(
                        code="MANAGED_CONTAINER_ENDPOINT_CHANGED",
                        message="Managed container endpoint changed.",
                        path="execution.cleanup",
                        source="managed_container_cleaner",
                    )
                ]
            )
        return Result.success(None)

    def cleanup(self, scope):
        self.calls.append(scope)
        if self.fail:
            return Result.failure(
                [
                    Issue(
                        code="MANAGED_CONTAINER_CLEANUP_FAILED",
                        message="Managed workflow containers could not be cleaned up.",
                        path="execution.cleanup",
                        source="managed_container_cleaner",
                    )
                ]
            )
        return Result.success(None)


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
    assert result.value.stdout.strip() == "[REDACTED]"


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
        "DOCKER_HOST",
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

    def track_termination(self, process):
        terminated.append(process)
        return original_terminate(self, process)

    monkeypatch.setattr(
        ProcessRunner,
        "_terminate",
        track_termination,
    )
    runner = _make_runner(timeout_seconds=0.1)

    result = runner.run(
        CommandSpec(argv=(sys.executable, "-c", "import time; time.sleep(10)"))
    )

    assert result.is_failure
    assert result.issues[0].code == "PROCESS_RUNNER_TIMEOUT"
    assert len(terminated) == 1
    assert terminated[0].poll() is not None


def test_process_runner_timeout_reaps_a_real_descendant_tree(tmp_path):
    pid_file = tmp_path / "grandchild.pid"
    grandchild = "import time; time.sleep(30)"
    parent = (
        "import pathlib, subprocess, sys, time; "
        f"child=subprocess.Popen([sys.executable, '-c', {grandchild!r}]); "
        f"pathlib.Path({str(pid_file)!r}).write_text(str(child.pid)); "
        "time.sleep(30)"
    )

    result = _make_runner(timeout_seconds=0.2).run(
        CommandSpec(argv=(sys.executable, "-c", parent))
    )

    assert result.is_failure
    assert result.issues[0].code == "PROCESS_RUNNER_TIMEOUT"
    grandchild_pid = int(pid_file.read_text())
    deadline = time.monotonic() + 2
    while time.monotonic() < deadline:
        try:
            state = Path(f"/proc/{grandchild_pid}/stat").read_text().split()[2]
        except (FileNotFoundError, ProcessLookupError):
            break
        if state == "Z":
            break
        time.sleep(0.02)
    else:
        os.kill(grandchild_pid, 0)
        pytest.fail("timed-out workflow descendant remained alive")


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
    assert chunks == ["[REDACTED]"]


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


def test_process_runner_redacts_unicode_private_prefix_at_byte_limit():
    private_reference = "/private/参考/genome.fa"
    visible_prefix = "visible:"
    private_bytes = private_reference.encode("utf-8")
    first_multibyte = private_bytes.index("参".encode("utf-8"))
    limit = len(visible_prefix.encode("utf-8")) + first_multibyte + 1
    streamed: list[str] = []
    script = (
        "import sys; "
        f"sys.stdout.buffer.write({(visible_prefix + private_reference).encode()!r})"
    )

    result = _make_runner(max_output_bytes=limit).run(
        CommandSpec(
            argv=(sys.executable, "-c", script),
            redaction_values=(private_reference,),
        ),
        output_callback=lambda _stream, lines: streamed.extend(lines),
    )

    assert result.is_success
    assert result.value.stdout.endswith("[REDACTED]")
    assert "".join(streamed) == result.value.stdout
    assert "/private/" not in result.value.stdout
    assert "/private/" not in "".join(streamed)


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

    def track_termination(self, process):
        terminated.append(process)
        return original_terminate(self, process)

    monkeypatch.setattr(
        ProcessRunner,
        "_terminate",
        track_termination,
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


def test_process_runner_redacts_adapter_private_values_from_capture_and_streams():
    private_reference = "/private/references/genome.fa"
    streamed: list[str] = []
    script = (
        "import sys; "
        f"print({private_reference!r}); "
        f"sys.stderr.write({(private_reference + '/index')!r} + '\\n')"
    )

    result = _make_runner().run(
        CommandSpec(
            argv=(sys.executable, "-c", script),
            redaction_values=(private_reference,),
        ),
        output_callback=lambda _stream, lines: streamed.extend(lines),
    )

    assert result.is_success
    assert private_reference not in result.value.stdout
    assert private_reference not in result.value.stderr
    assert private_reference not in "\n".join(streamed)
    assert "[REDACTED]" in result.value.stdout
    assert "[REDACTED]/index" in result.value.stderr


def test_process_runner_redacts_private_prefix_at_output_limit():
    private_reference = "/private/references/genome.fa"
    visible_prefix = "visible:"
    accepted = visible_prefix + private_reference[:-1]
    streamed: list[str] = []
    script = f"print({(visible_prefix + private_reference)!r}, end='')"

    result = _make_runner(max_output_bytes=len(accepted)).run(
        CommandSpec(
            argv=(sys.executable, "-c", script),
            redaction_values=(private_reference,),
        ),
        output_callback=lambda _stream, lines: streamed.extend(lines),
    )

    assert result.is_success
    assert result.value.stdout.endswith("[REDACTED]")
    assert "".join(streamed) == result.value.stdout
    assert private_reference[:-1] not in result.value.stdout
    assert private_reference[:-1] not in "".join(streamed)
    assert any(
        issue.code == "PROCESS_RUNNER_OUTPUT_LIMIT_EXCEEDED" for issue in result.issues
    )


@pytest.mark.parametrize("trailer", (b"\n", b""), ids=("newline", "eof"))
def test_process_runner_redacts_secret_split_at_stream_chunk_boundary(trailer):
    secret = "private-boundary-token"
    split = 7
    streamed: list[str] = []
    script = (
        "import os, time; "
        f"secret = {secret.encode()!r}; "
        f"os.write(1, b'x' * ({64 * 1024} - {split}) + secret[:{split}]); "
        "time.sleep(0.05); "
        f"os.write(1, secret[{split}:] + {trailer!r})"
    )

    result = _make_runner().run(
        CommandSpec(
            argv=(sys.executable, "-c", script),
            redaction_values=(secret,),
        ),
        output_callback=lambda _stream, lines: streamed.extend(lines),
    )

    assert result.is_success
    rendered_stream = "".join(streamed)
    assert secret not in rendered_stream
    assert secret not in result.value.stdout
    assert rendered_stream == "x" * (64 * 1024 - split) + "[REDACTED]"
    assert result.value.stdout.rstrip("\n") == rendered_stream


def test_streaming_redactor_masks_overlapping_private_values_across_chunks():
    private_values = ("bdb", "bb", "a")
    chunks = (
        "cccababbacb",
        "abbaab",
        "bcdd",
        "bbdbb",
        "dbdc",
        "dbad",
        "cbd",
        "c",
        "bbdaaa",
    )
    redactor = process_runner_module._StreamingRedactor(
        process_runner_module._BoundedLiteralMatcher(private_values)
    )

    rendered = "".join(redactor.feed(chunk) for chunk in chunks)
    rendered += redactor.finish()

    assert all(private not in rendered for private in private_values)


def test_process_runner_redacts_more_than_256_literals_across_stream_boundary():
    private_values = tuple(
        f"/private/run/sample-{index:04d}/reads.fastq.gz" for index in range(2_205)
    )
    secret = private_values[-1]
    split = 11
    streamed: list[str] = []
    script = (
        "import os, time; "
        f"secret = {secret.encode()!r}; "
        f"os.write(1, b'x' * ({64 * 1024} - {split}) + secret[:{split}]); "
        "time.sleep(0.05); "
        f"os.write(1, secret[{split}:] + b'\\n')"
    )

    result = _make_runner().run(
        CommandSpec(
            argv=(sys.executable, "-c", script),
            redaction_values=private_values,
        ),
        output_callback=lambda _stream, lines: streamed.extend(lines),
    )

    assert result.is_success
    rendered_stream = "".join(streamed)
    assert secret not in rendered_stream
    assert secret not in result.value.stdout
    assert rendered_stream == "x" * (64 * 1024 - split) + "[REDACTED]"
    assert result.value.stdout.rstrip("\n") == rendered_stream


def test_process_runner_rejects_one_oversized_redaction_literal_before_launch(
    tmp_path,
):
    marker = tmp_path / "started"
    oversized = "x" * (process_runner_module._MAX_STREAM_REDACTION_PATTERN_LENGTH + 1)

    result = _make_runner().run(
        CommandSpec(
            argv=(sys.executable, "-c", f"open({str(marker)!r}, 'w').close()"),
            redaction_values=(oversized,),
        )
    )

    assert result.is_failure
    assert result.issues[0].code == "PROCESS_RUNNER_REDACTION_BOUND_EXCEEDED"
    assert not marker.exists()


def test_process_runner_rejects_oversized_total_redaction_contract_before_launch(
    tmp_path,
):
    marker = tmp_path / "started"
    literal_length = process_runner_module._MAX_STREAM_REDACTION_PATTERN_LENGTH
    count = (
        process_runner_module._MAX_STREAM_REDACTION_TOTAL_CHARACTERS // literal_length
    )
    private_values = tuple(
        f"{index:08x}" + "x" * (literal_length - 8) for index in range(count + 1)
    )

    result = _make_runner().run(
        CommandSpec(
            argv=(sys.executable, "-c", f"open({str(marker)!r}, 'w').close()"),
            redaction_values=private_values,
        )
    )

    assert result.is_failure
    assert result.issues[0].code == "PROCESS_RUNNER_REDACTION_BOUND_EXCEEDED"
    assert not marker.exists()


def test_process_runner_accepts_near_maximum_redaction_contract_with_bounded_heap():
    literal_length = process_runner_module._MAX_STREAM_REDACTION_PATTERN_LENGTH
    reserve_for_runner_values = 16_384
    count = (
        process_runner_module._MAX_STREAM_REDACTION_TOTAL_CHARACTERS
        - reserve_for_runner_values
    ) // literal_length
    private_values = tuple(
        f"{index:08x}" + "x" * (literal_length - 8) for index in range(count)
    )
    secret = private_values[-1]

    tracemalloc.start()
    try:
        result = _make_runner().run(
            CommandSpec(
                argv=(sys.executable, "-c", f"print({secret!r})"),
                redaction_values=private_values,
            )
        )
        _current, peak = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()

    assert sum(map(len, private_values)) > 15 * 1024 * 1024
    assert result.is_success
    assert result.value.stdout == "[REDACTED]\n"
    assert peak < process_runner_module._MAX_STREAM_REDACTION_TOTAL_CHARACTERS


def test_process_runner_fails_closed_on_adversarial_redaction_candidate_density():
    literal_length = process_runner_module._MAX_STREAM_REDACTION_PATTERN_LENGTH
    private_values = tuple(
        "a" * 100 + f"{index:08x}" + "a" * (literal_length - 108)
        for index in range(300)
    )
    streamed: list[str] = []

    result = _make_runner().run(
        CommandSpec(
            argv=(sys.executable, "-c", f"print({'a' * literal_length!r})"),
            redaction_values=private_values,
        ),
        output_callback=lambda _stream, lines: streamed.extend(lines),
    )

    assert result.is_failure
    assert result.issues[0].code == "PROCESS_RUNNER_REDACTION_WORK_EXCEEDED"
    assert streamed == []


def test_process_runner_does_not_call_docker_without_declared_scope():
    cleaner = _RecordingContainerCleaner()

    result = _make_runner(managed_container_cleaner=cleaner).run(
        CommandSpec(argv=(sys.executable, "-c", "pass"))
    )

    assert result.is_success
    assert cleaner.calls == []


def test_process_runner_injects_only_server_owned_local_docker_boundary(monkeypatch):
    monkeypatch.setenv("DOCKER_HOST", "tcp://remote.example:2375")
    cleaner = _RecordingContainerCleaner()
    scope = "0" * 64
    script = (
        "import os; "
        "host = os.environ.get('DOCKER_HOST'); "
        "path = os.environ.get('PATH', '').split(os.pathsep)[0]; "
        "assert host == 'unix:///server/run/docker.sock'; "
        "assert path == '/server/bin'; "
        "print(host); print(path)"
    )

    result = _make_runner(managed_container_cleaner=cleaner).run(
        CommandSpec(
            argv=(sys.executable, "-c", script),
            managed_container_scope=scope,
            managed_container_endpoint_identity=_MANAGED_ENDPOINT_IDENTITY,
        )
    )

    assert result.is_success
    assert result.value.stdout.splitlines() == ["[REDACTED]", "[REDACTED]"]
    assert "/server" not in result.value.stdout
    assert "remote.example" not in result.value.stdout
    assert cleaner.calls == [scope]


def test_process_runner_refuses_managed_scope_without_cleaner(tmp_path):
    marker = tmp_path / "started"

    result = _make_runner().run(
        CommandSpec(
            argv=(sys.executable, "-c", f"open({str(marker)!r}, 'w').close()"),
            managed_container_scope="a" * 64,
            managed_container_endpoint_identity=_MANAGED_ENDPOINT_IDENTITY,
        )
    )

    assert result.is_failure
    assert result.issues[0].code == "PROCESS_RUNNER_CONTAINER_CLEANER_UNAVAILABLE"
    assert not marker.exists()


def test_process_runner_refuses_mismatched_managed_container_endpoint(tmp_path):
    cleaner = _RecordingContainerCleaner()
    marker = tmp_path / "started"
    mismatched_identity = "0" * 64
    assert mismatched_identity != cleaner.endpoint_identity

    result = _make_runner(managed_container_cleaner=cleaner).run(
        CommandSpec(
            argv=(sys.executable, "-c", f"open({str(marker)!r}, 'w').close()"),
            managed_container_scope="a" * 64,
            managed_container_endpoint_identity=mismatched_identity,
        )
    )

    assert result.is_failure
    assert result.issues[0].code == "PROCESS_RUNNER_CONTAINER_ENDPOINT_MISMATCH"
    assert not marker.exists()
    assert cleaner.calls == []
    rendered = str(result.issues[0].to_dict())
    assert mismatched_identity not in rendered
    assert "/server" not in rendered


def test_process_runner_rechecks_managed_endpoint_before_launch(tmp_path):
    cleaner = _RecordingContainerCleaner(endpoint_changed=True)
    marker = tmp_path / "started"

    result = _make_runner(managed_container_cleaner=cleaner).run(
        CommandSpec(
            argv=(sys.executable, "-c", f"open({str(marker)!r}, 'w').close()"),
            managed_container_scope="a" * 64,
            managed_container_endpoint_identity=_MANAGED_ENDPOINT_IDENTITY,
        )
    )

    assert result.is_failure
    assert result.issues[0].code == "PROCESS_RUNNER_CONTAINER_ENDPOINT_CHANGED"
    assert not marker.exists()
    assert cleaner.calls == ["a" * 64]
    assert "/server" not in str(result.issues[0].to_dict())


def test_process_runner_timeout_cleans_declared_container_scope():
    cleaner = _RecordingContainerCleaner()
    scope = "b" * 64

    result = _make_runner(
        timeout_seconds=0.1,
        managed_container_cleaner=cleaner,
    ).run(
        CommandSpec(
            argv=(sys.executable, "-c", "import time; time.sleep(10)"),
            managed_container_scope=scope,
            managed_container_endpoint_identity=_MANAGED_ENDPOINT_IDENTITY,
        )
    )

    assert result.is_failure
    assert result.issues[0].code == "PROCESS_RUNNER_TIMEOUT"
    assert cleaner.calls == [scope]


def test_process_runner_callback_failure_cleans_declared_container_scope():
    cleaner = _RecordingContainerCleaner()
    scope = "c" * 64

    result = _make_runner(managed_container_cleaner=cleaner).run(
        CommandSpec(
            argv=(sys.executable, "-c", "print('ready', flush=True)"),
            managed_container_scope=scope,
            managed_container_endpoint_identity=_MANAGED_ENDPOINT_IDENTITY,
        ),
        output_callback=lambda _stream, _lines: (_ for _ in ()).throw(
            RuntimeError("private callback failure")
        ),
    )

    assert result.is_failure
    assert result.issues[0].code == "PROCESS_RUNNER_OUTPUT_CALLBACK_FAILED"
    assert cleaner.calls == [scope]
    assert "private callback failure" not in str(result.issues[0].to_dict())


def test_process_runner_base_exception_cleans_and_preserves_control_flow():
    class AbortExecution(BaseException):
        pass

    cleaner = _RecordingContainerCleaner()
    scope = "d" * 64

    with pytest.raises(AbortExecution):
        _make_runner(managed_container_cleaner=cleaner).run(
            CommandSpec(
                argv=(
                    sys.executable,
                    "-c",
                    "import time; print('ready', flush=True); time.sleep(10)",
                ),
                managed_container_scope=scope,
                managed_container_endpoint_identity=_MANAGED_ENDPOINT_IDENTITY,
            ),
            output_callback=lambda _stream, _lines: (_ for _ in ()).throw(
                AbortExecution()
            ),
        )

    assert cleaner.calls == [scope]


def test_process_runner_base_exception_cleanup_failure_is_fail_closed():
    class AbortExecution(BaseException):
        pass

    cleaner = _RecordingContainerCleaner(fail=True)

    with pytest.raises(ProcessRunnerCleanupError):
        _make_runner(managed_container_cleaner=cleaner).run(
            CommandSpec(
                argv=(
                    sys.executable,
                    "-c",
                    "import time; print('ready', flush=True); time.sleep(10)",
                ),
                managed_container_scope="e" * 64,
                managed_container_endpoint_identity=_MANAGED_ENDPOINT_IDENTITY,
            ),
            output_callback=lambda _stream, _lines: (_ for _ in ()).throw(
                AbortExecution()
            ),
        )


def test_container_cleanup_failure_is_primary_after_nonzero_exit():
    cleaner = _RecordingContainerCleaner(fail=True)

    result = _make_runner(managed_container_cleaner=cleaner).run(
        CommandSpec(
            argv=(sys.executable, "-c", "raise SystemExit(3)"),
            managed_container_scope="1" * 64,
            managed_container_endpoint_identity=_MANAGED_ENDPOINT_IDENTITY,
        )
    )

    assert result.is_failure
    assert [issue.code for issue in result.issues] == [
        "MANAGED_CONTAINER_CLEANUP_FAILED",
        "PROCESS_RUNNER_NONZERO_EXIT",
    ]


@pytest.mark.skipif(not sys.platform.startswith("linux"), reason="Linux /proc gate")
def test_timeout_reaps_descendant_forked_late_from_term_handler(tmp_path):
    cleaner = _RecordingContainerCleaner()
    child_pid_path = tmp_path / "late-child.pid"
    script = f"""
import os
from pathlib import Path
import signal
import time

child_pid_path = Path({str(child_pid_path)!r})
forked = False

def terminate(_signum, _frame):
    global forked
    if forked:
        return
    forked = True
    pid = os.fork()
    if pid == 0:
        os.setsid()
        signal.signal(signal.SIGTERM, signal.SIG_IGN)
        os.close(1)
        os.close(2)
        while True:
            time.sleep(1)
    child_pid_path.write_text(str(pid), encoding="utf-8")

signal.signal(signal.SIGTERM, terminate)
print("ready", flush=True)
while True:
    time.sleep(1)
"""

    result = _make_runner(
        timeout_seconds=0.15,
        managed_container_cleaner=cleaner,
    ).run(
        CommandSpec(
            argv=(sys.executable, "-c", script),
            managed_container_scope="f" * 64,
            managed_container_endpoint_identity=_MANAGED_ENDPOINT_IDENTITY,
        )
    )

    assert result.is_failure
    assert result.issues[0].code == "PROCESS_RUNNER_TIMEOUT"
    child_pid = int(child_pid_path.read_text(encoding="utf-8"))
    with pytest.raises(ProcessLookupError):
        os.kill(child_pid, 0)
    assert cleaner.calls == ["f" * 64]


@pytest.mark.skipif(not sys.platform.startswith("linux"), reason="Linux /proc gate")
def test_bind_failure_after_fast_root_exit_reaps_adopted_setsid_child(
    tmp_path,
    monkeypatch,
):
    cleaner = _RecordingContainerCleaner()
    child_pid_path = tmp_path / "adopted-child.pid"
    script = f"""
import os
from pathlib import Path
import time

pid = os.fork()
if pid == 0:
    os.setsid()
    os.close(1)
    os.close(2)
    while True:
        time.sleep(1)
Path({str(child_pid_path)!r}).write_text(str(pid), encoding="utf-8")
os._exit(0)
"""

    def fail_bind_after_root_exit(self, root_pid):
        deadline = time.monotonic() + 2
        while time.monotonic() < deadline:
            if child_pid_path.exists() and not Path(f"/proc/{root_pid}").exists():
                break
            time.sleep(0.01)
        raise process_runner_module._SubreaperUnavailable(
            "simulated fast-root bind failure"
        )

    monkeypatch.setattr(
        process_runner_module._LinuxSubreaper,
        "bind",
        fail_bind_after_root_exit,
    )

    result = _make_runner(managed_container_cleaner=cleaner).run(
        CommandSpec(
            argv=(sys.executable, "-c", script),
            managed_container_scope="9" * 64,
            managed_container_endpoint_identity=_MANAGED_ENDPOINT_IDENTITY,
        )
    )

    assert result.is_failure
    assert result.issues[0].code == "PROCESS_RUNNER_SUBREAPER_UNAVAILABLE"
    child_pid = int(child_pid_path.read_text(encoding="utf-8"))
    with pytest.raises(ProcessLookupError):
        os.kill(child_pid, 0)
    assert cleaner.calls == ["9" * 64]


@pytest.mark.skipif(not sys.platform.startswith("linux"), reason="Linux prctl gate")
def test_process_runner_restores_host_subreaper_state():
    libc = ctypes.CDLL(None, use_errno=True)

    def current_state():
        value = ctypes.c_int()
        assert libc.prctl(37, ctypes.byref(value), 0, 0, 0) == 0
        return value.value

    before = current_state()

    result = _make_runner().run(CommandSpec(argv=(sys.executable, "-c", "pass")))

    assert result.is_success
    assert current_state() == before


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
