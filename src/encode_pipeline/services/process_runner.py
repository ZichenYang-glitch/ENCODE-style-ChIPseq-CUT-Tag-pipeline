"""Safe subprocess boundary for executing already-built CommandSpec instances."""

from __future__ import annotations

import codecs
import os
import selectors
import subprocess
import time
from collections.abc import Callable, Mapping
from dataclasses import dataclass
from typing import TYPE_CHECKING

from encode_pipeline.platform.results import Issue, Result

if TYPE_CHECKING:
    from encode_pipeline.platform.adapters import CommandSpec


OutputCallback = Callable[[str, tuple[str, ...]], None]


@dataclass(frozen=True)
class ProcessResult:
    """The outcome of a successfully-run subprocess."""

    exit_code: int
    stdout: str
    stderr: str


# Only runtime variables needed to locate and run the scientific toolchain are
# inherited. In particular, worker/API configuration and provider credentials
# must never become part of a Snakemake process (or one of its rule children).
_INHERITED_ENVIRONMENT_NAMES = frozenset(
    {
        "CONDA_DEFAULT_ENV",
        "CONDA_EXE",
        "CONDA_PREFIX",
        "CONDA_PYTHON_EXE",
        "CONDA_SHLVL",
        "CURL_CA_BUNDLE",
        "HOME",
        "JAVA_HOME",
        "LANG",
        "LANGUAGE",
        "LD_LIBRARY_PATH",
        "LC_ADDRESS",
        "LC_ALL",
        "LC_COLLATE",
        "LC_CTYPE",
        "LC_IDENTIFICATION",
        "LC_MEASUREMENT",
        "LC_MESSAGES",
        "LC_MONETARY",
        "LC_NAME",
        "LC_NUMERIC",
        "LC_PAPER",
        "LC_TELEPHONE",
        "LC_TIME",
        "LOGNAME",
        "MAMBA_EXE",
        "MAMBA_ROOT_PREFIX",
        "NO_COLOR",
        "PATH",
        "REQUESTS_CA_BUNDLE",
        "SAMTOOLS",
        "SHELL",
        "SSL_CERT_DIR",
        "SSL_CERT_FILE",
        "TEMP",
        "TERM",
        "TMP",
        "TMPDIR",
        "TZ",
        "USER",
        "XDG_CACHE_HOME",
        "XDG_DATA_DIRS",
        "_CE_CONDA",
        "_CE_M",
        "_CONDA_EXE",
        "_CONDA_ROOT",
    }
)
_PROTECTED_ENVIRONMENT_PREFIXES = (
    "ENCODE_PIPELINE_",
    "ENCODE_AGENT_",
    "ANTHROPIC_",
    "AWS_",
    "AZURE_",
    "CLAUDE_",
    "CODEX_",
    "DEEPSEEK_",
    "GEMINI_",
    "GITHUB_",
    "GOOGLE_",
    "OPENAI_",
)
_PROTECTED_ENVIRONMENT_SUFFIXES = (
    "_API_KEY",
    "_AUTH",
    "_CREDENTIAL",
    "_CREDENTIALS",
    "_PASSWORD",
    "_SECRET",
    "_TOKEN",
)
_PROTECTED_ENVIRONMENT_EXACT_NAMES = frozenset(
    {
        "DATABASE_URL",
        "REDIS_URL",
    }
)
_PROTECTED_ENVIRONMENT_FRAGMENTS = (
    "ACCESS_KEY",
    "API_KEY",
    "AUTH_TOKEN",
    "AUTHORIZATION",
    "CLIENT_SECRET",
    "PRIVATE_KEY",
)
_PROTECTED_ENVIRONMENT_SEGMENTS = frozenset(
    {
        "AUTH",
        "CREDENTIAL",
        "CREDENTIALS",
        "PASSWORD",
        "SECRET",
        "TOKEN",
    }
)
_READ_SIZE = 64 * 1024
_POST_TERMINATION_DRAIN_MAX_BYTES = 256 * 1024
_POST_TERMINATION_DRAIN_MAX_ITERATIONS = 8
_POST_TERMINATION_DRAIN_MAX_SECONDS = 0.05


def _is_protected_environment_name(name: str) -> bool:
    normalized = name.upper()
    segments = frozenset(normalized.split("_"))
    return (
        normalized in _PROTECTED_ENVIRONMENT_EXACT_NAMES
        or normalized.endswith(("_DATABASE_URL", "_REDIS_URL"))
        or normalized.startswith(_PROTECTED_ENVIRONMENT_PREFIXES)
        or normalized.endswith(_PROTECTED_ENVIRONMENT_SUFFIXES)
        or any(fragment in normalized for fragment in _PROTECTED_ENVIRONMENT_FRAGMENTS)
        or bool(segments & _PROTECTED_ENVIRONMENT_SEGMENTS)
    )


def _subprocess_environment(overrides: Mapping[str, str]) -> Result[dict[str, str]]:
    if any(_is_protected_environment_name(name) for name in overrides):
        return Result.failure(
            [
                Issue(
                    code="PROCESS_RUNNER_PROTECTED_ENVIRONMENT",
                    message="CommandSpec attempted to set a protected environment variable.",
                    severity="error",
                    path="command_spec.env",
                    source="process_runner",
                )
            ]
        )

    inherited = {
        name: value
        for name, value in os.environ.items()
        if name in _INHERITED_ENVIRONMENT_NAMES
    }
    inherited.update(overrides)
    return Result.success(inherited)


def _append_capture(
    capture: bytearray,
    data: bytes,
    max_bytes: int,
) -> tuple[bytes, bool]:
    """Append as much as fits and return the accepted bytes plus truncation."""
    remaining = max_bytes - len(capture)
    accepted = data[: max(remaining, 0)]
    if remaining > 0:
        capture.extend(accepted)
    return accepted, len(data) > max(remaining, 0)


def _decode(data: bytes) -> str:
    return data.decode("utf-8", errors="replace")


class _OutputLineDecoder:
    """Incrementally decode one stream and emit bounded, complete text chunks."""

    def __init__(self) -> None:
        decoder_factory = codecs.getincrementaldecoder("utf-8")
        self._decoder = decoder_factory(errors="replace")
        self._pending = ""
        self._finished = False

    def feed(self, data: bytes) -> tuple[str, ...]:
        if self._finished:
            raise RuntimeError("output decoder is already finished")
        self._pending += self._decoder.decode(data, final=False)
        return self._take_chunks(final=False)

    def finish(self) -> tuple[str, ...]:
        if self._finished:
            return ()
        self._finished = True
        self._pending += self._decoder.decode(b"", final=True)
        return self._take_chunks(final=True)

    def _take_chunks(self, *, final: bool) -> tuple[str, ...]:
        chunks: list[str] = []
        while self._pending:
            newline = self._pending.find("\n")
            if 0 <= newline <= _READ_SIZE:
                line = self._pending[:newline]
                self._pending = self._pending[newline + 1 :]
                if line.endswith("\r"):
                    line = line[:-1]
                chunks.append(line)
                continue
            if len(self._pending) >= _READ_SIZE:
                chunks.append(self._pending[:_READ_SIZE])
                self._pending = self._pending[_READ_SIZE:]
                continue
            break
        if final and self._pending:
            chunks.append(self._pending)
            self._pending = ""
        return tuple(chunks)


class ProcessRunner:
    """Execute a controlled command while draining and optionally streaming logs.

    The executable is allowlisted, inherited environment variables are reduced
    to a non-secret runtime allowlist, captured output is bounded, and stdout and
    stderr are drained concurrently. Non-zero exits remain structured successes
    with a warning so the lifecycle owner can map them to a durable run failure.
    """

    def __init__(
        self,
        *,
        allowed_executables: tuple[str, ...] = ("snakemake",),
        timeout_seconds: float = 300.0,
        max_output_bytes: int = 10_000_000,
        passthrough_exceptions: tuple[type[BaseException], ...] = (),
    ) -> None:
        if not isinstance(allowed_executables, tuple):
            raise ValueError("allowed_executables must be a tuple")
        if len(allowed_executables) == 0:
            raise ValueError("allowed_executables must be non-empty")
        for entry in allowed_executables:
            if not isinstance(entry, str) or not entry.strip():
                raise ValueError(
                    "allowed_executables entries must be non-empty strings"
                )
        if isinstance(timeout_seconds, bool) or not isinstance(
            timeout_seconds, (int, float)
        ):
            raise ValueError("timeout_seconds must be a number")
        if timeout_seconds <= 0:
            raise ValueError("timeout_seconds must be positive")
        if isinstance(max_output_bytes, bool) or not isinstance(max_output_bytes, int):
            raise ValueError("max_output_bytes must be an int")
        if max_output_bytes <= 0:
            raise ValueError("max_output_bytes must be positive")
        if not isinstance(passthrough_exceptions, tuple):
            raise ValueError("passthrough_exceptions must be a tuple")
        for exception_type in passthrough_exceptions:
            if not isinstance(exception_type, type) or not issubclass(
                exception_type, BaseException
            ):
                raise ValueError(
                    "passthrough_exceptions entries must be exception types"
                )

        self._allowed_executables = allowed_executables
        self._timeout_seconds = float(timeout_seconds)
        self._max_output_bytes = max_output_bytes
        self._passthrough_exceptions = passthrough_exceptions

    def run(
        self,
        spec: "CommandSpec",
        *,
        output_callback: OutputCallback | None = None,
    ) -> "Result[ProcessResult]":
        """Execute *spec* and return a structured, bounded result.

        When supplied, ``output_callback`` receives bounded line chunks as the
        child writes them. A callback failure terminates the direct child and is
        returned as a public-safe infrastructure failure.
        """
        from encode_pipeline.platform.adapters import CommandSpec

        if not isinstance(spec, CommandSpec):
            return Result.failure(
                [
                    Issue(
                        code="PROCESS_RUNNER_INVALID_COMMAND_SPEC",
                        message="spec must be a CommandSpec.",
                        severity="error",
                        path="command_spec",
                        source="process_runner",
                    )
                ]
            )
        if output_callback is not None and not callable(output_callback):
            return Result.failure(
                [
                    Issue(
                        code="PROCESS_RUNNER_INVALID_OUTPUT_CALLBACK",
                        message="output_callback must be callable.",
                        severity="error",
                        path="output_callback",
                        source="process_runner",
                    )
                ]
            )

        if spec.argv[0] not in self._allowed_executables:
            return Result.failure(
                [
                    Issue(
                        code="PROCESS_RUNNER_EXECUTABLE_NOT_ALLOWED",
                        message="Executable is not in the allowlist.",
                        severity="error",
                        path="command_spec.argv[0]",
                        source="process_runner",
                    )
                ]
            )

        environment_result = _subprocess_environment(spec.env)
        if environment_result.is_failure:
            return environment_result

        try:
            process = subprocess.Popen(
                spec.argv,
                cwd=spec.cwd,
                env=environment_result.value,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
        except FileNotFoundError:
            return Result.failure(
                [
                    Issue(
                        code="PROCESS_RUNNER_EXECUTABLE_NOT_FOUND",
                        message="Executable not found.",
                        severity="error",
                        path="command_spec.argv[0]",
                        source="process_runner",
                    )
                ]
            )
        except OSError:
            return Result.failure(
                [
                    Issue(
                        code="PROCESS_RUNNER_EXECUTION_ERROR",
                        message="Subprocess execution failed with an OS error.",
                        severity="error",
                        path="command_spec",
                        source="process_runner",
                    )
                ]
            )

        assert process.stdout is not None
        assert process.stderr is not None
        streams = {
            process.stdout.fileno(): ("stdout", process.stdout),
            process.stderr.fileno(): ("stderr", process.stderr),
        }
        captures = {"stdout": bytearray(), "stderr": bytearray()}
        decoders = {
            "stdout": _OutputLineDecoder(),
            "stderr": _OutputLineDecoder(),
        }
        truncated = {"stdout": False, "stderr": False}
        selector = selectors.DefaultSelector()
        started_at = time.monotonic()

        def consume(stream_name: str, data: bytes) -> bool:
            accepted, was_truncated = _append_capture(
                captures[stream_name], data, self._max_output_bytes
            )
            truncated[stream_name] = was_truncated or truncated[stream_name]
            if output_callback is None or not accepted:
                return True
            chunks = decoders[stream_name].feed(accepted)
            return not chunks or self._emit(
                process,
                output_callback,
                stream_name,
                chunks,
            )

        def drain_available_after_termination() -> bool:
            """Drain bytes ready now without waiting for descendant-held pipes."""
            drained_bytes = 0
            iterations = 0
            deadline = time.monotonic() + _POST_TERMINATION_DRAIN_MAX_SECONDS
            while (
                selector.get_map()
                and drained_bytes < _POST_TERMINATION_DRAIN_MAX_BYTES
                and iterations < _POST_TERMINATION_DRAIN_MAX_ITERATIONS
                and time.monotonic() < deadline
            ):
                iterations += 1
                try:
                    ready = selector.select(timeout=0)
                except OSError:
                    return True
                if not ready:
                    return True
                progressed = False
                for key, _ in ready:
                    remaining_drain_bytes = (
                        _POST_TERMINATION_DRAIN_MAX_BYTES - drained_bytes
                    )
                    if remaining_drain_bytes <= 0 or time.monotonic() >= deadline:
                        return True
                    descriptor = key.fd
                    stream_name, _pipe = streams[descriptor]
                    try:
                        data = os.read(
                            descriptor,
                            min(_READ_SIZE, remaining_drain_bytes),
                        )
                    except BlockingIOError:
                        continue
                    except OSError:
                        selector.unregister(descriptor)
                        progressed = True
                        continue
                    if not data:
                        selector.unregister(descriptor)
                        progressed = True
                        continue
                    progressed = True
                    drained_bytes += len(data)
                    if not consume(stream_name, data):
                        return False
                if not progressed:
                    return True
            return True

        try:
            for descriptor in streams:
                os.set_blocking(descriptor, False)
                selector.register(descriptor, selectors.EVENT_READ)

            while selector.get_map() or process.poll() is None:
                remaining = self._timeout_seconds - (time.monotonic() - started_at)
                if remaining <= 0:
                    self._terminate(process)
                    if not drain_available_after_termination():
                        return self._callback_failure(truncated)
                    if output_callback is not None:
                        for stream_name in ("stdout", "stderr"):
                            chunks = decoders[stream_name].finish()
                            if chunks and not self._emit(
                                process,
                                output_callback,
                                stream_name,
                                chunks,
                            ):
                                return self._callback_failure(truncated)
                    return self._failure_with_output_limit(
                        Issue(
                            code="PROCESS_RUNNER_TIMEOUT",
                            message="Subprocess timed out.",
                            severity="error",
                            path="command_spec",
                            source="process_runner",
                            context={"timeout_seconds": self._timeout_seconds},
                        ),
                        truncated,
                    )

                for key, _ in selector.select(timeout=min(remaining, 0.1)):
                    descriptor = key.fd
                    stream_name, _pipe = streams[descriptor]
                    try:
                        data = os.read(descriptor, _READ_SIZE)
                    except BlockingIOError:
                        continue
                    if not data:
                        selector.unregister(descriptor)
                        if output_callback is not None:
                            chunks = decoders[stream_name].finish()
                            if chunks and not self._emit(
                                process,
                                output_callback,
                                stream_name,
                                chunks,
                            ):
                                return self._callback_failure(truncated)
                        continue

                    if not consume(stream_name, data):
                        return self._callback_failure(truncated)

            exit_code = process.wait()
        except OSError:
            self._terminate(process)
            return self._failure_with_output_limit(
                Issue(
                    code="PROCESS_RUNNER_EXECUTION_ERROR",
                    message="Subprocess execution failed with an OS error.",
                    severity="error",
                    path="command_spec",
                    source="process_runner",
                ),
                truncated,
            )
        finally:
            # RQ's timeout is delivered asynchronously. Any exception, including
            # a BaseException raised between reads, must not leave the directly
            # spawned workflow process alive. Process-group cleanup is PR125.
            if process.poll() is None:
                self._terminate(process)
            selector.close()
            process.stdout.close()
            process.stderr.close()

        issues: list[Issue] = []
        if any(truncated.values()):
            issues.append(self._output_limit_issue())
        if exit_code != 0:
            issues.append(
                Issue(
                    code="PROCESS_RUNNER_NONZERO_EXIT",
                    message="Subprocess exited with a non-zero code.",
                    severity="warning",
                    path="command_spec",
                    source="process_runner",
                    context={"exit_code": exit_code},
                )
            )

        return Result.success(
            ProcessResult(
                exit_code=exit_code,
                stdout=_decode(bytes(captures["stdout"])),
                stderr=_decode(bytes(captures["stderr"])),
            ),
            issues=issues,
        )

    def _emit(
        self,
        process: subprocess.Popen[bytes],
        callback: OutputCallback,
        stream_name: str,
        lines: tuple[str, ...],
    ) -> bool:
        try:
            callback(stream_name, lines)
        except self._passthrough_exceptions:
            raise
        except Exception:
            ProcessRunner._terminate(process)
            return False
        return True

    @staticmethod
    def _terminate(process: subprocess.Popen[bytes]) -> None:
        if process.poll() is None:
            process.kill()
        try:
            process.wait(timeout=5)
        except (
            subprocess.TimeoutExpired
        ):  # pragma: no cover - kill is immediate on POSIX
            pass

    @classmethod
    def _callback_failure(
        cls,
        truncated: Mapping[str, bool],
    ) -> Result[ProcessResult]:
        return cls._failure_with_output_limit(
            Issue(
                code="PROCESS_RUNNER_OUTPUT_CALLBACK_FAILED",
                message="Subprocess output could not be persisted.",
                severity="error",
                path="output_callback",
                source="process_runner",
            ),
            truncated,
        )

    @classmethod
    def _failure_with_output_limit(
        cls,
        primary_issue: Issue,
        truncated: Mapping[str, bool],
    ) -> Result[ProcessResult]:
        issues = [primary_issue]
        if any(truncated.values()):
            issues.append(cls._output_limit_issue())
        return Result.failure(issues)

    @staticmethod
    def _output_limit_issue() -> Issue:
        return Issue(
            code="PROCESS_RUNNER_OUTPUT_LIMIT_EXCEEDED",
            message="Subprocess output exceeded the size limit and was truncated.",
            severity="warning",
            path="command_spec",
            source="process_runner",
        )
