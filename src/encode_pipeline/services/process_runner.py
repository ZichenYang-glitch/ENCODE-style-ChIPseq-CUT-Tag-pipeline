"""Safe subprocess boundary for executing already-built CommandSpec instances."""

from __future__ import annotations

import os
import subprocess
from dataclasses import dataclass
from typing import TYPE_CHECKING

from encode_pipeline.platform.results import Issue, Result

if TYPE_CHECKING:
    from encode_pipeline.platform.adapters import CommandSpec


@dataclass(frozen=True)
class ProcessResult:
    """The outcome of a successfully-run subprocess."""

    exit_code: int
    stdout: str
    stderr: str


def _decode_and_truncate(data: bytes, max_bytes: int) -> tuple[str, bool]:
    """Decode *data* as UTF-8, truncating to *max_bytes* at the byte level.

    Returns ``(decoded_string, was_truncated)``.
    """
    truncated = len(data) > max_bytes
    if truncated:
        data = data[:max_bytes]
    return data.decode("utf-8", errors="replace"), truncated


class ProcessRunner:
    """Execute a CommandSpec via subprocess.run with safety guardrails.

    Enforces an executable allowlist, output size limits, a wall-clock
    timeout, and inherits ``os.environ`` with optional overrides from the
    command spec.  Returns ``Result[ProcessResult]`` — nonzero exits are
    successes with a warning issue; infrastructure failures (timeout,
    executable not found, OSError) are error failures.
    """

    def __init__(
        self,
        *,
        allowed_executables: tuple[str, ...] = ("snakemake",),
        timeout_seconds: float = 300.0,
        max_output_bytes: int = 10_000_000,
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
        if isinstance(timeout_seconds, bool) or not isinstance(timeout_seconds, (int, float)):
            raise ValueError("timeout_seconds must be a number")
        if timeout_seconds <= 0:
            raise ValueError("timeout_seconds must be positive")
        if isinstance(max_output_bytes, bool) or not isinstance(max_output_bytes, int):
            raise ValueError("max_output_bytes must be an int")
        if max_output_bytes <= 0:
            raise ValueError("max_output_bytes must be positive")

        self._allowed_executables = allowed_executables
        self._timeout_seconds = float(timeout_seconds)
        self._max_output_bytes = max_output_bytes

    def run(self, spec: "CommandSpec") -> "Result[ProcessResult]":
        """Execute *spec* and return a structured result."""
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

        env = dict(os.environ)
        env.update(spec.env)

        try:
            completed = subprocess.run(
                spec.argv,
                cwd=spec.cwd,
                env=env,
                capture_output=True,
                timeout=self._timeout_seconds,
            )
        except subprocess.TimeoutExpired:
            return Result.failure(
                [
                    Issue(
                        code="PROCESS_RUNNER_TIMEOUT",
                        message="Subprocess timed out.",
                        severity="error",
                        path="command_spec",
                        source="process_runner",
                        context={"timeout_seconds": self._timeout_seconds},
                    )
                ]
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

        stdout_str, stdout_truncated = _decode_and_truncate(
            completed.stdout, self._max_output_bytes
        )
        stderr_str, stderr_truncated = _decode_and_truncate(
            completed.stderr, self._max_output_bytes
        )

        issues: list[Issue] = []
        if stdout_truncated or stderr_truncated:
            issues.append(
                Issue(
                    code="PROCESS_RUNNER_OUTPUT_LIMIT_EXCEEDED",
                    message="Subprocess output exceeded the size limit and was truncated.",
                    severity="warning",
                    path="command_spec",
                    source="process_runner",
                )
            )

        exit_code = completed.returncode
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
                stdout=stdout_str,
                stderr=stderr_str,
            ),
            issues=issues,
        )
