"""Server-owned cleanup for workflow containers on one local Docker daemon."""

from __future__ import annotations

from dataclasses import dataclass, field
import os
from pathlib import Path
import re
import selectors
import stat
import subprocess
import time

from encode_pipeline.platform.managed_containers import (
    MANAGED_CONTAINER_SCOPE_LABEL,
    managed_container_endpoint_identity,
)
from encode_pipeline.platform.results import Issue, Result


_CONTAINER_ID = re.compile(r"^[0-9a-f]{64}$")
_SCOPE = re.compile(r"^[0-9a-f]{64}$")
_MAX_CONTAINERS = 128
_MAX_OUTPUT_BYTES = 64 * 1024
_EMPTY_OBSERVATION_INTERVAL_SECONDS = 0.1


@dataclass(frozen=True)
class ManagedContainerCleaner:
    """Clean containers under a bounded, server-owned daemon quiescence policy.

    ``cleanup_timeout_seconds`` is the continuous late-registration observation
    budget after the process tree has been reaped. Activity restarts that budget,
    up to ``maximum_cleanup_timeout_seconds``. A separate final cleanup/audit
    budget handles registrations at the observation boundary. Docker exposes no
    in-flight-create fence; a request materializing after both server-owned
    budgets is an external-runtime residual outside this bounded contract.
    """

    executable: Path = field(repr=False)
    unix_socket: Path = field(default=Path("/var/run/docker.sock"), repr=False)
    cleanup_timeout_seconds: float = 2.0
    maximum_cleanup_timeout_seconds: float = 20.0
    final_cleanup_timeout_seconds: float = 5.0
    stop_timeout_seconds: int = 2
    _executable_identity: tuple[int, ...] = field(init=False, repr=False)
    _socket_identity: tuple[int, ...] = field(init=False, repr=False)

    def __post_init__(self) -> None:
        for name in ("executable", "unix_socket"):
            value = getattr(self, name)
            if not isinstance(value, Path) or not value.is_absolute():
                raise ValueError(f"{name} must be an absolute pathlib.Path")
            rendered = str(value)
            if rendered != str(Path(rendered)) or any(
                character in rendered for character in ("\x00", "\n", "\r")
            ):
                raise ValueError(f"{name} must be a canonical absolute path")
        if self.executable.name != "docker" or ":" in str(self.executable.parent):
            raise ValueError("executable must be a fixed Docker CLI path")
        for name, upper_bound in (
            ("cleanup_timeout_seconds", 60),
            ("maximum_cleanup_timeout_seconds", 120),
            ("final_cleanup_timeout_seconds", 60),
        ):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                raise ValueError(f"{name} must be a number")
            if not 0 < float(value) <= upper_bound:
                raise ValueError(f"{name} must be between 0 and {upper_bound}")
            object.__setattr__(self, name, float(value))
        if self.maximum_cleanup_timeout_seconds < self.cleanup_timeout_seconds:
            raise ValueError(
                "maximum_cleanup_timeout_seconds must cover cleanup_timeout_seconds"
            )
        if (
            isinstance(self.stop_timeout_seconds, bool)
            or not isinstance(self.stop_timeout_seconds, int)
            or not 1 <= self.stop_timeout_seconds <= 30
        ):
            raise ValueError("stop_timeout_seconds must be between 1 and 30")
        executable_identity, socket_identity = self._endpoint_identities()
        object.__setattr__(self, "_executable_identity", executable_identity)
        object.__setattr__(self, "_socket_identity", socket_identity)

    def cleanup(self, scope: str) -> Result[None]:
        """Observe, clean, and finally audit one bounded container scope."""
        if not isinstance(scope, str) or _SCOPE.fullmatch(scope) is None:
            return self._failure("MANAGED_CONTAINER_SCOPE_INVALID")
        started_at = time.monotonic()
        hard_deadline = started_at + self.maximum_cleanup_timeout_seconds
        quiet_deadline = started_at + self.cleanup_timeout_seconds
        try:
            while True:
                phase_deadline = min(quiet_deadline, hard_deadline)
                all_containers = self._list(
                    scope,
                    all_containers=True,
                    deadline=phase_deadline,
                )
                if all_containers:
                    self._clean_observed(
                        scope,
                        all_containers,
                        deadline=phase_deadline,
                    )
                    now = time.monotonic()
                    quiet_deadline = now + self.cleanup_timeout_seconds
                    if now >= hard_deadline:
                        break
                    continue
                if time.monotonic() >= phase_deadline:
                    break
                self._wait_for_next_observation(phase_deadline)
        except subprocess.TimeoutExpired:
            # The independent final phase gets a fresh budget and fresh daemon
            # observations even when an observation-phase CLI call timed out.
            pass
        except (OSError, subprocess.SubprocessError, ValueError):
            return self._failure("MANAGED_CONTAINER_CLEANUP_FAILED")

        final_deadline = time.monotonic() + self.final_cleanup_timeout_seconds
        try:
            while True:
                all_containers = self._list(
                    scope,
                    all_containers=True,
                    deadline=final_deadline,
                )
                if not all_containers:
                    # A distinct final audit is deliberately inside the fresh
                    # budget rather than reusing the last observation.
                    all_containers = self._list(
                        scope,
                        all_containers=True,
                        deadline=final_deadline,
                    )
                    if not all_containers:
                        return Result.success(None)
                self._clean_observed(
                    scope,
                    all_containers,
                    deadline=final_deadline,
                )
        except (OSError, subprocess.SubprocessError, ValueError):
            return self._failure("MANAGED_CONTAINER_CLEANUP_FAILED")

    def _clean_observed(
        self,
        scope: str,
        all_containers: tuple[str, ...],
        *,
        deadline: float,
    ) -> None:
        self._command(
            (
                "stop",
                "--time",
                str(self.stop_timeout_seconds),
                *all_containers,
            ),
            deadline=deadline,
        )
        running = self._list(
            scope,
            all_containers=False,
            deadline=deadline,
        )
        if running:
            self._command(("kill", *running), deadline=deadline)

        remaining = self._list(
            scope,
            all_containers=True,
            deadline=deadline,
        )
        if remaining:
            self._command(("rm", "--force", *remaining), deadline=deadline)

    @staticmethod
    def _wait_for_next_observation(deadline: float) -> None:
        remaining = deadline - time.monotonic()
        if remaining <= 0:
            raise subprocess.TimeoutExpired("docker", 0)
        time.sleep(min(_EMPTY_OBSERVATION_INTERVAL_SECONDS, remaining))

    @property
    def local_docker_host(self) -> str:
        """Return the server-owned local daemon URL injected into the workflow."""
        return f"unix://{self.unix_socket}"

    @property
    def endpoint_identity(self) -> str:
        """Return the private path identity required by a scoped command."""
        return managed_container_endpoint_identity(self.executable, self.unix_socket)

    def verify_endpoint(self) -> Result[None]:
        """Fail closed unless the captured executable/socket identities persist."""
        try:
            self._require_stable_endpoint()
        except (OSError, ValueError):
            return self._failure("MANAGED_CONTAINER_ENDPOINT_CHANGED")
        return Result.success(None)

    def _list(
        self,
        scope: str,
        *,
        all_containers: bool,
        deadline: float,
    ) -> tuple[str, ...]:
        argv = ["ps"]
        if all_containers:
            argv.append("--all")
        argv.extend(
            (
                "--quiet",
                "--no-trunc",
                "--filter",
                f"label={MANAGED_CONTAINER_SCOPE_LABEL}={scope}",
            )
        )
        stdout = self._command(tuple(argv), deadline=deadline)
        lines = tuple(line.strip() for line in stdout.splitlines() if line.strip())
        if (
            len(lines) > _MAX_CONTAINERS
            or len(lines) != len(set(lines))
            or any(_CONTAINER_ID.fullmatch(line) is None for line in lines)
        ):
            raise ValueError("docker returned invalid container identities")
        return lines

    def _command(self, argv: tuple[str, ...], *, deadline: float) -> str:
        self._require_stable_endpoint()
        if deadline <= time.monotonic():
            raise subprocess.TimeoutExpired(str(self.executable), 0)
        process = subprocess.Popen(
            (
                str(self.executable),
                "--host",
                f"unix://{self.unix_socket}",
                *argv,
            ),
            shell=False,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            env={"PATH": "/usr/bin:/bin"},
        )
        assert process.stdout is not None
        output = bytearray()
        selector = selectors.DefaultSelector()
        descriptor = process.stdout.fileno()
        try:
            os.set_blocking(descriptor, False)
            selector.register(descriptor, selectors.EVENT_READ)
            while selector.get_map() or process.poll() is None:
                remaining = deadline - time.monotonic()
                if remaining <= 0:
                    raise subprocess.TimeoutExpired(str(self.executable), 0)
                for _key, _events in selector.select(timeout=min(remaining, 0.05)):
                    remaining_bytes = _MAX_OUTPUT_BYTES - len(output) + 1
                    chunk = os.read(descriptor, min(8192, remaining_bytes))
                    if not chunk:
                        selector.unregister(descriptor)
                        continue
                    output.extend(chunk)
                    if len(output) > _MAX_OUTPUT_BYTES:
                        raise OSError("docker cleanup output exceeded the bound")
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                raise subprocess.TimeoutExpired(str(self.executable), 0)
            return_code = process.wait(timeout=remaining)
            if return_code != 0:
                raise OSError("docker cleanup command failed")
            return bytes(output).decode("utf-8", errors="replace")
        except BaseException:
            if process.poll() is None:
                process.kill()
            try:
                process.wait(timeout=0.1)
            except subprocess.TimeoutExpired:
                pass
            raise
        finally:
            selector.close()
            process.stdout.close()

    def _require_stable_endpoint(self) -> None:
        executable_identity, socket_identity = self._endpoint_identities()
        if (
            executable_identity != self._executable_identity
            or socket_identity != self._socket_identity
        ):
            raise OSError("managed Docker endpoints changed")

    def _endpoint_identities(self) -> tuple[tuple[int, ...], tuple[int, ...]]:
        executable = os.lstat(self.executable)
        socket_info = os.lstat(self.unix_socket)
        if (
            not stat.S_ISREG(executable.st_mode)
            or executable.st_nlink != 1
            or executable.st_mode & 0o111 == 0
            or executable.st_mode & 0o022 != 0
            or executable.st_uid not in {0, os.geteuid()}
            or not stat.S_ISSOCK(socket_info.st_mode)
            or socket_info.st_nlink != 1
        ):
            raise ValueError("managed Docker endpoints are invalid")
        return (
            (
                executable.st_dev,
                executable.st_ino,
                executable.st_mode,
                executable.st_size,
                executable.st_mtime_ns,
                executable.st_ctime_ns,
            ),
            (socket_info.st_dev, socket_info.st_ino, socket_info.st_mode),
        )

    @staticmethod
    def _failure(reason_code: str) -> Result[None]:
        return Result.failure(
            [
                Issue(
                    code=reason_code,
                    message="Managed workflow containers could not be cleaned up.",
                    severity="error",
                    path="execution.cleanup",
                    source="managed_container_cleaner",
                )
            ]
        )
