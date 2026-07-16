"""Safe subprocess boundary for executing already-built CommandSpec instances."""

from __future__ import annotations

import codecs
import ctypes
import os
import selectors
import signal
import subprocess
import sys
import threading
import time
from collections.abc import Callable, Mapping
from dataclasses import dataclass
from typing import TYPE_CHECKING

from encode_pipeline.platform.results import Issue, Result

if TYPE_CHECKING:
    from encode_pipeline.platform.adapters import CommandSpec
    from encode_pipeline.services.managed_containers import ManagedContainerCleaner


OutputCallback = Callable[[str, tuple[str, ...]], None]


@dataclass(frozen=True)
class ProcessResult:
    """The outcome of a successfully-run subprocess."""

    exit_code: int
    stdout: str
    stderr: str


# Only runtime variables needed to locate and run the scientific toolchain are
# inherited. In particular, worker/API configuration and provider credentials
# must never become part of a workflow process (or one of its task children).
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
        "DOCKER_HOST",
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
_OUTPUT_REDACTION = "[REDACTED]"
_MAX_STREAM_REDACTION_PATTERNS = 8_192
_MAX_STREAM_REDACTION_PATTERN_LENGTH = 4096
_MAX_STREAM_REDACTION_TOTAL_CHARACTERS = 16 * 1024 * 1024
_MAX_STREAM_REDACTION_TOTAL_UTF8_BYTES = 32 * 1024 * 1024
_LITERAL_ANCHOR_LENGTH = 8
_MIN_LITERAL_MATCH_WORK = 1_000_000
_LITERAL_MATCH_WORK_PER_CHARACTER = 64
_TERMINATION_GRACE_SECONDS = 2.0
_TERMINATION_POLL_SECONDS = 0.02
_FORCED_CLEANUP_SECONDS = 2.0
_SUBREAPER_LOCK = threading.RLock()
_PR_SET_CHILD_SUBREAPER = 36
_PR_GET_CHILD_SUBREAPER = 37


class ProcessRunnerCleanupError(BaseException):
    """An asynchronous abort could not prove that owned resources stopped."""


class _SubreaperUnavailable(OSError):
    """The strict Linux child-reaping boundary could not be established."""


class _RedactionWorkExceeded(RuntimeError):
    """A bounded literal search encountered adversarial candidate density."""


class _RedactionBoundExceeded(RuntimeError):
    """A literal-redaction contract exceeded its explicit resource bounds."""


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


def _redaction_values(spec: "CommandSpec") -> tuple[str, ...]:
    """Return longest-first private values that may be echoed by a child.

    Adapters can explicitly name submitted resource paths. The runner also
    treats its controlled cwd plus absolute argv/environment paths as private
    execution details. Values are never copied into public Issues.
    """
    values = set(getattr(spec, "redaction_values", ()))
    if spec.cwd:
        values.add(spec.cwd)
    values.update(value for value in spec.argv if value.startswith("/"))
    values.update(value for value in spec.env.values() if value.startswith("/"))
    return tuple(sorted((value for value in values if value), key=len, reverse=True))


class _BoundedLiteralMatcher:
    """Linear-space fixed-anchor matcher shared by every output path.

    Long literals select the rarest of five fixed-size anchors. Scanning is
    linear in output plus verified candidates, while an explicit weighted work
    limit fails closed for adversarial shared anchors. Short literals use at
    most seven exact-length scans. The index stores references and bounded
    anchors, never a second copy of the aggregate literal text.
    """

    def __init__(self, private_values: tuple[str, ...]) -> None:
        self.max_length = max((len(value) for value in private_values), default=0)
        self.max_utf8_length = max(
            (len(value.encode("utf-8", errors="replace")) for value in private_values),
            default=0,
        )
        self.newline_is_boundary = all("\n" not in value for value in private_values)
        short_by_length: dict[int, set[str]] = {}
        long_patterns: list[tuple[str, tuple[int, ...]]] = []
        anchor_frequency: dict[str, int] = {}
        for value in private_values:
            if len(value) < _LITERAL_ANCHOR_LENGTH:
                short_by_length.setdefault(len(value), set()).add(value)
                continue
            last = len(value) - _LITERAL_ANCHOR_LENGTH
            offsets = tuple(sorted({0, last // 4, last // 2, (last * 3) // 4, last}))
            long_patterns.append((value, offsets))
            for offset in offsets:
                anchor = value[offset : offset + _LITERAL_ANCHOR_LENGTH]
                anchor_frequency[anchor] = anchor_frequency.get(anchor, 0) + 1

        anchor_index: dict[str, list[tuple[str, int]]] = {}
        for value, offsets in long_patterns:
            offset = min(
                offsets,
                key=lambda candidate: (
                    anchor_frequency[
                        value[candidate : candidate + _LITERAL_ANCHOR_LENGTH]
                    ],
                    candidate,
                ),
            )
            anchor = value[offset : offset + _LITERAL_ANCHOR_LENGTH]
            anchor_index.setdefault(anchor, []).append((value, offset))
        self._short_by_length = tuple(
            (length, frozenset(values))
            for length, values in sorted(short_by_length.items())
        )
        self._anchor_index = {
            anchor: tuple(candidates) for anchor, candidates in anchor_index.items()
        }

    def mask(self, value: str) -> bytearray:
        mask = bytearray(len(value))
        for length, patterns in self._short_by_length:
            for start in range(0, len(value) - length + 1):
                if value[start : start + length] in patterns:
                    mask[start : start + length] = b"\x01" * length

        maximum_work = max(
            _MIN_LITERAL_MATCH_WORK,
            len(value) * _LITERAL_MATCH_WORK_PER_CHARACTER,
        )
        work = 0
        for anchor_start in range(
            0,
            len(value) - _LITERAL_ANCHOR_LENGTH + 1,
        ):
            anchor = value[anchor_start : anchor_start + _LITERAL_ANCHOR_LENGTH]
            for pattern, offset in self._anchor_index.get(anchor, ()):
                start = anchor_start - offset
                if start < 0 or start + len(pattern) > len(value):
                    continue
                work += len(pattern)
                if work > maximum_work:
                    raise _RedactionWorkExceeded
                if value.startswith(pattern, start):
                    mask[start : start + len(pattern)] = b"\x01" * len(pattern)
        return mask

    def redact(self, value: str) -> str:
        if not value or self.max_length == 0:
            return value
        return _render_redacted(value, self.mask(value))


def _build_bounded_literal_matcher(
    private_values: tuple[str, ...],
) -> _BoundedLiteralMatcher:
    """Build one matcher only after enforcing the shared redaction contract."""
    values = tuple(
        sorted(
            {value for value in private_values if value},
            key=lambda value: (-len(value), value),
        )
    )
    if (
        len(values) > _MAX_STREAM_REDACTION_PATTERNS
        or any(len(value) > _MAX_STREAM_REDACTION_PATTERN_LENGTH for value in values)
        or sum(len(value) for value in values) > _MAX_STREAM_REDACTION_TOTAL_CHARACTERS
        or sum(len(value.encode("utf-8", errors="replace")) for value in values)
        > _MAX_STREAM_REDACTION_TOTAL_UTF8_BYTES
    ):
        raise _RedactionBoundExceeded
    try:
        return _BoundedLiteralMatcher(values)
    except (MemoryError, OverflowError) as exc:
        raise _RedactionBoundExceeded from exc


def redact_bounded_literals(value: str, private_values: tuple[str, ...]) -> str:
    """Redact literals with the same size and candidate-work bounds as streams."""
    return _build_bounded_literal_matcher(private_values).redact(value)


def _redact_captured_output(
    value: bytes,
    matcher: _BoundedLiteralMatcher,
    *,
    truncated: bool,
) -> str:
    if truncated and value:
        # Capture is byte-bounded and may end inside UTF-8. Conservatively hide
        # one whole maximum-literal tail, then redact complete literals in the
        # safe prefix. At most the already-bounded tail is sacrificed.
        tail_length = min(len(value), matcher.max_utf8_length)
        safe = _decode(value[:-tail_length])
        return matcher.redact(safe) + _OUTPUT_REDACTION
    return matcher.redact(_decode(value))


def _render_redacted(value: str, mask: bytearray) -> str:
    segments: list[str] = []
    cursor = 0
    while cursor < len(value):
        redacted = bool(mask[cursor])
        end = cursor + 1
        while end < len(value) and bool(mask[end]) == redacted:
            end += 1
        segments.append(_OUTPUT_REDACTION if redacted else value[cursor:end])
        cursor = end
    return "".join(segments)


def _linux_process_table() -> dict[int, tuple[int, int]]:
    """Return pid -> (ppid, starttime) for the current Linux process table."""
    table: dict[int, tuple[int, int]] = {}
    try:
        entries = os.scandir("/proc")
    except OSError:
        return table
    with entries:
        for entry in entries:
            if not entry.name.isdigit():
                continue
            try:
                stat = _ProcessStat.read(entry.name)
            except (OSError, ValueError):
                continue
            if stat is not None:
                table[int(entry.name)] = stat
    return table


class _ProcessStat:
    """Small parser isolated for defensive `/proc/<pid>/stat` reads."""

    @staticmethod
    def read(pid: str) -> tuple[int, int] | None:
        with open(f"/proc/{pid}/stat", encoding="utf-8") as handle:
            value = handle.read(4096)
        closing = value.rfind(")")
        if closing < 0:
            return None
        fields = value[closing + 2 :].split()
        if len(fields) < 20:
            return None
        return int(fields[1]), int(fields[19])


def _linux_descendants(root_pid: int) -> tuple[tuple[int, int], ...]:
    table = _linux_process_table()
    frontier = [root_pid]
    descendants: list[tuple[int, int]] = []
    while frontier:
        parent = frontier.pop(0)
        children = sorted(
            pid for pid, (ppid, _starttime) in table.items() if ppid == parent
        )
        for pid in children:
            descendants.append((pid, table[pid][1]))
            frontier.append(pid)
    return tuple(descendants)


def _linux_process_matches(pid: int, starttime: int) -> bool:
    try:
        current = _ProcessStat.read(str(pid))
    except (OSError, ValueError):
        return False
    return current is not None and current[1] == starttime


class _LinuxSubreaper:
    """Temporarily adopt and track workflow descendants across reparenting."""

    def __init__(self, *, strict: bool) -> None:
        self._strict = strict
        self._enabled = False
        self._previous = 0
        self._libc = None
        self._root: tuple[int, int] | None = None
        self._spawned_pid: int | None = None
        self._tracked: dict[int, int] = {}
        self._baseline_children: set[tuple[int, int]] = set()
        self._signal_failed = False

    def __enter__(self) -> "_LinuxSubreaper":
        _SUBREAPER_LOCK.acquire()
        try:
            if not sys.platform.startswith("linux"):
                if self._strict:
                    raise _SubreaperUnavailable("Linux subreaper is unavailable")
                return self
            libc = ctypes.CDLL(None, use_errno=True)
            previous = ctypes.c_int()
            if (
                libc.prctl(
                    _PR_GET_CHILD_SUBREAPER,
                    ctypes.byref(previous),
                    0,
                    0,
                    0,
                )
                != 0
            ):
                raise _SubreaperUnavailable("could not read subreaper state")
            if libc.prctl(_PR_SET_CHILD_SUBREAPER, 1, 0, 0, 0) != 0:
                raise _SubreaperUnavailable("could not enable subreaper state")
            self._libc = libc
            self._previous = previous.value
            self._enabled = True
            table = _linux_process_table()
            current_pid = os.getpid()
            self._baseline_children = {
                (pid, starttime)
                for pid, (ppid, starttime) in table.items()
                if ppid == current_pid
            }
            return self
        except _SubreaperUnavailable:
            if not self._strict:
                self._enabled = False
                return self
            _SUBREAPER_LOCK.release()
            raise
        except BaseException:
            _SUBREAPER_LOCK.release()
            raise

    def __exit__(self, exc_type, _exc_value, _traceback) -> None:
        restore_failed = False
        try:
            self._reap_tracked()
            if self._enabled and self._libc is not None:
                restore_failed = (
                    self._libc.prctl(
                        _PR_SET_CHILD_SUBREAPER,
                        self._previous,
                        0,
                        0,
                        0,
                    )
                    != 0
                )
        finally:
            _SUBREAPER_LOCK.release()
        if restore_failed:
            if exc_type is not None:
                raise ProcessRunnerCleanupError(
                    "Workflow cleanup state could not be restored."
                ) from None
            raise _SubreaperUnavailable("could not restore subreaper state")

    def bind(self, root_pid: int) -> None:
        self._spawned_pid = root_pid
        try:
            root = _ProcessStat.read(str(root_pid))
        except (OSError, ValueError) as exc:
            if not self._strict:
                self._root = None
                return
            raise _SubreaperUnavailable("could not identify workflow root") from exc
        if root is None:
            if not self._strict:
                self._root = None
                return
            raise _SubreaperUnavailable("could not identify workflow root")
        self._root = (root_pid, root[1])

    def register_spawn(self, root_pid: int) -> None:
        """Record the Popen PID before a fast root can exit and reparent children."""
        self._spawned_pid = root_pid

    @property
    def has_unbound_spawn(self) -> bool:
        """Return whether Popen succeeded but its root identity was not bound."""
        return self._spawned_pid is not None and self._root is None

    def observe(self) -> tuple[tuple[int, int], ...]:
        """Remember descendants plus newly adopted children by pid/starttime."""
        table = _linux_process_table()
        root_pid = self._spawned_pid if self._root is None else self._root[0]
        if self._root is not None:
            root_starttime: int | None = self._root[1]
        elif root_pid is not None and (root_stat := table.get(root_pid)) is not None:
            root_starttime = root_stat[1]
        else:
            root_starttime = None
        identities = tuple(self._tracked.items())
        if root_pid is not None and root_starttime is not None:
            identities = ((root_pid, root_starttime), *identities)
        live_roots = {
            pid
            for pid, starttime in identities
            if table.get(pid, (None, None))[1] == starttime
        }
        frontier = list(live_roots)
        while frontier:
            parent = frontier.pop()
            for pid, (ppid, starttime) in table.items():
                if ppid != parent or pid in self._tracked or pid == root_pid:
                    continue
                self._tracked[pid] = starttime
                frontier.append(pid)

        if self._enabled:
            current_pid = os.getpid()
            for pid, (ppid, starttime) in table.items():
                identity = (pid, starttime)
                if (
                    ppid == current_pid
                    and pid != root_pid
                    and identity not in self._baseline_children
                    and (root_starttime is None or starttime >= root_starttime)
                ):
                    self._tracked[pid] = starttime

        return tuple(
            sorted(
                (
                    (pid, starttime)
                    for pid, starttime in self._tracked.items()
                    if table.get(pid, (None, None))[1] == starttime
                ),
                reverse=True,
            )
        )

    def terminate(self, process: subprocess.Popen[bytes]) -> bool:
        """TERM, then repeatedly discover/KILL/reap the bounded owned tree."""
        if self._root is None:
            try:
                if process.poll() is None:
                    process.terminate()
                    process.wait(timeout=_TERMINATION_GRACE_SECONDS)
            except ProcessLookupError:
                pass
            except subprocess.TimeoutExpired:
                process.kill()
                try:
                    process.wait(timeout=_FORCED_CLEANUP_SECONDS)
                except subprocess.TimeoutExpired:
                    return False
            return process.poll() is not None and self.terminate_adopted()
        term_deadline = time.monotonic() + _TERMINATION_GRACE_SECONDS
        term_signalled: set[tuple[int, int]] = set()
        while time.monotonic() < term_deadline:
            descendants = self.observe()
            identities = (
                (process.pid, self._root_starttime()),
                *descendants,
            )
            for pid, starttime in identities:
                if (pid, starttime) in term_signalled:
                    continue
                self._signal(pid, starttime, signal.SIGTERM)
                term_signalled.add((pid, starttime))
            process.poll()
            self._reap_tracked()
            if process.poll() is not None and not self.observe():
                return not self._signal_failed
            time.sleep(_TERMINATION_POLL_SECONDS)

        kill_deadline = time.monotonic() + _FORCED_CLEANUP_SECONDS
        stable_empty_rounds = 0
        while time.monotonic() < kill_deadline:
            descendants = self.observe()
            for pid, starttime in descendants:
                self._signal(pid, starttime, signal.SIGKILL)
            self._signal(process.pid, self._root_starttime(), signal.SIGKILL)
            process.poll()
            self._reap_tracked()
            live = self.observe()
            if process.poll() is not None and not live:
                stable_empty_rounds += 1
                if stable_empty_rounds >= 3:
                    return not self._signal_failed
            else:
                stable_empty_rounds = 0
            time.sleep(_TERMINATION_POLL_SECONDS)

        process.poll()
        self._reap_tracked()
        return (
            process.poll() is not None
            and not self.observe()
            and not self._signal_failed
        )

    def terminate_adopted(self) -> bool:
        """Clean children adopted after Popen even when the root identity vanished."""
        term_deadline = time.monotonic() + _TERMINATION_GRACE_SECONDS
        term_signalled: set[tuple[int, int]] = set()
        stable_empty_rounds = 0
        while time.monotonic() < term_deadline:
            descendants = self.observe()
            if not descendants:
                stable_empty_rounds += 1
                if stable_empty_rounds >= 3:
                    return not self._signal_failed
            else:
                stable_empty_rounds = 0
                for pid, starttime in descendants:
                    if (pid, starttime) not in term_signalled:
                        self._signal(pid, starttime, signal.SIGTERM)
                        term_signalled.add((pid, starttime))
            self._reap_tracked()
            time.sleep(_TERMINATION_POLL_SECONDS)

        kill_deadline = time.monotonic() + _FORCED_CLEANUP_SECONDS
        stable_empty_rounds = 0
        while time.monotonic() < kill_deadline:
            descendants = self.observe()
            for pid, starttime in descendants:
                self._signal(pid, starttime, signal.SIGKILL)
            self._reap_tracked()
            if not self.observe():
                stable_empty_rounds += 1
                if stable_empty_rounds >= 3:
                    return not self._signal_failed
            else:
                stable_empty_rounds = 0
            time.sleep(_TERMINATION_POLL_SECONDS)
        self._reap_tracked()
        return not self.observe() and not self._signal_failed

    def _root_starttime(self) -> int:
        return -1 if self._root is None else self._root[1]

    def _signal(self, pid: int, starttime: int, signum: signal.Signals) -> None:
        if starttime < 0 or not _linux_process_matches(pid, starttime):
            return
        try:
            os.kill(pid, signum)
        except ProcessLookupError:
            return
        except OSError:
            self._signal_failed = True

    def _reap_tracked(self) -> None:
        for pid, starttime in tuple(self._tracked.items()):
            if _linux_process_matches(pid, starttime):
                try:
                    os.waitpid(pid, os.WNOHANG)
                except (ChildProcessError, OSError):
                    continue
            else:
                self._tracked.pop(pid, None)


class _StreamingRedactor:
    """Redact literals while retaining only a bounded cross-chunk suffix."""

    def __init__(self, matcher: _BoundedLiteralMatcher) -> None:
        self._matcher = matcher
        self._pending = ""
        self._finished = False

    def feed(self, value: str) -> str:
        if self._finished:
            raise RuntimeError("stream redactor is already finished")
        if self._matcher.max_length == 0:
            return value
        self._pending += value
        cut = max(0, len(self._pending) - self._matcher.max_length + 1)
        if (
            self._matcher.newline_is_boundary
            and (newline := self._pending.rfind("\n")) >= 0
        ):
            cut = max(cut, newline + 1)
        if cut == 0:
            return ""

        mask = self._matcher.mask(self._pending)
        while cut < len(self._pending) and mask[cut - 1] and mask[cut]:
            cut += 1

        safe = self._pending[:cut]
        self._pending = self._pending[cut:]
        return _render_redacted(safe, mask[:cut])

    def finish(self, *, truncated: bool = False) -> str:
        if self._finished:
            return ""
        self._finished = True
        pending = self._pending
        self._pending = ""
        if truncated and pending:
            return _OUTPUT_REDACTION
        return self._matcher.redact(pending)


class _OutputLineDecoder:
    """Incrementally decode one stream and emit bounded, complete text chunks."""

    def __init__(self, matcher: _BoundedLiteralMatcher) -> None:
        decoder_factory = codecs.getincrementaldecoder("utf-8")
        self._decoder = decoder_factory(errors="replace")
        self._redactor = _StreamingRedactor(matcher)
        self._pending = ""
        self._finished = False

    def feed(self, data: bytes) -> tuple[str, ...]:
        if self._finished:
            raise RuntimeError("output decoder is already finished")
        decoded = self._decoder.decode(data, final=False)
        self._pending += self._redactor.feed(decoded)
        return self._take_chunks(final=False)

    def finish(self, *, truncated: bool = False) -> tuple[str, ...]:
        if self._finished:
            return ()
        self._finished = True
        decoded = "" if truncated else self._decoder.decode(b"", final=True)
        self._pending += self._redactor.feed(decoded)
        self._pending += self._redactor.finish(truncated=truncated)
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
        managed_container_cleaner: "ManagedContainerCleaner | None" = None,
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
        if managed_container_cleaner is not None:
            from encode_pipeline.services.managed_containers import (
                ManagedContainerCleaner,
            )

            if not isinstance(managed_container_cleaner, ManagedContainerCleaner):
                raise ValueError(
                    "managed_container_cleaner must be a ManagedContainerCleaner or None"
                )

        self._allowed_executables = allowed_executables
        self._timeout_seconds = float(timeout_seconds)
        self._max_output_bytes = max_output_bytes
        self._passthrough_exceptions = passthrough_exceptions
        self._managed_container_cleaner = managed_container_cleaner

    def run(
        self,
        spec: "CommandSpec",
        *,
        output_callback: OutputCallback | None = None,
    ) -> "Result[ProcessResult]":
        """Run one command inside a temporary descendant-cleanup boundary."""
        from encode_pipeline.platform.adapters import CommandSpec

        scope = spec.managed_container_scope if isinstance(spec, CommandSpec) else None
        if scope is not None and self._managed_container_cleaner is None:
            return Result.failure(
                [
                    Issue(
                        code="PROCESS_RUNNER_CONTAINER_CLEANER_UNAVAILABLE",
                        message="Managed container cleanup is unavailable.",
                        severity="error",
                        path="command_spec",
                        source="process_runner",
                    )
                ]
            )
        if (
            scope is not None
            and spec.managed_container_endpoint_identity
            != self._managed_container_cleaner.endpoint_identity
        ):
            return Result.failure(
                [
                    Issue(
                        code="PROCESS_RUNNER_CONTAINER_ENDPOINT_MISMATCH",
                        message="Managed container endpoint identity does not match.",
                        severity="error",
                        path="command_spec",
                        source="process_runner",
                    )
                ]
            )

        try:
            with _LinuxSubreaper(strict=scope is not None) as subreaper:
                self._active_process: subprocess.Popen[bytes] | None = None
                self._active_subreaper = subreaper
                self._tree_cleanup_failed = False
                try:
                    try:
                        result = self._run_inner(
                            spec,
                            output_callback=output_callback,
                        )
                    except BaseException:
                        cleanup_issues = self._cleanup_after_abort(scope)
                        if cleanup_issues:
                            raise ProcessRunnerCleanupError(
                                "Workflow cleanup could not be confirmed."
                            ) from None
                        raise

                    cleanup_issues: list[Issue] = []
                    if self._active_process is not None and (
                        self._active_process.poll() is None or subreaper.observe()
                    ):
                        if not self._terminate(self._active_process):
                            self._tree_cleanup_failed = True
                    if self._tree_cleanup_failed:
                        cleanup_issues.append(self._tree_cleanup_issue())
                    if scope is not None:
                        cleanup_issues.extend(self._cleanup_managed_containers(scope))
                    if cleanup_issues:
                        combined = (
                            [*result.issues, *cleanup_issues]
                            if result.is_failure
                            else [*cleanup_issues, *result.issues]
                        )
                        return Result.failure(combined)
                    return result
                finally:
                    self._active_process = None
                    self._active_subreaper = None
        except _SubreaperUnavailable:
            return Result.failure(
                [
                    Issue(
                        code="PROCESS_RUNNER_SUBREAPER_UNAVAILABLE",
                        message="Managed process cleanup is unavailable.",
                        severity="error",
                        path="command_spec",
                        source="process_runner",
                    )
                ]
            )

    def _run_inner(
        self,
        spec: "CommandSpec",
        *,
        output_callback: OutputCallback | None = None,
    ) -> "Result[ProcessResult]":
        """Execute *spec* and return a structured, bounded result.

        When supplied, ``output_callback`` receives bounded line chunks as the
        child writes them. A callback failure terminates the workflow tree and is
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
        if spec.managed_container_scope is not None:
            cleaner = self._managed_container_cleaner
            if cleaner is None:  # guarded by run(), retained as a local invariant
                return Result.failure(
                    [
                        Issue(
                            code="PROCESS_RUNNER_CONTAINER_CLEANER_UNAVAILABLE",
                            message="Managed container cleanup is unavailable.",
                            severity="error",
                            path="command_spec",
                            source="process_runner",
                        )
                    ]
                )
            environment_result.value["DOCKER_HOST"] = cleaner.local_docker_host
            current_path = environment_result.value.get("PATH", "")
            environment_result.value["PATH"] = os.pathsep.join(
                value
                for value in (str(cleaner.executable.parent), current_path)
                if value
            )
        private_output_values = _redaction_values(spec)
        if spec.managed_container_scope is not None:
            assert self._managed_container_cleaner is not None
            cleaner_values = {
                str(self._managed_container_cleaner.executable),
                str(self._managed_container_cleaner.executable.parent),
                str(self._managed_container_cleaner.unix_socket),
                self._managed_container_cleaner.local_docker_host,
            }
            private_output_values = tuple(
                sorted(
                    {*private_output_values, *cleaner_values},
                    key=len,
                    reverse=True,
                )
            )
        try:
            private_output_matcher = _build_bounded_literal_matcher(
                private_output_values
            )
        except _RedactionBoundExceeded:
            return Result.failure(
                [
                    Issue(
                        code="PROCESS_RUNNER_REDACTION_BOUND_EXCEEDED",
                        message="Command output redaction exceeds the supported bound.",
                        severity="error",
                        path="command_spec",
                        source="process_runner",
                    )
                ]
            )

        if spec.managed_container_scope is not None:
            assert self._managed_container_cleaner is not None
            endpoint_result = self._managed_container_cleaner.verify_endpoint()
            if endpoint_result.is_failure:
                return Result.failure(
                    [
                        Issue(
                            code="PROCESS_RUNNER_CONTAINER_ENDPOINT_CHANGED",
                            message="Managed container endpoint identity changed.",
                            severity="error",
                            path="command_spec",
                            source="process_runner",
                        )
                    ]
                )

        try:
            process = subprocess.Popen(
                spec.argv,
                cwd=spec.cwd,
                env=environment_result.value,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            self._active_process = process
            assert self._active_subreaper is not None
            self._active_subreaper.register_spawn(process.pid)
            self._active_subreaper.bind(process.pid)
        except _SubreaperUnavailable:
            raise
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
            "stdout": _OutputLineDecoder(private_output_matcher),
            "stderr": _OutputLineDecoder(private_output_matcher),
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
                assert self._active_subreaper is not None
                self._active_subreaper.observe()
                remaining = self._timeout_seconds - (time.monotonic() - started_at)
                if remaining <= 0:
                    self._terminate(process)
                    if not drain_available_after_termination():
                        return self._callback_failure(truncated)
                    if output_callback is not None:
                        for stream_name in ("stdout", "stderr"):
                            chunks = decoders[stream_name].finish(
                                truncated=truncated[stream_name]
                            )
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
                            chunks = decoders[stream_name].finish(
                                truncated=truncated[stream_name]
                            )
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
        except _RedactionWorkExceeded:
            self._terminate(process)
            return self._redaction_work_failure(truncated)
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
            # a BaseException raised between reads, must not leave the spawned
            # workflow tree alive.
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

        try:
            stdout = _redact_captured_output(
                bytes(captures["stdout"]),
                private_output_matcher,
                truncated=truncated["stdout"],
            )
            stderr = _redact_captured_output(
                bytes(captures["stderr"]),
                private_output_matcher,
                truncated=truncated["stderr"],
            )
        except _RedactionWorkExceeded:
            return self._redaction_work_failure(truncated)

        return Result.success(
            ProcessResult(
                exit_code=exit_code,
                stdout=stdout,
                stderr=stderr,
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
            self._terminate(process)
            return False
        return True

    def _terminate(self, process: subprocess.Popen[bytes]) -> bool:
        subreaper = self._active_subreaper
        if subreaper is not None:
            cleaned = subreaper.terminate(process)
            self._tree_cleanup_failed = self._tree_cleanup_failed or not cleaned
            try:
                process.wait(timeout=0.1)
            except subprocess.TimeoutExpired:
                pass
            return cleaned

        descendants = _linux_descendants(process.pid)
        if process.poll() is None:
            try:
                process.terminate()
            except ProcessLookupError:
                pass
            deadline = time.monotonic() + _TERMINATION_GRACE_SECONDS
            while process.poll() is None and time.monotonic() < deadline:
                time.sleep(_TERMINATION_POLL_SECONDS)

        # Nextflow receives TERM first so its pinned runtime can stop task
        # wrappers and containers. Any descendants that survive that bounded
        # grace are killed by their original pid/start-time identity. This does
        # not create a new process group, so RQ's existing horse-group
        # cancellation continues to cover the entire workflow tree.
        for pid, starttime in reversed(descendants):
            if not _linux_process_matches(pid, starttime):
                continue
            try:
                os.kill(pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
        if process.poll() is None:
            process.kill()
        try:
            process.wait(timeout=5)
        except (
            subprocess.TimeoutExpired
        ):  # pragma: no cover - kill is immediate on POSIX
            return False
        return process.poll() is not None

    def _cleanup_after_abort(self, scope: str | None) -> tuple[Issue, ...]:
        issues: list[Issue] = []
        process = self._active_process
        subreaper = self._active_subreaper
        has_descendants = bool(subreaper.observe()) if subreaper is not None else False
        if process is None and has_descendants and subreaper is not None:
            if not subreaper.terminate_adopted():
                issues.append(self._tree_cleanup_issue())
        elif process is not None and (
            process.poll() is None
            or has_descendants
            or (subreaper is not None and subreaper.has_unbound_spawn)
        ):
            if not self._terminate(process):
                issues.append(self._tree_cleanup_issue())
        if scope is not None:
            issues.extend(self._cleanup_managed_containers(scope))
        return tuple(issues)

    def _cleanup_managed_containers(self, scope: str) -> tuple[Issue, ...]:
        cleaner = self._managed_container_cleaner
        if cleaner is None:
            return (
                Issue(
                    code="PROCESS_RUNNER_CONTAINER_CLEANER_UNAVAILABLE",
                    message="Managed container cleanup is unavailable.",
                    severity="error",
                    path="command_spec",
                    source="process_runner",
                ),
            )
        result = cleaner.cleanup(scope)
        return result.issues if result.is_failure else ()

    @staticmethod
    def _tree_cleanup_issue() -> Issue:
        return Issue(
            code="PROCESS_RUNNER_TREE_CLEANUP_FAILED",
            message="Workflow process cleanup could not be confirmed.",
            severity="error",
            path="command_spec",
            source="process_runner",
        )

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
    def _redaction_work_failure(
        cls,
        truncated: Mapping[str, bool],
    ) -> Result[ProcessResult]:
        return cls._failure_with_output_limit(
            Issue(
                code="PROCESS_RUNNER_REDACTION_WORK_EXCEEDED",
                message="Command output redaction exceeded its work bound.",
                severity="error",
                path="command_spec",
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
