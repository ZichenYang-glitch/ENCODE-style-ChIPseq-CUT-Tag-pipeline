"""Worker-only hard timeout control flow for RQ job execution."""

from __future__ import annotations

from collections.abc import Callable
from contextvars import ContextVar
from dataclasses import dataclass
import errno
import os
from pathlib import Path
import signal
import threading
import time

from rq import Worker
from rq.timeouts import (
    BaseTimeoutException,
    JobTimeoutException,
    UnixSignalDeathPenalty,
)


from encode_pipeline.services.process_runner import (
    ProcessRunnerCleanupError,
    _LinuxSubreaper,
    _ProcessStat,
)


def _owned_stat(pid: int) -> tuple[int, int] | None:
    """Only an absent proc entry proves absence; unreadable is not absent."""
    try:
        current = _ProcessStat.read(str(pid))
    except FileNotFoundError:
        return None
    except (OSError, ValueError):
        raise ProcessRunnerCleanupError(
            "Owned process identity is unavailable."
        ) from None
    if (
        not isinstance(current, tuple)
        or len(current) != 2
        or any(type(value) is not int or value < 0 for value in current)
    ):
        raise ProcessRunnerCleanupError("Owned process identity is invalid.")
    return current


def _linux_process_matches(pid: int, starttime: int) -> bool:
    current = _owned_stat(pid)
    return current is not None and current[1] == starttime


def _owned_state(pid: int, starttime: int) -> str | None:
    if not _linux_process_matches(pid, starttime):
        return None
    try:
        value = Path(f"/proc/{pid}/stat").read_text()
    except FileNotFoundError:
        return None
    except (OSError, UnicodeError):
        raise ProcessRunnerCleanupError("Owned process state is unavailable.") from None
    try:
        closing = value.rfind(")")
        fields = value[closing + 2 :].split()
        if closing < 0 or len(fields) < 20 or len(fields[0]) != 1:
            raise ValueError
        if int(fields[19]) != starttime:
            return None
        return fields[0]
    except (ValueError, IndexError):
        raise ProcessRunnerCleanupError("Owned process state is invalid.") from None


def _require_frozen_member(pid: int, starttime: int) -> bool:
    state = _owned_state(pid, starttime)
    if state is None or state in {"Z", "X", "x"}:
        # A dead parent can already have reparented unseen children. Neither a
        # matching starttime on its zombie nor an empty tree proves ownership.
        raise ProcessRunnerCleanupError("Workflow tree was lost before freezing.")
    return state in {"T", "t"}


def _linux_descendants(root_pid: int) -> tuple[tuple[int, int], ...]:
    """Return one generation of descendants of a confirmed frozen parent.

    The outer loop registers and freezes this generation before exploring its
    children. Reading only owned thread children lists avoids treating failures
    in unrelated proc entries as a failure of this tree (or as an empty tree).
    """
    root = _owned_stat(root_pid)
    if root is None or not _require_frozen_member(root_pid, root[1]):
        raise ProcessRunnerCleanupError("Workflow parent is not frozen.")
    try:
        with os.scandir(f"/proc/{root_pid}/task") as entries:
            tids = sorted(int(entry.name) for entry in entries if entry.name.isdigit())
        if not tids:
            raise ProcessRunnerCleanupError("Owned process threads are unavailable.")
        children = set()
        for tid in tids:
            words = Path(f"/proc/{root_pid}/task/{tid}/children").read_text().split()
            for word in words:
                if not word.isdecimal() or int(word) <= 0:
                    raise ValueError
                children.add(int(word))
    except (OSError, ValueError):
        raise ProcessRunnerCleanupError(
            "Owned process children are unavailable."
        ) from None
    found = []
    for child in sorted(children):
        current = _owned_stat(child)
        if current is None or current[0] != root_pid:
            raise ProcessRunnerCleanupError("Owned child changed during discovery.")
        found.append((child, current[1]))
    if not _require_frozen_member(root_pid, root[1]):
        raise ProcessRunnerCleanupError("Workflow parent resumed during discovery.")
    return tuple(found)


def _signal_owned(pid: int, starttime: int, sig: signal.Signals) -> None:
    """Never signal a reused PID, including during the final cleanup pass."""
    if not _linux_process_matches(pid, starttime):
        return
    try:
        os.kill(pid, sig)
    except ProcessLookupError:
        return


def _stopped_or_gone(pid: int, starttime: int) -> bool:
    state = _owned_state(pid, starttime)
    return state is None or state in {"T", "t", "Z", "X", "x"}


def _kill_owned_horse_tree(
    horse_pid: int,
    root_starttime: int,
    kill_horse: Callable[[], None],
) -> None:
    """Freeze an owned tree before killing it, including nested sessions.

    The RQ parent retains its direct horse wait status. Only captured descendant
    identities are adopted/reaped; unrelated children are never collected via
    the subreaper's generic newly-adopted-child discovery.
    """
    tracked: dict[int, int] = {}
    cleaned = False
    try:
        with _LinuxSubreaper(strict=True):
            root = _owned_stat(horse_pid)
            if root is None or root[1] != root_starttime:
                raise ProcessRunnerCleanupError("Nested workflow root was lost.")
            try:
                _signal_owned(horse_pid, root_starttime, signal.SIGSTOP)
                deadline = time.monotonic() + 2.0
                previous = None
                while time.monotonic() < deadline:
                    # Parent-before-child stopping closes the fork window. A
                    # second full scan must see the same stopped identities.
                    members = ((horse_pid, root_starttime), *tracked.items())
                    frozen = [
                        _require_frozen_member(pid, starttime)
                        for pid, starttime in members
                    ]
                    if not all(frozen):
                        time.sleep(0.02)
                        continue
                    for pid, _starttime in members:
                        for child, starttime in _linux_descendants(pid):
                            if child in tracked and tracked[child] != starttime:
                                raise ProcessRunnerCleanupError(
                                    "Owned child identity changed."
                                )
                            tracked[child] = starttime
                    for pid, starttime in tracked.items():
                        _signal_owned(pid, starttime, signal.SIGSTOP)
                    identities = tuple(sorted(tracked.items()))
                    members = ((horse_pid, root_starttime), *identities)
                    # Inspect every member even if an earlier one is not yet
                    # stopped: a lost later member cannot be hidden by all().
                    frozen = [
                        _require_frozen_member(pid, starttime)
                        for pid, starttime in members
                    ]
                    stopped = all(frozen) and all(
                        _stopped_or_gone(pid, starttime) for pid, starttime in members
                    )
                    if stopped and identities == previous:
                        # Reconfirm live frozen identities after the last reads;
                        # a zombie root is never a substitute for this proof.
                        confirmed = [
                            _require_frozen_member(pid, starttime)
                            for pid, starttime in members
                        ]
                        if all(confirmed):
                            break
                    previous = identities
                    time.sleep(0.02)
                else:
                    raise ProcessRunnerCleanupError(
                        "Nested workflow tree did not stabilize."
                    )
            finally:
                # Best-effort termination is limited to still-confirmed owned
                # identities. Unknown identities may remain alive or stopped;
                # that failure must never become a cleanup acknowledgement.
                failed = False
                for pid, starttime in reversed(tuple(tracked.items())):
                    try:
                        _signal_owned(pid, starttime, signal.SIGKILL)
                    except (OSError, ProcessRunnerCleanupError):
                        failed = True
                try:
                    if _linux_process_matches(horse_pid, root_starttime):
                        kill_horse()
                except (OSError, ProcessRunnerCleanupError):
                    failed = True
                deadline = time.monotonic() + 2.0
                while time.monotonic() < deadline:
                    live = []
                    for pid, starttime in tracked.items():
                        try:
                            if not _linux_process_matches(pid, starttime):
                                continue
                            try:
                                os.waitpid(pid, os.WNOHANG)
                            except ChildProcessError:
                                pass
                            if _linux_process_matches(pid, starttime):
                                live.append(pid)
                        except (OSError, ProcessRunnerCleanupError):
                            failed = True
                            live.append(pid)
                    if failed:
                        break
                    if not live:
                        cleaned = not failed
                        break
                    time.sleep(0.02)
                if not cleaned:
                    raise ProcessRunnerCleanupError(
                        "Nested workflow cleanup could not be confirmed."
                    )
    except ProcessRunnerCleanupError:
        raise
    except OSError:
        raise ProcessRunnerCleanupError(
            "Nested workflow cleanup could not be confirmed."
        ) from None


@dataclass
class _ExecutionPhaseState:
    phase: str


_EXECUTION_PHASE: ContextVar[_ExecutionPhaseState | None] = ContextVar(
    "encode_pipeline_worker_execution_phase",
    default=None,
)
_MISSING = object()
_STOP_MONITOR: ContextVar[tuple[object, str, int, int | None] | None] = ContextVar(
    "encode_pipeline_stop_monitor",
    default=None,
)


class WorkerHardTimeout(BaseException):
    """An RQ job deadline that must bypass application ``Exception`` handlers."""


class WorkerUnixSignalDeathPenalty(UnixSignalDeathPenalty):
    """Use hard control flow only for RQ's main job timeout exception."""

    def __init__(
        self,
        timeout,
        exception=BaseTimeoutException,
        **kwargs,
    ) -> None:
        super().__init__(timeout, exception, **kwargs)

    def handle_death_penalty(self, signum, frame):
        state = _EXECUTION_PHASE.get()
        if (
            self._exception is JobTimeoutException
            and state is not None
            and state.phase == "main_job"
        ):
            raise WorkerHardTimeout(
                f"Task exceeded maximum timeout value ({self._timeout} seconds)"
            )
        return super().handle_death_penalty(signum, frame)


class DurableWorker(Worker):
    """RQ worker with hard main-job deadlines and native callback deadlines."""

    death_penalty_class = WorkerUnixSignalDeathPenalty

    def _nested_cleanup_lock(self):
        return self.__dict__.setdefault("_nested_horse_cleanup_lock", threading.RLock())

    @property
    def _stopped_job_id(self):
        context = _STOP_MONITOR.get()
        if context is None or context[0] is not self:
            # RQ's command thread and diagnostic readers do not acknowledge a
            # stop. In particular, they must not block behind a paused kill.
            return self.__dict__.get("_nested_stopped_job_id")
        with self._nested_cleanup_lock():
            stopped = self.__dict__.get("_nested_stopped_job_id")
            if stopped is not None:
                _, job_id, horse_pid, starttime = context
                if (
                    stopped != job_id
                    or starttime is None
                    or getattr(self, "_nested_cleanup_completed", None)
                    != (job_id, horse_pid, starttime)
                    or getattr(self, "_nested_cleanup_failed_pid", None) == horse_pid
                ):
                    raise ProcessRunnerCleanupError(
                        "Nested workflow cleanup could not be confirmed."
                    )
            return stopped

    @_stopped_job_id.setter
    def _stopped_job_id(self, value):
        with self._nested_cleanup_lock():
            self.__dict__["_nested_stopped_job_id"] = value

    def monitor_work_horse(self, job, queue):
        horse_pid = self.horse_pid
        try:
            root = _owned_stat(horse_pid) if horse_pid > 0 else None
        except ProcessRunnerCleanupError:
            root = None
        token = _STOP_MONITOR.set(
            (self, job.id, horse_pid, None if root is None else root[1])
        )
        try:
            # Keep the original RQ wait4, branch, callbacks and failure handler.
            # Every stop-marker consumption through the end of this monitor
            # requires the same proof, including a stop arriving after wait4.
            return super().monitor_work_horse(job, queue)
        finally:
            _STOP_MONITOR.reset(token)

    def wait_for_horse(self):
        horse_pid = self.horse_pid
        try:
            root = _owned_stat(horse_pid) if horse_pid > 0 else None
        except ProcessRunnerCleanupError:
            root = None
        # Never hold the stop lock across wait4: the stop thread must be able
        # to terminate the horse that RQ is waiting for.
        result = super().wait_for_horse()
        with self._nested_cleanup_lock():
            failed_pid = getattr(self, "_nested_cleanup_failed_pid", None)
            stopped_job = self.__dict__.get("_nested_stopped_job_id")
            completed = getattr(self, "_nested_cleanup_completed", None)
            # RQ records the stop request before invoking kill_horse. A wait
            # can reach this point even before that method obtains the lock.
            # Only a matching completed cleanup may acknowledge that request.
            unconfirmed_stop = stopped_job is not None and (
                root is None
                or result[0] != horse_pid
                or completed != (stopped_job, horse_pid, root[1])
            )
            if unconfirmed_stop or (failed_pid is not None and failed_pid == result[0]):
                raise ProcessRunnerCleanupError(
                    "Nested workflow cleanup could not be confirmed."
                )
        return result

    def kill_horse(self, sig: signal.Signals = signal.SIGKILL):
        """Kill only the horse group, never its pre-``setpgrp`` parent group."""
        with self._nested_cleanup_lock():
            horse_pid = self.horse_pid
            if horse_pid <= 0:
                self.log.debug("No live RQ work horse to kill")
                return None

            if sig == signal.SIGKILL:
                # Preparation is part of cleanup. Every early return/error
                # remains unconfirmed until the complete owned tree succeeds.
                self._nested_cleanup_failed_pid = horse_pid
                self._nested_cleanup_completed = None
            process_group = None
            for attempt in range(5):
                if self.horse_pid != horse_pid:
                    self.log.debug("RQ work horse changed before process-group kill")
                    return None
                try:
                    process_group = os.getpgid(horse_pid)
                except OSError as exc:
                    if exc.errno == errno.ESRCH:
                        self.log.debug("RQ work horse is already gone")
                        return None
                    raise
                if process_group == horse_pid:
                    break
                if attempt < 4:
                    time.sleep(0.01)

            if self.horse_pid != horse_pid:
                self.log.debug("RQ work horse changed before termination")
                return None

            def kill_owned_horse() -> None:
                try:
                    if process_group == horse_pid:
                        os.killpg(horse_pid, sig)
                        self.log.info(
                            "Killed RQ work horse process group %s", horse_pid
                        )
                    else:
                        # Never signal the worker's pre-setpgrp parent group.
                        os.kill(horse_pid, sig)
                        self.log.info("Killed pre-group RQ work horse %s", horse_pid)
                except OSError as exc:
                    if exc.errno == errno.ESRCH:
                        self.log.debug("RQ work horse is already gone")
                        return
                    raise

            if sig == signal.SIGKILL:
                owned_root = _owned_stat(horse_pid)
                if owned_root is None:
                    raise ProcessRunnerCleanupError(
                        "Workflow root was lost before freezing."
                    )
                _kill_owned_horse_tree(horse_pid, owned_root[1], kill_owned_horse)
                self._nested_cleanup_completed = (
                    self.__dict__.get("_nested_stopped_job_id"),
                    horse_pid,
                    owned_root[1],
                )
                self._nested_cleanup_failed_pid = None
            else:
                kill_owned_horse()
        return None

    def perform_job(self, job, queue) -> bool:
        job_attributes = getattr(job, "__dict__", None)
        original_override = (
            job_attributes.get("_execute", _MISSING)
            if isinstance(job_attributes, dict)
            else _MISSING
        )
        original_execute = job._execute

        def execute_with_hard_timeout():
            state = _ExecutionPhaseState(phase="main_job")
            token = _EXECUTION_PHASE.set(state)
            try:
                return original_execute()
            finally:
                _EXECUTION_PHASE.reset(token)

        job._execute = execute_with_hard_timeout
        try:
            return super().perform_job(job, queue)
        finally:
            if original_override is _MISSING:
                try:
                    del job._execute
                except AttributeError:  # pragma: no cover - slot-only RQ variant
                    job._execute = original_execute
            else:
                job._execute = original_override
