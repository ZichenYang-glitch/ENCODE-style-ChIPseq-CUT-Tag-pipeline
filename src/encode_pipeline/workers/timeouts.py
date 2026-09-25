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


def _process_stat_fields(pid: int) -> tuple[str, int, int, int, int] | None:
    """Read state/ppid/pgrp/session/starttime; ``None`` only when truly absent."""
    try:
        value = Path(f"/proc/{pid}/stat").read_text()
    except FileNotFoundError:
        return None
    except (OSError, UnicodeError):
        raise ProcessRunnerCleanupError("Process state is unavailable.") from None
    closing = value.rfind(")")
    if closing < 0:
        raise ProcessRunnerCleanupError("Process state is invalid.")
    fields = value[closing + 2 :].split()
    if len(fields) < 20 or len(fields[0]) != 1:
        raise ProcessRunnerCleanupError("Process state is invalid.")
    try:
        return (
            fields[0],
            int(fields[1]),
            int(fields[2]),
            int(fields[3]),
            int(fields[19]),
        )
    except ValueError:
        raise ProcessRunnerCleanupError("Process state is invalid.") from None


def _member_visibility(pid: int, starttime: int) -> str:
    """Classify a *registered* member: ``frozen``, ``live``, ``gone``, ``changed``.

    The identity decision goes through the injectable owned-stat seam, so a
    denied or malformed identity read is never silently treated as an exit; only
    ``gone`` may retire a member from the freeze proof.
    """
    current = _owned_stat(pid)
    if current is None:
        return "gone"
    if current[1] != starttime:
        return "changed"
    fields = _process_stat_fields(pid)
    if fields is None:
        # The entry vanished between the two reads: an exit, not a change.
        return "gone"
    state, _ppid, _pgrp, _session, observed = fields
    if observed != starttime:
        return "changed"
    if state in {"Z", "X", "x"}:
        return "gone"
    if state in {"T", "t"}:
        return "frozen"
    return "live"


def _live_process_census() -> tuple[tuple[int, int, str, int, int, int], ...]:
    """Snapshot live processes as (pid, ppid, state, pgrp, session, starttime).

    An entry that cannot be read is skipped instead of failing this tree: F1
    established that a failure in an unrelated proc entry is not a failure of
    the owned tree. A live process in the worker's session was created by this
    worker's user and is readable, so an unreadable entry cannot be one of the
    horse's descendants -- the only shape that can share its process group.
    """
    try:
        with os.scandir("/proc") as entries:
            pids = sorted(
                int(entry.name) for entry in entries if entry.name.isdecimal()
            )
    except OSError:
        raise ProcessRunnerCleanupError("Process table is unavailable.") from None
    rows = []
    for pid in pids:
        try:
            fields = _process_stat_fields(pid)
        except ProcessRunnerCleanupError:
            continue
        if fields is None:
            continue
        state, ppid, pgrp, session, starttime = fields
        if state in {"Z", "X", "x"}:
            continue
        rows.append((pid, ppid, state, pgrp, session, starttime))
    return tuple(rows)


def _live_group_members(group_id: int) -> tuple[tuple[int, int], ...]:
    """Live processes whose kernel-maintained process group is ``group_id``.

    Group membership cannot be changed by a third party and survives
    reparenting, so an empty group is a positive completeness proof for every
    descendant that never started its own session -- including the orphan of a
    member that exited inside the freeze window.
    """
    return tuple(
        (pid, starttime)
        for pid, _ppid, _state, pgrp, _session, starttime in _live_process_census()
        if pgrp == group_id
    )


def _uncovered_processes(
    horse_pid: int,
    registered: tuple[tuple[int, int], ...],
) -> tuple[dict[str, object], ...]:
    """Report live processes the ownership proof cannot account for.

    Diagnostic evidence only, never a substitute for a proof. A member that
    exited before the subreaper scope adopted the tree leaves no kernel process
    table entry behind, so an independent session subtree it started cannot be
    enumerated here at all; this census is the best available signal and must
    not be read as proof of absence.
    """
    known = set(registered)
    worker_pid = os.getpid()
    try:
        worker_session: int | None = os.getsid(0)
    except OSError:  # pragma: no cover - diagnostic best effort
        worker_session = None
    seen = []
    for pid, ppid, state, pgrp, session, starttime in _live_process_census():
        if pid == horse_pid or (pid, starttime) in known:
            continue
        if pgrp == horse_pid:
            reason = "process-group member outside the registered tree"
        elif ppid == worker_pid:
            reason = "reparented live child outside the registered tree"
        else:
            continue
        seen.append(
            {
                "pid": pid,
                "starttime": starttime,
                "ppid": ppid,
                "state": state,
                "pgrp": pgrp,
                "session": session,
                "worker_session": worker_session,
                "reason": reason,
            }
        )
    return tuple(seen)


@dataclass(frozen=True)
class _OwnedCleanupReport:
    """Ownership evidence from one confirmed frozen-tree cleanup.

    ``registered`` is complete strictly before the first kill signal: every
    identity the cleanup relies on was observed alive with that exact
    starttime. ``exited`` lists members that terminated inside the freeze
    window -- their own identity is proven, but any descendant they reparented
    away before this scope adopted the tree is not covered by the proof.
    ``uncovered`` is diagnostic evidence, never part of the proof.

    ``confirmed`` carries the outcome the worker acts on. A cleanup whose proof
    did not complete records ``False`` with a fixed ``unconfirmed_reason``
    instead of raising, and ``horse_starttime`` stays ``None`` when the root
    identity was never established.
    """

    horse_pid: int
    horse_starttime: int | None
    registered: tuple[tuple[int, int], ...]
    exited: tuple[tuple[int, int], ...]
    uncovered: tuple[dict[str, object], ...]
    confirmed: bool = True
    unconfirmed_reason: str | None = None
    unconfirmed_detail: str | None = None


# Fixed code per proof-failure message this module can produce, so a worker log
# names a stable reason instead of relaying a free-form message.
_OWNED_CLEANUP_UNCONFIRMED_CODES = {
    "Owned process identity is unavailable.": "OWNED_PROCESS_IDENTITY_UNAVAILABLE",
    "Owned process identity is invalid.": "OWNED_PROCESS_IDENTITY_INVALID",
    "Owned process state is unavailable.": "OWNED_PROCESS_STATE_UNAVAILABLE",
    "Owned process state is invalid.": "OWNED_PROCESS_STATE_INVALID",
    "Process state is unavailable.": "PROCESS_STATE_UNAVAILABLE",
    "Process state is invalid.": "PROCESS_STATE_INVALID",
    "Process table is unavailable.": "PROCESS_TABLE_UNAVAILABLE",
    "Workflow tree was lost before freezing.": "OWNED_TREE_LOST_BEFORE_FREEZE",
    "Workflow parent is not frozen.": "OWNED_TREE_PARENT_NOT_FROZEN",
    "Owned process threads are unavailable.": "OWNED_TREE_THREADS_UNAVAILABLE",
    "Owned process children are unavailable.": "OWNED_TREE_CHILDREN_UNAVAILABLE",
    "Owned child changed during discovery.": "OWNED_TREE_CHILD_CHANGED",
    "Workflow parent resumed during discovery.": "OWNED_TREE_PARENT_RESUMED",
    "Owned child identity changed.": "OWNED_TREE_CHILD_IDENTITY_CHANGED",
    "Nested workflow root was lost.": "OWNED_TREE_ROOT_LOST",
    "Nested workflow tree did not stabilize.": "OWNED_TREE_DID_NOT_STABILIZE",
    "Nested workflow cleanup could not be confirmed.": "OWNED_TREE_CLEANUP_UNCONFIRMED",
    "Workflow root was lost before freezing.": "OWNED_TREE_ROOT_LOST_BEFORE_FREEZE",
}
_OWNED_CLEANUP_UNCONFIRMED_FALLBACK = "OWNED_TREE_CLEANUP_UNCONFIRMED"


def _unconfirmed_reason_code(error: BaseException) -> str:
    """Map a cleanup failure to a fixed code; never relay a raw message."""
    message = str(error)
    if message in _OWNED_CLEANUP_UNCONFIRMED_CODES:
        return _OWNED_CLEANUP_UNCONFIRMED_CODES[message]
    if getattr(error, "errno", None) == errno.EPERM:
        return "OWNED_TREE_SIGNAL_DENIED"
    return _OWNED_CLEANUP_UNCONFIRMED_FALLBACK


def _kill_owned_horse_tree(
    horse_pid: int,
    root_starttime: int,
    kill_horse: Callable[[], None],
) -> _OwnedCleanupReport:
    """Freeze an owned tree before killing it, including nested sessions.

    The RQ parent retains its direct horse wait status. Only captured descendant
    identities are adopted/reaped; unrelated children are never collected via
    the subreaper's generic newly-adopted-child discovery.

    A registered descendant that exits *inside* the freeze window retires from
    the freeze proof instead of failing it. Its parent is stopped and cannot
    ``wait`` for it, so a Nextflow submission burst makes such exits routine;
    treating them as a lost tree strands the run in ``running`` forever. The
    retirement is bounded:

    * the root must stay live and frozen -- a lost root remains fatal;
    * a member whose identity changed remains fatal;
    * every identity the cleanup relies on is registered before any kill;
    * after the group kill no live process may remain in the horse's
      kernel-maintained process group.

    Residual, deliberately outside the proof: a member that exited *before* this
    scope adopted the tree, having already started its own session. Its
    independent session subtree has no kernel entry left to enumerate, so it
    cannot be attributed or killed. ``_uncovered_processes`` reports whatever is
    still observable; it is diagnostic evidence and never a substitute for the
    proof.
    """
    tracked: dict[int, int] = {}
    exited: dict[int, int] = {}
    cleaned = False
    uncovered: tuple[dict[str, object], ...] = ()
    registered: tuple[tuple[int, int], ...] = ((horse_pid, root_starttime),)
    try:
        with _LinuxSubreaper(strict=True):
            root = _owned_stat(horse_pid)
            if root is None or root[1] != root_starttime:
                raise ProcessRunnerCleanupError("Nested workflow root was lost.")
            try:
                _signal_owned(horse_pid, root_starttime, signal.SIGSTOP)
                deadline = time.monotonic() + 2.0
                previous = None

                def retire_exited_members() -> tuple[bool, bool]:
                    """Retire members that exited; report (retired, all frozen).

                    Every registered non-root identity is read through the
                    injectable owned-stat seam, so a denied or malformed identity
                    read stays fatal instead of passing as an exit. A member that
                    exited is retired with its identity recorded; only the root
                    keeps the strict live-and-frozen requirement.
                    """
                    retired_any = False
                    all_frozen = True
                    for pid, starttime in list(tracked.items()):
                        visibility = _member_visibility(pid, starttime)
                        if visibility == "changed":
                            raise ProcessRunnerCleanupError(
                                "Owned child identity changed."
                            )
                        if visibility == "gone":
                            exited[pid] = starttime
                            tracked.pop(pid, None)
                            retired_any = True
                            continue
                        if visibility != "frozen":
                            all_frozen = False
                    return retired_any, all_frozen

                while time.monotonic() < deadline:
                    # Parent-before-child stopping closes the fork window. A
                    # registered descendant that exited inside it is retired here
                    # with its identity recorded; only the root, and any identity
                    # that changed underneath us, stays fatal.
                    _require_frozen_member(horse_pid, root_starttime)
                    if retire_exited_members()[0]:
                        previous = None
                        time.sleep(0.02)
                        continue
                    for pid, starttime in tracked.items():
                        _signal_owned(pid, starttime, signal.SIGSTOP)
                    retired, all_frozen = retire_exited_members()
                    if retired or not all_frozen:
                        time.sleep(0.02)
                        continue
                    for pid in (horse_pid, *tracked):
                        for child, starttime in _linux_descendants(pid):
                            known = tracked.get(child, exited.get(child))
                            if known is None:
                                tracked[child] = starttime
                            elif known != starttime:
                                raise ProcessRunnerCleanupError(
                                    "Owned child identity changed."
                                )
                    # A member discovered in this pass may have exited already --
                    # its stopped parent cannot reap it -- so retire it before
                    # the next pass turns it into a fatal zombie requirement.
                    if retire_exited_members()[0]:
                        previous = None
                        time.sleep(0.02)
                        continue
                    identities = tuple(sorted(tracked.items()))
                    members = ((horse_pid, root_starttime), *identities)
                    # Inspect every member even if an earlier one is not yet
                    # stopped: a lost later member cannot be hidden by all().
                    stopped = [
                        _stopped_or_gone(pid, starttime) for pid, starttime in members
                    ]
                    retired, all_frozen = retire_exited_members()
                    if retired:
                        previous = None
                        time.sleep(0.02)
                        continue
                    if all_frozen and all(stopped) and identities == previous:
                        # Reconfirm after the last reads; a zombie root is never
                        # a substitute for this proof.
                        retired, all_frozen = retire_exited_members()
                        if (
                            not retired
                            and all_frozen
                            and _require_frozen_member(horse_pid, root_starttime)
                        ):
                            break
                    previous = identities
                    time.sleep(0.02)
                else:
                    raise ProcessRunnerCleanupError(
                        "Nested workflow tree did not stabilize."
                    )
                # Registration is complete here, strictly before the first kill
                # below: every identity this cleanup signals or vouches for was
                # observed alive with this exact starttime.
                registered = (
                    (horse_pid, root_starttime),
                    *tuple(sorted(tracked.items())),
                    *tuple(sorted(exited.items())),
                )
            finally:
                # Best-effort termination is limited to still-confirmed owned
                # identities. Unknown identities may remain alive or stopped;
                # that failure must never become a cleanup acknowledgement.
                failed = False
                for pid, starttime in reversed((*tracked.items(), *exited.items())):
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
                    for pid in exited:
                        # Reaping an already-exited member is hygiene only: a
                        # non-child wait is ignored and it can never acknowledge
                        # the cleanup.
                        try:
                            os.waitpid(pid, os.WNOHANG)
                        except (ChildProcessError, OSError):
                            pass
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
                    # Group completeness: membership survives reparenting, so an
                    # empty group covers every same-group descendant, including
                    # the orphan of a member that exited in the freeze window.
                    try:
                        group_live = _live_group_members(horse_pid)
                    except ProcessRunnerCleanupError:
                        failed = True
                        group_live = ()
                    if failed:
                        break
                    if not live and not group_live:
                        cleaned = True
                        break
                    time.sleep(0.02)
                try:
                    uncovered = _uncovered_processes(horse_pid, registered)
                except ProcessRunnerCleanupError:  # pragma: no cover - diagnostics
                    uncovered = (
                        {"reason": "uncovered-process census was unavailable"},
                    )
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
    return _OwnedCleanupReport(
        horse_pid=horse_pid,
        horse_starttime=root_starttime,
        registered=registered,
        exited=tuple(sorted(exited.items())),
        uncovered=uncovered,
    )


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
            if stopped is None:
                return None
            _, job_id, horse_pid, starttime = context
            # Fail closed *without* raising. An unconfirmed cleanup must not
            # acknowledge the stop, and it must not end the worker either: RQ
            # reads this twice per job -- once to choose the stopped branch and
            # once inside ``handle_job_failure`` -- and a raise from either
            # unwinds ``work()``'s bare ``except``, which ends the worker while
            # the run is still ``running``. Returning ``None`` instead sends RQ
            # down its unexpected-termination branch, which reports the job as
            # failed. Taking the lock keeps the refusal behind an in-progress
            # ``kill_horse``, so no result is decided before the proof returns.
            if (
                stopped != job_id
                or starttime is None
                or getattr(self, "_nested_cleanup_completed", None)
                != (job_id, horse_pid, starttime)
                or getattr(self, "_nested_cleanup_failed_pid", None) == horse_pid
            ):
                return None
            return stopped

    @_stopped_job_id.setter
    def _stopped_job_id(self, value):
        with self._nested_cleanup_lock():
            self.__dict__["_nested_stopped_job_id"] = value

    def handle_payload(self, message):
        """Run an external command without letting it end the pubsub thread.

        RQ's ``handle_payload`` calls ``handle_command``, which has no exception
        handling, and the pubsub thread's own handler re-raises everything that
        is not a Redis connection error. One escaping cleanup failure therefore
        ends the thread and, because ``workers/cli.py`` discards ``work()``'s
        return value, still exits the worker with status ``0``. Redis transport
        failures arrive as ``redis.exceptions.RedisError`` rather than
        ``OSError``, so RQ's connection-retry contract is left untouched.
        """
        try:
            return super().handle_payload(message)
        except (ProcessRunnerCleanupError, OSError) as error:
            self._record_unconfirmed_cleanup(self.horse_pid, None, error)
            return None

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
            self._release_consumed_stop(job.id)

    def kill_horse(self, sig: signal.Signals = signal.SIGKILL):
        """Kill only the horse group, never its pre-``setpgrp`` parent group.

        The cleanup outcome is data, never an escaping exception. RQ calls this
        from its pubsub thread -- where an exception ends the thread and the
        worker -- and from its own monitor deadline path, where it unwinds
        ``work()``. Either way the worker would die *without* acknowledging the
        stop: the run stays ``running`` and, under a fail-fast supervisor, the
        worker's exit tears down every other platform service.
        """
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
            try:
                self._kill_owned_horse(horse_pid, sig)
            except (ProcessRunnerCleanupError, OSError) as error:
                self._record_unconfirmed_cleanup(horse_pid, sig, error)
            return None

    def _kill_owned_horse(self, horse_pid: int, sig: signal.Signals) -> None:
        process_group = None
        for attempt in range(5):
            if self.horse_pid != horse_pid:
                self.log.debug("RQ work horse changed before process-group kill")
                return
            try:
                process_group = os.getpgid(horse_pid)
            except OSError as exc:
                if exc.errno == errno.ESRCH:
                    self.log.debug("RQ work horse is already gone")
                    return
                # The group is unreadable, so no group signal may be sent.
                # Fall back to the identity-guarded single-process signal
                # below: leaving the horse alive would block RQ's wait4 and
                # strand the run instead of reporting it failed.
                self.log.debug("RQ work horse process group is unavailable: %s", exc)
                process_group = None
                break
            if process_group == horse_pid:
                break
            if attempt < 4:
                time.sleep(0.01)

        if self.horse_pid != horse_pid:
            self.log.debug("RQ work horse changed before termination")
            return

        def kill_owned_horse() -> None:
            try:
                if process_group == horse_pid:
                    os.killpg(horse_pid, sig)
                    self.log.info("Killed RQ work horse process group %s", horse_pid)
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
            report = _kill_owned_horse_tree(horse_pid, owned_root[1], kill_owned_horse)
            # Evidence before acknowledgement: the completion marker below
            # is only written once the relaxation and every process the
            # proof could not account for have been recorded.
            self._record_owned_cleanup(report)
            self._nested_cleanup_completed = (
                self.__dict__.get("_nested_stopped_job_id"),
                horse_pid,
                owned_root[1],
            )
            self._nested_cleanup_failed_pid = None
        else:
            kill_owned_horse()

    def _record_unconfirmed_cleanup(
        self,
        horse_pid: int,
        sig: signal.Signals | None,
        error: BaseException,
    ) -> None:
        """Record a cleanup whose ownership proof did not complete.

        The refusal is data: the stop stays unacknowledged, so RQ reports the
        job as failed, and the worker keeps running. ``sig`` is ``None`` when
        the failure is reported from the command boundary rather than a kill.
        Repeats for the same stop and horse are collapsed into one record.
        """
        marker = (
            self.__dict__.get("_nested_stopped_job_id"),
            horse_pid,
            None if sig is None else int(sig),
        )
        if self.__dict__.get("_nested_cleanup_unconfirmed") == marker:
            return
        self.__dict__["_nested_cleanup_unconfirmed"] = marker
        self._record_owned_cleanup(
            _OwnedCleanupReport(
                horse_pid=horse_pid,
                horse_starttime=None,
                registered=(),
                exited=(),
                uncovered=(),
                confirmed=False,
                unconfirmed_reason=_unconfirmed_reason_code(error),
                unconfirmed_detail=f"{type(error).__name__}: {error}",
            )
        )

    def _release_consumed_stop(self, job_id: str) -> None:
        """Drop this job's stop and cleanup markers once its monitor returns.

        The failure marker is keyed by horse pid alone and a pid is reusable, so
        a marker left behind could refuse the acknowledgement of a later,
        unrelated job.
        """
        with self._nested_cleanup_lock():
            if self.__dict__.get("_nested_stopped_job_id") == job_id:
                self.__dict__["_nested_stopped_job_id"] = None
            for name in (
                "_nested_cleanup_failed_pid",
                "_nested_cleanup_completed",
                "_nested_cleanup_unconfirmed",
            ):
                self.__dict__.pop(name, None)

    def _record_owned_cleanup(self, report: _OwnedCleanupReport | None) -> None:
        """Surface a relaxed or incompletely attributed cleanup; never silently."""
        self.__dict__["_nested_cleanup_report"] = report
        if report is None:  # pragma: no cover - injected cleanup stubs
            return
        if not report.confirmed:
            self.log.warning(
                "Owned cleanup could not be confirmed (%s) for horse %s: %s. The "
                "stop stays unacknowledged and the run is reported as an "
                "unexpected execution failure; the worker keeps running.",
                report.unconfirmed_reason,
                report.horse_pid,
                report.unconfirmed_detail,
            )
            return
        self.log.info(
            "Owned cleanup confirmed: %d identity(ies) registered before the "
            "group signal (horse=%s starttime=%s)",
            len(report.registered),
            report.horse_pid,
            report.horse_starttime,
        )
        if report.exited:
            self.log.warning(
                "Owned cleanup retired %d member(s) that exited inside the "
                "freeze window: %s. A member that exited before this cleanup "
                "adopted the tree, having already started its own session, "
                "leaves an independent session subtree outside the ownership "
                "proof.",
                len(report.exited),
                report.exited,
            )
        if report.uncovered:
            self.log.warning(
                "Owned cleanup could not account for %d live process(es): %s",
                len(report.uncovered),
                report.uncovered,
            )

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
