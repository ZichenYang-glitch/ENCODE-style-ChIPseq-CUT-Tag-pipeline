"""Worker-only hard timeout control flow for RQ job execution."""

from __future__ import annotations

from contextvars import ContextVar
from dataclasses import dataclass
import errno
import os
import signal
import time

from rq import Worker
from rq.timeouts import (
    BaseTimeoutException,
    JobTimeoutException,
    UnixSignalDeathPenalty,
)


@dataclass
class _ExecutionPhaseState:
    phase: str


_EXECUTION_PHASE: ContextVar[_ExecutionPhaseState | None] = ContextVar(
    "encode_pipeline_worker_execution_phase",
    default=None,
)
_MISSING = object()


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

    def kill_horse(self, sig: signal.Signals = signal.SIGKILL):
        """Kill only the horse group, never its pre-``setpgrp`` parent group."""
        horse_pid = self.horse_pid
        if horse_pid <= 0:
            self.log.debug("No live RQ work horse to kill")
            return None

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
        try:
            if process_group == horse_pid:
                os.killpg(horse_pid, sig)
                self.log.info("Killed RQ work horse process group %s", horse_pid)
            else:
                # The child calls setpgrp before entering the job. If that tiny
                # fork window persists, kill only the child PID; never signal
                # the worker's shared parent process group.
                os.kill(horse_pid, sig)
                self.log.info("Killed pre-group RQ work horse %s", horse_pid)
        except OSError as exc:
            if exc.errno == errno.ESRCH:
                self.log.debug("RQ work horse is already gone")
                return None
            raise
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
