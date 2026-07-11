"""Worker-only hard timeout control flow for RQ job execution."""

from __future__ import annotations

from contextvars import ContextVar
from dataclasses import dataclass

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
