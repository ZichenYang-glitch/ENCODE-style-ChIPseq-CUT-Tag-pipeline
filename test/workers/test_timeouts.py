"""Tests for phase-aware worker-only RQ timeout control flow."""

from __future__ import annotations

import pytest
from rq import Worker
from rq.timeouts import (
    BaseTimeoutException,
    HorseMonitorTimeoutException,
    JobTimeoutException,
)

from encode_pipeline.workers.timeouts import (
    DurableWorker,
    WorkerHardTimeout,
    WorkerUnixSignalDeathPenalty,
    _EXECUTION_PHASE,
    _ExecutionPhaseState,
)


def _raised_timeout_type(penalty) -> type[BaseException]:
    try:
        penalty.handle_death_penalty(None, None)
    except BaseException as exc:
        return type(exc)
    raise AssertionError("death penalty did not raise")


def test_worker_hard_timeout_bypasses_application_exception_handlers():
    assert issubclass(WorkerHardTimeout, BaseException)
    assert not issubclass(WorkerHardTimeout, Exception)


def test_main_job_signal_raises_worker_hard_timeout_without_rewriting_penalty():
    penalty = WorkerUnixSignalDeathPenalty(17, JobTimeoutException)
    token = _EXECUTION_PHASE.set(_ExecutionPhaseState(phase="main_job"))
    try:
        assert penalty._exception is JobTimeoutException
        assert _raised_timeout_type(penalty) is WorkerHardTimeout
    finally:
        _EXECUTION_PHASE.reset(token)


@pytest.mark.parametrize("callback_phase", ("success", "failure", "stopped"))
def test_callback_job_timeout_remains_native(callback_phase):
    penalty = WorkerUnixSignalDeathPenalty(17, JobTimeoutException)

    assert callback_phase  # Names the three RQ callback paths explicitly.
    assert _EXECUTION_PHASE.get() is None
    assert _raised_timeout_type(penalty) is JobTimeoutException


def test_non_job_timeout_control_flow_is_preserved_inside_main_job():
    penalty = WorkerUnixSignalDeathPenalty(17, HorseMonitorTimeoutException)
    token = _EXECUTION_PHASE.set(_ExecutionPhaseState(phase="main_job"))
    try:
        assert _raised_timeout_type(penalty) is HorseMonitorTimeoutException
    finally:
        _EXECUTION_PHASE.reset(token)


def test_default_death_penalty_exception_is_preserved():
    penalty = WorkerUnixSignalDeathPenalty(17)

    assert penalty._exception is BaseTimeoutException


def test_durable_worker_limits_main_phase_to_user_execute(monkeypatch):
    phases: list[str | None] = []

    class FakeJob:
        def _execute(self):
            state = _EXECUTION_PHASE.get()
            phases.append(None if state is None else state.phase)
            return "result"

    job = FakeJob()
    original_execute = job._execute

    def fake_perform_job(_self, current_job, _queue):
        assert current_job._execute() == "result"
        phases.append(
            None if _EXECUTION_PHASE.get() is None else _EXECUTION_PHASE.get().phase
        )
        assert (
            _raised_timeout_type(WorkerUnixSignalDeathPenalty(17, JobTimeoutException))
            is JobTimeoutException
        )
        return True

    monkeypatch.setattr(Worker, "perform_job", fake_perform_job)

    assert DurableWorker.perform_job(object.__new__(DurableWorker), job, object())
    assert phases == ["main_job", None]
    assert job._execute == original_execute
    assert _EXECUTION_PHASE.get() is None


def test_durable_worker_resets_phase_after_user_exception(monkeypatch):
    class FakeJob:
        def _execute(self):
            assert _EXECUTION_PHASE.get().phase == "main_job"
            raise RuntimeError("user function failed")

    job = FakeJob()
    original_execute = job._execute

    def fake_perform_job(_self, current_job, _queue):
        try:
            current_job._execute()
        except RuntimeError:
            assert _EXECUTION_PHASE.get() is None
            assert (
                _raised_timeout_type(
                    WorkerUnixSignalDeathPenalty(17, JobTimeoutException)
                )
                is JobTimeoutException
            )
            return False
        raise AssertionError("user exception did not propagate")

    monkeypatch.setattr(Worker, "perform_job", fake_perform_job)

    assert not DurableWorker.perform_job(
        object.__new__(DurableWorker),
        job,
        object(),
    )
    assert job._execute == original_execute
    assert _EXECUTION_PHASE.get() is None


def test_prepare_failure_callback_timeout_remains_native(monkeypatch):
    class FakeJob:
        def _execute(self):  # pragma: no cover - prepare failure skips execution
            raise AssertionError("user function should not run")

    job = FakeJob()

    def fail_before_execute(_self, _job, _queue):
        assert _EXECUTION_PHASE.get() is None
        assert (
            _raised_timeout_type(WorkerUnixSignalDeathPenalty(17, JobTimeoutException))
            is JobTimeoutException
        )
        return False

    monkeypatch.setattr(Worker, "perform_job", fail_before_execute)

    assert not DurableWorker.perform_job(
        object.__new__(DurableWorker),
        job,
        object(),
    )
    assert _EXECUTION_PHASE.get() is None
