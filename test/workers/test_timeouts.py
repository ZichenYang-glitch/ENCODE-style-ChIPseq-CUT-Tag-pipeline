"""Tests for phase-aware worker-only RQ timeout control flow."""

from __future__ import annotations

import inspect
import os
import signal
from types import SimpleNamespace

import fakeredis
import pytest
import rq
from rq import Queue, Worker
from rq.command import send_stop_job_command
from rq.job import Job
from rq.serializers import JSONSerializer
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
from encode_pipeline.workers.jobs import handle_execution_stopped


def _raised_timeout_type(penalty) -> type[BaseException]:
    try:
        penalty.handle_death_penalty(None, None)
    except BaseException as exc:
        return type(exc)
    raise AssertionError("death penalty did not raise")


def _rq_execution_canary() -> str:
    return "original execution unexpectedly ran"


def test_supported_rq_minor_preserves_private_execute_contract():
    """Fail loudly before an RQ upgrade can bypass phase-aware timeouts."""
    major_minor = tuple(int(part) for part in rq.__version__.split(".")[:2])
    assert major_minor == (2, 10)
    assert tuple(inspect.signature(Job._execute).parameters) == ("self",)

    connection = fakeredis.FakeRedis()
    queue = Queue(
        "rq-compatibility-canary", connection=connection, serializer=JSONSerializer
    )
    job = queue.enqueue(_rq_execution_canary)
    calls: list[str] = []

    def replacement_execute():
        calls.append("patched")
        return "replacement execution ran"

    job._execute = replacement_execute

    assert job.perform() == "replacement execution ran"
    assert calls == ["patched"]


def test_rq_210_stop_callback_and_command_signatures_are_compatible():
    assert tuple(inspect.signature(send_stop_job_command).parameters) == (
        "connection",
        "job_id",
        "serializer",
    )
    assert tuple(inspect.signature(Job.execute_stopped_callback).parameters) == (
        "self",
        "death_penalty_class",
    )
    assert tuple(inspect.signature(handle_execution_stopped).parameters) == (
        "job",
        "_connection",
    )


def test_rq_worker_kill_horse_uses_process_group_sigkill(monkeypatch):
    worker = object.__new__(Worker)
    worker._horse_pid = 123
    worker.name = "canary-worker"
    worker.log = SimpleNamespace(info=lambda *_args: None, debug=lambda *_args: None)
    calls = []
    monkeypatch.setattr(os, "getpgid", lambda pid: 123 if pid == 123 else None)
    monkeypatch.setattr(
        os,
        "killpg",
        lambda process_group, signum: calls.append((process_group, signum)),
    )

    Worker.kill_horse(worker)

    assert calls == [(123, signal.SIGKILL)]


def test_rq_worker_horse_still_creates_its_own_process_group():
    assert "os.setpgrp()" in inspect.getsource(Worker.fork_work_horse)


def test_durable_worker_never_delegates_a_zero_horse_pid(monkeypatch):
    worker = object.__new__(DurableWorker)
    worker._horse_pid = 0
    worker.log = SimpleNamespace(debug=lambda *_args: None)
    delegated = []
    monkeypatch.setattr(
        Worker,
        "kill_horse",
        lambda *_args, **_kwargs: delegated.append(True),
    )

    DurableWorker.kill_horse(worker)

    assert delegated == []


def test_durable_worker_uses_positive_pid_snapshot_during_completion_race(
    monkeypatch,
):
    worker = object.__new__(DurableWorker)
    worker._horse_pid = 123
    worker.log = SimpleNamespace(
        info=lambda *_args: None,
        debug=lambda *_args: None,
    )
    calls = []

    def get_process_group(pid):
        assert pid == 123
        worker._horse_pid = 0
        return 123

    monkeypatch.setattr(os, "getpgid", get_process_group)
    monkeypatch.setattr(
        os,
        "killpg",
        lambda process_group, signum: calls.append((process_group, signum)),
    )

    DurableWorker.kill_horse(worker)

    assert calls == []


def test_durable_worker_never_kills_the_parent_process_group(monkeypatch):
    worker = object.__new__(DurableWorker)
    worker._horse_pid = 123
    worker.log = SimpleNamespace(
        info=lambda *_args: None,
        debug=lambda *_args: None,
    )
    group_calls = []
    process_calls = []
    monkeypatch.setattr(os, "getpgid", lambda pid: 456 if pid == 123 else None)
    monkeypatch.setattr(
        os,
        "killpg",
        lambda process_group, signum: group_calls.append((process_group, signum)),
    )
    monkeypatch.setattr(
        os,
        "kill",
        lambda pid, signum: process_calls.append((pid, signum)),
    )

    DurableWorker.kill_horse(worker)

    assert group_calls == []
    assert process_calls == [(123, signal.SIGKILL)]


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
