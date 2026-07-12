"""Integration tests for the PR123 worker ownership handshake."""

from __future__ import annotations

from dataclasses import replace
from datetime import datetime, timezone
import os
from pathlib import Path
import shutil
from types import SimpleNamespace

import fakeredis
import pytest
from rq import SimpleWorker
from rq.serializers import JSONSerializer
from rq.timeouts import JobTimeoutException

import encode_pipeline.workers.jobs as worker_jobs
from encode_pipeline.persistence.runtime import open_run_persistence
from encode_pipeline.platform.execution import RunExecutionAssignment
from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.platform.results import Issue, Result
from encode_pipeline.platform.runs import RunStatus
from encode_pipeline.services.defaults import (
    create_default_run_service,
    create_default_workflow_registry,
)
from encode_pipeline.services.local_execution import LocalExecutionService
from encode_pipeline.services.process_runner import ProcessRunner
from encode_pipeline.services.workflow_builds import WorkflowBuildIdentityProvider
from encode_pipeline.workers.jobs import (
    handle_execution_stopped,
    handle_work_horse_killed,
)
from encode_pipeline.workers.rq_queue import RqRunQueue
from encode_pipeline.workers.settings import (
    QUEUE_NAME_ENV,
    REDIS_URL_ENV,
    WORKSPACE_ROOT_ENV,
)
from encode_pipeline.workers.timeouts import WorkerHardTimeout

from .conftest import create_planned_run, worker_settings


def _configure_worker_environment(monkeypatch, configured):
    monkeypatch.setenv("ENCODE_PIPELINE_DATABASE_URL", configured.database_url)
    monkeypatch.setenv(REDIS_URL_ENV, configured.redis_url)
    monkeypatch.setenv(QUEUE_NAME_ENV, configured.queue_name)
    monkeypatch.setenv(WORKSPACE_ROOT_ENV, str(configured.workspace_root))
    executable_dir = configured.workspace_root.parent / "test-bin"
    executable_dir.mkdir(parents=True, exist_ok=True)
    snakemake = executable_dir / "snakemake"
    snakemake.write_text(
        "#!/bin/sh\n"
        'while [ "$#" -gt 0 ]; do\n'
        '  if [ "$1" = "--directory" ]; then shift; cd "$1"; break; fi\n'
        "  shift\n"
        "done\n"
        "mkdir -p results/multiqc\n"
        "printf 'output_type\\tpath\\n' > results/multiqc/result_manifest.tsv\n"
        "printf 'worker stdout\\n'\n"
        "printf 'worker stderr\\n' >&2\n",
        encoding="utf-8",
    )
    snakemake.chmod(0o755)
    monkeypatch.setenv("PATH", f"{executable_dir}{os.pathsep}{os.environ['PATH']}")


def _run_burst(connection, run_queue):
    worker = SimpleWorker(
        [run_queue._queue],
        connection=connection,
        serializer=JSONSerializer,
    )
    return worker.work(burst=True, logging_level="WARNING")


def _read_run(configured, run_id):
    persistence = open_run_persistence(configured.database_url)
    service = create_default_run_service(
        registry=create_default_workflow_registry(),
        repository=persistence.repository,
    )
    try:
        return (
            service.get_run(run_id),
            service.list_events(run_id, limit=100),
            service.get_execution_assignment(run_id),
        )
    finally:
        persistence.close()


def _copy_controlled_project(destination: Path) -> None:
    source = Path(__file__).resolve().parents[2]
    destination.mkdir()
    shutil.copy2(source / "pyproject.toml", destination / "pyproject.toml")
    for relative in (
        Path("docs/architecture/artifact-inventory.yaml"),
        Path("src/encode_pipeline"),
        Path("workflow"),
        Path("profiles/default"),
        Path("scripts"),
    ):
        if (source / relative).is_dir():
            shutil.copytree(source / relative, destination / relative)
        else:
            (destination / relative).parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source / relative, destination / relative)


def _fail_if_process_starts(*_args, **_kwargs):
    raise AssertionError("ProcessRunner must not run before build verification")


def test_rq_worker_rebuilds_dependencies_and_persists_handshake_event(
    tmp_path, monkeypatch
):
    configured = worker_settings(tmp_path)
    assignment = create_planned_run(
        configured,
        "handshake-run",
        assign_queue=configured.queue_name,
    )
    _configure_worker_environment(monkeypatch, configured)
    connection = fakeredis.FakeRedis()
    run_queue = RqRunQueue(configured, connection=connection)
    assert assignment is not None
    run_queue.enqueue_execution(assignment)

    assert _run_burst(connection, run_queue) is True

    record, events, persisted_assignment = _read_run(configured, "handshake-run")
    handshake = [
        event for event in events if event.event_type == "worker_dependencies_rebuilt"
    ]
    assert record.status is RunStatus.SUCCEEDED
    assert persisted_assignment is not None
    assert persisted_assignment.dispatched_at is not None
    assert persisted_assignment.claimed_at is not None
    assert len(handshake) == 1
    assert handshake[0].context == {
        "backend": "rq",
        "job_id": assignment.job_id,
        "queue_name": configured.queue_name,
        "workflow_id": record.workflow_id,
    }
    persistence = open_run_persistence(configured.database_url)
    try:
        service = create_default_run_service(
            registry=create_default_workflow_registry(),
            repository=persistence.repository,
        )
        assert service.list_logs("handshake-run", "stdout")[0].lines == (
            "worker stdout",
        )
        assert service.list_logs("handshake-run", "stderr")[0].lines == (
            "worker stderr",
        )
        artifacts = service.list_artifacts("handshake-run")
        assert len(artifacts) == 1
        assert artifacts[0].metadata["output_type"] == "result_manifest"
        assert artifacts[0].metadata["relative_path"] == (
            "results/multiqc/result_manifest.tsv"
        )
    finally:
        persistence.close()

    job = run_queue._queue.fetch_job(assignment.job_id)
    assert job is not None
    job.delete()
    run_queue.enqueue_execution(assignment)
    assert _run_burst(connection, run_queue) is True

    _, events_after_retry, assignment_after_retry = _read_run(
        configured,
        "handshake-run",
    )
    assert assignment_after_retry == persisted_assignment
    assert [event.event_type for event in events_after_retry].count(
        "worker_dependencies_rebuilt"
    ) == 1


def test_rq_worker_rejects_stale_job_identity_without_writing_handshake(
    tmp_path, monkeypatch
):
    configured = worker_settings(tmp_path)
    assignment = create_planned_run(
        configured,
        "stale-run",
        assign_queue=configured.queue_name,
    )
    _configure_worker_environment(monkeypatch, configured)
    connection = fakeredis.FakeRedis()
    run_queue = RqRunQueue(configured, connection=connection)
    assert assignment is not None
    stale_assignment = replace(assignment, job_id="stale-job-id")
    run_queue.enqueue_execution(stale_assignment)

    assert _run_burst(connection, run_queue) is True

    job = run_queue._queue.fetch_job("stale-job-id")
    _, events, _ = _read_run(configured, "stale-run")
    assert job is not None
    assert job.is_failed
    assert all(event.event_type != "worker_dependencies_rebuilt" for event in events)


def test_rq_worker_rejects_missing_durable_assignment(tmp_path, monkeypatch):
    configured = worker_settings(tmp_path)
    create_planned_run(configured, "unassigned-run")
    _configure_worker_environment(monkeypatch, configured)
    connection = fakeredis.FakeRedis()
    run_queue = RqRunQueue(configured, connection=connection)
    orphan_assignment = RunExecutionAssignment(
        run_id="unassigned-run",
        job_id="orphan-job-id",
        backend="rq",
        queue_name=configured.queue_name,
        created_at=datetime.now(timezone.utc),
    )
    run_queue.enqueue_execution(orphan_assignment)

    assert _run_burst(connection, run_queue) is True

    job = run_queue._queue.fetch_job("orphan-job-id")
    _, events, _ = _read_run(configured, "unassigned-run")
    assert job is not None
    assert job.is_failed
    assert all(event.event_type != "worker_dependencies_rebuilt" for event in events)


def test_rq_worker_fails_legacy_planned_run_without_build_before_claim(
    tmp_path,
    monkeypatch,
):
    configured = worker_settings(tmp_path)
    assignment = create_planned_run(
        configured,
        "missing-build-run",
        assign_queue=configured.queue_name,
        bind_build_identity=False,
    )
    assert assignment is not None
    _configure_worker_environment(monkeypatch, configured)
    monkeypatch.setattr(ProcessRunner, "run", _fail_if_process_starts)
    connection = fakeredis.FakeRedis()
    run_queue = RqRunQueue(configured, connection=connection)
    run_queue.enqueue_execution(assignment)

    assert _run_burst(connection, run_queue) is True

    job = run_queue._queue.fetch_job(assignment.job_id)
    record, events, persisted_assignment = _read_run(
        configured,
        assignment.run_id,
    )
    assert job is not None
    assert job.is_failed
    assert record.status is RunStatus.FAILED
    assert record.error is not None
    assert record.error.code == "RUN_WORKFLOW_BUILD_IDENTITY_MISSING"
    assert events[-1].issue == record.error
    assert events[-1].context["reason_code"] == record.error.code
    assert persisted_assignment is not None
    assert persisted_assignment.dispatched_at is not None
    assert persisted_assignment.claimed_at is None
    assert all(event.event_type != "worker_dependencies_rebuilt" for event in events)


def test_rq_worker_rejects_project_a_build_on_project_b_before_process(
    tmp_path,
    monkeypatch,
):
    configured = worker_settings(tmp_path)
    project_a = tmp_path / "project-a"
    _copy_controlled_project(project_a)
    snakefile_a = project_a / "workflow" / "Snakefile"
    snakefile_a.write_text(
        snakefile_a.read_text(encoding="utf-8") + "\n# project A drift\n",
        encoding="utf-8",
    )
    registry = create_default_workflow_registry()
    identity_result = WorkflowBuildIdentityProvider(
        registry,
        project_root=project_a,
    ).capture("encode-style-chipseq-cuttag-atac-mnase")
    assert identity_result.is_success
    assignment = create_planned_run(
        configured,
        "mismatched-build-run",
        assign_queue=configured.queue_name,
        build_identity=identity_result.value,
    )
    assert assignment is not None
    _configure_worker_environment(monkeypatch, configured)
    monkeypatch.setattr(ProcessRunner, "run", _fail_if_process_starts)
    connection = fakeredis.FakeRedis()
    run_queue = RqRunQueue(configured, connection=connection)
    run_queue.enqueue_execution(assignment)

    assert _run_burst(connection, run_queue) is True

    job = run_queue._queue.fetch_job(assignment.job_id)
    record, events, persisted_assignment = _read_run(
        configured,
        assignment.run_id,
    )
    assert job is not None
    assert job.is_failed
    assert record.status is RunStatus.FAILED
    assert record.error is not None
    assert record.error.code == "RUN_WORKFLOW_BUILD_IDENTITY_MISMATCH"
    assert events[-1].issue == record.error
    assert events[-1].context["reason_code"] == record.error.code
    assert persisted_assignment is not None
    assert persisted_assignment.claimed_at is None
    assert all(event.event_type != "worker_dependencies_rebuilt" for event in events)


def test_rq_worker_fails_closed_when_local_build_cannot_be_fingerprinted(
    tmp_path,
    monkeypatch,
):
    configured = worker_settings(tmp_path)
    assignment = create_planned_run(
        configured,
        "unavailable-build-run",
        assign_queue=configured.queue_name,
    )
    assert assignment is not None
    _configure_worker_environment(monkeypatch, configured)

    def unavailable(_self, _workflow_id):
        return Result.failure(
            [
                Issue(
                    code="WORKFLOW_BUILD_SOURCE_UNAVAILABLE",
                    message="Controlled source is unavailable.",
                    severity="error",
                    source="test",
                )
            ]
        )

    monkeypatch.setattr(WorkflowBuildIdentityProvider, "capture", unavailable)
    monkeypatch.setattr(ProcessRunner, "run", _fail_if_process_starts)
    connection = fakeredis.FakeRedis()
    run_queue = RqRunQueue(configured, connection=connection)
    run_queue.enqueue_execution(assignment)

    assert _run_burst(connection, run_queue) is True

    job = run_queue._queue.fetch_job(assignment.job_id)
    record, events, persisted_assignment = _read_run(
        configured,
        assignment.run_id,
    )
    assert job is not None
    assert job.is_failed
    assert record.status is RunStatus.FAILED
    assert record.error is not None
    assert record.error.code == "RUN_WORKFLOW_BUILD_IDENTITY_UNAVAILABLE"
    assert events[-1].issue == record.error
    assert events[-1].context["reason_code"] == record.error.code
    assert persisted_assignment is not None
    assert persisted_assignment.claimed_at is None
    assert all(event.event_type != "worker_dependencies_rebuilt" for event in events)


def test_rq_worker_treats_queued_cancellation_as_clean_noop_without_process(
    tmp_path, monkeypatch
):
    configured = worker_settings(tmp_path)
    assignment = create_planned_run(
        configured,
        "cancelled-run",
        assign_queue=configured.queue_name,
    )
    assert assignment is not None
    _configure_worker_environment(monkeypatch, configured)
    connection = fakeredis.FakeRedis()
    run_queue = RqRunQueue(configured, connection=connection)
    run_queue.enqueue_execution(assignment)

    persistence = open_run_persistence(configured.database_url)
    try:
        run_service = create_default_run_service(
            registry=create_default_workflow_registry(),
            repository=persistence.repository,
        )
        run_service.mark_execution_dispatched(
            assignment.run_id,
            job_id=assignment.job_id,
        )
        queued = run_service.queue_dispatched_run(
            assignment.run_id,
            job_id=assignment.job_id,
            backend=assignment.backend,
            queue_name=assignment.queue_name,
        )
        assert queued.status is RunStatus.QUEUED
        run_service.cancel_run("cancelled-run", reason="Cancelled before worker claim.")
    finally:
        persistence.close()

    def fail_if_process_starts(*_args, **_kwargs):
        raise AssertionError("ProcessRunner must not run for a cancelled queued job")

    monkeypatch.setattr(ProcessRunner, "run", fail_if_process_starts)

    assert _run_burst(connection, run_queue) is True

    job = run_queue._queue.fetch_job(assignment.job_id)
    record, events, _ = _read_run(configured, "cancelled-run")
    assert job is not None
    assert job.is_finished
    assert record.status is RunStatus.CANCELLED
    assert all(event.event_type != "worker_dependencies_rebuilt" for event in events)
    persistence = open_run_persistence(configured.database_url)
    try:
        run_service = create_default_run_service(
            registry=create_default_workflow_registry(),
            repository=persistence.repository,
        )
        assert run_service.list_logs("cancelled-run", "stdout") == ()
        assert run_service.list_logs("cancelled-run", "stderr") == ()
    finally:
        persistence.close()


def test_rq_worker_persists_nonzero_execution_failure(tmp_path, monkeypatch):
    configured = worker_settings(tmp_path)
    assignment = create_planned_run(
        configured,
        "nonzero-run",
        assign_queue=configured.queue_name,
    )
    assert assignment is not None
    _configure_worker_environment(monkeypatch, configured)
    snakemake = configured.workspace_root.parent / "test-bin" / "snakemake"
    snakemake.write_text(
        "#!/bin/sh\nprintf 'before failure\\n'\nprintf 'safe error\\n' >&2\nexit 9\n",
        encoding="utf-8",
    )
    snakemake.chmod(0o755)
    connection = fakeredis.FakeRedis()
    run_queue = RqRunQueue(configured, connection=connection)
    run_queue.enqueue_execution(assignment)

    assert _run_burst(connection, run_queue) is True

    job = run_queue._queue.fetch_job(assignment.job_id)
    record, _events, persisted_assignment = _read_run(configured, "nonzero-run")
    assert job is not None
    assert job.is_failed
    assert record.status is RunStatus.FAILED
    assert record.error is not None
    assert record.error.code == "RUN_EXECUTION_FAILED"
    assert record.error.context == {"reason_code": "LOCAL_RUN_EXECUTION_FAILED"}
    assert persisted_assignment is not None
    assert persisted_assignment.claimed_at is not None


def test_rq_worker_sanitizes_unexpected_execution_exception(tmp_path, monkeypatch):
    configured = worker_settings(tmp_path)
    assignment = create_planned_run(
        configured,
        "unexpected-run",
        assign_queue=configured.queue_name,
    )
    assert assignment is not None
    _configure_worker_environment(monkeypatch, configured)

    def explode(_self, _run_id):
        raise RuntimeError("redis://private-password@internal:6379/0")

    monkeypatch.setattr(LocalExecutionService, "execute", explode)
    connection = fakeredis.FakeRedis()
    run_queue = RqRunQueue(configured, connection=connection)
    run_queue.enqueue_execution(assignment)

    assert _run_burst(connection, run_queue) is True

    job = run_queue._queue.fetch_job(assignment.job_id)
    record, _events, _assignment = _read_run(configured, "unexpected-run")
    assert job is not None
    assert job.is_failed
    failure = job.latest_result()
    assert failure is not None
    assert "private-password" not in (failure.exc_string or "")
    assert record.status is RunStatus.FAILED
    assert record.error is not None
    assert record.error.code == "RUN_WORKER_FAILED"
    assert "private-password" not in str(record.error.to_dict())


def test_rq_worker_persists_and_reraises_job_timeout(tmp_path, monkeypatch):
    configured = worker_settings(tmp_path)
    assignment = create_planned_run(
        configured,
        "timeout-run",
        assign_queue=configured.queue_name,
    )
    assert assignment is not None
    _configure_worker_environment(monkeypatch, configured)

    def timeout(_self, _run_id):
        raise WorkerHardTimeout("RQ deadline reached")

    monkeypatch.setattr(LocalExecutionService, "execute", timeout)
    connection = fakeredis.FakeRedis()
    run_queue = RqRunQueue(configured, connection=connection)
    run_queue.enqueue_execution(assignment)

    assert _run_burst(connection, run_queue) is True

    job = run_queue._queue.fetch_job(assignment.job_id)
    record, _events, _assignment = _read_run(configured, "timeout-run")
    assert job is not None
    assert job.is_failed
    failure = job.latest_result()
    assert failure is not None
    assert "WorkerHardTimeout" in (failure.exc_string or "")
    assert record.status is RunStatus.FAILED
    assert record.error is not None
    assert record.error.code == "RUN_WORKER_FAILED"
    assert record.error.context == {"reason_code": "WORKER_JOB_TIMEOUT"}


def test_rq_worker_generalizes_failure_when_durable_mapping_itself_fails(
    tmp_path,
    monkeypatch,
):
    configured = worker_settings(tmp_path)
    assignment = create_planned_run(
        configured,
        "mapping-failure-run",
        assign_queue=configured.queue_name,
    )
    assert assignment is not None
    _configure_worker_environment(monkeypatch, configured)

    def execution_failure(_self, _run_id):
        raise RuntimeError("redis://execution-secret@internal:6379/0")

    def mapping_failure(*_args, **_kwargs):
        raise RuntimeError("sqlite:///mapping-secret/private.db")

    monkeypatch.setattr(LocalExecutionService, "execute", execution_failure)
    monkeypatch.setattr(
        worker_jobs,
        "_record_unexpected_failure",
        mapping_failure,
    )
    connection = fakeredis.FakeRedis()
    run_queue = RqRunQueue(configured, connection=connection)
    run_queue.enqueue_execution(assignment)

    assert _run_burst(connection, run_queue) is True

    job = run_queue._queue.fetch_job(assignment.job_id)
    record, _events, _assignment = _read_run(configured, "mapping-failure-run")
    assert job is not None
    assert job.is_failed
    failure = job.latest_result()
    assert failure is not None
    assert "execution-secret" not in (failure.exc_string or "")
    assert "mapping-secret" not in (failure.exc_string or "")
    assert "WorkerExecutionError" in (failure.exc_string or "")
    assert record.status is RunStatus.QUEUED


def test_unexpected_failure_mapping_swallows_repository_errors():
    class BrokenRunService:
        def get_run(self, _run_id):
            raise RuntimeError("sqlite:///private-backend.db")

    worker_jobs._record_unexpected_failure(BrokenRunService(), "run-1")


def test_post_success_artifact_exception_does_not_fail_execution_job():
    recorded: list[str] = []

    class Extraction:
        def extract(self, _run_id):
            raise RuntimeError("/private/workspace")

        def record_unexpected_failure(self, run_id):
            recorded.append(run_id)

    runtime = SimpleNamespace(
        local_execution_service=SimpleNamespace(
            execute=lambda _run_id: Result.success(object())
        ),
        artifact_extraction_service=Extraction(),
    )

    worker_jobs._execute_claimed_run(runtime, "run-1")

    assert recorded == ["run-1"]


def test_failure_mapping_does_not_swallow_rq_timeout():
    class TimedOutRunService:
        def get_run(self, _run_id):
            raise WorkerHardTimeout("RQ deadline reached during failure mapping")

    with pytest.raises(WorkerHardTimeout, match="during failure mapping"):
        worker_jobs._record_unexpected_failure_safely(
            TimedOutRunService(),
            "run-1",
        )


def test_worker_composition_failure_is_public_safe_and_durably_failed(
    tmp_path,
    monkeypatch,
):
    configured = worker_settings(tmp_path)
    assignment = create_planned_run(
        configured,
        "composition-failure-run",
        assign_queue=configured.queue_name,
    )
    assert assignment is not None
    _configure_worker_environment(monkeypatch, configured)
    monkeypatch.setattr(
        worker_jobs,
        "get_current_job",
        lambda: SimpleNamespace(
            id=assignment.job_id,
            origin=configured.queue_name,
        ),
    )

    def fail_composition():
        raise RuntimeError("sqlite:///private/path/platform.db")

    monkeypatch.setattr(worker_jobs, "open_worker_runtime", fail_composition)

    with pytest.raises(worker_jobs.WorkerExecutionError) as raised:
        worker_jobs.run_execution_job("composition-failure-run")

    assert "private/path" not in str(raised.value)
    assert raised.value.__cause__ is None
    record, _events, persisted_assignment = _read_run(
        configured,
        "composition-failure-run",
    )
    assert record.status is RunStatus.FAILED
    assert record.error is not None
    assert record.error.context == {"reason_code": "WORKER_INITIALIZATION_FAILED"}
    assert persisted_assignment is not None
    assert persisted_assignment.dispatched_at is not None


@pytest.mark.parametrize("drift_field", ("job_id", "queue_name"))
def test_worker_composition_fallback_refuses_identity_drift_without_mutation(
    tmp_path,
    monkeypatch,
    drift_field,
):
    configured = worker_settings(tmp_path)
    assignment = create_planned_run(
        configured,
        f"identity-drift-{drift_field}",
        assign_queue=configured.queue_name,
    )
    assert assignment is not None
    _configure_worker_environment(monkeypatch, configured)
    current_job = SimpleNamespace(
        id="wrong-job" if drift_field == "job_id" else assignment.job_id,
        origin="wrong-queue" if drift_field == "queue_name" else configured.queue_name,
    )
    monkeypatch.setattr(worker_jobs, "get_current_job", lambda: current_job)
    monkeypatch.setattr(
        worker_jobs,
        "open_worker_runtime",
        lambda **_kwargs: (_ for _ in ()).throw(RuntimeError("composition failed")),
    )

    with pytest.raises(worker_jobs.WorkerExecutionError):
        worker_jobs.run_execution_job(assignment.run_id)

    record, events, persisted_assignment = _read_run(configured, assignment.run_id)
    assert record.status is RunStatus.PLANNED
    assert persisted_assignment == assignment
    assert all(event.status is not RunStatus.QUEUED for event in events)
    assert all(event.status is not RunStatus.FAILED for event in events)


def test_worker_hard_timeout_during_composition_uses_durable_fallback(
    tmp_path,
    monkeypatch,
):
    configured = worker_settings(tmp_path)
    assignment = create_planned_run(
        configured,
        "composition-timeout-run",
        assign_queue=configured.queue_name,
    )
    assert assignment is not None
    _configure_worker_environment(monkeypatch, configured)
    monkeypatch.setattr(
        worker_jobs,
        "get_current_job",
        lambda: SimpleNamespace(
            id=assignment.job_id,
            origin=configured.queue_name,
        ),
    )
    monkeypatch.setattr(
        worker_jobs,
        "open_worker_runtime",
        lambda **_kwargs: (_ for _ in ()).throw(
            WorkerHardTimeout("RQ deadline reached")
        ),
    )

    with pytest.raises(WorkerHardTimeout, match="RQ deadline reached"):
        worker_jobs.run_execution_job(assignment.run_id)

    record, _events, persisted_assignment = _read_run(configured, assignment.run_id)
    assert record.status is RunStatus.FAILED
    assert record.error is not None
    assert record.error.context == {"reason_code": "WORKER_JOB_TIMEOUT"}
    assert persisted_assignment is not None
    assert persisted_assignment.dispatched_at is not None


def test_missing_worker_adapter_is_initialization_failure_not_identity_error(
    tmp_path,
    monkeypatch,
):
    configured = worker_settings(tmp_path)
    assignment = create_planned_run(
        configured,
        "missing-adapter-run",
        assign_queue=configured.queue_name,
    )
    assert assignment is not None
    _configure_worker_environment(monkeypatch, configured)
    monkeypatch.setattr(
        worker_jobs,
        "get_current_job",
        lambda: SimpleNamespace(
            id=assignment.job_id,
            origin=configured.queue_name,
        ),
    )
    monkeypatch.setattr(
        "encode_pipeline.services.defaults.create_default_workflow_registry",
        lambda **_kwargs: WorkflowRegistry(),
    )

    with pytest.raises(worker_jobs.WorkerExecutionError):
        worker_jobs.run_execution_job(assignment.run_id)

    record, _events, _assignment = _read_run(configured, assignment.run_id)
    assert record.status is RunStatus.FAILED
    assert record.error is not None
    assert record.error.context == {"reason_code": "WORKER_INITIALIZATION_FAILED"}


def test_work_horse_death_is_mapped_to_durable_failure(tmp_path, monkeypatch):
    configured = worker_settings(tmp_path)
    assignment = create_planned_run(
        configured,
        "killed-run",
        assign_queue=configured.queue_name,
    )
    assert assignment is not None
    _configure_worker_environment(monkeypatch, configured)

    class KilledJob:
        id = assignment.job_id
        func_name = "encode_pipeline.workers.jobs.run_execution_job"
        args = (assignment.run_id,)
        kwargs = {}
        origin = configured.queue_name

    handle_work_horse_killed(KilledJob(), 123, 9, object())

    record, _events, persisted_assignment = _read_run(configured, "killed-run")
    assert record.status is RunStatus.FAILED
    assert record.error is not None
    assert record.error.code == "RUN_WORKER_FAILED"
    assert record.error.context == {"reason_code": "WORKER_PROCESS_TERMINATED"}
    assert persisted_assignment is not None
    assert persisted_assignment.dispatched_at is not None


def test_stopped_callback_acknowledges_user_intent_after_horse_exit(
    tmp_path,
    monkeypatch,
):
    configured = worker_settings(tmp_path)
    assignment = create_planned_run(
        configured,
        "cancelled-running-run",
        assign_queue=configured.queue_name,
    )
    assert assignment is not None
    _configure_worker_environment(monkeypatch, configured)
    claimed = _prepare_running_assignment(configured, assignment, request_cancel=True)
    job = _execution_job(claimed)

    handle_execution_stopped(job, fakeredis.FakeRedis())

    record, events, persisted = _read_run(configured, assignment.run_id)
    assert record.status is RunStatus.CANCELLED
    assert record.ended_at is not None
    assert record.cancellation_reason == "User requested cancellation."
    assert persisted is not None
    assert persisted.cancellation_acknowledged_at is not None
    assert [event.event_type for event in events].count(
        "cancellation_acknowledged"
    ) == 1

    handle_execution_stopped(job, fakeredis.FakeRedis())
    _record, repeated_events, repeated_assignment = _read_run(
        configured,
        assignment.run_id,
    )
    assert repeated_events == events
    assert repeated_assignment == persisted


def test_stopped_callback_without_user_intent_fails_truthfully(
    tmp_path,
    monkeypatch,
):
    configured = worker_settings(tmp_path)
    assignment = create_planned_run(
        configured,
        "unexpected-stopped-run",
        assign_queue=configured.queue_name,
    )
    assert assignment is not None
    _configure_worker_environment(monkeypatch, configured)
    claimed = _prepare_running_assignment(configured, assignment)

    handle_execution_stopped(_execution_job(claimed), fakeredis.FakeRedis())

    record, events, persisted = _read_run(configured, assignment.run_id)
    assert record.status is RunStatus.FAILED
    assert record.cancellation_reason is None
    assert record.error is not None
    assert record.error.code == "RUN_WORKER_STOPPED_UNEXPECTEDLY"
    assert record.error.context == {"reason_code": "WORKER_STOP_WITHOUT_CANCELLATION"}
    assert persisted is not None
    assert persisted.cancellation_acknowledged_at is None
    assert events[-1].event_type == "execution_stopped_unexpectedly"


@pytest.mark.parametrize("claimed", (False, True))
def test_stopped_callback_before_running_fails_truthfully(
    tmp_path,
    monkeypatch,
    claimed,
):
    configured = worker_settings(tmp_path)
    assignment = create_planned_run(
        configured,
        f"unexpected-queued-stop-{claimed}",
        assign_queue=configured.queue_name,
    )
    assert assignment is not None
    _configure_worker_environment(monkeypatch, configured)
    persistence = open_run_persistence(configured.database_url)
    try:
        service = create_default_run_service(
            create_default_workflow_registry(),
            repository=persistence.repository,
        )
        service.mark_execution_dispatched(
            assignment.run_id,
            job_id=assignment.job_id,
        )
        service.queue_dispatched_run(
            assignment.run_id,
            job_id=assignment.job_id,
            backend=assignment.backend,
            queue_name=assignment.queue_name,
        )
        if claimed:
            service.claim_execution_assignment(
                assignment.run_id,
                job_id=assignment.job_id,
                backend=assignment.backend,
                queue_name=assignment.queue_name,
            )
        queued_assignment = service.get_execution_assignment(assignment.run_id)
        assert queued_assignment is not None
    finally:
        persistence.close()

    handle_execution_stopped(
        _execution_job(queued_assignment),
        fakeredis.FakeRedis(),
    )

    record, events, persisted = _read_run(configured, assignment.run_id)
    assert record.status is RunStatus.FAILED
    assert record.error is not None
    assert record.error.code == "RUN_WORKER_STOPPED_UNEXPECTEDLY"
    assert persisted is not None
    assert persisted.cancellation_acknowledged_at is None
    assert events[-1].event_type == "execution_stopped_unexpectedly"
    assert events[-1].context["previous_status"] == "queued"


def test_stopped_callback_before_dispatch_fails_planned_truthfully(
    tmp_path,
    monkeypatch,
):
    configured = worker_settings(tmp_path)
    assignment = create_planned_run(
        configured,
        "unexpected-planned-stop",
        assign_queue=configured.queue_name,
    )
    assert assignment is not None
    _configure_worker_environment(monkeypatch, configured)

    handle_execution_stopped(
        _execution_job(assignment),
        fakeredis.FakeRedis(),
    )

    record, events, persisted = _read_run(configured, assignment.run_id)
    assert record.status is RunStatus.FAILED
    assert record.error is not None
    assert record.error.code == "RUN_WORKER_STOPPED_UNEXPECTEDLY"
    assert persisted is not None
    assert persisted.dispatched_at is None
    assert persisted.cancellation_acknowledged_at is None
    assert events[-1].event_type == "execution_stopped_unexpectedly"
    assert events[-1].context["previous_status"] == "planned"


def test_stopped_callback_preserves_natural_terminal_race(tmp_path, monkeypatch):
    configured = worker_settings(tmp_path)
    assignment = create_planned_run(
        configured,
        "naturally-finished-run",
        assign_queue=configured.queue_name,
    )
    assert assignment is not None
    _configure_worker_environment(monkeypatch, configured)
    claimed = _prepare_running_assignment(configured, assignment, request_cancel=True)
    persistence = open_run_persistence(configured.database_url)
    try:
        service = create_default_run_service(
            create_default_workflow_registry(),
            repository=persistence.repository,
        )
        succeeded = service.transition_run(assignment.run_id, RunStatus.SUCCEEDED)
        events_before = service.list_events(assignment.run_id, limit=100)
    finally:
        persistence.close()

    handle_execution_stopped(_execution_job(claimed), fakeredis.FakeRedis())

    record, events, persisted = _read_run(configured, assignment.run_id)
    assert record == succeeded
    assert events == events_before
    assert persisted is not None
    assert persisted.cancellation_acknowledged_at is None


@pytest.mark.parametrize("field", ["job_id", "queue_name", "run_id", "func_name"])
def test_stopped_callback_identity_drift_never_mutates_another_run(
    tmp_path,
    monkeypatch,
    field,
):
    configured = worker_settings(tmp_path)
    assignment = create_planned_run(
        configured,
        f"stopped-drift-{field}",
        assign_queue=configured.queue_name,
    )
    assert assignment is not None
    _configure_worker_environment(monkeypatch, configured)
    claimed = _prepare_running_assignment(configured, assignment, request_cancel=True)
    job = _execution_job(claimed)
    if field == "job_id":
        job.id = "wrong-job"
    elif field == "queue_name":
        job.origin = "wrong-queue"
    elif field == "run_id":
        job.args = ("another-run",)
    else:
        job.func_name = "other.worker"
    record_before, events_before, persisted_before = _read_run(
        configured,
        assignment.run_id,
    )

    handle_execution_stopped(job, fakeredis.FakeRedis())

    record, events, persisted = _read_run(configured, assignment.run_id)
    assert record == record_before
    assert events == events_before
    assert persisted == persisted_before


def test_stopped_callback_does_not_swallow_rq_callback_timeout(monkeypatch):
    monkeypatch.setattr(
        worker_jobs,
        "load_worker_settings",
        lambda: (_ for _ in ()).throw(JobTimeoutException("callback deadline")),
    )

    with pytest.raises(JobTimeoutException, match="callback deadline"):
        handle_execution_stopped(
            SimpleNamespace(
                id="job-1",
                origin="runs",
                func_name="encode_pipeline.workers.jobs.run_execution_job",
                args=("run-1",),
                kwargs={},
            ),
            fakeredis.FakeRedis(),
        )


def _prepare_running_assignment(configured, assignment, *, request_cancel=False):
    persistence = open_run_persistence(configured.database_url)
    try:
        service = create_default_run_service(
            create_default_workflow_registry(),
            repository=persistence.repository,
        )
        service.mark_execution_dispatched(
            assignment.run_id,
            job_id=assignment.job_id,
        )
        service.queue_dispatched_run(
            assignment.run_id,
            job_id=assignment.job_id,
            backend=assignment.backend,
            queue_name=assignment.queue_name,
        )
        service.claim_execution_assignment(
            assignment.run_id,
            job_id=assignment.job_id,
            backend=assignment.backend,
            queue_name=assignment.queue_name,
        )
        service.transition_run(assignment.run_id, RunStatus.RUNNING)
        claimed = service.get_execution_assignment(assignment.run_id)
        assert claimed is not None
        if request_cancel:
            requested = service.request_execution_cancellation(
                assignment.run_id,
                job_id=claimed.job_id,
                backend=claimed.backend,
                queue_name=claimed.queue_name,
                reason="User requested cancellation.",
            )
            claimed = requested.assignment
        return claimed
    finally:
        persistence.close()


def _execution_job(assignment):
    return SimpleNamespace(
        id=assignment.job_id,
        origin=assignment.queue_name,
        func_name="encode_pipeline.workers.jobs.run_execution_job",
        args=(assignment.run_id,),
        kwargs={},
    )
