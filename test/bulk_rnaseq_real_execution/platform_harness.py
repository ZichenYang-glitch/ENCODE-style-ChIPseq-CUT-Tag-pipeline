"""Real SQLite/Redis/RQ orchestration for full bulk RNA-seq acceptance."""

from __future__ import annotations

from dataclasses import dataclass
import os
from pathlib import Path
import signal
import subprocess
import sys
import time
from uuid import uuid4

from redis import Redis
from rq.job import JobStatus

from encode_pipeline.persistence import DATABASE_URL_ENV, open_run_persistence
from encode_pipeline.platform.managed_containers import managed_container_scope
from encode_pipeline.platform.runs import RunStatus
from encode_pipeline.services.managed_containers import ManagedContainerCleaner
from encode_pipeline.services.run_cancellation import (
    RunCancellationResult,
    RunCancellationService,
)
from encode_pipeline.services.run_submission import RunSubmissionService
from encode_pipeline.services.validated_inputs import (
    ValidatedInputService,
    ValidatedRunCreationService,
)
from encode_pipeline.services.validation import ValidationService
from encode_pipeline.workers.rq_queue import (
    RqRunQueue,
    create_api_redis_connection,
)
from encode_pipeline.workers.runtime import open_worker_runtime
from encode_pipeline.workers.settings import (
    JOB_TIMEOUT_SECONDS_ENV,
    MANAGED_DOCKER_EXECUTABLE_ENV,
    MANAGED_DOCKER_SOCKET_ENV,
    QUEUE_NAME_ENV,
    REDIS_URL_ENV,
    WORKSPACE_ROOT_ENV,
    load_worker_settings,
)
from encode_pipeline.workers.timeouts import WorkerHardTimeout

from workers.process_helpers import terminate_rq_worker

from .support import (
    FIXTURE_MANIFEST_ENV,
    MANAGED_DOCKER_EXECUTABLE_ENV as GATE_DOCKER_EXECUTABLE_ENV,
    MANAGED_DOCKER_SOCKET_ENV as GATE_DOCKER_SOCKET_ENV,
    REQUIRE_REAL_EXECUTION_ENV,
    RUNTIME_ROOT_ENV,
    TEST_REDIS_URL_ENV,
    AcceptanceEvidence,
    AcceptanceFixture,
    GateSettings,
    build_acceptance_process_runner,
    build_results_composition,
    collect_success_evidence,
    assert_no_managed_containers,
    load_acceptance_fixture,
    managed_container_ids,
)


_DEFAULT_ACCEPTANCE_TIMEOUT_SECONDS = 7_200
_MAX_ACCEPTANCE_TIMEOUT_SECONDS = 14_400
_ACTIVITY_POLL_SECONDS = 0.05


@dataclass(frozen=True)
class SubmittedAcceptanceRun:
    """Public identities retained across the independent worker boundary."""

    run_id: str
    job_id: str
    validated_snapshot_id: str
    fixture_acceptance_manifest_sha256: str


@dataclass(frozen=True)
class ExecutionActivityEvidence:
    """Observed real workflow activity before cancellation or timeout."""

    worker_session_id: int
    process_group_count: int
    nextflow_observed: bool
    managed_container_observed: bool


@dataclass(frozen=True)
class TerminalLifecycleEvidence:
    """Public-safe durable state after one non-successful real execution."""

    run_id: str
    job_id: str
    lifecycle_status: str
    lifecycle_history: tuple[str, ...]
    event_types: tuple[str, ...]
    assignment_dispatched: bool
    assignment_claimed: bool
    cancellation_requested: bool
    cancellation_acknowledged: bool
    cancellation_reason: str | None
    error_code: str | None
    error_reason_code: str | None
    artifact_revision: int
    artifact_attempt_id: str | None
    artifact_attempt_status: str | None
    qc_revision: int
    qc_attempt_id: str | None
    qc_attempt_status: str | None
    artifact_count: int
    qc_metric_count: int
    rq_status: str
    rq_failed: bool
    rq_stopped: bool
    rq_finished: bool
    cleanup_confirmed: bool


class PlatformAcceptanceHarness:
    """Own one isolated queue, SQLite database, workspace root, and worker set."""

    def __init__(
        self,
        *,
        gate_settings: GateSettings,
        repository_root: Path,
        temporary_root: Path,
        job_timeout_seconds: int | None = None,
    ) -> None:
        if not isinstance(gate_settings, GateSettings):
            raise ValueError("gate_settings must be GateSettings")
        for name, value in (
            ("repository_root", repository_root),
            ("temporary_root", temporary_root),
        ):
            if not isinstance(value, Path) or not value.is_absolute():
                raise ValueError(f"{name} must be an absolute Path")
        self.gate_settings = gate_settings
        self.repository_root = repository_root
        self.temporary_root = temporary_root
        self.workspace_root = temporary_root / "workspaces"
        self.database_url = f"sqlite:///{temporary_root / 'platform.db'}"
        self.queue_name = f"bulk-rnaseq-acceptance-{uuid4().hex}"
        if job_timeout_seconds is not None and (
            isinstance(job_timeout_seconds, bool)
            or not isinstance(job_timeout_seconds, int)
            or not 1 <= job_timeout_seconds <= _MAX_ACCEPTANCE_TIMEOUT_SECONDS
        ):
            raise ValueError(
                "job_timeout_seconds must be between 1 and the acceptance bound"
            )
        settings_environment = dict(os.environ)
        settings_environment.update(
            {
                DATABASE_URL_ENV: self.database_url,
                REDIS_URL_ENV: gate_settings.redis_url,
                QUEUE_NAME_ENV: self.queue_name,
                WORKSPACE_ROOT_ENV: str(self.workspace_root),
                MANAGED_DOCKER_EXECUTABLE_ENV: str(gate_settings.docker_executable),
                MANAGED_DOCKER_SOCKET_ENV: str(gate_settings.docker_socket),
                JOB_TIMEOUT_SECONDS_ENV: (
                    str(job_timeout_seconds)
                    if job_timeout_seconds is not None
                    else os.environ.get(
                        JOB_TIMEOUT_SECONDS_ENV,
                        str(_DEFAULT_ACCEPTANCE_TIMEOUT_SECONDS),
                    )
                ),
            }
        )
        self.worker_settings = load_worker_settings(settings_environment)
        if self.worker_settings.job_timeout_seconds > _MAX_ACCEPTANCE_TIMEOUT_SECONDS:
            raise ValueError("real acceptance job timeout exceeds the fixed bound")
        self.composition = build_results_composition(
            gate_settings,
            project_root=repository_root,
        )
        self._connection: Redis | None = None
        self._run_queue: RqRunQueue | None = None
        self._submitted: list[SubmittedAcceptanceRun] = []
        self._worker_processes: list[subprocess.Popen[str]] = []

    def __enter__(self) -> PlatformAcceptanceHarness:
        self.temporary_root.mkdir(parents=True, exist_ok=True)
        connection = create_api_redis_connection(self.worker_settings)
        try:
            if connection.ping() is not True:
                raise AssertionError("real Redis did not answer PING")
        except Exception:
            connection.close()
            raise AssertionError("real Redis is unavailable for acceptance") from None
        self._connection = connection
        self._run_queue = RqRunQueue(
            self.worker_settings,
            connection=connection,
        )
        return self

    def __exit__(self, _exc_type, _exc_value, _traceback) -> None:
        try:
            self.close()
        except Exception:
            if _exc_type is None:
                raise

    def execute(self, fixture: AcceptanceFixture) -> AcceptanceEvidence:
        """Validate, preflight, enqueue, execute, reopen, and collect one run."""
        submitted = self.submit(fixture)
        self.run_worker()
        return self.collect(submitted)

    def submit(self, fixture: AcceptanceFixture) -> SubmittedAcceptanceRun:
        """Persist and enqueue one run without starting a worker.

        Keeping submission separate gives cancellation tests a durable point at
        which to start and observe a real worker before requesting a stop.
        """
        return self._submit(fixture)

    def run_worker(self) -> None:
        """Run and reap one bounded burst worker for the isolated queue."""
        self.wait_worker(self.start_worker())

    def start_worker(self) -> subprocess.Popen[str]:
        """Start a real worker session that a cancellation test may observe."""
        process = subprocess.Popen(
            (
                sys.executable,
                "-m",
                "bulk_rnaseq_real_execution.worker_entry",
            ),
            cwd=self.repository_root,
            env=self._worker_environment(),
            stdin=subprocess.DEVNULL,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            text=True,
            start_new_session=True,
        )
        self._worker_processes.append(process)
        return process

    def wait_for_execution_activity(
        self,
        submitted: SubmittedAcceptanceRun,
        process: subprocess.Popen[str],
        *,
        require_managed_container: bool,
        timeout_seconds: float = 900,
    ) -> ExecutionActivityEvidence:
        """Observe a claimed RQ job, real Nextflow, and optionally one container."""
        if submitted not in self._submitted:
            raise ValueError("submitted run is not owned by this acceptance harness")
        if process not in self._worker_processes:
            raise ValueError("process is not owned by this acceptance harness")
        if (
            isinstance(timeout_seconds, bool)
            or not isinstance(timeout_seconds, (int, float))
            or timeout_seconds <= 0
            or timeout_seconds > _MAX_ACCEPTANCE_TIMEOUT_SECONDS
        ):
            raise ValueError("activity timeout is outside the acceptance bound")

        run_queue = self._require_queue()
        job = run_queue._queue.fetch_job(submitted.job_id)
        if job is None:
            raise AssertionError("accepted RQ job disappeared")
        cleaner = ManagedContainerCleaner(
            executable=self.gate_settings.docker_executable,
            unix_socket=self.gate_settings.docker_socket,
        )
        if cleaner.verify_endpoint().is_failure:
            raise AssertionError("managed Docker endpoint changed during acceptance")
        scope = managed_container_scope(self.workspace_root / submitted.run_id)
        deadline = time.monotonic() + float(timeout_seconds)
        nextflow_observed = False
        container_observed = False

        persistence = open_run_persistence(self.database_url)
        try:
            from encode_pipeline.services.runs import RunService

            run_service = RunService(
                self.composition.registry,
                repository=persistence.repository,
            )
            while time.monotonic() < deadline:
                if process.poll() is not None:
                    raise AssertionError(
                        "DurableWorker exited before execution activity was observed"
                    )
                record = run_service.get_run(submitted.run_id)
                assignment = run_service.get_execution_assignment(submitted.run_id)
                job.refresh()
                if record.status.is_terminal:
                    raise AssertionError(
                        "bulk RNA-seq run became terminal before required activity"
                    )
                if (
                    record.status is RunStatus.RUNNING
                    and assignment is not None
                    and assignment.dispatched_at is not None
                    and assignment.claimed_at is not None
                    and job.get_status(refresh=True) is JobStatus.STARTED
                    and bool(job.worker_name)
                ):
                    nextflow_observed = nextflow_observed or bool(
                        _worker_session_nextflow_processes(
                            process.pid,
                            runtime_root=self.gate_settings.runtime_root,
                        )
                    )
                    if nextflow_observed and require_managed_container:
                        container_observed = container_observed or bool(
                            managed_container_ids(cleaner, scope, all_containers=False)
                        )
                    if nextflow_observed and (
                        container_observed or not require_managed_container
                    ):
                        process_groups = _worker_session_process_groups(process.pid)
                        if not process_groups:
                            raise AssertionError(
                                "real execution has no auditable worker process group"
                            )
                        return ExecutionActivityEvidence(
                            worker_session_id=process.pid,
                            process_group_count=len(process_groups),
                            nextflow_observed=True,
                            managed_container_observed=container_observed,
                        )
                time.sleep(_ACTIVITY_POLL_SECONDS)
        finally:
            persistence.close()
        raise AssertionError("required real execution activity was not observed")

    def request_cancellation(
        self,
        submitted: SubmittedAcceptanceRun,
        *,
        reason: str,
    ) -> RunCancellationResult:
        """Persist cancellation intent and publish the real RQ stop command."""
        if submitted not in self._submitted:
            raise ValueError("submitted run is not owned by this acceptance harness")
        persistence = open_run_persistence(self.database_url)
        try:
            from encode_pipeline.services.runs import RunService

            run_service = RunService(
                self.composition.registry,
                repository=persistence.repository,
            )
            return RunCancellationService(
                run_service,
                self._require_queue(),
            ).cancel_run(submitted.run_id, reason=reason)
        finally:
            persistence.close()

    def wait_worker(
        self,
        process: subprocess.Popen[str],
        *,
        timeout_seconds: float | None = None,
    ) -> None:
        """Wait for a harness-owned worker and enforce its hard time bound."""
        if process not in self._worker_processes:
            raise ValueError("process is not owned by this acceptance harness")
        if timeout_seconds is not None and (
            isinstance(timeout_seconds, bool)
            or not isinstance(timeout_seconds, (int, float))
            or timeout_seconds <= 0
            or timeout_seconds > _MAX_ACCEPTANCE_TIMEOUT_SECONDS
        ):
            raise ValueError("worker wait timeout is outside the acceptance bound")
        wait_timeout = (
            self.worker_settings.job_timeout_seconds + 90
            if timeout_seconds is None
            else float(timeout_seconds)
        )
        try:
            process.communicate(timeout=wait_timeout)
        except subprocess.TimeoutExpired:
            terminate_rq_worker(process)
            raise AssertionError("bulk RNA-seq DurableWorker timed out") from None
        finally:
            if process.poll() is None:
                terminate_rq_worker(process)
            self._worker_processes.remove(process)
        _assert_worker_session_reaped(process.pid)
        if process.returncode != 0:
            raise AssertionError("bulk RNA-seq DurableWorker exited unsuccessfully")

    def collect(
        self,
        submitted: SubmittedAcceptanceRun,
    ) -> AcceptanceEvidence:
        """Reopen SQLite and collect path-free evidence for one finished job."""
        if submitted not in self._submitted:
            raise ValueError("submitted run is not owned by this acceptance harness")
        run_queue = self._require_queue()
        job = run_queue._queue.fetch_job(submitted.job_id)
        if job is None:
            raise AssertionError("accepted RQ job disappeared")
        job.refresh()
        if not job.is_finished:
            raise AssertionError("accepted RQ job did not finish successfully")
        fixture = load_acceptance_fixture(self.gate_settings.fixture_manifest)
        if (
            fixture.acceptance_manifest_sha256
            != submitted.fixture_acceptance_manifest_sha256
        ):
            raise AssertionError("accepted fixture changed after durable submission")

        persistence = open_run_persistence(self.database_url)
        try:
            from encode_pipeline.services.runs import RunService

            run_service = RunService(
                self.composition.registry,
                repository=persistence.repository,
            )
            cleaner = ManagedContainerCleaner(
                executable=self.gate_settings.docker_executable,
                unix_socket=self.gate_settings.docker_socket,
            )
            return collect_success_evidence(
                run_service=run_service,
                run_id=submitted.run_id,
                expected_job_id=submitted.job_id,
                validated_snapshot_id=submitted.validated_snapshot_id,
                fixture=fixture,
                workspace_root=self.workspace_root,
                repository_root=self.repository_root,
                cleaner=cleaner,
            )
        finally:
            persistence.close()

    def collect_terminal(
        self,
        submitted: SubmittedAcceptanceRun,
    ) -> TerminalLifecycleEvidence:
        """Reopen SQLite and audit a cancelled or failed run and its cleanup."""
        if submitted not in self._submitted:
            raise ValueError("submitted run is not owned by this acceptance harness")
        run_queue = self._require_queue()
        job = run_queue._queue.fetch_job(submitted.job_id)
        if job is None:
            raise AssertionError("accepted RQ job disappeared")
        job.refresh()
        rq_status = job.get_status(refresh=True)
        if rq_status not in {JobStatus.FAILED, JobStatus.STOPPED}:
            raise AssertionError("accepted RQ job lacks a non-success terminal state")

        persistence = open_run_persistence(self.database_url)
        try:
            from encode_pipeline.services.runs import RunService

            run_service = RunService(
                self.composition.registry,
                repository=persistence.repository,
            )
            record = run_service.get_run(submitted.run_id)
            if record.status not in {RunStatus.CANCELLED, RunStatus.FAILED}:
                raise AssertionError("accepted lifecycle is not cancelled or failed")
            assignment = run_service.get_execution_assignment(submitted.run_id)
            if assignment is None or assignment.job_id != submitted.job_id:
                raise AssertionError("accepted lifecycle lost durable RQ ownership")
            events = run_service.list_events(submitted.run_id, limit=1000)
            state = run_service.get_result_state(submitted.run_id)
            artifacts = run_service.list_artifacts(submitted.run_id)
            metrics = run_service.list_qc_metrics(submitted.run_id)

            cleaner = ManagedContainerCleaner(
                executable=self.gate_settings.docker_executable,
                unix_socket=self.gate_settings.docker_socket,
            )
            scope = managed_container_scope(self.workspace_root / submitted.run_id)
            assert_no_managed_containers(cleaner, scope)

            lifecycle_history = tuple(
                event.status.value
                for event in events
                if event.status is not None
                and (event.event_type == "status_changed" or event.status.is_terminal)
            )
            error_reason_code = None
            if record.error is not None:
                candidate = record.error.context.get("reason_code")
                if isinstance(candidate, str):
                    error_reason_code = candidate
            return TerminalLifecycleEvidence(
                run_id=submitted.run_id,
                job_id=submitted.job_id,
                lifecycle_status=record.status.value,
                lifecycle_history=lifecycle_history,
                event_types=tuple(event.event_type for event in events),
                assignment_dispatched=assignment.dispatched_at is not None,
                assignment_claimed=assignment.claimed_at is not None,
                cancellation_requested=(
                    assignment.cancellation_requested_at is not None
                ),
                cancellation_acknowledged=(
                    assignment.cancellation_acknowledged_at is not None
                ),
                cancellation_reason=record.cancellation_reason,
                error_code=(record.error.code if record.error is not None else None),
                error_reason_code=error_reason_code,
                artifact_revision=state.artifact_revision,
                artifact_attempt_id=state.artifact_attempt_id,
                artifact_attempt_status=state.artifact_attempt_status,
                qc_revision=state.qc_revision,
                qc_attempt_id=state.qc_attempt_id,
                qc_attempt_status=state.qc_attempt_status,
                artifact_count=len(artifacts),
                qc_metric_count=len(metrics),
                rq_status=rq_status.value,
                rq_failed=job.is_failed,
                rq_stopped=job.is_stopped,
                rq_finished=job.is_finished,
                cleanup_confirmed=True,
            )
        finally:
            persistence.close()

    def close(self) -> None:
        """Clean only this harness's run scopes, jobs, and unique queue."""
        for process in tuple(self._worker_processes):
            terminate_rq_worker(process)
            self._worker_processes.remove(process)
        cleaner = ManagedContainerCleaner(
            executable=self.gate_settings.docker_executable,
            unix_socket=self.gate_settings.docker_socket,
        )
        cleanup_confirmed = True
        for submitted in self._submitted:
            scope = managed_container_scope(self.workspace_root / submitted.run_id)
            if cleaner.cleanup(scope).is_failure:
                cleanup_confirmed = False
        if self._run_queue is not None:
            for submitted in self._submitted:
                job = self._run_queue._queue.fetch_job(submitted.job_id)
                if job is not None:
                    job.delete()
            self._run_queue._queue.delete()
        if self._connection is not None:
            self._connection.close()
        self._run_queue = None
        self._connection = None
        if not cleanup_confirmed:
            raise AssertionError("managed container cleanup could not be confirmed")

    def _submit(self, fixture: AcceptanceFixture) -> SubmittedAcceptanceRun:
        if not isinstance(fixture, AcceptanceFixture):
            raise ValueError("fixture must be AcceptanceFixture")
        current_fixture = load_acceptance_fixture(self.gate_settings.fixture_manifest)
        if current_fixture != fixture:
            raise AssertionError("submitted fixture differs from its canonical closure")
        run_queue = self._require_queue()
        process_runner = build_acceptance_process_runner(
            settings=self.gate_settings,
            binding=self.composition.binding,
            timeout_seconds=self.worker_settings.job_timeout_seconds,
            passthrough_exceptions=(WorkerHardTimeout,),
        )
        with open_worker_runtime(
            self.worker_settings,
            registry=self.composition.registry,
            build_identity_provider=self.composition.build_identity_provider,
            process_runner=process_runner,
        ) as runtime:
            validation_service = ValidationService(registry=self.composition.registry)
            snapshot_result = ValidatedInputService(
                registry=self.composition.registry,
                validation_service=validation_service,
                build_identity_provider=self.composition.build_identity_provider,
                repository=runtime.persistence.repository,
            ).validate("bulk-rnaseq", fixture.workflow_inputs)
            if snapshot_result.is_failure or snapshot_result.value is None:
                raise AssertionError(
                    "bulk RNA-seq validation failed: "
                    + _issue_codes(snapshot_result.issues)
                )
            snapshot = snapshot_result.value
            created = ValidatedRunCreationService(
                run_service=runtime.run_service,
                build_identity_provider=self.composition.build_identity_provider,
            ).create_run("bulk-rnaseq", snapshot.snapshot_id)
            run_id = created.record.run_id
            preflight = runtime.preflight_service.preflight(run_id)
            if preflight.is_failure:
                raise AssertionError(
                    "bulk RNA-seq preflight failed: " + _issue_codes(preflight.issues)
                )
            queued = RunSubmissionService(
                run_service=runtime.run_service,
                run_queue=run_queue,
                build_identity_provider=runtime.build_identity_provider,
            ).start_run(run_id)
            if queued.status.value != "queued":
                raise AssertionError("bulk RNA-seq run was not durably queued")
            assignment = runtime.run_service.get_execution_assignment(run_id)
            if assignment is None:
                raise AssertionError("bulk RNA-seq run lacks an RQ assignment")
        submitted = SubmittedAcceptanceRun(
            run_id=run_id,
            job_id=assignment.job_id,
            validated_snapshot_id=snapshot.snapshot_id,
            fixture_acceptance_manifest_sha256=(fixture.acceptance_manifest_sha256),
        )
        self._submitted.append(submitted)
        return submitted

    def _worker_environment(self) -> dict[str, str]:
        environment = dict(os.environ)
        environment.update(
            {
                REQUIRE_REAL_EXECUTION_ENV: "1",
                RUNTIME_ROOT_ENV: str(self.gate_settings.runtime_root),
                FIXTURE_MANIFEST_ENV: str(self.gate_settings.fixture_manifest),
                TEST_REDIS_URL_ENV: self.gate_settings.redis_url,
                GATE_DOCKER_EXECUTABLE_ENV: str(self.gate_settings.docker_executable),
                GATE_DOCKER_SOCKET_ENV: str(self.gate_settings.docker_socket),
                DATABASE_URL_ENV: self.worker_settings.database_url,
                REDIS_URL_ENV: self.worker_settings.redis_url,
                QUEUE_NAME_ENV: self.worker_settings.queue_name,
                WORKSPACE_ROOT_ENV: str(self.worker_settings.workspace_root),
                JOB_TIMEOUT_SECONDS_ENV: str(self.worker_settings.job_timeout_seconds),
                MANAGED_DOCKER_EXECUTABLE_ENV: str(
                    self.worker_settings.managed_docker_executable
                ),
                MANAGED_DOCKER_SOCKET_ENV: str(
                    self.worker_settings.managed_docker_socket
                ),
                "PYTHONDONTWRITEBYTECODE": "1",
                "PYTHONPATH": os.pathsep.join(
                    (
                        str(self.repository_root / "src"),
                        str(self.repository_root / "test"),
                    )
                ),
            }
        )
        return environment

    def _require_queue(self) -> RqRunQueue:
        if self._run_queue is None:
            raise RuntimeError("acceptance harness is not open")
        return self._run_queue


def _issue_codes(issues) -> str:
    codes = sorted({getattr(issue, "code", "UNKNOWN") for issue in issues})
    return ",".join(codes) if codes else "UNKNOWN"


def _assert_worker_session_reaped(session_id: int) -> None:
    """Fail after cleaning any process group left by this harness's worker."""
    residual = _worker_session_process_groups(session_id)
    if not residual:
        return
    for process_group in residual:
        try:
            os.killpg(process_group, signal.SIGKILL)
        except ProcessLookupError:
            pass
    deadline = time.monotonic() + 5
    while time.monotonic() < deadline:
        if not _worker_session_process_groups(session_id):
            break
        time.sleep(0.05)
    raise AssertionError("bulk RNA-seq worker left a residual process group")


def _worker_session_process_groups(session_id: int) -> tuple[int, ...]:
    process_groups: set[int] = set()
    try:
        candidates = tuple(Path("/proc").iterdir())
    except OSError:
        raise AssertionError("worker process cleanup cannot be audited") from None
    for candidate in candidates:
        if not candidate.name.isdigit():
            continue
        try:
            pid = int(candidate.name)
            if os.getsid(pid) != session_id:
                continue
            process_group = os.getpgid(pid)
        except (OSError, ValueError):
            continue
        if process_group > 0 and process_group != os.getpgrp():
            process_groups.add(process_group)
    return tuple(sorted(process_groups))


def _worker_session_nextflow_processes(
    session_id: int,
    *,
    runtime_root: Path,
) -> tuple[int, ...]:
    if not runtime_root.is_absolute():
        raise ValueError("runtime_root must be absolute")
    runtime_token = os.fsencode(str(runtime_root))
    process_ids: list[int] = []
    try:
        candidates = tuple(Path("/proc").iterdir())
    except OSError:
        raise AssertionError("Nextflow process activity cannot be audited") from None
    for candidate in candidates:
        if not candidate.name.isdigit():
            continue
        try:
            pid = int(candidate.name)
            if os.getsid(pid) != session_id:
                continue
            with (candidate / "cmdline").open("rb") as stream:
                command = stream.read(1_048_577)
            if len(command) > 1_048_576:
                raise AssertionError("Nextflow process command exceeded audit bound")
        except (OSError, ValueError):
            continue
        if runtime_token in command and b"nextflow" in command.lower():
            process_ids.append(pid)
    return tuple(sorted(process_ids))
