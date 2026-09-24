"""Real SQLite/Redis/RQ orchestration for full bulk RNA-seq acceptance."""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
import json
import os
from pathlib import Path
import re
import signal
import stat
import subprocess
import sys
import time
from uuid import uuid4

from redis import Redis
from rq.job import JobStatus

from encode_pipeline.adapters.bulk_rnaseq.reference_profiles import (
    BULK_RNASEQ_REFERENCE_BINDING_CONTRACT,
)
from encode_pipeline.persistence import (
    DATABASE_URL_ENV,
    open_existing_run_persistence,
    open_run_persistence,
)
from encode_pipeline.platform.adapters import WorkflowInputs
from encode_pipeline.platform.managed_containers import managed_container_scope
from encode_pipeline.platform.runs import RunStatus
from encode_pipeline.services.managed_containers import ManagedContainerCleaner
from encode_pipeline.services.private_reference_profiles import (
    PRIVATE_REFERENCE_PROFILE_SCHEMA_VERSION,
    load_private_reference_profile_config,
)
from encode_pipeline.services.reference_profiles import ReferenceProfileService
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
    rq_job_timeout_seconds,
)
from encode_pipeline.workers.runtime import open_worker_runtime
from encode_pipeline.workers.settings import (
    JOB_TIMEOUT_SECONDS_ENV,
    MANAGED_DOCKER_EXECUTABLE_ENV,
    MANAGED_DOCKER_SOCKET_ENV,
    QUEUE_NAME_ENV,
    REFERENCE_PROFILE_CONFIG_ENV,
    REDIS_URL_ENV,
    WORKSPACE_ROOT_ENV,
    load_worker_settings,
)
from encode_pipeline.workers.timeouts import WorkerHardTimeout

from workers.process_helpers import terminate_rq_worker

from .cancellation_diagnostics import (
    SNAPSHOT_SCHEMA_VERSION,
    WorkerStreamCapture,
    resolve_private_diagnostics_root,
    write_cancellation_snapshot,
)
from .failure_diagnostics import preserve_execution_failure
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
    _write_canonical_evidence_document,
)


_DEFAULT_ACCEPTANCE_TIMEOUT_SECONDS = 7_200
_MAX_ACCEPTANCE_TIMEOUT_SECONDS = 14_400
_ACTIVITY_POLL_SECONDS = 0.05
_RQ_TERMINAL_STABILIZATION_SECONDS = 5.0
_RQ_TERMINAL_POLL_SECONDS = 0.05
_EARLY_TERMINAL_WORKER_WAIT_SECONDS = 120.0
_TERMINAL_LIFECYCLE_EVIDENCE_SCHEMA_VERSION = "1.0.0"
_TERMINAL_BEFORE_REQUIRED_ACTIVITY = "TERMINAL_BEFORE_REQUIRED_ACTIVITY"
_PATH_FREE_EVIDENCE_TOKEN = re.compile(r"[A-Za-z0-9_.:-]+")
_REFERENCE_CONFIG_DIRECTORY = "operator-reference-profile"
_REFERENCE_CONFIG_FILENAME = "reference-profiles.json"
_REFERENCE_CONFIG_KEY = "bulk-rnaseq-gate-tiny-private"
_REFERENCE_PROFILE_SAFE_KEY = "bulk-rnaseq-gate-tiny"
_REFERENCE_PROFILE_DISPLAY_NAME = "Bulk RNA-seq protected tiny"
_REFERENCE_PROFILE_ORGANISM = "Synthetic organism"
_REFERENCE_PROFILE_ASSEMBLY = "tiny"
_DIAGNOSTICS_REDIS_KEY_LIMIT = 64
_DIAGNOSTICS_PROC_SCAN_LIMIT = 8192
_DIAGNOSTICS_CMDLINE_LIMIT = 4096


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
    assertion_reason_code: str | None = None

    def to_dict(self) -> dict[str, object]:
        """Return the deliberately path-free Protected Gate projection."""
        return {
            "schema_version": _TERMINAL_LIFECYCLE_EVIDENCE_SCHEMA_VERSION,
            "assertion_reason_code": _path_free_evidence_token(
                self.assertion_reason_code
            ),
            "run_id": _path_free_evidence_token(self.run_id),
            "job_id": _path_free_evidence_token(self.job_id),
            "lifecycle_status": _path_free_evidence_token(self.lifecycle_status),
            "lifecycle_history": [
                _path_free_evidence_token(value) for value in self.lifecycle_history
            ],
            "event_types": [
                _path_free_evidence_token(value) for value in self.event_types
            ],
            "assignment_dispatched": self.assignment_dispatched,
            "assignment_claimed": self.assignment_claimed,
            "cancellation_requested": self.cancellation_requested,
            "cancellation_acknowledged": self.cancellation_acknowledged,
            "error_code": _path_free_evidence_token(self.error_code),
            "error_reason_code": _path_free_evidence_token(self.error_reason_code),
            "artifact_revision": self.artifact_revision,
            "artifact_attempt_id": _path_free_evidence_token(self.artifact_attempt_id),
            "artifact_attempt_status": _path_free_evidence_token(
                self.artifact_attempt_status
            ),
            "qc_revision": self.qc_revision,
            "qc_attempt_id": _path_free_evidence_token(self.qc_attempt_id),
            "qc_attempt_status": _path_free_evidence_token(self.qc_attempt_status),
            "artifact_count": self.artifact_count,
            "qc_metric_count": self.qc_metric_count,
            "rq_status": _path_free_evidence_token(self.rq_status),
            "rq_failed": self.rq_failed,
            "rq_stopped": self.rq_stopped,
            "rq_finished": self.rq_finished,
            "cleanup_confirmed": self.cleanup_confirmed,
        }


def prepare_acceptance_database(database_url: str) -> None:
    """Migrate the harness-owned fresh database exactly once.

    This is the harness's operator-equivalent preparation step. It must run
    before any worker runtime or existing-only service composition opens the
    database; every later harness access uses ``open_existing_run_persistence``
    so a missing preparation fails closed instead of being masked by an
    implicit migration.
    """
    persistence = open_run_persistence(database_url)
    persistence.close()


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
        self.reference_profile_config_path = (
            temporary_root / _REFERENCE_CONFIG_DIRECTORY / _REFERENCE_CONFIG_FILENAME
        )
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
                REFERENCE_PROFILE_CONFIG_ENV: str(self.reference_profile_config_path),
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
        self._worker_streams: list[
            tuple[subprocess.Popen[str], WorkerStreamCapture]
        ] = []
        self._private_diagnostics_root: Path | None = None
        self._private_diagnostics_error: str | None = None
        self._reference_profile_revision_id: str | None = None
        self._reference_profile_public_identity_sha256: str | None = None
        self._reference_profile_directory_identity: tuple[int, int] | None = None
        self._reference_profile_config_identity: tuple[int, int] | None = None

    def __enter__(self) -> PlatformAcceptanceHarness:
        self.temporary_root.mkdir(parents=True, exist_ok=True)
        prepare_acceptance_database(self.database_url)
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
        if _exc_type is not None:
            for submitted in self._submitted:
                if not (
                    self.temporary_root / "evidence/execution-failure.json"
                ).exists():
                    self._preserve_failure(submitted)
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
        """Start a real worker session that a cancellation test may observe.

        The session's raw stdout/stderr are kept in the owner-only private
        diagnostics directory rather than discarded, because a cancellation
        that never reaches a terminal state must still be explainable. They stay
        out of the public evidence tree, the download bundle, and mail.
        """
        root = self._resolve_private_diagnostics_root()
        capture = (
            WorkerStreamCapture(root=root, ordinal=len(self._worker_streams))
            if root is not None
            else None
        )
        try:
            process = subprocess.Popen(
                (
                    sys.executable,
                    "-m",
                    "bulk_rnaseq_real_execution.worker_entry",
                ),
                cwd=self.repository_root,
                env=self._worker_environment(),
                stdin=subprocess.DEVNULL,
                stdout=(
                    capture.stdout_handle if capture is not None else subprocess.DEVNULL
                ),
                stderr=(
                    capture.stderr_handle if capture is not None else subprocess.DEVNULL
                ),
                text=True,
                start_new_session=True,
            )
        except Exception:
            if capture is not None:
                capture.close(returncode=None)
            raise
        if capture is not None:
            self._worker_streams.append((process, capture))
        self._worker_processes.append(process)
        return process

    def _close_worker_streams(
        self,
        process: subprocess.Popen[str],
        *,
        returncode: int | None,
    ) -> None:
        for index, (owner, capture) in enumerate(tuple(self._worker_streams)):
            if owner is process:
                del self._worker_streams[index]
                capture.close(returncode=returncode)
                return

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
        scope = managed_container_scope(self.workspace_root / submitted.run_id)
        deadline = time.monotonic() + float(timeout_seconds)
        nextflow_observed = False
        container_observed = False
        terminal_observed = False
        endpoint_verified = False
        job = None
        cleaner = None

        persistence = open_existing_run_persistence(self.database_url)
        try:
            from encode_pipeline.services.runs import RunService

            run_service = RunService(
                self.composition.registry,
                repository=persistence.repository,
            )
            while time.monotonic() < deadline:
                record = run_service.get_run(submitted.run_id)
                assignment = run_service.get_execution_assignment(submitted.run_id)
                if record.status.is_terminal:
                    terminal_observed = True
                    break
                if not endpoint_verified:
                    try:
                        cleaner = ManagedContainerCleaner(
                            executable=self.gate_settings.docker_executable,
                            unix_socket=self.gate_settings.docker_socket,
                        )
                        endpoint_result = cleaner.verify_endpoint()
                    except Exception:
                        raise AssertionError(
                            "managed Docker endpoint is unavailable"
                        ) from None
                    if endpoint_result.is_failure:
                        raise AssertionError(
                            "managed Docker endpoint changed during acceptance"
                        )
                    endpoint_verified = True
                if job is None:
                    try:
                        job = run_queue._queue.fetch_job(submitted.job_id)
                    except Exception:
                        raise AssertionError(
                            "accepted RQ job state is unavailable"
                        ) from None
                    if job is None:
                        raise AssertionError("accepted RQ job state is unavailable")
                try:
                    job.refresh()
                    worker_name = job.worker_name
                    rq_status = job.get_status(refresh=False)
                except Exception:
                    raise AssertionError(
                        "accepted RQ job state is unavailable"
                    ) from None
                if process.poll() is not None:
                    raise AssertionError(
                        "DurableWorker exited before execution activity was observed"
                    )
                if (
                    record.status is RunStatus.RUNNING
                    and assignment is not None
                    and assignment.dispatched_at is not None
                    and assignment.claimed_at is not None
                    and rq_status is JobStatus.STARTED
                    and bool(worker_name)
                ):
                    nextflow_observed = nextflow_observed or bool(
                        _worker_session_nextflow_processes(
                            process.pid,
                            runtime_root=self.gate_settings.runtime_root,
                        )
                    )
                    if nextflow_observed and require_managed_container:
                        assert cleaner is not None
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
        if terminal_observed:
            evidence = self.collect_terminal(
                submitted,
                assertion_reason_code=_TERMINAL_BEFORE_REQUIRED_ACTIVITY,
                allow_unstable_rq=True,
            )
            message = _terminal_before_activity_message(evidence)
            try:
                self.wait_worker(
                    process,
                    timeout_seconds=_EARLY_TERMINAL_WORKER_WAIT_SECONDS,
                )
            except AssertionError:
                raise AssertionError(f"{message} worker_reap=FAILED") from None
            raise AssertionError(message)
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
        persistence = open_existing_run_persistence(self.database_url)
        try:
            from encode_pipeline.services.runs import RunService

            run_service = RunService(
                self.composition.registry,
                repository=persistence.repository,
            )
            result = RunCancellationService(
                run_service,
                self._require_queue(),
            ).cancel_run(submitted.run_id, reason=reason)
            self._capture_diagnostics_snapshot(
                "cancellation-requested",
                note=f"stop_requested={result.stop_requested}",
            )
            return result
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
            rq_job_timeout_seconds(self.worker_settings.job_timeout_seconds) + 60
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
            self._close_worker_streams(process, returncode=process.returncode)
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
            # Collect while RQ metadata still exists, before close deletes it.
            try:
                self.collect_terminal(submitted)
            except Exception:
                self._preserve_failure(submitted)
            raise AssertionError("accepted RQ job did not finish successfully")
        fixture = load_acceptance_fixture(self.gate_settings.fixture_manifest)
        if (
            fixture.acceptance_manifest_sha256
            != submitted.fixture_acceptance_manifest_sha256
        ):
            raise AssertionError("accepted fixture changed after durable submission")

        persistence = open_existing_run_persistence(self.database_url)
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
        *,
        assertion_reason_code: str | None = None,
        allow_unstable_rq: bool = False,
    ) -> TerminalLifecycleEvidence:
        """Audit a non-success run or an early-success diagnostic and cleanup."""
        if submitted not in self._submitted:
            raise ValueError("submitted run is not owned by this acceptance harness")
        run_queue = self._require_queue()
        try:
            job = run_queue._queue.fetch_job(submitted.job_id)
        except Exception:
            if not allow_unstable_rq:
                raise AssertionError("accepted RQ job state is unavailable") from None
            job = None
        try:
            rq_status = _wait_for_rq_terminal_status(job) if job is not None else None
        except Exception:
            if not allow_unstable_rq:
                raise AssertionError("accepted RQ job state is unavailable") from None
            rq_status = None
        if (
            rq_status not in {JobStatus.FAILED, JobStatus.STOPPED}
            and not allow_unstable_rq
        ):
            self._capture_diagnostics_snapshot(
                "rq-not-terminal",
                note=(
                    "rq_status="
                    + (
                        rq_status.value
                        if isinstance(rq_status, JobStatus)
                        else "unavailable"
                    )
                ),
            )
            raise AssertionError("accepted RQ job lacks a non-success terminal state")

        persistence = open_existing_run_persistence(self.database_url)
        try:
            from encode_pipeline.services.runs import RunService

            run_service = RunService(
                self.composition.registry,
                repository=persistence.repository,
            )
            record = run_service.get_run(submitted.run_id)
            expected_terminal_statuses = {RunStatus.CANCELLED, RunStatus.FAILED}
            if assertion_reason_code == _TERMINAL_BEFORE_REQUIRED_ACTIVITY:
                expected_terminal_statuses.add(RunStatus.SUCCEEDED)
            if record.status not in expected_terminal_statuses:
                raise AssertionError("accepted lifecycle terminal state is invalid")
            assignment = run_service.get_execution_assignment(submitted.run_id)
            if assignment is None or assignment.job_id != submitted.job_id:
                raise AssertionError("accepted lifecycle lost durable RQ ownership")
            events = run_service.list_events(submitted.run_id, limit=1000)
            state = run_service.get_result_state(submitted.run_id)
            artifacts = run_service.list_artifacts(submitted.run_id)
            metrics = run_service.list_qc_metrics(submitted.run_id)

            cleanup_confirmed = True
            try:
                cleaner = ManagedContainerCleaner(
                    executable=self.gate_settings.docker_executable,
                    unix_socket=self.gate_settings.docker_socket,
                )
                scope = managed_container_scope(self.workspace_root / submitted.run_id)
                assert_no_managed_containers(cleaner, scope)
            except Exception:
                cleanup_confirmed = False

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
            evidence = TerminalLifecycleEvidence(
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
                rq_status=(
                    rq_status.value
                    if isinstance(rq_status, JobStatus)
                    else "unavailable"
                ),
                rq_failed=rq_status is JobStatus.FAILED,
                rq_stopped=rq_status is JobStatus.STOPPED,
                rq_finished=rq_status is JobStatus.FINISHED,
                cleanup_confirmed=cleanup_confirmed,
                assertion_reason_code=assertion_reason_code,
            )
            self._publish_terminal_lifecycle_evidence(evidence)
            if record.status == RunStatus.FAILED:
                self._preserve_failure(submitted, reason_code=error_reason_code)
            if not cleanup_confirmed and not allow_unstable_rq:
                self._capture_diagnostics_snapshot(
                    "cleanup-incomplete",
                    note="cleanup_confirmed=false",
                )
                raise AssertionError("accepted lifecycle cleanup is incomplete")
            return evidence
        finally:
            persistence.close()

    def _preserve_failure(
        self, submitted: SubmittedAcceptanceRun, *, reason_code: str | None = None
    ) -> None:
        try:
            job = self._require_queue()._queue.fetch_job(submitted.job_id)
            rq_exception = getattr(job, "exc_info", None)
        except Exception:
            rq_exception = None
        preserve_execution_failure(
            workspace=self.workspace_root / submitted.run_id,
            evidence_root=self.temporary_root / "evidence",
            stage="platform",
            reason_code=reason_code,
            rq_exception=rq_exception if isinstance(rq_exception, str) else None,
        )

    def _publish_terminal_lifecycle_evidence(
        self,
        evidence: TerminalLifecycleEvidence,
    ) -> None:
        run_id = _path_free_evidence_token(evidence.run_id)
        if run_id in {None, "REDACTED"}:
            raise AssertionError("terminal lifecycle evidence identity is invalid")
        _write_canonical_evidence_document(
            evidence.to_dict(),
            (self.temporary_root / "evidence" / f"terminal-lifecycle-{run_id}.json"),
            failure_message="terminal lifecycle evidence could not be published",
        )

    def close(self) -> None:
        """Clean only this harness's run scopes, jobs, and unique queue."""
        # Capture the decisive private state before cleanup destroys it.
        self._capture_diagnostics_snapshot("pre-close")
        for process in tuple(self._worker_processes):
            terminate_rq_worker(process)
            self._worker_processes.remove(process)
            self._close_worker_streams(process, returncode=process.returncode)
        for _owner, capture in tuple(self._worker_streams):
            capture.close(returncode=None)
        self._worker_streams.clear()
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
        reference_config_cleanup_confirmed = self._cleanup_reference_profile_config()
        if not cleanup_confirmed or not reference_config_cleanup_confirmed:
            raise AssertionError("acceptance cleanup could not be confirmed")

    def _submit(self, fixture: AcceptanceFixture) -> SubmittedAcceptanceRun:
        if not isinstance(fixture, AcceptanceFixture):
            raise ValueError("fixture must be AcceptanceFixture")
        current_fixture = load_acceptance_fixture(self.gate_settings.fixture_manifest)
        if current_fixture != fixture:
            raise AssertionError("submitted fixture differs from its canonical closure")
        self._prepare_reference_profile_config(fixture)
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
            reference_profiles = ReferenceProfileService(
                repository=runtime.persistence.reference_profile_repository,
                private_config_provider=self._load_private_reference_profile_config,
                adapter_provider=self.composition.registry.get,
            )
            revision_id = self._ensure_reference_profile(reference_profiles)
            public_inputs = _public_fixture_inputs(fixture)
            validation_service = ValidationService(registry=self.composition.registry)
            snapshot_result = ValidatedInputService(
                registry=self.composition.registry,
                validation_service=validation_service,
                build_identity_provider=self.composition.build_identity_provider,
                repository=runtime.persistence.repository,
                reference_profile_binding_service=(
                    runtime.reference_profile_binding_service
                ),
                reference_profile_catalog=reference_profiles,
            ).validate(
                "bulk-rnaseq",
                public_inputs,
                reference_profile_revision_id=revision_id,
            )
            if snapshot_result.is_failure or snapshot_result.value is None:
                raise AssertionError(
                    "bulk RNA-seq validation failed: "
                    + _issue_codes(snapshot_result.issues)
                )
            snapshot = snapshot_result.value
            created = ValidatedRunCreationService(
                run_service=runtime.run_service,
                build_identity_provider=self.composition.build_identity_provider,
                reference_profile_binding_service=(
                    runtime.reference_profile_binding_service
                ),
            ).create_run("bulk-rnaseq", snapshot.snapshot_id)
            run_id = created.record.run_id
            snapshot_binding = runtime.run_service.get_validated_reference_binding(
                snapshot.snapshot_id
            )
            if (
                snapshot_binding is None
                or snapshot_binding.revision_id != revision_id
                or snapshot_binding.revision_public_identity_sha256
                != self._reference_profile_public_identity_sha256
            ):
                raise AssertionError(
                    "bulk RNA-seq snapshot lost its exact reference identity"
                )
            run_binding = runtime.run_service.get_run_reference_binding(run_id)
            if run_binding != snapshot_binding:
                raise AssertionError(
                    "bulk RNA-seq run differs from its snapshot reference identity"
                )
            preflight = runtime.preflight_service.preflight(run_id)
            if preflight.is_failure:
                raise AssertionError(
                    "bulk RNA-seq preflight failed: " + _issue_codes(preflight.issues)
                )
            queued = RunSubmissionService(
                run_service=runtime.run_service,
                run_queue=run_queue,
                build_identity_provider=runtime.build_identity_provider,
                reference_profile_resolver=runtime.reference_profile_resolver,
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
                REFERENCE_PROFILE_CONFIG_ENV: str(
                    self._require_reference_profile_config_path()
                ),
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

    def _prepare_reference_profile_config(
        self,
        fixture: AcceptanceFixture,
    ) -> None:
        if self._reference_profile_config_identity is not None:
            if (
                self._reference_profile_directory_identity is None
                or not _matches_directory_identity(
                    self.reference_profile_config_path.parent,
                    self._reference_profile_directory_identity,
                )
                or not _matches_regular_file_identity(
                    self.reference_profile_config_path,
                    self._reference_profile_config_identity,
                    expected_mode=0o600,
                )
            ):
                raise AssertionError("private reference profile config changed")
            return
        config_directory = self.reference_profile_config_path.parent
        descriptor = -1
        try:
            self.temporary_root.mkdir(parents=True, exist_ok=True)
            config_directory.mkdir(mode=0o700)
            directory_stat = config_directory.stat(follow_symlinks=False)
            directory_identity = (directory_stat.st_dev, directory_stat.st_ino)
            self._reference_profile_directory_identity = directory_identity
            if not _matches_directory_identity(
                config_directory,
                directory_identity,
            ):
                raise AssertionError("private reference profile directory is invalid")
            document = _private_reference_profile_document(fixture)
            descriptor = os.open(
                self.reference_profile_config_path,
                os.O_WRONLY
                | os.O_CREAT
                | os.O_EXCL
                | getattr(os, "O_CLOEXEC", 0)
                | getattr(os, "O_NOFOLLOW", 0),
                0o600,
            )
            created_stat = os.fstat(descriptor)
            identity = (created_stat.st_dev, created_stat.st_ino)
            if (
                not stat.S_ISREG(created_stat.st_mode)
                or created_stat.st_uid != os.getuid()
                or created_stat.st_nlink != 1
            ):
                raise AssertionError("private reference profile config is invalid")
            self._reference_profile_config_identity = identity
            os.fchmod(descriptor, 0o600)
            with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
                descriptor = -1
                json.dump(
                    document,
                    handle,
                    ensure_ascii=True,
                    sort_keys=True,
                    separators=(",", ":"),
                )
                handle.write("\n")
                handle.flush()
                os.fsync(handle.fileno())
            if not _matches_regular_file_identity(
                self.reference_profile_config_path,
                identity,
                expected_mode=0o600,
            ):
                raise AssertionError("private reference profile config is invalid")
        except Exception:
            if descriptor >= 0:
                os.close(descriptor)
            if self._cleanup_reference_profile_config():
                raise AssertionError(
                    "private reference profile config could not be prepared"
                ) from None
            raise AssertionError(
                "private reference profile config cleanup could not be confirmed"
            ) from None

    def _load_private_reference_profile_config(self):
        return load_private_reference_profile_config(
            self._require_reference_profile_config_path()
        )

    def _ensure_reference_profile(
        self,
        reference_profiles: ReferenceProfileService,
    ) -> str:
        revision_id = self._reference_profile_revision_id
        if revision_id is None:
            registered = reference_profiles.register(
                safe_key=_REFERENCE_PROFILE_SAFE_KEY,
                display_name=_REFERENCE_PROFILE_DISPLAY_NAME,
                organism=_REFERENCE_PROFILE_ORGANISM,
                assembly=_REFERENCE_PROFILE_ASSEMBLY,
                config_key=_REFERENCE_CONFIG_KEY,
            )
            enabled = reference_profiles.enable(
                registered.profile_id,
                revision_id=registered.revision_id,
            )
            if (
                not enabled.enabled
                or enabled.revision_id != registered.revision_id
                or enabled.public_identity_sha256 != registered.public_identity_sha256
            ):
                raise AssertionError("private reference profile was not enabled")
            self._reference_profile_revision_id = enabled.revision_id
            self._reference_profile_public_identity_sha256 = (
                enabled.public_identity_sha256
            )
            return enabled.revision_id
        summary = reference_profiles.get_revision_summary(revision_id)
        if (
            not summary.enabled
            or summary.public_identity_sha256
            != self._reference_profile_public_identity_sha256
        ):
            raise AssertionError("private reference profile identity changed")
        return revision_id

    def _require_reference_profile_config_path(self) -> Path:
        path = self.worker_settings.reference_profile_config
        directory_identity = self._reference_profile_directory_identity
        identity = self._reference_profile_config_identity
        if (
            path is None
            or directory_identity is None
            or identity is None
            or path != self.reference_profile_config_path
            or not _matches_directory_identity(path.parent, directory_identity)
            or not _matches_regular_file_identity(
                path,
                identity,
                expected_mode=0o600,
            )
        ):
            raise AssertionError("private reference profile config is unavailable")
        return path

    def _cleanup_reference_profile_config(self) -> bool:
        config_directory = self.reference_profile_config_path.parent
        directory_identity = self._reference_profile_directory_identity
        identity = self._reference_profile_config_identity
        if directory_identity is None:
            try:
                config_directory.stat(follow_symlinks=False)
            except FileNotFoundError:
                return True
            except OSError:
                return False
            return False
        if not _matches_directory_identity(config_directory, directory_identity):
            return False
        if identity is None:
            try:
                self.reference_profile_config_path.stat(follow_symlinks=False)
            except FileNotFoundError:
                pass
            except OSError:
                return False
            else:
                return False
        elif not _matches_regular_file_identity(
            self.reference_profile_config_path,
            identity,
            expected_mode=0o600,
        ):
            return False
        try:
            if identity is not None:
                self.reference_profile_config_path.unlink()
            config_directory.rmdir()
        except OSError:
            return False
        self._reference_profile_directory_identity = None
        self._reference_profile_config_identity = None
        return True

    def _require_queue(self) -> RqRunQueue:
        if self._run_queue is None:
            raise RuntimeError("acceptance harness is not open")
        return self._run_queue

    def _resolve_private_diagnostics_root(self) -> Path | None:
        """Resolve the owner-only diagnostics root once; never mask a failure."""
        if self._private_diagnostics_root is not None:
            return self._private_diagnostics_root
        if self._private_diagnostics_error is not None:
            return None
        try:
            self._private_diagnostics_root = resolve_private_diagnostics_root(
                self.temporary_root / "evidence"
            )
        except Exception as error:
            self._private_diagnostics_error = type(error).__name__
            return None
        return self._private_diagnostics_root

    def _capture_diagnostics_snapshot(
        self,
        label: str,
        *,
        note: str | None = None,
    ) -> str | None:
        """Record private pre-destruction state; never raise, never publish."""
        root = self._resolve_private_diagnostics_root()
        if root is None:
            return None
        sections: dict[str, object] = {}
        for name, collector in (
            ("redis", self._diagnostics_redis),
            ("rq_job", self._diagnostics_rq_job),
            ("sqlite", self._diagnostics_sqlite),
            ("processes", self._diagnostics_processes),
            ("containers", self._diagnostics_containers),
        ):
            try:
                sections[name] = collector()
            except Exception as error:
                sections[name] = {"section_error": type(error).__name__}
        document: dict[str, object] = {
            "schema": SNAPSHOT_SCHEMA_VERSION,
            "note": note,
            "run_ids": [submitted.run_id for submitted in self._submitted],
            "job_ids": [submitted.job_id for submitted in self._submitted],
            "worker_session_ids": [process.pid for process in self._worker_processes],
            "worker_returncodes": [
                process.returncode for process in self._worker_processes
            ],
            "sections": sections,
        }
        if self._private_diagnostics_error is not None:
            document["root_error"] = self._private_diagnostics_error
        return write_cancellation_snapshot(root=root, label=label, document=document)

    def _diagnostics_redis(self) -> dict[str, object]:
        """Dump every Redis key that names this harness's queue or jobs."""
        connection = self._connection
        if connection is None:
            return {"unavailable": "harness is not open"}
        tokens = [self.queue_name, *(s.job_id for s in self._submitted)]
        keys: set[str] = set()
        for token in tokens:
            try:
                keys.update(
                    _redis_text(key)
                    for key in connection.scan_iter(match=f"*{token}*", count=200)
                )
            except Exception as error:
                return {"scan_error": type(error).__name__}
        ordered = sorted(keys)
        return {
            "key_count": len(ordered),
            "keys": {
                key: _redis_key_state(connection, key)
                for key in ordered[:_DIAGNOSTICS_REDIS_KEY_LIMIT]
            },
        }

    def _diagnostics_rq_job(self) -> dict[str, object]:
        """Read the live RQ job objects for each submitted run."""
        run_queue = self._run_queue
        if run_queue is None:
            return {"unavailable": "harness is not open"}
        document: dict[str, object] = {}
        for submitted in self._submitted:
            entry: dict[str, object] = {}
            try:
                job = run_queue._queue.fetch_job(submitted.job_id)
            except Exception as error:
                document[submitted.job_id] = {"fetch_error": type(error).__name__}
                continue
            if job is None:
                document[submitted.job_id] = {"present": False}
                continue
            entry["present"] = True
            try:
                entry["status_refreshed"] = _enum_value(job.get_status(refresh=True))
            except Exception as error:
                entry["status_error"] = type(error).__name__
            for attribute in (
                "id",
                "origin",
                "description",
                "worker_name",
                "created_at",
                "enqueued_at",
                "started_at",
                "ended_at",
                "timeout",
                "result_ttl",
                "failure_ttl",
                "retries_left",
                "is_failed",
                "is_finished",
                "is_stopped",
                "is_queued",
                "is_started",
                "is_deferred",
                "is_scheduled",
                "is_canceled",
            ):
                entry[attribute] = _job_attribute(job, attribute)
            try:
                entry["exc_info"] = _bounded_text(job.exc_info)
            except Exception as error:
                entry["exc_info_error"] = type(error).__name__
            document[submitted.job_id] = entry
        return document

    def _diagnostics_sqlite(self) -> dict[str, object]:
        """Read the durable run, assignment, event, and result state."""
        persistence = open_existing_run_persistence(self.database_url)
        try:
            from encode_pipeline.services.runs import RunService

            run_service = RunService(
                self.composition.registry,
                repository=persistence.repository,
            )
            document: dict[str, object] = {}
            for submitted in self._submitted:
                entry: dict[str, object] = {
                    "run": _run_record_projection(
                        run_service.get_run(submitted.run_id)
                    ),
                    "assignment": _assignment_projection(
                        run_service.get_execution_assignment(submitted.run_id)
                    ),
                    "events": [
                        _event_projection(event)
                        for event in run_service.list_events(
                            submitted.run_id, limit=1000
                        )
                    ],
                    "result_state": _result_state_projection(
                        run_service.get_result_state(submitted.run_id)
                    ),
                    "artifact_count": len(run_service.list_artifacts(submitted.run_id)),
                    "qc_metric_count": len(
                        run_service.list_qc_metrics(submitted.run_id)
                    ),
                }
                document[submitted.run_id] = entry
            return document
        finally:
            persistence.close()

    def _diagnostics_processes(self) -> dict[str, object]:
        """Record every process still alive in a harness-owned worker session."""
        document: dict[str, object] = {
            "harness_pid": os.getpid(),
            "harness_process_group": os.getpgrp(),
            "workers": [],
        }
        for process in self._worker_processes:
            entry: dict[str, object] = {
                "pid": process.pid,
                "returncode": process.returncode,
            }
            try:
                entry["poll"] = process.poll()
            except Exception as error:
                entry["poll_error"] = type(error).__name__
            for name, reader in (
                ("session_id", lambda: os.getsid(process.pid)),
                (
                    "process_groups",
                    lambda: _worker_session_process_groups(process.pid),
                ),
                (
                    "nextflow_pids",
                    lambda: _worker_session_nextflow_processes(
                        process.pid,
                        runtime_root=self.gate_settings.runtime_root,
                    ),
                ),
                ("table", lambda: _process_table(process.pid)),
            ):
                try:
                    entry[name] = reader()
                except Exception as error:
                    entry[name] = {"error": type(error).__name__}
            document["workers"].append(entry)
        return document

    def _diagnostics_containers(self) -> dict[str, object]:
        """Record the managed container scope before cleanup removes it."""
        document: dict[str, object] = {}
        try:
            cleaner = ManagedContainerCleaner(
                executable=self.gate_settings.docker_executable,
                unix_socket=self.gate_settings.docker_socket,
            )
            document["endpoint_is_failure"] = bool(
                getattr(cleaner.verify_endpoint(), "is_failure", None)
            )
        except Exception as error:
            return {"endpoint_error": type(error).__name__}
        for submitted in self._submitted:
            scope = managed_container_scope(self.workspace_root / submitted.run_id)
            entry: dict[str, object] = {}
            for key, all_containers in (
                ("scope_containers", False),
                ("managed_containers", True),
            ):
                try:
                    entry[key] = [
                        str(value)
                        for value in managed_container_ids(
                            cleaner, scope, all_containers=all_containers
                        )
                    ]
                except Exception as error:
                    entry[key] = {"error": type(error).__name__}
            document[submitted.run_id] = entry
        return document


def _redis_text(value: object) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8", "replace")
    return str(value)


def _redis_key_state(connection, key: str) -> object:
    try:
        kind = _redis_text(connection.type(key))
    except Exception as error:
        return {"type_error": type(error).__name__}
    try:
        if kind == "string":
            return {"type": kind, "value": _bounded_text(connection.get(key))}
        if kind == "hash":
            return {
                "type": kind,
                "fields": {
                    _redis_text(field): _bounded_text(value)
                    for field, value in connection.hgetall(key).items()
                },
            }
        if kind == "list":
            return {
                "type": kind,
                "length": connection.llen(key),
                "values": [
                    _bounded_text(value) for value in connection.lrange(key, 0, 63)
                ],
            }
        if kind == "set":
            return {
                "type": kind,
                "length": connection.scard(key),
                "values": sorted(
                    _bounded_text(value) for value in connection.smembers(key)
                )[:64],
            }
        if kind == "zset":
            return {
                "type": kind,
                "length": connection.zcard(key),
                "values": [
                    _bounded_text(value) for value in connection.zrange(key, 0, 63)
                ],
            }
        return {"type": kind}
    except Exception as error:
        return {"type": kind, "read_error": type(error).__name__}


def _enum_value(value: object) -> object:
    return getattr(value, "value", value)


def _job_attribute(job, name: str) -> object:
    try:
        value = getattr(job, name, None)
    except Exception as error:
        return {"error": type(error).__name__}
    if callable(value):
        try:
            value = value()
        except Exception as error:
            return {"error": type(error).__name__}
    if isinstance(value, (bool, int, float)) or value is None:
        return value
    return _bounded_text(value)


def _bounded_text(value: object) -> object:
    if value is None:
        return None
    if isinstance(value, bytes):
        return value.decode("utf-8", "replace")[:_DIAGNOSTICS_CMDLINE_LIMIT]
    if isinstance(value, str):
        return value[:65536]
    return str(value)[:65536]


def _bounded_value(value: object) -> object:
    if isinstance(value, str):
        return value[:8192]
    if isinstance(value, dict):
        return {
            str(key): _bounded_value(item) for key, item in list(value.items())[:64]
        }
    if isinstance(value, (list, tuple)):
        return [_bounded_value(item) for item in list(value)[:64]]
    if isinstance(value, (bool, int, float)) or value is None:
        return value
    return _bounded_text(value)


def _run_record_projection(record) -> dict[str, object]:
    error = getattr(record, "error", None)
    return {
        "run_id": getattr(record, "run_id", None),
        "workflow_id": getattr(record, "workflow_id", None),
        "status": _enum_value(getattr(record, "status", None)),
        "created_at": getattr(record, "created_at", None),
        "updated_at": getattr(record, "updated_at", None),
        "started_at": getattr(record, "started_at", None),
        "ended_at": getattr(record, "ended_at", None),
        "current_stage": getattr(record, "current_stage", None),
        "cancellation_reason": getattr(record, "cancellation_reason", None),
        "error_code": getattr(error, "code", None),
        "error_context": _bounded_value(getattr(error, "context", None)),
    }


def _assignment_projection(assignment) -> dict[str, object] | None:
    if assignment is None:
        return None
    return {
        name: getattr(assignment, name, None)
        for name in (
            "run_id",
            "job_id",
            "backend",
            "queue_name",
            "created_at",
            "managed_container_scope",
            "dispatched_at",
            "claimed_at",
            "cancellation_requested_at",
            "cancellation_reason",
            "cancellation_acknowledged_at",
            "requeue_requested_at",
            "requeue_confirmed_at",
        )
    }


def _event_projection(event) -> dict[str, object]:
    issue = getattr(event, "issue", None)
    return {
        "sequence": getattr(event, "sequence", None),
        "event_type": getattr(event, "event_type", None),
        "timestamp": getattr(event, "timestamp", None),
        "status": _enum_value(getattr(event, "status", None)),
        "stage": getattr(event, "stage", None),
        "message": _bounded_value(getattr(event, "message", None)),
        "context": _bounded_value(getattr(event, "context", None)),
        "issue_code": getattr(issue, "code", None),
        "issue_message": _bounded_value(getattr(issue, "message", None)),
    }


def _result_state_projection(state) -> dict[str, object] | None:
    if state is None:
        return None
    return {
        name: getattr(state, name, None)
        for name in (
            "artifact_revision",
            "artifact_attempt_id",
            "artifact_attempt_status",
            "qc_revision",
            "qc_attempt_id",
            "qc_attempt_status",
        )
    }


def _parse_proc_stat(raw: str) -> dict[str, object] | None:
    start = raw.find("(")
    end = raw.rfind(")")
    if start == -1 or end < start:
        return None
    pid_text = raw[:start].strip()
    remainder = raw[end + 1 :].split()
    if not pid_text.isdigit() or len(remainder) < 20:
        return None

    def field(index: int) -> int | None:
        return int(remainder[index]) if remainder[index].isdigit() else None

    return {
        "comm": raw[start + 1 : end][:128],
        "state": remainder[0],
        "ppid": field(1),
        "pgrp": field(2),
        "session": field(3),
        "starttime_ticks": field(19),
    }


def _process_table(session_id: int) -> object:
    """List every live process whose session matches one worker session id."""
    try:
        entries = tuple(Path("/proc").iterdir())
    except OSError:
        return {"error": "OSError"}
    rows: list[dict[str, object]] = []
    truncated = False
    for count, entry in enumerate(entries):
        if not entry.name.isdigit():
            continue
        if count > _DIAGNOSTICS_PROC_SCAN_LIMIT:
            truncated = True
            break
        parsed = None
        try:
            parsed = _parse_proc_stat(
                (entry / "stat").read_text(encoding="utf-8", errors="replace")
            )
        except OSError:
            continue
        if parsed is None:
            continue
        pid = int(entry.name)
        if parsed["session"] != session_id and pid != session_id:
            continue
        row: dict[str, object] = {"pid": pid, **parsed}
        try:
            raw = (entry / "cmdline").read_bytes()[:_DIAGNOSTICS_CMDLINE_LIMIT]
        except OSError:
            row["cmdline"] = None
        else:
            row["cmdline"] = raw.replace(b"\0", b" ").decode("utf-8", "replace").strip()
        rows.append(row)
    rows.sort(key=lambda row: int(row["pid"]))
    return {"truncated": truncated, "processes": rows}


def _wait_for_rq_terminal_status(
    job,
    *,
    timeout_seconds: float = _RQ_TERMINAL_STABILIZATION_SECONDS,
    monotonic=time.monotonic,
    sleep=time.sleep,
) -> JobStatus | None:
    """Read RQ until its metadata follows the already-durable SQLite terminal."""
    if (
        isinstance(timeout_seconds, bool)
        or not isinstance(timeout_seconds, (int, float))
        or timeout_seconds <= 0
        or timeout_seconds > _MAX_ACCEPTANCE_TIMEOUT_SECONDS
    ):
        raise ValueError("RQ stabilization timeout is outside the acceptance bound")
    deadline = monotonic() + float(timeout_seconds)
    while True:
        status = job.get_status(refresh=True)
        if status in {JobStatus.FAILED, JobStatus.FINISHED, JobStatus.STOPPED}:
            return status
        remaining = deadline - monotonic()
        if remaining <= 0:
            return status
        sleep(min(_RQ_TERMINAL_POLL_SECONDS, remaining))


def _path_free_evidence_token(value: str | None) -> str | None:
    if value is None:
        return None
    if isinstance(value, str) and _PATH_FREE_EVIDENCE_TOKEN.fullmatch(value):
        return value
    return "REDACTED"


def _terminal_before_activity_message(
    evidence: TerminalLifecycleEvidence,
) -> str:
    error_code = _path_free_evidence_token(evidence.error_code) or "UNAVAILABLE"
    error_reason_code = (
        _path_free_evidence_token(evidence.error_reason_code) or "UNAVAILABLE"
    )
    return (
        f"{_TERMINAL_BEFORE_REQUIRED_ACTIVITY} "
        f"error_code={error_code} error_reason_code={error_reason_code}"
    )


def _private_reference_profile_document(
    fixture: AcceptanceFixture,
) -> dict[str, object]:
    workflow_inputs = fixture.workflow_inputs.to_dict()
    config = workflow_inputs.get("config")
    standard = config.get("standard") if isinstance(config, dict) else None
    reference = standard.get("reference") if isinstance(standard, dict) else None
    if not isinstance(reference, dict):
        raise AssertionError("acceptance fixture reference binding is unavailable")
    transcriptome = fixture.transcriptome
    transcript_fasta = getattr(transcriptome, "transcript_fasta", None)
    if not isinstance(transcript_fasta, Path) or not transcript_fasta.is_absolute():
        raise AssertionError("acceptance fixture transcriptome binding is unavailable")
    private_binding = {
        "schema_version": BULK_RNASEQ_REFERENCE_BINDING_CONTRACT,
        "reference": deepcopy(reference),
        "transcriptome": {
            "reference_id": transcriptome.reference_id,
            "fasta_sha256": transcriptome.fasta_sha256,
            "gtf_sha256": transcriptome.gtf_sha256,
            "transcript_fasta": str(transcript_fasta),
            "transcript_fasta_sha256": transcriptome.transcript_fasta_sha256,
        },
    }
    return {
        "schema_version": PRIVATE_REFERENCE_PROFILE_SCHEMA_VERSION,
        "profiles": {
            _REFERENCE_CONFIG_KEY: {"bindings": {"bulk-rnaseq": private_binding}}
        },
    }


def _public_fixture_inputs(fixture: AcceptanceFixture) -> WorkflowInputs:
    document = fixture.workflow_inputs.to_dict()
    config = deepcopy(document["config"])
    standard = config.get("standard") if isinstance(config, dict) else None
    if not isinstance(standard, dict) or not isinstance(
        standard.pop("reference", None), dict
    ):
        raise AssertionError("acceptance fixture reference binding is unavailable")
    return WorkflowInputs(
        config=config,
        samples=deepcopy(document["samples"]),
        options=deepcopy(document["options"]),
    )


def _matches_regular_file_identity(
    path: Path,
    identity: tuple[int, int],
    *,
    expected_mode: int,
) -> bool:
    try:
        observed = path.stat(follow_symlinks=False)
    except OSError:
        return False
    return (
        stat.S_ISREG(observed.st_mode)
        and stat.S_IMODE(observed.st_mode) == expected_mode
        and observed.st_uid == os.getuid()
        and observed.st_nlink == 1
        and (observed.st_dev, observed.st_ino) == identity
    )


def _matches_directory_identity(
    path: Path,
    identity: tuple[int, int],
) -> bool:
    try:
        observed = path.stat(follow_symlinks=False)
    except OSError:
        return False
    return (
        stat.S_ISDIR(observed.st_mode)
        and stat.S_IMODE(observed.st_mode) == 0o700
        and observed.st_uid == os.getuid()
        and (observed.st_dev, observed.st_ino) == identity
    )


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
