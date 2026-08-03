"""Module-level RQ job entry points."""

from __future__ import annotations

from rq import get_current_job
from rq.timeouts import JobTimeoutException

from encode_pipeline.persistence.migration_admission import MigrationAdmissionError
from encode_pipeline.persistence.runtime import open_run_persistence
from encode_pipeline.platform.adapters import (
    ReferenceProfileBindingAdapter,
    WorkflowAdapter,
    WorkflowInputs,
)
from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.platform.reference_profiles import BoundWorkflowReference
from encode_pipeline.platform.managed_containers import managed_container_scope
from encode_pipeline.platform.planning import WorkspacePathError, WorkspacePathPolicy
from encode_pipeline.platform.results import Issue, Result
from encode_pipeline.platform.runs import RunStatus
from encode_pipeline.platform.result_generations import new_result_attempt_id
from encode_pipeline.services.run_repositories import ConcurrentRunUpdateError
from encode_pipeline.services.runs import RunService
from encode_pipeline.services.managed_containers import ManagedContainerCleaner
from encode_pipeline.workers.runtime import open_worker_runtime
from encode_pipeline.workers.settings import load_worker_settings
from encode_pipeline.workers.timeouts import WorkerHardTimeout


class WorkerJobIdentityError(RuntimeError):
    """The dequeued RQ job does not own the durable execution assignment."""


class WorkerExecutionError(RuntimeError):
    """A correctly-owned execution job failed after its durable claim."""

    def __init__(
        self,
        message: str,
        *,
        reason_code: str | None = None,
    ) -> None:
        super().__init__(message)
        self.reason_code = reason_code


def run_execution_job(run_id: str) -> None:
    """Claim one durable assignment and execute its reconstructed workspace."""
    if not isinstance(run_id, str) or not run_id.strip():
        raise WorkerJobIdentityError("run_id must be a non-empty string")
    run_id = run_id.strip()
    current_job = None

    try:
        current_job = get_current_job()
        if current_job is None or not current_job.id:
            raise WorkerJobIdentityError(
                "run_execution_job must run inside an RQ job context"
            )
        with open_worker_runtime() as runtime:
            try:
                acquired = _initialize_execution_with_runtime(
                    runtime,
                    current_job,
                    run_id,
                )
            except WorkerHardTimeout:
                _handle_runtime_hard_timeout_safely(runtime, run_id)
                raise
            except WorkerJobIdentityError:
                raise
            except WorkerExecutionError:
                raise
            except Exception:
                _record_initialization_failure_fallback(run_id, current_job)
                raise WorkerExecutionError(
                    "durable local workflow execution could not be initialized"
                ) from None
            if not acquired:
                return
            try:
                _execute_claimed_run(runtime, run_id)
            except WorkerHardTimeout:
                _handle_runtime_hard_timeout_safely(runtime, run_id)
                raise
            except WorkerExecutionError:
                raise
            except Exception:
                _record_unexpected_failure_safely(runtime.run_service, run_id)
                raise WorkerExecutionError(
                    "durable local workflow execution failed unexpectedly"
                ) from None
    except WorkerHardTimeout as timeout:
        try:
            _record_initialization_failure_fallback(
                run_id,
                current_job,
                reason_code="WORKER_JOB_TIMEOUT",
            )
        except WorkerHardTimeout:
            pass
        raise timeout from None
    except MigrationAdmissionError as error:
        raise WorkerExecutionError(
            str(error),
            reason_code=error.reason_code,
        ) from None
    except (WorkerExecutionError, WorkerJobIdentityError):
        raise
    except Exception:
        # Runtime composition and early durable-state reads may fail before a
        # RunService is usable. Keep RQ failure metadata public-safe regardless.
        _record_initialization_failure_fallback(run_id, current_job)
        raise WorkerExecutionError(
            "durable local workflow execution could not be initialized"
        ) from None


def _initialize_execution_with_runtime(runtime, current_job, run_id: str) -> bool:
    try:
        record = runtime.run_service.get_run(run_id)
    except KeyError:
        raise WorkerJobIdentityError(
            "execution job does not reference durable workflow state"
        ) from None
    # A persisted workflow whose adapter cannot be rebuilt is an initialization
    # failure, not an identity mismatch. Let it reach the strict SQLite fallback.
    runtime.registry.get(record.workflow_id)
    _require_identity(
        "RQ job origin",
        getattr(current_job, "origin", None),
        runtime.settings.queue_name,
    )
    _require_durable_assignment_identity(runtime, current_job, run_id)
    if not _require_matching_workflow_build(runtime, record):
        return False
    try:
        runtime.run_service.mark_execution_dispatched(
            run_id,
            job_id=current_job.id,
        )
        runtime.run_service.queue_dispatched_run(
            run_id,
            job_id=current_job.id,
            backend="rq",
            queue_name=runtime.settings.queue_name,
        )
        claim = runtime.run_service.claim_execution_assignment(
            run_id,
            job_id=current_job.id,
            backend="rq",
            queue_name=runtime.settings.queue_name,
            context={"workflow_id": record.workflow_id},
        )
    except (ConcurrentRunUpdateError, KeyError, ValueError):
        try:
            if runtime.run_service.get_run(run_id).status is RunStatus.CANCELLED:
                return False
        except KeyError:
            pass
        raise WorkerJobIdentityError(
            "stale execution job does not own durable workflow state"
        ) from None
    if not claim.acquired:
        return False
    return True


def _require_durable_assignment_identity(runtime, current_job, run_id: str) -> None:
    """Reject stale jobs before build admission can mutate durable state."""
    assignment = runtime.run_service.get_execution_assignment(run_id)
    if assignment is None:
        raise WorkerJobIdentityError(
            "stale execution job does not own durable workflow state"
        )
    for label, actual, expected in (
        ("run ID", assignment.run_id, run_id),
        ("job ID", getattr(current_job, "id", None), assignment.job_id),
        ("backend", assignment.backend, "rq"),
        ("queue", assignment.queue_name, runtime.settings.queue_name),
    ):
        _require_identity(label, actual, expected)


def _require_matching_workflow_build(runtime, record) -> bool:
    """Fail closed before claim when durable and local workflow builds differ."""
    persisted = runtime.run_service.get_workflow_build_identity(record.run_id)
    if persisted is None:
        return _fail_workflow_build_identity(
            runtime.run_service,
            record.run_id,
            code="RUN_WORKFLOW_BUILD_IDENTITY_MISSING",
            message="Run has no durable workflow build identity.",
        )

    current_result = _capture_current_workflow_build(runtime, record)
    if current_result.is_failure:
        return _fail_workflow_build_identity(
            runtime.run_service,
            record.run_id,
            code="RUN_WORKFLOW_BUILD_IDENTITY_UNAVAILABLE",
            message="Workflow build identity could not be verified.",
        )
    if not persisted.matches(current_result.value):
        return _fail_workflow_build_identity(
            runtime.run_service,
            record.run_id,
            code="RUN_WORKFLOW_BUILD_IDENTITY_MISMATCH",
            message="Workflow build differs from the preflighted build.",
        )
    return True


def _capture_current_workflow_build(runtime, record):
    try:
        adapter = runtime.registry.get(record.workflow_id)
    except (KeyError, ValueError):
        return _workflow_build_unavailable()
    if not isinstance(adapter, ReferenceProfileBindingAdapter):
        return runtime.build_identity_provider.capture_executable(record.workflow_id)
    config = record.inputs.get("config", {})
    options = record.inputs.get("options", {})
    if not isinstance(config, dict) or not isinstance(options, dict):
        return _workflow_build_unavailable()
    try:
        resolved = runtime.reference_profile_resolver.resolve_run(
            record.run_id,
            record.workflow_id,
            WorkflowInputs(
                config=config,
                samples=record.inputs.get("samples"),
                options=options,
            ),
            require_enabled=False,
        )
    except Exception:
        return _workflow_build_unavailable()
    if (
        not isinstance(resolved, Result)
        or resolved.is_failure
        or not isinstance(resolved.value, BoundWorkflowReference)
        or not isinstance(resolved.value.adapter, WorkflowAdapter)
    ):
        return _workflow_build_unavailable()
    return runtime.build_identity_provider.capture_resolved_executable(
        resolved.value.adapter
    )


def _workflow_build_unavailable() -> Result:
    return Result.failure(
        [
            Issue(
                code="RUN_WORKFLOW_BUILD_IDENTITY_UNAVAILABLE",
                message="Workflow build identity could not be verified.",
            )
        ]
    )


def _fail_workflow_build_identity(
    run_service,
    run_id: str,
    *,
    code: str,
    message: str,
) -> bool:
    issue = Issue(
        code=code,
        message=message,
        severity="error",
        path="workflow",
        source="worker",
        hint="Create a new run and complete preflight with the current build.",
        context={"reason_code": code},
    )
    current = run_service.get_run(run_id)
    if current.status is RunStatus.CANCELLED:
        return False
    if current.status is RunStatus.QUEUED:
        run_service.transition_run(
            run_id,
            RunStatus.FAILED,
            stage="execution",
            message="Workflow build identity verification failed.",
            issue=issue,
            context={"reason_code": code},
        )
    raise WorkerExecutionError("durable workflow build identity check failed")


def _execute_claimed_run(runtime, run_id: str) -> None:
    execution_result = runtime.local_execution_service.execute(run_id)
    if execution_result.is_success:
        artifact_attempt_id = new_result_attempt_id()
        try:
            artifact_result = runtime.artifact_extraction_service.extract(
                run_id,
                attempt_id=artifact_attempt_id,
            )
        except WorkerHardTimeout:
            try:
                runtime.artifact_extraction_service.record_unexpected_failure(
                    run_id,
                    attempt_id=artifact_attempt_id,
                )
            except WorkerHardTimeout as timeout:
                raise timeout from None
            except Exception:
                pass
            raise
        except Exception:
            try:
                runtime.artifact_extraction_service.record_unexpected_failure(
                    run_id,
                    attempt_id=artifact_attempt_id,
                )
            except WorkerHardTimeout as timeout:
                raise timeout from None
            except Exception:
                pass
            return
        if artifact_result.is_failure:
            return
        try:
            artifact_generation = runtime.run_service.get_result_state(
                run_id
            ).artifact_generation
        except WorkerHardTimeout:
            raise
        except Exception:
            return
        if artifact_generation is None:
            return
        qc_attempt_id = new_result_attempt_id()
        try:
            runtime.qc_summary_indexing_service.index(
                run_id,
                artifact_result.value,
                attempt_id=qc_attempt_id,
                expected_artifact_generation=artifact_generation,
            )
        except WorkerHardTimeout:
            try:
                runtime.qc_summary_indexing_service.record_unexpected_failure(
                    run_id,
                    attempt_id=qc_attempt_id,
                    expected_artifact_generation=artifact_generation,
                )
            except WorkerHardTimeout as timeout:
                raise timeout from None
            except Exception:
                pass
            raise
        except Exception:
            try:
                runtime.qc_summary_indexing_service.record_unexpected_failure(
                    run_id,
                    attempt_id=qc_attempt_id,
                    expected_artifact_generation=artifact_generation,
                )
            except WorkerHardTimeout as timeout:
                raise timeout from None
            except Exception:
                pass
        return
    current = runtime.run_service.get_run(run_id)
    if current.status is RunStatus.CANCELLED:
        return
    if current.status is not RunStatus.FAILED:
        _record_unexpected_failure_safely(runtime.run_service, run_id)
    raise WorkerExecutionError("durable local workflow execution failed") from None


def _record_initialization_failure_fallback(
    run_id: str,
    current_job,
    *,
    reason_code: str = "WORKER_INITIALIZATION_FAILED",
) -> None:
    """Best-effort failure mapping using only canonical SQLite state."""
    try:
        job_id = getattr(current_job, "id", None)
        queue_name = getattr(current_job, "origin", None)
    except Exception:
        return
    if (
        not isinstance(job_id, str)
        or not job_id.strip()
        or not isinstance(queue_name, str)
        or not queue_name.strip()
    ):
        return

    persistence = None
    try:
        persistence = open_run_persistence()
        run_service = RunService(
            WorkflowRegistry(),
            repository=persistence.repository,
        )
        assignment = run_service.get_execution_assignment(run_id)
        if (
            assignment is None
            or assignment.run_id != run_id
            or assignment.job_id != job_id
            or assignment.backend != "rq"
            or assignment.queue_name != queue_name
        ):
            return
        current = run_service.get_run(run_id)
        if current.status not in {
            RunStatus.PLANNED,
            RunStatus.QUEUED,
            RunStatus.RUNNING,
        }:
            return
        if assignment.dispatched_at is None:
            run_service.mark_execution_dispatched(run_id, job_id=job_id)
        run_service.queue_dispatched_run(
            run_id,
            job_id=job_id,
            backend="rq",
            queue_name=queue_name,
        )
        _record_unexpected_failure(
            run_service,
            run_id,
            reason_code=reason_code,
        )
    except WorkerHardTimeout:
        raise
    except Exception:
        return
    finally:
        if persistence is not None:
            try:
                persistence.close()
            except WorkerHardTimeout:
                raise
            except Exception:
                pass


def _require_identity(label: str, actual: object, expected: object) -> None:
    if actual != expected:
        raise WorkerJobIdentityError(
            f"stale execution job: {label} does not match durable ownership"
        )


def _handle_runtime_hard_timeout_safely(runtime, run_id: str) -> None:
    """Preserve the active P0 signal while recording only public-safe evidence."""
    try:
        _cleanup_runtime_managed_containers(runtime, run_id)
    except WorkerHardTimeout:
        pass
    except Exception:
        try:
            _record_cleanup_failure_safely(runtime.run_service, run_id)
        except WorkerHardTimeout:
            pass
    try:
        _record_unexpected_failure_safely(
            runtime.run_service,
            run_id,
            reason_code="WORKER_JOB_TIMEOUT",
        )
    except WorkerHardTimeout:
        pass


def _cleanup_runtime_managed_containers(runtime, run_id: str) -> bool:
    cleaner = runtime.managed_container_cleaner
    if cleaner is None:
        return True
    try:
        workspace = runtime.local_run_driver.derive_workspace_dir(run_id)
        result = cleaner.cleanup(managed_container_scope(workspace))
    except (OSError, ValueError, WorkspacePathError):
        return _record_cleanup_failure_safely(runtime.run_service, run_id)
    if result.is_failure:
        return _record_cleanup_failure_safely(runtime.run_service, run_id)
    return True


def _cleanup_settings_managed_containers(
    settings,
    run_id: str,
    run_service,
) -> bool:
    executable = settings.managed_docker_executable
    if executable is None:
        return True
    try:
        workspace = WorkspacePathPolicy(base_dir=settings.workspace_root).resolve(
            run_id
        )
        cleaner = ManagedContainerCleaner(
            executable=executable,
            unix_socket=settings.managed_docker_socket,
        )
        result = cleaner.cleanup(managed_container_scope(workspace))
    except (OSError, ValueError, WorkspacePathError):
        return _record_cleanup_failure_safely(run_service, run_id)
    if result.is_failure:
        return _record_cleanup_failure_safely(run_service, run_id)
    return True


def _record_cleanup_failure_safely(run_service, run_id: str) -> bool:
    """Keep lifecycle non-terminal and persist only a public-safe diagnostic."""
    try:
        run_service.add_event(
            run_id=run_id,
            event_type="execution_cleanup_failed",
            message="Managed workflow cleanup could not be confirmed.",
            status=None,
            stage="execution",
            context={"reason_code": "MANAGED_CONTAINER_CLEANUP_FAILED"},
            issue=Issue(
                code="MANAGED_CONTAINER_CLEANUP_FAILED",
                message="Managed workflow cleanup could not be confirmed.",
                severity="error",
                path="execution.cleanup",
                source="worker",
            ),
        )
    except WorkerHardTimeout:
        raise
    except Exception:
        pass
    return False


def handle_work_horse_killed(job, _retpid, _ret_val, _rusage) -> None:
    """Persist an unexpected RQ work-horse death without trusting Redis state."""
    try:
        if (
            job is None
            or job.func_name != "encode_pipeline.workers.jobs.run_execution_job"
            or len(job.args) != 1
            or not isinstance(job.args[0], str)
            or job.kwargs != {}
        ):
            return
        run_id = job.args[0].strip()
        if not run_id:
            return
        with open_worker_runtime() as runtime:
            if job.origin != runtime.settings.queue_name:
                return
            assignment = runtime.run_service.get_execution_assignment(run_id)
            if (
                assignment is None
                or assignment.job_id != job.id
                or assignment.backend != "rq"
                or assignment.queue_name != runtime.settings.queue_name
            ):
                return
            if not _cleanup_runtime_managed_containers(runtime, run_id):
                return
            runtime.run_service.mark_execution_dispatched(
                run_id,
                job_id=job.id,
            )
            runtime.run_service.queue_dispatched_run(
                run_id,
                job_id=job.id,
                backend="rq",
                queue_name=runtime.settings.queue_name,
            )
            _record_unexpected_failure_safely(
                runtime.run_service,
                run_id,
                reason_code="WORKER_PROCESS_TERMINATED",
            )
    except WorkerHardTimeout:
        raise
    except Exception:
        # RQ's parent worker must remain healthy even if durable failure mapping
        # cannot open SQLite. A later operational reconciliation can retry it.
        return


def handle_execution_stopped(job, _connection) -> None:
    """Acknowledge a stop only after RQ has reaped the execution process group."""
    persistence = None
    try:
        if (
            job is None
            or job.func_name != "encode_pipeline.workers.jobs.run_execution_job"
            or len(job.args) != 1
            or not isinstance(job.args[0], str)
            or job.kwargs != {}
        ):
            return
        run_id = job.args[0].strip()
        job_id = getattr(job, "id", None)
        queue_name = getattr(job, "origin", None)
        if (
            not run_id
            or not isinstance(job_id, str)
            or not job_id.strip()
            or not isinstance(queue_name, str)
            or not queue_name.strip()
        ):
            return

        settings = load_worker_settings()
        if queue_name != settings.queue_name:
            return
        persistence = open_run_persistence(settings.database_url)
        run_service = RunService(
            WorkflowRegistry(),
            repository=persistence.repository,
        )
        assignment = run_service.get_execution_assignment(run_id)
        if (
            assignment is None
            or assignment.run_id != run_id
            or assignment.job_id != job_id
            or assignment.backend != "rq"
            or assignment.queue_name != queue_name
        ):
            return
        if not _cleanup_settings_managed_containers(
            settings,
            run_id,
            run_service,
        ):
            return
        run_service.acknowledge_execution_stop(
            run_id,
            job_id=job_id,
            backend="rq",
            queue_name=queue_name,
        )
    except (JobTimeoutException, WorkerHardTimeout):
        raise
    except Exception:
        # The RQ parent must remain healthy and proceed to its STOPPED metadata.
        # SQLite remains canonical and deliberately stays diagnosable if unavailable.
        return
    finally:
        if persistence is not None:
            try:
                persistence.close()
            except (JobTimeoutException, WorkerHardTimeout):
                raise
            except Exception:
                pass


def _record_unexpected_failure_safely(
    run_service,
    run_id: str,
    *,
    reason_code: str = "WORKER_UNEXPECTED_ERROR",
) -> None:
    """Never allow durable failure mapping to replace the generic job error."""
    try:
        _record_unexpected_failure(
            run_service,
            run_id,
            reason_code=reason_code,
        )
    except WorkerHardTimeout:
        raise
    except Exception:
        return


def _record_unexpected_failure(
    run_service,
    run_id: str,
    *,
    reason_code: str = "WORKER_UNEXPECTED_ERROR",
) -> None:
    """Best-effort public-safe mapping for an unexpected worker exception."""
    try:
        issue = Issue(
            code="RUN_WORKER_FAILED",
            message="The local execution worker failed unexpectedly.",
            severity="error",
            path="execution",
            source="worker",
            hint="Review persisted run events and logs.",
            context={"reason_code": reason_code},
        )
    except WorkerHardTimeout:
        raise
    except Exception:
        return
    for _attempt in range(2):
        try:
            current = run_service.get_run(run_id)
            if current.status.is_terminal:
                return
            if current.status not in {RunStatus.QUEUED, RunStatus.RUNNING}:
                return
            run_service.transition_run(
                run_id,
                RunStatus.FAILED,
                stage="execution",
                message="The local execution worker failed unexpectedly.",
                issue=issue,
                context={"reason_code": reason_code},
            )
            return
        except (ConcurrentRunUpdateError, KeyError, ValueError):
            continue
        except WorkerHardTimeout:
            raise
        except Exception:
            # Durable failure mapping must never expose a backend exception through
            # RQ's failure metadata. Reconciliation can retry the canonical update.
            return
