"""Workflow run lifecycle service."""

from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping
from datetime import datetime, timezone
from threading import RLock
from typing import Any
from uuid import uuid4

from encode_pipeline.platform.adapters import WorkflowInputs
from encode_pipeline.platform.builds import WorkflowBuildIdentity
from encode_pipeline.platform.execution import (
    RunExecutionAssignment,
    RunExecutionCancellationRequest,
    RunExecutionClaim,
    RunExecutionOwnership,
    RunExecutionStopAcknowledgement,
    build_execution_job_id,
)
from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.platform.results import Issue
from encode_pipeline.platform.runs import (
    RunArtifactRef,
    RunEvent,
    RunLogChunk,
    RunRecord,
    RunStatus,
    require_transition,
)
from encode_pipeline.services.run_repositories import (
    ConcurrentRunUpdateError,
    InMemoryRunRepository,
    RunEventDraft,
    RunRepository,
)


class RunCancellationNotAvailableError(RuntimeError):
    """Cancellation cannot safely mutate the current durable run state."""

    def __init__(self, record: RunRecord) -> None:
        super().__init__("Run cancellation is not available in the current state.")
        self.record = record


class RunService:
    """Owner of workflow run lifecycle mechanics."""

    def __init__(
        self,
        registry: WorkflowRegistry,
        id_factory: Callable[[], str] | None = None,
        repository: RunRepository | None = None,
    ) -> None:
        if not isinstance(registry, WorkflowRegistry):
            raise ValueError("RunService registry must be WorkflowRegistry")
        self._registry = registry
        self._id_factory = (
            id_factory if id_factory is not None else lambda: str(uuid4())
        )
        self._repository = (
            repository if repository is not None else InMemoryRunRepository()
        )
        self._lock = RLock()

    def create_run(
        self,
        workflow_id: str,
        inputs: WorkflowInputs,
        tags: Mapping[str, str] | None = None,
    ) -> RunRecord:
        """Create a new run in the created state."""
        with self._lock:
            adapter = self._registry.get(workflow_id)
            run_id = self._id_factory()
            if self._repository.contains_run(run_id):
                raise ValueError(f"Duplicate run_id: {run_id!r}")

            now = datetime.now(timezone.utc)
            record = RunRecord(
                run_id=run_id,
                workflow_id=adapter.metadata.workflow_id,
                inputs=inputs.to_dict(),
                status=RunStatus.CREATED,
                created_at=now,
                updated_at=now,
                started_at=None,
                ended_at=None,
                current_stage=None,
                cancellation_reason=None,
                error=None,
                tags=tags or {},
            )
            self._repository.create_run(
                record,
                RunEventDraft(
                    event_type="status_changed",
                    message="Run created.",
                    status=RunStatus.CREATED,
                    context={
                        "previous_status": None,
                        "new_status": RunStatus.CREATED.value,
                    },
                ),
            )
            return record

    def add_event(
        self,
        run_id: str,
        event_type: str,
        message: str,
        *,
        status: RunStatus | str | None = None,
        stage: str | None = None,
        context: Mapping[str, Any] | None = None,
        issue: Issue | None = None,
    ) -> RunEvent:
        """Record an event for a run without changing run status."""
        with self._lock:
            if isinstance(status, str):
                status = RunStatus(status)
            return self._emit_event(
                run_id=run_id,
                event_type=event_type,
                message=message,
                status=status,
                stage=stage,
                context=context,
                issue=issue,
            )

    def list_events(
        self,
        run_id: str,
        *,
        after: str | None = None,
        limit: int = 50,
    ) -> tuple[RunEvent, ...]:
        """Return ordered events with optional cursor pagination."""
        with self._lock:
            if limit <= 0:
                raise ValueError("limit must be positive")
            return self._repository.list_events(run_id, after=after, limit=limit)

    def append_log(
        self,
        run_id: str,
        stream_name: str,
        lines: Iterable[str],
    ) -> RunLogChunk:
        """Append a log chunk to a run stream."""
        with self._lock:
            return self._repository.append_log(run_id, stream_name, lines)

    def list_logs(
        self,
        run_id: str,
        stream_name: str = "stdout",
        *,
        after: str | None = None,
        limit: int = 50,
    ) -> tuple[RunLogChunk, ...]:
        """Return log chunks for a stream with optional cursor pagination."""
        with self._lock:
            if limit <= 0:
                raise ValueError("limit must be positive")
            return self._repository.list_logs(
                run_id,
                stream_name,
                after=after,
                limit=limit,
            )

    def get_run(self, run_id: str) -> RunRecord:
        """Return the run record for a run ID."""
        with self._lock:
            return self._repository.get_run(run_id)

    def list_runs(self) -> tuple[RunRecord, ...]:
        """Return all runs in creation order."""
        with self._lock:
            return self._repository.list_runs()

    def ensure_execution_assignment(
        self,
        run_id: str,
        *,
        queue_name: str,
        backend: str = "rq",
    ) -> RunExecutionAssignment:
        """Return the durable worker assignment for a planned run.

        The job identity is derived only from ``run_id`` so retries made while
        the run remains planned resolve to the same backend job.  The
        repository returns the first persisted assignment as the canonical
        value when this method is called more than once.
        """
        with self._lock:
            current = self._repository.get_run(run_id)
            if current.status is not RunStatus.PLANNED:
                raise ValueError(
                    "Execution assignments may only be created for planned runs; "
                    f"run {run_id!r} is {current.status.value!r}."
                )
            assignment = RunExecutionAssignment(
                run_id=run_id,
                job_id=build_execution_job_id(run_id),
                backend=backend,
                queue_name=queue_name,
                created_at=datetime.now(timezone.utc),
            )
            persisted = self._repository.ensure_execution_assignment(
                assignment,
                expected_status=RunStatus.PLANNED,
            )
            if persisted.backend != backend or persisted.queue_name != queue_name:
                raise ValueError(
                    "Existing execution assignment does not match the configured "
                    "backend and queue."
                )
            return persisted

    def get_execution_assignment(
        self,
        run_id: str,
    ) -> RunExecutionAssignment | None:
        """Return the durable worker assignment for ``run_id``, if present."""
        with self._lock:
            return self._repository.get_execution_assignment(run_id)

    def mark_execution_dispatched(
        self,
        run_id: str,
        *,
        job_id: str,
    ) -> RunExecutionAssignment:
        """Persist that the configured backend accepted this execution job."""
        with self._lock:
            return self._repository.mark_execution_dispatched(
                run_id,
                job_id=job_id,
                dispatched_at=datetime.now(timezone.utc),
                allowed_statuses=frozenset(
                    {
                        RunStatus.PLANNED,
                        RunStatus.QUEUED,
                        RunStatus.CANCELLED,
                    }
                ),
            )

    def queue_dispatched_run(
        self,
        run_id: str,
        *,
        job_id: str,
        backend: str,
        queue_name: str,
    ) -> RunRecord:
        """Atomically move a dispatched planned run into the worker queue."""
        with self._lock:
            current = self._repository.get_run(run_id)
            assignment = self._repository.get_execution_assignment(run_id)
            if assignment is None:
                raise ValueError("Run has no durable execution assignment.")
            if (
                assignment.job_id != job_id
                or assignment.backend != backend
                or assignment.queue_name != queue_name
            ):
                raise ValueError("execution assignment identity does not match")
            if assignment.dispatched_at is None:
                raise ValueError("execution assignment has not been dispatched")
            if current.status is not RunStatus.PLANNED:
                if current.status in {
                    RunStatus.QUEUED,
                    RunStatus.RUNNING,
                    RunStatus.SUCCEEDED,
                    RunStatus.FAILED,
                    RunStatus.CANCELLED,
                }:
                    return current
                raise ValueError(
                    "Only planned runs may enter the execution queue; "
                    f"run {run_id!r} is {current.status.value!r}."
                )

            now = datetime.now(timezone.utc)
            updated = RunRecord(
                run_id=current.run_id,
                workflow_id=current.workflow_id,
                inputs=current.inputs,
                status=RunStatus.QUEUED,
                created_at=current.created_at,
                updated_at=now,
                started_at=current.started_at,
                ended_at=current.ended_at,
                current_stage="execution",
                cancellation_reason=current.cancellation_reason,
                error=None,
                tags=current.tags,
            )
            queued = self._repository.queue_dispatched_run(
                updated,
                expected_status=RunStatus.PLANNED,
                job_id=job_id,
                backend=backend,
                queue_name=queue_name,
                event=RunEventDraft(
                    event_type="status_changed",
                    message="Run submitted for worker execution.",
                    status=RunStatus.QUEUED,
                    stage="execution",
                    context={
                        "previous_status": RunStatus.PLANNED.value,
                        "new_status": RunStatus.QUEUED.value,
                        "backend": assignment.backend,
                        "job_id": assignment.job_id,
                        "queue_name": assignment.queue_name,
                    },
                ),
            )
            if queued:
                return updated

            raced = self._repository.get_run(run_id)
            if raced.status in {
                RunStatus.QUEUED,
                RunStatus.RUNNING,
                RunStatus.SUCCEEDED,
                RunStatus.FAILED,
                RunStatus.CANCELLED,
            }:
                return raced
            raise ConcurrentRunUpdateError(
                f"Run {run_id!r} changed while it was being queued."
            )

    def claim_execution_assignment(
        self,
        run_id: str,
        *,
        job_id: str,
        backend: str,
        queue_name: str,
        context: Mapping[str, Any] | None = None,
    ) -> RunExecutionClaim:
        """Atomically claim a dispatched job and write its first worker event."""
        with self._lock:
            event_context = dict(context or {})
            event_context.update(
                {
                    "backend": backend,
                    "job_id": job_id,
                    "queue_name": queue_name,
                }
            )
            return self._repository.claim_execution_assignment(
                run_id,
                job_id=job_id,
                backend=backend,
                queue_name=queue_name,
                claimed_at=datetime.now(timezone.utc),
                allowed_statuses=frozenset({RunStatus.QUEUED}),
                event=RunEventDraft(
                    event_type="worker_dependencies_rebuilt",
                    message=(
                        "Worker rebuilt execution dependencies from durable state."
                    ),
                    stage="execution",
                    context=event_context,
                ),
            )

    def recover_interrupted_runs(self) -> tuple[RunRecord, ...]:
        """Recover active runs according to their durable execution owner.

        Validation is API-owned and therefore cannot survive an API restart.
        Queued work requires a durable dispatch marker and running work requires
        an atomic worker claim; those runs remain untouched while active records
        without the required ownership evidence are failed as orphans.
        Quiescent and terminal states are always preserved.
        """
        with self._lock:
            recovered: list[RunRecord] = []
            for current in self._repository.list_runs():
                if current.status is RunStatus.VALIDATING:
                    required_ownership = None
                    issue = Issue(
                        code="RUN_INTERRUPTED_BY_API_RESTART",
                        message="Run was interrupted by an API restart.",
                        severity="error",
                        path="run_id",
                        source="run_service",
                        hint=(
                            "Review the run events and submit a new preflight if "
                            "needed."
                        ),
                    )
                    reason_code = "API_RESTART_INTERRUPTED"
                    event_message = (
                        "Run marked failed because API-owned validation was "
                        "interrupted by an API restart."
                    )
                elif current.status in {RunStatus.QUEUED, RunStatus.RUNNING}:
                    required_ownership = (
                        RunExecutionOwnership.DISPATCHED
                        if current.status is RunStatus.QUEUED
                        else RunExecutionOwnership.CLAIMED
                    )
                    reason_code = "WORKER_OWNERSHIP_NOT_CONFIRMED"
                    issue = Issue(
                        code="RUN_ORPHANED_AFTER_API_RESTART",
                        message=(
                            "Run has an active status but no confirmed durable worker "
                            "ownership."
                        ),
                        severity="error",
                        path="run_id",
                        source="run_service",
                        hint="Review the run events and submit the run again.",
                    )
                    event_message = (
                        "Run marked failed because durable worker ownership was not "
                        "confirmed after the API restart."
                    )
                else:
                    continue

                now = datetime.now(timezone.utc)
                updated = RunRecord(
                    run_id=current.run_id,
                    workflow_id=current.workflow_id,
                    inputs=current.inputs,
                    status=RunStatus.FAILED,
                    created_at=current.created_at,
                    updated_at=now,
                    started_at=current.started_at,
                    ended_at=now,
                    current_stage=current.current_stage,
                    cancellation_reason=current.cancellation_reason,
                    error=issue,
                    tags=current.tags,
                )
                try:
                    was_recovered = self._repository.fail_interrupted_run_if_unowned(
                        updated,
                        expected_status=current.status,
                        required_ownership=required_ownership,
                        event=RunEventDraft(
                            event_type="run_recovered_after_restart",
                            message=event_message,
                            status=RunStatus.FAILED,
                            stage=current.current_stage,
                            context={
                                "previous_status": current.status.value,
                                "new_status": RunStatus.FAILED.value,
                                "reason_code": reason_code,
                            },
                            issue=issue,
                        ),
                    )
                except (ConcurrentRunUpdateError, KeyError):
                    # A concurrent writer already changed or deleted the run.
                    # Re-reading it on the next startup keeps recovery idempotent.
                    continue
                if not was_recovered:
                    continue
                recovered.append(updated)
            return tuple(recovered)

    def _emit_event(
        self,
        run_id: str,
        event_type: str,
        message: str,
        *,
        status: RunStatus | None = None,
        stage: str | None = None,
        context: Mapping[str, Any] | None = None,
        issue: Issue | None = None,
    ) -> RunEvent:
        """Append an event under the existing service lock."""
        return self._repository.add_event(
            run_id,
            RunEventDraft(
                event_type=event_type,
                status=status,
                stage=stage,
                message=message,
                context=context or {},
                issue=issue,
            ),
        )

    def get_workflow_build_identity(
        self,
        run_id: str,
    ) -> WorkflowBuildIdentity | None:
        """Return the durable workflow build bound at successful preflight."""
        with self._lock:
            return self._repository.get_workflow_build_identity(run_id)

    def complete_preflight(
        self,
        run_id: str,
        identity: WorkflowBuildIdentity,
        *,
        stage: str = "preflight",
        message: str = "Local preflight completed; dry-run succeeded.",
        context: Mapping[str, Any] | None = None,
    ) -> RunRecord:
        """Atomically bind build identity and transition VALIDATING to PLANNED."""
        if not isinstance(identity, WorkflowBuildIdentity):
            raise ValueError("identity must be a WorkflowBuildIdentity")
        with self._lock:
            current = self._repository.get_run(run_id)
            require_transition(current.status, RunStatus.PLANNED)
            if identity.workflow_id != current.workflow_id:
                raise ValueError("workflow build identity does not match the run")

            now = datetime.now(timezone.utc)
            updated = RunRecord(
                run_id=current.run_id,
                workflow_id=current.workflow_id,
                inputs=current.inputs,
                status=RunStatus.PLANNED,
                created_at=current.created_at,
                updated_at=now,
                started_at=current.started_at,
                ended_at=current.ended_at,
                current_stage=stage,
                cancellation_reason=current.cancellation_reason,
                error=None,
                tags=current.tags,
            )
            event_context = dict(context or {})
            event_context.update(
                {
                    "previous_status": current.status.value,
                    "new_status": RunStatus.PLANNED.value,
                    "workflow_build_scheme": identity.scheme,
                    "workflow_build_digest": identity.digest,
                }
            )
            self._repository.complete_preflight(
                updated,
                identity,
                expected_status=current.status,
                event=RunEventDraft(
                    event_type="preflight_completed",
                    message=message,
                    status=RunStatus.PLANNED,
                    stage=stage,
                    context=event_context,
                ),
            )
            return updated

    def transition_run(
        self,
        run_id: str,
        to_status: RunStatus | str,
        *,
        stage: str | None = None,
        message: str | None = None,
        event_type: str = "status_changed",
        context: Mapping[str, Any] | None = None,
        issue: Issue | None = None,
    ) -> RunRecord:
        """Transition a run to a new status, enforcing the PR99 graph."""
        with self._lock:
            current = self._repository.get_run(run_id)
            if isinstance(to_status, str):
                to_status = RunStatus(to_status)
            if current.status is RunStatus.RUNNING and to_status is RunStatus.CANCELLED:
                raise ValueError(
                    "RUNNING cancellation requires execution-stop acknowledgement."
                )
            require_transition(current.status, to_status)

            now = datetime.now(timezone.utc)
            updated_stage = stage if stage is not None else current.current_stage
            started_at = current.started_at
            ended_at = current.ended_at
            error = None

            if to_status is RunStatus.RUNNING:
                started_at = now
            if to_status.is_terminal:
                ended_at = now
            if to_status is RunStatus.FAILED and issue is not None:
                error = issue

            updated = RunRecord(
                run_id=current.run_id,
                workflow_id=current.workflow_id,
                inputs=current.inputs,
                status=to_status,
                created_at=current.created_at,
                updated_at=now,
                started_at=started_at,
                ended_at=ended_at,
                current_stage=updated_stage,
                cancellation_reason=current.cancellation_reason,
                error=error,
                tags=current.tags,
            )
            event_message = (
                message
                if message is not None
                else f"Status changed to {to_status.value}."
            )
            event_context = dict(context or {})
            event_context["previous_status"] = current.status.value
            event_context["new_status"] = to_status.value
            self._repository.update_run(
                updated,
                expected_status=current.status,
                event=RunEventDraft(
                    event_type=event_type,
                    message=event_message,
                    status=to_status,
                    stage=updated_stage,
                    context=event_context,
                    issue=issue,
                ),
            )
            return updated

    def request_execution_cancellation(
        self,
        run_id: str,
        *,
        job_id: str,
        backend: str,
        queue_name: str,
        reason: str,
    ) -> RunExecutionCancellationRequest:
        """Atomically persist user intent without claiming process termination."""
        if not isinstance(reason, str) or not reason.strip():
            raise ValueError("cancellation reason must be a non-empty string")
        with self._lock:
            return self._repository.request_execution_cancellation(
                run_id,
                job_id=job_id,
                backend=backend,
                queue_name=queue_name,
                requested_at=datetime.now(timezone.utc),
                reason=reason.strip(),
                event=RunEventDraft(
                    event_type="cancellation_requested",
                    message="Run cancellation requested.",
                    status=RunStatus.RUNNING,
                    stage="execution",
                    context={
                        "job_id": job_id,
                        "backend": backend,
                        "queue_name": queue_name,
                    },
                ),
            )

    def acknowledge_execution_stop(
        self,
        run_id: str,
        *,
        job_id: str,
        backend: str,
        queue_name: str,
    ) -> RunExecutionStopAcknowledgement:
        """Commit the terminal state after RQ has reaped the work horse."""
        issue = Issue(
            code="RUN_WORKER_STOPPED_UNEXPECTEDLY",
            message="The local execution worker stopped unexpectedly.",
            severity="error",
            path="execution",
            source="worker",
            hint="Review persisted run events and logs.",
            context={"reason_code": "WORKER_STOP_WITHOUT_CANCELLATION"},
        )
        with self._lock:
            return self._repository.acknowledge_execution_stop(
                run_id,
                job_id=job_id,
                backend=backend,
                queue_name=queue_name,
                acknowledged_at=datetime.now(timezone.utc),
                cancellation_event=RunEventDraft(
                    event_type="cancellation_acknowledged",
                    message="Run cancellation acknowledged after process termination.",
                    status=RunStatus.CANCELLED,
                    stage="execution",
                    context={
                        "job_id": job_id,
                        "backend": backend,
                        "queue_name": queue_name,
                    },
                ),
                unexpected_stop_event=RunEventDraft(
                    event_type="execution_stopped_unexpectedly",
                    message="The local execution worker stopped unexpectedly.",
                    status=RunStatus.FAILED,
                    stage="execution",
                    context={
                        "previous_status": RunStatus.RUNNING.value,
                        "new_status": RunStatus.FAILED.value,
                        "job_id": job_id,
                        "backend": backend,
                        "queue_name": queue_name,
                        "reason_code": "WORKER_STOP_WITHOUT_CANCELLATION",
                    },
                    issue=issue,
                ),
            )

    def cancel_run(
        self,
        run_id: str,
        reason: str | None = None,
    ) -> RunRecord:
        """Cancel a pre-running run, or return an already-terminal run unchanged.

        ``RUNNING`` cancellation must use :class:`RunCancellationService`, which
        records intent and waits for the worker stop acknowledgement. Expected-
        status writes and a canonical re-read make this pre-running path race
        safely with worker startup: whichever SQLite transition wins determines
        the outcome.
        """
        with self._lock:
            while True:
                current = self._repository.get_run(run_id)
                if current.status.is_terminal:
                    return current
                if current.status is RunStatus.RUNNING:
                    raise RunCancellationNotAvailableError(current)

                now = datetime.now(timezone.utc)
                updated = RunRecord(
                    run_id=current.run_id,
                    workflow_id=current.workflow_id,
                    inputs=current.inputs,
                    status=RunStatus.CANCELLED,
                    created_at=current.created_at,
                    updated_at=now,
                    started_at=current.started_at,
                    ended_at=now,
                    current_stage=current.current_stage,
                    cancellation_reason=reason,
                    error=None,
                    tags=current.tags,
                )
                try:
                    self._repository.update_run(
                        updated,
                        expected_status=current.status,
                        event=RunEventDraft(
                            event_type="status_changed",
                            message="Run cancelled.",
                            status=RunStatus.CANCELLED,
                            stage=current.current_stage,
                            context={
                                "previous_status": current.status.value,
                                "new_status": RunStatus.CANCELLED.value,
                                "cancellation_reason": reason,
                            },
                            issue=None,
                        ),
                    )
                except ConcurrentRunUpdateError:
                    # Another API/worker process advanced the monotonic state.
                    # Re-read SQLite and re-apply the terminal/RUNNING policy.
                    continue
                return updated

    def record_artifact(
        self,
        run_id: str,
        artifact: RunArtifactRef,
    ) -> RunArtifactRef:
        """Record an artifact reference for a run."""
        with self._lock:
            if artifact.run_id != run_id:
                raise ValueError(
                    f"Artifact run_id {artifact.run_id!r} does not match {run_id!r}"
                )
            return self._repository.record_artifact(run_id, artifact)

    def list_artifacts(self, run_id: str) -> tuple[RunArtifactRef, ...]:
        """Return artifact references in insertion order."""
        with self._lock:
            return self._repository.list_artifacts(run_id)
