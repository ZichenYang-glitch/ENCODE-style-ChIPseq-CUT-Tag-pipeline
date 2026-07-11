"""Workflow run lifecycle service."""

from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping
from datetime import datetime, timezone
from threading import RLock
from typing import Any
from uuid import uuid4

from encode_pipeline.platform.adapters import WorkflowInputs
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

    def recover_interrupted_runs(self) -> tuple[RunRecord, ...]:
        """Mark work interrupted by an API restart as failed.

        ``CREATED`` and ``PLANNED`` are quiescent persisted states and remain
        available for user action. ``VALIDATING``, ``QUEUED``, and ``RUNNING``
        imply work that belonged to the prior API process, so this method
        records an explicit terminal failure instead of claiming that work is
        still active after that process has disappeared.
        """
        recoverable_statuses = {
            RunStatus.VALIDATING,
            RunStatus.QUEUED,
            RunStatus.RUNNING,
        }
        with self._lock:
            recovered: list[RunRecord] = []
            for current in self._repository.list_runs():
                if current.status not in recoverable_statuses:
                    continue

                now = datetime.now(timezone.utc)
                issue = Issue(
                    code="RUN_INTERRUPTED_BY_API_RESTART",
                    message="Run was interrupted by an API restart.",
                    severity="error",
                    path="run_id",
                    source="run_service",
                    hint="Review the run events and submit a new preflight if needed.",
                )
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
                    self._repository.update_run(
                        updated,
                        expected_status=current.status,
                        event=RunEventDraft(
                            event_type="run_recovered_after_restart",
                            message=(
                                "Run marked failed because the API restarted before "
                                "active work completed."
                            ),
                            status=RunStatus.FAILED,
                            stage=current.current_stage,
                            context={
                                "previous_status": current.status.value,
                                "new_status": RunStatus.FAILED.value,
                                "reason_code": "API_RESTART_INTERRUPTED",
                            },
                            issue=issue,
                        ),
                    )
                except (ConcurrentRunUpdateError, KeyError):
                    # A concurrent writer already changed or deleted the run.
                    # Re-reading it on the next startup keeps recovery idempotent.
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

    def transition_run(
        self,
        run_id: str,
        to_status: RunStatus | str,
        *,
        stage: str | None = None,
        message: str | None = None,
        context: Mapping[str, Any] | None = None,
        issue: Issue | None = None,
    ) -> RunRecord:
        """Transition a run to a new status, enforcing the PR99 graph."""
        with self._lock:
            current = self._repository.get_run(run_id)
            if isinstance(to_status, str):
                to_status = RunStatus(to_status)
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
                    event_type="status_changed",
                    message=event_message,
                    status=to_status,
                    stage=updated_stage,
                    context=event_context,
                    issue=issue,
                ),
            )
            return updated

    def cancel_run(
        self,
        run_id: str,
        reason: str | None = None,
    ) -> RunRecord:
        """Cancel an active run, or return an already-terminal run unchanged."""
        with self._lock:
            current = self._repository.get_run(run_id)
            if current.status.is_terminal:
                return current

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
