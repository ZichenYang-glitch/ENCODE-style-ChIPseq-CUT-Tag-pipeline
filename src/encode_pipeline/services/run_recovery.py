"""Explicit, fail-closed administrator recovery for durable workflow runs."""

from __future__ import annotations

from collections.abc import Callable
from datetime import datetime, timezone
import re
from typing import NoReturn

from encode_pipeline.platform.builds import WorkflowBuildIdentity
from encode_pipeline.platform.execution import RunExecutionAssignment
from encode_pipeline.platform.result_generations import RunResultState
from encode_pipeline.platform.results import Issue
from encode_pipeline.platform.run_recovery import (
    ExecutionQueueEvidenceState,
    RunExecutionQueueEvidence,
    RunRecoveryAction,
    RunRecoveryActionResult,
    RunRecoveryCleanupState,
    RunRecoveryDiagnostic,
    RunRecoveryGap,
    RunRecoverySummary,
    RunResultIndexingDiagnostic,
    RunResultIndexingState,
)
from encode_pipeline.platform.runs import RunEvent, RunRecord, RunStatus
from encode_pipeline.services.run_queue import (
    RunQueueError,
    RunQueueIdentityError,
    RunQueueJobUnavailableError,
    RunQueueUnavailableError,
    RunRecoveryQueue,
)
from encode_pipeline.services.run_repositories import (
    ConcurrentRunUpdateError,
    RunEventDraft,
    RunRepository,
)


RUN_RECOVERY_FAIL_EVENT = "run_failed_by_admin_recovery"
RUN_RECOVERY_FAIL_REASON_CODE = "RUN_FAILED_BY_ADMIN_RECOVERY"
RUN_RECOVERY_REQUEUE_REQUESTED_EVENT = "run_requeue_requested_by_admin"
RUN_RECOVERY_REQUEUE_CONFIRMED_EVENT = "run_requeue_confirmed_by_admin"
RUN_RECOVERY_REQUEUE_REASON_CODE = "RUN_REQUEUED_BY_ADMIN_RECOVERY"

_ACTIVE_STATUSES = frozenset({RunStatus.QUEUED, RunStatus.RUNNING})
_SCHEDULABLE_QUEUE_STATES = frozenset(
    {
        ExecutionQueueEvidenceState.QUEUED,
        ExecutionQueueEvidenceState.DEFERRED,
        ExecutionQueueEvidenceState.SCHEDULED,
    }
)
_CONFIRMABLE_QUEUE_STATES = _SCHEDULABLE_QUEUE_STATES | {
    ExecutionQueueEvidenceState.STARTED_LIVE,
    ExecutionQueueEvidenceState.FINISHED,
    ExecutionQueueEvidenceState.FAILED,
    ExecutionQueueEvidenceState.STOPPED,
    ExecutionQueueEvidenceState.CANCELED,
}
_IN_FLIGHT_RESULT_STATES = frozenset(
    {
        RunResultIndexingState.NOT_STARTED,
        RunResultIndexingState.ARTIFACT_INDEXING,
        RunResultIndexingState.ARTIFACT_INDEXED,
        RunResultIndexingState.QC_INDEXING,
    }
)
_ERROR_MESSAGES = {
    "RUN_RECOVERY_NOT_FOUND": "Run was not found.",
    "RUN_RECOVERY_CONFLICT": "Run recovery preconditions did not match.",
    "RUN_RECOVERY_NOT_SAFE": "Requested run recovery action is not safe.",
    "RUN_RECOVERY_QUEUE_UNAVAILABLE": "Run recovery queue is unavailable.",
    "RUN_RECOVERY_CLEANUP_FAILED": "Run recovery cleanup failed.",
    "RUN_RECOVERY_DATA_INVALID": "Run recovery data is invalid.",
}


class RunRecoveryError(RuntimeError):
    """Stable, disclosure-safe failure from the public recovery service."""

    def __init__(self, code: str) -> None:
        if code not in _ERROR_MESSAGES:
            raise ValueError("run recovery error code is invalid")
        self.code = code
        super().__init__(_ERROR_MESSAGES[code])


class RunRecoveryService:
    """Diagnose and explicitly close recoverable local execution gaps."""

    def __init__(
        self,
        repository: RunRepository,
        queue: RunRecoveryQueue,
        *,
        cleanup: Callable[[str], bool] | None = None,
        cleanup_endpoint_identity: str | None = None,
        clock: Callable[[], datetime] | None = None,
    ) -> None:
        if repository is None:
            raise ValueError("repository is required")
        if queue is None:
            raise ValueError("queue is required")
        if cleanup is not None and not callable(cleanup):
            raise ValueError("cleanup must be callable or None")
        if (cleanup is None) != (cleanup_endpoint_identity is None):
            raise ValueError(
                "cleanup and cleanup_endpoint_identity must be configured together"
            )
        if cleanup_endpoint_identity is not None and (
            not isinstance(cleanup_endpoint_identity, str)
            or re.fullmatch(r"[0-9a-f]{64}", cleanup_endpoint_identity) is None
        ):
            raise ValueError("cleanup_endpoint_identity is invalid")
        if clock is not None and not callable(clock):
            raise ValueError("clock must be callable or None")
        self._repository = repository
        self._queue = queue
        self._cleanup = cleanup
        self._cleanup_endpoint_identity = cleanup_endpoint_identity
        self._clock = clock or (lambda: datetime.now(timezone.utc))

    def diagnose(self, run_id: str) -> RunRecoveryDiagnostic:
        """Return one bounded diagnosis without mutating SQLite or the queue."""
        normalized_run_id = _run_id(run_id)
        record = self._read_run(normalized_run_id)
        assignment = self._read_assignment(normalized_run_id)
        identity = self._read_build_identity(normalized_run_id)
        result_state = self._read_result_state(normalized_run_id)
        result_indexing = _result_indexing_diagnostic(result_state)
        events = self._read_events(normalized_run_id)

        issue_codes: list[str] = []
        gaps: set[RunRecoveryGap] = set()
        actions: set[RunRecoveryAction] = set()
        cleanup = RunRecoveryCleanupState.NOT_REQUIRED
        queue_evidence = RunExecutionQueueEvidence(
            state=ExecutionQueueEvidenceState.UNKNOWN
        )

        if identity is None and (
            record.status
            in {
                RunStatus.PLANNED,
                RunStatus.QUEUED,
                RunStatus.RUNNING,
                RunStatus.SUCCEEDED,
            }
            or assignment is not None
        ):
            _database_issue(
                issue_codes,
                gaps,
                "RUN_RECOVERY_BUILD_IDENTITY_MISSING",
            )
        elif identity is not None and identity.workflow_id != record.workflow_id:
            _database_issue(
                issue_codes,
                gaps,
                "RUN_RECOVERY_BUILD_IDENTITY_MISMATCH",
            )

        if (
            record.status
            in {
                RunStatus.CREATED,
                RunStatus.VALIDATING,
            }
            and assignment is not None
        ):
            _database_issue(
                issue_codes,
                gaps,
                "RUN_RECOVERY_ASSIGNMENT_UNEXPECTED",
            )
        if record.status in {RunStatus.CREATED, RunStatus.VALIDATING} and (
            identity is not None
        ):
            _database_issue(
                issue_codes,
                gaps,
                "RUN_RECOVERY_BUILD_IDENTITY_UNEXPECTED",
            )

        if record.status in _ACTIVE_STATUSES:
            if record.ended_at is not None:
                _database_issue(
                    issue_codes,
                    gaps,
                    "RUN_RECOVERY_ACTIVE_RUN_ENDED",
                )
            if record.status is RunStatus.QUEUED and record.started_at is not None:
                _database_issue(
                    issue_codes,
                    gaps,
                    "RUN_RECOVERY_QUEUED_RUN_STARTED",
                )
            if record.status is RunStatus.RUNNING and record.started_at is None:
                _database_issue(
                    issue_codes,
                    gaps,
                    "RUN_RECOVERY_RUNNING_START_MISSING",
                )
            if assignment is None:
                _database_issue(
                    issue_codes,
                    gaps,
                    "RUN_RECOVERY_ASSIGNMENT_MISSING",
                )
                if record.status is RunStatus.RUNNING:
                    _database_issue(
                        issue_codes,
                        gaps,
                        "RUN_RECOVERY_CLAIM_MARKER_MISSING",
                    )
            else:
                if assignment.cancellation_acknowledged_at is not None:
                    _database_issue(
                        issue_codes,
                        gaps,
                        "RUN_RECOVERY_ACTIVE_CANCELLATION_ACKNOWLEDGED",
                    )
                queue_evidence = self._inspect_queue(assignment)
                self._diagnose_active_assignment(
                    record=record,
                    assignment=assignment,
                    evidence=queue_evidence,
                    issue_codes=issue_codes,
                    gaps=gaps,
                    actions=actions,
                )
        elif (
            record.status.is_terminal
            and assignment is not None
            and assignment.requeue_requested_at is not None
            and assignment.requeue_confirmed_at is None
            and not _requeue_request_abandoned(events)
        ):
            # The replacement may finish between Redis acceptance and the
            # confirmation transaction.  Inspect only this exceptional
            # terminal handshake; ordinary historical terminal runs never
            # consult Redis.
            queue_evidence = self._inspect_queue(assignment)
            gaps.add(RunRecoveryGap.CALLBACK)
            issue_codes.append("RUN_RECOVERY_REQUEUE_CONFIRMATION_PENDING")
            if assignment.claimed_at is not None or (
                queue_evidence.state in _CONFIRMABLE_QUEUE_STATES
                and (
                    not queue_evidence.is_exact_terminal
                    or queue_evidence.requeue_delivery_matches_request
                )
            ):
                actions.add(RunRecoveryAction.REQUEUE)
            elif queue_evidence.state is ExecutionQueueEvidenceState.UNAVAILABLE:
                gaps.add(RunRecoveryGap.QUEUE)
                issue_codes.append("RUN_RECOVERY_QUEUE_UNAVAILABLE")
            elif queue_evidence.state is ExecutionQueueEvidenceState.IDENTITY_DRIFT:
                gaps.add(RunRecoveryGap.QUEUE)
                issue_codes.append("RUN_RECOVERY_QUEUE_IDENTITY_DRIFT")
            else:
                gaps.add(RunRecoveryGap.QUEUE)
                issue_codes.append("RUN_RECOVERY_OWNER_UNPROVEN")

        if assignment is not None and assignment.claimed_at is not None:
            cleanup_binding_issue = self._cleanup_binding_issue(assignment)
            cleanup_ownership_unresolved = (
                record.status in _ACTIVE_STATUSES
                and queue_evidence.state
                in {
                    ExecutionQueueEvidenceState.MISSING,
                    ExecutionQueueEvidenceState.STARTED_UNPROVEN,
                    ExecutionQueueEvidenceState.UNKNOWN,
                    ExecutionQueueEvidenceState.IDENTITY_DRIFT,
                    ExecutionQueueEvidenceState.UNAVAILABLE,
                }
            )
            cleanup_action_requested = bool(
                actions.intersection(
                    {RunRecoveryAction.FAIL, RunRecoveryAction.REQUEUE}
                )
            )
            if queue_evidence.is_exact_terminal and record.status in _ACTIVE_STATUSES:
                gaps.add(RunRecoveryGap.CLEANUP)
                cleanup = (
                    RunRecoveryCleanupState.PENDING
                    if cleanup_binding_issue is None
                    else RunRecoveryCleanupState.BLOCKED
                )
            elif cleanup_ownership_unresolved:
                gaps.add(RunRecoveryGap.CLEANUP)
                cleanup = RunRecoveryCleanupState.BLOCKED
            if cleanup_binding_issue is not None and (
                cleanup_action_requested
                or cleanup_ownership_unresolved
                or (
                    queue_evidence.is_exact_terminal
                    and record.status in _ACTIVE_STATUSES
                )
            ):
                issue_codes.append(cleanup_binding_issue)
                gaps.add(RunRecoveryGap.CLEANUP)
                cleanup = RunRecoveryCleanupState.BLOCKED
                actions.difference_update(
                    {RunRecoveryAction.FAIL, RunRecoveryAction.REQUEUE}
                )

        if record.status is RunStatus.SUCCEEDED:
            if result_indexing.state is not RunResultIndexingState.QC_INDEXED:
                if (
                    result_indexing.state in _IN_FLIGHT_RESULT_STATES
                    and assignment is not None
                    and assignment.claimed_at is not None
                ):
                    if queue_evidence.state is ExecutionQueueEvidenceState.UNKNOWN:
                        queue_evidence = self._inspect_queue(assignment)
                    if (
                        queue_evidence.state
                        is not ExecutionQueueEvidenceState.STARTED_LIVE
                    ):
                        _record_result_owner_gap(
                            queue_evidence,
                            issue_codes=issue_codes,
                            gaps=gaps,
                        )
                        gaps.add(RunRecoveryGap.RESULT_INDEXING)
                        issue_codes.append("RUN_RECOVERY_RESULT_INDEXING_PENDING")
                else:
                    gaps.add(RunRecoveryGap.RESULT_INDEXING)
                    issue_codes.append("RUN_RECOVERY_RESULT_INDEXING_PENDING")
        elif _has_result_evidence(result_state):
            _database_issue(
                issue_codes,
                gaps,
                "RUN_RECOVERY_RESULT_STATE_INCONSISTENT",
            )

        if _cleanup_failed(events):
            if "RUN_RECOVERY_CLEANUP_NOT_CONFIRMED" not in issue_codes:
                issue_codes.append("RUN_RECOVERY_CLEANUP_NOT_CONFIRMED")
            cleanup = RunRecoveryCleanupState.BLOCKED
            gaps.add(RunRecoveryGap.CLEANUP)

        if RunRecoveryGap.DATABASE in gaps:
            actions.clear()
        if assignment is not None and (
            record.status in _ACTIVE_STATUSES
            and assignment.cancellation_requested_at is not None
            and assignment.cancellation_acknowledged_at is None
        ):
            actions.clear()
            gaps.add(RunRecoveryGap.CALLBACK)
            issue_codes.append("RUN_RECOVERY_CANCELLATION_ACK_PENDING")

        ordered_gaps = tuple(sorted(gaps, key=lambda item: item.value))
        ordered_actions = tuple(sorted(actions, key=lambda item: item.value))
        return RunRecoveryDiagnostic(
            run_id=record.run_id,
            workflow_id=record.workflow_id,
            status=record.status,
            diagnosis_code=_diagnosis_code(
                record=record,
                assignment=assignment,
                evidence=queue_evidence,
                result_indexing=result_indexing,
                gaps=ordered_gaps,
                issue_codes=tuple(issue_codes),
            ),
            assignment=assignment,
            queue_evidence=queue_evidence,
            build_identity_present=identity is not None,
            result_indexing=result_indexing,
            cleanup=cleanup,
            gaps=ordered_gaps,
            allowed_actions=ordered_actions,
            issue_codes=tuple(dict.fromkeys(issue_codes)),
        )

    def summarize(self) -> RunRecoverySummary:
        """Aggregate typed gap counts for machine-readable doctor output."""
        try:
            records = self._repository.list_runs()
        except Exception:
            raise RunRecoveryError("RUN_RECOVERY_DATA_INVALID") from None
        if not isinstance(records, tuple) or any(
            not isinstance(record, RunRecord) for record in records
        ):
            raise RunRecoveryError("RUN_RECOVERY_DATA_INVALID")

        counts = {gap: 0 for gap in RunRecoveryGap}
        queue_unavailable = False
        for record in records:
            try:
                diagnostic = self.diagnose(record.run_id)
            except RunRecoveryError:
                counts[RunRecoveryGap.DATABASE] += 1
                continue
            for gap in diagnostic.gaps:
                counts[gap] += 1
            queue_unavailable = queue_unavailable or (
                diagnostic.queue_evidence.state
                is ExecutionQueueEvidenceState.UNAVAILABLE
            )
        return RunRecoverySummary(
            runs_examined=len(records),
            database_gaps=counts[RunRecoveryGap.DATABASE],
            queue_gaps=counts[RunRecoveryGap.QUEUE],
            callback_gaps=counts[RunRecoveryGap.CALLBACK],
            result_indexing_gaps=counts[RunRecoveryGap.RESULT_INDEXING],
            cleanup_gaps=counts[RunRecoveryGap.CLEANUP],
            queue_unavailable=queue_unavailable,
        )

    def fail_run(
        self,
        run_id: str,
        *,
        expected_status: RunStatus,
        expected_assignment: RunExecutionAssignment,
    ) -> RunRecoveryActionResult:
        """Fail one exact stuck run after final ownership and cleanup checks."""
        normalized_run_id = _run_id(run_id)
        _expected_recovery_identity(
            normalized_run_id,
            expected_status,
            expected_assignment,
        )
        diagnostic = self.diagnose(normalized_run_id)
        self._require_expected_diagnostic(
            diagnostic,
            expected_status=expected_status,
            expected_assignment=expected_assignment,
        )
        if not diagnostic.can_fail:
            self._raise_unsafe_diagnostic(diagnostic)

        if expected_assignment.claimed_at is not None:
            if self._cleanup_binding_issue(expected_assignment) is not None:
                raise RunRecoveryError("RUN_RECOVERY_CLEANUP_FAILED")
            cleanup_scope = expected_assignment.managed_container_scope
            assert cleanup_scope is not None
            assert self._cleanup is not None
            try:
                cleanup_succeeded = self._cleanup(cleanup_scope)
            except Exception:
                raise RunRecoveryError("RUN_RECOVERY_CLEANUP_FAILED") from None
            if cleanup_succeeded is not True:
                raise RunRecoveryError("RUN_RECOVERY_CLEANUP_FAILED")

        final_evidence = self._inspect_queue(expected_assignment)
        if expected_assignment.claimed_at is None:
            if not (final_evidence.is_missing or final_evidence.is_exact_terminal):
                self._raise_unsafe_evidence(final_evidence)
        elif not final_evidence.is_exact_terminal:
            self._raise_unsafe_evidence(final_evidence)

        current = self._read_run(normalized_run_id)
        if current.status is not expected_status:
            raise RunRecoveryError("RUN_RECOVERY_CONFLICT")
        now = self._now()
        issue = Issue(
            code=RUN_RECOVERY_FAIL_REASON_CODE,
            message="Run was failed by an administrator after ownership loss.",
            severity="error",
            path="run_id",
            source="run_recovery",
            hint="Review the path-free recovery diagnosis before creating a new run.",
        )
        failed = RunRecord(
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
            changed = self._repository.fail_run_by_recovery(
                failed,
                expected_status=expected_status,
                expected_assignment=expected_assignment,
                event=RunEventDraft(
                    event_type=RUN_RECOVERY_FAIL_EVENT,
                    message="Run failed by explicit administrator recovery.",
                    status=RunStatus.FAILED,
                    stage=current.current_stage,
                    context={
                        "previous_status": expected_status.value,
                        "new_status": RunStatus.FAILED.value,
                        "reason_code": RUN_RECOVERY_FAIL_REASON_CODE,
                    },
                    issue=issue,
                ),
            )
        except (ConcurrentRunUpdateError, KeyError):
            raise RunRecoveryError("RUN_RECOVERY_CONFLICT") from None
        except (TypeError, ValueError):
            raise RunRecoveryError("RUN_RECOVERY_DATA_INVALID") from None
        return RunRecoveryActionResult(
            run_id=normalized_run_id,
            action=RunRecoveryAction.FAIL,
            previous_status=expected_status,
            status=RunStatus.FAILED,
            assignment=expected_assignment,
            reason_code=RUN_RECOVERY_FAIL_REASON_CODE,
            changed=changed,
        )

    def requeue_run(
        self,
        run_id: str,
        *,
        expected_status: RunStatus,
        expected_assignment: RunExecutionAssignment,
    ) -> RunRecoveryActionResult:
        """Requeue one never-claimed assignment or close its pending handshake."""
        normalized_run_id = _run_id(run_id)
        _expected_recovery_identity(
            normalized_run_id,
            expected_status,
            expected_assignment,
        )
        diagnostic = self.diagnose(normalized_run_id)
        self._require_expected_diagnostic(
            diagnostic,
            expected_status=expected_status,
            expected_assignment=expected_assignment,
        )
        if not diagnostic.can_requeue:
            self._raise_unsafe_diagnostic(diagnostic)
        if expected_assignment.requeue_confirmed_at is not None:
            raise RunRecoveryError("RUN_RECOVERY_CONFLICT")

        prepared = expected_assignment
        preparation_created = False
        if expected_assignment.requeue_requested_at is None:
            if expected_status is not RunStatus.QUEUED:
                raise RunRecoveryError("RUN_RECOVERY_NOT_SAFE")
            requested_at = self._now()
            try:
                preparation = self._repository.prepare_execution_requeue(
                    normalized_run_id,
                    expected_status=expected_status,
                    expected_assignment=expected_assignment,
                    requested_at=requested_at,
                    event=RunEventDraft(
                        event_type=RUN_RECOVERY_REQUEUE_REQUESTED_EVENT,
                        message="Run requeue requested by an administrator.",
                        status=RunStatus.QUEUED,
                        stage="execution",
                        context={"reason_code": RUN_RECOVERY_REQUEUE_REASON_CODE},
                    ),
                )
            except (ConcurrentRunUpdateError, KeyError):
                raise RunRecoveryError("RUN_RECOVERY_CONFLICT") from None
            except (TypeError, ValueError):
                raise RunRecoveryError("RUN_RECOVERY_DATA_INVALID") from None
            prepared = preparation.assignment
            preparation_created = preparation.created

        evidence = self._inspect_queue(prepared)
        backend_job_id: str
        if prepared.claimed_at is not None:
            # A monotonic claim after the request transaction is durable proof
            # that the stable replacement identity was delivered. Confirmation
            # never schedules another job, even if Redis evidence later expires.
            backend_job_id = prepared.job_id
        elif (
            evidence.state in _SCHEDULABLE_QUEUE_STATES
            or evidence.state is ExecutionQueueEvidenceState.STARTED_LIVE
        ) and evidence.requeue_delivery_matches_request:
            backend_job_id = prepared.job_id
        elif evidence.is_exact_terminal and evidence.requeue_delivery_matches_request:
            # The exact terminal record carries the request-bound delivery
            # marker written with the replacement job. Closing its SQLite audit
            # marker cannot schedule another worker.
            backend_job_id = prepared.job_id
        elif prepared.claimed_at is None and (
            evidence.permits_requeue
            or evidence.state in _SCHEDULABLE_QUEUE_STATES
            or evidence.state is ExecutionQueueEvidenceState.STARTED_LIVE
        ):
            # A durable claim is the final execution guard.  If the prior
            # replacement acquired ownership, this assignment is no longer an
            # eligible input here; if it never claimed, retrying the same RQ
            # identity cannot repeat scientific execution.
            try:
                backend_job_id = self._queue.requeue_execution(prepared)
            except RunQueueUnavailableError:
                raise RunRecoveryError("RUN_RECOVERY_QUEUE_UNAVAILABLE") from None
            except (RunQueueIdentityError, RunQueueJobUnavailableError, RunQueueError):
                raise RunRecoveryError("RUN_RECOVERY_NOT_SAFE") from None
            except Exception:
                raise RunRecoveryError("RUN_RECOVERY_QUEUE_UNAVAILABLE") from None
        else:
            self._raise_unsafe_evidence(evidence)
        if backend_job_id != prepared.job_id:
            raise RunRecoveryError("RUN_RECOVERY_NOT_SAFE")

        confirmed_at = self._now()
        try:
            confirmed = self._repository.confirm_execution_requeue(
                normalized_run_id,
                expected_status=expected_status,
                expected_assignment=prepared,
                confirmed_at=confirmed_at,
                event=RunEventDraft(
                    event_type=RUN_RECOVERY_REQUEUE_CONFIRMED_EVENT,
                    message="Run requeue confirmed by the execution queue.",
                    status=None,
                    stage="execution",
                    context={"reason_code": RUN_RECOVERY_REQUEUE_REASON_CODE},
                ),
            )
        except (ConcurrentRunUpdateError, KeyError):
            raise RunRecoveryError("RUN_RECOVERY_CONFLICT") from None
        except (TypeError, ValueError):
            raise RunRecoveryError("RUN_RECOVERY_DATA_INVALID") from None
        return RunRecoveryActionResult(
            run_id=normalized_run_id,
            action=RunRecoveryAction.REQUEUE,
            previous_status=expected_status,
            status=self._read_run(normalized_run_id).status,
            assignment=confirmed,
            reason_code=RUN_RECOVERY_REQUEUE_REASON_CODE,
            changed=preparation_created or prepared.requeue_confirmed_at is None,
        )

    def _diagnose_active_assignment(
        self,
        *,
        record: RunRecord,
        assignment: RunExecutionAssignment,
        evidence: RunExecutionQueueEvidence,
        issue_codes: list[str],
        gaps: set[RunRecoveryGap],
        actions: set[RunRecoveryAction],
    ) -> None:
        if assignment.dispatched_at is None:
            _database_issue(
                issue_codes,
                gaps,
                "RUN_RECOVERY_DISPATCH_MARKER_MISSING",
            )
        if record.status is RunStatus.RUNNING and assignment.claimed_at is None:
            _database_issue(
                issue_codes,
                gaps,
                "RUN_RECOVERY_CLAIM_MARKER_MISSING",
            )

        if evidence.state is ExecutionQueueEvidenceState.UNAVAILABLE:
            gaps.add(RunRecoveryGap.QUEUE)
            issue_codes.append("RUN_RECOVERY_QUEUE_UNAVAILABLE")
        elif evidence.state is ExecutionQueueEvidenceState.IDENTITY_DRIFT:
            gaps.add(RunRecoveryGap.QUEUE)
            issue_codes.append("RUN_RECOVERY_QUEUE_IDENTITY_DRIFT")
        elif evidence.state in {
            ExecutionQueueEvidenceState.STARTED_UNPROVEN,
            ExecutionQueueEvidenceState.UNKNOWN,
        }:
            gaps.add(RunRecoveryGap.QUEUE)
            issue_codes.append("RUN_RECOVERY_OWNER_UNPROVEN")

        cancellation_pending = (
            assignment.cancellation_requested_at is not None
            and assignment.cancellation_acknowledged_at is None
        )
        database_consistent = RunRecoveryGap.DATABASE not in gaps
        if assignment.claimed_at is None:
            if evidence.state is ExecutionQueueEvidenceState.MISSING:
                gaps.add(RunRecoveryGap.QUEUE)
                issue_codes.append("RUN_RECOVERY_JOB_MISSING_BEFORE_CLAIM")
            elif evidence.is_exact_terminal:
                gaps.add(RunRecoveryGap.QUEUE)
                issue_codes.append("RUN_RECOVERY_JOB_TERMINAL_BEFORE_CLAIM")
            if database_consistent and not cancellation_pending:
                if assignment.requeue_requested_at is None:
                    if evidence.permits_requeue:
                        actions.update(
                            {RunRecoveryAction.FAIL, RunRecoveryAction.REQUEUE}
                        )
                elif assignment.requeue_confirmed_at is None:
                    issue_codes.append("RUN_RECOVERY_REQUEUE_CONFIRMATION_PENDING")
                    if evidence.state in _CONFIRMABLE_QUEUE_STATES:
                        actions.add(RunRecoveryAction.REQUEUE)
                    elif evidence.permits_requeue:
                        actions.update(
                            {RunRecoveryAction.FAIL, RunRecoveryAction.REQUEUE}
                        )
                elif evidence.permits_requeue:
                    actions.add(RunRecoveryAction.FAIL)
        else:
            if evidence.is_exact_terminal:
                gaps.add(RunRecoveryGap.CALLBACK)
                issue_codes.append("RUN_RECOVERY_TERMINAL_CALLBACK_PENDING")
                if database_consistent and not cancellation_pending:
                    actions.add(RunRecoveryAction.FAIL)
            elif evidence.state in _SCHEDULABLE_QUEUE_STATES:
                gaps.add(RunRecoveryGap.QUEUE)
                issue_codes.append("RUN_RECOVERY_STATUS_QUEUE_MISMATCH")
            elif evidence.state is ExecutionQueueEvidenceState.MISSING:
                gaps.add(RunRecoveryGap.QUEUE)
                issue_codes.append("RUN_RECOVERY_CLAIMED_OWNER_MISSING")

            if (
                assignment.requeue_requested_at is not None
                and assignment.requeue_confirmed_at is None
                and (
                    assignment.claimed_at is not None
                    or evidence.state in _CONFIRMABLE_QUEUE_STATES
                )
            ):
                actions.add(RunRecoveryAction.REQUEUE)

    def _cleanup_binding_issue(
        self,
        assignment: RunExecutionAssignment,
    ) -> str | None:
        scope = assignment.managed_container_scope
        endpoint_identity = assignment.managed_container_endpoint_identity
        if scope is None or endpoint_identity is None:
            return "RUN_RECOVERY_CLEANUP_BINDING_MISSING"
        if self._cleanup is None or self._cleanup_endpoint_identity is None:
            return "RUN_RECOVERY_CLEANUP_UNAVAILABLE"
        if endpoint_identity != self._cleanup_endpoint_identity:
            return "RUN_RECOVERY_CLEANUP_ENDPOINT_MISMATCH"
        return None

    def _read_run(self, run_id: str) -> RunRecord:
        try:
            record = self._repository.get_run(run_id)
        except KeyError:
            raise RunRecoveryError("RUN_RECOVERY_NOT_FOUND") from None
        except Exception:
            raise RunRecoveryError("RUN_RECOVERY_DATA_INVALID") from None
        if not isinstance(record, RunRecord) or record.run_id != run_id:
            raise RunRecoveryError("RUN_RECOVERY_DATA_INVALID")
        return record

    def _read_assignment(self, run_id: str) -> RunExecutionAssignment | None:
        try:
            assignment = self._repository.get_execution_assignment(run_id)
        except Exception:
            raise RunRecoveryError("RUN_RECOVERY_DATA_INVALID") from None
        if assignment is not None and (
            not isinstance(assignment, RunExecutionAssignment)
            or assignment.run_id != run_id
        ):
            raise RunRecoveryError("RUN_RECOVERY_DATA_INVALID")
        return assignment

    def _read_build_identity(self, run_id: str) -> WorkflowBuildIdentity | None:
        try:
            identity = self._repository.get_workflow_build_identity(run_id)
        except Exception:
            raise RunRecoveryError("RUN_RECOVERY_DATA_INVALID") from None
        if identity is not None and not isinstance(identity, WorkflowBuildIdentity):
            raise RunRecoveryError("RUN_RECOVERY_DATA_INVALID")
        return identity

    def _read_result_state(self, run_id: str) -> RunResultState:
        try:
            state = self._repository.get_result_state(run_id)
        except Exception:
            raise RunRecoveryError("RUN_RECOVERY_DATA_INVALID") from None
        if not isinstance(state, RunResultState) or state.run_id != run_id:
            raise RunRecoveryError("RUN_RECOVERY_DATA_INVALID")
        return state

    def _read_events(self, run_id: str) -> tuple[RunEvent, ...]:
        try:
            events = self._repository.list_events(run_id, limit=10_000)
        except Exception:
            raise RunRecoveryError("RUN_RECOVERY_DATA_INVALID") from None
        if not isinstance(events, tuple) or any(
            not isinstance(event, RunEvent) or event.run_id != run_id
            for event in events
        ):
            raise RunRecoveryError("RUN_RECOVERY_DATA_INVALID")
        return events

    def _inspect_queue(
        self,
        assignment: RunExecutionAssignment,
    ) -> RunExecutionQueueEvidence:
        try:
            evidence = self._queue.inspect_execution(assignment)
        except Exception:
            return RunExecutionQueueEvidence(
                state=ExecutionQueueEvidenceState.UNAVAILABLE
            )
        if not isinstance(evidence, RunExecutionQueueEvidence):
            raise RunRecoveryError("RUN_RECOVERY_DATA_INVALID")
        return evidence

    @staticmethod
    def _require_expected_diagnostic(
        diagnostic: RunRecoveryDiagnostic,
        *,
        expected_status: RunStatus,
        expected_assignment: RunExecutionAssignment,
    ) -> None:
        if (
            diagnostic.status is not expected_status
            or diagnostic.assignment != expected_assignment
        ):
            raise RunRecoveryError("RUN_RECOVERY_CONFLICT")

    @staticmethod
    def _raise_unsafe_diagnostic(diagnostic: RunRecoveryDiagnostic) -> NoReturn:
        if diagnostic.queue_evidence.state is ExecutionQueueEvidenceState.UNAVAILABLE:
            raise RunRecoveryError("RUN_RECOVERY_QUEUE_UNAVAILABLE")
        raise RunRecoveryError("RUN_RECOVERY_NOT_SAFE")

    @staticmethod
    def _raise_unsafe_evidence(evidence: RunExecutionQueueEvidence) -> NoReturn:
        if evidence.state is ExecutionQueueEvidenceState.UNAVAILABLE:
            raise RunRecoveryError("RUN_RECOVERY_QUEUE_UNAVAILABLE")
        raise RunRecoveryError("RUN_RECOVERY_NOT_SAFE")

    def _now(self) -> datetime:
        value = self._clock()
        if not isinstance(value, datetime):
            raise RunRecoveryError("RUN_RECOVERY_DATA_INVALID")
        return value


def _run_id(value: object) -> str:
    if not isinstance(value, str) or not value.strip():
        raise RunRecoveryError("RUN_RECOVERY_DATA_INVALID")
    return value.strip()


def _expected_recovery_identity(
    run_id: str,
    expected_status: object,
    expected_assignment: object,
) -> None:
    if (
        not isinstance(expected_status, RunStatus)
        or not isinstance(expected_assignment, RunExecutionAssignment)
        or expected_assignment.run_id != run_id
    ):
        raise RunRecoveryError("RUN_RECOVERY_DATA_INVALID")


def _database_issue(
    issue_codes: list[str],
    gaps: set[RunRecoveryGap],
    code: str,
) -> None:
    issue_codes.append(code)
    gaps.add(RunRecoveryGap.DATABASE)


def _cleanup_failed(events: tuple[RunEvent, ...]) -> bool:
    failure_sequences = [
        event.sequence
        for event in events
        if event.event_type == "execution_cleanup_failed"
        or (
            event.issue is not None
            and event.issue.code == "MANAGED_CONTAINER_CLEANUP_FAILED"
        )
    ]
    if not failure_sequences:
        return False
    closure_sequences = [
        event.sequence
        for event in events
        if event.event_type
        in {
            RUN_RECOVERY_FAIL_EVENT,
            "cancellation_acknowledged",
            "execution_stopped_unexpectedly",
        }
    ]
    return max(failure_sequences) > max(closure_sequences, default=0)


def _requeue_request_abandoned(events: tuple[RunEvent, ...]) -> bool:
    request_sequences = [
        event.sequence
        for event in events
        if event.event_type == RUN_RECOVERY_REQUEUE_REQUESTED_EVENT
    ]
    if not request_sequences:
        return False
    abandonment_sequences = [
        event.sequence
        for event in events
        if event.event_type == RUN_RECOVERY_FAIL_EVENT
    ]
    return max(abandonment_sequences, default=0) > max(request_sequences)


def _record_result_owner_gap(
    evidence: RunExecutionQueueEvidence,
    *,
    issue_codes: list[str],
    gaps: set[RunRecoveryGap],
) -> None:
    if evidence.state is ExecutionQueueEvidenceState.UNAVAILABLE:
        gaps.add(RunRecoveryGap.QUEUE)
        issue_codes.append("RUN_RECOVERY_QUEUE_UNAVAILABLE")
    elif evidence.state is ExecutionQueueEvidenceState.IDENTITY_DRIFT:
        gaps.add(RunRecoveryGap.QUEUE)
        issue_codes.append("RUN_RECOVERY_QUEUE_IDENTITY_DRIFT")
    elif evidence.is_exact_terminal:
        gaps.add(RunRecoveryGap.CALLBACK)
        issue_codes.append("RUN_RECOVERY_RESULT_CALLBACK_PENDING")
    elif evidence.state is ExecutionQueueEvidenceState.MISSING:
        gaps.add(RunRecoveryGap.QUEUE)
        issue_codes.append("RUN_RECOVERY_RESULT_OWNER_MISSING")
    else:
        gaps.add(RunRecoveryGap.QUEUE)
        issue_codes.append("RUN_RECOVERY_OWNER_UNPROVEN")


def _has_result_evidence(state: RunResultState) -> bool:
    return any(
        (
            state.artifact_revision,
            state.artifact_generation,
            state.artifact_attempt_id,
            state.artifact_outcome,
            state.qc_revision,
            state.qc_generation,
            state.qc_attempt_id,
            state.qc_outcome,
        )
    )


def _result_indexing_diagnostic(
    state: RunResultState,
) -> RunResultIndexingDiagnostic:
    if state.artifact_attempt_status == "pending":
        classification = RunResultIndexingState.ARTIFACT_INDEXING
    elif state.artifact_attempt_status == "failed":
        classification = RunResultIndexingState.ARTIFACT_FAILED
    elif state.artifact_attempt_status not in {None, "succeeded"}:
        classification = RunResultIndexingState.INCONSISTENT
    elif state.artifact_generation is None:
        classification = (
            RunResultIndexingState.NOT_STARTED
            if not _has_result_evidence(state)
            else RunResultIndexingState.INCONSISTENT
        )
    elif state.qc_attempt_status == "pending":
        classification = RunResultIndexingState.QC_INDEXING
    elif state.qc_attempt_status == "failed":
        classification = RunResultIndexingState.QC_FAILED
    elif state.qc_attempt_status not in {None, "succeeded"}:
        classification = RunResultIndexingState.INCONSISTENT
    elif state.qc_outcome == "invalidated":
        classification = RunResultIndexingState.QC_INVALIDATED
    elif state.qc_generation is not None and state.qc_outcome == "succeeded":
        classification = RunResultIndexingState.QC_INDEXED
    elif state.qc_generation is None and state.qc_outcome is None:
        classification = RunResultIndexingState.ARTIFACT_INDEXED
    else:
        classification = RunResultIndexingState.INCONSISTENT
    return RunResultIndexingDiagnostic(
        state=classification,
        artifact_attempt_status=state.artifact_attempt_status,
        artifact_outcome=state.artifact_outcome,
        artifact_generation_present=state.artifact_generation is not None,
        qc_attempt_status=state.qc_attempt_status,
        qc_outcome=state.qc_outcome,
        qc_generation_present=state.qc_generation is not None,
    )


def _diagnosis_code(
    *,
    record: RunRecord,
    assignment: RunExecutionAssignment | None,
    evidence: RunExecutionQueueEvidence,
    result_indexing: RunResultIndexingDiagnostic,
    gaps: tuple[RunRecoveryGap, ...],
    issue_codes: tuple[str, ...],
) -> str:
    if RunRecoveryGap.DATABASE in gaps:
        return "RUN_RECOVERY_DATA_INCONSISTENT"
    if "RUN_RECOVERY_CANCELLATION_ACK_PENDING" in issue_codes:
        return "RUN_RECOVERY_CANCELLATION_PENDING"
    if evidence.state is ExecutionQueueEvidenceState.UNAVAILABLE:
        return "RUN_RECOVERY_QUEUE_UNAVAILABLE"
    if evidence.state is ExecutionQueueEvidenceState.IDENTITY_DRIFT:
        return "RUN_RECOVERY_IDENTITY_DRIFT"
    if (
        "RUN_RECOVERY_CLAIMED_OWNER_MISSING" in issue_codes
        or "RUN_RECOVERY_OWNER_UNPROVEN" in issue_codes
    ):
        return "RUN_RECOVERY_OWNER_UNPROVEN"
    for cleanup_code in (
        "RUN_RECOVERY_CLEANUP_BINDING_MISSING",
        "RUN_RECOVERY_CLEANUP_ENDPOINT_MISMATCH",
        "RUN_RECOVERY_CLEANUP_UNAVAILABLE",
    ):
        if cleanup_code in issue_codes:
            return cleanup_code
    if "RUN_RECOVERY_TERMINAL_CALLBACK_PENDING" in issue_codes:
        return "RUN_RECOVERY_TERMINAL_CALLBACK_PENDING"
    if "RUN_RECOVERY_JOB_MISSING_BEFORE_CLAIM" in issue_codes:
        return "RUN_RECOVERY_JOB_MISSING_BEFORE_CLAIM"
    if "RUN_RECOVERY_JOB_TERMINAL_BEFORE_CLAIM" in issue_codes:
        return "RUN_RECOVERY_JOB_TERMINAL_BEFORE_CLAIM"
    if "RUN_RECOVERY_REQUEUE_CONFIRMATION_PENDING" in issue_codes:
        return "RUN_RECOVERY_REQUEUE_CONFIRMATION_PENDING"
    if "RUN_RECOVERY_STATUS_QUEUE_MISMATCH" in issue_codes:
        return "RUN_RECOVERY_STATUS_QUEUE_MISMATCH"
    if "RUN_RECOVERY_CLEANUP_NOT_CONFIRMED" in issue_codes:
        return "RUN_RECOVERY_CLEANUP_NOT_CONFIRMED"
    if RunRecoveryGap.RESULT_INDEXING in gaps:
        return f"RUN_RECOVERY_{result_indexing.state.value.upper()}"
    if record.status not in _ACTIVE_STATUSES:
        return "RUN_RECOVERY_NOT_REQUIRED"
    if assignment is not None and evidence.is_live:
        return "RUN_RECOVERY_OWNER_CONFIRMED"
    return "RUN_RECOVERY_REVIEW_REQUIRED"


__all__ = [
    "RUN_RECOVERY_FAIL_EVENT",
    "RUN_RECOVERY_FAIL_REASON_CODE",
    "RUN_RECOVERY_REQUEUE_CONFIRMED_EVENT",
    "RUN_RECOVERY_REQUEUE_REASON_CODE",
    "RUN_RECOVERY_REQUEUE_REQUESTED_EVENT",
    "RunRecoveryError",
    "RunRecoveryService",
]
