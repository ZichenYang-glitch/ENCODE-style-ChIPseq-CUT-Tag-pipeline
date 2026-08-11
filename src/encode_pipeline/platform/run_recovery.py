"""Public-safe diagnostics for explicit administrator run recovery."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any

from encode_pipeline.platform.execution import RunExecutionAssignment
from encode_pipeline.platform.runs import RunStatus


class ExecutionQueueEvidenceState(str, Enum):
    """Bounded projection of execution evidence from a queue backend."""

    UNAVAILABLE = "unavailable"
    MISSING = "missing"
    IDENTITY_DRIFT = "identity_drift"
    QUEUED = "queued"
    DEFERRED = "deferred"
    SCHEDULED = "scheduled"
    STARTED_LIVE = "started_live"
    STARTED_UNPROVEN = "started_unproven"
    FINISHED = "finished"
    FAILED = "failed"
    STOPPED = "stopped"
    CANCELED = "canceled"
    UNKNOWN = "unknown"


_EXACT_TERMINAL_QUEUE_STATES = frozenset(
    {
        ExecutionQueueEvidenceState.FINISHED,
        ExecutionQueueEvidenceState.FAILED,
        ExecutionQueueEvidenceState.STOPPED,
        ExecutionQueueEvidenceState.CANCELED,
    }
)
_LIVE_QUEUE_STATES = frozenset(
    {
        ExecutionQueueEvidenceState.QUEUED,
        ExecutionQueueEvidenceState.DEFERRED,
        ExecutionQueueEvidenceState.SCHEDULED,
        ExecutionQueueEvidenceState.STARTED_LIVE,
    }
)


@dataclass(frozen=True)
class RunExecutionQueueEvidence:
    """One path-free queue observation after strict assignment inspection."""

    state: ExecutionQueueEvidenceState
    requeue_delivery_matches_request: bool = False

    def __post_init__(self) -> None:
        if not isinstance(self.state, ExecutionQueueEvidenceState):
            raise ValueError("queue evidence state is invalid")
        if not isinstance(self.requeue_delivery_matches_request, bool):
            raise ValueError("requeue delivery evidence must be a bool")

    @property
    def is_missing(self) -> bool:
        return self.state is ExecutionQueueEvidenceState.MISSING

    @property
    def is_exact_terminal(self) -> bool:
        return self.state in _EXACT_TERMINAL_QUEUE_STATES

    @property
    def is_live(self) -> bool:
        return self.state in _LIVE_QUEUE_STATES

    @property
    def permits_requeue(self) -> bool:
        return self.is_missing or self.is_exact_terminal

    def to_dict(self) -> dict[str, str | bool]:
        return {
            "state": self.state.value,
            "requeue_delivery_matches_request": (self.requeue_delivery_matches_request),
        }


class RunResultIndexingState(str, Enum):
    """Path-free projection of canonical artifact and QC indexing state."""

    NOT_STARTED = "not_started"
    ARTIFACT_INDEXING = "artifact_indexing"
    ARTIFACT_FAILED = "artifact_failed"
    ARTIFACT_INDEXED = "artifact_indexed"
    QC_INDEXING = "qc_indexing"
    QC_FAILED = "qc_failed"
    QC_INDEXED = "qc_indexed"
    QC_INVALIDATED = "qc_invalidated"
    INCONSISTENT = "inconsistent"


class RunRecoveryCleanupState(str, Enum):
    """Whether managed execution cleanup gates a terminal recovery action."""

    NOT_REQUIRED = "not_required"
    PENDING = "pending"
    BLOCKED = "blocked"


class RunRecoveryGap(str, Enum):
    """Stable doctor categories for one diagnosed run."""

    DATABASE = "database"
    QUEUE = "queue"
    CALLBACK = "callback"
    RESULT_INDEXING = "result_indexing"
    CLEANUP = "cleanup"


class RunRecoveryAction(str, Enum):
    """Explicit administrator mutations currently safe for one diagnosis."""

    FAIL = "fail"
    REQUEUE = "requeue"


@dataclass(frozen=True)
class RunResultIndexingDiagnostic:
    """Safe summary derived only from the persisted ``RunResultState``."""

    state: RunResultIndexingState
    artifact_attempt_status: str | None
    artifact_outcome: str | None
    artifact_generation_present: bool
    qc_attempt_status: str | None
    qc_outcome: str | None
    qc_generation_present: bool

    def __post_init__(self) -> None:
        if not isinstance(self.state, RunResultIndexingState):
            raise ValueError("result indexing state is invalid")
        for name in (
            "artifact_attempt_status",
            "artifact_outcome",
            "qc_attempt_status",
            "qc_outcome",
        ):
            value = getattr(self, name)
            if value is not None and (not isinstance(value, str) or not value):
                raise ValueError(f"{name} must be a non-empty string or None")
        for name in ("artifact_generation_present", "qc_generation_present"):
            if not isinstance(getattr(self, name), bool):
                raise ValueError(f"{name} must be a bool")

    def to_dict(self) -> dict[str, Any]:
        return {
            "state": self.state.value,
            "artifact_attempt_status": self.artifact_attempt_status,
            "artifact_outcome": self.artifact_outcome,
            "artifact_generation_present": self.artifact_generation_present,
            "qc_attempt_status": self.qc_attempt_status,
            "qc_outcome": self.qc_outcome,
            "qc_generation_present": self.qc_generation_present,
        }


@dataclass(frozen=True)
class RunRecoveryDiagnostic:
    """Administrator-facing diagnosis without inputs, paths, or raw errors."""

    run_id: str
    workflow_id: str
    status: RunStatus
    diagnosis_code: str
    assignment: RunExecutionAssignment | None = field(repr=False)
    queue_evidence: RunExecutionQueueEvidence
    build_identity_present: bool
    result_indexing: RunResultIndexingDiagnostic
    cleanup: RunRecoveryCleanupState
    gaps: tuple[RunRecoveryGap, ...] = ()
    allowed_actions: tuple[RunRecoveryAction, ...] = ()
    issue_codes: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if not isinstance(self.run_id, str) or not self.run_id.strip():
            raise ValueError("run_id must be a non-empty string")
        if not isinstance(self.workflow_id, str) or not self.workflow_id.strip():
            raise ValueError("workflow_id must be a non-empty string")
        if not isinstance(self.status, RunStatus):
            raise ValueError("status must be a RunStatus")
        if not isinstance(self.diagnosis_code, str) or not self.diagnosis_code:
            raise ValueError("diagnosis_code must be a non-empty string")
        if self.assignment is not None:
            if not isinstance(self.assignment, RunExecutionAssignment):
                raise ValueError("assignment is invalid")
            if self.assignment.run_id != self.run_id:
                raise ValueError("assignment does not match run_id")
        if not isinstance(self.queue_evidence, RunExecutionQueueEvidence):
            raise ValueError("queue_evidence is invalid")
        if not isinstance(self.build_identity_present, bool):
            raise ValueError("build_identity_present must be a bool")
        if not isinstance(self.result_indexing, RunResultIndexingDiagnostic):
            raise ValueError("result_indexing is invalid")
        if not isinstance(self.cleanup, RunRecoveryCleanupState):
            raise ValueError("cleanup state is invalid")
        normalized_gaps = tuple(self.gaps)
        if any(not isinstance(gap, RunRecoveryGap) for gap in normalized_gaps):
            raise ValueError("gaps contain an invalid value")
        if len(normalized_gaps) != len(set(normalized_gaps)):
            raise ValueError("gaps must not contain duplicates")
        object.__setattr__(self, "gaps", normalized_gaps)
        normalized_actions = tuple(self.allowed_actions)
        if any(
            not isinstance(action, RunRecoveryAction) for action in normalized_actions
        ):
            raise ValueError("allowed_actions contain an invalid value")
        if len(normalized_actions) != len(set(normalized_actions)):
            raise ValueError("allowed_actions must not contain duplicates")
        object.__setattr__(self, "allowed_actions", normalized_actions)
        normalized_codes = tuple(self.issue_codes)
        if any(not isinstance(code, str) or not code for code in normalized_codes):
            raise ValueError("issue_codes must contain non-empty strings")
        object.__setattr__(self, "issue_codes", normalized_codes)

    @property
    def can_fail(self) -> bool:
        return RunRecoveryAction.FAIL in self.allowed_actions

    @property
    def can_requeue(self) -> bool:
        return RunRecoveryAction.REQUEUE in self.allowed_actions

    def to_dict(self) -> dict[str, Any]:
        assignment = self.assignment
        assignment_payload = None
        if assignment is not None:
            assignment_payload = {
                "job_id": assignment.job_id,
                "backend": assignment.backend,
                "queue_name": assignment.queue_name,
                "dispatched": assignment.dispatched_at is not None,
                "claimed": assignment.claimed_at is not None,
                "cancellation_requested": (
                    assignment.cancellation_requested_at is not None
                ),
                "cancellation_acknowledged": (
                    assignment.cancellation_acknowledged_at is not None
                ),
                "requeue_requested": assignment.requeue_requested_at is not None,
                "requeue_confirmed": assignment.requeue_confirmed_at is not None,
            }
        return {
            "run_id": self.run_id,
            "workflow_id": self.workflow_id,
            "status": self.status.value,
            "diagnosis_code": self.diagnosis_code,
            "assignment": assignment_payload,
            "queue_evidence": self.queue_evidence.to_dict(),
            "build_identity_present": self.build_identity_present,
            "result_indexing": self.result_indexing.to_dict(),
            "cleanup": self.cleanup.value,
            "gaps": [gap.value for gap in self.gaps],
            "allowed_actions": [action.value for action in self.allowed_actions],
            "issue_codes": list(self.issue_codes),
        }


@dataclass(frozen=True)
class RunRecoveryActionResult:
    """Bounded result of one committed administrator recovery action."""

    run_id: str
    action: RunRecoveryAction
    previous_status: RunStatus
    status: RunStatus
    assignment: RunExecutionAssignment
    reason_code: str
    changed: bool

    def __post_init__(self) -> None:
        if not isinstance(self.run_id, str) or not self.run_id.strip():
            raise ValueError("run_id must be a non-empty string")
        if not isinstance(self.action, RunRecoveryAction):
            raise ValueError("action must be a RunRecoveryAction")
        if not isinstance(self.previous_status, RunStatus):
            raise ValueError("previous_status must be a RunStatus")
        if not isinstance(self.status, RunStatus):
            raise ValueError("status must be a RunStatus")
        if (
            not isinstance(self.assignment, RunExecutionAssignment)
            or self.assignment.run_id != self.run_id
        ):
            raise ValueError("assignment does not match run_id")
        if not isinstance(self.reason_code, str) or not self.reason_code:
            raise ValueError("reason_code must be a non-empty string")
        if not isinstance(self.changed, bool):
            raise ValueError("changed must be a bool")

    def to_dict(self) -> dict[str, Any]:
        return {
            "run_id": self.run_id,
            "action": self.action.value,
            "previous_status": self.previous_status.value,
            "status": self.status.value,
            "reason_code": self.reason_code,
            "changed": self.changed,
            "job_id": self.assignment.job_id,
            "requeue_requested_at": _optional_datetime(
                self.assignment.requeue_requested_at
            ),
            "requeue_confirmed_at": _optional_datetime(
                self.assignment.requeue_confirmed_at
            ),
        }


@dataclass(frozen=True)
class RunRecoverySummary:
    """Aggregate, path-free recovery gap counts for local doctor output."""

    runs_examined: int
    database_gaps: int
    queue_gaps: int
    callback_gaps: int
    result_indexing_gaps: int
    cleanup_gaps: int
    queue_unavailable: bool

    def __post_init__(self) -> None:
        for name in (
            "runs_examined",
            "database_gaps",
            "queue_gaps",
            "callback_gaps",
            "result_indexing_gaps",
            "cleanup_gaps",
        ):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value < 0:
                raise ValueError(f"{name} must be a nonnegative integer")
            if value > self.runs_examined:
                raise ValueError(f"{name} cannot exceed runs_examined")
        if not isinstance(self.queue_unavailable, bool):
            raise ValueError("queue_unavailable must be a bool")

    def to_dict(self) -> dict[str, int | bool]:
        return {
            "runs_examined": self.runs_examined,
            "database_gaps": self.database_gaps,
            "queue_gaps": self.queue_gaps,
            "callback_gaps": self.callback_gaps,
            "result_indexing_gaps": self.result_indexing_gaps,
            "cleanup_gaps": self.cleanup_gaps,
            "queue_unavailable": self.queue_unavailable,
        }


def _optional_datetime(value: object) -> str | None:
    return value.isoformat() if value is not None else None
