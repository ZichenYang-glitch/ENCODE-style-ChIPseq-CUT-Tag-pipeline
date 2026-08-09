"""Tests for public-safe administrator recovery diagnostic models."""

from __future__ import annotations

from datetime import datetime, timedelta, timezone
import json

import pytest

from encode_pipeline.platform.execution import RunExecutionAssignment
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
from encode_pipeline.platform.runs import RunStatus


NOW = datetime(2026, 8, 9, 12, 0, tzinfo=timezone.utc)


def _assignment(*, run_id: str = "run-1") -> RunExecutionAssignment:
    return RunExecutionAssignment(
        run_id=run_id,
        job_id="job-1",
        backend="rq",
        queue_name="runs",
        created_at=NOW,
        dispatched_at=NOW + timedelta(seconds=1),
        claimed_at=NOW + timedelta(seconds=2),
        cancellation_requested_at=NOW + timedelta(seconds=3),
        cancellation_reason="operator request /private/lab/run-1",
        cancellation_acknowledged_at=NOW + timedelta(seconds=4),
        requeue_requested_at=NOW + timedelta(seconds=5),
        requeue_confirmed_at=NOW + timedelta(seconds=6),
    )


def _result_indexing_kwargs() -> dict[str, object]:
    return {
        "state": RunResultIndexingState.QC_INDEXED,
        "artifact_attempt_status": "succeeded",
        "artifact_outcome": "indexed",
        "artifact_generation_present": True,
        "qc_attempt_status": "succeeded",
        "qc_outcome": "indexed",
        "qc_generation_present": True,
    }


def _result_indexing() -> RunResultIndexingDiagnostic:
    return RunResultIndexingDiagnostic(**_result_indexing_kwargs())


def _diagnostic_kwargs() -> dict[str, object]:
    return {
        "run_id": "run-1",
        "workflow_id": "workflow-1",
        "status": RunStatus.RUNNING,
        "diagnosis_code": "CANCELLATION_ACK_PENDING",
        "assignment": _assignment(),
        "queue_evidence": RunExecutionQueueEvidence(
            state=ExecutionQueueEvidenceState.STOPPED,
            requeue_delivery_matches_request=True,
        ),
        "build_identity_present": True,
        "result_indexing": _result_indexing(),
        "cleanup": RunRecoveryCleanupState.PENDING,
        "gaps": (RunRecoveryGap.CALLBACK, RunRecoveryGap.CLEANUP),
        "allowed_actions": (RunRecoveryAction.FAIL, RunRecoveryAction.REQUEUE),
        "issue_codes": ("CANCELLATION_ACK_PENDING",),
    }


def _action_result_kwargs() -> dict[str, object]:
    return {
        "run_id": "run-1",
        "action": RunRecoveryAction.REQUEUE,
        "previous_status": RunStatus.QUEUED,
        "status": RunStatus.QUEUED,
        "assignment": _assignment(),
        "reason_code": "RUN_REQUEUED_BY_ADMIN_RECOVERY",
        "changed": True,
    }


def _summary_kwargs() -> dict[str, object]:
    return {
        "runs_examined": 3,
        "database_gaps": 1,
        "queue_gaps": 2,
        "callback_gaps": 1,
        "result_indexing_gaps": 1,
        "cleanup_gaps": 1,
        "queue_unavailable": False,
    }


@pytest.mark.parametrize(
    "state",
    [
        ExecutionQueueEvidenceState.FINISHED,
        ExecutionQueueEvidenceState.FAILED,
        ExecutionQueueEvidenceState.STOPPED,
        ExecutionQueueEvidenceState.CANCELED,
    ],
)
def test_exact_terminal_queue_evidence_is_bounded_and_actionable(state):
    evidence = RunExecutionQueueEvidence(state=state)

    assert evidence.is_exact_terminal is True
    assert evidence.is_live is False
    assert evidence.permits_requeue is True


@pytest.mark.parametrize(
    "state",
    [
        ExecutionQueueEvidenceState.UNAVAILABLE,
        ExecutionQueueEvidenceState.IDENTITY_DRIFT,
        ExecutionQueueEvidenceState.STARTED_UNPROVEN,
        ExecutionQueueEvidenceState.UNKNOWN,
    ],
)
def test_ambiguous_queue_evidence_never_authorizes_mutation(state):
    evidence = RunExecutionQueueEvidence(state=state)

    assert evidence.is_exact_terminal is False
    assert evidence.permits_requeue is False


def test_missing_queue_evidence_permits_only_unclaimed_requeue_decision():
    evidence = RunExecutionQueueEvidence(state=ExecutionQueueEvidenceState.MISSING)

    assert evidence.is_missing is True
    assert evidence.is_exact_terminal is False
    assert evidence.permits_requeue is True


@pytest.mark.parametrize(
    "state",
    [
        ExecutionQueueEvidenceState.QUEUED,
        ExecutionQueueEvidenceState.DEFERRED,
        ExecutionQueueEvidenceState.SCHEDULED,
        ExecutionQueueEvidenceState.STARTED_LIVE,
    ],
)
def test_live_queue_evidence_never_authorizes_requeue(state):
    evidence = RunExecutionQueueEvidence(state=state)

    assert evidence.is_live is True
    assert evidence.permits_requeue is False


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"state": "queued"}, "state is invalid"),
        (
            {
                "state": ExecutionQueueEvidenceState.QUEUED,
                "requeue_delivery_matches_request": 1,
            },
            "must be a bool",
        ),
    ],
)
def test_queue_evidence_rejects_unbounded_repository_data(kwargs, message):
    with pytest.raises(ValueError, match=message):
        RunExecutionQueueEvidence(**kwargs)


def test_queue_evidence_serializes_only_bounded_fields():
    evidence = RunExecutionQueueEvidence(
        state=ExecutionQueueEvidenceState.FINISHED,
        requeue_delivery_matches_request=True,
    )

    assert evidence.to_dict() == {
        "state": "finished",
        "requeue_delivery_matches_request": True,
    }


def test_result_indexing_projection_preserves_only_structured_state():
    diagnostic = _result_indexing()

    assert diagnostic.to_dict() == {
        "state": "qc_indexed",
        "artifact_attempt_status": "succeeded",
        "artifact_outcome": "indexed",
        "artifact_generation_present": True,
        "qc_attempt_status": "succeeded",
        "qc_outcome": "indexed",
        "qc_generation_present": True,
    }


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("state", "qc_indexed", "state is invalid"),
        ("artifact_attempt_status", "", "artifact_attempt_status"),
        ("artifact_outcome", 1, "artifact_outcome"),
        ("qc_attempt_status", "", "qc_attempt_status"),
        ("qc_outcome", 1, "qc_outcome"),
        ("artifact_generation_present", 1, "artifact_generation_present"),
        ("qc_generation_present", None, "qc_generation_present"),
    ],
)
def test_result_indexing_projection_rejects_corrupted_repository_data(
    field, value, message
):
    kwargs = _result_indexing_kwargs()
    kwargs[field] = value

    with pytest.raises(ValueError, match=message):
        RunResultIndexingDiagnostic(**kwargs)


def test_recovery_diagnostic_is_path_safe_and_exposes_bounded_actions():
    diagnostic = RunRecoveryDiagnostic(**_diagnostic_kwargs())

    payload = diagnostic.to_dict()
    rendered = json.dumps(payload, sort_keys=True)
    assert diagnostic.can_fail is True
    assert diagnostic.can_requeue is True
    assert payload["assignment"] == {
        "job_id": "job-1",
        "backend": "rq",
        "queue_name": "runs",
        "dispatched": True,
        "claimed": True,
        "cancellation_requested": True,
        "cancellation_acknowledged": True,
        "requeue_requested": True,
        "requeue_confirmed": True,
    }
    assert payload["gaps"] == ["callback", "cleanup"]
    assert payload["allowed_actions"] == ["fail", "requeue"]
    assert "/private/lab" not in rendered
    assert "cancellation_reason" not in rendered
    assert "operator request" not in repr(diagnostic)


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("run_id", " ", "^run_id must be a non-empty string$"),
        ("workflow_id", "", "workflow_id"),
        ("status", "running", "RunStatus"),
        ("diagnosis_code", "", "diagnosis_code"),
        ("assignment", object(), "assignment is invalid"),
        (
            "assignment",
            _assignment(run_id="other-run"),
            "assignment does not match run_id",
        ),
        ("queue_evidence", object(), "queue_evidence"),
        ("build_identity_present", 1, "build_identity_present"),
        ("result_indexing", object(), "result_indexing"),
        ("cleanup", "pending", "cleanup state"),
        ("gaps", ("queue",), "gaps contain"),
        ("gaps", (RunRecoveryGap.QUEUE, RunRecoveryGap.QUEUE), "duplicates"),
        ("allowed_actions", ("fail",), "allowed_actions contain"),
        (
            "allowed_actions",
            (RunRecoveryAction.FAIL, RunRecoveryAction.FAIL),
            "duplicates",
        ),
        ("issue_codes", ("",), "issue_codes"),
    ],
)
def test_recovery_diagnostic_rejects_invalid_public_contract_data(
    field, value, message
):
    kwargs = _diagnostic_kwargs()
    kwargs[field] = value

    with pytest.raises(ValueError, match=message):
        RunRecoveryDiagnostic(**kwargs)


def test_recovery_action_result_serializes_opaque_identity_and_timestamps():
    result = RunRecoveryActionResult(**_action_result_kwargs())

    assert result.to_dict() == {
        "run_id": "run-1",
        "action": "requeue",
        "previous_status": "queued",
        "status": "queued",
        "reason_code": "RUN_REQUEUED_BY_ADMIN_RECOVERY",
        "changed": True,
        "job_id": "job-1",
        "requeue_requested_at": "2026-08-09T12:00:05+00:00",
        "requeue_confirmed_at": "2026-08-09T12:00:06+00:00",
    }
    assert "backend" not in result.to_dict()
    assert "queue_name" not in result.to_dict()


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("run_id", "", "^run_id must be a non-empty string$"),
        ("action", "requeue", "RunRecoveryAction"),
        ("previous_status", "queued", "previous_status"),
        ("status", "queued", "status"),
        ("assignment", object(), "assignment does not match"),
        (
            "assignment",
            _assignment(run_id="other-run"),
            "assignment does not match",
        ),
        ("reason_code", "", "reason_code"),
        ("changed", 1, "changed"),
    ],
)
def test_recovery_action_result_rejects_invalid_committed_data(field, value, message):
    kwargs = _action_result_kwargs()
    kwargs[field] = value

    with pytest.raises(ValueError, match=message):
        RunRecoveryActionResult(**kwargs)


def test_recovery_summary_serializes_machine_readable_gap_counts():
    summary = RunRecoverySummary(**_summary_kwargs())

    assert summary.to_dict() == _summary_kwargs()


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("runs_examined", True, "nonnegative integer"),
        ("database_gaps", -1, "nonnegative integer"),
        ("queue_gaps", 1.5, "nonnegative integer"),
        ("callback_gaps", 4, "cannot exceed"),
        ("result_indexing_gaps", 4, "cannot exceed"),
        ("cleanup_gaps", 4, "cannot exceed"),
        ("queue_unavailable", 0, "must be a bool"),
    ],
)
def test_recovery_summary_rejects_impossible_gap_counts(field, value, message):
    kwargs = _summary_kwargs()
    kwargs[field] = value

    with pytest.raises(ValueError, match=message):
        RunRecoverySummary(**kwargs)
