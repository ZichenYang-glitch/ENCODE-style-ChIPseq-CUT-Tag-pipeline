"""Real cancellation and timeout acceptance for bulk RNA-seq lifecycle safety."""

from __future__ import annotations

from pathlib import Path
import time

import pytest

from .platform_harness import PlatformAcceptanceHarness, TerminalLifecycleEvidence
from .support import load_acceptance_fixture, require_gate_settings


pytestmark = pytest.mark.bulk_rnaseq_real_execution
REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
_CANCELLATION_REASON = "User requested cancellation."
_WORKFLOW_TIMEOUT_SECONDS = 20


def test_real_container_execution_is_truthfully_cancelled_and_reaped(
    tmp_path: Path,
) -> None:
    """Cancel a live Nextflow/container run through SQLite, Redis, and RQ."""
    settings = require_gate_settings()
    fixture = load_acceptance_fixture(settings.fixture_manifest)

    with PlatformAcceptanceHarness(
        gate_settings=settings,
        repository_root=REPOSITORY_ROOT,
        temporary_root=(tmp_path / "cancelled-platform").resolve(),
    ) as harness:
        submitted = harness.submit(fixture)
        worker = harness.start_worker()
        activity = harness.wait_for_execution_activity(
            submitted,
            worker,
            require_managed_container=True,
        )
        cancellation = harness.request_cancellation(
            submitted,
            reason=_CANCELLATION_REASON,
        )
        assert cancellation.stop_requested or (
            cancellation.record.status.value == "cancelled"
        )
        harness.wait_worker(worker, timeout_seconds=120)
        evidence = harness.collect_terminal(submitted)

    assert activity.nextflow_observed is True
    assert activity.managed_container_observed is True
    assert activity.process_group_count >= 1
    assert evidence.lifecycle_status == "cancelled"
    assert evidence.lifecycle_history[-3:] == (
        "queued",
        "running",
        "cancelled",
    )
    assert evidence.assignment_dispatched is True
    assert evidence.assignment_claimed is True
    assert evidence.cancellation_requested is True
    assert evidence.cancellation_acknowledged is True
    assert evidence.cancellation_reason == _CANCELLATION_REASON
    assert evidence.event_types.count("cancellation_requested") == 1
    assert evidence.event_types.count("cancellation_acknowledged") == 1
    assert evidence.event_types.count("execution_cleanup_failed") == 0
    assert evidence.error_code is None
    assert evidence.error_reason_code is None
    _assert_no_result_attempts(evidence)
    assert evidence.rq_status == "stopped"
    assert evidence.rq_stopped is True
    assert evidence.rq_failed is False
    assert evidence.rq_finished is False
    assert evidence.cleanup_confirmed is True


def test_real_nextflow_timeout_persists_failure_and_reaps_execution_tree(
    tmp_path: Path,
) -> None:
    """Apply the platform workflow timeout to a live pinned Nextflow process."""
    settings = require_gate_settings()
    fixture = load_acceptance_fixture(settings.fixture_manifest)

    with PlatformAcceptanceHarness(
        gate_settings=settings,
        repository_root=REPOSITORY_ROOT,
        temporary_root=(tmp_path / "timed-out-platform").resolve(),
        job_timeout_seconds=_WORKFLOW_TIMEOUT_SECONDS,
    ) as harness:
        assert harness.worker_settings.job_timeout_seconds == _WORKFLOW_TIMEOUT_SECONDS
        submitted = harness.submit(fixture)
        worker = harness.start_worker()
        activity = harness.wait_for_execution_activity(
            submitted,
            worker,
            require_managed_container=False,
        )
        activity_observed_at = time.monotonic()
        harness.wait_worker(worker)
        elapsed_after_activity = time.monotonic() - activity_observed_at
        evidence = harness.collect_terminal(submitted)

    assert activity.nextflow_observed is True
    assert activity.process_group_count >= 1
    assert elapsed_after_activity >= _WORKFLOW_TIMEOUT_SECONDS - 2
    assert evidence.lifecycle_status == "failed"
    assert evidence.lifecycle_history[-3:] == ("queued", "running", "failed")
    assert evidence.assignment_dispatched is True
    assert evidence.assignment_claimed is True
    assert evidence.cancellation_requested is False
    assert evidence.cancellation_acknowledged is False
    assert evidence.cancellation_reason is None
    assert evidence.event_types.count("execution_cleanup_failed") == 0
    assert evidence.error_code == "RUN_EXECUTION_FAILED"
    assert evidence.error_reason_code == "PROCESS_RUNNER_TIMEOUT"
    _assert_no_result_attempts(evidence)
    assert evidence.rq_status == "failed"
    assert evidence.rq_failed is True
    assert evidence.rq_stopped is False
    assert evidence.rq_finished is False
    assert evidence.cleanup_confirmed is True


def _assert_no_result_attempts(evidence: TerminalLifecycleEvidence) -> None:
    assert evidence.artifact_revision == 0
    assert evidence.artifact_attempt_id is None
    assert evidence.artifact_attempt_status is None
    assert evidence.qc_revision == 0
    assert evidence.qc_attempt_id is None
    assert evidence.qc_attempt_status is None
    assert evidence.artifact_count == 0
    assert evidence.qc_metric_count == 0
