from __future__ import annotations

from dataclasses import replace
from datetime import datetime, timedelta, timezone
from threading import Barrier, Thread

import pytest

from encode_pipeline.platform.adapters import WorkflowInputs
from encode_pipeline.platform.builds import WorkflowBuildIdentity
from encode_pipeline.platform.runs import RunRecord, RunStatus
from encode_pipeline.platform.snapshots import (
    PAYLOAD_DIGEST_SCHEME,
    VALIDATION_EVIDENCE_OUTCOME,
    ValidatedInputSnapshot,
    build_workflow_inputs_digest,
    canonical_workflow_inputs_json,
)
from encode_pipeline.services.run_repositories import (
    InMemoryRunRepository,
    RunEventDraft,
    ValidatedSnapshotBuildMismatchError,
    ValidatedSnapshotExpiredError,
    ValidatedSnapshotReplayConflictError,
)


NOW = datetime(2026, 7, 14, 9, 0, tzinfo=timezone.utc)


def _identity(digest: str = "a" * 64) -> WorkflowBuildIdentity:
    return WorkflowBuildIdentity(
        workflow_id="workflow-a",
        adapter_version="1.0.0",
        scheme="sha256-tree-v1",
        logical_entrypoint="workflow/Snakefile",
        digest=digest,
        captured_at=NOW,
    )


def _snapshot() -> ValidatedInputSnapshot:
    payload = canonical_workflow_inputs_json(
        WorkflowInputs(
            config={"threads": 1},
            samples=[{"sample": "S1"}],
            options={},
        )
    )
    return ValidatedInputSnapshot(
        snapshot_id="vsnap_0123456789abcdef0123456789abcdef",
        workflow_id="workflow-a",
        adapter_version="1.0.0",
        schema_version="1.0.0",
        schema_dialect="https://json-schema.org/draft/2020-12/schema",
        workflow_build_identity=_identity(),
        canonical_payload=payload,
        payload_digest_scheme=PAYLOAD_DIGEST_SCHEME,
        payload_digest=build_workflow_inputs_digest(payload),
        validation_outcome=VALIDATION_EVIDENCE_OUTCOME,
        validation_issue_codes=(),
        validated_at=NOW,
        expires_at=NOW + timedelta(minutes=30),
    )


def _record(run_id: str, *, tags=None) -> RunRecord:
    return RunRecord(
        run_id=run_id,
        workflow_id="workflow-a",
        inputs=_snapshot().to_workflow_inputs().to_dict(),
        status=RunStatus.CREATED,
        created_at=NOW + timedelta(minutes=1),
        updated_at=NOW + timedelta(minutes=1),
        started_at=None,
        ended_at=None,
        current_stage=None,
        cancellation_reason=None,
        error=None,
        tags=tags or {},
    )


def _event() -> RunEventDraft:
    return RunEventDraft(
        event_type="status_changed",
        message="Run created.",
        status=RunStatus.CREATED,
        context={"previous_status": None, "new_status": "created"},
    )


def _consume(repository, record, **changes):
    arguments = {
        "workflow_id": "workflow-a",
        "expected_build_identity": _identity(),
        "record": record,
        "consumed_at": NOW + timedelta(minutes=1),
        "event": _event(),
    }
    arguments.update(changes)
    record = replace(
        record,
        created_at=arguments["consumed_at"],
        updated_at=arguments["consumed_at"],
    )
    arguments["record"] = record
    return repository.consume_validated_input_snapshot(
        _snapshot().snapshot_id,
        **arguments,
    )


def test_in_memory_snapshot_create_read_is_defensive() -> None:
    repository = InMemoryRunRepository()
    snapshot = _snapshot()

    repository.create_validated_input_snapshot(snapshot)
    loaded = repository.get_validated_input_snapshot(snapshot.snapshot_id)

    assert loaded == snapshot
    assert loaded is not snapshot
    with pytest.raises(ValueError, match="Duplicate"):
        repository.create_validated_input_snapshot(snapshot)


def test_in_memory_snapshot_first_use_atomically_creates_run_and_event() -> None:
    repository = InMemoryRunRepository()
    snapshot = _snapshot()
    repository.create_validated_input_snapshot(snapshot)

    creation = _consume(repository, _record("run-1"))

    assert creation.created is True
    assert creation.record.run_id == "run-1"
    assert repository.get_run("run-1") == creation.record
    assert len(repository.list_events("run-1")) == 1
    consumed = repository.get_validated_input_snapshot(snapshot.snapshot_id)
    assert consumed.consumed_run_id == "run-1"
    assert consumed.consumed_at == NOW + timedelta(minutes=1)


def test_in_memory_snapshot_identical_replay_returns_canonical_advanced_run() -> None:
    repository = InMemoryRunRepository()
    repository.create_validated_input_snapshot(_snapshot())
    first = _consume(repository, _record("run-1", tags={"owner": "lab"}))
    planned = replace(
        first.record,
        status=RunStatus.VALIDATING,
        updated_at=NOW + timedelta(minutes=2),
    )
    repository.update_run(
        planned,
        expected_status=RunStatus.CREATED,
        event=RunEventDraft(
            event_type="status_changed",
            message="Validation started.",
            status=RunStatus.VALIDATING,
        ),
    )

    replay = _consume(repository, _record("run-other", tags={"owner": "lab"}))

    assert replay.created is False
    assert replay.record == planned
    assert len(repository.list_runs()) == 1
    assert len(repository.list_events("run-1")) == 2


def test_in_memory_snapshot_replay_with_different_tags_is_rejected() -> None:
    repository = InMemoryRunRepository()
    repository.create_validated_input_snapshot(_snapshot())
    _consume(repository, _record("run-1", tags={"owner": "lab"}))

    with pytest.raises(ValidatedSnapshotReplayConflictError):
        _consume(repository, _record("run-2", tags={"owner": "other"}))

    assert [run.run_id for run in repository.list_runs()] == ["run-1"]


def test_in_memory_snapshot_expiry_and_build_mismatch_create_nothing() -> None:
    repository = InMemoryRunRepository()
    repository.create_validated_input_snapshot(_snapshot())

    with pytest.raises(ValidatedSnapshotExpiredError):
        _consume(
            repository,
            _record("run-expired"),
            consumed_at=NOW + timedelta(minutes=30),
        )
    with pytest.raises(ValidatedSnapshotBuildMismatchError):
        _consume(
            repository,
            _record("run-stale"),
            expected_build_identity=_identity("b" * 64),
        )

    assert repository.list_runs() == ()
    assert (
        repository.get_validated_input_snapshot(_snapshot().snapshot_id).consumed_run_id
        is None
    )


def test_in_memory_snapshot_cross_workflow_is_indistinguishable_from_missing() -> None:
    repository = InMemoryRunRepository()
    repository.create_validated_input_snapshot(_snapshot())

    with pytest.raises(KeyError):
        _consume(repository, _record("run-1"), workflow_id="workflow-b")

    assert repository.list_runs() == ()


def test_in_memory_snapshot_rejects_tampered_storage_before_run_creation() -> None:
    repository = InMemoryRunRepository()
    snapshot = _snapshot()
    repository.create_validated_input_snapshot(snapshot)
    stored = repository._validated_input_snapshots[snapshot.snapshot_id]
    object.__setattr__(stored, "payload_digest", "b" * 64)

    with pytest.raises(ValueError, match="digest"):
        _consume(repository, _record("run-1"))

    assert repository.list_runs() == ()


def test_in_memory_concurrent_snapshot_consumption_converges_on_one_run() -> None:
    repository = InMemoryRunRepository()
    repository.create_validated_input_snapshot(_snapshot())
    barrier = Barrier(3)
    results = []
    failures = []

    def consume(run_id: str) -> None:
        barrier.wait()
        try:
            results.append(_consume(repository, _record(run_id)))
        except BaseException as error:  # pragma: no cover - assertion reports it
            failures.append(error)

    threads = [Thread(target=consume, args=(f"run-{index}",)) for index in (1, 2)]
    for thread in threads:
        thread.start()
    barrier.wait()
    for thread in threads:
        thread.join(timeout=3)

    assert failures == []
    assert len(results) == 2
    assert sum(result.created for result in results) == 1
    assert len({result.record.run_id for result in results}) == 1
    assert len(repository.list_runs()) == 1
