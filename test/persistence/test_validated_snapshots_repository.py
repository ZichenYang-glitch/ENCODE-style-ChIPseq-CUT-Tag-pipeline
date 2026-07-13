from __future__ import annotations

from datetime import datetime, timedelta, timezone
from threading import Barrier, Thread

import pytest
from sqlalchemy import text

from encode_pipeline.persistence import open_run_persistence
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
    RunEventDraft,
    ValidatedSnapshotExpiredError,
    ValidatedSnapshotReplayConflictError,
)


NOW = datetime(2026, 7, 14, 10, 0, tzinfo=timezone.utc)


def _identity() -> WorkflowBuildIdentity:
    return WorkflowBuildIdentity(
        workflow_id="workflow-a",
        adapter_version="1.0.0",
        scheme="sha256-tree-v1",
        logical_entrypoint="workflow/Snakefile",
        digest="a" * 64,
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
        validation_issue_codes=("VALIDATION_WARNING",),
        validated_at=NOW,
        expires_at=NOW + timedelta(minutes=30),
    )


def _record(run_id: str, consumed_at: datetime, tags=None) -> RunRecord:
    return RunRecord(
        run_id=run_id,
        workflow_id="workflow-a",
        inputs=_snapshot().to_workflow_inputs().to_dict(),
        status=RunStatus.CREATED,
        created_at=consumed_at,
        updated_at=consumed_at,
        started_at=None,
        ended_at=None,
        current_stage=None,
        cancellation_reason=None,
        error=None,
        tags=tags or {},
    )


def _consume(repository, run_id: str, *, tags=None, consumed_at=None):
    timestamp = consumed_at or NOW + timedelta(minutes=1)
    return repository.consume_validated_input_snapshot(
        _snapshot().snapshot_id,
        workflow_id="workflow-a",
        expected_build_identity=_identity(),
        record=_record(run_id, timestamp, tags),
        consumed_at=timestamp,
        event=RunEventDraft(
            event_type="status_changed",
            message="Run created.",
            status=RunStatus.CREATED,
            context={"previous_status": None, "new_status": "created"},
        ),
    )


def test_sqlalchemy_snapshot_survives_restart_and_consumes_atomically(tmp_path) -> None:
    database_url = f"sqlite:///{tmp_path / 'platform.db'}"
    first = open_run_persistence(database_url)
    first.repository.create_validated_input_snapshot(_snapshot())
    first.close()

    second = open_run_persistence(database_url)
    loaded = second.repository.get_validated_input_snapshot(_snapshot().snapshot_id)
    assert loaded == _snapshot()

    creation = _consume(second.repository, "run-1", tags={"owner": "lab"})
    assert creation.created is True
    assert second.repository.get_run("run-1") == creation.record
    assert len(second.repository.list_events("run-1")) == 1
    consumed = second.repository.get_validated_input_snapshot(_snapshot().snapshot_id)
    assert consumed.consumed_run_id == "run-1"
    second.close()

    third = open_run_persistence(database_url)
    replay = _consume(third.repository, "run-other", tags={"owner": "lab"})
    assert replay.created is False
    assert replay.record.run_id == "run-1"
    assert len(third.repository.list_runs()) == 1
    third.close()


def test_sqlalchemy_snapshot_expiry_and_replay_conflict_roll_back(tmp_path) -> None:
    persistence = open_run_persistence(f"sqlite:///{tmp_path / 'platform.db'}")
    repository = persistence.repository
    repository.create_validated_input_snapshot(_snapshot())

    with pytest.raises(ValidatedSnapshotExpiredError):
        _consume(
            repository,
            "run-expired",
            consumed_at=NOW + timedelta(minutes=30),
        )
    assert repository.list_runs() == ()
    _consume(repository, "run-1", tags={"owner": "lab"})
    with pytest.raises(ValidatedSnapshotReplayConflictError):
        _consume(repository, "run-2", tags={"owner": "other"})
    assert [run.run_id for run in repository.list_runs()] == ["run-1"]
    persistence.close()


def test_sqlalchemy_snapshot_replay_rejects_corrupt_linked_run_timestamp(
    tmp_path,
) -> None:
    persistence = open_run_persistence(f"sqlite:///{tmp_path / 'platform.db'}")
    repository = persistence.repository
    repository.create_validated_input_snapshot(_snapshot())
    _consume(repository, "run-1", tags={"owner": "lab"})
    second_time = NOW + timedelta(minutes=2)
    repository.create_run(
        _record("run-2", second_time, tags={"owner": "lab"}),
        RunEventDraft(
            event_type="status_changed",
            message="Run created.",
            status=RunStatus.CREATED,
            context={"previous_status": None, "new_status": "created"},
        ),
    )
    with persistence.engine.begin() as connection:
        connection.execute(
            text(
                "UPDATE validated_input_snapshots SET consumed_run_id = :run_id "
                "WHERE snapshot_id = :snapshot_id"
            ),
            {"run_id": "run-2", "snapshot_id": _snapshot().snapshot_id},
        )

    with pytest.raises(ValueError, match="consumption time"):
        _consume(repository, "run-other", tags={"owner": "lab"})

    persistence.close()


def test_sqlalchemy_snapshot_corruption_fails_before_run_creation(tmp_path) -> None:
    persistence = open_run_persistence(f"sqlite:///{tmp_path / 'platform.db'}")
    repository = persistence.repository
    repository.create_validated_input_snapshot(_snapshot())
    with persistence.engine.begin() as connection:
        connection.execute(
            text(
                "UPDATE validated_input_snapshots SET payload_digest = :digest "
                "WHERE snapshot_id = :snapshot_id"
            ),
            {"digest": "b" * 64, "snapshot_id": _snapshot().snapshot_id},
        )

    with pytest.raises(ValueError, match="digest"):
        _consume(repository, "run-1")

    assert repository.list_runs() == ()
    persistence.close()


def test_two_sqlalchemy_repositories_converge_on_one_snapshot_run(tmp_path) -> None:
    database_url = f"sqlite:///{tmp_path / 'platform.db'}"
    first = open_run_persistence(database_url)
    second = open_run_persistence(database_url)
    first.repository.create_validated_input_snapshot(_snapshot())
    barrier = Barrier(3)
    results = []
    failures = []

    def consume(repository, run_id: str) -> None:
        barrier.wait()
        try:
            results.append(_consume(repository, run_id))
        except BaseException as error:  # pragma: no cover - assertion reports it
            failures.append(error)

    threads = [
        Thread(target=consume, args=(first.repository, "run-1")),
        Thread(target=consume, args=(second.repository, "run-2")),
    ]
    for thread in threads:
        thread.start()
    barrier.wait()
    for thread in threads:
        thread.join(timeout=5)

    assert failures == []
    assert sum(result.created for result in results) == 1
    assert len({result.record.run_id for result in results}) == 1
    assert len(first.repository.list_runs()) == 1
    first.close()
    second.close()
