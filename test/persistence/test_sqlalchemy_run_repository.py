"""Integration tests for SQLAlchemy-backed workflow run persistence."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from dataclasses import replace
from datetime import datetime, timezone

import pytest

from encode_pipeline.persistence import (
    SqlAlchemyRunRepository,
    create_database_engine,
    create_session_factory,
    upgrade_database,
)
from encode_pipeline.platform.adapters import WorkflowInputs
from encode_pipeline.platform.results import Issue
from encode_pipeline.platform.runs import RunArtifactRef, RunRecord, RunStatus
from encode_pipeline.services.defaults import create_default_workflow_registry
from encode_pipeline.services.run_repositories import (
    ConcurrentRunUpdateError,
    RunEventDraft,
)
from encode_pipeline.services.runs import RunService


WORKFLOW_ID = "encode-style-chipseq-cuttag-atac-mnase"


@pytest.fixture
def database_url(tmp_path):
    return f"sqlite:///{tmp_path / 'platform.db'}"


@pytest.fixture
def repository(database_url):
    upgrade_database(database_url)
    engine = create_database_engine(database_url)
    repo = SqlAlchemyRunRepository(create_session_factory(engine))
    yield repo
    engine.dispose()


def test_run_service_round_trip_survives_repository_reopen(database_url):
    upgrade_database(database_url)
    first_engine = create_database_engine(database_url)
    first_service = RunService(
        registry=create_default_workflow_registry(),
        id_factory=lambda: "run-1",
        repository=SqlAlchemyRunRepository(create_session_factory(first_engine)),
    )
    inputs = WorkflowInputs(
        config={"assay": "chipseq"},
        samples=[{"id": "S1", "fastq_1": "fixtures/S1.fastq.gz"}],
        options={"cores": 2},
    )
    created = first_service.create_run(WORKFLOW_ID, inputs, tags={"owner": "test"})
    validating = first_service.transition_run(
        created.run_id,
        RunStatus.VALIDATING,
        stage="validation",
    )
    issue = Issue(
        code="INPUT_WARNING",
        message="Input is usable with a warning.",
        severity="warning",
        source="adapter",
        technical_message="private detail",
        context={"sample": "S1"},
    )
    first_service.add_event(
        created.run_id,
        "issue_added",
        "Validation warning recorded.",
        status=RunStatus.VALIDATING,
        stage="validation",
        issue=issue,
    )
    stdout = first_service.append_log(created.run_id, "stdout", ["line 1"])
    stderr = first_service.append_log(created.run_id, "stderr", ["warning"])
    artifact = RunArtifactRef(
        artifact_id="artifact-1",
        run_id=created.run_id,
        artifact_type="file",
        name="summary.json",
        uri="run://runs/run-1/artifacts/summary.json",
        mime_type="application/json",
        produced_at=datetime.now(timezone.utc),
        metadata={"sample": "S1"},
    )
    first_service.record_artifact(created.run_id, artifact)
    first_engine.dispose()

    second_engine = create_database_engine(database_url)
    second_service = RunService(
        registry=create_default_workflow_registry(),
        repository=SqlAlchemyRunRepository(create_session_factory(second_engine)),
    )

    assert second_service.get_run("run-1") == validating
    events = second_service.list_events("run-1")
    assert [event.sequence for event in events] == [1, 2, 3]
    assert events[-1].issue == issue
    assert second_service.list_events("run-1", after="evt-1", limit=1) == (events[1],)
    assert second_service.list_logs("run-1", "stdout") == (stdout,)
    assert second_service.list_logs("run-1", "stderr") == (stderr,)
    assert second_service.list_artifacts("run-1") == (artifact,)
    second_engine.dispose()


def test_repository_preserves_creation_and_artifact_insertion_order(repository):
    repository.create_run(_record("run-2"), _created_event())
    repository.create_run(_record("run-1"), _created_event())
    artifact_2 = _artifact("run-2", "artifact-2")
    artifact_1 = _artifact("run-2", "artifact-1")
    repository.record_artifact("run-2", artifact_2)
    repository.record_artifact("run-2", artifact_1)

    assert [record.run_id for record in repository.list_runs()] == ["run-2", "run-1"]
    assert repository.list_artifacts("run-2") == (artifact_2, artifact_1)


def test_repository_preserves_key_and_cursor_errors(repository):
    with pytest.raises(KeyError):
        repository.get_run("missing")
    with pytest.raises(KeyError):
        repository.list_events("missing")
    with pytest.raises(KeyError):
        repository.list_logs("missing")
    with pytest.raises(KeyError):
        repository.list_artifacts("missing")

    repository.create_run(_record("run-1"), _created_event())
    with pytest.raises(KeyError):
        repository.list_events("run-1", after="evt-missing")
    with pytest.raises(KeyError):
        repository.list_logs("run-1", after="log-missing")


def test_repository_rejects_duplicate_public_ids(repository):
    record = _record("run-1")
    repository.create_run(record, _created_event())

    with pytest.raises(ValueError, match="Duplicate run_id"):
        repository.create_run(record, _created_event())

    artifact = _artifact("run-1", "artifact-1")
    repository.record_artifact("run-1", artifact)
    with pytest.raises(ValueError, match="Duplicate artifact_id"):
        repository.record_artifact("run-1", artifact)


def test_update_rejects_stale_status(repository):
    record = _record("run-1")
    repository.create_run(record, _created_event())
    updated = replace(record, status=RunStatus.VALIDATING)

    with pytest.raises(ConcurrentRunUpdateError):
        repository.update_run(
            updated,
            expected_status=RunStatus.PLANNED,
            event=RunEventDraft(
                event_type="status_changed",
                message="Stale update.",
                status=RunStatus.VALIDATING,
            ),
        )

    assert repository.get_run("run-1") == record
    assert len(repository.list_events("run-1")) == 1


def test_run_update_and_event_roll_back_together(repository, monkeypatch):
    record = _record("run-1")
    repository.create_run(record, _created_event())
    updated = replace(
        record,
        status=RunStatus.VALIDATING,
        updated_at=datetime.now(timezone.utc),
    )

    def fail_event(_session, _run_id, _draft):
        raise RuntimeError("event insert failed")

    monkeypatch.setattr(repository, "_insert_event", fail_event)
    with pytest.raises(RuntimeError, match="event insert failed"):
        repository.update_run(
            updated,
            expected_status=RunStatus.CREATED,
            event=RunEventDraft(
                event_type="status_changed",
                message="Run validating.",
                status=RunStatus.VALIDATING,
            ),
        )

    assert repository.get_run("run-1") == record
    assert len(repository.list_events("run-1")) == 1


def test_concurrent_appends_keep_unique_monotonic_sequences(repository):
    repository.create_run(_record("run-1"), _created_event())

    with ThreadPoolExecutor(max_workers=8) as pool:
        events = tuple(
            pool.map(
                lambda number: repository.add_event(
                    "run-1",
                    RunEventDraft(
                        event_type="progress",
                        message=f"event {number}",
                    ),
                ),
                range(24),
            )
        )
        logs = tuple(
            pool.map(
                lambda number: repository.append_log(
                    "run-1", "stdout", [f"line {number}"]
                ),
                range(24),
            )
        )

    assert sorted(event.sequence for event in events) == list(range(2, 26))
    assert [event.sequence for event in repository.list_events("run-1")] == list(
        range(1, 26)
    )
    assert sorted(chunk.sequence for chunk in logs) == list(range(1, 25))
    assert [chunk.sequence for chunk in repository.list_logs("run-1")] == list(
        range(1, 25)
    )


def _record(run_id: str) -> RunRecord:
    now = datetime.now(timezone.utc)
    return RunRecord(
        run_id=run_id,
        workflow_id=WORKFLOW_ID,
        inputs={"config": {}, "samples": None, "options": {}},
        status=RunStatus.CREATED,
        created_at=now,
        updated_at=now,
        started_at=None,
        ended_at=None,
        current_stage=None,
        cancellation_reason=None,
        error=None,
        tags={"test": "true"},
    )


def _created_event() -> RunEventDraft:
    return RunEventDraft(
        event_type="status_changed",
        message="Run created.",
        status=RunStatus.CREATED,
    )


def _artifact(run_id: str, artifact_id: str) -> RunArtifactRef:
    return RunArtifactRef(
        artifact_id=artifact_id,
        run_id=run_id,
        artifact_type="file",
        name=f"{artifact_id}.json",
        uri=f"run://runs/{run_id}/artifacts/{artifact_id}.json",
        mime_type="application/json",
        produced_at=datetime.now(timezone.utc),
        metadata={},
    )
