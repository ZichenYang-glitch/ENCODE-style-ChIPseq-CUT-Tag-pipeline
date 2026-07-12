"""Integration tests for SQLAlchemy-backed workflow run persistence."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from dataclasses import replace
from datetime import datetime, timedelta, timezone
from threading import Barrier

import pytest

from encode_pipeline.persistence import (
    SqlAlchemyRunRepository,
    create_database_engine,
    create_session_factory,
    upgrade_database,
)
from encode_pipeline.platform.adapters import WorkflowInputs
from encode_pipeline.platform.builds import WorkflowBuildIdentity
from encode_pipeline.platform.execution import RunExecutionAssignment
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


def test_complete_preflight_identity_run_and_event_round_trip(repository):
    created = _record("run-1")
    validating = replace(created, status=RunStatus.VALIDATING)
    repository.create_run(created, _created_event())
    repository.update_run(
        validating,
        expected_status=RunStatus.CREATED,
        event=RunEventDraft(
            event_type="status_changed",
            message="Run validating.",
            status=RunStatus.VALIDATING,
        ),
    )
    planned = replace(validating, status=RunStatus.PLANNED)
    identity = _build_identity()

    repository.complete_preflight(
        planned,
        identity,
        expected_status=RunStatus.VALIDATING,
        event=RunEventDraft(
            event_type="preflight_completed",
            message="Preflight complete.",
            status=RunStatus.PLANNED,
        ),
    )

    assert repository.get_run("run-1") == planned
    assert repository.get_workflow_build_identity("run-1") == identity
    assert repository.list_events("run-1")[-1].status is RunStatus.PLANNED


def test_complete_preflight_rolls_back_identity_and_run_when_event_fails(
    repository,
    monkeypatch,
):
    validating = replace(_record("run-1"), status=RunStatus.VALIDATING)
    repository.create_run(validating, _created_event())

    def fail_event(_session, _run_id, _draft):
        raise RuntimeError("event insert failed")

    monkeypatch.setattr(repository, "_insert_event", fail_event)
    with pytest.raises(RuntimeError, match="event insert failed"):
        repository.complete_preflight(
            replace(validating, status=RunStatus.PLANNED),
            _build_identity(),
            expected_status=RunStatus.VALIDATING,
            event=RunEventDraft(
                event_type="preflight_completed",
                message="Preflight complete.",
                status=RunStatus.PLANNED,
            ),
        )

    assert repository.get_run("run-1") == validating
    assert repository.get_workflow_build_identity("run-1") is None


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


def test_independent_repositories_serialize_event_and_log_sequences(database_url):
    upgrade_database(database_url)
    first_engine = create_database_engine(database_url)
    second_engine = create_database_engine(database_url)
    first_repository = SqlAlchemyRunRepository(create_session_factory(first_engine))
    second_repository = SqlAlchemyRunRepository(create_session_factory(second_engine))
    first_repository.create_run(_record("run-1"), _created_event())
    event_barrier = Barrier(2)

    def append_event(repository, label):
        event_barrier.wait(timeout=5)
        return repository.add_event(
            "run-1",
            RunEventDraft(event_type=label, message=label),
        )

    try:
        with ThreadPoolExecutor(max_workers=2) as pool:
            futures = (
                pool.submit(append_event, first_repository, "first"),
                pool.submit(append_event, second_repository, "second"),
            )
            for future in futures:
                future.result()

        log_barrier = Barrier(2)

        def append_log(repository, label):
            log_barrier.wait(timeout=5)
            return repository.append_log("run-1", "stdout", [label])

        with ThreadPoolExecutor(max_workers=2) as pool:
            futures = (
                pool.submit(append_log, first_repository, "first"),
                pool.submit(append_log, second_repository, "second"),
            )
            for future in futures:
                future.result()

        events = first_repository.list_events("run-1", limit=10)
        logs = second_repository.list_logs("run-1", "stdout", limit=10)
        assert [event.sequence for event in events] == [1, 2, 3]
        assert {event.event_type for event in events[1:]} == {"first", "second"}
        assert [chunk.sequence for chunk in logs] == [1, 2]
        assert {chunk.lines for chunk in logs} == {("first",), ("second",)}
    finally:
        first_engine.dispose()
        second_engine.dispose()


def test_execution_assignment_round_trips_and_returns_canonical_existing(
    repository,
):
    repository.create_run(_record("run-1"), _created_event())
    original = _assignment("run-1", "job-original")
    replacement = replace(
        original,
        job_id="job-replacement",
        queue_name="priority",
        created_at=original.created_at + timedelta(seconds=1),
    )

    assert repository.get_execution_assignment("run-1") is None
    assert (
        repository.ensure_execution_assignment(
            original, expected_status=RunStatus.CREATED
        )
        == original
    )
    assert (
        repository.ensure_execution_assignment(
            replacement, expected_status=RunStatus.CREATED
        )
        == original
    )
    assert repository.get_execution_assignment("run-1") == original


def test_dispatch_mark_is_idempotent_after_status_changes(repository):
    record = _record("run-1")
    repository.create_run(record, _created_event())
    assignment = repository.ensure_execution_assignment(
        _assignment(record.run_id, "job-1"),
        expected_status=RunStatus.CREATED,
    )
    dispatched_at = datetime.now(timezone.utc)
    dispatched = repository.mark_execution_dispatched(
        record.run_id,
        job_id=assignment.job_id,
        dispatched_at=dispatched_at,
        allowed_statuses=frozenset({RunStatus.CREATED}),
    )
    repository.update_run(
        replace(record, status=RunStatus.VALIDATING),
        expected_status=RunStatus.CREATED,
        event=RunEventDraft(
            event_type="status_changed",
            message="Run advanced.",
            status=RunStatus.VALIDATING,
        ),
    )

    retried = repository.mark_execution_dispatched(
        record.run_id,
        job_id=assignment.job_id,
        dispatched_at=dispatched_at + timedelta(seconds=1),
        allowed_statuses=frozenset({RunStatus.CREATED}),
    )

    assert retried == dispatched


def test_queue_dispatched_run_is_atomic_and_idempotent(repository):
    planned = replace(_record("run-1"), status=RunStatus.PLANNED)
    repository.create_run(planned, _created_event())
    assignment = repository.ensure_execution_assignment(
        _assignment(planned.run_id, "job-1"),
        expected_status=RunStatus.PLANNED,
    )
    queued = replace(
        planned,
        status=RunStatus.QUEUED,
        updated_at=datetime.now(timezone.utc),
        current_stage="execution",
    )
    event = RunEventDraft(
        event_type="status_changed",
        message="Run queued.",
        status=RunStatus.QUEUED,
    )

    with pytest.raises(ValueError, match="has not been dispatched"):
        repository.queue_dispatched_run(
            queued,
            expected_status=RunStatus.PLANNED,
            job_id=assignment.job_id,
            backend=assignment.backend,
            queue_name=assignment.queue_name,
            event=event,
        )

    repository.mark_execution_dispatched(
        planned.run_id,
        job_id=assignment.job_id,
        dispatched_at=datetime.now(timezone.utc),
        allowed_statuses=frozenset({RunStatus.PLANNED}),
    )
    for backend, queue_name in (
        ("other", assignment.queue_name),
        (assignment.backend, "other"),
    ):
        with pytest.raises(ValueError, match="identity does not match"):
            repository.queue_dispatched_run(
                queued,
                expected_status=RunStatus.PLANNED,
                job_id=assignment.job_id,
                backend=backend,
                queue_name=queue_name,
                event=event,
            )
    assert repository.get_run(planned.run_id).status is RunStatus.PLANNED

    assert repository.queue_dispatched_run(
        queued,
        expected_status=RunStatus.PLANNED,
        job_id=assignment.job_id,
        backend=assignment.backend,
        queue_name=assignment.queue_name,
        event=event,
    )
    assert not repository.queue_dispatched_run(
        queued,
        expected_status=RunStatus.PLANNED,
        job_id=assignment.job_id,
        backend=assignment.backend,
        queue_name=assignment.queue_name,
        event=event,
    )

    assert repository.get_run(planned.run_id) == queued
    events = repository.list_events(planned.run_id)
    assert [item.status for item in events].count(RunStatus.QUEUED) == 1


def test_execution_assignment_rejects_cross_run_job_reuse(repository):
    repository.create_run(_record("run-1"), _created_event())
    repository.create_run(_record("run-2"), _created_event())
    repository.ensure_execution_assignment(
        _assignment("run-1", "shared-job"),
        expected_status=RunStatus.CREATED,
    )

    with pytest.raises(ValueError, match="shared-job.*already assigned"):
        repository.ensure_execution_assignment(
            _assignment("run-2", "shared-job"),
            expected_status=RunStatus.CREATED,
        )


def test_execution_assignment_requires_a_persisted_run(repository):
    with pytest.raises(KeyError, match="missing"):
        repository.ensure_execution_assignment(
            _assignment("missing", "job-1"),
            expected_status=RunStatus.CREATED,
        )
    assert repository.get_execution_assignment("missing") is None


def test_independent_repositories_concurrently_choose_one_assignment(database_url):
    upgrade_database(database_url)
    first_engine = create_database_engine(database_url)
    second_engine = create_database_engine(database_url)
    first_repository = SqlAlchemyRunRepository(create_session_factory(first_engine))
    second_repository = SqlAlchemyRunRepository(create_session_factory(second_engine))
    first_repository.create_run(_record("run-1"), _created_event())
    first = _assignment("run-1", "job-first")
    second = replace(
        first,
        job_id="job-second",
        created_at=first.created_at + timedelta(seconds=1),
    )
    barrier = Barrier(2)

    def ensure(repository, assignment):
        barrier.wait()
        return repository.ensure_execution_assignment(
            assignment,
            expected_status=RunStatus.CREATED,
        )

    try:
        with ThreadPoolExecutor(max_workers=2) as pool:
            futures = (
                pool.submit(ensure, first_repository, first),
                pool.submit(ensure, second_repository, second),
            )
            resolved = tuple(future.result() for future in futures)

        assert resolved[0] == resolved[1]
        assert resolved[0] in {first, second}
        assert first_repository.get_execution_assignment("run-1") == resolved[0]
        assert second_repository.get_execution_assignment("run-1") == resolved[0]
    finally:
        first_engine.dispose()
        second_engine.dispose()


def test_independent_services_concurrently_queue_one_status_event(database_url):
    upgrade_database(database_url)
    first_engine = create_database_engine(database_url)
    second_engine = create_database_engine(database_url)
    first_service = RunService(
        create_default_workflow_registry(),
        id_factory=lambda: "run-1",
        repository=SqlAlchemyRunRepository(create_session_factory(first_engine)),
    )
    second_service = RunService(
        create_default_workflow_registry(),
        repository=SqlAlchemyRunRepository(create_session_factory(second_engine)),
    )
    first_service.create_run(WORKFLOW_ID, WorkflowInputs(config={}))
    first_service.transition_run("run-1", RunStatus.VALIDATING)
    first_service.transition_run("run-1", RunStatus.PLANNED)
    assignment = first_service.ensure_execution_assignment(
        "run-1",
        queue_name="default",
    )
    first_service.mark_execution_dispatched("run-1", job_id=assignment.job_id)
    barrier = Barrier(2)

    def queue(service):
        barrier.wait(timeout=5)
        return service.queue_dispatched_run(
            "run-1",
            job_id=assignment.job_id,
            backend=assignment.backend,
            queue_name=assignment.queue_name,
        )

    try:
        with ThreadPoolExecutor(max_workers=2) as pool:
            first = pool.submit(queue, first_service)
            second = pool.submit(queue, second_service)
            results = (first.result(), second.result())

        assert all(record.status is RunStatus.QUEUED for record in results)
        events = first_service.list_events("run-1")
        queued_events = [
            event
            for event in events
            if event.status is RunStatus.QUEUED
            and event.context.get("new_status") == RunStatus.QUEUED.value
        ]
        assert len(queued_events) == 1
    finally:
        first_engine.dispose()
        second_engine.dispose()


def test_recovery_and_dispatch_race_preserves_exactly_one_outcome(database_url):
    upgrade_database(database_url)
    first_engine = create_database_engine(database_url)
    second_engine = create_database_engine(database_url)
    first_service = RunService(
        create_default_workflow_registry(),
        id_factory=lambda: "run-1",
        repository=SqlAlchemyRunRepository(create_session_factory(first_engine)),
    )
    second_service = RunService(
        create_default_workflow_registry(),
        repository=SqlAlchemyRunRepository(create_session_factory(second_engine)),
    )
    first_service.create_run(WORKFLOW_ID, WorkflowInputs(config={}))
    first_service.transition_run("run-1", RunStatus.VALIDATING)
    first_service.transition_run("run-1", RunStatus.PLANNED)
    assignment = first_service.ensure_execution_assignment(
        "run-1",
        queue_name="default",
    )
    first_service.transition_run("run-1", RunStatus.QUEUED)
    barrier = Barrier(2)

    def recover():
        barrier.wait(timeout=5)
        return first_service.recover_interrupted_runs()

    def dispatch():
        barrier.wait(timeout=5)
        try:
            return second_service.mark_execution_dispatched(
                "run-1",
                job_id=assignment.job_id,
            )
        except ConcurrentRunUpdateError:
            return None

    try:
        with ThreadPoolExecutor(max_workers=2) as pool:
            recovery_future = pool.submit(recover)
            dispatch_future = pool.submit(dispatch)
            recovered = recovery_future.result()
            dispatched = dispatch_future.result()

        final_record = second_service.get_run("run-1")
        final_assignment = second_service.get_execution_assignment("run-1")
        assert final_assignment is not None
        if final_record.status is RunStatus.QUEUED:
            assert recovered == ()
            assert dispatched is not None
            assert final_assignment.dispatched_at is not None
        else:
            assert final_record.status is RunStatus.FAILED
            assert [record.run_id for record in recovered] == ["run-1"]
            assert dispatched is None
            assert final_assignment.dispatched_at is None
        assert second_service.list_events("run-1")[-1].status is final_record.status
    finally:
        first_engine.dispose()
        second_engine.dispose()


def test_worker_claim_and_cancel_race_never_writes_after_cancellation(database_url):
    upgrade_database(database_url)
    first_engine = create_database_engine(database_url)
    second_engine = create_database_engine(database_url)
    first_service = RunService(
        create_default_workflow_registry(),
        id_factory=lambda: "run-1",
        repository=SqlAlchemyRunRepository(create_session_factory(first_engine)),
    )
    second_service = RunService(
        create_default_workflow_registry(),
        repository=SqlAlchemyRunRepository(create_session_factory(second_engine)),
    )
    first_service.create_run(WORKFLOW_ID, WorkflowInputs(config={}))
    first_service.transition_run("run-1", RunStatus.VALIDATING)
    first_service.transition_run("run-1", RunStatus.PLANNED)
    assignment = first_service.ensure_execution_assignment(
        "run-1",
        queue_name="default",
    )
    first_service.mark_execution_dispatched("run-1", job_id=assignment.job_id)
    first_service.transition_run("run-1", RunStatus.QUEUED)
    barrier = Barrier(2)

    def claim():
        barrier.wait(timeout=5)
        try:
            return first_service.claim_execution_assignment(
                "run-1",
                job_id=assignment.job_id,
                backend=assignment.backend,
                queue_name=assignment.queue_name,
            )
        except ConcurrentRunUpdateError:
            return None

    def cancel():
        barrier.wait(timeout=5)
        return second_service.cancel_run("run-1", reason="race test")

    try:
        with ThreadPoolExecutor(max_workers=2) as pool:
            claim_future = pool.submit(claim)
            cancel_future = pool.submit(cancel)
            claim_result = claim_future.result()
            cancelled = cancel_future.result()

        assert cancelled.status is RunStatus.CANCELLED
        assert second_service.get_run("run-1").status is RunStatus.CANCELLED
        events = second_service.list_events("run-1", limit=100)
        assert events[-1].status is RunStatus.CANCELLED
        handshake_sequences = [
            event.sequence
            for event in events
            if event.event_type == "worker_dependencies_rebuilt"
        ]
        if claim_result is None:
            assert handshake_sequences == []
        else:
            assert claim_result.acquired is True
            assert handshake_sequences[-1] < events[-1].sequence
    finally:
        first_engine.dispose()
        second_engine.dispose()


def test_sqlalchemy_cancellation_intent_and_acknowledgement_survive_reopen(
    database_url,
):
    upgrade_database(database_url)
    first_engine = create_database_engine(database_url)
    first_service, assignment = _running_service(first_engine)

    requested = first_service.request_execution_cancellation(
        "run-1",
        job_id=assignment.job_id,
        backend=assignment.backend,
        queue_name=assignment.queue_name,
        reason="User requested cancellation.",
    )
    assert requested.record.status is RunStatus.RUNNING
    first_engine.dispose()

    second_engine = create_database_engine(database_url)
    second_service = RunService(
        create_default_workflow_registry(),
        repository=SqlAlchemyRunRepository(create_session_factory(second_engine)),
    )
    persisted = second_service.get_execution_assignment("run-1")
    assert persisted is not None
    assert persisted.cancellation_requested_at is not None
    assert persisted.cancellation_reason == "User requested cancellation."
    assert persisted.cancellation_acknowledged_at is None
    assert second_service.get_run("run-1").status is RunStatus.RUNNING

    acknowledged = second_service.acknowledge_execution_stop(
        "run-1",
        job_id=persisted.job_id,
        backend=persisted.backend,
        queue_name=persisted.queue_name,
    )
    assert acknowledged.record.status is RunStatus.CANCELLED
    assert acknowledged.record.cancellation_reason == persisted.cancellation_reason
    second_engine.dispose()

    third_engine = create_database_engine(database_url)
    third_service = RunService(
        create_default_workflow_registry(),
        repository=SqlAlchemyRunRepository(create_session_factory(third_engine)),
    )
    final_assignment = third_service.get_execution_assignment("run-1")
    assert final_assignment is not None
    assert final_assignment.cancellation_acknowledged_at is not None
    assert third_service.get_run("run-1") == acknowledged.record
    event_types = [
        event.event_type for event in third_service.list_events("run-1", limit=100)
    ]
    assert event_types.count("cancellation_requested") == 1
    assert event_types.count("cancellation_acknowledged") == 1
    third_engine.dispose()


def test_sqlalchemy_cancellation_intent_rolls_back_when_event_insert_fails(
    database_url,
    monkeypatch,
):
    upgrade_database(database_url)
    engine = create_database_engine(database_url)
    setup_service, assignment = _running_service(engine)
    repository = SqlAlchemyRunRepository(create_session_factory(engine))
    service = RunService(
        create_default_workflow_registry(),
        repository=repository,
    )
    events_before = setup_service.list_events("run-1", limit=100)

    def fail_event(_session, _run_id, _draft):
        raise RuntimeError("event insert failed")

    monkeypatch.setattr(repository, "_insert_event", fail_event)
    with pytest.raises(RuntimeError, match="event insert failed"):
        service.request_execution_cancellation(
            "run-1",
            job_id=assignment.job_id,
            backend=assignment.backend,
            queue_name=assignment.queue_name,
            reason="User requested cancellation.",
        )

    persisted = setup_service.get_execution_assignment("run-1")
    assert persisted is not None
    assert persisted.cancellation_requested_at is None
    assert setup_service.get_run("run-1").status is RunStatus.RUNNING
    assert setup_service.list_events("run-1", limit=100) == events_before
    engine.dispose()


def test_sqlalchemy_stop_ack_rolls_back_when_terminal_event_insert_fails(
    database_url,
    monkeypatch,
):
    upgrade_database(database_url)
    engine = create_database_engine(database_url)
    setup_service, assignment = _running_service(engine)
    requested = setup_service.request_execution_cancellation(
        "run-1",
        job_id=assignment.job_id,
        backend=assignment.backend,
        queue_name=assignment.queue_name,
        reason="User requested cancellation.",
    )
    repository = SqlAlchemyRunRepository(create_session_factory(engine))
    service = RunService(
        create_default_workflow_registry(),
        repository=repository,
    )
    events_before = setup_service.list_events("run-1", limit=100)

    def fail_event(_session, _run_id, _draft):
        raise RuntimeError("event insert failed")

    monkeypatch.setattr(repository, "_insert_event", fail_event)
    with pytest.raises(RuntimeError, match="event insert failed"):
        service.acknowledge_execution_stop(
            "run-1",
            job_id=assignment.job_id,
            backend=assignment.backend,
            queue_name=assignment.queue_name,
        )

    persisted = setup_service.get_execution_assignment("run-1")
    assert persisted is not None
    assert persisted.cancellation_requested_at == (
        requested.assignment.cancellation_requested_at
    )
    assert persisted.cancellation_acknowledged_at is None
    assert setup_service.get_run("run-1").status is RunStatus.RUNNING
    assert setup_service.list_events("run-1", limit=100) == events_before
    engine.dispose()


def test_sqlalchemy_concurrent_cancellation_requests_create_one_intent_event(
    database_url,
):
    upgrade_database(database_url)
    setup_engine = create_database_engine(database_url)
    _service, assignment = _running_service(setup_engine)
    setup_engine.dispose()
    first_engine = create_database_engine(database_url)
    second_engine = create_database_engine(database_url)
    first_service = RunService(
        create_default_workflow_registry(),
        repository=SqlAlchemyRunRepository(create_session_factory(first_engine)),
    )
    second_service = RunService(
        create_default_workflow_registry(),
        repository=SqlAlchemyRunRepository(create_session_factory(second_engine)),
    )
    barrier = Barrier(2)

    def request(service, reason):
        barrier.wait(timeout=5)
        return service.request_execution_cancellation(
            "run-1",
            job_id=assignment.job_id,
            backend=assignment.backend,
            queue_name=assignment.queue_name,
            reason=reason,
        )

    try:
        with ThreadPoolExecutor(max_workers=2) as pool:
            futures = (
                pool.submit(request, first_service, "first reason"),
                pool.submit(request, second_service, "second reason"),
            )
            completed = [future.result() for future in futures]

        assert sorted(result.created for result in completed) == [False, True]
        assert completed[0].assignment == completed[1].assignment
        events = first_service.list_events("run-1", limit=100)
        assert [event.event_type for event in events].count(
            "cancellation_requested"
        ) == 1
    finally:
        first_engine.dispose()
        second_engine.dispose()


def test_sqlalchemy_natural_completion_and_stop_ack_are_single_terminal_cas(
    database_url,
):
    upgrade_database(database_url)
    setup_engine = create_database_engine(database_url)
    setup_service, assignment = _running_service(setup_engine)
    setup_service.request_execution_cancellation(
        "run-1",
        job_id=assignment.job_id,
        backend=assignment.backend,
        queue_name=assignment.queue_name,
        reason="User requested cancellation.",
    )
    setup_engine.dispose()
    completion_engine = create_database_engine(database_url)
    acknowledgement_engine = create_database_engine(database_url)
    completion_service = RunService(
        create_default_workflow_registry(),
        repository=SqlAlchemyRunRepository(create_session_factory(completion_engine)),
    )
    acknowledgement_service = RunService(
        create_default_workflow_registry(),
        repository=SqlAlchemyRunRepository(
            create_session_factory(acknowledgement_engine)
        ),
    )
    barrier = Barrier(2)

    def complete():
        barrier.wait(timeout=5)
        try:
            return completion_service.transition_run("run-1", RunStatus.SUCCEEDED)
        except (ConcurrentRunUpdateError, ValueError):
            return completion_service.get_run("run-1")

    def acknowledge():
        barrier.wait(timeout=5)
        return acknowledgement_service.acknowledge_execution_stop(
            "run-1",
            job_id=assignment.job_id,
            backend=assignment.backend,
            queue_name=assignment.queue_name,
        ).record

    try:
        with ThreadPoolExecutor(max_workers=2) as pool:
            futures = (pool.submit(complete), pool.submit(acknowledge))
            [future.result() for future in futures]

        final = completion_service.get_run("run-1")
        assert final.status in {RunStatus.SUCCEEDED, RunStatus.CANCELLED}
        terminal_events = [
            event
            for event in completion_service.list_events("run-1", limit=100)
            if event.status is not None and event.status.is_terminal
        ]
        assert len(terminal_events) == 1
        assert terminal_events[0].status is final.status
    finally:
        completion_engine.dispose()
        acknowledgement_engine.dispose()


def _running_service(engine):
    service = RunService(
        create_default_workflow_registry(),
        id_factory=lambda: "run-1",
        repository=SqlAlchemyRunRepository(create_session_factory(engine)),
    )
    service.create_run(WORKFLOW_ID, WorkflowInputs(config={}))
    service.transition_run("run-1", RunStatus.VALIDATING)
    service.transition_run("run-1", RunStatus.PLANNED)
    assignment = service.ensure_execution_assignment(
        "run-1",
        queue_name="default",
    )
    service.mark_execution_dispatched("run-1", job_id=assignment.job_id)
    service.transition_run("run-1", RunStatus.QUEUED)
    service.claim_execution_assignment(
        "run-1",
        job_id=assignment.job_id,
        backend=assignment.backend,
        queue_name=assignment.queue_name,
    )
    service.transition_run("run-1", RunStatus.RUNNING)
    claimed = service.get_execution_assignment("run-1")
    assert claimed is not None
    return service, claimed


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


def _assignment(run_id: str, job_id: str) -> RunExecutionAssignment:
    return RunExecutionAssignment(
        run_id=run_id,
        job_id=job_id,
        backend="rq",
        queue_name="default",
        created_at=datetime.now(timezone.utc),
    )


def _build_identity() -> WorkflowBuildIdentity:
    return WorkflowBuildIdentity(
        workflow_id=WORKFLOW_ID,
        adapter_version="0.3.0",
        scheme="sha256-tree-v1",
        logical_entrypoint="workflow/Snakefile",
        digest="a" * 64,
        captured_at=datetime.now(timezone.utc),
    )
