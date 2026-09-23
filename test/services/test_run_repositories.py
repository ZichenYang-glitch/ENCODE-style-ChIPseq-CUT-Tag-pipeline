"""Contract tests shared by run repository implementations."""

from __future__ import annotations

from collections import UserDict
from dataclasses import replace
from datetime import datetime, timedelta, timezone
from decimal import Decimal
from hashlib import sha256

import pytest

from encode_pipeline.persistence import (
    SqlAlchemyRunRepository,
    create_database_engine,
    create_session_factory,
    upgrade_database,
)
from encode_pipeline.platform.execution import RunExecutionAssignment
from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.platform.results import Issue
from encode_pipeline.platform.builds import WorkflowBuildIdentity
from encode_pipeline.platform.runs import (
    RunArtifactRef,
    RunQcMetric,
    RunRecord,
    RunStatus,
    build_qc_metric_id,
)
from encode_pipeline.services.run_repositories import (
    ConcurrentRunUpdateError,
    InMemoryRunRepository,
    RunEventDraft,
    RunRepository,
)

from encode_pipeline.services.runs import RunService


@pytest.fixture(params=("memory", "sqlite"))
def repository(request, tmp_path):
    """Only construction and cleanup differ; all behavior assertions are shared."""
    if request.param == "memory":
        yield InMemoryRunRepository()
        return

    database_url = f"sqlite:///{tmp_path / 'runs.db'}"
    upgrade_database(database_url)
    engine = create_database_engine(database_url)
    try:
        yield SqlAlchemyRunRepository(create_session_factory(engine))
    finally:
        try:
            assert engine.pool.checkedout() == 0
        finally:
            engine.dispose()


@pytest.fixture(
    params=(object(), None, [("key", "value")]), ids=("object", "none", "pairs")
)
def invalid_context(request):
    return request.param


def test_repository_create_is_atomic_when_event_is_invalid(repository, invalid_context):
    # Draft construction stays permissive; rejection belongs to the write entry.
    draft = RunEventDraft(
        event_type="status_changed",
        message="Run created.",
        status=RunStatus.CREATED,
        context=invalid_context,
    )
    with pytest.raises(ValueError, match="context must be a mapping"):
        repository.create_run(_record(), draft)

    assert repository.list_runs() == ()
    assert not repository.contains_run("run-1")
    assert repository.get_workflow_build_identity("run-1") is None
    with pytest.raises(KeyError):
        repository.list_events("run-1")


def test_repository_update_is_atomic_when_event_is_invalid(repository, invalid_context):
    record = _record()
    repository.create_run(record, _created_event())
    before_events = repository.list_events(record.run_id)
    draft = RunEventDraft(
        event_type="status_changed",
        message="Run validating.",
        status=RunStatus.VALIDATING,
        context=invalid_context,
    )

    with pytest.raises(ValueError, match="context must be a mapping"):
        repository.update_run(
            replace(record, status=RunStatus.VALIDATING),
            expected_status=RunStatus.CREATED,
            event=draft,
        )

    assert repository.get_run(record.run_id) == record
    assert len(repository.list_events(record.run_id)) == 1
    assert repository.list_events(record.run_id) == before_events
    assert repository.get_workflow_build_identity(record.run_id) is None


def test_repository_replace_artifacts_is_atomic_and_idempotent(repository):
    succeeded = replace(
        _record(),
        status=RunStatus.SUCCEEDED,
        ended_at=datetime.now(timezone.utc),
    )
    repository.create_run(succeeded, _created_event())
    artifacts = (_artifact("run-1", "artifact-1"),)
    draft = RunEventDraft(
        event_type="artifacts_indexed",
        message="Workflow artifacts indexed.",
        status=RunStatus.SUCCEEDED,
        context={"artifact_count": 1},
    )
    attempt_id = "resultattempt-" + "a" * 64
    repository.begin_artifact_result_attempt(
        "run-1",
        attempt_id=attempt_id,
        expected_status=RunStatus.SUCCEEDED,
    )

    first = repository.replace_artifacts(
        "run-1",
        artifacts,
        attempt_id=attempt_id,
        expected_status=RunStatus.SUCCEEDED,
        event=draft,
    )
    second = repository.replace_artifacts(
        "run-1",
        artifacts,
        attempt_id=attempt_id,
        expected_status=RunStatus.SUCCEEDED,
        event=draft,
    )

    assert first is not None
    assert second is None
    assert repository.list_artifacts("run-1") == artifacts
    assert [event.event_type for event in repository.list_events("run-1")].count(
        "artifacts_indexed"
    ) == 1


def test_repository_artifact_replacement_rolls_back_when_event_is_invalid(repository):
    succeeded = replace(
        _record(),
        status=RunStatus.SUCCEEDED,
        ended_at=datetime.now(timezone.utc),
    )
    repository.create_run(succeeded, _created_event())
    original = (_artifact("run-1", "artifact-original"),)
    _replace_artifacts(repository, "run-1", original)
    artifact_generation = repository.get_result_state("run-1").artifact_generation
    assert artifact_generation is not None
    qc_attempt_id = "resultattempt-" + "b" * 64
    repository.begin_qc_result_attempt(
        "run-1",
        attempt_id=qc_attempt_id,
        expected_artifact_generation=artifact_generation,
        expected_artifacts=original,
        expected_status=RunStatus.SUCCEEDED,
    )
    repository.replace_qc_metrics(
        "run-1",
        (),
        attempt_id=qc_attempt_id,
        expected_artifact_generation=artifact_generation,
        expected_artifacts=original,
        expected_status=RunStatus.SUCCEEDED,
        event=RunEventDraft(
            event_type="qc_metrics_indexed",
            message="Workflow QC metrics indexed.",
            status=RunStatus.SUCCEEDED,
        ),
    )
    artifact_attempt_id = "resultattempt-" + "c" * 64
    repository.begin_artifact_result_attempt(
        "run-1",
        attempt_id=artifact_attempt_id,
        expected_status=RunStatus.SUCCEEDED,
    )
    before_state = repository.get_result_state("run-1")
    before_events = repository.list_events("run-1")

    with pytest.raises(TypeError):
        repository.replace_artifacts(
            "run-1",
            (_artifact("run-1", "artifact-replacement"),),
            attempt_id=artifact_attempt_id,
            expected_status=RunStatus.SUCCEEDED,
            event=RunEventDraft(
                event_type="artifacts_indexed",
                message="Workflow artifacts indexed.",
                status=RunStatus.SUCCEEDED,
                context=object(),
            ),
        )

    assert repository.list_artifacts("run-1") == original
    assert repository.list_qc_metrics("run-1") == ()
    assert repository.get_result_state("run-1") == before_state
    assert repository.list_events("run-1") == before_events


def test_repository_qc_replacement_rolls_back_when_event_is_invalid(repository):
    succeeded = replace(
        _record(),
        status=RunStatus.SUCCEEDED,
        ended_at=datetime.now(timezone.utc),
    )
    repository.create_run(succeeded, _created_event())
    artifacts = (_artifact("run-1", "artifact-1"),)
    _replace_artifacts(repository, "run-1", artifacts)
    artifact_generation = repository.get_result_state("run-1").artifact_generation
    assert artifact_generation is not None
    attempt_id = "resultattempt-" + "d" * 64
    repository.begin_qc_result_attempt(
        "run-1",
        attempt_id=attempt_id,
        expected_artifact_generation=artifact_generation,
        expected_artifacts=artifacts,
        expected_status=RunStatus.SUCCEEDED,
    )
    before_state = repository.get_result_state("run-1")
    before_events = repository.list_events("run-1")
    metric_key = "mapping.rate"
    metric = RunQcMetric(
        metric_id=build_qc_metric_id(metric_key, "run", None, None),
        run_id="run-1",
        metric_key=metric_key,
        display_name="Mapping rate",
        value=Decimal("0.5"),
        unit="fraction",
        scope="run",
        source_artifact_id=artifacts[0].artifact_id,
        produced_at=datetime.now(timezone.utc),
    )

    with pytest.raises(TypeError):
        repository.replace_qc_metrics(
            "run-1",
            (metric,),
            attempt_id=attempt_id,
            expected_artifact_generation=artifact_generation,
            expected_artifacts=artifacts,
            expected_status=RunStatus.SUCCEEDED,
            event=RunEventDraft(
                event_type="qc_metrics_indexed",
                message="Workflow QC metrics indexed.",
                status=RunStatus.SUCCEEDED,
                context=object(),
            ),
        )

    assert repository.list_artifacts("run-1") == artifacts
    assert repository.list_qc_metrics("run-1") == ()
    assert repository.get_result_state("run-1") == before_state
    assert repository.list_events("run-1") == before_events


def test_repository_artifact_failure_rolls_back_when_event_is_invalid(repository):
    succeeded = replace(
        _record(),
        status=RunStatus.SUCCEEDED,
        ended_at=datetime.now(timezone.utc),
    )
    repository.create_run(succeeded, _created_event())
    artifacts = (_artifact("run-1", "artifact-1"),)
    _replace_artifacts(repository, "run-1", artifacts)
    artifact_generation = repository.get_result_state("run-1").artifact_generation
    assert artifact_generation is not None
    qc_attempt_id = "resultattempt-" + "e" * 64
    repository.begin_qc_result_attempt(
        "run-1",
        attempt_id=qc_attempt_id,
        expected_artifact_generation=artifact_generation,
        expected_artifacts=artifacts,
        expected_status=RunStatus.SUCCEEDED,
    )
    metric_key = "mapping.rate"
    metric = RunQcMetric(
        metric_id=build_qc_metric_id(metric_key, "run", None, None),
        run_id="run-1",
        metric_key=metric_key,
        display_name="Mapping rate",
        value=Decimal("0.5"),
        unit="fraction",
        scope="run",
        source_artifact_id=artifacts[0].artifact_id,
        produced_at=datetime.now(timezone.utc),
    )
    repository.replace_qc_metrics(
        "run-1",
        (metric,),
        attempt_id=qc_attempt_id,
        expected_artifact_generation=artifact_generation,
        expected_artifacts=artifacts,
        expected_status=RunStatus.SUCCEEDED,
        event=RunEventDraft(
            event_type="qc_metrics_indexed",
            message="Workflow QC metrics indexed.",
            status=RunStatus.SUCCEEDED,
        ),
    )
    artifact_attempt_id = "resultattempt-" + "f" * 64
    repository.begin_artifact_result_attempt(
        "run-1",
        attempt_id=artifact_attempt_id,
        expected_status=RunStatus.SUCCEEDED,
    )
    before_artifacts = repository.list_artifacts("run-1")
    before_metrics = repository.list_qc_metrics("run-1")
    before_state = repository.get_result_state("run-1")
    before_events = repository.list_events("run-1")

    with pytest.raises(TypeError):
        repository.record_artifact_failure(
            "run-1",
            attempt_id=artifact_attempt_id,
            reason_code="ARTIFACT_EXTRACTION_FAILED",
            expected_status=RunStatus.SUCCEEDED,
            event=RunEventDraft(
                event_type="artifact_extraction_failed",
                message="Workflow artifacts could not be indexed.",
                status=RunStatus.SUCCEEDED,
                context=object(),
            ),
        )

    assert before_metrics
    assert repository.list_artifacts("run-1") == before_artifacts
    assert repository.list_qc_metrics("run-1") == before_metrics
    assert repository.get_result_state("run-1") == before_state
    assert repository.list_events("run-1") == before_events


def test_repository_qc_failure_rolls_back_when_event_is_invalid(repository):
    succeeded = replace(
        _record(),
        status=RunStatus.SUCCEEDED,
        ended_at=datetime.now(timezone.utc),
    )
    repository.create_run(succeeded, _created_event())
    artifacts = (_artifact("run-1", "artifact-1"),)
    _replace_artifacts(repository, "run-1", artifacts)
    artifact_generation = repository.get_result_state("run-1").artifact_generation
    assert artifact_generation is not None
    first_attempt_id = "resultattempt-" + "1" * 64
    repository.begin_qc_result_attempt(
        "run-1",
        attempt_id=first_attempt_id,
        expected_artifact_generation=artifact_generation,
        expected_artifacts=artifacts,
        expected_status=RunStatus.SUCCEEDED,
    )
    metric_key = "mapping.rate"
    metric = RunQcMetric(
        metric_id=build_qc_metric_id(metric_key, "run", None, None),
        run_id="run-1",
        metric_key=metric_key,
        display_name="Mapping rate",
        value=Decimal("0.5"),
        unit="fraction",
        scope="run",
        source_artifact_id=artifacts[0].artifact_id,
        produced_at=datetime.now(timezone.utc),
    )
    repository.replace_qc_metrics(
        "run-1",
        (metric,),
        attempt_id=first_attempt_id,
        expected_artifact_generation=artifact_generation,
        expected_artifacts=artifacts,
        expected_status=RunStatus.SUCCEEDED,
        event=RunEventDraft(
            event_type="qc_metrics_indexed",
            message="Workflow QC metrics indexed.",
            status=RunStatus.SUCCEEDED,
        ),
    )
    failure_attempt_id = "resultattempt-" + "2" * 64
    repository.begin_qc_result_attempt(
        "run-1",
        attempt_id=failure_attempt_id,
        expected_artifact_generation=artifact_generation,
        expected_artifacts=artifacts,
        expected_status=RunStatus.SUCCEEDED,
    )
    before_artifacts = repository.list_artifacts("run-1")
    before_metrics = repository.list_qc_metrics("run-1")
    before_state = repository.get_result_state("run-1")
    before_events = repository.list_events("run-1")

    with pytest.raises(TypeError):
        repository.record_qc_metrics_failure(
            "run-1",
            attempt_id=failure_attempt_id,
            expected_artifact_generation=artifact_generation,
            reason_code="QC_INDEXING_ADAPTER_FAILED",
            expected_status=RunStatus.SUCCEEDED,
            event=RunEventDraft(
                event_type="qc_metrics_indexing_failed",
                message="Workflow QC metrics could not be indexed.",
                status=RunStatus.SUCCEEDED,
                context=object(),
            ),
        )

    assert before_metrics
    assert repository.list_artifacts("run-1") == before_artifacts
    assert repository.list_qc_metrics("run-1") == before_metrics
    assert repository.get_result_state("run-1") == before_state
    assert repository.list_events("run-1") == before_events


def test_repository_replace_artifacts_rejects_non_succeeded_without_mutation(
    repository,
):
    repository.create_run(_record(), _created_event())

    with pytest.raises(ConcurrentRunUpdateError):
        repository.replace_artifacts(
            "run-1",
            (_artifact("run-1", "artifact-1"),),
            attempt_id="resultattempt-" + "a" * 64,
            expected_status=RunStatus.SUCCEEDED,
            event=RunEventDraft(
                event_type="artifacts_indexed",
                message="Workflow artifacts indexed.",
            ),
        )

    assert repository.list_artifacts("run-1") == ()


def test_repository_artifact_queries_are_sorted_paginated_and_run_scoped(repository):
    succeeded = replace(
        _record(),
        status=RunStatus.SUCCEEDED,
        ended_at=datetime.now(timezone.utc),
    )
    repository.create_run(succeeded, _created_event())
    repository.create_run(replace(succeeded, run_id="run-2"), _created_event())
    artifact_z = _artifact("run-1", "artifact-z")
    artifact_a = _artifact("run-1", "artifact-a")
    artifact_m = _artifact("run-1", "artifact-m")
    other = _artifact("run-2", "artifact-other")
    _replace_artifacts(repository, "run-1", (artifact_z, artifact_a, artifact_m))
    _replace_artifacts(repository, "run-2", (other,))

    assert repository.list_artifacts("run-1") == (
        artifact_a,
        artifact_m,
        artifact_z,
    )
    assert repository.list_artifacts("run-1", limit=2) == (
        artifact_a,
        artifact_m,
    )
    assert repository.list_artifacts("run-1", after="artifact-m", limit=2) == (
        artifact_z,
    )
    assert repository.get_artifact("run-1", "artifact-m") == artifact_m

    with pytest.raises(KeyError):
        repository.list_artifacts("run-1", after=other.artifact_id)
    with pytest.raises(KeyError):
        repository.get_artifact("run-1", other.artifact_id)


def test_repository_complete_preflight_atomically_binds_build_identity(repository):
    created = _record()
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

    event = repository.complete_preflight(
        planned,
        identity,
        expected_status=RunStatus.VALIDATING,
        event=RunEventDraft(
            event_type="preflight_completed",
            message="Preflight complete.",
            status=RunStatus.PLANNED,
        ),
    )

    assert event.status is RunStatus.PLANNED
    assert repository.get_run(created.run_id) == planned
    assert repository.get_workflow_build_identity(created.run_id) == identity


def test_repository_complete_preflight_rolls_back_invalid_event(
    repository, invalid_context
):
    validating = replace(_record(), status=RunStatus.VALIDATING)
    repository.create_run(validating, _created_event())
    before_events = repository.list_events(validating.run_id)
    draft = RunEventDraft(
        event_type="preflight_completed",
        message="Preflight complete.",
        status=RunStatus.PLANNED,
        context=invalid_context,
    )

    with pytest.raises(ValueError, match="context must be a mapping"):
        repository.complete_preflight(
            replace(validating, status=RunStatus.PLANNED),
            _build_identity(),
            expected_status=RunStatus.VALIDATING,
            event=draft,
        )

    assert repository.get_run(validating.run_id) == validating
    assert repository.get_workflow_build_identity(validating.run_id) is None
    assert repository.list_events(validating.run_id) == before_events


def test_repository_execution_assignment_is_idempotent_per_run(repository):
    record = _record()
    repository.create_run(record, _created_event())
    original = _assignment(record.run_id, "job-original")
    replacement = _assignment(record.run_id, "job-replacement")

    assert repository.get_execution_assignment(record.run_id) is None
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
    assert repository.get_execution_assignment(record.run_id) == original


def test_repository_dispatch_mark_is_idempotent_after_status_changes(repository):
    record = _record()
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
        dispatched_at=dispatched_at,
        allowed_statuses=frozenset({RunStatus.CREATED}),
    )

    assert retried == dispatched
    assert retried.dispatched_at == dispatched_at


def test_repository_queue_dispatched_run_is_atomic_and_idempotent(repository):
    planned = replace(_record(), status=RunStatus.PLANNED)
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


def test_repository_execution_assignment_rejects_cross_run_job_reuse(repository):
    first = _record()
    second = replace(first, run_id="run-2")
    repository.create_run(first, _created_event())
    repository.create_run(second, _created_event())
    repository.ensure_execution_assignment(
        _assignment(first.run_id, "shared-job"),
        expected_status=RunStatus.CREATED,
    )

    with pytest.raises(ValueError, match="shared-job.*already assigned"):
        repository.ensure_execution_assignment(
            _assignment(second.run_id, "shared-job"),
            expected_status=RunStatus.CREATED,
        )


def test_repository_execution_assignment_requires_a_persisted_run(repository):

    with pytest.raises(KeyError, match="missing"):
        repository.ensure_execution_assignment(
            _assignment("missing", "job-1"),
            expected_status=RunStatus.CREATED,
        )
    assert repository.get_execution_assignment("missing") is None


def test_repository_creation_round_trips_and_duplicate_id_preserves_state(repository):
    first = replace(
        _record(),
        inputs={"config": {"enabled": True}, "samples": [], "options": {}},
        tags={"purpose": "contract"},
    )
    second = replace(first, run_id="run-2", workflow_id="other-workflow")
    assert repository.list_runs() == ()
    assert not repository.contains_run(first.run_id)

    first_event = repository.create_run(first, _created_event())
    second_event = repository.create_run(second, _created_event())

    assert repository.contains_run(first.run_id)
    assert repository.get_run(first.run_id) == first
    assert repository.get_run(second.run_id) == second
    assert repository.list_runs() == (first, second)
    assert repository.get_run_requester_user_id(first.run_id) is None
    assert repository.list_events(first.run_id) == (first_event,)
    assert repository.list_events(second.run_id) == (second_event,)
    assert first_event.run_id == first.run_id
    assert first_event.status is RunStatus.CREATED
    assert first_event.sequence == second_event.sequence == 1

    with pytest.raises(ValueError, match="Duplicate run_id"):
        repository.create_run(
            replace(first, tags={"purpose": "replacement"}), _created_event()
        )

    assert repository.list_runs() == (first, second)
    assert repository.list_events(first.run_id) == (first_event,)
    assert repository.list_events(second.run_id) == (second_event,)


def test_repository_accepts_non_dict_mapping_at_event_write_entries(repository):
    context = UserDict({"sample": "S1", "counts": [1, 2]})
    created = _record()
    draft = RunEventDraft(
        "status_changed", "Created.", status=created.status, context=context
    )
    first = repository.create_run(created, draft)
    validating = replace(created, status=RunStatus.VALIDATING)
    second = repository.update_run(
        validating,
        expected_status=RunStatus.CREATED,
        event=replace(draft, status=RunStatus.VALIDATING),
    )
    planned = replace(validating, status=RunStatus.PLANNED)
    identity = _build_identity()
    third = repository.complete_preflight(
        planned,
        identity,
        expected_status=RunStatus.VALIDATING,
        event=replace(draft, status=RunStatus.PLANNED),
    )
    fourth = repository.add_event(
        created.run_id, replace(draft, event_type="note", status=None)
    )

    events = (first, second, third, fourth)
    assert repository.list_events(created.run_id) == events
    assert [event.sequence for event in events] == [1, 2, 3, 4]
    assert all(event.context == context for event in events)
    assert repository.get_run(created.run_id) == planned
    assert repository.get_workflow_build_identity(created.run_id) == identity


def test_repository_update_returns_persisted_event_and_rejects_stale_status(repository):
    record = _record()
    created_event = repository.create_run(record, _created_event())
    updated = replace(
        record,
        status=RunStatus.VALIDATING,
        updated_at=record.updated_at + timedelta(seconds=1),
        current_stage="validation",
    )
    draft = RunEventDraft(
        "status_changed",
        "Validating.",
        status=updated.status,
        stage="validation",
        context={"previous_status": "created", "new_status": "validating"},
    )
    event = repository.update_run(updated, expected_status=record.status, event=draft)
    assert repository.get_run(record.run_id) == updated
    assert repository.list_events(record.run_id) == (created_event, event)
    assert (event.sequence, event.status, event.stage, event.context) == (
        2,
        RunStatus.VALIDATING,
        "validation",
        draft.context,
    )

    # A delayed writer must not overwrite the winner or append a second event.
    with pytest.raises(ConcurrentRunUpdateError):
        repository.update_run(
            replace(updated, status=RunStatus.PLANNED),
            expected_status=RunStatus.CREATED,
            event=draft,
        )
    assert repository.get_run(record.run_id) == updated
    assert repository.list_events(record.run_id) == (created_event, event)


def test_repository_missing_run_operations_raise_key_error_without_creating_state(
    repository,
):
    missing = replace(_record(), run_id="missing")
    operations = (
        lambda: repository.get_run(missing.run_id),
        lambda: repository.get_run_requester_user_id(missing.run_id),
        lambda: repository.list_events(missing.run_id),
        lambda: repository.list_logs(missing.run_id),
        lambda: repository.add_event(missing.run_id, _created_event()),
        lambda: repository.append_log(missing.run_id, "stdout", ["line"]),
        lambda: repository.update_run(
            missing, expected_status=missing.status, event=_created_event()
        ),
        lambda: repository.complete_preflight(
            replace(missing, status=RunStatus.PLANNED),
            _build_identity(),
            expected_status=RunStatus.VALIDATING,
            event=_created_event(),
        ),
    )
    for operation in operations:
        with pytest.raises(KeyError) as raised:
            operation()
        assert raised.value.args == (missing.run_id,)
    assert not repository.contains_run(missing.run_id)
    assert repository.list_runs() == ()
    assert repository.get_workflow_build_identity(missing.run_id) is None


def test_repository_events_are_append_only_ordered_and_paginated_per_run(repository):
    record = _record()
    first = repository.create_run(record, _created_event())
    other = repository.create_run(replace(record, run_id="run-2"), _created_event())
    issue = Issue(
        code="TEST_WARNING",
        message="Synthetic warning.",
        severity="warning",
        source="test",
    )
    draft = RunEventDraft(
        "note",
        "Repeated note.",
        stage="validation",
        context={"attempt": 1},
        issue=issue,
    )
    second = repository.add_event(record.run_id, draft)
    third = repository.add_event(record.run_id, draft)
    events = (first, second, third)

    assert repository.list_events(record.run_id) == events
    assert [event.sequence for event in events] == [1, 2, 3]
    assert len({event.event_id for event in events}) == 3
    assert second.issue == third.issue == issue
    assert second.context == third.context == draft.context
    assert repository.get_run(record.run_id) == record
    assert repository.list_events(record.run_id, limit=1) == (first,)
    assert repository.list_events(record.run_id, after=first.event_id, limit=1) == (
        second,
    )
    assert repository.list_events(record.run_id, after=second.event_id, limit=10) == (
        third,
    )
    assert repository.list_events(record.run_id, after=third.event_id, limit=1) == ()
    assert repository.list_events("run-2") == (other,)
    for run_id, cursor in ((record.run_id, "evt-missing"), ("run-2", third.event_id)):
        with pytest.raises(KeyError):
            repository.list_events(run_id, after=cursor)
    assert repository.list_events(record.run_id) == events


def test_repository_logs_are_append_only_ordered_and_paginated_per_stream(repository):
    record = _record()
    repository.create_run(record, _created_event())
    repository.create_run(replace(record, run_id="run-2"), _created_event())
    before_events = repository.list_events(record.run_id)
    assert repository.list_logs(record.run_id) == ()
    lines = ("first line\n", "second line")
    first = repository.append_log(record.run_id, "stdout", iter(lines))
    error = repository.append_log(record.run_id, "stderr", ["warning"])
    second = repository.append_log(record.run_id, "stdout", iter(lines))
    other = repository.append_log("run-2", "stdout", ["other run"])

    assert first.lines == second.lines == lines
    assert (first.sequence, second.sequence, error.sequence, other.sequence) == (
        1,
        2,
        1,
        1,
    )
    assert first.chunk_id != second.chunk_id
    assert (first.run_id, first.stream_name) == (record.run_id, "stdout")
    assert repository.list_logs(record.run_id) == (first, second)
    assert repository.list_logs(record.run_id, limit=1) == (first,)
    assert repository.list_logs(record.run_id, after=first.chunk_id, limit=1) == (
        second,
    )
    assert repository.list_logs(record.run_id, after=second.chunk_id, limit=10) == ()
    assert repository.list_logs(record.run_id, "stderr") == (error,)
    assert repository.list_logs("run-2") == (other,)
    for run_id, stream, cursor in (
        (record.run_id, "stdout", "log-missing"),
        (record.run_id, "stderr", second.chunk_id),
        ("run-2", "stdout", second.chunk_id),
    ):
        with pytest.raises(KeyError):
            repository.list_logs(run_id, stream, after=cursor)
    assert repository.get_run(record.run_id) == record
    assert repository.list_events(record.run_id) == before_events


@pytest.mark.parametrize(
    "terminal", (RunStatus.SUCCEEDED, RunStatus.FAILED, RunStatus.CANCELLED)
)
def test_service_lifecycle_and_terminal_retry_preserve_repository_state(
    repository, terminal
):
    # Graph validation belongs to RunService, not to repository.update_run.
    record = _record()
    repository.create_run(record, _created_event())
    service = RunService(WorkflowRegistry(), repository=repository)
    before_events = repository.list_events(record.run_id)
    with pytest.raises(ValueError, match="Illegal transition"):
        service.transition_run(record.run_id, RunStatus.RUNNING)
    assert repository.get_run(record.run_id) == record
    assert repository.list_events(record.run_id) == before_events

    path = {
        RunStatus.SUCCEEDED: (
            RunStatus.VALIDATING,
            RunStatus.PLANNED,
            RunStatus.QUEUED,
            RunStatus.RUNNING,
            RunStatus.SUCCEEDED,
        ),
        RunStatus.FAILED: (RunStatus.VALIDATING, RunStatus.FAILED),
        RunStatus.CANCELLED: (RunStatus.CANCELLED,),
    }[terminal]
    for status in path:
        updated = service.transition_run(record.run_id, status)
        assert repository.get_run(record.run_id) == updated
    assert updated.ended_at is not None
    events = repository.list_events(record.run_id)
    assert tuple(event.status for event in events) == (RunStatus.CREATED, *path)
    assert [event.sequence for event in events] == list(range(1, len(events) + 1))
    assert service.cancel_run(record.run_id, reason="Late retry") == updated
    with pytest.raises(ValueError, match="Illegal transition"):
        service.transition_run(record.run_id, RunStatus.VALIDATING)
    assert repository.get_run(record.run_id) == updated
    assert repository.list_events(record.run_id) == events


def _record() -> RunRecord:
    now = datetime.now(timezone.utc)
    return RunRecord(
        run_id="run-1",
        workflow_id="fake",
        inputs={"config": {}, "samples": None, "options": {}},
        status=RunStatus.CREATED,
        created_at=now,
        updated_at=now,
        started_at=None,
        ended_at=None,
        current_stage=None,
        cancellation_reason=None,
        error=None,
        tags={},
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
        name=f"{artifact_id}.txt",
        uri=f"run://runs/{run_id}/artifacts/{artifact_id}",
        mime_type="text/plain",
        produced_at=datetime.now(timezone.utc),
        revision="artifactrev-" + sha256(artifact_id.encode()).hexdigest(),
        metadata={},
    )


def _replace_artifacts(
    repository: RunRepository,
    run_id: str,
    artifacts: tuple[RunArtifactRef, ...],
) -> None:
    attempt_id = "resultattempt-" + sha256(f"artifacts:{run_id}".encode()).hexdigest()
    repository.begin_artifact_result_attempt(
        run_id,
        attempt_id=attempt_id,
        expected_status=RunStatus.SUCCEEDED,
    )
    repository.replace_artifacts(
        run_id,
        artifacts,
        attempt_id=attempt_id,
        expected_status=RunStatus.SUCCEEDED,
        event=RunEventDraft(
            event_type="artifacts_indexed",
            message="Workflow artifacts indexed.",
            status=RunStatus.SUCCEEDED,
        ),
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
        workflow_id="fake",
        adapter_version="1.0.0",
        scheme="sha256-tree-v1",
        logical_entrypoint="workflow/Snakefile",
        digest="a" * 64,
        captured_at=datetime.now(timezone.utc),
    )
