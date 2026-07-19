"""Contract tests shared by run repository implementations."""

from __future__ import annotations

from dataclasses import replace
from datetime import datetime, timezone
from decimal import Decimal
from hashlib import sha256

import pytest

from encode_pipeline.platform.execution import RunExecutionAssignment
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
)


def test_in_memory_create_is_atomic_when_event_is_invalid():
    repository = InMemoryRunRepository()

    with pytest.raises(ValueError, match="context must be a mapping"):
        repository.create_run(
            _record(),
            RunEventDraft(
                event_type="status_changed",
                message="Run created.",
                status=RunStatus.CREATED,
                context=object(),
            ),
        )

    assert repository.list_runs() == ()


def test_in_memory_update_is_atomic_when_event_is_invalid():
    repository = InMemoryRunRepository()
    record = _record()
    repository.create_run(record, _created_event())

    with pytest.raises(ValueError, match="context must be a mapping"):
        repository.update_run(
            replace(record, status=RunStatus.VALIDATING),
            expected_status=RunStatus.CREATED,
            event=RunEventDraft(
                event_type="status_changed",
                message="Run validating.",
                status=RunStatus.VALIDATING,
                context=object(),
            ),
        )

    assert repository.get_run(record.run_id) == record
    assert len(repository.list_events(record.run_id)) == 1


def test_in_memory_replace_artifacts_is_atomic_and_idempotent():
    repository = InMemoryRunRepository()
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


def test_in_memory_artifact_replacement_rolls_back_when_event_is_invalid():
    repository = InMemoryRunRepository()
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


def test_in_memory_qc_replacement_rolls_back_when_event_is_invalid():
    repository = InMemoryRunRepository()
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


def test_in_memory_artifact_failure_rolls_back_when_event_is_invalid():
    repository = InMemoryRunRepository()
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


def test_in_memory_qc_failure_rolls_back_when_event_is_invalid():
    repository = InMemoryRunRepository()
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


def test_in_memory_replace_artifacts_rejects_non_succeeded_without_mutation():
    repository = InMemoryRunRepository()
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


def test_in_memory_artifact_queries_are_sorted_paginated_and_run_scoped():
    repository = InMemoryRunRepository()
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


def test_in_memory_complete_preflight_atomically_binds_build_identity():
    repository = InMemoryRunRepository()
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


def test_in_memory_complete_preflight_rolls_back_invalid_event():
    repository = InMemoryRunRepository()
    validating = replace(_record(), status=RunStatus.VALIDATING)
    repository.create_run(validating, _created_event())

    with pytest.raises(ValueError, match="context must be a mapping"):
        repository.complete_preflight(
            replace(validating, status=RunStatus.PLANNED),
            _build_identity(),
            expected_status=RunStatus.VALIDATING,
            event=RunEventDraft(
                event_type="preflight_completed",
                message="Preflight complete.",
                status=RunStatus.PLANNED,
                context=object(),
            ),
        )

    assert repository.get_run(validating.run_id) == validating
    assert repository.get_workflow_build_identity(validating.run_id) is None


def test_in_memory_execution_assignment_is_idempotent_per_run():
    repository = InMemoryRunRepository()
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


def test_in_memory_dispatch_mark_is_idempotent_after_status_changes():
    repository = InMemoryRunRepository()
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


def test_in_memory_queue_dispatched_run_is_atomic_and_idempotent():
    repository = InMemoryRunRepository()
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


def test_in_memory_execution_assignment_rejects_cross_run_job_reuse():
    repository = InMemoryRunRepository()
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


def test_in_memory_execution_assignment_requires_a_persisted_run():
    repository = InMemoryRunRepository()

    with pytest.raises(KeyError, match="missing"):
        repository.ensure_execution_assignment(
            _assignment("missing", "job-1"),
            expected_status=RunStatus.CREATED,
        )
    assert repository.get_execution_assignment("missing") is None


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
    repository: InMemoryRunRepository,
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
