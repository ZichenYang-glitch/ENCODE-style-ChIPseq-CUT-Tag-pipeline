"""Internal atomic artifact/QC publication against both real repositories."""

from dataclasses import replace
from datetime import datetime, timezone
from decimal import Decimal
from hashlib import sha256

import pytest
from sqlalchemy import func, select

from encode_pipeline.persistence import (
    SqlAlchemyRunRepository,
    create_database_engine,
    create_session_factory,
    upgrade_database,
)
from encode_pipeline.persistence.models import ArtifactPublicationRow
from encode_pipeline.platform.registry import WorkflowRegistry
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
    result_bundle_is_complete,
)
from encode_pipeline.services.runs import RunService

NOW = datetime(2026, 9, 24, tzinfo=timezone.utc)


@pytest.fixture(params=("memory", "sqlite"))
def store(request, tmp_path):
    if request.param == "memory":
        yield InMemoryRunRepository()
        return
    url = f"sqlite:///{tmp_path / 'bundle.sqlite'}"
    upgrade_database(url)
    engine = create_database_engine(url)
    try:
        yield SqlAlchemyRunRepository(create_session_factory(engine))
    finally:
        assert engine.pool.checkedout() == 0
        engine.dispose()


def service_for(store):
    store.create_run(
        RunRecord(
            run_id="run-1",
            workflow_id="bundle",
            inputs={"config": {}, "samples": None, "options": {}},
            status=RunStatus.SUCCEEDED,
            created_at=NOW,
            updated_at=NOW,
            ended_at=NOW,
            started_at=NOW,
            current_stage=None,
            cancellation_reason=None,
            error=None,
        ),
        RunEventDraft(
            event_type="status_changed",
            message="Succeeded.",
            status=RunStatus.SUCCEEDED,
        ),
    )
    return RunService(repository=store, registry=WorkflowRegistry())


def candidates(version="one"):
    artifact = RunArtifactRef(
        artifact_id="summary",
        run_id="run-1",
        artifact_type="file",
        name="summary.tsv",
        uri="run://runs/run-1/artifacts/summary",
        produced_at=NOW,
        revision="artifactrev-" + sha256(version.encode()).hexdigest(),
        mime_type="text/tab-separated-values",
        metadata={"output_type": "summary"},
    )
    metric = RunQcMetric(
        metric_id=build_qc_metric_id("pets", "sample", "sample-a", None),
        run_id="run-1",
        metric_key="pets",
        display_name="PETs",
        value=Decimal(8),
        unit="count",
        scope="sample",
        sample_id="sample-a",
        source_artifact_id="summary",
        produced_at=NOW,
    )
    return (artifact,), (metric,)


def publications(store):
    if isinstance(store, InMemoryRunRepository):
        return len(store._artifact_publications)
    with store._session_factory() as session:
        return session.scalar(select(func.count()).select_from(ArtifactPublicationRow))


def publish(service, artifacts, metrics, *, attempt=None, generation=None):
    if attempt is None:
        state = service.begin_artifact_result_attempt("run-1")
        attempt, generation = state.artifact_attempt_id, state.artifact_generation
    return service.publish_result_bundle(
        "run-1",
        artifacts,
        metrics,
        attempt_id=attempt,
        expected_artifact_generation=generation,
    )


def public_snapshot(store):
    state = store.get_result_state("run-1")
    return (
        store.list_artifacts("run-1"),
        store.list_qc_metrics("run-1"),
        state.artifact_generation,
        state.qc_generation,
        publications(store),
    )


def test_bundle_one_transaction_generations_idempotence_and_old_attempt(store):
    service = service_for(store)
    artifacts, metrics = candidates()
    state = publish(service, artifacts, metrics)
    assert result_bundle_is_complete(state, state.artifact_attempt_id)
    assert state.qc_artifact_generation == state.artifact_generation
    assert service.list_artifacts("run-1") == artifacts
    assert service.list_qc_metrics("run-1") == metrics
    assert publications(store) == 1
    events = store.list_events("run-1")
    assert [e.event_type for e in events][-2:] == [
        "artifacts_indexed",
        "qc_metrics_indexed",
    ]
    assert (
        publish(service, artifacts, metrics, attempt=state.artifact_attempt_id) == state
    )
    assert store.list_events("run-1") == events
    assert publications(store) == 1
    with pytest.raises(ConcurrentRunUpdateError):
        publish(service, *candidates("other"), attempt=state.artifact_attempt_id)
    new = publish(service, artifacts, metrics)
    assert new.artifact_generation == state.artifact_generation
    assert new.qc_generation == state.qc_generation
    assert store.list_events("run-1") == events
    assert publications(store) == 1
    with pytest.raises(ConcurrentRunUpdateError):
        publish(service, artifacts, metrics, attempt=state.artifact_attempt_id)


def test_bundle_invalid_qc_preserves_old_complete_results(store):
    service = service_for(store)
    publish(service, *candidates())
    before = public_snapshot(store)
    artifacts, metrics = candidates("next")
    with pytest.raises(ValueError):
        publish(service, artifacts, (replace(metrics[0], source_artifact_id="absent"),))
    assert public_snapshot(store) == before


def test_bundle_stale_generation_rejected_without_partial_rows(store):
    service = service_for(store)
    first = publish(service, *candidates())
    pending = service.begin_artifact_result_attempt("run-1")
    before = public_snapshot(store)
    with pytest.raises(ConcurrentRunUpdateError):
        publish(
            service,
            *candidates("next"),
            attempt=pending.artifact_attempt_id,
            generation="artifactgen-" + "f" * 64,
        )
    assert public_snapshot(store) == before
    assert not result_bundle_is_complete(pending, first.artifact_attempt_id)


def test_bundle_late_event_error_rolls_back_every_public_write(store, monkeypatch):
    service = service_for(store)
    publish(service, *candidates())
    before = public_snapshot(store)
    events = store.list_events("run-1")
    method = (
        "_make_event" if isinstance(store, InMemoryRunRepository) else "_insert_event"
    )
    original = getattr(store, method)

    def fail(*args, **kwargs):
        event = original(*args, **kwargs)
        if args[-1].event_type == "qc_metrics_indexed":
            raise RuntimeError("injected after QC event construction/SQL flush")
        return event

    monkeypatch.setattr(store, method, fail)
    with pytest.raises(RuntimeError, match="injected"):
        publish(service, *candidates("next"))
    assert public_snapshot(store) == before
    assert store.list_events("run-1") == events


def test_sql_reader_sees_old_complete_bundle_until_commit(tmp_path, monkeypatch):
    url = f"sqlite:///{tmp_path / 'readers.sqlite'}"
    upgrade_database(url)
    engine = create_database_engine(url)
    try:
        factory = create_session_factory(engine)
        writer, reader = (
            SqlAlchemyRunRepository(factory),
            SqlAlchemyRunRepository(factory),
        )
        service = service_for(writer)
        publish(service, *candidates())
        before = public_snapshot(reader)
        original = writer._insert_event
        observations = []

        def observe(*args, **kwargs):
            result = original(*args, **kwargs)
            observations.append(public_snapshot(reader))
            return result

        monkeypatch.setattr(writer, "_insert_event", observe)
        publish(service, *candidates("next"))
        assert observations and all(value == before for value in observations)
        after = public_snapshot(reader)
        assert after != before
        assert after[2] != before[2] and after[3] != before[3]
        assert after[4] == 2
        assert reader.get_result_state("run-1").qc_artifact_generation == after[2]
    finally:
        engine.dispose()
