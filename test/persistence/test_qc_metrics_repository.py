"""Repository parity tests for complete durable QC metric indexes."""

from __future__ import annotations

from dataclasses import replace
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from decimal import Decimal
from threading import Barrier

import pytest

from encode_pipeline.persistence import (
    SqlAlchemyRunRepository,
    create_database_engine,
    create_session_factory,
    upgrade_database,
)
from encode_pipeline.platform.runs import (
    RunArtifactRef,
    RunQcMetric,
    RunRecord,
    RunStatus,
)
from encode_pipeline.services.run_repositories import (
    ConcurrentRunUpdateError,
    InMemoryRunRepository,
    RunEventDraft,
)


NOW = datetime(2026, 7, 12, tzinfo=timezone.utc)


@pytest.fixture(params=("memory", "sql"))
def repository(request, tmp_path):
    if request.param == "memory":
        yield InMemoryRunRepository()
        return
    database_url = f"sqlite:///{tmp_path / 'platform.db'}"
    upgrade_database(database_url)
    engine = create_database_engine(database_url)
    yield SqlAlchemyRunRepository(create_session_factory(engine))
    engine.dispose()


def _run(run_id="run-1"):
    return RunRecord(
        run_id=run_id,
        workflow_id="workflow",
        inputs={"config": {}, "samples": None, "options": {}},
        status=RunStatus.SUCCEEDED,
        created_at=NOW,
        updated_at=NOW,
        started_at=NOW,
        ended_at=NOW,
        current_stage="execution",
        cancellation_reason=None,
        error=None,
    )


def _event(event_type, *, context=None):
    return RunEventDraft(
        event_type=event_type,
        message=event_type.replace("_", " "),
        status=RunStatus.SUCCEEDED,
        stage="qc_summary_indexing",
        context={} if context is None else context,
    )


def _artifact(artifact_id, output_type="qc_summary"):
    return RunArtifactRef(
        artifact_id=artifact_id,
        run_id="run-1",
        artifact_type="file",
        name=f"{artifact_id}.tsv",
        uri=f"run://runs/run-1/artifacts/{artifact_id}",
        mime_type="text/tab-separated-values",
        produced_at=NOW,
        metadata={
            "output_type": output_type,
            "relative_path": f"results/{artifact_id}.tsv",
            "size_bytes": 10,
        },
    )


def _metric(metric_id="qcmetric-1", value=Decimal("9007199254740993")):
    return RunQcMetric(
        metric_id=metric_id,
        run_id="run-1",
        metric_key="sequencing.total_reads",
        display_name="Total reads",
        value=value,
        unit="count",
        scope="sample",
        sample_id="S1",
        experiment_id="EXP1",
        assay="chipseq",
        qc_flag=None,
        source_artifact_id="artifact-1",
        produced_at=NOW,
    )


def _prepare(repository, artifacts):
    repository.create_run(_run(), _event("status_changed"))
    repository.replace_artifacts(
        "run-1",
        artifacts,
        expected_status=RunStatus.SUCCEEDED,
        event=_event("artifacts_indexed", context={"artifact_count": len(artifacts)}),
    )


def test_repositories_atomically_replace_qc_with_order_independent_sources(repository):
    artifacts = (_artifact("artifact-1"), _artifact("artifact-2", "other"))
    _prepare(repository, artifacts)
    metric = _metric()

    first = repository.replace_qc_metrics(
        "run-1",
        (metric,),
        expected_artifacts=tuple(reversed(artifacts)),
        expected_status=RunStatus.SUCCEEDED,
        event=_event("qc_metrics_indexed", context={"metric_count": 1}),
    )
    duplicate = repository.replace_qc_metrics(
        "run-1",
        (metric,),
        expected_artifacts=artifacts,
        expected_status=RunStatus.SUCCEEDED,
        event=_event("qc_metrics_indexed", context={"metric_count": 1}),
    )

    assert first is not None
    assert duplicate is None
    assert repository.list_qc_metrics("run-1") == (metric,)
    assert [event.event_type for event in repository.list_events("run-1")].count(
        "qc_metrics_indexed"
    ) == 1


def test_repositories_reject_stale_source_set_without_partial_qc_change(repository):
    artifacts = (_artifact("artifact-1"),)
    _prepare(repository, artifacts)
    original = _metric()
    repository.replace_qc_metrics(
        "run-1",
        (original,),
        expected_artifacts=artifacts,
        expected_status=RunStatus.SUCCEEDED,
        event=_event("qc_metrics_indexed", context={"metric_count": 1}),
    )
    stale = (_artifact("artifact-stale"),)
    changed = replace(original, value=Decimal("2"))

    with pytest.raises(ConcurrentRunUpdateError):
        repository.replace_qc_metrics(
            "run-1",
            (changed,),
            expected_artifacts=stale,
            expected_status=RunStatus.SUCCEEDED,
            event=_event("qc_metrics_indexed", context={"metric_count": 1}),
        )

    assert repository.list_qc_metrics("run-1") == (original,)


def test_repositories_reject_missing_metric_source_and_invalid_decimal(repository):
    artifacts = (_artifact("artifact-1"),)
    _prepare(repository, artifacts)

    with pytest.raises(ValueError):
        repository.replace_qc_metrics(
            "run-1",
            (replace(_metric(), source_artifact_id="artifact-other"),),
            expected_artifacts=artifacts,
            expected_status=RunStatus.SUCCEEDED,
            event=_event("qc_metrics_indexed"),
        )
    with pytest.raises(ValueError):
        repository.replace_qc_metrics(
            "run-1",
            (_metric(value=Decimal("0.1234567890123")),),
            expected_artifacts=artifacts,
            expected_status=RunStatus.SUCCEEDED,
            event=_event("qc_metrics_indexed"),
        )

    assert repository.list_qc_metrics("run-1") == ()


def test_changed_artifact_generation_invalidates_qc_but_reordering_does_not(repository):
    artifacts = (_artifact("artifact-1"), _artifact("artifact-2", "other"))
    _prepare(repository, artifacts)
    repository.replace_qc_metrics(
        "run-1",
        (_metric(),),
        expected_artifacts=artifacts,
        expected_status=RunStatus.SUCCEEDED,
        event=_event("qc_metrics_indexed", context={"metric_count": 1}),
    )

    reordered = repository.replace_artifacts(
        "run-1",
        tuple(reversed(artifacts)),
        expected_status=RunStatus.SUCCEEDED,
        event=_event("artifacts_indexed", context={"artifact_count": 2}),
    )
    assert reordered is None
    assert repository.list_qc_metrics("run-1") == (_metric(),)

    repository.replace_artifacts(
        "run-1",
        (_artifact("artifact-3", "other"),),
        expected_status=RunStatus.SUCCEEDED,
        event=_event("artifacts_indexed", context={"artifact_count": 1}),
    )

    assert repository.list_qc_metrics("run-1") == ()
    event_types = [event.event_type for event in repository.list_events("run-1")]
    assert event_types[-2:] == ["artifacts_indexed", "qc_metrics_invalidated"]


def test_first_artifact_index_does_not_claim_qc_was_invalidated(repository):
    _prepare(repository, (_artifact("artifact-1"),))

    assert "qc_metrics_invalidated" not in {
        event.event_type for event in repository.list_events("run-1")
    }


def test_changed_artifacts_invalidate_a_confirmed_empty_qc_index(repository):
    artifacts = (_artifact("artifact-1"),)
    _prepare(repository, artifacts)
    repository.replace_qc_metrics(
        "run-1",
        (),
        expected_artifacts=artifacts,
        expected_status=RunStatus.SUCCEEDED,
        event=_event("qc_metrics_indexed", context={"metric_count": 0}),
    )

    repository.replace_artifacts(
        "run-1",
        (_artifact("artifact-2", "other"),),
        expected_status=RunStatus.SUCCEEDED,
        event=_event("artifacts_indexed", context={"artifact_count": 1}),
    )

    assert [event.event_type for event in repository.list_events("run-1")][-2:] == [
        "artifacts_indexed",
        "qc_metrics_invalidated",
    ]


def test_repeated_identical_qc_failure_is_an_atomic_noop(repository):
    _prepare(repository, (_artifact("artifact-1"),))
    failure = _event(
        "qc_metrics_indexing_failed",
        context={"reason_code": "QC_INDEXING_VALIDATION_FAILED"},
    )

    first = repository.record_qc_metrics_failure(
        "run-1",
        reason_code="QC_INDEXING_VALIDATION_FAILED",
        expected_status=RunStatus.SUCCEEDED,
        event=failure,
    )
    second = repository.record_qc_metrics_failure(
        "run-1",
        reason_code="QC_INDEXING_VALIDATION_FAILED",
        expected_status=RunStatus.SUCCEEDED,
        event=failure,
    )

    assert first is not None
    assert second is None
    assert [event.event_type for event in repository.list_events("run-1")].count(
        "qc_metrics_indexing_failed"
    ) == 1


def test_sqlalchemy_qc_decimal_text_survives_reopen_without_float_rounding(tmp_path):
    database_url = f"sqlite:///{tmp_path / 'platform.db'}"
    upgrade_database(database_url)
    first_engine = create_database_engine(database_url)
    first = SqlAlchemyRunRepository(create_session_factory(first_engine))
    artifacts = (_artifact("artifact-1"),)
    _prepare(first, artifacts)
    metrics = (
        _metric("qcmetric-large", Decimal("9007199254740993")),
        replace(
            _metric("qcmetric-boundary"),
            metric_key="library.pbc2",
            value=Decimal("99999999999999999999999999.999999999999"),
        ),
    )
    first.replace_qc_metrics(
        "run-1",
        metrics,
        expected_artifacts=artifacts,
        expected_status=RunStatus.SUCCEEDED,
        event=_event("qc_metrics_indexed", context={"metric_count": 2}),
    )
    first_engine.dispose()

    second_engine = create_database_engine(database_url)
    reopened = SqlAlchemyRunRepository(create_session_factory(second_engine))
    try:
        assert reopened.list_qc_metrics("run-1") == tuple(
            sorted(metrics, key=lambda metric: metric.metric_id)
        )
    finally:
        second_engine.dispose()


def test_sqlalchemy_qc_replace_rolls_back_rows_when_event_insert_fails(
    tmp_path,
    monkeypatch,
):
    database_url = f"sqlite:///{tmp_path / 'platform.db'}"
    upgrade_database(database_url)
    engine = create_database_engine(database_url)
    repository = SqlAlchemyRunRepository(create_session_factory(engine))
    artifacts = (_artifact("artifact-1"),)
    _prepare(repository, artifacts)
    original = _metric()
    repository.replace_qc_metrics(
        "run-1",
        (original,),
        expected_artifacts=artifacts,
        expected_status=RunStatus.SUCCEEDED,
        event=_event("qc_metrics_indexed", context={"metric_count": 1}),
    )
    before_events = repository.list_events("run-1")

    def fail_event(*_args, **_kwargs):
        raise RuntimeError("event insert failed")

    monkeypatch.setattr(repository, "_insert_event", fail_event)
    with pytest.raises(RuntimeError, match="event insert failed"):
        repository.replace_qc_metrics(
            "run-1",
            (replace(original, value=Decimal("2")),),
            expected_artifacts=artifacts,
            expected_status=RunStatus.SUCCEEDED,
            event=_event("qc_metrics_indexed", context={"metric_count": 1}),
        )

    assert repository.list_qc_metrics("run-1") == (original,)
    assert repository.list_events("run-1") == before_events
    engine.dispose()


def test_sqlalchemy_artifact_and_qc_replace_race_cannot_leave_orphan_sources(
    tmp_path,
):
    database_url = f"sqlite:///{tmp_path / 'platform.db'}"
    upgrade_database(database_url)
    setup_engine = create_database_engine(database_url)
    setup = SqlAlchemyRunRepository(create_session_factory(setup_engine))
    original_artifacts = (_artifact("artifact-1"),)
    _prepare(setup, original_artifacts)
    setup_engine.dispose()
    qc_engine = create_database_engine(database_url)
    artifact_engine = create_database_engine(database_url)
    qc_repository = SqlAlchemyRunRepository(create_session_factory(qc_engine))
    artifact_repository = SqlAlchemyRunRepository(
        create_session_factory(artifact_engine)
    )
    barrier = Barrier(2)

    def replace_qc():
        barrier.wait(timeout=5)
        try:
            qc_repository.replace_qc_metrics(
                "run-1",
                (_metric(),),
                expected_artifacts=original_artifacts,
                expected_status=RunStatus.SUCCEEDED,
                event=_event("qc_metrics_indexed", context={"metric_count": 1}),
            )
        except ConcurrentRunUpdateError:
            return "stale"
        return "indexed"

    def replace_artifacts():
        barrier.wait(timeout=5)
        artifact_repository.replace_artifacts(
            "run-1",
            (_artifact("artifact-new", "other"),),
            expected_status=RunStatus.SUCCEEDED,
            event=_event("artifacts_indexed", context={"artifact_count": 1}),
        )

    try:
        with ThreadPoolExecutor(max_workers=2) as pool:
            qc_future = pool.submit(replace_qc)
            artifact_future = pool.submit(replace_artifacts)
            qc_future.result(timeout=10)
            artifact_future.result(timeout=10)

        assert qc_repository.list_qc_metrics("run-1") == ()
        assert qc_repository.list_artifacts("run-1") == (
            _artifact("artifact-new", "other"),
        )
    finally:
        qc_engine.dispose()
        artifact_engine.dispose()
