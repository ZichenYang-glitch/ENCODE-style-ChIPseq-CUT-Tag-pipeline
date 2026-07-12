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
from encode_pipeline.persistence.models import RunQcMetricRow
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


def _artifact(artifact_id, output_type="qc_summary", *, run_id="run-1"):
    return RunArtifactRef(
        artifact_id=artifact_id,
        run_id=run_id,
        artifact_type="file",
        name=f"{artifact_id}.tsv",
        uri=f"run://runs/{run_id}/artifacts/{artifact_id}",
        mime_type="text/tab-separated-values",
        produced_at=NOW,
        metadata={
            "output_type": output_type,
            "relative_path": f"results/{artifact_id}.tsv",
            "size_bytes": 10,
        },
    )


def _metric(
    metric_key="sequencing.total_reads",
    value=Decimal("9007199254740993"),
    *,
    scope="sample",
    sample_id="S1",
    experiment_id="EXP1",
    run_id="run-1",
    source_artifact_id="artifact-1",
):
    return RunQcMetric(
        metric_id=build_qc_metric_id(
            metric_key,
            scope,
            sample_id,
            experiment_id,
        ),
        run_id=run_id,
        metric_key=metric_key,
        display_name="Total reads",
        value=value,
        unit="count",
        scope=scope,
        sample_id=sample_id,
        experiment_id=experiment_id,
        assay="chipseq",
        qc_flag=None,
        source_artifact_id=source_artifact_id,
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


def test_repositories_page_qc_metrics_in_stable_run_scoped_id_order(repository):
    artifacts = (_artifact("artifact-1"),)
    _prepare(repository, artifacts)
    repository.create_run(_run("run-2"), _event("status_changed"))
    other_artifacts = (_artifact("artifact-other", run_id="run-2"),)
    repository.replace_artifacts(
        "run-2",
        other_artifacts,
        expected_status=RunStatus.SUCCEEDED,
        event=_event("artifacts_indexed", context={"artifact_count": 1}),
    )
    metrics = (
        _metric("sequencing.total_reads"),
        _metric("peaks.count"),
        _metric("library.pbc2"),
    )
    repository.replace_qc_metrics(
        "run-1",
        metrics,
        expected_artifacts=artifacts,
        expected_status=RunStatus.SUCCEEDED,
        event=_event("qc_metrics_indexed", context={"metric_count": 3}),
    )
    other_metric = _metric(
        "library.nrf",
        run_id="run-2",
        source_artifact_id="artifact-other",
    )
    repository.replace_qc_metrics(
        "run-2",
        (other_metric,),
        expected_artifacts=other_artifacts,
        expected_status=RunStatus.SUCCEEDED,
        event=_event("qc_metrics_indexed", context={"metric_count": 1}),
    )
    ordered = tuple(sorted(metrics, key=lambda metric: metric.metric_id))

    assert repository.list_qc_metrics("run-1") == ordered
    assert repository.list_qc_metrics("run-1", limit=2) == ordered[:2]
    assert (
        repository.list_qc_metrics(
            "run-1",
            after=ordered[1].metric_id,
            limit=2,
        )
        == ordered[2:]
    )
    assert other_metric not in repository.list_qc_metrics("run-1")

    with pytest.raises(KeyError):
        repository.list_qc_metrics("run-1", after=other_metric.metric_id)
    with pytest.raises(KeyError):
        repository.list_qc_metrics(
            "run-1",
            after=build_qc_metric_id("library.pbc1", "sample", "S1", "EXP1"),
        )


def test_repositories_validate_a_qc_cursor_row_before_skipping_it(repository):
    artifacts = (_artifact("artifact-1"),)
    _prepare(repository, artifacts)
    metrics = (
        _metric("sequencing.total_reads"),
        _metric("peaks.count"),
    )
    repository.replace_qc_metrics(
        "run-1",
        metrics,
        expected_artifacts=artifacts,
        expected_status=RunStatus.SUCCEEDED,
        event=_event("qc_metrics_indexed", context={"metric_count": 2}),
    )
    ordered = tuple(sorted(metrics, key=lambda metric: metric.metric_id))
    cursor = ordered[0]
    if isinstance(repository, InMemoryRunRepository):
        repository._qc_metrics["run-1"][cursor.metric_id] = replace(
            cursor,
            display_name="Private/path",
        )
    else:
        with repository._session_factory.begin() as session:
            row = (
                session.query(RunQcMetricRow)
                .filter_by(
                    run_id="run-1",
                    metric_id=cursor.metric_id,
                )
                .one()
            )
            row.display_name = "Private/path"

    with pytest.raises(ValueError):
        repository.list_qc_metrics("run-1", after=cursor.metric_id, limit=1)


@pytest.mark.parametrize("corruption", ("run_id", "storage_key"))
def test_in_memory_qc_cursor_rejects_corrupt_storage_identity(corruption):
    repository = InMemoryRunRepository()
    artifacts = (_artifact("artifact-1"),)
    _prepare(repository, artifacts)
    metrics = (
        _metric("sequencing.total_reads"),
        _metric("peaks.count"),
    )
    repository.replace_qc_metrics(
        "run-1",
        metrics,
        expected_artifacts=artifacts,
        expected_status=RunStatus.SUCCEEDED,
        event=_event("qc_metrics_indexed", context={"metric_count": 2}),
    )
    cursor = min(metrics, key=lambda metric: metric.metric_id)
    if corruption == "run_id":
        corrupted = replace(cursor, run_id="run-other")
    else:
        other_key = "library.nrf"
        corrupted = replace(
            cursor,
            metric_id=build_qc_metric_id(
                other_key,
                cursor.scope,
                cursor.sample_id,
                cursor.experiment_id,
            ),
            metric_key=other_key,
        )
    repository._qc_metrics["run-1"][cursor.metric_id] = corrupted

    with pytest.raises(ValueError):
        repository.list_qc_metrics("run-1", after=cursor.metric_id, limit=1)


@pytest.mark.parametrize("corruption", ("run_id", "storage_key"))
def test_in_memory_selected_qc_row_rejects_corrupt_storage_identity(corruption):
    repository = InMemoryRunRepository()
    artifacts = (_artifact("artifact-1"),)
    _prepare(repository, artifacts)
    metric = _metric()
    repository.replace_qc_metrics(
        "run-1",
        (metric,),
        expected_artifacts=artifacts,
        expected_status=RunStatus.SUCCEEDED,
        event=_event("qc_metrics_indexed", context={"metric_count": 1}),
    )
    if corruption == "run_id":
        corrupted = replace(metric, run_id="run-other")
    else:
        other_key = "library.nrf"
        corrupted = replace(
            metric,
            metric_id=build_qc_metric_id(
                other_key,
                metric.scope,
                metric.sample_id,
                metric.experiment_id,
            ),
            metric_key=other_key,
        )
    repository._qc_metrics["run-1"][metric.metric_id] = corrupted

    with pytest.raises(ValueError):
        repository.list_qc_metrics("run-1", limit=1)


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


def test_repositories_reject_well_formed_metric_id_for_different_semantics(
    repository,
):
    artifacts = (_artifact("artifact-1"),)
    _prepare(repository, artifacts)
    metric = _metric()
    wrong_id = build_qc_metric_id(
        "peaks.count",
        metric.scope,
        metric.sample_id,
        metric.experiment_id,
    )

    with pytest.raises(ValueError):
        repository.replace_qc_metrics(
            "run-1",
            (replace(metric, metric_id=wrong_id),),
            expected_artifacts=artifacts,
            expected_status=RunStatus.SUCCEEDED,
            event=_event("qc_metrics_indexed"),
        )

    assert repository.list_qc_metrics("run-1") == ()


def test_repositories_cannot_persist_same_semantics_under_two_legal_hashes(
    repository,
):
    artifacts = (_artifact("artifact-1"),)
    _prepare(repository, artifacts)
    metric = _metric()
    second_legal_hash = build_qc_metric_id(
        "library.pbc2",
        metric.scope,
        metric.sample_id,
        metric.experiment_id,
    )

    with pytest.raises(ValueError):
        repository.replace_qc_metrics(
            "run-1",
            (metric, replace(metric, metric_id=second_legal_hash)),
            expected_artifacts=artifacts,
            expected_status=RunStatus.SUCCEEDED,
            event=_event("qc_metrics_indexed"),
        )

    assert repository.list_qc_metrics("run-1") == ()


@pytest.mark.parametrize("sample_id", ("-S1", ".S1", "_S1"))
def test_repositories_accept_canonical_leading_punctuation_identifiers(
    repository,
    sample_id,
):
    artifacts = (_artifact("artifact-1"),)
    _prepare(repository, artifacts)
    metric = _metric(sample_id=sample_id)

    repository.replace_qc_metrics(
        "run-1",
        (metric,),
        expected_artifacts=artifacts,
        expected_status=RunStatus.SUCCEEDED,
        event=_event("qc_metrics_indexed"),
    )

    assert repository.list_qc_metrics("run-1") == (metric,)


@pytest.mark.parametrize("sample_id", ("/private", "..", "bad\nidentifier"))
def test_repositories_reject_unsafe_identifiers(repository, sample_id):
    artifacts = (_artifact("artifact-1"),)
    _prepare(repository, artifacts)
    metric = _metric(sample_id=sample_id)

    with pytest.raises(ValueError):
        repository.replace_qc_metrics(
            "run-1",
            (metric,),
            expected_artifacts=artifacts,
            expected_status=RunStatus.SUCCEEDED,
            event=_event("qc_metrics_indexed"),
        )

    assert repository.list_qc_metrics("run-1") == ()


@pytest.mark.parametrize(
    "changes",
    (
        {"metric_id": "metric-unsafe"},
        {"metric_id": None},
        {"metric_key": "../unsafe"},
        {"display_name": "Private/path"},
        {"unit": "percent"},
        {"unit": None},
        {"scope": "project"},
        {"scope": "run"},
        {"scope": "experiment", "sample_id": "S1", "experiment_id": "EXP1"},
        {"sample_id": "S1/private"},
        {"assay": "chipseq/private"},
        {"qc_flag": "unknown"},
        {"qc_flag": ["pass"]},
        {"produced_at": datetime(2026, 7, 12)},
    ),
)
def test_repositories_reject_invalid_durable_metric_fields_atomically(
    repository,
    changes,
):
    artifacts = (_artifact("artifact-1"),)
    _prepare(repository, artifacts)

    with pytest.raises(ValueError):
        repository.replace_qc_metrics(
            "run-1",
            (replace(_metric(), **changes),),
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
        _metric("sequencing.total_reads", Decimal("9007199254740993")),
        replace(
            _metric("library.pbc2"),
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
