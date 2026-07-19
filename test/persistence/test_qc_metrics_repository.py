"""Repository parity tests for complete durable QC metric indexes."""

from __future__ import annotations

from dataclasses import replace
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timedelta, timezone
from decimal import Decimal
from hashlib import sha256
from threading import Barrier, Event

import pytest
from sqlalchemy import event as sqlalchemy_event

from encode_pipeline.persistence import (
    SqlAlchemyRunRepository,
    create_database_engine,
    create_session_factory,
    upgrade_database,
)
from encode_pipeline.persistence.models import RunQcMetricRow
from encode_pipeline.platform.result_generations import (
    artifact_manifest_digest,
    new_result_attempt_id,
    qc_metric_manifest_digest,
)
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
    ResultGenerationChangedError,
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


def _artifact(
    artifact_id,
    output_type="qc_summary",
    *,
    run_id="run-1",
    revision_seed=None,
    produced_at=NOW,
):
    return RunArtifactRef(
        artifact_id=artifact_id,
        run_id=run_id,
        artifact_type="file",
        name=f"{artifact_id}.tsv",
        uri=f"run://runs/{run_id}/artifacts/{artifact_id}",
        mime_type="text/tab-separated-values",
        produced_at=produced_at,
        revision="artifactrev-"
        + sha256(
            f"{artifact_id}:{output_type}:{revision_seed or 'initial'}".encode()
        ).hexdigest(),
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
    unit="count",
    produced_at=NOW,
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
        unit=unit,
        scope=scope,
        sample_id=sample_id,
        experiment_id=experiment_id,
        assay="chipseq",
        qc_flag=None,
        source_artifact_id=source_artifact_id,
        produced_at=produced_at,
    )


def _prepare(repository, artifacts):
    repository.create_run(_run(), _event("status_changed"))
    attempt_id = new_result_attempt_id()
    repository.begin_artifact_result_attempt(
        "run-1",
        attempt_id=attempt_id,
        expected_status=RunStatus.SUCCEEDED,
    )
    repository.replace_artifacts(
        "run-1",
        artifacts,
        attempt_id=attempt_id,
        expected_status=RunStatus.SUCCEEDED,
        event=_event("artifacts_indexed", context={"artifact_count": len(artifacts)}),
    )


def test_persisted_artifact_metadata_cannot_mutate_a_generation_in_place(repository):
    artifact = RunArtifactRef(
        artifact_id="artifact-1",
        run_id="run-1",
        artifact_type="file",
        name="artifact-1.tsv",
        uri="run://runs/run-1/artifacts/artifact-1",
        mime_type="text/tab-separated-values",
        produced_at=NOW,
        revision="artifactrev-" + "a" * 64,
        metadata={
            "output_type": "qc_summary",
            "relative_path": "results/artifact-1.tsv",
            "size_bytes": 10,
            "identity": {"labels": ["initial"]},
        },
    )
    _prepare(repository, (artifact,))
    before = repository.get_result_state("run-1")
    persisted = repository.list_artifacts("run-1")[0]

    with pytest.raises(TypeError):
        persisted.metadata["size_bytes"] = 11
    with pytest.raises(AttributeError):
        persisted.metadata["identity"]["labels"].append("mutated")

    after = repository.get_result_state("run-1")
    reread = repository.list_artifacts("run-1")[0]
    assert after == before
    assert artifact_manifest_digest((reread,)) == before.artifact_manifest_digest
    assert reread.to_dict()["metadata"]["identity"]["labels"] == ["initial"]


def _begin_artifact_attempt(repository, *, run_id="run-1"):
    attempt_id = new_result_attempt_id()
    repository.begin_artifact_result_attempt(
        run_id,
        attempt_id=attempt_id,
        expected_status=RunStatus.SUCCEEDED,
    )
    return attempt_id


def _begin_qc_attempt(repository, *, run_id="run-1"):
    state = repository.get_result_state(run_id)
    assert state.artifact_generation is not None
    artifacts = repository.list_artifacts(run_id)
    attempt_id = new_result_attempt_id()
    repository.begin_qc_result_attempt(
        run_id,
        attempt_id=attempt_id,
        expected_artifact_generation=state.artifact_generation,
        expected_artifacts=artifacts,
        expected_status=RunStatus.SUCCEEDED,
    )
    return attempt_id, state.artifact_generation


def _replace_qc(
    repository,
    metrics,
    *,
    artifacts=None,
    expected_artifacts=None,
    run_id="run-1",
    expected_status=RunStatus.SUCCEEDED,
    event=None,
):
    if artifacts is None:
        artifacts = expected_artifacts
    assert artifacts is not None
    attempt_id, artifact_generation = _begin_qc_attempt(
        repository,
        run_id=run_id,
    )
    return repository.replace_qc_metrics(
        run_id,
        metrics,
        attempt_id=attempt_id,
        expected_artifact_generation=artifact_generation,
        expected_artifacts=artifacts,
        expected_status=expected_status,
        event=(
            _event("qc_metrics_indexed", context={"metric_count": len(metrics)})
            if event is None
            else event
        ),
    )


def test_same_length_artifact_revision_replacement_advances_generation_and_invalidates_qc(
    repository,
):
    first = (_artifact("artifact-1", revision_seed="same-size-a"),)
    _prepare(repository, first)
    _replace_qc(repository, (_metric(),), artifacts=first)
    before = repository.get_result_state("run-1")
    assert before.artifact_generation is not None
    assert before.qc_generation is not None

    second = (_artifact("artifact-1", revision_seed="same-size-b"),)
    attempt_id = new_result_attempt_id()
    repository.begin_artifact_result_attempt(
        "run-1",
        attempt_id=attempt_id,
        expected_status=RunStatus.SUCCEEDED,
    )
    repository.replace_artifacts(
        "run-1",
        second,
        attempt_id=attempt_id,
        expected_status=RunStatus.SUCCEEDED,
        event=_event("artifacts_indexed"),
    )

    after = repository.get_result_state("run-1")
    assert after.artifact_revision == before.artifact_revision + 1
    assert after.artifact_generation != before.artifact_generation
    assert after.qc_revision == before.qc_revision + 1
    assert after.qc_generation is None
    assert after.qc_outcome == "invalidated"
    assert repository.list_qc_metrics("run-1") == ()


def test_artifact_page_read_is_generation_bound_when_ids_and_count_are_reused(
    repository,
):
    first = (
        _artifact("artifact-a", revision_seed="generation-a"),
        _artifact("artifact-b", revision_seed="generation-a"),
    )
    _prepare(repository, first)
    generation_a = repository.get_result_state("run-1").artifact_generation
    assert generation_a is not None
    observed_generation, page = repository.list_artifacts_page(
        "run-1",
        expected_generation=generation_a,
        limit=1,
    )
    assert observed_generation == generation_a
    assert tuple(item.artifact_id for item in page) == ("artifact-a",)

    second = (
        _artifact("artifact-a", revision_seed="generation-b"),
        _artifact("artifact-b", revision_seed="generation-b"),
    )
    attempt_id = _begin_artifact_attempt(repository)
    repository.replace_artifacts(
        "run-1",
        second,
        attempt_id=attempt_id,
        expected_status=RunStatus.SUCCEEDED,
        event=_event("artifacts_indexed"),
    )

    with pytest.raises(ResultGenerationChangedError):
        repository.list_artifacts_page(
            "run-1",
            expected_generation=generation_a,
            after="artifact-a",
            limit=1,
        )


def test_artifact_detail_checks_generation_before_deleted_row(repository):
    first = (_artifact("artifact-a", revision_seed="generation-a"),)
    _prepare(repository, first)
    generation_a = repository.get_result_state("run-1").artifact_generation
    assert generation_a is not None

    attempt_id = _begin_artifact_attempt(repository)
    repository.replace_artifacts(
        "run-1",
        (),
        attempt_id=attempt_id,
        expected_status=RunStatus.SUCCEEDED,
        event=_event("artifacts_indexed"),
    )
    generation_b = repository.get_result_state("run-1").artifact_generation
    assert generation_b is not None and generation_b != generation_a

    with pytest.raises(ResultGenerationChangedError):
        repository.get_artifact_at_generation(
            "run-1",
            "artifact-a",
            expected_generation=generation_a,
        )
    with pytest.raises(KeyError):
        repository.get_artifact_at_generation(
            "run-1",
            "artifact-a",
            expected_generation=generation_b,
        )


def test_result_timestamps_are_utc_canonical_and_match_persisted_generations(
    repository,
):
    offset = timezone(timedelta(hours=8))
    local_time = datetime(2026, 1, 2, 3, 4, tzinfo=offset)
    expected_utc = datetime(2026, 1, 1, 19, 4, tzinfo=timezone.utc)
    artifacts = (_artifact("artifact-1", produced_at=local_time),)
    _prepare(repository, artifacts)

    persisted_artifacts = repository.list_artifacts("run-1")
    artifact_state = repository.get_result_state("run-1")
    assert persisted_artifacts[0].produced_at == expected_utc
    assert artifact_state.artifact_manifest_digest == artifact_manifest_digest(
        persisted_artifacts
    )

    metrics = (_metric(produced_at=local_time),)
    _replace_qc(repository, metrics, artifacts=persisted_artifacts)
    persisted_metrics = repository.list_qc_metrics("run-1")
    qc_state = repository.get_result_state("run-1")
    assert persisted_metrics[0].produced_at == expected_utc
    assert qc_state.qc_manifest_digest == qc_metric_manifest_digest(persisted_metrics)


def test_stale_qc_attempt_success_and_failure_cannot_override_newer_success(repository):
    artifacts = (_artifact("artifact-1"),)
    _prepare(repository, artifacts)
    generation = repository.get_result_state("run-1").artifact_generation
    assert generation is not None
    old_attempt = "resultattempt-" + "a" * 64
    new_attempt = "resultattempt-" + "b" * 64
    for attempt_id in (old_attempt, new_attempt):
        repository.begin_qc_result_attempt(
            "run-1",
            attempt_id=attempt_id,
            expected_artifact_generation=generation,
            expected_artifacts=artifacts,
            expected_status=RunStatus.SUCCEEDED,
        )
    metrics = (_metric(value=Decimal("11")),)
    repository.replace_qc_metrics(
        "run-1",
        metrics,
        attempt_id=new_attempt,
        expected_artifact_generation=generation,
        expected_artifacts=artifacts,
        expected_status=RunStatus.SUCCEEDED,
        event=_event("qc_metrics_indexed"),
    )
    events = repository.list_events("run-1")

    with pytest.raises(ConcurrentRunUpdateError):
        repository.begin_qc_result_attempt(
            "run-1",
            attempt_id=old_attempt,
            expected_artifact_generation=generation,
            expected_artifacts=artifacts,
            expected_status=RunStatus.SUCCEEDED,
        )
    with pytest.raises(ConcurrentRunUpdateError):
        repository.replace_qc_metrics(
            "run-1",
            (_metric(value=Decimal("10")),),
            attempt_id=old_attempt,
            expected_artifact_generation=generation,
            expected_artifacts=artifacts,
            expected_status=RunStatus.SUCCEEDED,
            event=_event("qc_metrics_indexed"),
        )
    with pytest.raises(ConcurrentRunUpdateError):
        repository.record_qc_metrics_failure(
            "run-1",
            attempt_id=old_attempt,
            expected_artifact_generation=generation,
            reason_code="QC_OLD_FAILURE",
            expected_status=RunStatus.SUCCEEDED,
            event=_event("qc_metrics_indexing_failed"),
        )
    assert repository.list_qc_metrics("run-1") == metrics
    assert repository.list_events("run-1") == events


def test_superseded_artifact_attempt_cannot_be_registered_again(repository):
    artifacts = (_artifact("artifact-1"),)
    _prepare(repository, artifacts)
    first_attempt = repository.get_result_state("run-1").artifact_attempt_id
    assert first_attempt is not None
    second_attempt = _begin_artifact_attempt(repository)
    repository.replace_artifacts(
        "run-1",
        artifacts,
        attempt_id=second_attempt,
        expected_status=RunStatus.SUCCEEDED,
        event=_event("artifacts_indexed"),
    )
    before = repository.get_result_state("run-1")

    with pytest.raises(ConcurrentRunUpdateError):
        repository.begin_artifact_result_attempt(
            "run-1",
            attempt_id=first_attempt,
            expected_status=RunStatus.SUCCEEDED,
        )

    assert repository.get_result_state("run-1") == before


def test_sqlite_reopen_preserves_superseded_qc_attempt_tombstone(tmp_path):
    database_url = f"sqlite:///{tmp_path / 'attempts.db'}"
    upgrade_database(database_url)
    first_engine = create_database_engine(database_url)
    first = SqlAlchemyRunRepository(create_session_factory(first_engine))
    artifacts = (_artifact("artifact-1"),)
    _prepare(first, artifacts)
    generation = first.get_result_state("run-1").artifact_generation
    assert generation is not None
    old_attempt = new_result_attempt_id()
    new_attempt = new_result_attempt_id()
    for attempt_id in (old_attempt, new_attempt):
        first.begin_qc_result_attempt(
            "run-1",
            attempt_id=attempt_id,
            expected_artifact_generation=generation,
            expected_artifacts=artifacts,
            expected_status=RunStatus.SUCCEEDED,
        )
    first.replace_qc_metrics(
        "run-1",
        (_metric(),),
        attempt_id=new_attempt,
        expected_artifact_generation=generation,
        expected_artifacts=artifacts,
        expected_status=RunStatus.SUCCEEDED,
        event=_event("qc_metrics_indexed"),
    )
    first_engine.dispose()

    reopened_engine = create_database_engine(database_url)
    reopened = SqlAlchemyRunRepository(create_session_factory(reopened_engine))
    with pytest.raises(ConcurrentRunUpdateError):
        reopened.begin_qc_result_attempt(
            "run-1",
            attempt_id=old_attempt,
            expected_artifact_generation=generation,
            expected_artifacts=artifacts,
            expected_status=RunStatus.SUCCEEDED,
        )
    assert reopened.get_result_state("run-1").qc_outcome == "succeeded"
    assert reopened.list_qc_metrics("run-1") == (_metric(),)
    reopened_engine.dispose()


def test_qc_attempt_cannot_bind_a_stale_manifest_to_the_current_generation(repository):
    first = (_artifact("artifact-1", revision_seed="first"),)
    _prepare(repository, first)
    second = (_artifact("artifact-1", revision_seed="second"),)
    artifact_attempt = _begin_artifact_attempt(repository)
    repository.replace_artifacts(
        "run-1",
        second,
        attempt_id=artifact_attempt,
        expected_status=RunStatus.SUCCEEDED,
        event=_event("artifacts_indexed"),
    )
    metrics = (_metric(value=Decimal("11")),)
    _replace_qc(repository, metrics, artifacts=second)
    before_state = repository.get_result_state("run-1")
    before_events = repository.list_events("run-1")

    with pytest.raises(ResultGenerationChangedError):
        repository.begin_qc_result_attempt(
            "run-1",
            attempt_id=new_result_attempt_id(),
            expected_artifact_generation=before_state.artifact_generation,
            expected_artifacts=first,
            expected_status=RunStatus.SUCCEEDED,
        )

    assert repository.get_result_state("run-1") == before_state
    assert repository.list_qc_metrics("run-1") == metrics
    assert repository.list_events("run-1") == before_events


def test_same_attempt_failure_after_qc_success_is_a_noop(repository):
    artifacts = (_artifact("artifact-1"),)
    _prepare(repository, artifacts)
    generation = repository.get_result_state("run-1").artifact_generation
    assert generation is not None
    attempt_id = new_result_attempt_id()
    repository.begin_qc_result_attempt(
        "run-1",
        attempt_id=attempt_id,
        expected_artifact_generation=generation,
        expected_artifacts=artifacts,
        expected_status=RunStatus.SUCCEEDED,
    )
    repository.replace_qc_metrics(
        "run-1",
        (_metric(),),
        attempt_id=attempt_id,
        expected_artifact_generation=generation,
        expected_artifacts=artifacts,
        expected_status=RunStatus.SUCCEEDED,
        event=_event("qc_metrics_indexed"),
    )
    before = repository.list_events("run-1")

    result = repository.record_qc_metrics_failure(
        "run-1",
        attempt_id=attempt_id,
        expected_artifact_generation=generation,
        reason_code="QC_HARD_TIMEOUT",
        expected_status=RunStatus.SUCCEEDED,
        event=_event("qc_metrics_indexing_failed"),
    )
    assert result is None
    assert repository.list_events("run-1") == before


def test_qc_page_read_is_generation_bound_when_metric_ids_and_count_are_reused(
    repository,
):
    artifacts = (_artifact("artifact-1"),)
    _prepare(repository, artifacts)
    first_metrics = (_metric(value=Decimal("10")),)
    _replace_qc(repository, first_metrics, artifacts=artifacts)
    first_generation = repository.get_result_state("run-1").qc_generation
    assert first_generation is not None
    assert repository.list_qc_metrics_page(
        "run-1",
        expected_generation=first_generation,
        limit=10,
    ) == (first_generation, first_metrics)

    second_metrics = (_metric(value=Decimal("11")),)
    _replace_qc(repository, second_metrics, artifacts=artifacts)
    second_generation = repository.get_result_state("run-1").qc_generation
    assert second_generation is not None
    assert second_generation != first_generation
    with pytest.raises(ResultGenerationChangedError):
        repository.list_qc_metrics_page(
            "run-1",
            expected_generation=first_generation,
            limit=10,
        )
    assert repository.list_qc_metrics_page(
        "run-1",
        expected_generation=second_generation,
        limit=10,
    ) == (second_generation, second_metrics)


def test_repositories_page_qc_metrics_in_stable_run_scoped_id_order(repository):
    artifacts = (_artifact("artifact-1"),)
    _prepare(repository, artifacts)
    repository.create_run(_run("run-2"), _event("status_changed"))
    other_artifacts = (_artifact("artifact-other", run_id="run-2"),)
    other_artifact_attempt = _begin_artifact_attempt(repository, run_id="run-2")
    repository.replace_artifacts(
        "run-2",
        other_artifacts,
        attempt_id=other_artifact_attempt,
        expected_status=RunStatus.SUCCEEDED,
        event=_event("artifacts_indexed", context={"artifact_count": 1}),
    )
    metrics = (
        _metric("sequencing.total_reads"),
        _metric("peaks.count"),
        _metric("library.pbc2"),
    )
    _replace_qc(
        repository,
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
    _replace_qc(
        repository,
        (other_metric,),
        expected_artifacts=other_artifacts,
        run_id="run-2",
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


def test_repositories_round_trip_a_workflow_neutral_score_unit(repository):
    artifacts = (_artifact("artifact-1"),)
    _prepare(repository, artifacts)
    metric = _metric(
        "rseqc.tin.mean_score",
        value=Decimal("72.125"),
        unit="score",
    )

    _replace_qc(
        repository,
        (metric,),
        expected_artifacts=artifacts,
        expected_status=RunStatus.SUCCEEDED,
        event=_event("qc_metrics_indexed", context={"metric_count": 1}),
    )

    assert repository.list_qc_metrics("run-1") == (metric,)


def test_repositories_validate_a_qc_cursor_row_before_skipping_it(repository):
    artifacts = (_artifact("artifact-1"),)
    _prepare(repository, artifacts)
    metrics = (
        _metric("sequencing.total_reads"),
        _metric("peaks.count"),
    )
    _replace_qc(
        repository,
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
    _replace_qc(
        repository,
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
    _replace_qc(
        repository,
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

    first = _replace_qc(
        repository,
        (metric,),
        expected_artifacts=tuple(reversed(artifacts)),
        expected_status=RunStatus.SUCCEEDED,
        event=_event("qc_metrics_indexed", context={"metric_count": 1}),
    )
    duplicate = _replace_qc(
        repository,
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
    _replace_qc(
        repository,
        (original,),
        expected_artifacts=artifacts,
        expected_status=RunStatus.SUCCEEDED,
        event=_event("qc_metrics_indexed", context={"metric_count": 1}),
    )
    stale = (_artifact("artifact-stale"),)
    changed = replace(original, value=Decimal("2"))

    with pytest.raises(ConcurrentRunUpdateError):
        _replace_qc(
            repository,
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
        _replace_qc(
            repository,
            (replace(_metric(), source_artifact_id="artifact-other"),),
            expected_artifacts=artifacts,
            expected_status=RunStatus.SUCCEEDED,
            event=_event("qc_metrics_indexed"),
        )
    with pytest.raises(ValueError):
        _replace_qc(
            repository,
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
        _replace_qc(
            repository,
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
        _replace_qc(
            repository,
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

    _replace_qc(
        repository,
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
        _replace_qc(
            repository,
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
        _replace_qc(
            repository,
            (replace(_metric(), **changes),),
            expected_artifacts=artifacts,
            expected_status=RunStatus.SUCCEEDED,
            event=_event("qc_metrics_indexed"),
        )

    assert repository.list_qc_metrics("run-1") == ()


def test_changed_artifact_generation_invalidates_qc_but_reordering_does_not(repository):
    artifacts = (_artifact("artifact-1"), _artifact("artifact-2", "other"))
    _prepare(repository, artifacts)
    _replace_qc(
        repository,
        (_metric(),),
        expected_artifacts=artifacts,
        expected_status=RunStatus.SUCCEEDED,
        event=_event("qc_metrics_indexed", context={"metric_count": 1}),
    )

    reorder_attempt = _begin_artifact_attempt(repository)
    reordered = repository.replace_artifacts(
        "run-1",
        tuple(reversed(artifacts)),
        attempt_id=reorder_attempt,
        expected_status=RunStatus.SUCCEEDED,
        event=_event("artifacts_indexed", context={"artifact_count": 2}),
    )
    assert reordered is None
    assert repository.list_qc_metrics("run-1") == (_metric(),)

    replacement_attempt = _begin_artifact_attempt(repository)
    repository.replace_artifacts(
        "run-1",
        (_artifact("artifact-3", "other"),),
        attempt_id=replacement_attempt,
        expected_status=RunStatus.SUCCEEDED,
        event=_event("artifacts_indexed", context={"artifact_count": 1}),
    )

    assert repository.list_qc_metrics("run-1") == ()
    final_events = repository.list_events("run-1")[-2:]
    assert [event.event_type for event in final_events] == [
        "qc_metrics_invalidated",
        "artifacts_indexed",
    ]
    assert {event.context["artifact_generation"] for event in final_events} == {
        repository.get_result_state("run-1").artifact_generation
    }


def test_first_artifact_index_does_not_claim_qc_was_invalidated(repository):
    _prepare(repository, (_artifact("artifact-1"),))

    assert "qc_metrics_invalidated" not in {
        event.event_type for event in repository.list_events("run-1")
    }


def test_changed_artifacts_invalidate_a_confirmed_empty_qc_index(repository):
    artifacts = (_artifact("artifact-1"),)
    _prepare(repository, artifacts)
    _replace_qc(
        repository,
        (),
        expected_artifacts=artifacts,
        expected_status=RunStatus.SUCCEEDED,
        event=_event("qc_metrics_indexed", context={"metric_count": 0}),
    )

    replacement_attempt = _begin_artifact_attempt(repository)
    repository.replace_artifacts(
        "run-1",
        (_artifact("artifact-2", "other"),),
        attempt_id=replacement_attempt,
        expected_status=RunStatus.SUCCEEDED,
        event=_event("artifacts_indexed", context={"artifact_count": 1}),
    )

    final_events = repository.list_events("run-1")[-2:]
    assert [event.event_type for event in final_events] == [
        "qc_metrics_invalidated",
        "artifacts_indexed",
    ]
    assert {event.context["artifact_generation"] for event in final_events} == {
        repository.get_result_state("run-1").artifact_generation
    }


def test_repeated_identical_qc_failure_is_an_atomic_noop(repository):
    _prepare(repository, (_artifact("artifact-1"),))
    failure = _event(
        "qc_metrics_indexing_failed",
        context={"reason_code": "QC_INDEXING_VALIDATION_FAILED"},
    )
    attempt_id, artifact_generation = _begin_qc_attempt(repository)

    first = repository.record_qc_metrics_failure(
        "run-1",
        attempt_id=attempt_id,
        expected_artifact_generation=artifact_generation,
        reason_code="QC_INDEXING_VALIDATION_FAILED",
        expected_status=RunStatus.SUCCEEDED,
        event=failure,
    )
    second = repository.record_qc_metrics_failure(
        "run-1",
        attempt_id=attempt_id,
        expected_artifact_generation=artifact_generation,
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
    _replace_qc(
        first,
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
    _replace_qc(
        repository,
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
        _replace_qc(
            repository,
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
    qc_attempt, artifact_generation = _begin_qc_attempt(qc_repository)
    artifact_attempt = _begin_artifact_attempt(artifact_repository)
    barrier = Barrier(2)

    def replace_qc():
        barrier.wait(timeout=5)
        try:
            qc_repository.replace_qc_metrics(
                "run-1",
                (_metric(),),
                attempt_id=qc_attempt,
                expected_artifact_generation=artifact_generation,
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
            attempt_id=artifact_attempt,
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


def test_sqlalchemy_generation_bound_page_read_is_atomic_during_replacement(
    tmp_path,
):
    database_url = f"sqlite:///{tmp_path / 'platform.db'}"
    upgrade_database(database_url)
    setup_engine = create_database_engine(database_url)
    setup = SqlAlchemyRunRepository(create_session_factory(setup_engine))
    artifacts = (_artifact("artifact-1"),)
    _prepare(setup, artifacts)
    original = tuple(
        sorted(
            (
                _metric("sequencing.total_reads", value=Decimal("10")),
                _metric("peaks.count", value=Decimal("20")),
            ),
            key=lambda metric: metric.metric_id,
        )
    )
    _replace_qc(setup, original, artifacts=artifacts)
    original_generation = setup.get_result_state("run-1").qc_generation
    assert original_generation is not None
    setup_engine.dispose()

    reader_engine = create_database_engine(database_url)
    writer_engine = create_database_engine(database_url)
    reader = SqlAlchemyRunRepository(create_session_factory(reader_engine))
    writer = SqlAlchemyRunRepository(create_session_factory(writer_engine))
    replacement = tuple(
        replace(metric, value=metric.value + Decimal("1")) for metric in original
    )
    writer_attempt, artifact_generation = _begin_qc_attempt(writer)
    state_selected = Event()
    replacement_committed = Event()

    def pause_after_state_select(
        _connection,
        _cursor,
        statement,
        _parameters,
        _context,
        _executemany,
    ):
        if "FROM run_result_states" not in statement or state_selected.is_set():
            return
        state_selected.set()
        assert replacement_committed.wait(timeout=5)

    sqlalchemy_event.listen(
        reader_engine,
        "after_cursor_execute",
        pause_after_state_select,
    )

    def read_second_page():
        return reader.list_qc_metrics_page(
            "run-1",
            expected_generation=original_generation,
            after=original[0].metric_id,
            limit=1,
        )

    def replace_generation():
        assert state_selected.wait(timeout=5)
        try:
            writer.replace_qc_metrics(
                "run-1",
                replacement,
                attempt_id=writer_attempt,
                expected_artifact_generation=artifact_generation,
                expected_artifacts=artifacts,
                expected_status=RunStatus.SUCCEEDED,
                event=_event("qc_metrics_indexed", context={"metric_count": 2}),
            )
        finally:
            replacement_committed.set()

    try:
        with ThreadPoolExecutor(max_workers=2) as pool:
            read_future = pool.submit(read_second_page)
            replace_future = pool.submit(replace_generation)
            page = read_future.result(timeout=10)
            replace_future.result(timeout=10)

        assert page == (
            original_generation,
            (original[1],),
        )
        final_generation, final_metrics = reader.list_qc_metrics_page(
            "run-1",
            expected_generation=None,
        )
        assert final_generation != original_generation
        assert final_metrics == replacement
    finally:
        sqlalchemy_event.remove(
            reader_engine,
            "after_cursor_execute",
            pause_after_state_select,
        )
        reader_engine.dispose()
        writer_engine.dispose()


def test_sqlalchemy_artifact_generation_bound_page_read_is_atomic_during_replacement(
    tmp_path,
):
    database_url = f"sqlite:///{tmp_path / 'platform.db'}"
    upgrade_database(database_url)
    setup_engine = create_database_engine(database_url)
    setup = SqlAlchemyRunRepository(create_session_factory(setup_engine))
    original = (
        _artifact("artifact-a", revision_seed="generation-a"),
        _artifact("artifact-b", revision_seed="generation-a"),
    )
    _prepare(setup, original)
    original_generation = setup.get_result_state("run-1").artifact_generation
    assert original_generation is not None
    setup_engine.dispose()

    reader_engine = create_database_engine(database_url)
    writer_engine = create_database_engine(database_url)
    reader = SqlAlchemyRunRepository(create_session_factory(reader_engine))
    writer = SqlAlchemyRunRepository(create_session_factory(writer_engine))
    replacement = (
        _artifact("artifact-a", revision_seed="generation-b"),
        _artifact("artifact-b", revision_seed="generation-b"),
    )
    writer_attempt = _begin_artifact_attempt(writer)
    state_selected = Event()
    replacement_committed = Event()

    def pause_after_state_select(
        _connection,
        _cursor,
        statement,
        _parameters,
        _context,
        _executemany,
    ):
        if "FROM run_result_states" not in statement or state_selected.is_set():
            return
        state_selected.set()
        assert replacement_committed.wait(timeout=5)

    sqlalchemy_event.listen(
        reader_engine,
        "after_cursor_execute",
        pause_after_state_select,
    )

    def read_second_page():
        return reader.list_artifacts_page(
            "run-1",
            expected_generation=original_generation,
            after="artifact-a",
            limit=1,
        )

    def replace_generation():
        assert state_selected.wait(timeout=5)
        try:
            writer.replace_artifacts(
                "run-1",
                replacement,
                attempt_id=writer_attempt,
                expected_status=RunStatus.SUCCEEDED,
                event=_event("artifacts_indexed", context={"artifact_count": 2}),
            )
        finally:
            replacement_committed.set()

    try:
        with ThreadPoolExecutor(max_workers=2) as pool:
            read_future = pool.submit(read_second_page)
            replace_future = pool.submit(replace_generation)
            page = read_future.result(timeout=10)
            replace_future.result(timeout=10)

        assert page == (
            original_generation,
            (original[1],),
        )
        final_generation, final_artifacts = reader.list_artifacts_page(
            "run-1",
            expected_generation=None,
        )
        assert final_generation != original_generation
        assert final_artifacts == replacement
    finally:
        sqlalchemy_event.remove(
            reader_engine,
            "after_cursor_execute",
            pause_after_state_select,
        )
        reader_engine.dispose()
        writer_engine.dispose()
