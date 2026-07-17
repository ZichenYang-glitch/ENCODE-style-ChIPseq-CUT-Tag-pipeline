"""Tests for the read-only persisted QC metrics API."""

from __future__ import annotations

import asyncio
from collections.abc import Iterator
from dataclasses import replace
from datetime import datetime, timezone
from decimal import Decimal
import inspect
import re
from threading import Event, Thread

import fastapi.routing
import httpx
import pytest
from sqlalchemy import update

from encode_pipeline.api.main import create_app
from encode_pipeline.api.models import QcMetricResponse
from encode_pipeline.api.routes.qc_metrics import list_run_qc_metrics
from encode_pipeline.persistence.models import RunQcMetricRow
from encode_pipeline.platform.adapters import WorkflowInputs
from encode_pipeline.platform.result_generations import (
    QcMetricCursor,
    encode_qc_metric_cursor,
)
from encode_pipeline.platform.runs import (
    RunArtifactRef,
    RunQcMetric,
    RunStatus,
    build_qc_metric_id,
)
from api_test_client import ApiTestClient


WORKFLOW_ID = "encode-style-chipseq-cuttag-atac-mnase"
PRODUCED_AT = datetime(2026, 7, 13, 8, 30, tzinfo=timezone.utc)
QC_GENERATION_PATTERN = re.compile(r"^qcgen-[0-9a-f]{64}$")
ARTIFACT_REVISION = f"artifactrev-{'a' * 64}"


async def _run_in_joined_test_thread(function, *args, **kwargs):
    """Exercise sync FastAPI routes without leaking Python 3.13 threads."""
    completed = Event()
    results: list[object] = []
    exceptions: list[BaseException] = []

    def invoke() -> None:
        try:
            results.append(function(*args, **kwargs))
        except BaseException as exc:
            exceptions.append(exc)
        finally:
            completed.set()

    thread = Thread(target=invoke)
    thread.start()
    try:
        while not completed.is_set():
            await asyncio.sleep(0.001)
    finally:
        thread.join(timeout=3)
        if thread.is_alive():  # pragma: no cover - test deadlock guard
            raise RuntimeError("test threadpool call did not terminate")
    if exceptions:
        raise exceptions[0]
    return results[0]


@pytest.fixture(autouse=True)
def joined_test_threadpool(monkeypatch):
    monkeypatch.setattr(
        fastapi.routing,
        "run_in_threadpool",
        _run_in_joined_test_thread,
    )


@pytest.fixture
def client(tmp_path) -> Iterator[ApiTestClient]:
    app = create_app(
        database_url=f"sqlite:///{tmp_path / 'platform.db'}",
        workspace_root=tmp_path / "workspaces",
    )
    with ApiTestClient(app) as test_client:
        yield test_client


def _create_succeeded_run(client: ApiTestClient) -> str:
    record = client.app.state.run_service.create_run(
        WORKFLOW_ID,
        WorkflowInputs(config={}),
    )
    run_id = record.run_id
    for status in (
        RunStatus.VALIDATING,
        RunStatus.PLANNED,
        RunStatus.QUEUED,
        RunStatus.RUNNING,
        RunStatus.SUCCEEDED,
    ):
        client.app.state.run_service.transition_run(run_id, status)
    return run_id


def _artifact(run_id: str) -> RunArtifactRef:
    return RunArtifactRef(
        artifact_id="artifact-qc-summary",
        run_id=run_id,
        artifact_type="file",
        name="qc_summary.tsv",
        uri=f"run://runs/{run_id}/artifacts/artifact-qc-summary",
        mime_type="text/tab-separated-values",
        produced_at=PRODUCED_AT,
        revision=ARTIFACT_REVISION,
        metadata={
            "relative_path": "results/multiqc/qc_summary.tsv",
            "output_type": "qc_summary",
            "size_bytes": 12,
        },
    )


def _metric(
    run_id: str,
    metric_key: str = "sequencing.total_reads",
    *,
    value: Decimal = Decimal("9007199254740993.123456789012"),
    scope: str = "sample",
    sample_id: str | None = "-S1",
    experiment_id: str | None = "EXP1",
    assay: str | None = "chipseq",
    qc_flag: str | None = None,
    unit: str = "count",
) -> RunQcMetric:
    return RunQcMetric(
        metric_id=build_qc_metric_id(
            metric_key,
            scope,
            sample_id,
            experiment_id,
        ),
        run_id=run_id,
        metric_key=metric_key,
        display_name=metric_key.replace(".", " ").title(),
        value=value,
        unit=unit,
        scope=scope,
        sample_id=sample_id,
        experiment_id=experiment_id,
        assay=assay,
        qc_flag=qc_flag,
        source_artifact_id="artifact-qc-summary",
        produced_at=PRODUCED_AT,
    )


def _record_metrics(
    client: ApiTestClient,
    run_id: str,
    metrics: tuple[RunQcMetric, ...],
) -> None:
    artifact = _artifact(run_id)
    client.app.state.run_service.replace_artifacts(run_id, (artifact,))
    client.app.state.run_service.replace_qc_metrics(
        run_id,
        metrics,
        expected_artifacts=(artifact,),
    )


def test_qc_metric_route_uses_fastapi_threadpool_for_sync_database_io():
    assert not inspect.iscoroutinefunction(list_run_qc_metrics)


def test_existing_run_without_qc_metrics_returns_empty_page(client):
    run_id = _create_succeeded_run(client)

    response = client.get(f"/api/v1/runs/{run_id}/qc-metrics")

    assert response.status_code == 200
    assert response.json() == {
        "ok": True,
        "run_id": run_id,
        "qc_generation": None,
        "qc_metrics": [],
        "next_cursor": None,
        "issues": [],
    }


def test_qc_metrics_are_stably_paginated_and_decimal_is_lossless(client):
    run_id = _create_succeeded_run(client)
    metrics = (
        _metric(run_id, "sequencing.total_reads"),
        _metric(run_id, "peaks.count", value=Decimal("42")),
        _metric(run_id, "library.pbc2", value=Decimal("1.2500")),
    )
    _record_metrics(client, run_id, metrics)
    ordered = sorted(metrics, key=lambda metric: metric.metric_id)

    first = client.get(
        f"/api/v1/runs/{run_id}/qc-metrics",
        params={"limit": 2},
    )
    generation = first.json()["qc_generation"]
    cursor = first.json()["next_cursor"]
    second = client.get(
        f"/api/v1/runs/{run_id}/qc-metrics",
        params={"after": cursor, "generation": generation, "limit": 2},
    )

    assert first.status_code == second.status_code == 200
    assert [item["metric_id"] for item in first.json()["qc_metrics"]] == [
        item.metric_id for item in ordered[:2]
    ]
    assert cursor.startswith("qccur_")
    assert [item["metric_id"] for item in second.json()["qc_metrics"]] == [
        ordered[2].metric_id
    ]
    values = {
        item["metric_key"]: item["value"]
        for item in first.json()["qc_metrics"] + second.json()["qc_metrics"]
    }
    assert values == {
        "sequencing.total_reads": "9007199254740993.123456789012",
        "peaks.count": "42",
        "library.pbc2": "1.25",
    }


def test_qc_pages_expose_generation_and_use_an_opaque_generation_bound_cursor(client):
    run_id = _create_succeeded_run(client)
    metrics = (
        _metric(run_id, "sequencing.total_reads"),
        _metric(run_id, "peaks.count", value=Decimal("42")),
    )
    _record_metrics(client, run_id, metrics)

    first = client.get(
        f"/api/v1/runs/{run_id}/qc-metrics",
        params={"limit": 1},
    )
    generation = first.json()["qc_generation"]
    cursor = first.json()["next_cursor"]
    second = client.get(
        f"/api/v1/runs/{run_id}/qc-metrics",
        params={
            "after": cursor,
            "generation": generation,
            "limit": 1,
        },
    )

    assert first.status_code == second.status_code == 200
    assert QC_GENERATION_PATTERN.fullmatch(generation)
    assert cursor.startswith("qccur_")
    assert cursor not in {metric.metric_id for metric in metrics}
    assert second.json()["qc_generation"] == generation


def test_qc_cursor_rejects_same_ids_and_count_from_a_new_value_generation(client):
    run_id = _create_succeeded_run(client)
    artifact = _artifact(run_id)
    first_generation = (
        _metric(run_id, "sequencing.total_reads", value=Decimal("100")),
        _metric(run_id, "peaks.count", value=Decimal("10")),
    )
    client.app.state.run_service.replace_artifacts(run_id, (artifact,))
    client.app.state.run_service.replace_qc_metrics(
        run_id,
        first_generation,
        expected_artifacts=(artifact,),
    )
    first_page = client.get(
        f"/api/v1/runs/{run_id}/qc-metrics",
        params={"limit": 1},
    )
    old_cursor = first_page.json()["next_cursor"]
    old_generation = first_page.json()["qc_generation"]

    second_generation = tuple(
        replace(metric, value=metric.value + 1) for metric in first_generation
    )
    client.app.state.run_service.replace_qc_metrics(
        run_id,
        second_generation,
        expected_artifacts=(artifact,),
    )
    stale = client.get(
        f"/api/v1/runs/{run_id}/qc-metrics",
        params={
            "after": old_cursor,
            "generation": old_generation,
            "limit": 1,
        },
    )
    current = client.get(f"/api/v1/runs/{run_id}/qc-metrics")

    assert stale.status_code == 409
    assert stale.json()["issues"][0]["code"] == ("RUN_QC_METRIC_GENERATION_CHANGED")
    assert stale.json()["qc_generation"] is None
    assert stale.json()["qc_metrics"] == []
    assert current.status_code == 200
    assert current.json()["qc_generation"] != old_generation
    assert {item["value"] for item in current.json()["qc_metrics"]} == {"11", "101"}


def test_qc_metrics_default_page_is_bounded_to_fifty(client):
    run_id = _create_succeeded_run(client)
    metrics = tuple(_metric(run_id, f"metric{i}") for i in range(51))
    _record_metrics(client, run_id, metrics)

    response = client.get(f"/api/v1/runs/{run_id}/qc-metrics")

    assert response.status_code == 200
    assert len(response.json()["qc_metrics"]) == 50
    assert response.json()["next_cursor"].startswith("qccur_")


def test_qc_metric_projection_preserves_nullable_fields(client):
    run_id = _create_succeeded_run(client)
    metric = _metric(
        run_id,
        "run.elapsed_seconds",
        value=Decimal("0"),
        scope="run",
        sample_id=None,
        experiment_id=None,
        assay=None,
    )
    _record_metrics(client, run_id, (metric,))

    response = client.get(f"/api/v1/runs/{run_id}/qc-metrics")

    assert response.status_code == 200
    assert response.json()["qc_metrics"] == [
        {
            "metric_id": metric.metric_id,
            "metric_key": "run.elapsed_seconds",
            "display_name": "Run Elapsed_Seconds",
            "value": "0",
            "unit": "count",
            "scope": "run",
            "sample_id": None,
            "experiment_id": None,
            "assay": None,
            "qc_flag": None,
            "source_artifact_id": "artifact-qc-summary",
            "produced_at": "2026-07-13T08:30:00Z",
        }
    ]
    assert "run_id" not in response.json()["qc_metrics"][0]


def test_qc_metric_projection_preserves_a_workflow_neutral_score_unit(client):
    run_id = _create_succeeded_run(client)
    metric = _metric(
        run_id,
        "rseqc.tin.mean_score",
        value=Decimal("72.125"),
        unit="score",
    )
    _record_metrics(client, run_id, (metric,))

    response = client.get(f"/api/v1/runs/{run_id}/qc-metrics")

    assert response.status_code == 200
    assert response.json()["qc_metrics"][0]["unit"] == "score"
    assert response.json()["qc_metrics"][0]["value"] == "72.125"


def test_unknown_run_returns_stable_redacted_404(client):
    response = client.get("/api/v1/runs/run-missing/qc-metrics")

    assert response.status_code == 404
    assert response.json()["qc_generation"] is None
    assert response.json()["issues"][0]["code"] == "RUN_NOT_FOUND"
    assert response.json()["issues"][0]["context"] == {"run_id": "run-missing"}
    assert "technical_message" not in response.text


@pytest.mark.parametrize(
    "change",
    (
        {"metric_id": build_qc_metric_id("peaks.count", "sample", "-S1", "EXP1")},
        {"display_name": "DATABASE_URL=sqlite"},
        {"source_artifact_id": "../private"},
    ),
)
def test_qc_metric_projection_rejects_invalid_or_disclosure_unsafe_values(change):
    metric = replace(_metric("run-1"), **change)

    with pytest.raises(ValueError):
        QcMetricResponse.from_metric(metric, expected_run_id="run-1")


def test_qc_metric_domain_rejects_naive_produced_time_before_projection():
    with pytest.raises(ValueError, match="timezone-aware"):
        replace(_metric("run-1"), produced_at=datetime(2026, 7, 13, 8, 30))


def test_qc_metric_projection_rejects_cross_run_domain_value():
    with pytest.raises(ValueError):
        QcMetricResponse.from_metric(
            _metric("run-other"),
            expected_run_id="run-1",
        )


def test_unknown_and_cross_run_cursors_share_the_same_400(client):
    first_run = _create_succeeded_run(client)
    second_run = _create_succeeded_run(client)
    current = _metric(first_run, "library.nrf")
    other = _metric(second_run, "library.nrf")
    _record_metrics(client, first_run, (current,))
    _record_metrics(client, second_run, (other,))
    first_generation = client.get(f"/api/v1/runs/{first_run}/qc-metrics").json()[
        "qc_generation"
    ]
    second_generation = client.get(f"/api/v1/runs/{second_run}/qc-metrics").json()[
        "qc_generation"
    ]
    missing_metric_id = build_qc_metric_id(
        "library.pbc1",
        "sample",
        "-S1",
        "EXP1",
    )
    missing_cursor = encode_qc_metric_cursor(
        QcMetricCursor(
            run_id=first_run,
            qc_generation=first_generation,
            after_metric_id=missing_metric_id,
        )
    )
    cross_run_cursor = encode_qc_metric_cursor(
        QcMetricCursor(
            run_id=second_run,
            qc_generation=second_generation,
            after_metric_id=other.metric_id,
        )
    )

    missing = client.get(
        f"/api/v1/runs/{first_run}/qc-metrics",
        params={"after": missing_cursor, "generation": first_generation},
    )
    cross_run = client.get(
        f"/api/v1/runs/{first_run}/qc-metrics",
        params={"after": cross_run_cursor, "generation": first_generation},
    )

    assert missing.status_code == cross_run.status_code == 400
    assert missing.json()["issues"][0]["code"] == ("RUN_QC_METRIC_CURSOR_NOT_FOUND")
    assert cross_run.json()["issues"][0]["code"] == ("RUN_QC_METRIC_CURSOR_NOT_FOUND")
    assert other.metric_id not in cross_run.text


@pytest.mark.parametrize("params", ({"limit": 0}, {"limit": 101}, {"after": "bad"}))
def test_invalid_query_parameters_use_qc_specific_envelope(client, params):
    run_id = _create_succeeded_run(client)

    response = client.get(
        f"/api/v1/runs/{run_id}/qc-metrics",
        params=params,
    )

    assert response.status_code == 400
    assert response.json() == {
        "ok": False,
        "run_id": run_id,
        "qc_generation": None,
        "qc_metrics": [],
        "issues": [
            {
                "code": "API_REQUEST_INVALID",
                "message": "QC metric query parameters are invalid.",
                "severity": "error",
                "path": "query",
                "source": "api",
                "hint": (
                    "Use a limit between 1 and 100 and a cursor returned "
                    "by this endpoint."
                ),
                "context": {},
            }
        ],
    }
    assert "technical_message" not in response.text


def test_damaged_persisted_metric_fails_closed_without_disclosure(client):
    run_id = _create_succeeded_run(client)
    metric = _metric(run_id)
    _record_metrics(client, run_id, (metric,))
    private_value = "WORKSPACE_SECRET=hidden"
    with client.app.state.persistence.engine.begin() as connection:
        connection.execute(
            update(RunQcMetricRow)
            .where(
                RunQcMetricRow.run_id == run_id,
                RunQcMetricRow.metric_id == metric.metric_id,
            )
            .values(display_name=private_value)
        )

    response = client.get(f"/api/v1/runs/{run_id}/qc-metrics")

    assert response.status_code == 500
    assert response.json()["qc_generation"] is None
    assert response.json()["qc_metrics"] == []
    assert response.json()["issues"][0]["code"] == "RUN_QC_METRIC_DATA_INVALID"
    assert private_value not in response.text
    assert "technical_message" not in response.text


def test_unexpected_qc_repository_failure_uses_declared_safe_envelope(
    client,
    monkeypatch,
):
    run_id = _create_succeeded_run(client)
    private_text = "sqlite:///private/qc.db SECRET_TOKEN=value"

    def fail(*_args, **_kwargs):
        raise RuntimeError(private_text)

    monkeypatch.setattr(client.app.state.run_service, "list_qc_metrics_page", fail)

    async def request() -> httpx.Response:
        transport = httpx.ASGITransport(
            app=client.app,
            raise_app_exceptions=False,
        )
        async with httpx.AsyncClient(
            transport=transport,
            base_url="http://testserver",
        ) as async_client:
            return await async_client.get(f"/api/v1/runs/{run_id}/qc-metrics")

    response = asyncio.run(request())

    assert response.status_code == 500
    assert response.json() == {
        "ok": False,
        "run_id": run_id,
        "qc_generation": None,
        "qc_metrics": [],
        "issues": [
            {
                "code": "INTERNAL_SERVER_ERROR",
                "message": "Run QC metrics are temporarily unavailable.",
                "severity": "error",
                "path": "qc_metrics",
                "source": "runtime",
                "context": {},
            }
        ],
    }
    assert private_text not in response.text
    assert "technical_message" not in response.text


def test_damaged_cursor_row_cannot_be_skipped(client):
    run_id = _create_succeeded_run(client)
    metrics = (
        _metric(run_id, "sequencing.total_reads"),
        _metric(run_id, "peaks.count"),
    )
    _record_metrics(client, run_id, metrics)
    cursor = min(metrics, key=lambda metric: metric.metric_id)
    generation = client.get(f"/api/v1/runs/{run_id}/qc-metrics").json()["qc_generation"]
    encoded_cursor = encode_qc_metric_cursor(
        QcMetricCursor(
            run_id=run_id,
            qc_generation=generation,
            after_metric_id=cursor.metric_id,
        )
    )
    with client.app.state.persistence.engine.begin() as connection:
        connection.execute(
            update(RunQcMetricRow)
            .where(
                RunQcMetricRow.run_id == run_id,
                RunQcMetricRow.metric_id == cursor.metric_id,
            )
            .values(display_name="Private/path")
        )

    response = client.get(
        f"/api/v1/runs/{run_id}/qc-metrics",
        params={"after": encoded_cursor, "generation": generation},
    )

    assert response.status_code == 500
    assert response.json()["qc_generation"] is None
    assert response.json()["issues"][0]["code"] == "RUN_QC_METRIC_DATA_INVALID"
    assert "Private/path" not in response.text


def test_qc_metrics_survive_sqlite_reopen(tmp_path):
    database_url = f"sqlite:///{tmp_path / 'platform.db'}"
    workspace_root = tmp_path / "workspaces"
    first_app = create_app(database_url=database_url, workspace_root=workspace_root)
    with ApiTestClient(first_app) as first_client:
        run_id = _create_succeeded_run(first_client)
        metric = _metric(run_id)
        _record_metrics(first_client, run_id, (metric,))
    first_app.state.persistence.close()

    second_app = create_app(database_url=database_url, workspace_root=workspace_root)
    with ApiTestClient(second_app) as second_client:
        response = second_client.get(f"/api/v1/runs/{run_id}/qc-metrics")
    second_app.state.persistence.close()

    assert response.status_code == 200
    assert response.json()["qc_metrics"][0]["metric_id"] == metric.metric_id
