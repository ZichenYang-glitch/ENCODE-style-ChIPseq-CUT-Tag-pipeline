"""Read-only run-history API contract and disclosure tests."""

from __future__ import annotations

from datetime import datetime, timedelta, timezone
import asyncio
import inspect as python_inspect
from threading import get_ident

import pytest
import httpx
from sqlalchemy import delete, update

from api_test_client import ApiTestClient
from encode_pipeline.api.main import create_app
from encode_pipeline.api.routes.runs import list_runs
from encode_pipeline.persistence.models import RunRow
from encode_pipeline.platform.runs import RunRecord, RunStatus
from encode_pipeline.services.run_repositories import RunEventDraft


NOW = datetime(2026, 7, 14, 8, 0, tzinfo=timezone.utc)


@pytest.fixture
def client():
    app = create_app()
    with ApiTestClient(app) as test_client:
        yield test_client


def _record(
    run_id: str,
    *,
    workflow_id: str = "workflow-a",
    status: RunStatus = RunStatus.SUCCEEDED,
    created_at: datetime = NOW,
) -> RunRecord:
    ended_at = created_at + timedelta(minutes=1) if status.is_terminal else None
    started_at = (
        created_at + timedelta(seconds=1) if status is not RunStatus.CREATED else None
    )
    return RunRecord(
        run_id=run_id,
        workflow_id=workflow_id,
        inputs={"workspace": "/private/workspace", "token": "SECRET"},
        status=status,
        created_at=created_at,
        updated_at=ended_at or started_at or created_at,
        started_at=started_at,
        ended_at=ended_at,
        current_stage="execution" if started_at is not None else None,
        cancellation_reason="private cancellation"
        if status is RunStatus.CANCELLED
        else None,
        error=None,
        tags={"owner": "private"},
    )


def _create(client: ApiTestClient, record: RunRecord) -> None:
    client.app.state.persistence.repository.create_run(
        record,
        RunEventDraft(
            event_type="status_changed",
            message="Run created.",
            status=record.status,
        ),
    )


def test_list_runs_is_sync_with_stable_operation_id_and_empty_envelope(client) -> None:
    assert not python_inspect.iscoroutinefunction(list_runs)
    operation = client.get("/openapi.json").json()["paths"]["/api/v1/runs"]["get"]
    assert operation["operationId"] == "listRuns"

    response = client.get("/api/v1/runs")

    assert response.status_code == 200
    assert response.json() == {
        "ok": True,
        "runs": [],
        "next_cursor": None,
        "issues": [],
    }


def test_list_runs_projects_only_safe_summary_fields_in_stable_order(client) -> None:
    _create(client, _record("run-a"))
    _create(client, _record("run-b"))

    response = client.get("/api/v1/runs")

    assert response.status_code == 200
    body = response.json()
    assert [run["run_id"] for run in body["runs"]] == ["run-b", "run-a"]
    assert set(body["runs"][0]) == {
        "run_id",
        "workflow_id",
        "status",
        "created_at",
        "updated_at",
        "started_at",
        "ended_at",
        "current_stage",
    }
    assert "/private" not in response.text
    assert "SECRET" not in response.text
    assert "technical_message" not in response.text


def test_list_runs_pages_and_binds_cursor_to_filters(client) -> None:
    _create(client, _record("run-a", workflow_id="workflow-a"))
    _create(client, _record("run-c", workflow_id="workflow-a"))
    _create(
        client,
        _record("run-b", workflow_id="workflow-b", status=RunStatus.FAILED),
    )

    first = client.get(
        "/api/v1/runs",
        params={"limit": 1, "workflow_id": "workflow-a"},
    )
    assert first.status_code == 200
    assert [run["run_id"] for run in first.json()["runs"]] == ["run-c"]
    cursor = first.json()["next_cursor"]
    assert cursor.startswith("runhist_")

    second = client.get(
        "/api/v1/runs",
        params={
            "limit": 1,
            "workflow_id": "workflow-a",
            "after": cursor,
        },
    )
    assert second.status_code == 200
    assert [run["run_id"] for run in second.json()["runs"]] == ["run-a"]
    assert second.json()["next_cursor"] is None

    cross_filter = client.get(
        "/api/v1/runs",
        params={"workflow_id": "workflow-b", "after": cursor},
    )
    assert cross_filter.status_code == 400
    assert cross_filter.json()["issues"][0]["code"] == ("RUN_HISTORY_CURSOR_INVALID")
    assert cursor not in cross_filter.text


@pytest.mark.parametrize(
    "params",
    [
        {"limit": 0},
        {"limit": 101},
        {"status": "unknown"},
        {"workflow_id": "x" * 256},
        {"after": "not-a-cursor"},
    ],
)
def test_invalid_query_uses_stable_history_envelope(client, params) -> None:
    response = client.get("/api/v1/runs", params=params)

    assert response.status_code == 400
    assert response.json()["ok"] is False
    assert response.json()["runs"] == []
    assert response.json()["next_cursor"] is None
    assert response.json()["issues"][0]["code"] == "API_REQUEST_INVALID"
    assert "technical_message" not in response.text


def test_deleted_cursor_boundary_returns_safe_not_found(client) -> None:
    _create(client, _record("run-a"))
    _create(client, _record("run-b"))
    first = client.get("/api/v1/runs", params={"limit": 1})
    cursor = first.json()["next_cursor"]
    with client.app.state.persistence.engine.begin() as connection:
        connection.execute(delete(RunRow).where(RunRow.run_id == "run-b"))

    response = client.get(
        "/api/v1/runs",
        params={"limit": 1, "after": cursor},
    )

    assert response.status_code == 400
    assert response.json()["issues"][0]["code"] == ("RUN_HISTORY_CURSOR_NOT_FOUND")
    assert cursor not in response.text


def test_corrupt_persisted_summary_fails_closed_without_sensitive_text(client) -> None:
    _create(client, _record("run-a"))
    private_value = "/private/workspace/stage"
    with client.app.state.persistence.engine.begin() as connection:
        connection.execute(
            update(RunRow)
            .where(RunRow.run_id == "run-a")
            .values(current_stage=private_value)
        )

    response = client.get("/api/v1/runs")

    assert response.status_code == 500
    assert response.json()["runs"] == []
    assert response.json()["issues"][0]["code"] == "RUN_HISTORY_DATA_INVALID"
    assert private_value not in response.text
    assert "technical_message" not in response.text


@pytest.mark.parametrize(
    "changes",
    [
        {"status": "succeeded", "started_at": NOW, "ended_at": None},
        {"status": "running", "started_at": None, "ended_at": None},
        {"status": "planned", "started_at": NOW, "ended_at": None},
        {"status": "running", "started_at": NOW, "ended_at": NOW},
    ],
)
def test_lifecycle_corruption_fails_closed_as_invalid_history_data(
    client,
    changes,
) -> None:
    _create(client, _record("run-a", status=RunStatus.CREATED))
    with client.app.state.persistence.engine.begin() as connection:
        connection.execute(
            update(RunRow).where(RunRow.run_id == "run-a").values(**changes)
        )

    response = client.get("/api/v1/runs")

    assert response.status_code == 500
    assert response.json()["runs"] == []
    assert response.json()["issues"][0]["code"] == "RUN_HISTORY_DATA_INVALID"
    assert "technical_message" not in response.text


def test_list_runs_repository_work_executes_off_request_event_loop(
    client, monkeypatch
) -> None:
    request_thread = get_ident()
    called_threads: list[int] = []
    repository = client.app.state.persistence.repository
    original = repository.list_run_summaries

    def observed(**kwargs):
        called_threads.append(get_ident())
        return original(**kwargs)

    monkeypatch.setattr(repository, "list_run_summaries", observed)

    response = client.get("/api/v1/runs")

    assert response.status_code == 200
    assert called_threads and called_threads[0] != request_thread


def test_unexpected_failure_uses_history_specific_safe_envelope(
    client, monkeypatch
) -> None:
    private_text = "sqlite:///private/history.db SECRET_TOKEN=value"

    def fail(**_kwargs):
        raise RuntimeError(private_text)

    monkeypatch.setattr(client.app.state.run_service, "list_run_history", fail)

    async def request() -> httpx.Response:
        transport = httpx.ASGITransport(
            app=client.app,
            raise_app_exceptions=False,
        )
        async with httpx.AsyncClient(
            transport=transport,
            base_url="http://testserver",
        ) as async_client:
            return await async_client.get("/api/v1/runs")

    response = asyncio.run(request())

    assert response.status_code == 500
    assert response.json()["ok"] is False
    assert response.json()["runs"] == []
    assert response.json()["next_cursor"] is None
    assert response.json()["issues"][0]["code"] == "INTERNAL_SERVER_ERROR"
    assert private_text not in response.text
    assert "technical_message" not in response.text
