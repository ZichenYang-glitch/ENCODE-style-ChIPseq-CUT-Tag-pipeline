"""Shared HTTP pagination contract over real memory and file-backed SQL stores."""

from __future__ import annotations

import asyncio
from dataclasses import replace
from datetime import datetime, timezone
import json
from types import SimpleNamespace

import httpx
import pytest

from api_test_client import seeded_auth_cookies
from conftest import seed_test_authentication
from encode_pipeline.api.main import create_app
from encode_pipeline.persistence import open_run_persistence
from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.platform.runs import RunRecord, RunStatus
from encode_pipeline.services.run_repositories import (
    InMemoryRunRepository,
    RunEventDraft,
)


@pytest.fixture(params=("sqlite", "memory"))
def progress_api(request, tmp_path, monkeypatch):
    """Only repository construction differs; auth always uses a real test session."""
    database = tmp_path / "pagination.db"

    def opener(url):
        persistence = open_run_persistence(url)
        if request.param == "memory":
            persistence = replace(persistence, repository=InMemoryRunRepository())
        return persistence

    app = create_app(
        database_url=f"sqlite:///{database}",
        workspace_root=tmp_path / "workspaces",
        project_root=tmp_path,
        registry=WorkflowRegistry(),
        persistence_opener=opener,
    )
    seed_test_authentication(app)
    repository = app.state.persistence.repository
    calls = []
    responses = []

    def observe(name):
        original = getattr(repository, name)

        def invoke(*args, **kwargs):
            calls.append({"method": name, "args": args, "kwargs": kwargs})
            return original(*args, **kwargs)  # Observe, never replace the query result.

        monkeypatch.setattr(repository, name, invoke)

    observe("list_events")
    observe("list_logs")

    def seed(count=3, run_id="run-pagination"):
        now = datetime(2026, 9, 23, tzinfo=timezone.utc)
        repository.create_run(
            RunRecord(
                run_id=run_id,
                workflow_id="pagination-test",
                inputs={},
                status=RunStatus.CREATED,
                created_at=now,
                updated_at=now,
                started_at=None,
                ended_at=None,
                current_stage=None,
                cancellation_reason=None,
                error=None,
            ),
            RunEventDraft(event_type="created", message="event-1"),
        )
        for index in range(2, count + 1):
            repository.add_event(
                run_id, RunEventDraft(event_type="test", message=f"event-{index}")
            )
        for index in range(1, count + 1):
            for stream in ("stdout", "stderr"):
                repository.append_log(run_id, stream, [f"{stream}-{index}"])
        return run_id

    async def get_async(
        path, *, params=None, authenticated=True, raise_exceptions=False
    ):
        transport = httpx.ASGITransport(app=app, raise_app_exceptions=raise_exceptions)
        async with httpx.AsyncClient(
            transport=transport,
            base_url="http://testserver",
            cookies=seeded_auth_cookies(app) if authenticated else None,
        ) as client:
            response = await client.get(path, params=params)
        responses.append(
            {
                "path": path,
                "params": params,
                "authenticated": authenticated,
                "status": response.status_code,
                "body": response.json(),
                "queries": list(calls),
            }
        )
        return response

    def get(path, **kwargs):
        return asyncio.run(get_async(path, **kwargs))

    try:
        yield SimpleNamespace(
            app=app, repository=repository, seed=seed, get=get, calls=calls
        )
    finally:
        (tmp_path / "http-evidence.json").write_text(
            json.dumps(
                {
                    "backend": request.param,
                    "database": str(database),
                    "responses": responses,
                },
                indent=2,
            )
            + "\n"
        )
        app.state.run_queue.close()
        try:
            assert app.state.persistence.engine.pool.checkedout() == 0
        finally:
            app.state.persistence.close()


@pytest.fixture(params=("events", "logs"))
def endpoint(request):
    return request.param


def items(body, endpoint):
    return body["events" if endpoint == "events" else "chunks"]


def identifier(item, endpoint):
    return item["event_id" if endpoint == "events" else "chunk_id"]


@pytest.mark.parametrize(
    "limit", [None, 1, 50, 100, 101, 2**63 - 2, 2**63 - 1, 2**63, 2**100, 0, -1, "text"]
)
def test_http_limit_bounds(progress_api, endpoint, limit):
    run_id = progress_api.seed()
    params = {} if limit is None else {"limit": str(limit)}
    response = progress_api.get(f"/api/v1/runs/{run_id}/{endpoint}", params=params)
    valid = limit is None or isinstance(limit, int) and 1 <= limit <= 100
    assert response.status_code == (200 if valid else 400)
    body = response.json()
    if valid:
        size = 50 if limit is None else limit
        assert body["ok"] is True and body["issues"] == []
        rows = items(body, endpoint)
        assert [row["sequence"] for row in rows] == list(range(1, min(size, 3) + 1))
        assert body["next_cursor"] == (
            identifier(rows[-1], endpoint) if size < 3 else None
        )
        assert progress_api.calls == [
            {
                "method": f"list_{endpoint}",
                "args": (run_id,) if endpoint == "events" else (run_id, "stdout"),
                "kwargs": {"after": None, "limit": size + 1},
            }
        ]
    else:
        assert body["ok"] is False
        assert [issue["code"] for issue in body["issues"]] == ["API_REQUEST_INVALID"]
        assert body["issues"][0]["technical_message"] is None
        assert progress_api.calls == []


@pytest.mark.parametrize("count", (99, 100, 101))
def test_page_lookahead_cursor_order_and_final_page(progress_api, endpoint, count):
    run_id = progress_api.seed(count)
    # Logs exercise both streams independently through the very same route/query.
    for stream in ("stdout", "stderr") if endpoint == "logs" else (None,):
        params = {"limit": 100}
        if stream:
            params["stream_name"] = stream
        rows = []
        cursors = []
        while True:
            response = progress_api.get(
                f"/api/v1/runs/{run_id}/{endpoint}", params=params
            )
            assert response.status_code == 200
            body = response.json()
            page = items(body, endpoint)
            assert len(page) == min(100, count - len(rows))
            rows.extend(page)
            cursor = body["next_cursor"]
            if len(rows) == count:
                assert cursor is None
                break
            assert cursor == identifier(page[-1], endpoint)
            assert cursor not in cursors
            cursors.append(cursor)
            params["after"] = cursor
        assert [row["sequence"] for row in rows] == list(range(1, count + 1))
        assert len({identifier(row, endpoint) for row in rows}) == count
        if stream:
            assert all(row["stream_name"] == stream for row in rows)
            assert [row["lines"] for row in rows] == [
                [f"{stream}-{i}"] for i in range(1, count + 1)
            ]
        else:
            assert [row["message"] for row in rows] == [
                f"event-{i}" for i in range(1, count + 1)
            ]
        params["after"] = identifier(rows[-1], endpoint)
        empty = progress_api.get(f"/api/v1/runs/{run_id}/{endpoint}", params=params)
        assert empty.status_code == 200
        assert items(empty.json(), endpoint) == []
        assert empty.json()["next_cursor"] is None
    assert all(call["kwargs"]["limit"] == 101 for call in progress_api.calls)


def test_missing_run_cursor_and_authentication_keep_existing_responses(
    progress_api, endpoint
):
    run_id = progress_api.seed()
    missing = progress_api.get(f"/api/v1/runs/missing/{endpoint}")
    assert missing.status_code == 404
    assert missing.json()["issues"][0]["code"] == "RUN_NOT_FOUND"
    invalid = progress_api.get(
        f"/api/v1/runs/{run_id}/{endpoint}", params={"after": "absent"}
    )
    assert invalid.status_code == 400
    assert invalid.json()["issues"][0]["code"] == "RUN_CURSOR_NOT_FOUND"
    progress_api.calls.clear()
    anonymous = progress_api.get(
        f"/api/v1/runs/{run_id}/{endpoint}", authenticated=False
    )
    assert anonymous.status_code == 401
    assert progress_api.calls == []


def test_log_cursor_cannot_cross_stream(progress_api):
    run_id = progress_api.seed(3)
    chunk = progress_api.repository.append_log(run_id, "stderr", ["stderr-only"])
    response = progress_api.get(
        f"/api/v1/runs/{run_id}/logs",
        params={"stream_name": "stdout", "after": chunk.chunk_id},
    )
    assert response.status_code == 400
    assert response.json()["issues"][0]["code"] == "RUN_CURSOR_NOT_FOUND"
    # No logs in an otherwise valid stream remains an empty, successful page.
    empty = progress_api.get(
        f"/api/v1/runs/{run_id}/logs", params={"stream_name": "empty"}
    )
    assert empty.status_code == 200
    assert empty.json()["chunks"] == [] and empty.json()["next_cursor"] is None


def test_openapi_progress_limits(progress_api, endpoint):
    operation = progress_api.app.openapi()["paths"][
        f"/api/v1/runs/{{run_id}}/{endpoint}"
    ]["get"]
    schema = next(p["schema"] for p in operation["parameters"] if p["name"] == "limit")
    assert schema["type"] == "integer"
    assert {key: schema.get(key) for key in ("default", "minimum", "maximum")} == {
        "default": 50,
        "minimum": 1,
        "maximum": 100,
    }
