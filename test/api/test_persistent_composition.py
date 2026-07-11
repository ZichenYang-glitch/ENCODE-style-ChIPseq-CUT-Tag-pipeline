"""API-level durable SQLite composition and restart recovery tests."""

from __future__ import annotations

import asyncio

import pytest

from api_test_client import ApiTestClient
from encode_pipeline.api.main import create_app
from encode_pipeline.persistence import RunPersistence, SqlAlchemyRunRepository
from encode_pipeline.platform.adapters import WorkflowInputs
from encode_pipeline.platform.runs import RunStatus


WORKFLOW_ID = "encode-style-chipseq-cuttag-atac-mnase"


def _database_url(tmp_path) -> str:
    return f"sqlite:///{tmp_path / 'platform.db'}"


def _create_run(client: ApiTestClient) -> str:
    response = client.post(
        f"/api/v1/workflows/{WORKFLOW_ID}/runs",
        json={"config": {"samples": "samples.tsv"}},
    )
    assert response.status_code == 201
    return response.json()["run"]["run_id"]


def test_default_api_composition_uses_environment_configured_sqlite(
    tmp_path, monkeypatch
):
    database_url = _database_url(tmp_path)
    monkeypatch.setenv("ENCODE_PIPELINE_DATABASE_URL", database_url)

    app = create_app()

    assert app.state.database_url == database_url
    assert isinstance(app.state.persistence.repository, SqlAlchemyRunRepository)
    app.state.persistence.close()


def test_api_lifespan_closes_owned_persistence(tmp_path, monkeypatch):
    app = create_app(database_url=_database_url(tmp_path))
    closed_urls: list[str] = []
    original_close = RunPersistence.close

    def tracked_close(persistence: RunPersistence) -> None:
        closed_urls.append(persistence.database_url)
        original_close(persistence)

    monkeypatch.setattr(RunPersistence, "close", tracked_close)

    async def exercise_lifespan() -> None:
        async with app.router.lifespan_context(app):
            pass

    asyncio.run(exercise_lifespan())

    assert closed_urls == [_database_url(tmp_path)]


def test_api_run_is_visible_after_app_factory_restart(tmp_path):
    database_url = _database_url(tmp_path)
    first_app = create_app(database_url=database_url)
    with ApiTestClient(first_app) as client:
        run_id = _create_run(client)
    first_app.state.run_service.add_event(run_id, "checkpoint", "Before restart.")
    first_app.state.run_service.append_log(run_id, "stdout", ["before restart"])
    first_app.state.persistence.close()

    second_app = create_app(database_url=database_url)
    with ApiTestClient(second_app) as client:
        response = client.get(f"/api/v1/runs/{run_id}")
        events_response = client.get(f"/api/v1/runs/{run_id}/events")
        logs_response = client.get(f"/api/v1/runs/{run_id}/logs")

    assert response.status_code == 200
    body = response.json()
    assert body["ok"] is True
    assert body["run"]["run_id"] == run_id
    assert body["run"]["status"] == "created"
    assert events_response.status_code == 200
    assert [event["event_type"] for event in events_response.json()["events"]] == [
        "status_changed",
        "checkpoint",
    ]
    assert logs_response.status_code == 200
    assert logs_response.json()["chunks"][0]["lines"] == ["before restart"]
    second_app.state.persistence.close()


@pytest.mark.parametrize(
    ("active_status", "transitions"),
    [
        (RunStatus.VALIDATING, (RunStatus.VALIDATING,)),
        (
            RunStatus.QUEUED,
            (RunStatus.VALIDATING, RunStatus.PLANNED, RunStatus.QUEUED),
        ),
        (
            RunStatus.RUNNING,
            (
                RunStatus.VALIDATING,
                RunStatus.PLANNED,
                RunStatus.QUEUED,
                RunStatus.RUNNING,
            ),
        ),
    ],
)
def test_api_restart_recovers_active_run_with_public_safe_failure(
    tmp_path,
    active_status: RunStatus,
    transitions: tuple[RunStatus, ...],
):
    database_url = _database_url(tmp_path)
    first_app = create_app(database_url=database_url)
    created = first_app.state.run_service.create_run(
        WORKFLOW_ID,
        WorkflowInputs(config={"samples": "private/input.tsv"}),
    )
    for status in transitions:
        first_app.state.run_service.transition_run(
            created.run_id,
            status,
            stage="preflight" if status is RunStatus.VALIDATING else None,
        )
    first_app.state.persistence.close()

    restarted_app = create_app(database_url=database_url)
    with ApiTestClient(restarted_app) as client:
        run_response = client.get(f"/api/v1/runs/{created.run_id}")
        events_response = client.get(f"/api/v1/runs/{created.run_id}/events")

    assert restarted_app.state.recovered_run_ids == (created.run_id,)
    assert run_response.status_code == 200
    run = run_response.json()["run"]
    assert run["status"] == "failed"
    assert run["error"]["code"] == "RUN_INTERRUPTED_BY_API_RESTART"
    assert "private/input.tsv" not in str(run["error"])
    assert events_response.status_code == 200
    recovery_event = events_response.json()["events"][-1]
    assert recovery_event["event_type"] == "run_recovered_after_restart"
    assert recovery_event["status"] == "failed"
    assert recovery_event["context"] == {
        "previous_status": active_status.value,
        "new_status": "failed",
        "reason_code": "API_RESTART_INTERRUPTED",
    }
    assert "private/input.tsv" not in str(recovery_event)
    restarted_app.state.persistence.close()


def test_api_restart_preserves_planned_run_for_future_execution(tmp_path):
    database_url = _database_url(tmp_path)
    first_app = create_app(database_url=database_url)
    created = first_app.state.run_service.create_run(
        WORKFLOW_ID,
        WorkflowInputs(config={"samples": "samples.tsv"}),
    )
    first_app.state.run_service.transition_run(created.run_id, RunStatus.VALIDATING)
    planned = first_app.state.run_service.transition_run(
        created.run_id, RunStatus.PLANNED
    )
    first_app.state.persistence.close()

    restarted_app = create_app(database_url=database_url)
    with ApiTestClient(restarted_app) as client:
        response = client.get(f"/api/v1/runs/{created.run_id}")

    assert restarted_app.state.recovered_run_ids == ()
    assert response.status_code == 200
    assert response.json()["run"]["status"] == planned.status.value
    restarted_app.state.persistence.close()
