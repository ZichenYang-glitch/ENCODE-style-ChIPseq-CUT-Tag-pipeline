"""API-level durable SQLite composition and restart recovery tests."""

from __future__ import annotations

import asyncio

import pytest

from api_test_client import ApiTestClient
from encode_pipeline.api.main import create_app
from encode_pipeline.persistence import RunPersistence, SqlAlchemyRunRepository
from encode_pipeline.platform.adapters import WorkflowInputs
from encode_pipeline.platform.runs import RunStatus
from encode_pipeline.workers.rq_queue import RqRunQueue


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
    closed_queues: list[str] = []
    original_close = RunPersistence.close
    original_queue_close = RqRunQueue.close

    def tracked_close(persistence: RunPersistence) -> None:
        closed_urls.append(persistence.database_url)
        original_close(persistence)

    def tracked_queue_close(run_queue: RqRunQueue) -> None:
        closed_queues.append(run_queue.queue_name)
        original_queue_close(run_queue)

    monkeypatch.setattr(RunPersistence, "close", tracked_close)
    monkeypatch.setattr(RqRunQueue, "close", tracked_queue_close)

    async def exercise_lifespan() -> None:
        async with app.router.lifespan_context(app):
            pass

    asyncio.run(exercise_lifespan())

    assert closed_urls == [_database_url(tmp_path)]
    assert closed_queues == [app.state.worker_settings.queue_name]


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


def test_api_restart_fails_api_owned_validation_with_public_safe_failure(tmp_path):
    database_url = _database_url(tmp_path)
    first_app = create_app(database_url=database_url)
    created = first_app.state.run_service.create_run(
        WORKFLOW_ID,
        WorkflowInputs(config={"samples": "private/input.tsv"}),
    )
    first_app.state.run_service.transition_run(
        created.run_id,
        RunStatus.VALIDATING,
        stage="preflight",
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
        "previous_status": RunStatus.VALIDATING.value,
        "new_status": "failed",
        "reason_code": "API_RESTART_INTERRUPTED",
    }
    assert "private/input.tsv" not in str(recovery_event)
    restarted_app.state.persistence.close()


@pytest.mark.parametrize("active_status", [RunStatus.QUEUED, RunStatus.RUNNING])
def test_api_restart_fails_unassigned_worker_state_as_orphan(
    tmp_path,
    active_status: RunStatus,
):
    database_url = _database_url(tmp_path)
    first_app = create_app(database_url=database_url)
    created = first_app.state.run_service.create_run(
        WORKFLOW_ID,
        WorkflowInputs(config={"samples": "private/input.tsv"}),
    )
    for status in (RunStatus.VALIDATING, RunStatus.PLANNED, RunStatus.QUEUED):
        first_app.state.run_service.transition_run(created.run_id, status)
    if active_status is RunStatus.RUNNING:
        first_app.state.run_service.transition_run(created.run_id, RunStatus.RUNNING)
    first_app.state.persistence.close()

    restarted_app = create_app(database_url=database_url)
    with ApiTestClient(restarted_app) as client:
        run_response = client.get(f"/api/v1/runs/{created.run_id}")
        events_response = client.get(f"/api/v1/runs/{created.run_id}/events")

    assert restarted_app.state.recovered_run_ids == (created.run_id,)
    assert run_response.status_code == 200
    run = run_response.json()["run"]
    assert run["status"] == "failed"
    assert run["error"]["code"] == "RUN_ORPHANED_AFTER_API_RESTART"
    assert "private/input.tsv" not in str(run["error"])
    recovery_event = events_response.json()["events"][-1]
    assert recovery_event["context"] == {
        "previous_status": active_status.value,
        "new_status": "failed",
        "reason_code": "WORKER_OWNERSHIP_NOT_CONFIRMED",
    }
    restarted_app.state.persistence.close()


@pytest.mark.parametrize("active_status", [RunStatus.QUEUED, RunStatus.RUNNING])
def test_api_restarts_preserve_worker_owned_active_run_without_recovery_noise(
    tmp_path,
    active_status: RunStatus,
):
    database_url = _database_url(tmp_path)
    first_app = create_app(database_url=database_url)
    created = first_app.state.run_service.create_run(
        WORKFLOW_ID,
        WorkflowInputs(config={"samples": "samples.tsv"}),
    )
    first_app.state.run_service.transition_run(created.run_id, RunStatus.VALIDATING)
    first_app.state.run_service.transition_run(created.run_id, RunStatus.PLANNED)
    assignment = first_app.state.run_service.ensure_execution_assignment(
        created.run_id,
        queue_name="runs",
    )
    assignment = first_app.state.run_service.mark_execution_dispatched(
        created.run_id,
        job_id=assignment.job_id,
    )
    first_app.state.run_service.transition_run(created.run_id, RunStatus.QUEUED)
    if active_status is RunStatus.RUNNING:
        claim = first_app.state.run_service.claim_execution_assignment(
            created.run_id,
            job_id=assignment.job_id,
            backend=assignment.backend,
            queue_name=assignment.queue_name,
        )
        assignment = claim.assignment
        first_app.state.run_service.transition_run(created.run_id, RunStatus.RUNNING)
    event_count = len(first_app.state.run_service.list_events(created.run_id))
    first_app.state.persistence.close()

    second_app = create_app(database_url=database_url)
    with ApiTestClient(second_app) as client:
        second_run_response = client.get(f"/api/v1/runs/{created.run_id}")
        second_events_response = client.get(f"/api/v1/runs/{created.run_id}/events")

    assert second_app.state.recovered_run_ids == ()
    assert second_run_response.json()["run"]["status"] == active_status.value
    assert len(second_events_response.json()["events"]) == event_count
    assert second_app.state.run_service.get_execution_assignment(created.run_id) == (
        assignment
    )
    second_app.state.persistence.close()

    third_app = create_app(database_url=database_url)
    with ApiTestClient(third_app) as client:
        third_run_response = client.get(f"/api/v1/runs/{created.run_id}")
        third_events_response = client.get(f"/api/v1/runs/{created.run_id}/events")

    assert third_app.state.recovered_run_ids == ()
    assert third_run_response.json()["run"]["status"] == active_status.value
    assert len(third_events_response.json()["events"]) == event_count
    third_app.state.persistence.close()


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
