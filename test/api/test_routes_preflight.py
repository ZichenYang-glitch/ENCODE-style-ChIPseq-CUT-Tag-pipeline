"""Tests for the preflight trigger route."""

from __future__ import annotations

import asyncio
from pathlib import Path

import httpx

from encode_pipeline.api.main import create_app
from encode_pipeline.platform.results import Result
from encode_pipeline.platform.runs import RunStatus
from encode_pipeline.services.command_builder import CommandBuilder
from encode_pipeline.services.defaults import create_default_workspace_planner
from encode_pipeline.services.local_run_driver import LocalRunDriver
from encode_pipeline.services.materialization import WorkspaceMaterializer
from encode_pipeline.services.planning import ExecutionPlanner
from encode_pipeline.services.preflight import LocalPreflightService
from encode_pipeline.services.process_runner import ProcessResult, ProcessRunner
from encode_pipeline.services.run_repositories import ConcurrentRunUpdateError


class _FakeProcessRunner(ProcessRunner):
    def __init__(self, *, exit_code: int = 0):
        super().__init__(allowed_executables=("snakemake",))
        self._exit_code = exit_code

    def run(self, spec):
        return Result.success(
            ProcessResult(
                exit_code=self._exit_code,
                stdout="",
                stderr="",
            )
        )


def _make_valid_samples_tsv(tmp_path: Path) -> Path:
    samples_tsv = tmp_path / "samples.tsv"
    samples_tsv.write_text(
        "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\tpeak_mode\tgenome\tbowtie2_index\n"
        "S1\t/abs/S1_1.fq.gz\t/abs/S1_2.fq.gz\tPE\tchipseq\tH3K27ac\tnarrow\ths\t/abs/bt2/GRCh38\n",
        encoding="utf-8",
    )
    return samples_tsv


def _test_app(tmp_path: Path):
    app = create_app()
    registry = app.state.registry
    run_service = app.state.run_service
    local_run_driver = LocalRunDriver(
        run_service=run_service,
        materializer=WorkspaceMaterializer(),
        command_builder=CommandBuilder(registry=registry),
        workspace_root=tmp_path / "workspaces",
        process_runner=_FakeProcessRunner(),
    )
    app.state.local_run_driver = local_run_driver
    app.state.preflight_service = LocalPreflightService(
        run_service=run_service,
        execution_planner=ExecutionPlanner(run_service=run_service),
        workspace_planner=create_default_workspace_planner(registry=registry),
        local_run_driver=local_run_driver,
    )
    return app


def _client(tmp_path: Path) -> httpx.AsyncClient:
    return httpx.AsyncClient(
        transport=httpx.ASGITransport(app=_test_app(tmp_path)),
        base_url="http://testserver",
        follow_redirects=True,
    )


async def _create_run(client: httpx.AsyncClient, tmp_path: Path):
    samples_tsv = _make_valid_samples_tsv(tmp_path)
    response = await client.post(
        "/api/v1/workflows/encode-style-chipseq-cuttag-atac-mnase/runs",
        json={
            "config": {
                "samples": str(samples_tsv),
                "threads": 1,
                "genome_resources": {"hs": {"effective_genome_size": "hs"}},
            },
            "samples": str(samples_tsv),
            "options": {},
        },
    )
    assert response.status_code == 201
    return response.json()["run"]


def test_trigger_preflight_returns_202_and_validating_status(tmp_path):
    async def scenario() -> None:
        async with _client(tmp_path) as client:
            run = await _create_run(client, tmp_path)

            response = await client.post(f"/api/v1/runs/{run['run_id']}/preflight")

            assert response.status_code == 202
            body = response.json()
            assert body["ok"] is True
            # Response body is serialized before BackgroundTasks execute.
            assert body["run"]["status"] == "validating"

    asyncio.run(scenario())


def test_trigger_preflight_reaches_planned_on_get(tmp_path):
    async def scenario() -> None:
        async with _client(tmp_path) as client:
            run = await _create_run(client, tmp_path)

            await client.post(f"/api/v1/runs/{run['run_id']}/preflight")

            # httpx ASGITransport waits for the background task before returning control.
            response = await client.get(f"/api/v1/runs/{run['run_id']}")
            assert response.status_code == 200
            body = response.json()
            assert body["run"]["status"] == "planned"

            events = (await client.get(f"/api/v1/runs/{run['run_id']}/events")).json()[
                "events"
            ]
            event_types = [event["event_type"] for event in events]
            assert "workspace_materialized" in event_types
            assert "command_built" in event_types
            assert "dry_run_completed" in event_types
            assert "preflight_completed" in event_types

    asyncio.run(scenario())


def test_trigger_preflight_returns_404_for_unknown_run(tmp_path):
    async def scenario() -> None:
        async with _client(tmp_path) as client:
            response = await client.post("/api/v1/runs/does-not-exist/preflight")
            assert response.status_code == 404

    asyncio.run(scenario())


def test_trigger_preflight_returns_409_when_already_triggered(tmp_path):
    async def scenario() -> None:
        async with _client(tmp_path) as client:
            run = await _create_run(client, tmp_path)
            await client.post(f"/api/v1/runs/{run['run_id']}/preflight")

            response = await client.post(f"/api/v1/runs/{run['run_id']}/preflight")

            assert response.status_code == 409
            assert response.json()["issues"][0]["code"] == "PREFLIGHT_ALREADY_TRIGGERED"

    asyncio.run(scenario())


def test_trigger_preflight_returns_409_when_another_request_wins_transition(
    tmp_path,
    monkeypatch,
):
    async def scenario() -> None:
        app = _test_app(tmp_path)
        run_service = app.state.run_service
        original_transition = run_service.transition_run

        def concurrent_transition(run_id, to_status, **kwargs):
            original_transition(run_id, to_status, **kwargs)
            raise ConcurrentRunUpdateError("another request won")

        monkeypatch.setattr(run_service, "transition_run", concurrent_transition)
        async with httpx.AsyncClient(
            transport=httpx.ASGITransport(app=app),
            base_url="http://testserver",
        ) as client:
            run = await _create_run(client, tmp_path)

            response = await client.post(f"/api/v1/runs/{run['run_id']}/preflight")

            assert response.status_code == 409
            assert response.json()["issues"][0]["code"] == (
                "PREFLIGHT_ALREADY_TRIGGERED"
            )
            assert response.json()["issues"][0]["context"] == {
                "current_status": "validating"
            }
            assert run_service.get_run(run["run_id"]).status is RunStatus.VALIDATING

    asyncio.run(scenario())


def test_trigger_preflight_remains_responsive_after_prior_requests(
    tmp_path,
):
    """Preflight completes after earlier API calls in one ASGI client session."""

    async def scenario() -> None:
        async with _client(tmp_path) as client:
            run = await _create_run(client, tmp_path)

            for _ in range(3):
                response = await client.get("/api/v1/workflows")
                assert response.status_code == 200
                response = await client.get(
                    "/api/v1/workflows/encode-style-chipseq-cuttag-atac-mnase/schema"
                )
                assert response.status_code == 200

            response = await client.post(f"/api/v1/runs/{run['run_id']}/preflight")
            assert response.status_code == 202
            assert response.json()["run"]["status"] == "validating"

            response = await client.get(f"/api/v1/runs/{run['run_id']}")
            assert response.status_code == 200
            assert response.json()["run"]["status"] == "planned"

    asyncio.run(scenario())
