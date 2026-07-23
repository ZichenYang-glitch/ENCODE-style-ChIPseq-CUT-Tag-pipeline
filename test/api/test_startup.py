"""Lightweight startup smoke test for the API app factory."""

from __future__ import annotations

import hashlib
import inspect
import json
from pathlib import Path

import pytest

fastapi = pytest.importorskip("fastapi")

from encode_pipeline.api import dependencies  # noqa: E402
from encode_pipeline.api.main import create_app  # noqa: E402
from encode_pipeline.adapters.bulk_rnaseq.deployment import (  # noqa: E402
    MANAGED_DOCKER_EXECUTABLE_ENV,
    MANAGED_DOCKER_SOCKET_ENV,
    RUNTIME_ROOT_ENV,
    TRANSCRIPTOME_BINDING_MANIFEST_ENV,
)
from encode_pipeline.services.workflow_info import WorkflowInfoService  # noqa: E402
from encode_pipeline.workers.rq_queue import RqRunQueue  # noqa: E402
from encode_pipeline.workers.settings import (  # noqa: E402
    QUEUE_NAME_ENV,
    REDIS_URL_ENV,
    WORKSPACE_ROOT_ENV,
)


def test_create_app_builds_expected_app() -> None:
    app = create_app()
    assert app.title == "HelixWeave API"
    assert app.description == "Reproducible omics workflows, from inputs to evidence."


def _collect_route_paths(routes, prefix: str = "") -> set[str]:
    paths: set[str] = set()
    for route in routes:
        if type(route).__name__ == "_IncludedRouter":
            include_prefix = getattr(route.include_context, "prefix", "")
            paths.update(
                _collect_route_paths(
                    route.original_router.routes, prefix + include_prefix
                )
            )
        elif hasattr(route, "path"):
            paths.add(prefix + route.path)
    return paths


def _collect_route_endpoints(routes, prefix: str = "") -> list[tuple[str, object]]:
    endpoints: list[tuple[str, object]] = []
    for route in routes:
        if type(route).__name__ == "_IncludedRouter":
            include_prefix = getattr(route.include_context, "prefix", "")
            endpoints.extend(
                _collect_route_endpoints(
                    route.original_router.routes, prefix + include_prefix
                )
            )
        elif hasattr(route, "path") and hasattr(route, "endpoint"):
            endpoints.append((prefix + route.path, route.endpoint))
    return endpoints


def test_expected_routes_are_registered() -> None:
    app = create_app()
    paths = {path.rstrip("/") for path in _collect_route_paths(app.routes)}
    assert "/api/v1/workflows" in paths
    assert "/api/v1/workflows/{workflow_id}/schema" in paths
    assert "/api/v1/workflows/{workflow_id}/validate" in paths
    assert "/api/v1/workflows/{workflow_id}/agent/chat" in paths
    assert "/api/v1/workflows/{workflow_id}/runs" in paths
    assert "/api/v1/runs/{run_id}" in paths
    assert "/api/v1/runs/{run_id}/events" in paths
    assert "/api/v1/runs/{run_id}/logs" in paths
    assert "/api/v1/runs/{run_id}/start" in paths
    assert "/api/v1/runs/{run_id}/cancel" in paths
    assert "/api/v1/runs/{run_id}/preflight" in paths
    assert "/api/v1/runs/{run_id}/artifacts/{artifact_id}/download" in paths


def test_api_dependencies_are_async_to_avoid_testclient_threadpool_hang() -> None:
    dependency_names = [
        "get_registry",
        "get_validation_service",
        "get_agent_service",
        "get_run_service",
        "get_artifact_download_service",
        "get_run_submission_service",
        "get_run_cancellation_service",
        "get_preflight_service",
    ]

    for name in dependency_names:
        assert inspect.iscoroutinefunction(getattr(dependencies, name))


def test_create_app_exposes_preflight_service_and_local_run_driver() -> None:
    app = create_app()
    assert hasattr(app.state, "run_submission_service")
    assert hasattr(app.state, "run_cancellation_service")
    assert hasattr(app.state, "preflight_service")
    assert hasattr(app.state, "artifact_download_service")
    assert hasattr(app.state, "local_run_driver")
    assert not hasattr(app.state, "stub_execution_driver")


def test_create_app_aligns_command_and_identity_project_roots(tmp_path: Path) -> None:
    project_root = (tmp_path / "source").resolve()
    app = create_app(
        database_url=f"sqlite:///{tmp_path / 'platform.db'}",
        project_root=project_root,
    )

    assert app.state.build_identity_provider.project_root == project_root
    assert app.state.local_run_driver._command_builder._project_root == project_root

    app.state.run_queue.close()
    app.state.persistence.close()


def test_create_app_composes_shared_api_and_worker_settings(
    tmp_path: Path,
    monkeypatch,
) -> None:
    database_url = f"sqlite:///{tmp_path / 'platform.db'}"
    workspace_root = tmp_path / "shared-workspaces"
    monkeypatch.setenv(REDIS_URL_ENV, "redis://redis.internal:6380/4")
    monkeypatch.setenv(QUEUE_NAME_ENV, "epigenomics")
    monkeypatch.setenv(WORKSPACE_ROOT_ENV, str(workspace_root))

    app = create_app(database_url=database_url)

    assert app.state.database_url == database_url
    assert app.state.workspace_root == workspace_root
    assert app.state.worker_settings.database_url == database_url
    assert app.state.worker_settings.workspace_root == workspace_root
    assert app.state.worker_settings.redis_url == "redis://redis.internal:6380/4"
    assert app.state.worker_settings.queue_name == "epigenomics"
    assert isinstance(app.state.run_queue, RqRunQueue)
    assert app.state.run_queue.queue_name == "epigenomics"
    assert app.state.local_run_driver._workspace_root == workspace_root

    app.state.run_queue.close()
    app.state.persistence.close()


def test_create_app_keeps_bulk_authoring_when_configured_docker_is_unavailable(
    tmp_path: Path,
    monkeypatch,
) -> None:
    transcript = (tmp_path / "fixture/transcripts.fa").resolve()
    transcript.parent.mkdir()
    transcript.write_text(">TX1\nACGT\n", encoding="utf-8")
    transcript_sha256 = hashlib.sha256(transcript.read_bytes()).hexdigest()
    binding_manifest = (tmp_path / "transcriptome-binding.json").resolve()
    binding_manifest.write_text(
        json.dumps(
            {
                "schema_version": "1.0.0",
                "reference_id": "tiny",
                "fasta_sha256": "a" * 64,
                "gtf_sha256": "b" * 64,
                "transcript_fasta": str(transcript),
                "transcript_fasta_sha256": transcript_sha256,
            }
        ),
        encoding="utf-8",
    )
    docker = (tmp_path / "bin/docker").resolve()
    docker.parent.mkdir()
    docker.write_text("#!/bin/sh\nexit 1\n", encoding="utf-8")
    docker.chmod(0o755)
    missing_socket = (tmp_path / "run/docker.sock").resolve()
    monkeypatch.setenv(RUNTIME_ROOT_ENV, str((tmp_path / "runtime").resolve()))
    monkeypatch.setenv(
        TRANSCRIPTOME_BINDING_MANIFEST_ENV,
        str(binding_manifest),
    )
    monkeypatch.setenv(MANAGED_DOCKER_EXECUTABLE_ENV, str(docker))
    monkeypatch.setenv(MANAGED_DOCKER_SOCKET_ENV, str(missing_socket))

    app = create_app(database_url=f"sqlite:///{tmp_path / 'platform.db'}")

    descriptor = WorkflowInfoService(app.state.registry).get_descriptor("bulk-rnaseq")
    assert descriptor.is_success
    assert descriptor.value.availability.execution == "unavailable"
    assert descriptor.value.capabilities.supports == (
        "validation",
        "input_authoring",
    )
    assert app.state.local_run_driver._process_runner._allowed_executables == (
        "snakemake",
    )
    assert app.state.local_run_driver._process_runner._managed_container_cleaner is None

    app.state.run_queue.close()
    app.state.persistence.close()


def test_only_explicit_blocking_routes_use_fastapi_threadpool() -> None:
    app = create_app()
    route_handlers = [
        (path.rstrip("/"), endpoint)
        for path, endpoint in _collect_route_endpoints(app.routes)
        if path.startswith("/api/v1/")
    ]

    assert route_handlers
    for path, endpoint in route_handlers:
        if path in {
            "/api/v1/workflows",
            "/api/v1/workflows/{workflow_id}",
            "/api/v1/workflows/{workflow_id}/validate",
            "/api/v1/workflows/{workflow_id}/runs",
            "/api/v1/runs",
            "/api/v1/runs/{run_id}/start",
            "/api/v1/runs/{run_id}/cancel",
            "/api/v1/runs/{run_id}/qc-metrics",
            "/api/v1/runs/{run_id}/artifacts/{artifact_id}/download",
        }:
            assert not inspect.iscoroutinefunction(endpoint), endpoint
        else:
            assert inspect.iscoroutinefunction(endpoint), endpoint
