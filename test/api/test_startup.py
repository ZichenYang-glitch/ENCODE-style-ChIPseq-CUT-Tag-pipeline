"""Lightweight startup smoke test for the API app factory."""

from __future__ import annotations

import inspect

import pytest

fastapi = pytest.importorskip("fastapi")

from encode_pipeline.api import dependencies  # noqa: E402
from encode_pipeline.api.main import create_app  # noqa: E402


def test_create_app_builds_expected_app() -> None:
    app = create_app()
    assert app.title == "Workflow Platform API"


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
    assert "/api/v1/runs/{run_id}/cancel" in paths
    assert "/api/v1/runs/{run_id}/preflight" in paths


def test_api_dependencies_are_async_to_avoid_testclient_threadpool_hang() -> None:
    dependency_names = [
        "get_registry",
        "get_validation_service",
        "get_agent_service",
        "get_run_service",
        "get_preflight_service",
    ]

    for name in dependency_names:
        assert inspect.iscoroutinefunction(getattr(dependencies, name))


def test_create_app_exposes_preflight_service_and_local_run_driver() -> None:
    app = create_app()
    assert hasattr(app.state, "preflight_service")
    assert hasattr(app.state, "local_run_driver")
    assert not hasattr(app.state, "stub_execution_driver")


def test_api_route_handlers_are_async_to_avoid_testclient_threadpool_hang() -> None:
    app = create_app()
    route_handlers = [
        endpoint
        for path, endpoint in _collect_route_endpoints(app.routes)
        if path.startswith("/api/v1/")
    ]

    assert route_handlers
    for endpoint in route_handlers:
        assert inspect.iscoroutinefunction(endpoint), endpoint
