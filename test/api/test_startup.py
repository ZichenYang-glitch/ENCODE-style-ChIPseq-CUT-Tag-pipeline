"""Lightweight startup smoke test for the API app factory."""

from __future__ import annotations

import pytest

fastapi = pytest.importorskip("fastapi")

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
                _collect_route_paths(route.original_router.routes, prefix + include_prefix)
            )
        elif hasattr(route, "path"):
            paths.add(prefix + route.path)
    return paths


def test_expected_routes_are_registered() -> None:
    app = create_app()
    paths = {path.rstrip("/") for path in _collect_route_paths(app.routes)}
    assert "/api/v1/workflows" in paths
    assert "/api/v1/workflows/{workflow_id}/schema" in paths
    assert "/api/v1/workflows/{workflow_id}/validate" in paths
    assert "/api/v1/workflows/{workflow_id}/agent/chat" in paths
