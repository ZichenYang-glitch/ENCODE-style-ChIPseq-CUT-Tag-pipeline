"""Contract and drift tests for the OpenAPI export."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

from encode_pipeline.api.main import create_app


EXPECTED_OPERATIONS = {
    ("GET", "/api/v1/workflows/"): "listWorkflows",
    ("GET", "/api/v1/workflows/{workflow_id}/schema"): "getWorkflowSchema",
    ("POST", "/api/v1/workflows/{workflow_id}/validate"): "validateWorkflow",
    ("POST", "/api/v1/workflows/{workflow_id}/runs"): "createRun",
    ("GET", "/api/v1/runs/{run_id}"): "getRun",
    ("POST", "/api/v1/runs/{run_id}/cancel"): "cancelRun",
    ("GET", "/api/v1/runs/{run_id}/events"): "listRunEvents",
    ("GET", "/api/v1/runs/{run_id}/logs"): "listRunLogs",
    ("POST", "/api/v1/runs/{run_id}/preflight"): "triggerPreflight",
    ("POST", "/api/v1/workflows/{workflow_id}/agent/chat"): "chatWithWorkflowAgent",
}

HTTP_METHODS = {"get", "put", "post", "delete", "options", "head", "patch", "trace"}

REPO_ROOT = Path(__file__).resolve().parents[2]
COMMITTED_OPENAPI_PATH = REPO_ROOT / "frontend" / "openapi.json"


def _collect_operations(schema: dict) -> dict[tuple[str, str], str]:
    operations: dict[tuple[str, str], str] = {}
    for path, path_item in schema.get("paths", {}).items():
        for method, operation in path_item.items():
            if method.lower() not in HTTP_METHODS:
                continue
            if not isinstance(operation, dict):
                continue
            op_id = operation.get("operationId")
            if op_id is None:
                raise AssertionError(f"Missing operationId for {method.upper()} {path}")
            operations[(method.upper(), path)] = op_id
    return operations


def test_exported_openapi_matches_app_openapi(tmp_path):
    output = tmp_path / "openapi.json"
    result = subprocess.run(
        [sys.executable, "scripts/export_openapi.py", "--output", str(output)],
        check=True,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0

    app_schema = create_app().openapi()
    exported = json.loads(output.read_text(encoding="utf-8"))

    assert exported["info"]["title"] == app_schema["info"]["title"]
    assert exported["info"]["version"] == app_schema["info"]["version"]
    assert exported["paths"].keys() == app_schema["paths"].keys()
    assert "/api/v1/workflows/" in exported["paths"]
    assert "/api/v1/workflows/{workflow_id}/runs" in exported["paths"]
    assert "/api/v1/runs/{run_id}" in exported["paths"]
    assert "/api/v1/runs/{run_id}/preflight" in exported["paths"]
    assert "/api/v1/workflows/{workflow_id}/agent/chat" in exported["paths"]


def test_openapi_operations_match_expected_contract():
    schema = create_app().openapi()
    operations = _collect_operations(schema)

    assert len(operations) == len(EXPECTED_OPERATIONS)

    operation_ids = list(operations.values())
    assert len(operation_ids) == len(set(operation_ids)), "duplicate operation_id found"

    assert operations == EXPECTED_OPERATIONS


def test_committed_openapi_matches_current_app():
    """Drift gate: the committed frontend/openapi.json must equal create_app().openapi().

    This prevents frontend code from running against a stale contract.
    """
    assert COMMITTED_OPENAPI_PATH.exists(), "committed openapi.json not found"
    committed = json.loads(COMMITTED_OPENAPI_PATH.read_text(encoding="utf-8"))
    current = create_app().openapi()

    assert committed == current, (
        "committed openapi.json drifts from create_app().openapi()"
    )
