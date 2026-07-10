"""Contract and drift tests for the OpenAPI export."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

from encode_pipeline.api.main import create_app


EXPECTED_OPERATION_IDS = frozenset(
    {
        "listWorkflows",
        "getWorkflowSchema",
        "validateWorkflow",
        "createRun",
        "getRun",
        "cancelRun",
        "listRunEvents",
        "listRunLogs",
        "triggerPreflight",
        "chatWithWorkflowAgent",
    }
)

REPO_ROOT = Path(__file__).resolve().parents[2]
COMMITTED_OPENAPI_PATH = REPO_ROOT / "frontend" / "openapi.json"


def _operation_ids(schema: dict) -> set[str]:
    ids: set[str] = set()
    for path_item in schema.get("paths", {}).values():
        for operation in path_item.values():
            if isinstance(operation, dict):
                op_id = operation.get("operationId")
                if op_id is not None:
                    ids.add(op_id)
    return ids


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


def test_openapi_operation_ids_are_unique_and_stable():
    schema = create_app().openapi()
    operation_ids = _operation_ids(schema)

    assert len(operation_ids) == len(EXPECTED_OPERATION_IDS)
    assert operation_ids == EXPECTED_OPERATION_IDS


def test_committed_openapi_matches_current_app():
    """Drift gate: the committed frontend/openapi.json must equal create_app().openapi().

    This prevents frontend code from running against a stale contract.
    """
    assert COMMITTED_OPENAPI_PATH.exists(), "committed openapi.json not found"
    committed = json.loads(COMMITTED_OPENAPI_PATH.read_text(encoding="utf-8"))
    current = create_app().openapi()

    assert committed == current, "committed openapi.json drifts from create_app().openapi()"
