"""Contract tests for the OpenAPI export script."""

from __future__ import annotations

import json
import subprocess
import sys

from encode_pipeline.api.main import create_app


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
