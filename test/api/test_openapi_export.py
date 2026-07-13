"""Contract and drift tests for the OpenAPI export."""

from __future__ import annotations

import json
import os
import site
import subprocess
import sys
from pathlib import Path

import pytest
from pydantic import ValidationError

from encode_pipeline.api.main import create_app
from encode_pipeline.api.models import WorkflowInputLimitsResponse
from encode_pipeline.persistence import open_run_persistence
from encode_pipeline.platform.adapters import (
    MAX_AUTHORING_REQUEST_BYTES,
    MAX_SAMPLE_CELL_LENGTH,
    MAX_SAMPLE_COLUMN_NAME_LENGTH,
    MAX_SAMPLE_COLUMNS,
    MAX_SAMPLE_ROWS,
    WorkflowInputs,
)
from encode_pipeline.platform.runs import RunStatus
from encode_pipeline.services.defaults import (
    create_default_run_service,
    create_default_workflow_registry,
)
from encode_pipeline.workers.settings import (
    JOB_TIMEOUT_SECONDS_ENV,
    QUEUE_NAME_ENV,
    REDIS_API_READ_TIMEOUT_SECONDS_ENV,
    REDIS_CONNECT_TIMEOUT_SECONDS_ENV,
    REDIS_URL_ENV,
    WORKSPACE_ROOT_ENV,
)


EXPECTED_OPERATIONS = {
    ("GET", "/api/v1/workflows/"): "listWorkflows",
    ("GET", "/api/v1/workflows/{workflow_id}/schema"): "getWorkflowSchema",
    ("POST", "/api/v1/workflows/{workflow_id}/validate"): "validateWorkflow",
    ("POST", "/api/v1/workflows/{workflow_id}/runs"): "createRun",
    ("GET", "/api/v1/runs/{run_id}"): "getRun",
    ("POST", "/api/v1/runs/{run_id}/start"): "startRun",
    ("POST", "/api/v1/runs/{run_id}/cancel"): "cancelRun",
    ("GET", "/api/v1/runs/{run_id}/events"): "listRunEvents",
    ("GET", "/api/v1/runs/{run_id}/logs"): "listRunLogs",
    ("GET", "/api/v1/runs/{run_id}/artifacts"): "listRunArtifacts",
    ("GET", "/api/v1/runs/{run_id}/qc-metrics"): "listRunQcMetrics",
    (
        "GET",
        "/api/v1/runs/{run_id}/artifacts/{artifact_id}",
    ): "getRunArtifact",
    (
        "GET",
        "/api/v1/runs/{run_id}/artifacts/{artifact_id}/download",
    ): "downloadRunArtifact",
    ("POST", "/api/v1/runs/{run_id}/preflight"): "triggerPreflight",
    ("POST", "/api/v1/workflows/{workflow_id}/agent/chat"): "chatWithWorkflowAgent",
}

HTTP_METHODS = {"get", "put", "post", "delete", "options", "head", "patch", "trace"}

REPO_ROOT = Path(__file__).resolve().parents[2]
COMMITTED_OPENAPI_PATH = REPO_ROOT / "frontend" / "openapi.json"
WORKFLOW_ID = "encode-style-chipseq-cuttag-atac-mnase"
INPUT_LIMIT_FIELDS = (
    ("max_request_bytes", MAX_AUTHORING_REQUEST_BYTES),
    ("max_sample_rows", MAX_SAMPLE_ROWS),
    ("max_sample_columns", MAX_SAMPLE_COLUMNS),
    ("max_sample_column_name_length", MAX_SAMPLE_COLUMN_NAME_LENGTH),
    ("max_sample_cell_length", MAX_SAMPLE_CELL_LENGTH),
)


def _isolated_app_schema(tmp_path: Path, name: str) -> dict:
    runtime_root = tmp_path / name
    app = create_app(
        database_url=f"sqlite:///{runtime_root / 'platform.db'}",
        workspace_root=runtime_root / "workspaces",
    )
    try:
        return app.openapi()
    finally:
        try:
            app.state.run_queue.close()
        finally:
            app.state.persistence.close()


def _export_schema(output: Path, *, environment: dict[str, str] | None = None) -> None:
    result = subprocess.run(
        [sys.executable, "scripts/export_openapi.py", "--output", str(output)],
        check=False,
        capture_output=True,
        text=True,
        env=environment,
    )
    assert result.returncode == 0, (
        "OpenAPI export subprocess failed.\n"
        f"stdout:\n{result.stdout}\n"
        f"stderr:\n{result.stderr}"
    )


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
    _export_schema(output)

    app_schema = _isolated_app_schema(tmp_path, "app-schema")
    exported = json.loads(output.read_text(encoding="utf-8"))

    assert exported["info"]["title"] == app_schema["info"]["title"]
    assert exported["info"]["version"] == app_schema["info"]["version"]
    assert exported["paths"].keys() == app_schema["paths"].keys()
    assert "/api/v1/workflows/" in exported["paths"]
    assert "/api/v1/workflows/{workflow_id}/runs" in exported["paths"]
    assert "/api/v1/runs/{run_id}" in exported["paths"]
    assert "/api/v1/runs/{run_id}/preflight" in exported["paths"]
    assert "/api/v1/workflows/{workflow_id}/agent/chat" in exported["paths"]


def test_openapi_operations_match_expected_contract(tmp_path):
    schema = _isolated_app_schema(tmp_path, "operations")
    operations = _collect_operations(schema)

    assert len(operations) == len(EXPECTED_OPERATIONS)

    operation_ids = list(operations.values())
    assert len(operation_ids) == len(set(operation_ids)), "duplicate operation_id found"

    assert operations == EXPECTED_OPERATIONS


def test_workflow_authoring_contract_is_typed_versioned_and_bounded(tmp_path):
    schema = _isolated_app_schema(tmp_path, "workflow-authoring-contract")
    schema_response = schema["components"]["schemas"]["SchemaResponse"]
    assert "schema_hints" not in schema_response["properties"]
    assert "schema" in schema_response["required"]
    assert schema_response["properties"]["schema"] == {
        "anyOf": [
            {"$ref": "#/components/schemas/WorkflowSchemaResponse"},
            {"type": "null"},
        ]
    }

    for component_name in ("JsonValue-Input", "JsonValue-Output"):
        component_ref = f"#/components/schemas/{component_name}"
        json_value = schema["components"]["schemas"][component_name]
        assert json_value == {
            "anyOf": [
                {
                    "additionalProperties": {"$ref": component_ref},
                    "type": "object",
                },
                {
                    "items": {"$ref": component_ref},
                    "type": "array",
                },
                {"type": "string"},
                {"type": "integer"},
                {"type": "number"},
                {"type": "boolean"},
                {"type": "null"},
            ]
        }

    contract = schema["components"]["schemas"]["WorkflowSchemaResponse"]
    assert contract["additionalProperties"] is False
    assert set(contract["required"]) == {
        "schema_version",
        "schema_dialect",
        "coverage",
        "authoring_modes",
        "input_modes",
        "limits",
        "config_schema",
        "sample_schema",
        "option_schema",
    }
    assert contract["properties"]["schema_dialect"]["const"] == (
        "https://json-schema.org/draft/2020-12/schema"
    )

    limits = schema["components"]["schemas"]["WorkflowInputLimitsResponse"]
    assert set(limits["required"]) == {name for name, _ in INPUT_LIMIT_FIELDS}
    for field_name, ceiling in INPUT_LIMIT_FIELDS:
        assert limits["properties"][field_name]["minimum"] == ceiling
        assert limits["properties"][field_name]["maximum"] == ceiling

    for request_name in ("ValidationRequest", "RunCreateRequest"):
        request = schema["components"]["schemas"][request_name]
        assert request["additionalProperties"] is False
        samples = request["properties"]["samples"]["anyOf"]
        inline = next(item for item in samples if item.get("type") == "array")
        assert inline["minItems"] == 1
        assert inline["maxItems"] == 1000
        assert inline["items"]["minProperties"] == 1
        assert inline["items"]["maxProperties"] == 64
        assert inline["items"]["propertyNames"]["maxLength"] == 128
        assert inline["items"]["additionalProperties"]["maxLength"] == 4096


@pytest.mark.parametrize(("field_name", "ceiling"), INPUT_LIMIT_FIELDS)
@pytest.mark.parametrize("offset", (-1, 1), ids=("lower", "upper"))
def test_workflow_input_limit_response_requires_exact_platform_ceiling(
    field_name: str,
    ceiling: int,
    offset: int,
):
    values = dict(INPUT_LIMIT_FIELDS)
    values[field_name] = ceiling + offset

    with pytest.raises(ValidationError):
        WorkflowInputLimitsResponse(**values)


def test_authoring_operations_declare_stable_too_large_envelopes(tmp_path):
    schema = _isolated_app_schema(tmp_path, "authoring-too-large")
    validate = schema["paths"]["/api/v1/workflows/{workflow_id}/validate"]["post"]
    create = schema["paths"]["/api/v1/workflows/{workflow_id}/runs"]["post"]

    assert validate["responses"]["413"]["description"] == "Request body too large."
    assert create["responses"]["413"]["description"] == "Request body too large."
    assert validate["responses"]["413"]["content"]["application/json"]["schema"] == {
        "$ref": "#/components/schemas/ValidationResponse"
    }
    assert create["responses"]["413"]["content"]["application/json"]["schema"] == {
        "$ref": "#/components/schemas/RunResponse"
    }


def test_artifact_operations_do_not_publish_unreachable_422_responses(tmp_path):
    schema = _isolated_app_schema(tmp_path, "artifact-responses")

    for path in (
        "/api/v1/runs/{run_id}/artifacts",
        "/api/v1/runs/{run_id}/artifacts/{artifact_id}",
        "/api/v1/runs/{run_id}/artifacts/{artifact_id}/download",
    ):
        responses = schema["paths"][path]["get"]["responses"]
        assert "422" not in responses
        assert "4XX" in responses


def test_artifact_download_contract_is_binary_and_has_bounded_error_envelopes(
    tmp_path,
):
    schema = _isolated_app_schema(tmp_path, "artifact-download-contract")
    operation = schema["paths"][
        "/api/v1/runs/{run_id}/artifacts/{artifact_id}/download"
    ]["get"]

    assert operation["operationId"] == "downloadRunArtifact"
    assert set(operation["responses"]) == {"200", "404", "409", "4XX", "500"}
    binary = operation["responses"]["200"]["content"]["application/octet-stream"][
        "schema"
    ]
    assert binary == {"type": "string", "format": "binary"}
    for status in ("404", "409", "4XX", "500"):
        response_schema = operation["responses"][status]["content"]["application/json"][
            "schema"
        ]
        assert response_schema == {
            "$ref": "#/components/schemas/RunArtifactDownloadErrorResponse"
        }
    issue_properties = schema["components"]["schemas"]["ArtifactDownloadIssueResponse"][
        "properties"
    ]
    assert "technical_message" not in issue_properties
    context_ref = issue_properties["context"]["$ref"]
    assert context_ref == ("#/components/schemas/ArtifactDownloadIssueContextResponse")
    context_schema = schema["components"]["schemas"][
        "ArtifactDownloadIssueContextResponse"
    ]
    assert context_schema["additionalProperties"] is False
    assert set(context_schema["properties"]) == {"reason_code"}


def test_qc_metric_operation_is_lossless_and_has_declared_error_envelopes(tmp_path):
    schema = _isolated_app_schema(tmp_path, "qc-metric-contract")
    operation = schema["paths"]["/api/v1/runs/{run_id}/qc-metrics"]["get"]

    assert operation["operationId"] == "listRunQcMetrics"
    assert "422" not in operation["responses"]
    assert {"200", "400", "404", "4XX", "500"} <= set(operation["responses"])
    value_schema = schema["components"]["schemas"]["QcMetricResponse"]["properties"][
        "value"
    ]
    assert value_schema["type"] == "string"
    required = set(schema["components"]["schemas"]["QcMetricResponse"]["required"])
    assert {"sample_id", "experiment_id", "assay", "qc_flag"} <= required


def test_committed_openapi_matches_current_app(tmp_path):
    """Drift gate: the committed frontend/openapi.json must equal create_app().openapi().

    This prevents frontend code from running against a stale contract.
    """
    assert COMMITTED_OPENAPI_PATH.exists(), "committed openapi.json not found"
    committed = json.loads(COMMITTED_OPENAPI_PATH.read_text(encoding="utf-8"))
    current = _isolated_app_schema(tmp_path, "committed-contract")

    assert committed == current, (
        "committed openapi.json drifts from create_app().openapi()"
    )


def test_export_does_not_recover_environment_configured_runtime_database(
    tmp_path,
):
    runtime_database = tmp_path / "configured-runtime" / "platform.db"
    database_url = f"sqlite:///{runtime_database}"
    output = tmp_path / "configured-openapi.json"
    environment = dict(os.environ)
    environment["ENCODE_PIPELINE_DATABASE_URL"] = database_url

    persistence = open_run_persistence(database_url)
    try:
        run_service = create_default_run_service(
            registry=create_default_workflow_registry(),
            repository=persistence.repository,
        )
        created = run_service.create_run(WORKFLOW_ID, WorkflowInputs(config={}))
        run_service.transition_run(created.run_id, RunStatus.VALIDATING)
        events_before_export = run_service.list_events(created.run_id)
    finally:
        persistence.close()

    _export_schema(output, environment=environment)

    assert output.is_file()
    persistence = open_run_persistence(database_url)
    try:
        run_service = create_default_run_service(
            registry=create_default_workflow_registry(),
            repository=persistence.repository,
        )
        assert run_service.get_run(created.run_id).status is RunStatus.VALIDATING
        assert run_service.list_events(created.run_id) == events_before_export
    finally:
        persistence.close()


def test_export_does_not_create_default_runtime_database(tmp_path):
    fake_home = tmp_path / "isolated-home"
    fake_home.mkdir()
    output = tmp_path / "default-openapi.json"
    environment = dict(os.environ)
    environment.pop("ENCODE_PIPELINE_DATABASE_URL", None)
    # The CI environment may install editable dependencies in the user site.
    # Preserve that import location while HOME independently exercises runtime
    # default-path isolation in the child process.
    environment["PYTHONUSERBASE"] = site.getuserbase()
    environment["HOME"] = str(fake_home)

    _export_schema(output, environment=environment)

    assert output.is_file()
    assert not (fake_home / ".encode-pipeline" / "platform.db").exists()


def test_export_ignores_invalid_environment_worker_settings(tmp_path):
    output = tmp_path / "isolated-worker-openapi.json"
    environment = dict(os.environ)
    environment[REDIS_URL_ENV] = "not-a-redis-url"
    environment[QUEUE_NAME_ENV] = " "
    environment[JOB_TIMEOUT_SECONDS_ENV] = "not-an-integer"
    environment[REDIS_CONNECT_TIMEOUT_SECONDS_ENV] = "not-a-timeout"
    environment[REDIS_API_READ_TIMEOUT_SECONDS_ENV] = "not-a-timeout"
    environment[WORKSPACE_ROOT_ENV] = "relative/runtime-workspaces"

    _export_schema(output, environment=environment)

    assert output.is_file()
    exported = json.loads(output.read_text(encoding="utf-8"))
    assert exported["info"]["title"] == "Workflow Platform API"
