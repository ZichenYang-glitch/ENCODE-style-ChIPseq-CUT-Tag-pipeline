"""Tests for the FastAPI validation MVP routes."""

from __future__ import annotations

import os
import tempfile
from collections.abc import Iterator

import pytest

from encode_pipeline.api.main import create_app
from encode_pipeline.platform.adapters import WorkflowCapabilities, WorkflowMetadata
from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.services.validation import ValidationService
from api_test_client import ApiTestClient

fastapi = pytest.importorskip("fastapi")


@pytest.fixture
def client() -> Iterator[ApiTestClient]:
    """Default app wired to the bundled ENCODE-style adapter."""
    app = create_app()
    with ApiTestClient(app) as tc:
        yield tc


@pytest.fixture
def minimal_config_and_samples() -> Iterator[tuple[str, str]]:
    """Provide a minimal valid config file and sample sheet."""
    with tempfile.TemporaryDirectory() as tmpdir:
        samples_path = os.path.join(tmpdir, "samples.tsv")
        config_path = os.path.join(tmpdir, "config.yaml")
        with open(samples_path, "w") as f:
            f.write(
                "sample\tfastq_1\tlayout\tassay\ttarget\tpeak_mode\tgenome\tbowtie2_index\n"
                "ctl1\tctl1.fastq.gz\tSE\tchipseq\tH3K4me3\tnarrow\thg38\t/path/to/bowtie2\n"
            )
        with open(config_path, "w") as f:
            f.write(f"samples: {samples_path}\n")
        yield config_path, samples_path


@pytest.fixture
def unsupported_capability_client() -> Iterator[ApiTestClient]:
    """App with a workflow that does not support validation."""

    class StubAdapter:
        metadata = WorkflowMetadata(
            workflow_id="no-validation",
            name="No Validation",
            version="0.0.1",
            description="Adapter without validation capability.",
            engines=("snakemake",),
            tags=("stub",),
        )
        capabilities = WorkflowCapabilities(supports=())

        def schema(self):
            from encode_pipeline.platform.adapters import WorkflowSchema

            return WorkflowSchema()

        def validate(self, inputs):
            from encode_pipeline.platform.results import Result

            return Result.success(None)

        def preview_dag(self, inputs):
            from encode_pipeline.platform.results import Result

            return Result.success(None)

        def plan_workspace(self, inputs, workspace):
            from encode_pipeline.platform.results import Result

            return Result.success(None)

        def build_command(self, plan):
            from encode_pipeline.platform.results import Result

            return Result.success(None)

        def extract_artifacts(self, inputs, workspace):
            from encode_pipeline.platform.results import Result

            return Result.success(())

    registry = WorkflowRegistry(adapters=[StubAdapter()])
    service = ValidationService(registry=registry)

    app = create_app()
    app.state.registry = registry
    app.state.validation_service = service
    with ApiTestClient(app) as tc:
        yield tc


def test_list_workflows_returns_encode_workflow(client: ApiTestClient) -> None:
    response = client.get("/api/v1/workflows")
    assert response.status_code == 200
    data = response.json()
    assert data["ok"] is True
    assert len(data["workflows"]) == 1
    item = data["workflows"][0]
    assert item["metadata"]["workflow_id"] == "encode-style-chipseq-cuttag-atac-mnase"
    assert "validation" in item["capabilities"]["supports"]
    assert data["issues"] == []


def test_get_schema_returns_hints(client: ApiTestClient) -> None:
    workflow_id = "encode-style-chipseq-cuttag-atac-mnase"
    response = client.get(f"/api/v1/workflows/{workflow_id}/schema")
    assert response.status_code == 200
    data = response.json()
    assert data["ok"] is True
    assert data["workflow_id"] == workflow_id
    assert "config_schema" in data["schema_hints"]
    assert "sample_schema" in data["schema_hints"]
    assert "option_schema" in data["schema_hints"]
    assert data["issues"] == []


def test_get_schema_unknown_workflow_returns_404(client: ApiTestClient) -> None:
    response = client.get("/api/v1/workflows/missing-workflow/schema")
    assert response.status_code == 404
    data = response.json()
    assert data["ok"] is False
    assert data["workflow_id"] == "missing-workflow"
    assert data["schema_hints"] == {}
    assert data["issues"][0]["code"] == "WORKFLOW_NOT_FOUND"


def test_validate_success(
    client: ApiTestClient,
    minimal_config_and_samples: tuple[str, str],
) -> None:
    config_path, samples_path = minimal_config_and_samples
    workflow_id = "encode-style-chipseq-cuttag-atac-mnase"
    response = client.post(
        f"/api/v1/workflows/{workflow_id}/validate",
        json={
            "config": {"samples": samples_path},
            "samples": samples_path,
            "options": {"strict_inputs": False},
        },
    )
    assert response.status_code == 200
    data = response.json()
    assert data["ok"] is True
    assert data["workflow_id"] == workflow_id
    assert data["issues"] == []


def test_validate_failure_with_structured_issues(client: ApiTestClient) -> None:
    workflow_id = "encode-style-chipseq-cuttag-atac-mnase"
    response = client.post(
        f"/api/v1/workflows/{workflow_id}/validate",
        json={"config": {}, "samples": "nonexistent.tsv"},
    )
    assert response.status_code == 200
    data = response.json()
    assert data["ok"] is False
    assert data["workflow_id"] == workflow_id
    assert len(data["issues"]) >= 1
    codes = {issue["code"] for issue in data["issues"]}
    assert "ENCODE_CONFIG_INVALID" in codes or "ENCODE_SAMPLES_INVALID" in codes


def test_validate_unknown_workflow_returns_404(client: ApiTestClient) -> None:
    response = client.post(
        "/api/v1/workflows/missing-workflow/validate",
        json={"config": {}, "samples": "samples.tsv"},
    )
    assert response.status_code == 404
    data = response.json()
    assert data["ok"] is False
    assert data["issues"][0]["code"] == "WORKFLOW_NOT_FOUND"


def test_validate_unsupported_capability_returns_409(
    unsupported_capability_client: ApiTestClient,
) -> None:
    response = unsupported_capability_client.post(
        "/api/v1/workflows/no-validation/validate",
        json={"config": {}, "samples": "samples.tsv"},
    )
    assert response.status_code == 409
    data = response.json()
    assert data["ok"] is False
    assert data["issues"][0]["code"] == "WORKFLOW_CAPABILITY_UNSUPPORTED"


def test_validate_malformed_json_returns_400(client: ApiTestClient) -> None:
    workflow_id = "encode-style-chipseq-cuttag-atac-mnase"
    response = client.post(
        f"/api/v1/workflows/{workflow_id}/validate",
        content="not-json",
        headers={"content-type": "application/json"},
    )
    assert response.status_code == 400
    data = response.json()
    assert data["ok"] is False
    assert data["issues"][0]["code"] == "API_REQUEST_INVALID"


def test_validate_wrong_field_types_returns_400(client: ApiTestClient) -> None:
    workflow_id = "encode-style-chipseq-cuttag-atac-mnase"
    response = client.post(
        f"/api/v1/workflows/{workflow_id}/validate",
        json={"config": "must-be-object", "samples": "samples.tsv"},
    )
    assert response.status_code == 400
    data = response.json()
    assert data["ok"] is False
    assert data["issues"][0]["code"] == "API_REQUEST_INVALID"


def test_api_routes_import_boundary() -> None:
    """API routes must not import validator/sample-loader/workflow-engine internals."""
    from encode_pipeline.api import routes
    from encode_pipeline.api.routes import workflows

    assert "encode_pipeline.config.validator" not in str(routes)
    assert "encode_pipeline.samples.load" not in str(workflows)
    assert "snakemake" not in str(workflows)
    assert "Snakefile" not in str(workflows)
