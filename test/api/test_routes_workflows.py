"""Tests for the FastAPI validation MVP routes."""

from __future__ import annotations

import os
import tempfile
from collections.abc import Iterator

import pytest

from encode_pipeline.api.main import create_app
from encode_pipeline.platform.adapters import (
    MAX_AUTHORING_REQUEST_BYTES,
    MAX_SAMPLE_CELL_LENGTH,
    MAX_SAMPLE_COLUMN_NAME_LENGTH,
    MAX_SAMPLE_COLUMNS,
    MAX_SAMPLE_ROWS,
    WorkflowCapabilities,
    WorkflowMetadata,
)
from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.services.validation import ValidationService
from encode_pipeline.services.validated_inputs import ValidatedInputService
from encode_pipeline.services.workflow_builds import WorkflowBuildIdentityProvider
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
        fastq_path = os.path.join(tmpdir, "ctl1.fastq.gz")
        bowtie2_path = os.path.join(tmpdir, "bowtie2", "index")
        with open(samples_path, "w") as f:
            f.write(
                "sample\tfastq_1\tlayout\tassay\ttarget\tpeak_mode\tgenome\tbowtie2_index\n"
                f"ctl1\t{fastq_path}\tSE\tchipseq\tH3K4me3\tnarrow\thg38\t{bowtie2_path}\n"
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
    app.state.validated_input_service = ValidatedInputService(
        registry=registry,
        validation_service=service,
        build_identity_provider=WorkflowBuildIdentityProvider(registry),
        repository=app.state.persistence.repository,
    )
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


def test_get_schema_returns_versioned_renderable_contract(
    client: ApiTestClient,
) -> None:
    workflow_id = "encode-style-chipseq-cuttag-atac-mnase"
    response = client.get(f"/api/v1/workflows/{workflow_id}/schema")
    assert response.status_code == 200
    data = response.json()
    assert data["ok"] is True
    assert data["workflow_id"] == workflow_id
    assert data["schema"]["schema_version"] == "1.0.0"
    assert data["schema"]["schema_dialect"] == (
        "https://json-schema.org/draft/2020-12/schema"
    )
    assert data["schema"]["coverage"] == {
        "config": "partial",
        "samples": "complete",
        "options": "complete",
    }
    assert data["schema"]["sample_schema"]["type"] == "array"
    assert data["schema"]["limits"] == {
        "max_request_bytes": MAX_AUTHORING_REQUEST_BYTES,
        "max_sample_rows": MAX_SAMPLE_ROWS,
        "max_sample_columns": MAX_SAMPLE_COLUMNS,
        "max_sample_column_name_length": MAX_SAMPLE_COLUMN_NAME_LENGTH,
        "max_sample_cell_length": MAX_SAMPLE_CELL_LENGTH,
    }
    assert data["issues"] == []


def test_get_schema_unknown_workflow_returns_404(client: ApiTestClient) -> None:
    response = client.get("/api/v1/workflows/missing-workflow/schema")
    assert response.status_code == 404
    data = response.json()
    assert data["ok"] is False
    assert data["workflow_id"] == "missing-workflow"
    assert data["schema"] is None
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
    assert data["value"] is None
    assert data["snapshot"]["workflow_id"] == workflow_id
    assert data["snapshot"]["schema_version"] == "1.0.0"
    assert len(data["snapshot"]["payload_digest"]) == 64
    assert data["issues"] == []


def test_validate_inline_rows_without_config_samples_is_successful(
    client: ApiTestClient,
    tmp_path,
) -> None:
    workflow_id = "encode-style-chipseq-cuttag-atac-mnase"
    row = {
        "sample": "S1",
        "fastq_1": str((tmp_path / "S1.R1.fastq.gz").resolve()),
        "fastq_2": str((tmp_path / "S1.R2.fastq.gz").resolve()),
        "layout": "PE",
        "assay": "chipseq",
        "target": "CTCF",
        "peak_mode": "narrow",
        "genome": "hs",
        "bowtie2_index": str((tmp_path / "indices/hs").resolve()),
    }

    response = client.post(
        f"/api/v1/workflows/{workflow_id}/validate",
        json={"config": {}, "samples": [row], "options": {}},
    )

    assert response.status_code == 200
    data = response.json()
    assert data["ok"] is True
    assert data["workflow_id"] == workflow_id
    assert data["value"] is None
    assert data["snapshot"]["snapshot_id"].startswith("vsnap_")
    assert data["issues"] == []
    assert tempfile.gettempdir() not in response.text


def test_validate_rejects_oversized_sample_cell_with_safe_400(
    client: ApiTestClient,
) -> None:
    workflow_id = "encode-style-chipseq-cuttag-atac-mnase"
    response = client.post(
        f"/api/v1/workflows/{workflow_id}/validate",
        json={"config": {}, "samples": [{"sample": "x" * 4097}]},
    )

    assert response.status_code == 400
    data = response.json()
    assert data["ok"] is False
    assert data["value"] is None
    assert data["issues"][0]["code"] == "API_REQUEST_INVALID"
    assert "x" * 128 not in response.text


def test_validate_rejects_control_characters_in_sample_path_with_safe_400(
    client: ApiTestClient,
) -> None:
    workflow_id = "encode-style-chipseq-cuttag-atac-mnase"
    response = client.post(
        f"/api/v1/workflows/{workflow_id}/validate",
        json={"config": {}, "samples": "private\npath.tsv"},
    )

    assert response.status_code == 400
    data = response.json()
    assert data["issues"][0]["code"] == "API_REQUEST_INVALID"
    assert "private" not in response.text


@pytest.mark.parametrize(
    "unsafe_value",
    [2**53, -(2**53)],
)
def test_validate_rejects_unsafe_json_integers_without_echo(
    client: ApiTestClient,
    unsafe_value: int,
) -> None:
    workflow_id = "encode-style-chipseq-cuttag-atac-mnase"
    response = client.post(
        f"/api/v1/workflows/{workflow_id}/validate",
        json={"config": {"threads": unsafe_value}, "samples": None, "options": {}},
    )

    assert response.status_code == 400
    body = response.json()
    assert body["issues"][0]["code"] == "API_REQUEST_INVALID"
    assert body["issues"][0]["technical_message"] is None
    assert str(unsafe_value) not in response.text


def test_separate_successful_validations_create_distinct_opaque_snapshots(
    client: ApiTestClient,
) -> None:
    workflow_id = "encode-style-chipseq-cuttag-atac-mnase"
    request = {
        "config": {},
        "samples": [
            {
                "sample": "S1",
                "fastq_1": "/tmp/S1.R1.fastq.gz",
                "layout": "SE",
                "assay": "chipseq",
                "target": "CTCF",
                "peak_mode": "narrow",
                "genome": "hs",
                "bowtie2_index": "/tmp/indices/hs",
            }
        ],
        "options": {},
    }

    first = client.post(
        f"/api/v1/workflows/{workflow_id}/validate", json=request
    ).json()["snapshot"]
    second = client.post(
        f"/api/v1/workflows/{workflow_id}/validate", json=request
    ).json()["snapshot"]

    assert first["snapshot_id"] != second["snapshot_id"]
    assert first["payload_digest"] == second["payload_digest"]


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
