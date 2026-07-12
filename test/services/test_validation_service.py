"""Tests for the workflow validation service."""

import os
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

from encode_pipeline.platform.adapters import (
    CommandSpec,
    DagPreview,
    WorkflowCapabilities,
    WorkflowInputs,
    WorkflowMetadata,
    WorkflowSchema,
    WorkspacePlan,
)
from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.platform.results import Issue, Result
from encode_pipeline.services.validation import ValidationService


SRC_ROOT = Path(__file__).resolve().parents[2] / "src"


class FakeAdapter:
    def __init__(
        self,
        workflow_id: str = "fake",
        *,
        supports: tuple[str, ...] = ("validation",),
        result: Result[object] | None = None,
        exception: Exception | None = None,
    ) -> None:
        self.metadata = WorkflowMetadata(
            workflow_id=workflow_id,
            name="Fake Workflow",
            version="1.0.0",
        )
        self.capabilities = WorkflowCapabilities(supports=supports)
        self.result = result or Result.success({"validated": True})
        self.exception = exception
        self.validate_called = False
        self.received_inputs: WorkflowInputs | None = None

    def schema(self) -> WorkflowSchema:
        return WorkflowSchema(config_schema={"type": "object"})

    def validate(self, inputs: WorkflowInputs) -> Result[object]:
        self.validate_called = True
        self.received_inputs = inputs
        if self.exception is not None:
            raise self.exception
        return self.result

    def preview_dag(self, inputs: WorkflowInputs) -> Result[DagPreview]:
        return Result.success(DagPreview())

    def plan_workspace(
        self,
        inputs: WorkflowInputs,
        workspace: str | Path,
    ) -> Result[WorkspacePlan]:
        return Result.success(WorkspacePlan(directories=[str(workspace)]))

    def build_command(self, plan: WorkspacePlan) -> Result[CommandSpec]:
        return Result.success(CommandSpec(argv=["run-workflow"]))

    def extract_artifacts(self, inputs, workspace):
        return Result.success(())


def test_service_rejects_non_workflow_registry_registry():
    with pytest.raises(ValueError, match="WorkflowRegistry"):
        ValidationService(registry=object())


def test_service_delegates_to_registered_validation_capable_adapter():
    adapter_result = Result.success({"sample_count": 1})
    adapter = FakeAdapter(result=adapter_result)
    inputs = WorkflowInputs(config={"samples": "samples.tsv"})
    service = ValidationService(registry=WorkflowRegistry(adapters=[adapter]))

    result = service.validate("fake", inputs)

    assert result is adapter_result
    assert adapter.validate_called
    assert adapter.received_inputs is inputs


def test_unknown_workflow_returns_workflow_not_found_issue():
    service = ValidationService(registry=WorkflowRegistry())

    result = service.validate("missing", WorkflowInputs(config={}))

    assert result.is_failure
    assert result.value is None
    assert result.errors == (
        Issue(
            code="WORKFLOW_NOT_FOUND",
            message="Workflow was not found.",
            source="registry",
            path="workflow_id",
            context={"workflow_id": "missing"},
        ),
    )


def test_invalid_workflow_id_from_registry_get_propagates_value_error():
    service = ValidationService(registry=WorkflowRegistry())

    with pytest.raises(ValueError):
        service.validate("", WorkflowInputs(config={}))


def test_adapter_without_validation_support_returns_capability_issue():
    adapter = FakeAdapter(supports=("dag_preview",))
    service = ValidationService(registry=WorkflowRegistry(adapters=[adapter]))

    result = service.validate("fake", WorkflowInputs(config={}))

    assert result.is_failure
    assert result.errors == (
        Issue(
            code="WORKFLOW_CAPABILITY_UNSUPPORTED",
            message="Workflow does not support validation.",
            source="registry",
            path="workflow.capabilities",
            context={"workflow_id": "fake", "capability": "validation"},
        ),
    )


def test_unsupported_capability_adapter_does_not_call_validate():
    adapter = FakeAdapter(supports=("dag_preview",))
    service = ValidationService(registry=WorkflowRegistry(adapters=[adapter]))

    service.validate("fake", WorkflowInputs(config={}))

    assert not adapter.validate_called


def test_unexpected_adapter_exception_propagates():
    adapter = FakeAdapter(exception=RuntimeError("adapter failed"))
    service = ValidationService(registry=WorkflowRegistry(adapters=[adapter]))

    with pytest.raises(RuntimeError, match="adapter failed"):
        service.validate("fake", WorkflowInputs(config={}))


def test_services_package_exports_validation_service():
    code = """
        from encode_pipeline.services import ValidationService
        print(ValidationService.__name__)
    """
    env = dict(os.environ)
    env["PYTHONPATH"] = str(SRC_ROOT)
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    proc = subprocess.run(
        [sys.executable, "-c", textwrap.dedent(code)],
        capture_output=True,
        text=True,
        env=env,
    )

    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.strip() == "ValidationService"


def test_importing_validation_service_does_not_import_workflow_specific_modules():
    code = """
        import sys
        import encode_pipeline.services.validation

        forbidden = [
            "encode_pipeline.adapters.encode",
            "encode_pipeline.config.validator",
            "encode_pipeline.samples",
            "encode_pipeline.cli.dag",
            "fastapi",
            "pydantic",
            "snakemake",
        ]
        for name in forbidden:
            print(f"{name}={name in sys.modules}")
    """
    env = dict(os.environ)
    env["PYTHONPATH"] = str(SRC_ROOT)
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    proc = subprocess.run(
        [sys.executable, "-c", textwrap.dedent(code)],
        capture_output=True,
        text=True,
        env=env,
    )

    assert proc.returncode == 0, proc.stderr
    assert set(proc.stdout.splitlines()) == {
        "encode_pipeline.adapters.encode=False",
        "encode_pipeline.config.validator=False",
        "encode_pipeline.samples=False",
        "encode_pipeline.cli.dag=False",
        "fastapi=False",
        "pydantic=False",
        "snakemake=False",
    }
