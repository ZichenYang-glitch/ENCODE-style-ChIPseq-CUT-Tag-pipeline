"""Tests for workflow-platform adapter registry primitives."""

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
from encode_pipeline.platform.results import Result


SRC_ROOT = Path(__file__).resolve().parents[2] / "src"


class FakeAdapter:
    def __init__(self, workflow_id: str, name: str = "Fake Workflow") -> None:
        self.metadata = WorkflowMetadata(
            workflow_id=workflow_id,
            name=name,
            version="1.0.0",
        )
        self.capabilities = WorkflowCapabilities(supports=("validation",))

    def schema(self) -> WorkflowSchema:
        return WorkflowSchema(config_schema={"type": "object"})

    def validate(self, inputs: WorkflowInputs) -> Result[object]:
        return Result.success({"validated": dict(inputs.config)})

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


def test_fake_adapter_can_be_registered_and_resolved_by_workflow_id():
    adapter = FakeAdapter("fake")
    registry = WorkflowRegistry(adapters=[adapter])

    assert registry.get("fake") is adapter


def test_has_returns_true_or_false_for_valid_workflow_ids():
    registry = WorkflowRegistry(adapters=[FakeAdapter("fake")])

    assert registry.has("fake")
    assert not registry.has("missing")


def test_list_metadata_returns_metadata_tuple_in_registration_order():
    first = FakeAdapter("first", name="First Workflow")
    second = FakeAdapter("second", name="Second Workflow")
    registry = WorkflowRegistry(adapters=[first, second])

    assert registry.list_metadata() == (first.metadata, second.metadata)


def test_duplicate_workflow_id_raises_value_error():
    with pytest.raises(ValueError, match="Duplicate workflow_id"):
        WorkflowRegistry(adapters=[FakeAdapter("fake"), FakeAdapter("fake")])


def test_missing_workflow_id_raises_key_error():
    registry = WorkflowRegistry(adapters=[FakeAdapter("fake")])

    with pytest.raises(KeyError):
        registry.get("missing")


@pytest.mark.parametrize("workflow_id", ["", "   ", None, 123])
def test_invalid_lookup_ids_raise_value_error(workflow_id):
    registry = WorkflowRegistry(adapters=[FakeAdapter("fake")])

    with pytest.raises(ValueError):
        registry.get(workflow_id)
    with pytest.raises(ValueError):
        registry.has(workflow_id)


def test_non_adapter_entries_raise_value_error():
    with pytest.raises(ValueError, match="WorkflowRegistry adapters"):
        WorkflowRegistry(adapters=[object()])


def test_registry_does_not_expose_mutable_adapter_mapping():
    registry = WorkflowRegistry(adapters=[FakeAdapter("fake")])

    assert not hasattr(registry, "adapters")
    assert not hasattr(registry, "_adapters_by_id")


def test_platform_exports_workflow_registry():
    code = """
        from encode_pipeline.platform import WorkflowRegistry
        print(WorkflowRegistry.__name__)
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
    assert proc.stdout.strip() == "WorkflowRegistry"


def test_importing_platform_registry_does_not_import_workflow_specific_modules():
    code = """
        import sys
        import encode_pipeline.platform.registry

        forbidden = [
            "encode_pipeline.adapters.encode",
            "encode_pipeline.config.validator",
            "encode_pipeline.samples",
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
        "fastapi=False",
        "pydantic=False",
        "snakemake=False",
    }
