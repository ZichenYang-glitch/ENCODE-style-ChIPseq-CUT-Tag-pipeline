"""Tests for default workflow service composition."""

import asyncio
import os
import subprocess
import sys
import textwrap
from pathlib import Path

from encode_pipeline.adapters.encode import EncodeStyleWorkflowAdapter
from encode_pipeline.api.models import AgentRequest
from encode_pipeline.platform.adapters import WorkflowInputs
from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.platform.results import Issue
from encode_pipeline.services.runs import RunService
from encode_pipeline.services.validation import ValidationService

from encode_pipeline.services.defaults import (
    create_default_agent_service,
    create_default_validation_service,
    create_default_workflow_registry,
)


SRC_ROOT = Path(__file__).resolve().parents[2] / "src"
WORKFLOW_ID = "encode-style-chipseq-cuttag-atac-mnase"
VALID_SAMPLES = (
    "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\tpeak_mode\tgenome\tbowtie2_index\n"
    "S1\tR1.fq\tR2.fq\tPE\tchipseq\tCTCF\tnarrow\ths\tidx\n"
)


def _write_samples(path: Path) -> Path:
    path.write_text(VALID_SAMPLES, encoding="utf-8")
    return path


def test_create_default_workflow_registry_returns_workflow_registry():
    registry = create_default_workflow_registry()

    assert isinstance(registry, WorkflowRegistry)


def test_default_registry_lists_encode_metadata():
    registry = create_default_workflow_registry()

    metadata = registry.list_metadata()

    assert len(metadata) == 1
    assert metadata[0].workflow_id == WORKFLOW_ID
    assert metadata[0].name == "ENCODE-style ChIP-seq/CUT&Tag/ATAC/MNase"


def test_default_registry_resolves_encode_adapter():
    registry = create_default_workflow_registry()

    adapter = registry.get(WORKFLOW_ID)

    assert isinstance(adapter, EncodeStyleWorkflowAdapter)


def test_repeated_default_registry_calls_return_distinct_objects():
    first = create_default_workflow_registry()
    second = create_default_workflow_registry()

    assert first is not second


def test_repeated_default_registry_calls_return_distinct_adapter_instances():
    first = create_default_workflow_registry()
    second = create_default_workflow_registry()

    assert first.get(WORKFLOW_ID) is not second.get(WORKFLOW_ID)


def test_create_default_validation_service_returns_validation_service():
    service = create_default_validation_service()

    assert isinstance(service, ValidationService)


def test_default_validation_service_validates_encode_inputs(tmp_path):
    samples_path = _write_samples(tmp_path / "samples.tsv")
    service = create_default_validation_service()

    result = service.validate(
        WORKFLOW_ID,
        WorkflowInputs(config={"samples": str(samples_path), "use_control": False}),
    )

    assert result.is_success
    assert isinstance(result.value, dict)
    assert set(result.value) == {"config", "samples"}
    assert result.value["config"]["samples"] == str(samples_path)
    assert result.value["samples"][0]["id"] == "S1"


def test_unknown_workflow_through_default_service_returns_workflow_not_found():
    service = create_default_validation_service()

    result = service.validate("missing", WorkflowInputs(config={}))

    assert result.is_failure
    assert result.errors == (
        Issue(
            code="WORKFLOW_NOT_FOUND",
            message="Workflow was not found.",
            source="registry",
            path="workflow_id",
            context={"workflow_id": "missing"},
        ),
    )


def test_create_default_agent_service_returns_agent_service():
    service = create_default_agent_service()

    assert service.__class__.__name__ == "AgentService"


def test_create_default_agent_service_wires_safety_components():
    service = create_default_agent_service()
    request = AgentRequest(message="Explain.")

    response = asyncio.run(service.chat(WORKFLOW_ID, request))

    assert response.ok is True
    assert "mock explanation" in response.message.lower()


def test_create_default_agent_service_redacts_paths_in_response():
    service = create_default_agent_service()
    request = AgentRequest(message="Where is my data?")

    response = asyncio.run(service.chat(WORKFLOW_ID, request))

    assert "/home" not in response.message.lower()


def test_services_package_exports_default_factory_functions():
    code = """
        from encode_pipeline.services import (
            create_default_agent_service,
            create_default_validation_service,
            create_default_workflow_registry,
        )
        print(create_default_workflow_registry.__name__)
        print(create_default_validation_service.__name__)
        print(create_default_agent_service.__name__)
    """
    proc = _run_python(code)

    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.splitlines() == [
        "create_default_workflow_registry",
        "create_default_validation_service",
        "create_default_agent_service",
    ]


def test_importing_services_defaults_does_not_import_heavy_workflow_modules():
    code = """
        import sys
        import encode_pipeline.services.defaults

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
    proc = _run_python(code)

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


def test_importing_services_package_does_not_import_heavy_workflow_modules():
    code = """
        import sys
        import encode_pipeline.services

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
    proc = _run_python(code)

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


def test_calling_default_registry_factory_intentionally_imports_encode_adapter():
    code = """
        import sys
        from encode_pipeline.services.defaults import create_default_workflow_registry

        print(f"before={ 'encode_pipeline.adapters.encode' in sys.modules }")
        create_default_workflow_registry()
        print(f"after={ 'encode_pipeline.adapters.encode' in sys.modules }")
    """
    proc = _run_python(code)

    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.splitlines() == [
        "before=False",
        "after=True",
    ]


def test_create_default_local_run_driver_returns_instance():
    from encode_pipeline.services import create_default_local_run_driver

    registry = create_default_workflow_registry()
    run_service = RunService(registry=registry)
    driver = create_default_local_run_driver(run_service=run_service)

    assert driver.__class__.__name__ == "LocalRunDriver"


def test_create_default_workspace_materializer_returns_instance():
    from encode_pipeline.services import create_default_workspace_materializer

    materializer = create_default_workspace_materializer()
    assert materializer.__class__.__name__ == "WorkspaceMaterializer"


def test_create_default_command_builder_returns_instance():
    from encode_pipeline.services import create_default_command_builder

    builder = create_default_command_builder()

    assert builder.__class__.__name__ == "CommandBuilder"


def test_create_default_command_builder_forwards_project_root(tmp_path):
    from encode_pipeline.services import create_default_command_builder

    project_root = (tmp_path / "source").resolve()
    builder = create_default_command_builder(project_root=project_root)

    assert builder._project_root == project_root


def _run_python(code: str) -> subprocess.CompletedProcess[str]:
    env = dict(os.environ)
    env["PYTHONPATH"] = str(SRC_ROOT)
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    return subprocess.run(
        [sys.executable, "-c", textwrap.dedent(code)],
        capture_output=True,
        text=True,
        env=env,
    )
