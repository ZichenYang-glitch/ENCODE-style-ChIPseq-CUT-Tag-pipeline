"""Tests for default workflow service composition."""

import asyncio
import os
import stat
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

from encode_pipeline.adapters.encode import EncodeStyleWorkflowAdapter
from encode_pipeline.adapters.bulk_rnaseq import (
    BulkRnaSeqExecutionBinding,
    BulkRnaSeqResultsWorkflowAdapter,
    BulkRnaSeqTranscriptomeBinding,
    BulkRnaSeqWorkflowAdapter,
    RuntimeAssetBinding,
)
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
    create_default_process_runner,
)
from encode_pipeline.workers.settings import WorkerSettings


SRC_ROOT = Path(__file__).resolve().parents[2] / "src"
WORKFLOW_ID = "encode-style-chipseq-cuttag-atac-mnase"
BULK_WORKFLOW_ID = "bulk-rnaseq"


def _write_samples(path: Path) -> Path:
    root = path.parent.resolve()
    path.write_text(
        "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\tpeak_mode\tgenome\tbowtie2_index\n"
        f"S1\t{root / 'R1.fq'}\t{root / 'R2.fq'}\tPE\tchipseq\tCTCF\tnarrow\ths\t{root / 'idx'}\n",
        encoding="utf-8",
    )
    return path


def test_create_default_workflow_registry_returns_workflow_registry():
    registry = create_default_workflow_registry()

    assert isinstance(registry, WorkflowRegistry)


def test_default_registry_lists_encode_and_bulk_metadata_without_runtime():
    registry = create_default_workflow_registry(environ={})

    metadata = registry.list_metadata()

    assert len(metadata) == 2
    assert metadata[0].workflow_id == WORKFLOW_ID
    assert metadata[0].name == "ENCODE-style ChIP-seq/CUT&Tag/ATAC/MNase"
    assert metadata[1].workflow_id == BULK_WORKFLOW_ID
    assert metadata[1].name == "Bulk RNA-seq"


def test_default_registry_resolves_encode_adapter():
    registry = create_default_workflow_registry()

    adapter = registry.get(WORKFLOW_ID)

    assert isinstance(adapter, EncodeStyleWorkflowAdapter)


def test_default_registry_resolves_authoring_only_bulk_adapter_without_runtime():
    registry = create_default_workflow_registry(environ={})

    adapter = registry.get(BULK_WORKFLOW_ID)

    assert isinstance(adapter, BulkRnaSeqWorkflowAdapter)
    assert adapter.capabilities.supports == ("validation", "input_authoring")
    availability = adapter.execution_availability()
    assert availability.execution == "not_configured"
    assert availability.reason_code == "WORKFLOW_EXECUTION_NOT_CONFIGURED"


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


def test_default_process_runner_allows_only_bundled_engine_without_bulk_runtime(
    tmp_path,
):
    settings = WorkerSettings(
        database_url=f"sqlite:///{tmp_path / 'db.sqlite'}",
        redis_url="redis://localhost:6379/0",
        queue_name="tests",
        workspace_root=(tmp_path / "workspaces").resolve(),
    )

    runner = create_default_process_runner(
        registry=create_default_workflow_registry(environ={}),
        settings=settings,
    )

    assert runner._allowed_executables == ("snakemake",)


def test_default_process_runner_uses_exact_admitted_bulk_server_coordinates(
    tmp_path,
    monkeypatch,
):
    import encode_pipeline.services.managed_containers as managed_containers_module

    docker = (tmp_path / "bin/docker").resolve()
    socket = (tmp_path / "docker.sock").resolve()
    docker.parent.mkdir()
    docker.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
    docker.chmod(0o755)
    socket.write_bytes(b"")
    real_lstat = os.lstat
    socket_stat = real_lstat(socket)
    socket_values = list(socket_stat)
    socket_values[0] = stat.S_IFSOCK | 0o660
    simulated_socket_stat = os.stat_result(socket_values)

    def lstat_with_simulated_socket(path):
        if Path(path) == socket:
            return simulated_socket_stat
        return real_lstat(path)

    monkeypatch.setattr(
        managed_containers_module.os,
        "lstat",
        lstat_with_simulated_socket,
    )
    binding = BulkRnaSeqExecutionBinding(
        assets=RuntimeAssetBinding(
            root=(tmp_path / "runtime").resolve(),
            docker_executable=docker,
            docker_socket=socket,
        ),
        transcriptome=BulkRnaSeqTranscriptomeBinding(
            reference_id="tiny",
            fasta_sha256="a" * 64,
            gtf_sha256="b" * 64,
            transcript_fasta=(tmp_path / "transcripts.fa").resolve(),
            transcript_fasta_sha256="c" * 64,
        ),
    )
    registry = WorkflowRegistry(
        [
            EncodeStyleWorkflowAdapter(),
            BulkRnaSeqResultsWorkflowAdapter(execution=binding),
        ]
    )
    settings = WorkerSettings(
        database_url=f"sqlite:///{tmp_path / 'db.sqlite'}",
        redis_url="redis://localhost:6379/0",
        queue_name="tests",
        workspace_root=(tmp_path / "workspaces").resolve(),
        managed_docker_executable=docker,
        managed_docker_socket=socket,
    )

    runner = create_default_process_runner(registry=registry, settings=settings)

    assert runner._allowed_executables == ("snakemake", "/usr/bin/unshare")

    mismatched = WorkerSettings(
        database_url=settings.database_url,
        redis_url="redis://localhost:6379/0",
        queue_name="tests",
        workspace_root=settings.workspace_root,
        managed_docker_executable=(tmp_path / "other/docker").resolve(),
        managed_docker_socket=socket,
    )
    with pytest.raises(ValueError, match="do not match"):
        create_default_process_runner(registry=registry, settings=mismatched)


def test_default_process_runner_disables_bulk_when_cleaner_cannot_bind(
    tmp_path,
    monkeypatch,
):
    import encode_pipeline.services.defaults as defaults

    docker = (tmp_path / "bin/docker").resolve()
    socket = (tmp_path / "missing/docker.sock").resolve()
    binding = BulkRnaSeqExecutionBinding(
        assets=RuntimeAssetBinding(
            root=(tmp_path / "runtime").resolve(),
            docker_executable=docker,
            docker_socket=socket,
        ),
        transcriptome=BulkRnaSeqTranscriptomeBinding(
            reference_id="tiny",
            fasta_sha256="a" * 64,
            gtf_sha256="b" * 64,
            transcript_fasta=(tmp_path / "transcripts.fa").resolve(),
            transcript_fasta_sha256="c" * 64,
        ),
    )
    adapter = BulkRnaSeqResultsWorkflowAdapter(execution=binding)
    registry = WorkflowRegistry([EncodeStyleWorkflowAdapter(), adapter])
    settings = WorkerSettings(
        database_url=f"sqlite:///{tmp_path / 'db.sqlite'}",
        redis_url="redis://localhost:6379/0",
        queue_name="tests",
        workspace_root=(tmp_path / "workspaces").resolve(),
        managed_docker_executable=docker,
        managed_docker_socket=socket,
    )

    def unavailable_cleaner(_settings):
        raise FileNotFoundError("private Docker endpoint")

    monkeypatch.setattr(
        defaults,
        "create_default_managed_container_cleaner",
        unavailable_cleaner,
    )

    runner = create_default_process_runner(registry=registry, settings=settings)

    assert runner._allowed_executables == ("snakemake",)
    assert runner._managed_container_cleaner is None
    assert adapter.execution_binding is None
    assert adapter.execution_availability().execution == "unavailable"


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
