"""Tests for worker-side dependency reconstruction."""

from __future__ import annotations

from pathlib import Path
import shutil

from encode_pipeline.services.local_run_driver import LocalRunDriver
from encode_pipeline.services.local_execution import LocalExecutionService
from encode_pipeline.services.artifact_extraction import ArtifactExtractionService
from encode_pipeline.services.preflight import LocalPreflightService
from encode_pipeline.services.qc_summary_indexing import QcSummaryIndexingService
from encode_pipeline.workers.runtime import open_worker_runtime
from encode_pipeline.workers.timeouts import WorkerHardTimeout

from .conftest import WORKFLOW_ID, create_planned_run, worker_settings


def test_open_worker_runtime_reopens_sqlite_and_full_execution_dependencies(tmp_path):
    configured = worker_settings(tmp_path)
    create_planned_run(configured, "runtime-run")

    with open_worker_runtime(configured) as runtime:
        record = runtime.run_service.get_run("runtime-run")

        assert runtime.settings is configured
        assert runtime.persistence.database_url == configured.database_url
        assert record.workflow_id == WORKFLOW_ID
        assert runtime.registry.get(WORKFLOW_ID).metadata.workflow_id == WORKFLOW_ID
        assert isinstance(runtime.local_run_driver, LocalRunDriver)
        assert isinstance(runtime.local_execution_service, LocalExecutionService)
        assert isinstance(
            runtime.artifact_extraction_service,
            ArtifactExtractionService,
        )
        assert isinstance(
            runtime.qc_summary_indexing_service,
            QcSummaryIndexingService,
        )
        assert isinstance(runtime.preflight_service, LocalPreflightService)
        assert runtime.local_run_driver._workspace_root == configured.workspace_root
        process_timeout = runtime.local_run_driver._process_runner._timeout_seconds
        assert process_timeout == configured.job_timeout_seconds
        assert runtime.local_run_driver._process_runner._passthrough_exceptions == (
            WorkerHardTimeout,
        )


def test_open_worker_runtime_aligns_command_and_identity_project_roots(tmp_path):
    from encode_pipeline.services.defaults import (
        create_default_workflow_build_identity_provider,
    )

    configured = worker_settings(tmp_path)
    project_root = (tmp_path / "source").resolve()
    inventory = project_root / "docs/architecture/artifact-inventory.yaml"
    inventory.parent.mkdir(parents=True)
    shutil.copy2(
        Path(__file__).resolve().parents[2]
        / "docs/architecture/artifact-inventory.yaml",
        inventory,
    )
    provider = create_default_workflow_build_identity_provider(
        project_root=project_root,
    )

    with open_worker_runtime(
        configured,
        build_identity_provider=provider,
    ) as runtime:
        assert runtime.build_identity_provider is provider
        assert runtime.command_builder._project_root == project_root
        assert runtime.local_run_driver._command_builder is runtime.command_builder


def test_open_worker_runtime_closes_persistence_if_composition_fails(
    tmp_path, monkeypatch
):
    configured = worker_settings(tmp_path)
    closed = False

    class BrokenPersistence:
        repository = object()

        def close(self):
            nonlocal closed
            closed = True

    monkeypatch.setattr(
        "encode_pipeline.workers.runtime.open_run_persistence",
        lambda _database_url: BrokenPersistence(),
    )
    monkeypatch.setattr(
        "encode_pipeline.services.defaults.create_default_workflow_registry",
        lambda **_kwargs: (_ for _ in ()).throw(RuntimeError("composition failed")),
    )

    try:
        open_worker_runtime(configured)
    except RuntimeError as exc:
        assert str(exc) == "composition failed"
    else:  # pragma: no cover - protects the resource cleanup assertion
        raise AssertionError("composition unexpectedly succeeded")

    assert closed is True
