"""Tests for worker-side dependency reconstruction."""

from __future__ import annotations

from pathlib import Path
import shutil
import traceback
from dataclasses import replace

import pytest

from encode_pipeline.persistence.runtime import open_run_persistence
from encode_pipeline.persistence.repositories import SqlAlchemyRunRepository
from encode_pipeline.persistence.input_registry import (
    SqlAlchemyInputRegistryRepository,
)
from encode_pipeline.services.local_run_driver import LocalRunDriver
from encode_pipeline.services.local_execution import LocalExecutionService
from encode_pipeline.services.managed_input_verification import (
    ManagedInputVerificationService,
)
from encode_pipeline.services.artifact_extraction import ArtifactExtractionService
from encode_pipeline.services.preflight import LocalPreflightService
from encode_pipeline.services.process_runner import ProcessRunner
from encode_pipeline.services.qc_summary_indexing import QcSummaryIndexingService
from encode_pipeline.services.reference_profile_runtime import (
    ReferenceProfileBindingService,
    ReferenceProfileRuntimeResolver,
)
from encode_pipeline.workers.runtime import WorkerRuntime, open_worker_runtime
from encode_pipeline.workers.timeouts import WorkerHardTimeout

from .conftest import WORKFLOW_ID, create_planned_run, worker_settings


def _prepare_database(configured) -> None:
    persistence = open_run_persistence(configured.database_url)
    persistence.close()


def test_open_worker_runtime_reopens_sqlite_and_full_execution_dependencies(tmp_path):
    configured = worker_settings(tmp_path)
    create_planned_run(configured, "runtime-run")

    with open_worker_runtime(configured) as runtime:
        record = runtime.run_service.get_run("runtime-run")

        assert runtime.settings is configured
        assert runtime.persistence.database_url == configured.database_url
        assert isinstance(runtime.persistence.repository, SqlAlchemyRunRepository)
        assert isinstance(
            runtime.persistence.input_registry_repository,
            SqlAlchemyInputRegistryRepository,
        )
        assert runtime.run_service._repository is runtime.persistence.repository
        assert record.workflow_id == WORKFLOW_ID
        assert runtime.registry.get(WORKFLOW_ID).metadata.workflow_id == WORKFLOW_ID
        assert isinstance(runtime.local_run_driver, LocalRunDriver)
        assert isinstance(runtime.local_execution_service, LocalExecutionService)
        assert isinstance(
            runtime.managed_input_verifier,
            ManagedInputVerificationService,
        )
        assert (
            runtime.managed_input_verifier._input_registry_repository
            is runtime.persistence.input_registry_repository
        )
        assert runtime.managed_input_verifier._storage_pool_config_path is None
        assert (
            runtime.local_execution_service._managed_input_verifier
            is runtime.managed_input_verifier
        )
        assert isinstance(
            runtime.artifact_extraction_service,
            ArtifactExtractionService,
        )
        assert isinstance(
            runtime.qc_summary_indexing_service,
            QcSummaryIndexingService,
        )
        assert isinstance(runtime.preflight_service, LocalPreflightService)
        assert isinstance(
            runtime.reference_profile_binding_service,
            ReferenceProfileBindingService,
        )
        assert isinstance(
            runtime.reference_profile_resolver,
            ReferenceProfileRuntimeResolver,
        )
        assert (
            runtime.workspace_planner._reference_profile_resolver
            is runtime.reference_profile_resolver
        )
        assert (
            runtime.command_builder._reference_profile_resolver
            is runtime.reference_profile_resolver
        )
        assert (
            runtime.artifact_extraction_service._reference_profile_resolver
            is runtime.reference_profile_resolver
        )
        assert (
            runtime.qc_summary_indexing_service._reference_profile_resolver
            is runtime.reference_profile_resolver
        )
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
    _prepare_database(configured)
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


def test_open_worker_runtime_uses_settings_admitted_encode_runtime_root(tmp_path):
    configured = worker_settings(tmp_path)
    _prepare_database(configured)
    runtime_root = (
        tmp_path
        / "runtimes"
        / "encode"
        / ("sha256-" + "a" * 64)
        / "payload"
        / "contracts"
        / "encode-runtime"
    ).resolve()
    inventory = runtime_root / "docs" / "architecture" / "artifact-inventory.yaml"
    inventory.parent.mkdir(parents=True)
    shutil.copy2(
        Path(__file__).resolve().parents[2]
        / "docs/architecture/artifact-inventory.yaml",
        inventory,
    )
    (runtime_root / "workflow").mkdir()
    (runtime_root / "workflow" / "Snakefile").write_text("rule all:\n    input: []\n")
    (runtime_root / "profiles" / "default").mkdir(parents=True)
    (runtime_root / "profiles" / "default" / "config.yaml").write_text("cores: 1\n")
    (runtime_root / "scripts").mkdir()
    (runtime_root / "scripts" / "runtime.py").write_text("# controlled\n")
    configured = replace(configured, encode_runtime_root=runtime_root)

    with open_worker_runtime(configured) as runtime:
        assert runtime.build_identity_provider.project_root == runtime_root
        assert runtime.command_builder._project_root == runtime_root


def test_open_worker_runtime_accepts_only_deployment_owned_registry_and_runner(
    tmp_path,
):
    from encode_pipeline.services.defaults import create_default_workflow_registry

    configured = worker_settings(tmp_path)
    _prepare_database(configured)
    registry = create_default_workflow_registry()
    runner = ProcessRunner(
        allowed_executables=("/opt/helixweave/nextflow",),
        timeout_seconds=17,
    )

    with open_worker_runtime(
        configured,
        registry=registry,
        process_runner=runner,
    ) as runtime:
        assert runtime.registry is registry
        assert runtime.local_run_driver._process_runner is runner


def test_open_worker_runtime_uses_provider_registry_and_rejects_mismatch(tmp_path):
    from encode_pipeline.services.workflow_builds import (
        WorkflowBuildIdentityProvider,
    )
    from encode_pipeline.platform.registry import WorkflowRegistry

    configured = worker_settings(tmp_path)
    _prepare_database(configured)
    selected = WorkflowRegistry()
    provider = WorkflowBuildIdentityProvider(
        selected,
        project_root=Path(__file__).resolve().parents[2],
    )

    with open_worker_runtime(
        configured,
        build_identity_provider=provider,
    ) as runtime:
        assert runtime.registry is selected
        assert runtime.build_identity_provider is provider

    mismatched = WorkflowRegistry()
    try:
        open_worker_runtime(
            configured,
            registry=mismatched,
            build_identity_provider=provider,
        )
    except ValueError as exc:
        assert str(exc) == "registry must be the build_identity_provider registry"
    else:  # pragma: no cover - protects the identity/execution binding invariant
        raise AssertionError("mismatched registry/provider unexpectedly accepted")


def test_open_worker_runtime_rejects_invalid_composition_overrides(tmp_path):
    configured = worker_settings(tmp_path)

    for kwargs in ({"registry": object()}, {"process_runner": object()}):
        try:
            open_worker_runtime(configured, **kwargs)
        except ValueError:
            pass
        else:  # pragma: no cover - guards the deployment-owned seam
            raise AssertionError("invalid worker composition unexpectedly succeeded")


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
        "encode_pipeline.workers.runtime.open_existing_run_persistence",
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


@pytest.mark.parametrize(
    "close_error",
    (
        RuntimeError("/private/database"),
        WorkerHardTimeout("second deadline during close"),
    ),
)
def test_runtime_close_cannot_replace_active_worker_hard_timeout(close_error):
    class RuntimeWithBrokenClose:
        __enter__ = WorkerRuntime.__enter__
        __exit__ = WorkerRuntime.__exit__

        def close(self):
            raise close_error

    with pytest.raises(WorkerHardTimeout, match="original RQ deadline") as raised:
        with RuntimeWithBrokenClose():
            raise WorkerHardTimeout("original RQ deadline")

    assert raised.value is not close_error
    formatted = "".join(traceback.format_exception(raised.value))
    assert "/private/database" not in formatted
    assert "second deadline during close" not in formatted


def test_runtime_close_failure_without_active_timeout_still_propagates():
    class RuntimeWithBrokenClose:
        __enter__ = WorkerRuntime.__enter__
        __exit__ = WorkerRuntime.__exit__

        def close(self):
            raise RuntimeError("close failed")

    with pytest.raises(RuntimeError, match="close failed"):
        with RuntimeWithBrokenClose():
            pass


def test_runtime_close_hard_timeout_suppresses_private_body_exception():
    class RuntimeWithTimedOutClose:
        __enter__ = WorkerRuntime.__enter__
        __exit__ = WorkerRuntime.__exit__

        def close(self):
            raise WorkerHardTimeout("deadline during runtime close")

    with pytest.raises(
        WorkerHardTimeout,
        match="deadline during runtime close",
    ) as raised:
        with RuntimeWithTimedOutClose():
            raise RuntimeError("/private/workspace")

    formatted = "".join(traceback.format_exception(raised.value))
    assert "/private/workspace" not in formatted
