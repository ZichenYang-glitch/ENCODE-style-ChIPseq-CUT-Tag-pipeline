"""Tests for worker-side dependency reconstruction."""

from __future__ import annotations

from encode_pipeline.services.local_run_driver import LocalRunDriver
from encode_pipeline.services.preflight import LocalPreflightService
from encode_pipeline.workers.runtime import open_worker_runtime

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
        assert isinstance(runtime.preflight_service, LocalPreflightService)
        assert runtime.local_run_driver._workspace_root == configured.workspace_root


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
        lambda: (_ for _ in ()).throw(RuntimeError("composition failed")),
    )

    try:
        open_worker_runtime(configured)
    except RuntimeError as exc:
        assert str(exc) == "composition failed"
    else:  # pragma: no cover - protects the resource cleanup assertion
        raise AssertionError("composition unexpectedly succeeded")

    assert closed is True
