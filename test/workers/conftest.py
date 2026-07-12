"""Shared worker-boundary test helpers."""

from __future__ import annotations

from pathlib import Path

import yaml

from encode_pipeline.persistence.runtime import open_run_persistence
from encode_pipeline.platform.adapters import WorkflowInputs
from encode_pipeline.platform.builds import WorkflowBuildIdentity
from encode_pipeline.platform.execution import RunExecutionAssignment
from encode_pipeline.platform.runs import RunStatus
from encode_pipeline.services.defaults import (
    create_default_workflow_build_identity_provider,
    create_default_workflow_registry,
)
from encode_pipeline.services.materialization import WorkspaceMaterializer
from encode_pipeline.services.planning import ExecutionPlanner, WorkspacePlanner
from encode_pipeline.services.runs import RunService
from encode_pipeline.workers.settings import WorkerSettings


WORKFLOW_ID = "encode-style-chipseq-cuttag-atac-mnase"
PROFILE_ROOT = Path(__file__).resolve().parents[1] / "profiles" / "platform_worker_tiny"


def worker_settings(tmp_path: Path, queue_name: str = "worker-tests") -> WorkerSettings:
    """Return isolated file-backed settings for one worker test."""
    return WorkerSettings(
        database_url=f"sqlite:///{tmp_path / 'platform.db'}",
        redis_url="redis://unused.test/0",
        queue_name=queue_name,
        workspace_root=tmp_path / "workspaces",
    )


def create_planned_run(
    settings: WorkerSettings,
    run_id: str,
    *,
    assign_queue: str | None = None,
    bind_build_identity: bool = True,
    build_identity: WorkflowBuildIdentity | None = None,
    enable_qc_summary: bool = False,
    samples_path: Path | None = None,
) -> RunExecutionAssignment | None:
    """Persist a PLANNED run and optionally its canonical execution job ID."""
    persistence = open_run_persistence(settings.database_url)
    try:
        registry = create_default_workflow_registry()
        service = RunService(
            registry,
            id_factory=lambda: run_id,
            repository=persistence.repository,
        )
        config = yaml.safe_load(
            (PROFILE_ROOT / "config.yaml").read_text(encoding="utf-8")
        )
        configured_samples = (
            (PROFILE_ROOT / "samples.tsv").resolve()
            if samples_path is None
            else samples_path.resolve()
        )
        config["samples"] = str(configured_samples)
        if enable_qc_summary:
            config["qc"]["summary"] = True
        service.create_run(WORKFLOW_ID, WorkflowInputs(config=config))
        service.transition_run(run_id, RunStatus.VALIDATING, stage="preflight")
        base_result = ExecutionPlanner(service).plan_run(run_id)
        assert base_result.is_success
        base_plan = base_result.value
        workspace_dir = settings.workspace_root / run_id
        workspace_result = WorkspacePlanner(registry).plan_workspace(
            base_plan,
            workspace_dir,
        )
        assert workspace_result.is_success
        workspace_plan = workspace_result.value
        materialized = WorkspaceMaterializer().materialize(
            workspace_plan.workspace_plan,
            workspace_dir,
        )
        assert materialized.is_success
        if bind_build_identity:
            if build_identity is None:
                identity_result = create_default_workflow_build_identity_provider(
                    registry=registry
                ).capture(WORKFLOW_ID)
                assert identity_result.is_success
                build_identity = identity_result.value
            service.complete_preflight(
                run_id,
                build_identity,
                stage="preflight",
            )
        else:
            service.transition_run(run_id, RunStatus.PLANNED, stage="preflight")
        if assign_queue is None:
            return None
        assignment = service.ensure_execution_assignment(
            run_id,
            queue_name=assign_queue,
        )
        return assignment
    finally:
        persistence.close()
