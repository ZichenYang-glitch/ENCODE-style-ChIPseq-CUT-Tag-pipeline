"""Shared worker-boundary test helpers."""

from __future__ import annotations

from pathlib import Path

from encode_pipeline.persistence.runtime import open_run_persistence
from encode_pipeline.platform.adapters import WorkflowInputs
from encode_pipeline.platform.execution import RunExecutionAssignment
from encode_pipeline.platform.runs import RunStatus
from encode_pipeline.services.defaults import (
    create_default_workflow_registry,
)
from encode_pipeline.services.runs import RunService
from encode_pipeline.workers.settings import WorkerSettings


WORKFLOW_ID = "encode-style-chipseq-cuttag-atac-mnase"


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
) -> RunExecutionAssignment | None:
    """Persist a PLANNED run and optionally its canonical execution job ID."""
    persistence = open_run_persistence(settings.database_url)
    try:
        service = RunService(
            create_default_workflow_registry(),
            id_factory=lambda: run_id,
            repository=persistence.repository,
        )
        service.create_run(WORKFLOW_ID, WorkflowInputs(config={}))
        service.transition_run(run_id, RunStatus.VALIDATING, stage="preflight")
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
