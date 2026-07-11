"""Module-level RQ job entry points."""

from __future__ import annotations

from rq import get_current_job

from encode_pipeline.services.run_repositories import ConcurrentRunUpdateError
from encode_pipeline.workers.runtime import open_worker_runtime


class WorkerJobIdentityError(RuntimeError):
    """The dequeued RQ job does not own the durable execution assignment."""


def run_execution_job(run_id: str) -> None:
    """Rebuild dependencies and prove ownership of a durable assignment.

    PR123 deliberately stops at the worker-boundary handshake. Real Snakemake
    execution and lifecycle transitions are introduced in PR124.
    """
    if not isinstance(run_id, str) or not run_id.strip():
        raise WorkerJobIdentityError("run_id must be a non-empty string")
    run_id = run_id.strip()

    current_job = get_current_job()
    if current_job is None or not current_job.id:
        raise WorkerJobIdentityError(
            "run_execution_job must run inside an RQ job context"
        )

    with open_worker_runtime() as runtime:
        try:
            record = runtime.run_service.get_run(run_id)
            runtime.registry.get(record.workflow_id)
        except KeyError:
            raise WorkerJobIdentityError(
                "execution job does not reference durable workflow state"
            ) from None
        _require_identity(
            "RQ job origin",
            getattr(current_job, "origin", None),
            runtime.settings.queue_name,
        )
        try:
            claim = runtime.run_service.claim_execution_assignment(
                run_id,
                job_id=current_job.id,
                backend="rq",
                queue_name=runtime.settings.queue_name,
                context={"workflow_id": record.workflow_id},
            )
        except (ConcurrentRunUpdateError, KeyError, ValueError):
            raise WorkerJobIdentityError(
                "stale execution job does not own durable workflow state"
            ) from None
        if not claim.acquired:
            return


def _require_identity(label: str, actual: object, expected: object) -> None:
    if actual != expected:
        raise WorkerJobIdentityError(
            f"stale execution job: {label} does not match durable ownership"
        )
