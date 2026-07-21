"""Explicit result-capable DurableWorker bootstrap for the private gate."""

from __future__ import annotations

from pathlib import Path

import encode_pipeline.workers.jobs as worker_jobs
from encode_pipeline.workers.cli import main as worker_main
from encode_pipeline.workers.runtime import open_worker_runtime
from encode_pipeline.workers.settings import load_worker_settings
from encode_pipeline.workers.timeouts import WorkerHardTimeout

from .support import (
    build_acceptance_process_runner,
    build_results_composition,
    require_gate_settings,
)


def main() -> int:
    """Run one burst worker with the real result-capable deployment registry."""
    repository_root = Path(__file__).resolve().parents[2]
    gate_settings = require_gate_settings()
    composition = build_results_composition(
        gate_settings,
        project_root=repository_root,
    )
    worker_settings = load_worker_settings()
    original = worker_jobs.open_worker_runtime

    def open_results_runtime():
        return open_worker_runtime(
            worker_settings,
            registry=composition.registry,
            build_identity_provider=composition.build_identity_provider,
            process_runner=build_acceptance_process_runner(
                settings=gate_settings,
                binding=composition.binding,
                timeout_seconds=worker_settings.job_timeout_seconds,
                passthrough_exceptions=(WorkerHardTimeout,),
            ),
        )

    worker_jobs.open_worker_runtime = open_results_runtime
    try:
        return worker_main(("--burst",))
    finally:
        worker_jobs.open_worker_runtime = original


if __name__ == "__main__":
    raise SystemExit(main())
