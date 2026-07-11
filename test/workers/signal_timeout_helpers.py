"""Importable helpers for real RQ signal-timeout integration tests."""

from __future__ import annotations

import os
from pathlib import Path
import time

from encode_pipeline.platform.results import Result
from encode_pipeline.services.local_execution import LocalExecutionService
from encode_pipeline.workers.jobs import run_execution_job


ENTERED_MARKER_ENV = "ENCODE_PIPELINE_TEST_TIMEOUT_ENTERED_MARKER"
CAUGHT_MARKER_ENV = "ENCODE_PIPELINE_TEST_TIMEOUT_CAUGHT_MARKER"


def run_execution_with_blocking_exception_handler(run_id: str) -> None:
    """Run the real durable job with an ``Exception``-swallowing test seam.

    The worker subprocess imports this helper. The temporary class replacement is
    isolated to its RQ work horse and lets the integration test prove that a real
    SIGALRM raises ``WorkerHardTimeout`` rather than the catchable native RQ
    timeout exception.
    """
    original_execute = LocalExecutionService.execute

    def block_until_rq_timeout(self, current_run_id: str):
        if current_run_id != run_id:
            raise RuntimeError("timeout helper received the wrong run identity")
        _marker_path(ENTERED_MARKER_ENV).write_text("entered\n", encoding="utf-8")
        try:
            time.sleep(60)
        except Exception:
            _marker_path(CAUGHT_MARKER_ENV).write_text("caught\n", encoding="utf-8")
            return Result.success(self._run_service.get_run(current_run_id))
        raise RuntimeError("timeout helper unexpectedly reached the sleep deadline")

    LocalExecutionService.execute = block_until_rq_timeout
    try:
        run_execution_job(run_id)
    finally:
        LocalExecutionService.execute = original_execute


def _marker_path(environment_name: str) -> Path:
    configured = os.environ.get(environment_name)
    if configured is None:
        raise RuntimeError("timeout integration marker is not configured")
    return Path(configured)
