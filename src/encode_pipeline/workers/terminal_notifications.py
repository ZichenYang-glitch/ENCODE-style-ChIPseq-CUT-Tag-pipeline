"""Worker-only isolation for terminal notification deadline interrupts."""

from __future__ import annotations

from dataclasses import dataclass

from encode_pipeline.platform.notifications import TerminalRunNotifier
from encode_pipeline.platform.runs import RunStatus
from encode_pipeline.workers.timeouts import WorkerHardTimeout


@dataclass(frozen=True)
class WorkerTerminalRunNotifier:
    """Prevent best-effort SMTP from changing an RQ job's durable outcome."""

    delegate: TerminalRunNotifier

    def __post_init__(self) -> None:
        if not callable(getattr(self.delegate, "notify_terminal_run", None)):
            raise ValueError("delegate must implement TerminalRunNotifier")

    def notify_terminal_run(
        self,
        run_id: str,
        status: RunStatus,
        *,
        include_qc: bool = False,
    ) -> None:
        try:
            self.delegate.notify_terminal_run(
                run_id,
                status,
                include_qc=include_qc,
            )
        except WorkerHardTimeout:
            return
