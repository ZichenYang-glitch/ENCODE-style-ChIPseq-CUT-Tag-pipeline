"""Queue boundary for durable workflow execution."""

from __future__ import annotations

from typing import Protocol, runtime_checkable

from encode_pipeline.platform.execution import RunExecutionAssignment


@runtime_checkable
class RunQueue(Protocol):
    """Submit stable run identifiers to an out-of-process worker backend."""

    def enqueue_execution(self, assignment: RunExecutionAssignment) -> str:
        """Enqueue one execution and return the backend job identifier."""
