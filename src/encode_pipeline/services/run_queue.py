"""Queue boundary for durable workflow execution."""

from __future__ import annotations

from typing import Protocol, runtime_checkable

from encode_pipeline.platform.execution import RunExecutionAssignment


class RunQueueError(RuntimeError):
    """Base error raised by durable execution queue adapters."""


class RunQueueIdentityError(RunQueueError):
    """A backend job does not match its durable run assignment."""


class RunQueueJobUnavailableError(RunQueueError):
    """A durable job identity cannot be reused for scheduling."""


class RunQueueUnavailableError(RunQueueError):
    """The queue backend could not confirm submission."""


class RunQueueStopUnavailableError(RunQueueError):
    """The queue backend could not confirm delivery of a stop command."""


@runtime_checkable
class RunQueue(Protocol):
    """Submit stable run identifiers to an out-of-process worker backend."""

    @property
    def backend(self) -> str:
        """Return the durable backend identity."""

    @property
    def queue_name(self) -> str:
        """Return the durable queue identity."""

    def enqueue_execution(self, assignment: RunExecutionAssignment) -> str:
        """Enqueue one execution and return the backend job identifier."""


@runtime_checkable
class RunStopQueue(Protocol):
    """Send a stop command for one strictly matched execution assignment."""

    @property
    def backend(self) -> str:
        """Return the durable backend identity."""

    @property
    def queue_name(self) -> str:
        """Return the durable queue identity."""

    def request_stop(self, assignment: RunExecutionAssignment) -> None:
        """Send the backend stop command or raise a bounded queue error."""
