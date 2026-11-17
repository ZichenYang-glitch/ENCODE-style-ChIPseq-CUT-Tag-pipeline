"""The Peri session factory used by the single authoring owner.

Platform grants, capability tokens, and generation fences stay on the owner.
This module only opens the fixed Peri process. ENCODE field projection and
the repository protocol remain P1b and are not part of this seam.

A substitute runtime can prove that the owner calls this seam. It does not
qualify another agent implementation.
"""

from __future__ import annotations

from pathlib import Path
from typing import Protocol

from collections.abc import Callable

from encode_pipeline.services.peri_runtime.session import NativeSession
from encode_pipeline.services.peri_runtime.stub import (
    open_qualification_model_transport,
)


ToolHandler = Callable[..., dict[str, object]]


class RuntimeSession(Protocol):
    """The session surface the owner drives. It does not hold platform grants."""

    instance: str
    service_generation: int
    state: str
    handler_calls: int
    directory: Path
    tool_arrivals: int
    errors: list[str]
    approvals: list[dict[str, object]]
    denied_callbacks: list[str]
    authorization_check: Callable[[], bool] | None

    def initialize(self) -> None:
        """Qualify the process before the owner marks it ready."""

    def run_authoring_turn(self, final_text: str) -> str:
        """Run one turn. Egress still asks the owner's callback."""

    def close(self) -> None:
        """Stop the process. A second call is harmless."""


class AuthoringRuntime(Protocol):
    """Open one runtime session for a generation the owner has already fenced."""

    def open_session(
        self,
        binary: Path,
        directory: Path,
        *,
        service_generation: int,
        tool_handler: ToolHandler,
    ) -> RuntimeSession:
        """Start one session. The owner retains its platform grant."""


class PeriAuthoringRuntime:
    """The accepted Peri implementation and its labeled qualification stub."""

    def open_session(
        self,
        binary: Path,
        directory: Path,
        *,
        service_generation: int,
        tool_handler: ToolHandler,
    ) -> RuntimeSession:
        return NativeSession(
            binary,
            directory,
            service_generation=service_generation,
            tool_handler=tool_handler,
            model_transport=open_qualification_model_transport,
        )
