"""In-process port for the restricted authoring bridge.

The default implementation admits nothing. A qualified runtime is installed by
the bridge owner; manual draft editing does not depend on it.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

from encode_pipeline.platform.authoring import BRIDGE_UNAVAILABLE


@dataclass(frozen=True)
class BridgeStatus:
    available: bool
    reason_code: str | None = None


@dataclass(frozen=True)
class BridgeAdmission:
    instance_id: str
    state: str


@dataclass(frozen=True)
class BridgeSessionRequest:
    sid: str
    generation: int
    owner_epoch: str
    capability_token: str
    schema_digest: str
    projection_digest: str
    auth_session_digest: str = ""
    user_id: str = ""
    instance_id: str | None = None


@dataclass(frozen=True)
class BridgeTurnRequest:
    sid: str
    turn_id: str
    generation: int
    instance_id: str
    owner_epoch: str
    capability_token: str
    question: dict[str, object]
    issue_codes: tuple[str, ...]
    schema_digest: str
    projection_digest: str
    base_revision: int


class BridgeDisconnected(Exception):
    """The native session closed before producing a final answer."""


class AuthoringBridge(Protocol):
    """One bridge owner. Each active authoring session has one native."""

    def status(self) -> BridgeStatus:
        """Return whether a qualified runtime is admitted."""

    def owner_epoch(self) -> str:
        """Return the epoch minted when this owner became ready to accept work."""

    def admit(self, request: BridgeSessionRequest) -> BridgeAdmission:
        """Start one native for a session that is already fenced to this epoch."""

    def execute_turn(self, request: BridgeTurnRequest) -> str:
        """Return one strict JSON object, or raise ``BridgeDisconnected``."""

    def close(self, request: BridgeSessionRequest) -> None:
        """Stop the native for this generation. A missing process is success."""


class UnavailableAuthoringBridge:
    """Explicit disabled bridge. Draft save and validate stay available."""

    def status(self) -> BridgeStatus:
        return BridgeStatus(available=False, reason_code=BRIDGE_UNAVAILABLE)

    def owner_epoch(self) -> str:
        raise BridgeDisconnected

    def admit(self, request: BridgeSessionRequest) -> BridgeAdmission:
        del request
        raise BridgeDisconnected

    def execute_turn(self, request: BridgeTurnRequest) -> str:
        del request
        raise BridgeDisconnected

    def close(self, request: BridgeSessionRequest) -> None:
        del request
