"""One bridge owner. Each authoring session has at most one live native."""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass
import json
from pathlib import Path
import shutil
from secrets import token_hex
import sys
import tempfile
import threading

from encode_pipeline.platform.authoring import BRIDGE_UNAVAILABLE, sha256_text
from encode_pipeline.services.authoring_bridge import (
    BridgeAdmission,
    BridgeDisconnected,
    BridgeSessionRequest,
    BridgeStatus,
    BridgeTurnRequest,
)
from encode_pipeline.services.peri_runtime.identity import binary_matches
from encode_pipeline.services.peri_runtime.isolation import (
    isolation_paths_exist,
    native_launcher,
)
from encode_pipeline.services.peri_runtime.runtime import (
    AuthoringRuntime,
    PeriAuthoringRuntime,
    RuntimeSession,
)
from encode_pipeline.services.peri_runtime.session import NativeRejected


ToolDispatch = Callable[..., dict[str, object]]


@dataclass
class _PlatformGrant:
    """HelixWeave authority for one live session. The runtime does not own it."""

    capability_token: str
    user_id: str
    auth_session_digest: str
    forwarding: bool = True
    turn: BridgeTurnRequest | None = None


def _tool_denied() -> dict[str, object]:
    return {
        "contract_version": "hw-a1/1",
        "data": None,
        "issues": [
            {
                "code": "SCOPE_DENIED",
                "retry_class": "fatal",
                "severity": "error",
            }
        ],
        "ok": False,
    }


class NativeAuthoringBridge:
    """Explicitly enabled Peri bridge. A failed admission leaves manual editing usable."""

    def __init__(
        self,
        *,
        binary: Path,
        state_root: Path | None = None,
        repository: object | None = None,
        runtime: AuthoringRuntime | None = None,
    ) -> None:
        self._binary = binary
        self._repository = repository
        self._runtime = runtime if runtime is not None else PeriAuthoringRuntime()
        self._epoch = token_hex(16)
        self._lock = threading.Lock()
        self._sessions: dict[str, RuntimeSession] = {}
        self._grants: dict[str, _PlatformGrant] = {}
        self._handler: ToolDispatch | None = None
        self.closure_notes: dict[str, dict[str, object]] = {}
        self._owned_root = state_root is None
        self._root = state_root or Path(tempfile.mkdtemp(prefix="helixweave-peri-"))
        self._root.mkdir(parents=True, exist_ok=True)
        self._available = binary_matches(binary) and isolation_paths_exist()

    def bind_tool_handler(self, handler: ToolDispatch) -> None:
        self._handler = handler

    def status(self) -> BridgeStatus:
        if self._available:
            return BridgeStatus(available=True)
        return BridgeStatus(available=False, reason_code=BRIDGE_UNAVAILABLE)

    def owner_epoch(self) -> str:
        if not self.status().available:
            raise BridgeDisconnected
        return self._epoch

    def admit(self, request: BridgeSessionRequest) -> BridgeAdmission:
        if not self.status().available or request.owner_epoch != self._epoch:
            raise BridgeDisconnected
        with self._lock:
            current = self._sessions.get(request.sid)
            if (
                current is not None
                and current.service_generation == request.generation
                and current.state == "ready"
            ):
                return BridgeAdmission(instance_id=current.instance, state="ready")
            if current is not None:
                self._close_locked(request.sid)
            directory = (
                self._root / f"{request.sid}-{request.generation}-{token_hex(4)}"
            )
            session = self._runtime.open_session(
                self._binary,
                directory,
                service_generation=request.generation,
                tool_handler=self._dispatch,
            )
            self._sessions[request.sid] = session
            self._grants[request.sid] = _PlatformGrant(
                capability_token=request.capability_token,
                user_id=request.user_id,
                auth_session_digest=request.auth_session_digest,
            )
        try:
            session.initialize()
        except (NativeRejected, OSError):
            self.close(request)
            raise BridgeDisconnected from None
        if session.state != "ready" or session.handler_calls != 0:
            self.close(request)
            raise BridgeDisconnected
        return BridgeAdmission(instance_id=session.instance, state="ready")

    def execute_turn(self, request: BridgeTurnRequest) -> str:
        session = self._live(request)
        if not self._accepts_forwarding(request):
            raise BridgeDisconnected
        answer = _proposal_text(request)
        with self._lock:
            grant = self._grants.get(request.sid)
            if grant is None:
                raise BridgeDisconnected
            grant.turn = request
        session.authorization_check = lambda: self._egress_allowed(request)
        try:
            text = session.run_authoring_turn(answer)
        except (NativeRejected, OSError):
            self._clear_turn(request)
            self.close(
                BridgeSessionRequest(
                    sid=request.sid,
                    generation=request.generation,
                    owner_epoch=request.owner_epoch,
                    capability_token="",
                    schema_digest=request.schema_digest,
                    projection_digest=request.projection_digest,
                    instance_id=request.instance_id,
                )
            )
            raise BridgeDisconnected from None
        finally:
            self._clear_turn(request)
        if not self._accepts_forwarding(request) or not text:
            raise BridgeDisconnected
        return text

    def suppress_forwarding(self, request: BridgeSessionRequest) -> None:
        """Stop model and tool egress for one live generation before it is fenced."""

        with self._lock:
            session = self._sessions.get(request.sid)
            if session is None or session.service_generation != request.generation:
                return
            if (
                request.instance_id is not None
                and session.instance != request.instance_id
            ):
                return
            grant = self._grants.get(request.sid)
            if grant is not None:
                grant.forwarding = False
            session.authorization_check = lambda: False

    def close(self, request: BridgeSessionRequest) -> None:
        with self._lock:
            current = self._sessions.get(request.sid)
            if current is None or current.service_generation != request.generation:
                return
            if request.owner_epoch != self._epoch:
                return
            if (
                request.instance_id is not None
                and current.instance != request.instance_id
            ):
                return
            self._close_locked(request.sid)

    def drop_authorizations(self, digests: Sequence[str]) -> int:
        wanted = {digest for digest in digests if digest}
        if not wanted:
            return 0
        with self._lock:
            matched = [
                sid
                for sid, grant in self._grants.items()
                if grant.auth_session_digest in wanted and sid in self._sessions
            ]
            for sid in matched:
                self._close_locked(sid)
        return len(matched)

    def sweep(self, allowed: Callable[[RuntimeSession], bool] | None = None) -> int:
        checker = self._session_allowed if allowed is None else allowed
        with self._lock:
            stale = []
            for sid, session in self._sessions.items():
                try:
                    keep = checker(session)
                except Exception:
                    keep = False
                if not keep:
                    stale.append(sid)
            for sid in stale:
                self._close_locked(sid)
        return len(stale)

    def diagnostics(self, sid: str) -> dict[str, object]:
        session = self._sessions.get(sid)
        if session is None:
            return {}
        stub = getattr(session, "stub", None)
        requests = [] if stub is None else list(stub.requests)
        return {
            "pid": getattr(session, "pid", None),
            "instance": session.instance,
            "generation": session.service_generation,
            "state": session.state,
            "handler_calls": session.handler_calls,
            "tool_arrivals": session.tool_arrivals,
            "approvals": list(session.approvals),
            "denied_callbacks": list(session.denied_callbacks),
            "argv": getattr(session, "argv", None),
            "model_requests": len(requests),
            "provider_requests": requests,
            "tool_payloads": list(getattr(session, "tool_payloads", ())),
            "native_pid": getattr(session, "native_pid", None),
            "leader_pid": getattr(session, "leader_pid", None),
            "scope_unit": getattr(session, "scope_unit", None),
            "scope_cgroup": getattr(session, "scope_cgroup", None),
            "stub_failures": [] if stub is None else list(stub.failures),
        }

    def _stored_capability_token(self, request: BridgeTurnRequest) -> str:
        with self._lock:
            session = self._sessions.get(request.sid)
            grant = self._grants.get(request.sid)
            if (
                session is None
                or grant is None
                or session.service_generation != request.generation
                or session.instance != request.instance_id
                or session.state != "ready"
            ):
                return ""
            return grant.capability_token

    def _accepts_forwarding(self, request: BridgeTurnRequest) -> bool:
        with self._lock:
            session = self._sessions.get(request.sid)
            grant = self._grants.get(request.sid)
            if (
                session is None
                or grant is None
                or session.service_generation != request.generation
                or session.instance != request.instance_id
                or session.state != "ready"
            ):
                return False
            return grant.forwarding is not False

    def _egress_allowed(self, request: BridgeTurnRequest) -> bool:
        if not self._accepts_forwarding(request):
            return False
        token = self._stored_capability_token(request)
        checker = getattr(self._repository, "tool_read_state", None)
        if checker is None or not token:
            return False
        reason, _record, _turn = checker(
            sid=request.sid,
            turn_id=request.turn_id,
            generation=request.generation,
            instance_id=request.instance_id,
            owner_epoch=request.owner_epoch,
            capability_digest=sha256_text(token),
        )
        return reason == "ok"

    def _session_allowed(self, session: RuntimeSession) -> bool:
        sid = next(
            (key for key, item in self._sessions.items() if item is session),
            None,
        )
        grant = None if sid is None else self._grants.get(sid)
        if grant is None or not grant.user_id or not grant.auth_session_digest:
            return False
        if session.state == "blocked":
            live = getattr(self._repository, "login_is_live", None)
            return bool(
                live is not None and live(grant.user_id, grant.auth_session_digest)
            )
        checker = getattr(self._repository, "connection_still_authorized", None)
        if checker is None or sid is None:
            return False
        return bool(
            checker(
                sid=sid,
                generation=session.service_generation,
                instance_id=session.instance,
                owner_epoch=self._epoch,
                user_id=grant.user_id,
                auth_session_digest=grant.auth_session_digest,
            )
        )

    def _live(self, request: BridgeTurnRequest) -> RuntimeSession:
        with self._lock:
            session = self._sessions.get(request.sid)
        if (
            session is None
            or session.state != "ready"
            or session.instance != request.instance_id
            or session.service_generation != request.generation
            or self._epoch != request.owner_epoch
        ):
            raise BridgeDisconnected
        return session

    def _dispatch(
        self,
        name: str,
        arguments: dict[str, object],
        session: RuntimeSession,
    ) -> dict[str, object]:
        del arguments
        with self._lock:
            grant = next(
                (
                    self._grants[sid]
                    for sid, item in self._sessions.items()
                    if item is session and sid in self._grants
                ),
                None,
            )
            request = None if grant is None else grant.turn
        handler = self._handler
        token = (
            self._stored_capability_token(request)
            if isinstance(request, BridgeTurnRequest)
            else ""
        )
        if handler is None or not isinstance(request, BridgeTurnRequest) or not token:
            return _tool_denied()
        if not self._accepts_forwarding(request):
            return _tool_denied()
        result = handler(
            sid=request.sid,
            turn_id=request.turn_id,
            generation=request.generation,
            instance_id=request.instance_id,
            owner_epoch=request.owner_epoch,
            capability_token=token,
            name=name,
        )
        if not self._accepts_forwarding(request):
            return _tool_denied()
        if (
            name in {"get_authoring_schema", "get_draft_summary"}
            and isinstance(result, dict)
            and result.get("ok") is True
        ):
            print(f"peri-authoring-tool {name}", file=sys.stderr, flush=True)
        return result

    def _close_locked(self, sid: str) -> None:
        session = self._sessions.pop(sid, None)
        if session is None:
            return
        self._grants.pop(sid, None)
        stub = getattr(session, "stub", None)
        note: dict[str, object] = {
            "errors": list(session.errors),
            "handler_calls": session.handler_calls,
            "tool_arrivals": session.tool_arrivals,
            "approvals": list(session.approvals),
            "denied_callbacks": list(session.denied_callbacks),
            "stub_failures": [] if stub is None else list(stub.failures),
        }
        log_path = session.directory / "stderr.log"
        session.close()
        if log_path.exists():
            note["stderr_tail"] = log_path.read_text(
                encoding="utf-8", errors="replace"
            )[-1500:]
        self.closure_notes[sid] = note
        shutil.rmtree(session.directory, ignore_errors=True)

    def _clear_turn(self, request: BridgeTurnRequest) -> None:
        with self._lock:
            grant = self._grants.get(request.sid)
            if grant is not None and grant.turn is request:
                grant.turn = None

    def shutdown(self) -> None:
        native_launcher().abandon_pending()
        with self._lock:
            for sid in list(self._sessions):
                self._close_locked(sid)
        if self._owned_root:
            shutil.rmtree(self._root, ignore_errors=True)


def _proposal_text(request: BridgeTurnRequest) -> str:
    question = request.question
    tokens = question.get("field_tokens")
    preferences = question.get("preferences")
    token = tokens[0] if isinstance(tokens, list) and tokens else None
    value = None
    if isinstance(preferences, list):
        for item in preferences:
            if isinstance(item, dict) and item.get("field_token") == token:
                value = item.get("value")
    if (
        isinstance(token, str)
        and isinstance(value, int)
        and not isinstance(value, bool)
    ):
        changes = [{"field_token": token, "value": value}]
        kind = "proposal"
        reason_codes: list[str] = []
    else:
        changes = []
        kind = "explanation"
        reason_codes = [code for code in request.issue_codes if isinstance(code, str)]
    return json.dumps(
        {
            "base_revision": request.base_revision,
            "changes": changes,
            "kind": kind,
            "reason_codes": reason_codes,
            "schema_digest": request.schema_digest,
        },
        separators=(",", ":"),
        sort_keys=True,
    )
