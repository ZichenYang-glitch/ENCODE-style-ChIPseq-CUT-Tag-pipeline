"""One bridge owner. Each authoring session has at most one live native."""

from __future__ import annotations

from collections.abc import Callable, Sequence
import json
from pathlib import Path
import shutil
from secrets import token_hex
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
from encode_pipeline.services.peri_runtime.isolation import isolation_paths_exist
from encode_pipeline.services.peri_runtime.session import NativeRejected, NativeSession


ToolDispatch = Callable[..., dict[str, object]]


class NativeAuthoringBridge:
    """Explicitly enabled Peri bridge. A failed admission leaves manual editing usable."""

    def __init__(
        self,
        *,
        binary: Path,
        state_root: Path | None = None,
        repository: object | None = None,
    ) -> None:
        self._binary = binary
        self._repository = repository
        self._epoch = token_hex(16)
        self._lock = threading.Lock()
        self._sessions: dict[str, NativeSession] = {}
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
            session = NativeSession(
                self._binary,
                directory,
                service_generation=request.generation,
                tool_handler=self._dispatch,
            )
            session.bound_sid = request.sid
            session.user_id = request.user_id
            session.auth_session_digest = request.auth_session_digest
            self._sessions[request.sid] = session
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
        answer = _proposal_text(request)
        session.turn_context = request
        session.authorization_check = lambda: self._egress_allowed(request)
        try:
            text = session.run_authoring_turn(answer)
        except (NativeRejected, OSError):
            session.turn_context = None
            self.close(
                BridgeSessionRequest(
                    sid=request.sid,
                    generation=request.generation,
                    owner_epoch=request.owner_epoch,
                    capability_token=request.capability_token,
                    schema_digest=request.schema_digest,
                    projection_digest=request.projection_digest,
                )
            )
            raise BridgeDisconnected from None
        finally:
            session.turn_context = None
        if not text:
            raise BridgeDisconnected
        return text

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
                for sid, session in self._sessions.items()
                if session.auth_session_digest in wanted
            ]
            for sid in matched:
                self._close_locked(sid)
        return len(matched)

    def sweep(self, allowed: Callable[[NativeSession], bool] | None = None) -> int:
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
        return {
            "pid": session.pid,
            "instance": session.instance,
            "generation": session.service_generation,
            "state": session.state,
            "handler_calls": session.handler_calls,
            "tool_arrivals": session.tool_arrivals,
            "approvals": list(session.approvals),
            "denied_callbacks": list(session.denied_callbacks),
            "argv": session.argv,
            "model_requests": 0 if session.stub is None else len(session.stub.requests),
            "provider_requests": []
            if session.stub is None
            else list(session.stub.requests),
            "tool_payloads": list(session.tool_payloads),
            "native_pid": session.native_pid,
            "leader_pid": session.leader_pid,
            "scope_unit": session.scope_unit,
            "scope_cgroup": session.scope_cgroup,
            "stub_failures": []
            if session.stub is None
            else list(session.stub.failures),
        }

    def _egress_allowed(self, request: BridgeTurnRequest) -> bool:
        checker = getattr(self._repository, "tool_read_state", None)
        if checker is None:
            return False
        reason, _record, _turn = checker(
            sid=request.sid,
            turn_id=request.turn_id,
            generation=request.generation,
            instance_id=request.instance_id,
            owner_epoch=request.owner_epoch,
            capability_digest=sha256_text(request.capability_token),
        )
        return reason == "ok"

    def _session_allowed(self, session: NativeSession) -> bool:
        if not session.user_id or not session.auth_session_digest:
            return False
        if session.state == "blocked":
            live = getattr(self._repository, "login_is_live", None)
            return bool(
                live is not None and live(session.user_id, session.auth_session_digest)
            )
        checker = getattr(self._repository, "connection_still_authorized", None)
        if checker is None or not session.bound_sid:
            return False
        return bool(
            checker(
                sid=session.bound_sid,
                generation=session.service_generation,
                instance_id=session.instance,
                owner_epoch=self._epoch,
                user_id=session.user_id,
                auth_session_digest=session.auth_session_digest,
            )
        )

    def _live(self, request: BridgeTurnRequest) -> NativeSession:
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
        session: NativeSession,
    ) -> dict[str, object]:
        del arguments
        request = session.turn_context
        handler = self._handler
        if handler is None or not isinstance(request, BridgeTurnRequest):
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
        return handler(
            sid=request.sid,
            turn_id=request.turn_id,
            generation=request.generation,
            instance_id=request.instance_id,
            owner_epoch=request.owner_epoch,
            capability_token=request.capability_token,
            name=name,
        )

    def _close_locked(self, sid: str) -> None:
        session = self._sessions.pop(sid, None)
        if session is None:
            return
        note: dict[str, object] = {
            "errors": list(session.errors),
            "handler_calls": session.handler_calls,
            "tool_arrivals": session.tool_arrivals,
            "approvals": list(session.approvals),
            "denied_callbacks": list(session.denied_callbacks),
            "stub_failures": []
            if session.stub is None
            else list(session.stub.failures),
        }
        log_path = session.directory / "stderr.log"
        session.close()
        if log_path.exists():
            note["stderr_tail"] = log_path.read_text(
                encoding="utf-8", errors="replace"
            )[-1500:]
        self.closure_notes[sid] = note
        shutil.rmtree(session.directory, ignore_errors=True)

    def shutdown(self) -> None:
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
