"""One bridge owner and the narrow API connection to it.

API workers do not spawn Peri. They send open, turn, and close to the owner
over a private Unix socket. The owner fences foreign epochs before it listens.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
import socket
import socketserver
import threading
from typing import Any

from encode_pipeline.platform.authoring import BRIDGE_UNAVAILABLE
from encode_pipeline.services.authoring_bridge import (
    BridgeAdmission,
    BridgeDisconnected,
    BridgeSessionRequest,
    BridgeStatus,
    BridgeTurnRequest,
)


_MONITOR_SECONDS = 5.0


class RemoteAuthoringBridge:
    """API-side client. A missing owner is unavailable, not a spawned native."""

    def __init__(self, socket_path: Path) -> None:
        self._path = socket_path

    def status(self) -> BridgeStatus:
        try:
            payload = self._call("status")
        except BridgeDisconnected:
            return BridgeStatus(available=False, reason_code=BRIDGE_UNAVAILABLE)
        if payload.get("available") is True and isinstance(
            payload.get("owner_epoch"), str
        ):
            return BridgeStatus(available=True)
        return BridgeStatus(available=False, reason_code=BRIDGE_UNAVAILABLE)

    def owner_epoch(self) -> str:
        payload = self._call("status")
        epoch = payload.get("owner_epoch")
        if payload.get("available") is not True or not isinstance(epoch, str):
            raise BridgeDisconnected
        return epoch

    def admit(self, request: BridgeSessionRequest) -> BridgeAdmission:
        payload = self._call("open", **_session_payload(request))
        instance_id = payload.get("instance_id")
        state = payload.get("state")
        if not isinstance(instance_id, str) or state != "ready":
            raise BridgeDisconnected
        return BridgeAdmission(instance_id=instance_id, state="ready")

    def execute_turn(self, request: BridgeTurnRequest) -> str:
        payload = self._call("turn", **_turn_payload(request))
        answer = payload.get("answer")
        if not isinstance(answer, str) or not answer:
            raise BridgeDisconnected
        return answer

    def close(self, request: BridgeSessionRequest) -> None:
        try:
            self._call("close", **_session_payload(request))
        except BridgeDisconnected:
            return

    def drop_authorizations(self, digests: tuple[str, ...]) -> int:
        payload = self._call("drop_auth", session_digests=list(digests))
        closed = payload.get("closed")
        return closed if isinstance(closed, int) else 0

    def _call(self, method: str, **fields: object) -> dict[str, Any]:
        message = {"method": method, **fields}
        try:
            with socket.socket(socket.AF_UNIX, socket.SOCK_STREAM) as client:
                client.settimeout(60)
                client.connect(str(self._path))
                client.sendall(json.dumps(message, sort_keys=True).encode() + b"\n")
                chunks = []
                while True:
                    piece = client.recv(65536)
                    if not piece:
                        break
                    chunks.append(piece)
                    if b"\n" in piece:
                        break
        except OSError as error:
            raise BridgeDisconnected from error
        raw = b"".join(chunks).split(b"\n", 1)[0]
        try:
            payload = json.loads(raw.decode("utf-8"))
        except (UnicodeError, json.JSONDecodeError) as error:
            raise BridgeDisconnected from error
        if not isinstance(payload, dict) or payload.get("ok") is not True:
            raise BridgeDisconnected
        return payload


class OwnerServer:
    """Accept trusted API calls only after the previous epoch is fenced."""

    def __init__(
        self,
        bridge: Any,
        repository: Any,
        socket_path: Path,
        *,
        before_listen: Any = None,
    ) -> None:
        self._bridge = bridge
        self._repository = repository
        self._path = socket_path
        self._before_listen = before_listen
        self._server: socketserver.UnixStreamServer | None = None
        self._thread: threading.Thread | None = None
        self._stop = threading.Event()
        self._monitor: threading.Thread | None = None
        self.retired = 0

    def start(self) -> None:
        epoch = self._bridge.owner_epoch()
        self.retired = int(self._repository.retire_foreign_epochs(epoch))
        if self._before_listen is not None:
            self._before_listen()
        if self._path.exists():
            self._path.unlink()
        self._path.parent.mkdir(parents=True, exist_ok=True)
        self._server = socketserver.ThreadingUnixStreamServer(
            str(self._path),
            _handler(self),
        )
        os.chmod(self._path, 0o600)
        self._thread = threading.Thread(
            target=self._server.serve_forever,
            daemon=True,
        )
        self._thread.start()
        self._monitor = threading.Thread(target=self._monitor_loop, daemon=True)
        self._monitor.start()

    def sweep(self) -> int:
        sweeper = getattr(self._bridge, "sweep", None)
        if sweeper is None:
            return 0
        return int(sweeper())

    def shutdown(self) -> None:
        self._stop.set()
        if self._server is not None:
            self._server.shutdown()
            self._server.server_close()
        closer = getattr(self._bridge, "shutdown", None)
        if closer is not None:
            closer()
        if self._path.exists():
            self._path.unlink()

    def handle(self, payload: dict[str, Any]) -> dict[str, Any]:
        method = payload.get("method")
        if method == "status":
            status = self._bridge.status()
            epoch = None
            try:
                epoch = self._bridge.owner_epoch()
            except BridgeDisconnected:
                epoch = None
            return {
                "ok": True,
                "available": status.available,
                "owner_epoch": epoch,
                "reason_code": status.reason_code,
            }
        epoch = payload.get("owner_epoch")
        try:
            current = self._bridge.owner_epoch()
        except BridgeDisconnected:
            return {"ok": False, "error": "unavailable"}
        if method in {"open", "turn", "close"} and epoch != current:
            return {"ok": False, "error": "epoch_rejected"}
        try:
            if method == "open":
                admission = self._bridge.admit(_session_request(payload))
                return {
                    "ok": True,
                    "instance_id": admission.instance_id,
                    "state": admission.state,
                }
            if method == "turn":
                answer = self._bridge.execute_turn(_turn_request(payload))
                return {"ok": True, "answer": answer}
            if method == "close":
                self._bridge.close(_session_request(payload))
                return {"ok": True, "closed": True}
            if method == "drop_auth":
                digests = payload.get("session_digests")
                if not isinstance(digests, list):
                    return {"ok": False, "error": "invalid_request"}
                drop = getattr(self._bridge, "drop_authorizations", None)
                closed = 0 if drop is None else int(drop(tuple(digests)))
                return {"ok": True, "closed": closed}
            if method == "sweep":
                return {"ok": True, "closed": self.sweep()}
        except BridgeDisconnected:
            return {"ok": False, "error": "disconnected"}
        return {"ok": False, "error": "unknown_method"}

    def _monitor_loop(self) -> None:
        while not self._stop.wait(_MONITOR_SECONDS):
            try:
                self.sweep()
            except Exception:
                continue


def _handler(owner: OwnerServer) -> type[socketserver.StreamRequestHandler]:
    class Handler(socketserver.StreamRequestHandler):
        def handle(self) -> None:
            raw = self.rfile.readline(65536)
            try:
                payload = json.loads(raw.decode("utf-8"))
            except (UnicodeError, json.JSONDecodeError):
                payload = {}
            if not isinstance(payload, dict):
                payload = {}
            response = owner.handle(payload)
            self.wfile.write(json.dumps(response, sort_keys=True).encode() + b"\n")

    return Handler


def _session_payload(request: BridgeSessionRequest) -> dict[str, object]:
    return {
        "auth_session_digest": request.auth_session_digest,
        "capability_token": request.capability_token,
        "generation": request.generation,
        "instance_id": request.instance_id,
        "owner_epoch": request.owner_epoch,
        "projection_digest": request.projection_digest,
        "schema_digest": request.schema_digest,
        "sid": request.sid,
        "user_id": request.user_id,
    }


def _turn_payload(request: BridgeTurnRequest) -> dict[str, object]:
    return {
        "base_revision": request.base_revision,
        "capability_token": request.capability_token,
        "generation": request.generation,
        "instance_id": request.instance_id,
        "issue_codes": list(request.issue_codes),
        "owner_epoch": request.owner_epoch,
        "projection_digest": request.projection_digest,
        "question": request.question,
        "schema_digest": request.schema_digest,
        "sid": request.sid,
        "turn_id": request.turn_id,
    }


def _session_request(payload: dict[str, Any]) -> BridgeSessionRequest:
    return BridgeSessionRequest(
        sid=str(payload.get("sid", "")),
        generation=int(payload.get("generation", 0)),
        owner_epoch=str(payload.get("owner_epoch", "")),
        capability_token=str(payload.get("capability_token", "")),
        schema_digest=str(payload.get("schema_digest", "")),
        projection_digest=str(payload.get("projection_digest", "")),
        auth_session_digest=str(payload.get("auth_session_digest", "")),
        user_id=str(payload.get("user_id", "")),
        instance_id=(
            None
            if payload.get("instance_id") in {None, ""}
            else str(payload.get("instance_id"))
        ),
    )


def _turn_request(payload: dict[str, Any]) -> BridgeTurnRequest:
    question = payload.get("question")
    codes = payload.get("issue_codes")
    return BridgeTurnRequest(
        sid=str(payload.get("sid", "")),
        turn_id=str(payload.get("turn_id", "")),
        generation=int(payload.get("generation", 0)),
        instance_id=str(payload.get("instance_id", "")),
        owner_epoch=str(payload.get("owner_epoch", "")),
        capability_token=str(payload.get("capability_token", "")),
        question=question if isinstance(question, dict) else {},
        issue_codes=tuple(code for code in codes if isinstance(code, str))
        if isinstance(codes, list)
        else (),
        schema_digest=str(payload.get("schema_digest", "")),
        projection_digest=str(payload.get("projection_digest", "")),
        base_revision=int(payload.get("base_revision", 0)),
    )


def main() -> None:
    """Run the owner from the operator's environment. Not an API worker."""

    import os as _os
    import sys

    from encode_pipeline.persistence import (
        DATABASE_URL_ENV,
        create_session_factory,
        open_run_persistence,
    )
    from encode_pipeline.persistence.authoring import SqlAlchemyAuthoringRepository
    from encode_pipeline.services.authoring import AuthoringService
    from encode_pipeline.services.data_registry import DataRegistryService
    from encode_pipeline.services.defaults import (
        create_default_validation_service,
        create_default_workflow_registry,
    )
    from encode_pipeline.services.peri_runtime.owner import NativeAuthoringBridge
    from encode_pipeline.services.validated_inputs import ValidatedInputService

    binary = _os.environ.get("HELIXWEAVE_PERI_BINARY")
    socket_path = _os.environ.get("HELIXWEAVE_PERI_OWNER_SOCKET")
    if not binary or not socket_path:
        print(
            "HELIXWEAVE_PERI_BINARY and HELIXWEAVE_PERI_OWNER_SOCKET are required",
            file=sys.stderr,
        )
        raise SystemExit(2)
    persistence = open_run_persistence(_os.environ.get(DATABASE_URL_ENV))
    repository = SqlAlchemyAuthoringRepository(
        create_session_factory(persistence.engine)
    )
    registry = create_default_workflow_registry()
    bridge = NativeAuthoringBridge(
        binary=Path(binary),
        repository=repository,
    )
    if not bridge.status().available:
        print("fixed Peri runtime is not admitted", file=sys.stderr)
        raise SystemExit(1)
    service = AuthoringService(
        repository=repository,
        registry=registry,
        data_registry=DataRegistryService(
            repository=persistence.data_registry_repository
        ),
        validated_inputs=ValidatedInputService(
            registry=registry,
            validation_service=create_default_validation_service(registry=registry),
            build_identity_provider=None,
            repository=persistence.repository,
        ),
        snapshots=persistence.repository,
        bridge=bridge,
    )
    bridge.bind_tool_handler(service.read_tool)
    server = OwnerServer(bridge, repository, Path(socket_path))
    server.start()
    try:
        threading.Event().wait()
    finally:
        server.shutdown()
        persistence.close()


if __name__ == "__main__":
    main()
