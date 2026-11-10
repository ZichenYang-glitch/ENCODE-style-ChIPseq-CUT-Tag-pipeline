"""One bridge owner and the narrow API connection to it.

API workers do not spawn Peri. They send open, turn, and close to the owner
over a private Unix socket. The owner takes the endpoint lock before it fences
foreign epochs or listens.
"""

from __future__ import annotations

import fcntl
import json
import os
from pathlib import Path
import select
import socket
import socketserver
import threading
import time
from typing import Any

from encode_pipeline.platform.authoring import BRIDGE_UNAVAILABLE, CONNECTION_LOST
from encode_pipeline.services.authoring_bridge import (
    BridgeAdmission,
    BridgeDisconnected,
    BridgeSessionRequest,
    BridgeStatus,
    BridgeTurnRequest,
)


_MONITOR_SECONDS = 5.0
_CONTROL_TIMEOUT_SECONDS = 180.0
_LONG_CONTROL_METHODS = frozenset({"open", "turn"})


class OwnerEndpointBusy(RuntimeError):
    """Another live owner already holds this endpoint."""


class _AbandonedCall:
    """One caller-loss fence that is in progress, waiting, or already durable."""

    def __init__(self, payload: dict[str, Any]) -> None:
        self.payload = payload
        self.state = "in_progress"


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
        control_timeout: float = _CONTROL_TIMEOUT_SECONDS,
    ) -> None:
        self._bridge = bridge
        self._repository = repository
        self._path = socket_path
        self._before_listen = before_listen
        self._control_timeout = control_timeout
        self._server: socketserver.UnixStreamServer | None = None
        self._thread: threading.Thread | None = None
        self._stop = threading.Event()
        self._monitor: threading.Thread | None = None
        self._lock_file = None
        self._owns_endpoint = False
        self._socket_inode: int | None = None
        self._monitor_interval = _MONITOR_SECONDS
        self._abandoned: dict[tuple[str, str, str], _AbandonedCall] = {}
        self._abandon_lock = threading.Lock()
        self.retired = 0

    def start(self) -> None:
        if self._owns_endpoint:
            raise OwnerEndpointBusy(str(self._path))
        self._acquire_endpoint()
        try:
            epoch = self._bridge.owner_epoch()
            self.retired = int(self._repository.retire_foreign_epochs(epoch))
            if self._before_listen is not None:
                self._before_listen()
            if self._path.exists():
                self._path.unlink()
            self._server = socketserver.ThreadingUnixStreamServer(
                str(self._path),
                _handler(self),
            )
            os.chmod(self._path, 0o600)
            self._socket_inode = self._path.stat().st_ino
            self._thread = threading.Thread(
                target=self._server.serve_forever,
                daemon=True,
            )
            self._thread.start()
            self._monitor = threading.Thread(target=self._monitor_loop, daemon=True)
            self._monitor.start()
        except Exception:
            self._release_endpoint()
            raise

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
            self._server = None
        closer = getattr(self._bridge, "shutdown", None)
        if closer is not None:
            closer()
        self._release_endpoint()

    def _acquire_endpoint(self) -> None:
        self._path.parent.mkdir(parents=True, exist_ok=True)
        lock_path = Path(str(self._path) + ".lock")
        lock_file = lock_path.open("a+")
        try:
            fcntl.flock(lock_file.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError:
            lock_file.close()
            raise OwnerEndpointBusy(str(self._path)) from None
        self._lock_file = lock_file
        self._owns_endpoint = True

    def _release_endpoint(self) -> None:
        if not self._owns_endpoint:
            return
        inode = self._socket_inode
        if inode is not None and self._path.exists():
            try:
                current = self._path.stat().st_ino
            except OSError:
                current = None
            if current == inode:
                self._path.unlink()
        lock_file = self._lock_file
        self._owns_endpoint = False
        self._socket_inode = None
        self._lock_file = None
        if lock_file is not None and not lock_file.closed:
            try:
                fcntl.flock(lock_file.fileno(), fcntl.LOCK_UN)
            except OSError:
                pass
            lock_file.close()

    def abandon_caller(self, payload: dict[str, Any]) -> None:
        """Fence one in-flight control call after its API caller disappears."""

        if payload.get("method") not in _LONG_CONTROL_METHODS:
            return
        sid = str(payload.get("sid") or "")
        generation = str(payload.get("generation") or "")
        key = (str(payload.get("method")), sid, generation)
        work = self._claim_abandon(key, payload)
        if work is None:
            return
        self._run_abandon(key, work)

    def _claim_abandon(
        self,
        key: tuple[str, str, str],
        payload: dict[str, Any],
    ) -> _AbandonedCall | None:
        with self._abandon_lock:
            current = self._abandoned.get(key)
            if current is None:
                current = _AbandonedCall(payload)
                self._abandoned[key] = current
                return current
            if current.state == "pending":
                current.state = "in_progress"
                return current
            return None

    def _run_abandon(
        self,
        key: tuple[str, str, str],
        work: _AbandonedCall,
    ) -> None:
        payload = work.payload
        sid = str(payload.get("sid") or "")
        generation = str(payload.get("generation") or "")
        instance = payload.get("instance_id")
        request = BridgeSessionRequest(
            sid=sid,
            generation=int(generation or "0"),
            owner_epoch=str(payload.get("owner_epoch") or ""),
            capability_token="",
            schema_digest=str(payload.get("schema_digest") or ""),
            projection_digest=str(payload.get("projection_digest") or ""),
            instance_id=None if instance in {None, ""} else str(instance),
        )
        failed = False
        stop = getattr(self._bridge, "suppress_forwarding", None)
        if stop is not None:
            try:
                stop(request)
            except Exception:
                failed = True
        fail = getattr(self._repository, "fail_connection", None)
        if fail is not None and sid and generation:
            try:
                fail(
                    sid,
                    generation=int(generation),
                    owner_epoch=request.owner_epoch,
                    reason_code=CONNECTION_LOST,
                )
            except Exception:
                failed = True
        try:
            self._bridge.close(request)
        except Exception:
            failed = True
        with self._abandon_lock:
            if self._abandoned.get(key) is work:
                work.state = "pending" if failed else "done"

    def finish_request(self, connection: socket.socket, raw: bytes, writer) -> None:
        try:
            payload = json.loads(raw.decode("utf-8"))
        except (UnicodeError, json.JSONDecodeError):
            payload = {}
        if not isinstance(payload, dict):
            payload = {}
        method = payload.get("method")
        lost = threading.Event()
        finished = threading.Event()
        if method in _LONG_CONTROL_METHODS:
            deadline = time.monotonic() + self._control_timeout
            watcher = threading.Thread(
                target=_watch_caller,
                args=(connection, lost, finished, self, payload, deadline),
                daemon=True,
            )
            watcher.start()
        try:
            response = self.handle(payload)
        finally:
            finished.set()
        if lost.is_set():
            self.abandon_caller(payload)
            return
        body = json.dumps(response, sort_keys=True).encode() + b"\n"
        try:
            writer(body)
        except (BrokenPipeError, ConnectionResetError, OSError):
            if method in _LONG_CONTROL_METHODS:
                self.abandon_caller(payload)

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

    def pending_abandon_count(self) -> int:
        with self._abandon_lock:
            return sum(
                1 for work in self._abandoned.values() if work.state == "pending"
            )

    def retry_abandoned(self) -> int:
        """Retry persistent fences that have not succeeded. Caller loss is enough."""

        with self._abandon_lock:
            claimed = []
            for key, work in self._abandoned.items():
                if work.state == "pending":
                    work.state = "in_progress"
                    claimed.append((key, work))
        for key, work in claimed:
            self._run_abandon(key, work)
        return len(claimed)

    def _monitor_loop(self) -> None:
        while not self._stop.wait(self._monitor_interval):
            try:
                self.retry_abandoned()
            except Exception:
                pass
            try:
                self.sweep()
            except Exception:
                continue


def _watch_caller(
    connection,
    lost: threading.Event,
    finished: threading.Event,
    owner,
    payload,
    deadline: float,
) -> None:
    while not finished.is_set():
        if time.monotonic() >= deadline:
            lost.set()
            owner.abandon_caller(payload)
            return
        remaining = min(0.05, max(0.0, deadline - time.monotonic()))
        readable, _writable, _errors = select.select([connection], [], [], remaining)
        if finished.is_set():
            return
        if not readable:
            continue
        try:
            peeked = connection.recv(1, socket.MSG_PEEK | socket.MSG_DONTWAIT)
        except BlockingIOError:
            continue
        except OSError:
            lost.set()
            owner.abandon_caller(payload)
            return
        if peeked == b"":
            lost.set()
            owner.abandon_caller(payload)
            return


def _handler(owner: OwnerServer) -> type[socketserver.StreamRequestHandler]:
    class Handler(socketserver.StreamRequestHandler):
        def handle(self) -> None:
            raw = self.rfile.readline(65536)

            def send(body: bytes) -> None:
                self.wfile.write(body)
                self.wfile.flush()

            owner.finish_request(self.connection, raw, send)

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
