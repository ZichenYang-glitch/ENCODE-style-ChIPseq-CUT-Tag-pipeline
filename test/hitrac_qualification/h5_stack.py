"""Shared H5 real-platform harness: dedicated Redis, capture SMTP, original stack.

This support module is imported only by the explicit ``real_execution`` H5
qualification tests. It starts a private Redis server, the original FastAPI
application under Uvicorn and the original RQ worker process, all bound to
loopback or a task-owned directory. It never touches a shared Redis, the
production database, another queue, or an external SMTP server. The capture
SMTP sink is an in-process loopback listener used as the notification
transport boundary; no real email leaves the host.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import email
import json
import os
from pathlib import Path
import signal
import socket
import subprocess
import threading
import time

import httpx

from encode_pipeline.persistence import upgrade_database
from encode_pipeline.persistence.database import (
    create_database_engine,
    create_session_factory,
)
from encode_pipeline.persistence.authentication import (
    SqlAlchemyAuthenticationRepository,
)
from encode_pipeline.services.authentication_service import AccountAdministrationService

WORKFLOW_ID = "hitrac-preprocess"
ADMIN_USERNAME = "h5-admin"
ADMIN_PASSWORD = "h5-admin-password-0000"
MEMBER_USERNAME = "h5-member"
MEMBER_PASSWORD = "h5-member-password-000"


class HarnessError(RuntimeError):
    """The H5 harness itself failed; never confuse this with product behaviour."""


def _free_port() -> int:
    with socket.socket() as probe:
        probe.bind(("127.0.0.1", 0))
        return probe.getsockname()[1]


def _wait_http(url: str, process: subprocess.Popen, timeout: float) -> None:
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if process.poll() is not None:
            raise HarnessError(f"process {process.pid} exited {process.returncode}")
        try:
            with httpx.Client(trust_env=False, timeout=1.0) as client:
                if client.get(url).status_code == 200:
                    return
        except httpx.TransportError:
            pass
        time.sleep(0.1)
    raise HarnessError(f"timed out waiting for {url}")


def read_proc_identity(pid: int) -> dict | None:
    """Return one process's durable identity; None only when truly absent."""
    try:
        text = Path(f"/proc/{pid}/stat").read_text()
    except FileNotFoundError:
        return None
    fields = text[text.rfind(")") + 2 :].split()
    return {
        "pid": pid,
        "starttime": int(fields[19]),
        "state": fields[0],
        "group": int(fields[2]),
    }


def _peak_rss(pid: int) -> int | None:
    try:
        for line in Path(f"/proc/{pid}/status").read_text().splitlines():
            if line.startswith("VmHWM:"):
                return int(line.split()[1]) * 1024
    except (FileNotFoundError, IndexError, ValueError):
        return None
    return None


def identity_alive(identity: dict) -> bool:
    current = read_proc_identity(identity["pid"])
    return current is not None and current["starttime"] == identity["starttime"]


def wait_gone(identity: dict, timeout: float) -> bool:
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if not identity_alive(identity):
            return True
        time.sleep(0.02)
    return not identity_alive(identity)


def process_children(pid: int) -> list[int]:
    try:
        text = Path(f"/proc/{pid}/task/{pid}/children").read_text()
    except FileNotFoundError:
        return []
    return [int(word) for word in text.split()]


def process_tree(root_pid: int) -> list[dict]:
    found = {}
    pending = [root_pid]
    while pending:
        current = pending.pop()
        identity = read_proc_identity(current)
        if identity is None or current in found:
            continue
        found[current] = identity
        pending.extend(process_children(current))
    return list(found.values())


class RedisServer:
    """One dedicated throwaway Redis server inside the task root."""

    def __init__(self, root: Path, redis_server: str, ld_library_path: str | None):
        self.root = root
        self.port = _free_port()
        self.url = f"redis://127.0.0.1:{self.port}/0"
        data = root / "redis-data"
        data.mkdir(parents=True)
        self.log_path = root / "redis.log"
        environment = os.environ.copy()
        if ld_library_path:
            environment["LD_LIBRARY_PATH"] = ld_library_path
        self.log_handle = self.log_path.open("ab")
        self.process = subprocess.Popen(
            [
                redis_server,
                "--bind",
                "127.0.0.1",
                "--port",
                str(self.port),
                "--dir",
                str(data),
                "--save",
                "",
                "--appendonly",
                "no",
            ],
            stdout=self.log_handle,
            stderr=subprocess.STDOUT,
            env=environment,
            start_new_session=True,
        )
        deadline = time.monotonic() + 20
        while time.monotonic() < deadline:
            if self.process.poll() is not None:
                raise HarnessError("dedicated Redis exited during startup")
            try:
                import redis

                client = redis.Redis.from_url(
                    self.url, socket_connect_timeout=0.5, socket_timeout=0.5
                )
                try:
                    if client.ping():
                        self.server_version = client.info(section="server")[
                            "redis_version"
                        ]
                        return
                finally:
                    client.close()
            except Exception:
                time.sleep(0.1)
        raise HarnessError("dedicated Redis did not become ready")

    def stop(self) -> dict:
        result = {"pid": self.process.pid, "url": self.url}
        if self.process.poll() is None:
            self.process.terminate()
            try:
                self.process.wait(timeout=10)
            except subprocess.TimeoutExpired:
                self.process.kill()
                self.process.wait(timeout=5)
        result["exit_code"] = self.process.returncode
        self.log_handle.close()
        with socket.socket() as probe:
            result["port_closed"] = probe.connect_ex(("127.0.0.1", self.port)) != 0
        return result


class CaptureSmtp:
    """Minimal loopback SMTP capture sink; the only notification transport used."""

    def __init__(self, root: Path):
        self.root = root
        self.messages: list[dict] = []
        self._lock = threading.Lock()
        self._listener = socket.socket()
        self._listener.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        self._listener.bind(("127.0.0.1", 0))
        self._listener.listen(8)
        self._listener.settimeout(0.2)
        self.port = self._listener.getsockname()[1]
        self._stopping = False
        self._thread = threading.Thread(target=self._serve, daemon=True)
        self._thread.start()

    def _serve(self) -> None:
        while not self._stopping:
            try:
                connection, _ = self._listener.accept()
            except socket.timeout:
                continue
            except OSError:
                return
            threading.Thread(
                target=self._handle, args=(connection,), daemon=True
            ).start()

    @staticmethod
    def _read_line(connection) -> bytes:
        data = b""
        while not data.endswith(b"\r\n"):
            chunk = connection.recv(1)
            if not chunk:
                raise ConnectionError
            data += chunk
            if len(data) > 4096:
                raise ConnectionError
        return data

    def _handle(self, connection) -> None:
        try:
            connection.settimeout(30)
            connection.sendall(b"220 h5-capture.local ESMTP\r\n")
            recipients: list[str] = []
            sender = ""
            while True:
                line = self._read_line(connection)
                verb = line[:4].upper()
                if verb in (b"EHLO", b"HELO"):
                    connection.sendall(b"250-h5-capture.local\r\n250 8BITMIME\r\n")
                elif verb == b"MAIL":
                    sender = line.decode("ascii", "replace").strip()
                    connection.sendall(b"250 OK\r\n")
                elif verb == b"RCPT":
                    recipients.append(line.decode("ascii", "replace").strip())
                    connection.sendall(b"250 OK\r\n")
                elif verb == b"DATA":
                    connection.sendall(b"354 End data with <CR><LF>.<CR><LF>\r\n")
                    body = b""
                    while not body.endswith(b"\r\n.\r\n"):
                        chunk = connection.recv(65536)
                        if not chunk:
                            raise ConnectionError
                        body += chunk
                        if len(body) > 4 * 1024 * 1024:
                            raise ConnectionError
                    with self._lock:
                        self.messages.append(
                            {
                                "sender": sender,
                                "recipients": recipients,
                                "data": body.decode("utf-8", "replace"),
                            }
                        )
                    connection.sendall(b"250 OK\r\n")
                    recipients = []
                    sender = ""
                elif verb == b"RSET":
                    recipients = []
                    sender = ""
                    connection.sendall(b"250 OK\r\n")
                elif verb == b"NOOP":
                    connection.sendall(b"250 OK\r\n")
                elif verb == b"QUIT":
                    connection.sendall(b"221 Bye\r\n")
                    return
                else:
                    connection.sendall(b"502 Command not implemented\r\n")
        except (ConnectionError, OSError, socket.timeout):
            return
        finally:
            try:
                connection.close()
            except OSError:
                pass

    def captured(self) -> list[dict]:
        with self._lock:
            return [dict(message) for message in self.messages]

    def parsed_subjects(self) -> list[str]:
        subjects = []
        for message in self.captured():
            parsed = email.message_from_string(
                message["data"].replace("\r\n.\r\n", "\r\n")
            )
            subjects.append(str(parsed.get("Subject", "")))
        return subjects

    def stop(self) -> dict:
        self._stopping = True
        self._thread.join(timeout=5)
        self._listener.close()
        return {"port": self.port, "messages": len(self.captured())}


@dataclass
class ManagedService:
    name: str
    process: subprocess.Popen
    log_path: Path
    log_handle: object = None
    identity: dict | None = None


@dataclass
class Stack:
    """One isolated original API + worker deployment over a dedicated SQLite."""

    root: Path
    redis_url: str
    queue_name: str
    job_timeout_seconds: int | None = None
    smtp: CaptureSmtp | None = None
    database_path: Path | None = None
    workspace_root_override: Path | None = None
    reference_config_path: Path | None = None
    services: list[ManagedService] = field(default_factory=list)
    api_port: int = 0
    admin_user_id: str = ""

    def __post_init__(self) -> None:
        database = self.database_path or self.root / "platform.db"
        self.database_url = f"sqlite:///{database}"
        self.workspace_root = self.workspace_root_override or (self.root / "workspaces")
        (self.root / "tmp").mkdir(mode=0o700, exist_ok=True)
        self.workspace_root.mkdir(mode=0o700, exist_ok=True)
        (self.root / "logs").mkdir(mode=0o700, exist_ok=True)

    def prepare_database(self, coordinates: dict) -> dict:
        """Original migration, admin bootstrap and reference registration."""
        upgrade_database(self.database_url)
        engine = create_database_engine(self.database_url)
        try:
            repository = SqlAlchemyAuthenticationRepository(
                create_session_factory(engine)
            )
            account = AccountAdministrationService(
                repository=repository
            ).bootstrap_initial_administrator(ADMIN_USERNAME, ADMIN_PASSWORD)
            self.admin_user_id = account.user_id
        finally:
            engine.dispose()

        reference_config = self.root / "reference-profiles.json"
        reference_config.write_text(
            json.dumps(
                {
                    "schema_version": "helixweave-reference-profiles-v1",
                    "profiles": {
                        "h5-tiny": {
                            "bindings": {
                                WORKFLOW_ID: {
                                    "schema_version": "hitrac-reference-profile-v1",
                                    "binding": coordinates["REFERENCE_BINDING"],
                                    "sha256": coordinates["REFERENCE_SHA256"],
                                }
                            }
                        }
                    },
                }
            )
        )
        reference_config.chmod(0o600)
        registration = self._admin_cli(
            [
                "reference-profile",
                "register",
                "--safe-key",
                "h5-tiny",
                "--display-name",
                "H5 tiny reference",
                "--organism",
                "synthetic",
                "--assembly",
                "tiny",
                "--config-key",
                "h5-tiny",
            ],
            reference_config,
        )
        self._admin_cli(
            ["reference-profile", "verify", registration["revision_id"]],
            reference_config,
        )
        self._admin_cli(
            [
                "reference-profile",
                "enable",
                registration["profile_id"],
                "--revision-id",
                registration["revision_id"],
            ],
            reference_config,
        )
        return registration

    def _admin_cli(self, arguments: list[str], reference_config: Path) -> dict:
        argv = [
            _task_python(),
            "-I",
            "-B",
            "-m",
            "encode_pipeline.cli.admin",
            "--database-url",
            self.database_url,
            "--reference-profile-config",
            str(reference_config),
            *arguments,
        ]
        environment = {
            name: value
            for name, value in os.environ.items()
            if not name.startswith(("ENCODE_PIPELINE_", "HELIXWEAVE_"))
        }
        environment["PYTHONDONTWRITEBYTECODE"] = "1"
        environment["TMPDIR"] = str(self.root / "tmp")
        completed = subprocess.run(
            argv,
            capture_output=True,
            text=True,
            timeout=60,
            check=False,
            env=environment,
        )
        (
            self.root
            / "logs"
            / f"admin-{len(self.services)}-{int(time.time() * 1000)}.log"
        ).write_text(
            json.dumps(
                {
                    "argv": argv,
                    "exit_code": completed.returncode,
                    "stdout": completed.stdout,
                    "stderr": completed.stderr,
                },
                indent=2,
            )
        )
        if completed.returncode != 0:
            raise HarnessError(f"admin CLI failed: {completed.stderr.strip()}")
        return json.loads(completed.stdout)

    def environment(self) -> dict:
        environment = {}
        for name, value in os.environ.items():
            if not name.startswith(("ENCODE_PIPELINE_", "HELIXWEAVE_")):
                environment[name] = value
        environment.update(
            {
                "ENCODE_PIPELINE_DATABASE_URL": self.database_url,
                "ENCODE_PIPELINE_WORKSPACE_ROOT": str(self.workspace_root),
                "ENCODE_PIPELINE_REDIS_URL": self.redis_url,
                "ENCODE_PIPELINE_QUEUE_NAME": self.queue_name,
                "ENCODE_PIPELINE_REFERENCE_PROFILE_CONFIG": str(
                    self.reference_config_path
                    or (self.root / "reference-profiles.json")
                ),
                "ENCODE_PIPELINE_REDIS_CONNECT_TIMEOUT_SECONDS": "2",
                "ENCODE_PIPELINE_REDIS_API_READ_TIMEOUT_SECONDS": "5",
                "PYTHONDONTWRITEBYTECODE": "1",
                "TMPDIR": str(self.root / "tmp"),
                "NO_PROXY": "127.0.0.1,localhost,::1",
                "no_proxy": "127.0.0.1,localhost,::1",
            }
        )
        if self.job_timeout_seconds is not None:
            environment["ENCODE_PIPELINE_JOB_TIMEOUT_SECONDS"] = str(
                self.job_timeout_seconds
            )
        if self.smtp is None:
            environment["HELIXWEAVE_TERMINAL_EMAIL_ENABLED"] = "false"
        else:
            environment.update(
                {
                    "HELIXWEAVE_TERMINAL_EMAIL_ENABLED": "true",
                    "HELIXWEAVE_TERMINAL_EMAIL_ADMIN_RECIPIENTS": (
                        "h5-review@example.test"
                    ),
                    "HELIXWEAVE_TERMINAL_EMAIL_FROM": "helixweave-h5@example.test",
                    "HELIXWEAVE_TERMINAL_EMAIL_APPLICATION_BASE_URL": (
                        "http://127.0.0.1"
                    ),
                    "HELIXWEAVE_SMTP_HOST": "127.0.0.1",
                    "HELIXWEAVE_SMTP_PORT": str(self.smtp.port),
                    "HELIXWEAVE_SMTP_TLS_MODE": "local_plaintext",
                }
            )
        return environment

    def start(self, coordinates: dict) -> None:
        environment = self.environment()
        environment["HELIXWEAVE_HITRAC_RUNTIME_BINDING"] = coordinates[
            "RUNTIME_BINDING"
        ]
        environment["HELIXWEAVE_HITRAC_RUNTIME_SHA256"] = coordinates[
            "RUNTIME_BINDING_SHA256"
        ]
        self.api_port = _free_port()
        self._start(
            "api",
            [
                _task_python(),
                "-I",
                "-B",
                "-m",
                "uvicorn",
                "encode_pipeline.api.main:create_app",
                "--factory",
                "--host",
                "127.0.0.1",
                "--port",
                str(self.api_port),
            ],
            environment,
        )
        self._start(
            "worker",
            [_task_python(), "-I", "-B", "-m", "encode_pipeline.workers.cli"],
            environment,
        )
        _wait_http(
            f"http://127.0.0.1:{self.api_port}/api/v1/auth/session",
            self.services[0].process,
            60,
        )
        self._wait_worker(60)

    def start_api_only(self, coordinates: dict, name: str = "api-reopen") -> None:
        """Start one fresh API process on this stack's existing database."""
        environment = self.environment()
        environment["HELIXWEAVE_HITRAC_RUNTIME_BINDING"] = coordinates[
            "RUNTIME_BINDING"
        ]
        environment["HELIXWEAVE_HITRAC_RUNTIME_SHA256"] = coordinates[
            "RUNTIME_BINDING_SHA256"
        ]
        self.api_port = _free_port()
        self._start(
            name,
            [
                _task_python(),
                "-I",
                "-B",
                "-m",
                "uvicorn",
                "encode_pipeline.api.main:create_app",
                "--factory",
                "--host",
                "127.0.0.1",
                "--port",
                str(self.api_port),
            ],
            environment,
        )
        _wait_http(
            f"http://127.0.0.1:{self.api_port}/api/v1/auth/session",
            self.services[-1].process,
            60,
        )

    def worker_process(self, name: str = "worker") -> ManagedService:
        for service in self.services:
            if service.name == name:
                return service
        raise HarnessError(f"service {name} was never started")

    def _start(self, name: str, argv: list[str], environment: dict) -> None:
        log_path = self.root / "logs" / f"{name}.log"
        handle = log_path.open("ab")
        process = subprocess.Popen(
            argv,
            cwd=self.root,
            env=environment,
            stdout=handle,
            stderr=subprocess.STDOUT,
            start_new_session=True,
        )
        self.services.append(
            ManagedService(
                name, process, log_path, handle, read_proc_identity(process.pid)
            )
        )

    def _wait_worker(self, timeout: float) -> None:
        from redis import Redis
        from rq import Worker

        deadline = time.monotonic() + timeout
        while time.monotonic() < deadline:
            if self.services[-1].process.poll() is not None:
                raise HarnessError("worker exited during startup")
            connection = Redis.from_url(
                self.redis_url, socket_connect_timeout=1, socket_timeout=1
            )
            try:
                if any(
                    self.queue_name in worker.queue_names()
                    for worker in Worker.all(connection=connection)
                ):
                    return
            except Exception:
                pass
            finally:
                connection.close()
            time.sleep(0.2)
        raise HarnessError("RQ worker did not register")

    def client(self) -> httpx.Client:
        return httpx.Client(
            base_url=f"http://127.0.0.1:{self.api_port}",
            trust_env=False,
            timeout=30.0,
        )

    def stop(self) -> list[dict]:
        stopped = []
        for service in reversed(self.services):
            record = {
                "name": service.name,
                "pid": service.process.pid,
                "identity": service.identity,
            }
            if service.process.poll() is None:
                record["peak_rss_bytes"] = _peak_rss(service.process.pid)
                os.killpg(service.process.pid, signal.SIGTERM)
                try:
                    service.process.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    os.killpg(service.process.pid, signal.SIGKILL)
                    service.process.wait(timeout=5)
            record["exit_code"] = service.process.returncode
            stopped.append(record)
        for service in self.services:
            if service.log_handle is not None:
                try:
                    service.log_handle.close()
                except Exception:
                    pass
        return stopped


_TASK_PYTHON: str | None = None


def set_task_python(executable: str) -> None:
    global _TASK_PYTHON
    _TASK_PYTHON = executable


def _task_python() -> str:
    if _TASK_PYTHON is None:
        raise HarnessError("task python was not configured")
    return _TASK_PYTHON
