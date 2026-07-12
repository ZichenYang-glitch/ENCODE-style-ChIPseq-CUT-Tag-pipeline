#!/usr/bin/env python3
"""Run the durable local platform stack as one supervised foreground process."""

from __future__ import annotations

import argparse
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
import json
import os
from pathlib import Path
import signal
import socket
import subprocess
import sys
import time
from urllib.parse import urlparse
from urllib.request import ProxyHandler, build_opener

from redis import Redis
from rq import Worker


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]


@dataclass(frozen=True)
class RuntimeConfig:
    project_root: Path
    frontend_root: Path
    runtime_root: Path
    redis_url: str
    queue_name: str
    api_host: str
    api_port: int
    frontend_host: str
    frontend_port: int
    readiness_timeout: float

    @property
    def database_url(self) -> str:
        return f"sqlite:///{self.runtime_root / 'platform.db'}"

    @property
    def workspace_root(self) -> Path:
        return self.runtime_root / "workspaces"


@dataclass
class ManagedProcess:
    name: str
    process: subprocess.Popen[str]
    log_handle: object


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Start Redis, FastAPI, one RQ worker, and the React frontend."
    )
    parser.add_argument("--project-root", type=Path, default=REPOSITORY_ROOT)
    parser.add_argument(
        "--frontend-root", type=Path, default=REPOSITORY_ROOT / "frontend"
    )
    parser.add_argument(
        "--runtime-root",
        type=Path,
        default=REPOSITORY_ROOT / ".local" / "platform-demo",
    )
    parser.add_argument("--redis-url", default="redis://127.0.0.1:6379/0")
    parser.add_argument("--queue-name", default="encode-pipeline-demo")
    parser.add_argument("--api-host", default="127.0.0.1")
    parser.add_argument("--api-port", type=int, default=8000)
    parser.add_argument("--frontend-host", default="127.0.0.1")
    parser.add_argument("--frontend-port", type=int, default=5173)
    parser.add_argument("--readiness-timeout", type=float, default=60.0)
    return parser


def config_from_args(args: argparse.Namespace) -> RuntimeConfig:
    return RuntimeConfig(
        project_root=args.project_root.expanduser().resolve(),
        frontend_root=args.frontend_root.expanduser().resolve(),
        runtime_root=args.runtime_root.expanduser().resolve(),
        redis_url=args.redis_url,
        queue_name=args.queue_name,
        api_host=args.api_host,
        api_port=args.api_port,
        frontend_host=args.frontend_host,
        frontend_port=args.frontend_port,
        readiness_timeout=args.readiness_timeout,
    )


def build_shared_environment(
    config: RuntimeConfig, environ: Mapping[str, str] | None = None
) -> dict[str, str]:
    environment = dict(os.environ if environ is None else environ)
    existing_pythonpath = environment.get("PYTHONPATH")
    source_root = str(config.project_root / "src")
    no_proxy = environment.get("NO_PROXY", environment.get("no_proxy", ""))
    no_proxy_entries = [entry for entry in no_proxy.split(",") if entry]
    for entry in ("127.0.0.1", "localhost", "::1"):
        if entry not in no_proxy_entries:
            no_proxy_entries.append(entry)
    normalized_no_proxy = ",".join(no_proxy_entries)
    environment.update(
        {
            "ENCODE_PIPELINE_DATABASE_URL": config.database_url,
            "ENCODE_PIPELINE_WORKSPACE_ROOT": str(config.workspace_root),
            "ENCODE_PIPELINE_REDIS_URL": config.redis_url,
            "ENCODE_PIPELINE_QUEUE_NAME": config.queue_name,
            "ENCODE_PIPELINE_REDIS_CONNECT_TIMEOUT_SECONDS": "2",
            "ENCODE_PIPELINE_REDIS_API_READ_TIMEOUT_SECONDS": "5",
            "PYTHONDONTWRITEBYTECODE": "1",
            "PYTHONPATH": (
                source_root
                if not existing_pythonpath
                else os.pathsep.join((source_root, existing_pythonpath))
            ),
            "VITE_API_PROXY_TARGET": f"http://{config.api_host}:{config.api_port}",
            "TMPDIR": str(config.runtime_root / "tmp"),
            "NO_PROXY": normalized_no_proxy,
            "no_proxy": normalized_no_proxy,
        }
    )
    return environment


class PlatformSupervisor:
    def __init__(self, config: RuntimeConfig) -> None:
        self.config = config
        self.environment = build_shared_environment(config)
        self.processes: list[ManagedProcess] = []
        self._stopping = False

    def run(self) -> int:
        self._prepare()
        previous_handlers = {
            signum: signal.signal(signum, self._handle_signal)
            for signum in (signal.SIGINT, signal.SIGTERM)
        }
        try:
            self._ensure_redis()
            self._start_services()
            self._wait_ready()
            print(
                f"Platform ready: http://{self.config.frontend_host}:"
                f"{self.config.frontend_port}"
            )
            print(f"Runtime data: {self.config.runtime_root}")
            print("Press Ctrl-C to stop every service.")
            while not self._stopping:
                for managed in self.processes:
                    return_code = managed.process.poll()
                    if return_code is not None:
                        raise RuntimeError(
                            f"{managed.name} exited unexpectedly with status {return_code}; "
                            f"see {self.config.runtime_root / 'logs' / (managed.name + '.log')}"
                        )
                time.sleep(0.2)
            return 0
        except KeyboardInterrupt:
            return 0
        finally:
            self.stop()
            for signum, handler in previous_handlers.items():
                signal.signal(signum, handler)

    def _prepare(self) -> None:
        if not (self.config.project_root / "src" / "encode_pipeline").is_dir():
            raise ValueError("project root must contain src/encode_pipeline")
        if not (self.config.frontend_root / "package.json").is_file():
            raise ValueError("frontend root must contain package.json")
        if self.config.readiness_timeout <= 0:
            raise ValueError("readiness timeout must be positive")
        self.config.runtime_root.mkdir(parents=True, exist_ok=True)
        self.config.workspace_root.mkdir(parents=True, exist_ok=True)
        (self.config.runtime_root / "tmp").mkdir(parents=True, exist_ok=True)
        (self.config.runtime_root / "logs").mkdir(parents=True, exist_ok=True)
        (self.config.runtime_root / "supervisor.pid").write_text(
            str(os.getpid()), encoding="utf-8"
        )

    def _ensure_redis(self) -> None:
        if _redis_ready(self.config.redis_url):
            return
        parsed = urlparse(self.config.redis_url)
        if parsed.scheme != "redis" or parsed.hostname not in {
            "127.0.0.1",
            "localhost",
        }:
            raise RuntimeError(
                "configured Redis is unavailable and is not a local endpoint"
            )
        if parsed.password or parsed.username or parsed.path not in {"", "/", "/0"}:
            raise RuntimeError(
                "automatic Redis startup supports only an unauthenticated local DB 0"
            )
        port = parsed.port or 6379
        redis_root = self.config.runtime_root / "redis"
        redis_root.mkdir(parents=True, exist_ok=True)
        self._start(
            "redis",
            [
                "redis-server",
                "--bind",
                parsed.hostname or "127.0.0.1",
                "--port",
                str(port),
                "--dir",
                str(redis_root),
                "--save",
                "",
                "--appendonly",
                "no",
            ],
            cwd=self.config.runtime_root,
        )
        _wait_until(
            lambda: _redis_ready(self.config.redis_url),
            self.config.readiness_timeout,
            "Redis",
            self._assert_processes_alive,
        )

    def _start_services(self) -> None:
        self._start(
            "api",
            [
                sys.executable,
                "-m",
                "uvicorn",
                "encode_pipeline.api.main:create_app",
                "--factory",
                "--host",
                self.config.api_host,
                "--port",
                str(self.config.api_port),
            ],
            cwd=self.config.project_root,
        )
        self._start(
            "worker",
            [sys.executable, "-m", "encode_pipeline.workers.cli"],
            cwd=self.config.project_root,
        )
        self._start(
            "frontend",
            [
                "npm",
                "run",
                "dev",
                "--",
                "--host",
                self.config.frontend_host,
                "--port",
                str(self.config.frontend_port),
                "--strictPort",
            ],
            cwd=self.config.frontend_root,
        )

    def _start(self, name: str, argv: Sequence[str], *, cwd: Path) -> None:
        log_path = self.config.runtime_root / "logs" / f"{name}.log"
        log_handle = log_path.open("a", encoding="utf-8")
        try:
            process = subprocess.Popen(
                list(argv),
                cwd=cwd,
                env=self.environment,
                stdout=log_handle,
                stderr=subprocess.STDOUT,
                text=True,
                start_new_session=True,
            )
        except Exception:
            log_handle.close()
            raise
        self.processes.append(ManagedProcess(name, process, log_handle))
        (self.config.runtime_root / "service-pids.json").write_text(
            json.dumps({item.name: item.process.pid for item in self.processes}),
            encoding="utf-8",
        )

    def _wait_ready(self) -> None:
        _wait_until(
            lambda: _http_ready(
                f"http://{self.config.api_host}:{self.config.api_port}/api/v1/workflows/"
            ),
            self.config.readiness_timeout,
            "FastAPI",
            self._assert_processes_alive,
        )
        _wait_until(
            self._worker_ready,
            self.config.readiness_timeout,
            "RQ worker",
            self._assert_processes_alive,
        )
        _wait_until(
            lambda: _http_ready(
                f"http://{self.config.frontend_host}:{self.config.frontend_port}"
            ),
            self.config.readiness_timeout,
            "frontend",
            self._assert_processes_alive,
        )

    def _worker_ready(self) -> bool:
        connection = Redis.from_url(
            self.config.redis_url,
            socket_connect_timeout=1,
            socket_timeout=1,
        )
        try:
            return any(
                self.config.queue_name in worker.queue_names()
                for worker in Worker.all(connection=connection)
            )
        except Exception:
            return False
        finally:
            connection.close()

    def _assert_processes_alive(self) -> None:
        if self._stopping:
            raise KeyboardInterrupt
        for managed in self.processes:
            return_code = managed.process.poll()
            if return_code is not None:
                raise RuntimeError(
                    f"{managed.name} exited during startup with status {return_code}"
                )

    def _handle_signal(self, _signum: int, _frame: object) -> None:
        self._stopping = True

    def stop(self) -> None:
        try:
            if not self.processes:
                return
            self._stopping = True
            for managed in reversed(self.processes):
                _signal_session(managed.process.pid, signal.SIGTERM)
            deadline = time.monotonic() + 3
            while time.monotonic() < deadline and any(
                managed.process.poll() is None for managed in self.processes
            ):
                time.sleep(0.05)
            for managed in reversed(self.processes):
                _signal_session(managed.process.pid, signal.SIGKILL)
                try:
                    managed.process.wait(timeout=2)
                except subprocess.TimeoutExpired:
                    try:
                        managed.process.kill()
                    except ProcessLookupError:
                        pass
                managed.log_handle.close()
            self.processes.clear()
        finally:
            (self.config.runtime_root / "supervisor.pid").unlink(missing_ok=True)


def _redis_ready(redis_url: str) -> bool:
    connection = Redis.from_url(
        redis_url, socket_connect_timeout=0.5, socket_timeout=0.5
    )
    try:
        return bool(connection.ping())
    except Exception:
        return False
    finally:
        connection.close()


def _http_ready(url: str) -> bool:
    try:
        opener = build_opener(ProxyHandler({}))
        with opener.open(url, timeout=1) as response:  # noqa: S310 - fixed local URLs
            return 200 <= response.status < 500
    except Exception:
        return False


def _wait_until(
    predicate,
    timeout_seconds: float,
    description: str,
    health_check=lambda: None,
) -> None:
    deadline = time.monotonic() + timeout_seconds
    while time.monotonic() < deadline:
        health_check()
        if predicate():
            return
        time.sleep(0.1)
    raise RuntimeError(f"Timed out waiting for {description} readiness")


def _session_process_groups(session_id: int) -> tuple[int, ...]:
    groups: set[int] = set()
    try:
        candidates = tuple(Path("/proc").iterdir())
    except OSError:
        return (session_id,)
    for candidate in candidates:
        if not candidate.name.isdigit():
            continue
        try:
            pid = int(candidate.name)
            if os.getsid(pid) != session_id:
                continue
            group = os.getpgid(pid)
        except (OSError, ValueError):
            continue
        if group > 0 and group != os.getpgrp():
            groups.add(group)
    if session_id > 0 and session_id != os.getpgrp():
        groups.add(session_id)
    return tuple(sorted(groups, reverse=True))


def _signal_session(session_id: int, signum: signal.Signals) -> None:
    for group in _session_process_groups(session_id):
        try:
            os.killpg(group, signum)
        except ProcessLookupError:
            pass


def _port_available(host: str, port: int) -> bool:
    with socket.socket() as probe:
        probe.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        try:
            probe.bind((host, port))
        except OSError:
            return False
    return True


def main(argv: Sequence[str] | None = None) -> int:
    config = config_from_args(build_parser().parse_args(argv))
    if not _port_available(config.api_host, config.api_port):
        raise SystemExit(f"API port {config.api_port} is already in use")
    if not _port_available(config.frontend_host, config.frontend_port):
        raise SystemExit(f"frontend port {config.frontend_port} is already in use")
    try:
        return PlatformSupervisor(config).run()
    except (RuntimeError, ValueError, OSError) as error:
        print(f"platform startup failed: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
