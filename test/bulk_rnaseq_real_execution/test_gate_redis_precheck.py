"""Redis format and reachability judged before protected Bulk science.

These cases stay in the deterministic suite. They use a local socket or a
stand-in connection and do not contact a production Redis service.
"""

from __future__ import annotations

from contextlib import contextmanager
import os
from pathlib import Path
import socket
import tempfile
import threading

import pytest

from encode_pipeline.workers.rq_queue import create_api_redis_connection
from encode_pipeline.workers.settings import (
    DEFAULT_REDIS_API_READ_TIMEOUT_SECONDS,
    DEFAULT_REDIS_CONNECT_TIMEOUT_SECONDS,
    REDIS_API_READ_TIMEOUT_SECONDS_ENV,
    REDIS_CONNECT_TIMEOUT_SECONDS_ENV,
    WorkerSettings,
)

from .support import GateSettings, require_gate_settings, verify_gate_redis


_PATH_TOKEN = "private-redis-path-token"
_SECRET = "s3cret-token"
_FORMAT_ERROR = "^redis_url must be a valid Redis URL$"
_UNREACHABLE_ERROR = "^Redis is unreachable$"


class _StandInRedis:
    def __init__(
        self,
        *,
        result: object = True,
        error: BaseException | None = None,
    ) -> None:
        self._result = result
        self._error = error
        self.closed = False

    def ping(self) -> object:
        if self._error is not None:
            raise self._error
        return self._result

    def close(self) -> None:
        self.closed = True


def _canonical_unix_url(path: Path) -> str:
    rendered = path.as_posix()
    if not rendered.startswith("/"):
        raise RuntimeError("Unix socket coordinate must be absolute")
    return "unix://" + rendered


def _single_slash_unix_url(path: Path) -> str:
    return "unix:" + path.as_posix()


def _assert_hidden(text: str, *secrets: str) -> None:
    for secret in secrets:
        assert secret not in text


def _gate_environment(root: Path, redis_url: str) -> dict[str, str]:
    runtime_root = root / "runtime"
    runtime_root.mkdir()
    fixture_manifest = root / "fixture.json"
    fixture_manifest.write_text("{}\n", encoding="utf-8")
    docker_executable = root / "bin" / "docker"
    docker_executable.parent.mkdir()
    docker_executable.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
    docker_executable.chmod(0o755)
    docker_socket = root / "docker.sock"
    docker_socket.write_text("coordinate\n", encoding="utf-8")
    return {
        "HELIXWEAVE_REQUIRE_BULK_RNASEQ_REAL_EXECUTION": "1",
        "HELIXWEAVE_BULK_RNASEQ_RUNTIME_ROOT": str(runtime_root),
        "HELIXWEAVE_BULK_RNASEQ_FIXTURE_MANIFEST": str(fixture_manifest),
        "ENCODE_PIPELINE_TEST_REDIS_URL": redis_url,
        "ENCODE_PIPELINE_MANAGED_DOCKER_EXECUTABLE": str(docker_executable),
        "ENCODE_PIPELINE_MANAGED_DOCKER_SOCKET": str(docker_socket),
    }


def _pop_resp(buffer: bytes) -> tuple[list[bytes] | None, bytes]:
    if not buffer.startswith(b"*"):
        return None, buffer
    header, separator, rest = buffer.partition(b"\r\n")
    if not separator:
        return None, buffer
    try:
        count = int(header[1:])
    except ValueError:
        return None, buffer
    args: list[bytes] = []
    remaining = rest
    for _ in range(count):
        if not remaining.startswith(b"$"):
            return None, buffer
        length_line, separator, remaining = remaining.partition(b"\r\n")
        if not separator:
            return None, buffer
        try:
            length = int(length_line[1:])
        except ValueError:
            return None, buffer
        if len(remaining) < length + 2:
            return None, buffer
        args.append(remaining[:length])
        if remaining[length : length + 2] != b"\r\n":
            return None, buffer
        remaining = remaining[length + 2 :]
    return args, remaining


def _redis_reply(command: list[bytes]) -> bytes:
    name = command[0].upper() if command else b""
    if name == b"PING":
        return b"+PONG\r\n"
    if name == b"HELLO":
        # redis-py 8.1 asks for RESP3 before PING. The map only carries the
        # protocol number; it has no address, credential, or path.
        return b"%1\r\n$5\r\nproto\r\n:3\r\n"
    return b"-ERR unsupported\r\n"


def _serve_ping(path: Path, ready: threading.Event, stop: threading.Event) -> None:
    server = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
    try:
        server.bind(str(path))
        server.listen(1)
        server.settimeout(0.2)
        ready.set()
        while not stop.is_set():
            try:
                connection, _address = server.accept()
            except socket.timeout:
                continue
            with connection:
                connection.settimeout(0.5)
                buffer = b""
                while not stop.is_set():
                    command, rest = _pop_resp(buffer)
                    if command is None:
                        try:
                            chunk = connection.recv(4096)
                        except socket.timeout:
                            continue
                        if not chunk:
                            break
                        buffer += chunk
                        continue
                    buffer = rest
                    connection.sendall(_redis_reply(command))
    finally:
        server.close()


@contextmanager
def _pingable_unix_redis():
    directory = Path(tempfile.mkdtemp(prefix="hwg", dir="/tmp"))
    path = directory / "r.sock"
    ready = threading.Event()
    stop = threading.Event()
    thread = threading.Thread(
        target=_serve_ping,
        args=(path, ready, stop),
        daemon=True,
    )
    thread.start()
    try:
        if not ready.wait(2):
            raise RuntimeError("local Redis stand-in did not become ready")
        yield path
    finally:
        stop.set()
        thread.join(timeout=2)
        path.unlink(missing_ok=True)
        os.rmdir(directory)


def test_single_slash_unix_url_is_invalid_before_any_connection(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def refuse_connection(*_args: object, **_kwargs: object) -> None:
        raise AssertionError("connection attempted")

    monkeypatch.setattr(
        "encode_pipeline.workers.rq_queue.create_api_redis_connection",
        refuse_connection,
    )
    path = Path("/tmp") / _PATH_TOKEN / "redis.sock"
    url = " " + _single_slash_unix_url(path) + " "

    with pytest.raises(AssertionError, match=_FORMAT_ERROR) as caught:
        verify_gate_redis(url, {})

    _assert_hidden(str(caught.value), _PATH_TOKEN, url.strip())
    assert caught.value.__cause__ is None
    assert caught.value.__suppress_context__ is True


def test_gate_rejects_single_slash_redis_before_returning_settings(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def refuse_connection(*_args: object, **_kwargs: object) -> None:
        raise AssertionError("connection attempted")

    monkeypatch.setattr(
        "encode_pipeline.workers.rq_queue.create_api_redis_connection",
        refuse_connection,
    )
    path = Path("/tmp") / _PATH_TOKEN / "redis.sock"
    url = _single_slash_unix_url(path)

    with pytest.raises(AssertionError, match=_FORMAT_ERROR) as caught:
        require_gate_settings(
            _gate_environment(tmp_path, url),
            _socket_probe=lambda _path: True,
        )

    _assert_hidden(str(caught.value), _PATH_TOKEN)
    assert caught.value.__cause__ is None
    assert caught.value.__suppress_context__ is True


def test_canonical_unix_socket_answers_ping_through_the_gate(tmp_path: Path) -> None:
    with _pingable_unix_redis() as path:
        url = _canonical_unix_url(path)
        settings = require_gate_settings(
            _gate_environment(tmp_path, url),
            _socket_probe=lambda _path: True,
        )

    assert isinstance(settings, GateSettings)
    rendered = repr(settings).lower()
    assert "redis" not in rendered
    assert path.name not in rendered
    assert str(path.parent.name) not in rendered


def test_valid_unix_socket_that_is_absent_is_unreachable(tmp_path: Path) -> None:
    directory = Path(tempfile.mkdtemp(prefix="hwg", dir="/tmp"))
    path = directory / _PATH_TOKEN
    url = _canonical_unix_url(path)
    try:
        with pytest.raises(AssertionError, match=_UNREACHABLE_ERROR) as caught:
            require_gate_settings(
                _gate_environment(tmp_path, url),
                _socket_probe=lambda _socket: True,
            )
    finally:
        os.rmdir(directory)

    text = str(caught.value)
    assert "valid Redis URL" not in text
    _assert_hidden(text, _PATH_TOKEN, url)
    assert caught.value.__context__ is None
    assert caught.value.__cause__ is None


@pytest.mark.parametrize(
    "redis_url",
    (
        "redis://127.0.0.1:6379/15",
        "rediss://127.0.0.1:6379/0",
    ),
)
def test_tcp_and_tls_addresses_parse_without_connecting(redis_url: str) -> None:
    observed: list[WorkerSettings] = []

    def connector(settings: WorkerSettings) -> None:
        observed.append(settings)

    settings = verify_gate_redis(redis_url, {}, _connector=connector)

    assert observed == [settings]
    assert settings.redis_url == redis_url
    assert settings.redis_connect_timeout_seconds == (
        DEFAULT_REDIS_CONNECT_TIMEOUT_SECONDS
    )
    assert settings.redis_api_read_timeout_seconds == (
        DEFAULT_REDIS_API_READ_TIMEOUT_SECONDS
    )


def test_production_connection_receives_gate_timeouts_without_ping() -> None:
    settings = verify_gate_redis(
        "redis://127.0.0.1:6379/15",
        {
            REDIS_CONNECT_TIMEOUT_SECONDS_ENV: "1.25",
            REDIS_API_READ_TIMEOUT_SECONDS_ENV: "4.5",
        },
        _connector=lambda _settings: None,
    )

    client = create_api_redis_connection(settings)
    try:
        kwargs = client.connection_pool.connection_kwargs
        assert kwargs["socket_connect_timeout"] == 1.25
        assert kwargs["socket_timeout"] == 4.5
        assert kwargs["retry_on_timeout"] is False
        assert client.connection_pool._created_connections == 0
    finally:
        client.close()


def test_client_exception_does_not_reveal_credentials_or_paths(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    url = f"rediss://operator:{_SECRET}@127.0.0.1:6379/0"

    def leak(_settings: WorkerSettings) -> None:
        raise ConnectionError(f"refused {_SECRET} {_PATH_TOKEN} {url}")

    monkeypatch.setattr(
        "encode_pipeline.workers.rq_queue.create_api_redis_connection",
        leak,
    )

    with pytest.raises(AssertionError, match=_UNREACHABLE_ERROR) as caught:
        verify_gate_redis(url, {})

    _assert_hidden(str(caught.value), _SECRET, _PATH_TOKEN, url, "rediss://")
    assert caught.value.__context__ is None
    assert caught.value.__cause__ is None


def test_negative_ping_is_unreachable_not_a_format_error(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    stand_in = _StandInRedis(result=False)

    def factory(_settings: WorkerSettings) -> _StandInRedis:
        return stand_in

    monkeypatch.setattr(
        "encode_pipeline.workers.rq_queue.create_api_redis_connection",
        factory,
    )

    with pytest.raises(AssertionError, match=_UNREACHABLE_ERROR):
        verify_gate_redis("redis://127.0.0.1:6379/15", {})

    assert stand_in.closed is True


def test_invalid_timeout_fails_before_connection() -> None:
    called = False

    def connector(_settings: WorkerSettings) -> None:
        nonlocal called
        called = True

    with pytest.raises(AssertionError, match="positive finite number") as caught:
        verify_gate_redis(
            "redis://127.0.0.1:6379/15",
            {REDIS_CONNECT_TIMEOUT_SECONDS_ENV: _PATH_TOKEN},
            _connector=connector,
        )

    assert called is False
    _assert_hidden(str(caught.value), _PATH_TOKEN)
