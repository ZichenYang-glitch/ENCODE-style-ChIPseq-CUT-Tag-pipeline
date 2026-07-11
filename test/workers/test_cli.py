"""Tests for the minimal worker command-line entry point."""

from __future__ import annotations

from encode_pipeline.workers import cli

from .conftest import worker_settings


def test_worker_cli_starts_named_json_worker_in_burst_mode(tmp_path, monkeypatch):
    configured = worker_settings(tmp_path)
    connection_closed = False

    class FakeConnection:
        def close(self):
            nonlocal connection_closed
            connection_closed = True

    connection = FakeConnection()
    queue = object()
    captured = {}

    class FakeWorker:
        def __init__(self, queues, *, connection, serializer):
            captured["queues"] = queues
            captured["connection"] = connection
            captured["serializer"] = serializer

        def work(self, *, burst):
            captured["burst"] = burst

    monkeypatch.setattr(cli, "load_worker_settings", lambda: configured)
    monkeypatch.setattr(cli, "create_redis_connection", lambda _settings: connection)
    monkeypatch.setattr(
        cli,
        "create_rq_queue",
        lambda _settings, *, connection: queue,
    )
    monkeypatch.setattr(cli, "Worker", FakeWorker)

    assert cli.main(["--burst"]) == 0
    assert captured == {
        "queues": [queue],
        "connection": connection,
        "serializer": cli.JSONSerializer,
        "burst": True,
    }
    assert connection_closed is True


def test_worker_cli_closes_redis_if_worker_construction_fails(tmp_path, monkeypatch):
    configured = worker_settings(tmp_path)
    connection_closed = False

    class FakeConnection:
        def close(self):
            nonlocal connection_closed
            connection_closed = True

    monkeypatch.setattr(cli, "load_worker_settings", lambda: configured)
    monkeypatch.setattr(
        cli, "create_redis_connection", lambda _settings: FakeConnection()
    )
    monkeypatch.setattr(
        cli,
        "create_rq_queue",
        lambda _settings, *, connection: object(),
    )
    monkeypatch.setattr(
        cli,
        "Worker",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(RuntimeError("worker failed")),
    )

    try:
        cli.main(["--burst"])
    except RuntimeError as exc:
        assert str(exc) == "worker failed"
    else:  # pragma: no cover - protects the cleanup assertion
        raise AssertionError("worker construction unexpectedly succeeded")

    assert connection_closed is True
