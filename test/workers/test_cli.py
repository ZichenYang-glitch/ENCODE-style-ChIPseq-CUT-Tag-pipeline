"""Tests for the minimal worker command-line entry point."""

from __future__ import annotations

from rq import Worker

from encode_pipeline.persistence.migration_admission import MigrationAdmissionError
from encode_pipeline.persistence.runtime import DatabaseSchemaNotReadyError
from encode_pipeline.workers import cli
from encode_pipeline.workers.timeouts import WorkerUnixSignalDeathPenalty

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
        def __init__(
            self,
            queues,
            *,
            connection,
            serializer,
            work_horse_killed_handler,
        ):
            captured["queues"] = queues
            captured["connection"] = connection
            captured["serializer"] = serializer
            captured["work_horse_killed_handler"] = work_horse_killed_handler

        def work(self, *, burst):
            captured["burst"] = burst

    monkeypatch.setattr(cli, "load_worker_settings", lambda: configured)
    monkeypatch.setattr(
        cli,
        "open_existing_run_persistence",
        lambda _database_url: type("Prepared", (), {"close": lambda self: None})(),
    )
    monkeypatch.setattr(
        cli,
        "create_worker_redis_connection",
        lambda _settings: connection,
    )
    monkeypatch.setattr(
        cli,
        "create_rq_queue",
        lambda _settings, *, connection: queue,
    )
    monkeypatch.setattr(cli, "DurableWorker", FakeWorker)

    assert cli.main(["--burst"]) == 0
    assert captured == {
        "queues": [queue],
        "connection": connection,
        "serializer": cli.JSONSerializer,
        "work_horse_killed_handler": cli.handle_work_horse_killed,
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
        cli,
        "open_existing_run_persistence",
        lambda _database_url: type("Prepared", (), {"close": lambda self: None})(),
    )
    monkeypatch.setattr(
        cli,
        "create_worker_redis_connection",
        lambda _settings: FakeConnection(),
    )
    monkeypatch.setattr(
        cli,
        "create_rq_queue",
        lambda _settings, *, connection: object(),
    )
    monkeypatch.setattr(
        cli,
        "DurableWorker",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(RuntimeError("worker failed")),
    )

    try:
        cli.main(["--burst"])
    except RuntimeError as exc:
        assert str(exc) == "worker failed"
    else:  # pragma: no cover - protects the cleanup assertion
        raise AssertionError("worker construction unexpectedly succeeded")

    assert connection_closed is True


def test_worker_cli_rejects_migration_inventory_before_settings_or_redis(
    monkeypatch,
    capsys,
):
    settings_loaded = False

    def reject_inventory():
        raise MigrationAdmissionError("MIGRATION_REVISION_UNKNOWN")

    def load_settings():
        nonlocal settings_loaded
        settings_loaded = True
        raise AssertionError("settings must not load after migration rejection")

    monkeypatch.setattr(
        cli,
        "verify_migration_execution_inventory",
        reject_inventory,
    )
    monkeypatch.setattr(cli, "load_worker_settings", load_settings)

    assert cli.main(["--burst"]) == 2
    assert settings_loaded is False
    assert capsys.readouterr().err == (
        "migration execution admission failed [MIGRATION_REVISION_UNKNOWN]\n"
    )


def test_worker_cli_rejects_unprepared_schema_before_redis(
    tmp_path,
    monkeypatch,
    capsys,
):
    configured = worker_settings(tmp_path)
    redis_opened = False

    monkeypatch.setattr(cli, "load_worker_settings", lambda: configured)
    monkeypatch.setattr(
        cli,
        "open_existing_run_persistence",
        lambda _database_url: (_ for _ in ()).throw(
            DatabaseSchemaNotReadyError("DATABASE_SCHEMA_NOT_CURRENT")
        ),
    )

    def open_redis(_settings):
        nonlocal redis_opened
        redis_opened = True
        raise AssertionError("Redis must not open for an unprepared database")

    monkeypatch.setattr(cli, "create_worker_redis_connection", open_redis)

    assert cli.main(["--burst"]) == 2
    assert redis_opened is False
    assert "DATABASE_SCHEMA_NOT_CURRENT" in capsys.readouterr().err


def test_worker_cli_rejects_invalid_enabled_notifications_before_redis(
    tmp_path,
    monkeypatch,
    capsys,
):
    configured = worker_settings(tmp_path)
    persistence_opened = False
    redis_opened = False
    monkeypatch.setattr(cli, "load_worker_settings", lambda: configured)
    monkeypatch.setenv("HELIXWEAVE_TERMINAL_EMAIL_ENABLED", "true")

    def open_persistence(_database_url):
        nonlocal persistence_opened
        persistence_opened = True
        raise AssertionError("SQLite must not open for invalid notification config")

    def open_redis(_settings):
        nonlocal redis_opened
        redis_opened = True
        raise AssertionError("Redis must not open for invalid notification config")

    monkeypatch.setattr(cli, "open_existing_run_persistence", open_persistence)
    monkeypatch.setattr(cli, "create_worker_redis_connection", open_redis)

    assert cli.main(["--burst"]) == 2
    assert persistence_opened is False
    assert redis_opened is False
    assert capsys.readouterr().err == (
        "Worker notification configuration is invalid.\n"
    )


def test_durable_worker_uses_hard_timeout_death_penalty():
    assert issubclass(cli.DurableWorker, Worker)
    assert cli.DurableWorker.death_penalty_class is WorkerUnixSignalDeathPenalty
