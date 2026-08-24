"""Tests for durable local platform runtime composition."""

from __future__ import annotations

import sqlite3

import pytest

from encode_pipeline.persistence import runtime
from encode_pipeline.persistence.migration_admission import MigrationAdmissionError
from encode_pipeline.persistence.runtime import (
    DATABASE_URL_ENV,
    DatabaseSchemaNotReadyError,
    open_existing_run_persistence,
    open_run_persistence,
    resolve_database_url,
)


def test_resolve_database_url_prefers_explicit_file_url(tmp_path, monkeypatch):
    monkeypatch.setenv(DATABASE_URL_ENV, "sqlite:////ignored/platform.db")
    explicit = f"sqlite:///{tmp_path / 'explicit.db'}"

    assert resolve_database_url(explicit) == explicit


def test_resolve_database_url_uses_environment_override(tmp_path, monkeypatch):
    configured = f"sqlite:///{tmp_path / 'configured.db'}"
    monkeypatch.setenv(DATABASE_URL_ENV, configured)

    assert resolve_database_url() == configured


def test_resolve_database_url_defaults_under_platform_home(tmp_path, monkeypatch):
    monkeypatch.delenv(DATABASE_URL_ENV, raising=False)
    monkeypatch.setattr(runtime.Path, "home", classmethod(lambda _cls: tmp_path))

    assert resolve_database_url() == (
        f"sqlite:///{tmp_path / '.encode-pipeline' / 'platform.db'}"
    )


@pytest.mark.parametrize(
    "database_url",
    [
        "",
        "sqlite://",
        "sqlite:///:memory:",
        "sqlite:///relative/platform.db",
        "sqlite:///~/platform.db",
        "postgresql://db",
    ],
)
def test_resolve_database_url_rejects_non_durable_backends(database_url):
    with pytest.raises(ValueError):
        resolve_database_url(database_url)


def test_open_run_persistence_creates_missing_parent_and_migrates(tmp_path):
    database_path = tmp_path / "nested" / "state" / "platform.db"
    persistence = open_run_persistence(f"sqlite:///{database_path}")

    assert database_path.exists()
    assert persistence.database_url == f"sqlite:///{database_path}"
    assert persistence.repository.list_runs() == ()
    persistence.close()


def test_open_existing_persistence_accepts_exact_admitted_head(tmp_path, monkeypatch):
    database_path = tmp_path / "platform.db"
    database_url = f"sqlite:///{database_path}"
    prepared = open_run_persistence(database_url)
    prepared.close()
    monkeypatch.setattr(
        runtime,
        "upgrade_database",
        lambda _database_url: (_ for _ in ()).throw(
            AssertionError("service open must not run Alembic")
        ),
    )

    persistence = open_existing_run_persistence(database_url)

    assert persistence.repository.list_runs() == ()
    persistence.close()


def test_open_existing_persistence_rejects_missing_database_without_creating_it(
    tmp_path,
):
    database_path = tmp_path / "missing" / "platform.db"

    with pytest.raises(DatabaseSchemaNotReadyError) as raised:
        open_existing_run_persistence(f"sqlite:///{database_path}")

    assert raised.value.reason_code == "DATABASE_SCHEMA_UNAVAILABLE"
    assert not database_path.parent.exists()


def test_open_existing_persistence_rejects_prior_head_without_migrating(
    tmp_path,
    monkeypatch,
):
    database_path = tmp_path / "platform.db"
    database_url = f"sqlite:///{database_path}"
    prepared = open_run_persistence(database_url)
    prepared.close()
    with sqlite3.connect(database_path) as connection:
        connection.execute("UPDATE alembic_version SET version_num = ?", ("old",))
    monkeypatch.setattr(
        runtime,
        "upgrade_database",
        lambda _database_url: (_ for _ in ()).throw(
            AssertionError("service open must not run Alembic")
        ),
    )

    with pytest.raises(DatabaseSchemaNotReadyError) as raised:
        open_existing_run_persistence(database_url)

    assert raised.value.reason_code == "DATABASE_SCHEMA_NOT_CURRENT"
    with sqlite3.connect(database_path) as connection:
        assert connection.execute(
            "SELECT version_num FROM alembic_version"
        ).fetchone() == ("old",)


def test_resolve_database_url_rejects_unparseable_url():
    with pytest.raises(ValueError, match="valid SQLAlchemy URL"):
        resolve_database_url("nonsense-no-scheme")


def test_open_existing_persistence_rejects_engine_head_drift_and_closes(
    tmp_path,
    monkeypatch,
):
    database_path = tmp_path / "platform.db"
    database_url = f"sqlite:///{database_path}"
    prepared = open_run_persistence(database_url)
    prepared.close()
    opened = []
    real_open = runtime._open_persistence

    def recording_open(resolved_url):
        persistence = real_open(resolved_url)
        opened.append(persistence)
        return persistence

    closed = []
    real_close = runtime.RunPersistence.close

    def recording_close(self):
        closed.append(True)
        real_close(self)

    monkeypatch.setattr(runtime, "_open_persistence", recording_open)
    monkeypatch.setattr(runtime.RunPersistence, "close", recording_close)
    monkeypatch.setattr(runtime, "_read_engine_heads", lambda _engine: ("drifted",))

    with pytest.raises(DatabaseSchemaNotReadyError) as raised:
        open_existing_run_persistence(database_url)

    assert raised.value.reason_code == "DATABASE_SCHEMA_NOT_CURRENT"
    assert len(opened) == 1
    assert closed == [True]


def test_verify_existing_schema_rejects_in_memory_database_defense():
    with pytest.raises(DatabaseSchemaNotReadyError) as raised:
        runtime._verify_existing_schema("sqlite:///:memory:")

    assert raised.value.reason_code == "DATABASE_SCHEMA_UNAVAILABLE"


def test_open_existing_persistence_maps_admission_failure(tmp_path, monkeypatch):
    database_path = tmp_path / "platform.db"
    open_run_persistence(f"sqlite:///{database_path}").close()

    def denied():
        raise MigrationAdmissionError("DENIED_FOR_TEST")

    monkeypatch.setattr(runtime, "verify_migration_execution_inventory", denied)

    with pytest.raises(DatabaseSchemaNotReadyError) as raised:
        open_existing_run_persistence(f"sqlite:///{database_path}")

    assert raised.value.reason_code == "DATABASE_SCHEMA_ADMISSION_FAILED"


def test_open_existing_persistence_maps_unreadable_database_file(tmp_path):
    database_path = tmp_path / "platform.db"
    database_path.write_bytes(b"not a sqlite database")

    with pytest.raises(DatabaseSchemaNotReadyError) as raised:
        open_existing_run_persistence(f"sqlite:///{database_path}")

    assert raised.value.reason_code == "DATABASE_SCHEMA_UNAVAILABLE"


class _Rows(list):
    def all(self):
        return list(self)


class _StubEngine:
    def __init__(self, rows=(), error=None):
        self._rows = rows
        self._error = error

    def connect(self):
        return self

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def exec_driver_sql(self, _sql):
        if self._error is not None:
            raise self._error
        return _Rows(self._rows)


def test_read_engine_heads_rejects_non_string_head():
    engine = _StubEngine(rows=[(123,)])

    with pytest.raises(DatabaseSchemaNotReadyError) as raised:
        runtime._read_engine_heads(engine)

    assert raised.value.reason_code == "DATABASE_SCHEMA_UNAVAILABLE"


def test_read_engine_heads_maps_driver_failure():
    engine = _StubEngine(error=RuntimeError("driver exploded"))

    with pytest.raises(DatabaseSchemaNotReadyError) as raised:
        runtime._read_engine_heads(engine)

    assert raised.value.reason_code == "DATABASE_SCHEMA_UNAVAILABLE"


def test_read_engine_heads_reraises_not_ready():
    engine = _StubEngine(error=DatabaseSchemaNotReadyError("DATABASE_SCHEMA_GONE"))

    with pytest.raises(DatabaseSchemaNotReadyError) as raised:
        runtime._read_engine_heads(engine)

    assert raised.value.reason_code == "DATABASE_SCHEMA_GONE"
