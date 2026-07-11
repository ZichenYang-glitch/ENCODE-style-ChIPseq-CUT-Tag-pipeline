"""Tests for durable local platform runtime composition."""

from __future__ import annotations

import pytest

from encode_pipeline.persistence import runtime
from encode_pipeline.persistence.runtime import (
    DATABASE_URL_ENV,
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
