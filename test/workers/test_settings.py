"""Tests for shared worker configuration."""

from __future__ import annotations

from pathlib import Path

import pytest

from encode_pipeline.persistence.runtime import DATABASE_URL_ENV
from encode_pipeline.workers import settings
from encode_pipeline.workers.settings import (
    DEFAULT_JOB_TIMEOUT_SECONDS,
    DEFAULT_QUEUE_NAME,
    DEFAULT_REDIS_URL,
    JOB_TIMEOUT_SECONDS_ENV,
    QUEUE_NAME_ENV,
    REDIS_URL_ENV,
    WORKSPACE_ROOT_ENV,
    WorkerSettings,
    load_worker_settings,
)


def test_load_worker_settings_reads_shared_environment(tmp_path):
    database_url = f"sqlite:///{tmp_path / 'platform.db'}"
    workspace_root = tmp_path / "workspaces"

    configured = load_worker_settings(
        {
            DATABASE_URL_ENV: database_url,
            REDIS_URL_ENV: "redis://redis.internal:6380/4",
            QUEUE_NAME_ENV: "epigenomics",
            WORKSPACE_ROOT_ENV: str(workspace_root),
            JOB_TIMEOUT_SECONDS_ENV: "3600",
        }
    )

    assert configured == WorkerSettings(
        database_url=database_url,
        redis_url="redis://redis.internal:6380/4",
        queue_name="epigenomics",
        workspace_root=workspace_root,
        job_timeout_seconds=3600,
    )


def test_load_worker_settings_uses_local_defaults(tmp_path, monkeypatch):
    monkeypatch.setattr(settings.Path, "home", classmethod(lambda _cls: tmp_path))

    configured = load_worker_settings({})

    assert configured.database_url == (
        f"sqlite:///{tmp_path / '.encode-pipeline' / 'platform.db'}"
    )
    assert configured.redis_url == DEFAULT_REDIS_URL
    assert configured.queue_name == DEFAULT_QUEUE_NAME
    assert configured.workspace_root == tmp_path / ".encode-pipeline" / "workspaces"
    assert configured.job_timeout_seconds == DEFAULT_JOB_TIMEOUT_SECONDS


def test_worker_settings_requires_absolute_workspace(tmp_path):
    with pytest.raises(ValueError, match="workspace_root must be an absolute path"):
        WorkerSettings(
            database_url=f"sqlite:///{tmp_path / 'platform.db'}",
            redis_url=DEFAULT_REDIS_URL,
            queue_name=DEFAULT_QUEUE_NAME,
            workspace_root=Path("relative/workspaces"),
        )


@pytest.mark.parametrize(
    ("field_name", "value"),
    [("redis_url", "  "), ("queue_name", "")],
)
def test_worker_settings_rejects_blank_transport_values(tmp_path, field_name, value):
    values = {
        "database_url": f"sqlite:///{tmp_path / 'platform.db'}",
        "redis_url": DEFAULT_REDIS_URL,
        "queue_name": DEFAULT_QUEUE_NAME,
        "workspace_root": tmp_path / "workspaces",
    }
    values[field_name] = value

    with pytest.raises(ValueError, match=field_name):
        WorkerSettings(**values)


@pytest.mark.parametrize("redis_url", ["http://redis.example", "not-a-url"])
def test_worker_settings_rejects_non_redis_urls_without_echoing_value(
    tmp_path, redis_url
):
    with pytest.raises(ValueError, match="valid Redis URL") as exc_info:
        WorkerSettings(
            database_url=f"sqlite:///{tmp_path / 'platform.db'}",
            redis_url=redis_url,
            queue_name=DEFAULT_QUEUE_NAME,
            workspace_root=tmp_path / "workspaces",
        )

    assert redis_url not in str(exc_info.value)


@pytest.mark.parametrize("value", ["0", "-1", "not-an-integer", ""])
def test_load_worker_settings_rejects_invalid_job_timeout(value):
    with pytest.raises(ValueError, match="job_timeout_seconds"):
        load_worker_settings({JOB_TIMEOUT_SECONDS_ENV: value})


def test_worker_settings_repr_hides_redis_credentials(tmp_path):
    secret = "redis-password"
    configured = WorkerSettings(
        database_url=f"sqlite:///{tmp_path / 'platform.db'}",
        redis_url=f"redis://worker:{secret}@redis.internal:6379/0",
        queue_name=DEFAULT_QUEUE_NAME,
        workspace_root=tmp_path / "workspaces",
    )

    assert secret not in repr(configured)
    assert configured.redis_url.endswith("@redis.internal:6379/0")
