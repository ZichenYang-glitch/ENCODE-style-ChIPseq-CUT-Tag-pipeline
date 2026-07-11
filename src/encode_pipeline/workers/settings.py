"""Environment-backed settings shared by API queue clients and workers."""

from __future__ import annotations

import math
import os
from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path

from redis.connection import parse_url

from encode_pipeline.persistence.runtime import (
    DATABASE_URL_ENV,
    resolve_database_url,
)


REDIS_URL_ENV = "ENCODE_PIPELINE_REDIS_URL"
REDIS_CONNECT_TIMEOUT_SECONDS_ENV = "ENCODE_PIPELINE_REDIS_CONNECT_TIMEOUT_SECONDS"
REDIS_API_READ_TIMEOUT_SECONDS_ENV = "ENCODE_PIPELINE_REDIS_API_READ_TIMEOUT_SECONDS"
QUEUE_NAME_ENV = "ENCODE_PIPELINE_QUEUE_NAME"
WORKSPACE_ROOT_ENV = "ENCODE_PIPELINE_WORKSPACE_ROOT"
JOB_TIMEOUT_SECONDS_ENV = "ENCODE_PIPELINE_JOB_TIMEOUT_SECONDS"

DEFAULT_REDIS_URL = "redis://localhost:6379/0"
DEFAULT_REDIS_CONNECT_TIMEOUT_SECONDS = 2.0
DEFAULT_REDIS_API_READ_TIMEOUT_SECONDS = 5.0
DEFAULT_QUEUE_NAME = "encode-pipeline"
DEFAULT_JOB_TIMEOUT_SECONDS = 604_800


@dataclass(frozen=True)
class WorkerSettings:
    """Validated process-independent local worker configuration."""

    database_url: str
    redis_url: str = field(repr=False)
    queue_name: str
    workspace_root: Path
    job_timeout_seconds: int = DEFAULT_JOB_TIMEOUT_SECONDS
    redis_connect_timeout_seconds: float = DEFAULT_REDIS_CONNECT_TIMEOUT_SECONDS
    redis_api_read_timeout_seconds: float = DEFAULT_REDIS_API_READ_TIMEOUT_SECONDS

    def __post_init__(self) -> None:
        database_url = resolve_database_url(self.database_url)
        redis_url = _non_empty(self.redis_url, "redis_url")
        try:
            parse_url(redis_url)
        except (TypeError, ValueError):
            # Do not echo a configured URL because it may contain credentials.
            raise ValueError("redis_url must be a valid Redis URL") from None
        queue_name = _non_empty(self.queue_name, "queue_name")
        job_timeout_seconds = _positive_int(
            self.job_timeout_seconds,
            "job_timeout_seconds",
        )
        redis_connect_timeout_seconds = _positive_float(
            self.redis_connect_timeout_seconds,
            "redis_connect_timeout_seconds",
        )
        redis_api_read_timeout_seconds = _positive_float(
            self.redis_api_read_timeout_seconds,
            "redis_api_read_timeout_seconds",
        )
        if not isinstance(self.workspace_root, Path):
            raise ValueError("workspace_root must be a pathlib.Path")
        workspace_root = self.workspace_root.expanduser()
        if not workspace_root.is_absolute():
            raise ValueError("workspace_root must be an absolute path")

        object.__setattr__(self, "database_url", database_url)
        object.__setattr__(self, "redis_url", redis_url)
        object.__setattr__(self, "queue_name", queue_name)
        object.__setattr__(self, "workspace_root", workspace_root)
        object.__setattr__(self, "job_timeout_seconds", job_timeout_seconds)
        object.__setattr__(
            self,
            "redis_connect_timeout_seconds",
            redis_connect_timeout_seconds,
        )
        object.__setattr__(
            self,
            "redis_api_read_timeout_seconds",
            redis_api_read_timeout_seconds,
        )


def load_worker_settings(
    environ: Mapping[str, str] | None = None,
) -> WorkerSettings:
    """Load worker settings from *environ*, applying local-only defaults."""
    source = os.environ if environ is None else environ
    database_url = source.get(DATABASE_URL_ENV)
    if database_url is None:
        database_url = f"sqlite:///{Path.home() / '.encode-pipeline' / 'platform.db'}"

    workspace_value = source.get(WORKSPACE_ROOT_ENV)
    workspace_root = (
        Path.home() / ".encode-pipeline" / "workspaces"
        if workspace_value is None
        else Path(workspace_value)
    )

    return WorkerSettings(
        database_url=database_url,
        redis_url=source.get(REDIS_URL_ENV, DEFAULT_REDIS_URL),
        queue_name=source.get(QUEUE_NAME_ENV, DEFAULT_QUEUE_NAME),
        workspace_root=workspace_root,
        job_timeout_seconds=_positive_int(
            source.get(
                JOB_TIMEOUT_SECONDS_ENV,
                str(DEFAULT_JOB_TIMEOUT_SECONDS),
            ),
            "job_timeout_seconds",
        ),
        redis_connect_timeout_seconds=_positive_float(
            source.get(
                REDIS_CONNECT_TIMEOUT_SECONDS_ENV,
                str(DEFAULT_REDIS_CONNECT_TIMEOUT_SECONDS),
            ),
            "redis_connect_timeout_seconds",
        ),
        redis_api_read_timeout_seconds=_positive_float(
            source.get(
                REDIS_API_READ_TIMEOUT_SECONDS_ENV,
                str(DEFAULT_REDIS_API_READ_TIMEOUT_SECONDS),
            ),
            "redis_api_read_timeout_seconds",
        ),
    )


def _non_empty(value: str, field_name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{field_name} must be a non-empty string")
    return value.strip()


def _positive_int(value: object, field_name: str) -> int:
    if isinstance(value, bool):
        raise ValueError(f"{field_name} must be a positive integer")
    try:
        normalized = int(value)
    except (TypeError, ValueError):
        raise ValueError(f"{field_name} must be a positive integer") from None
    if normalized <= 0:
        raise ValueError(f"{field_name} must be a positive integer")
    return normalized


def _positive_float(value: object, field_name: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{field_name} must be a positive finite number")
    try:
        normalized = float(value)
    except (TypeError, ValueError):
        raise ValueError(f"{field_name} must be a positive finite number") from None
    if normalized <= 0 or not math.isfinite(normalized):
        raise ValueError(f"{field_name} must be a positive finite number")
    return normalized
