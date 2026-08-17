"""Environment-backed settings shared by API queue clients and workers."""

from __future__ import annotations

import math
import os
import re
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
STORAGE_POOL_CONFIG_ENV = "ENCODE_PIPELINE_STORAGE_POOL_CONFIG"
REFERENCE_PROFILE_CONFIG_ENV = "ENCODE_PIPELINE_REFERENCE_PROFILE_CONFIG"
JOB_TIMEOUT_SECONDS_ENV = "ENCODE_PIPELINE_JOB_TIMEOUT_SECONDS"
MANAGED_DOCKER_EXECUTABLE_ENV = "ENCODE_PIPELINE_MANAGED_DOCKER_EXECUTABLE"
MANAGED_DOCKER_SOCKET_ENV = "ENCODE_PIPELINE_MANAGED_DOCKER_SOCKET"
ENCODE_RUNTIME_ROOT_ENV = "HELIXWEAVE_ENCODE_RUNTIME_ROOT"
ENCODE_RUNNER_ROOT_ENV = "HELIXWEAVE_ENCODE_RUNNER_ROOT"
ENCODE_CONDA_PREFIX_ENV = "HELIXWEAVE_ENCODE_CONDA_PREFIX"

DEFAULT_REDIS_URL = "redis://localhost:6379/0"
DEFAULT_REDIS_CONNECT_TIMEOUT_SECONDS = 2.0
DEFAULT_REDIS_API_READ_TIMEOUT_SECONDS = 5.0
DEFAULT_QUEUE_NAME = "encode-pipeline"
DEFAULT_JOB_TIMEOUT_SECONDS = 604_800
DEFAULT_MANAGED_DOCKER_SOCKET = Path("/var/run/docker.sock")
_DEPLOYMENT_IDENTITY = re.compile(r"^sha256-[0-9a-f]{64}$")


@dataclass(frozen=True)
class WorkerSettings:
    """Validated process-independent local worker configuration."""

    database_url: str
    redis_url: str = field(repr=False)
    queue_name: str
    workspace_root: Path
    storage_pool_config: Path | None = field(default=None, repr=False)
    reference_profile_config: Path | None = field(default=None, repr=False)
    encode_runtime_root: Path | None = field(default=None, repr=False)
    encode_runner_root: Path | None = field(default=None, repr=False)
    encode_conda_prefix: Path | None = field(default=None, repr=False)
    job_timeout_seconds: int = DEFAULT_JOB_TIMEOUT_SECONDS
    redis_connect_timeout_seconds: float = DEFAULT_REDIS_CONNECT_TIMEOUT_SECONDS
    redis_api_read_timeout_seconds: float = DEFAULT_REDIS_API_READ_TIMEOUT_SECONDS
    managed_docker_executable: Path | None = field(default=None, repr=False)
    managed_docker_socket: Path = field(
        default=DEFAULT_MANAGED_DOCKER_SOCKET,
        repr=False,
    )

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
        storage_pool_config = self.storage_pool_config
        if storage_pool_config is not None:
            storage_pool_config = _absolute_path(
                storage_pool_config,
                "storage_pool_config",
            )
        reference_profile_config = self.reference_profile_config
        if reference_profile_config is not None:
            reference_profile_config = _absolute_path(
                reference_profile_config,
                "reference_profile_config",
            )
        encode_runtime_root = self.encode_runtime_root
        if encode_runtime_root is not None:
            encode_runtime_root = _absolute_path(
                encode_runtime_root,
                "encode_runtime_root",
            )
            runtime_parents = encode_runtime_root.parents
            if (
                len(runtime_parents) < 5
                or encode_runtime_root.name != "encode-runtime"
                or encode_runtime_root.parent.name != "contracts"
                or encode_runtime_root.parent.parent.name != "payload"
                or _DEPLOYMENT_IDENTITY.fullmatch(
                    encode_runtime_root.parent.parent.parent.name
                )
                is None
                or encode_runtime_root.parent.parent.parent.parent.name != "encode"
                or encode_runtime_root.parent.parent.parent.parent.parent.name
                != "runtimes"
            ):
                raise ValueError(
                    "encode_runtime_root must identify an immutable deployment"
                )
        encode_runner_root = self.encode_runner_root
        if encode_runner_root is not None:
            encode_runner_root = _absolute_path(
                encode_runner_root,
                "encode_runner_root",
            )
        encode_conda_prefix = self.encode_conda_prefix
        if encode_conda_prefix is not None:
            encode_conda_prefix = _absolute_path(
                encode_conda_prefix,
                "encode_conda_prefix",
            )
        if (encode_runner_root is None) != (encode_conda_prefix is None):
            raise ValueError(
                "encode_runner_root and encode_conda_prefix must be configured together"
            )
        if encode_runner_root is not None and (
            encode_runtime_root is None
            or encode_runner_root.name != "runner"
            or encode_conda_prefix is None
            or encode_conda_prefix.name != "conda-envs"
            or encode_runner_root.parent != encode_conda_prefix.parent
            or encode_runner_root.parent.name
            != encode_runtime_root.parent.parent.parent.name
        ):
            raise ValueError(
                "scientific runtime coordinates must match encode_runtime_root"
            )
        managed_docker_executable = self.managed_docker_executable
        if managed_docker_executable is not None:
            managed_docker_executable = _absolute_path(
                managed_docker_executable,
                "managed_docker_executable",
            )
        managed_docker_socket = _absolute_path(
            self.managed_docker_socket,
            "managed_docker_socket",
        )

        object.__setattr__(self, "database_url", database_url)
        object.__setattr__(self, "redis_url", redis_url)
        object.__setattr__(self, "queue_name", queue_name)
        object.__setattr__(self, "workspace_root", workspace_root)
        object.__setattr__(self, "storage_pool_config", storage_pool_config)
        object.__setattr__(
            self,
            "reference_profile_config",
            reference_profile_config,
        )
        object.__setattr__(self, "encode_runtime_root", encode_runtime_root)
        object.__setattr__(self, "encode_runner_root", encode_runner_root)
        object.__setattr__(self, "encode_conda_prefix", encode_conda_prefix)
        object.__setattr__(
            self,
            "managed_docker_executable",
            managed_docker_executable,
        )
        object.__setattr__(self, "managed_docker_socket", managed_docker_socket)
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
        storage_pool_config=(
            None
            if source.get(STORAGE_POOL_CONFIG_ENV) is None
            else Path(source[STORAGE_POOL_CONFIG_ENV])
        ),
        reference_profile_config=(
            None
            if source.get(REFERENCE_PROFILE_CONFIG_ENV) is None
            else Path(source[REFERENCE_PROFILE_CONFIG_ENV])
        ),
        encode_runtime_root=(
            None
            if source.get(ENCODE_RUNTIME_ROOT_ENV) is None
            else Path(source[ENCODE_RUNTIME_ROOT_ENV])
        ),
        encode_runner_root=(
            None
            if source.get(ENCODE_RUNNER_ROOT_ENV) is None
            else Path(source[ENCODE_RUNNER_ROOT_ENV])
        ),
        encode_conda_prefix=(
            None
            if source.get(ENCODE_CONDA_PREFIX_ENV) is None
            else Path(source[ENCODE_CONDA_PREFIX_ENV])
        ),
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
        managed_docker_executable=(
            None
            if source.get(MANAGED_DOCKER_EXECUTABLE_ENV) is None
            else Path(source[MANAGED_DOCKER_EXECUTABLE_ENV])
        ),
        managed_docker_socket=Path(
            source.get(
                MANAGED_DOCKER_SOCKET_ENV,
                str(DEFAULT_MANAGED_DOCKER_SOCKET),
            )
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


def _absolute_path(value: object, field_name: str) -> Path:
    if not isinstance(value, Path):
        raise ValueError(f"{field_name} must be a pathlib.Path")
    expanded = value.expanduser()
    rendered = str(expanded)
    if (
        not expanded.is_absolute()
        or rendered != str(Path(rendered))
        or any(character in rendered for character in ("\x00", "\n", "\r"))
    ):
        raise ValueError(f"{field_name} must be a canonical absolute path")
    return expanded
