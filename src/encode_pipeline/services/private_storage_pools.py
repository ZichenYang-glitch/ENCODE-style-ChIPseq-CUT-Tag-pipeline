"""Load operator-private StoragePool roots without exposing their topology."""

from __future__ import annotations

from contextlib import suppress
import json
import os
from pathlib import Path
import re
import stat
from types import MappingProxyType
from typing import Mapping


PRIVATE_STORAGE_POOL_SCHEMA_VERSION = "storage-pool-roots-v1"

_MAX_CONFIG_BYTES = 64 * 1024
_MAX_STORAGE_POOLS = 256
_MAX_ROOT_BYTES = 4096
_SAFE_CONFIG_KEY = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._:-]{0,254}$")
_PUBLIC_ERROR = "private storage pool configuration is invalid"


class PrivateStoragePoolConfigError(RuntimeError):
    """Stable, redacted failure for private root configuration."""

    def __init__(self) -> None:
        super().__init__(_PUBLIC_ERROR)


class PrivateStoragePoolConfig:
    """Bounded immutable mapping from opaque database keys to private roots."""

    __slots__ = ("_roots",)

    def __init__(self, roots: Mapping[str, Path]) -> None:
        self._roots = MappingProxyType(dict(roots))

    @classmethod
    def from_file(cls, config_path: Path) -> PrivateStoragePoolConfig:
        """Load the closed versioned JSON document with redacted failures."""
        try:
            document = _read_bounded_regular_file(config_path)
            if len(document) > _MAX_CONFIG_BYTES:
                raise ValueError("oversized")
            payload = _strict_json_loads(document)
            roots = _validate_document(payload)
        except Exception:
            raise PrivateStoragePoolConfigError() from None
        return cls(roots)

    @property
    def config_keys(self) -> tuple[str, ...]:
        """Return only the opaque configured keys, never their root values."""
        return tuple(sorted(self._roots))

    def root_for(self, config_key: str) -> Path:
        """Resolve an exact opaque key for private process composition."""
        if (
            not isinstance(config_key, str)
            or _SAFE_CONFIG_KEY.fullmatch(config_key) is None
        ):
            raise PrivateStoragePoolConfigError()
        try:
            return self._roots[config_key]
        except KeyError:
            raise PrivateStoragePoolConfigError() from None

    def __repr__(self) -> str:
        return f"{type(self).__name__}(configured_pool_count={len(self._roots)})"


def load_private_storage_pool_config(
    config_path: Path,
) -> PrivateStoragePoolConfig:
    """Load an operator-private StoragePool root document."""
    return PrivateStoragePoolConfig.from_file(config_path)


def _strict_json_loads(document: bytes) -> object:
    def object_pairs(pairs: list[tuple[str, object]]) -> dict[str, object]:
        result: dict[str, object] = {}
        for key, value in pairs:
            if key in result:
                raise ValueError("duplicate")
            result[key] = value
        return result

    def reject_number(_value: str) -> float:
        raise ValueError("number")

    try:
        return json.loads(
            document.decode("utf-8"),
            object_pairs_hook=object_pairs,
            parse_float=reject_number,
            parse_constant=reject_number,
        )
    except (UnicodeError, ValueError, json.JSONDecodeError, RecursionError):
        raise ValueError("document") from None


def _validate_document(payload: object) -> dict[str, Path]:
    if not isinstance(payload, dict) or set(payload) != {
        "schema_version",
        "storage_pool_roots",
    }:
        raise ValueError("shape")
    if payload["schema_version"] != PRIVATE_STORAGE_POOL_SCHEMA_VERSION:
        raise ValueError("version")
    configured = payload["storage_pool_roots"]
    if not isinstance(configured, dict) or len(configured) > _MAX_STORAGE_POOLS:
        raise ValueError("roots")

    roots: dict[str, Path] = {}
    for config_key, raw_root in configured.items():
        if (
            not isinstance(config_key, str)
            or _SAFE_CONFIG_KEY.fullmatch(config_key) is None
            or not isinstance(raw_root, str)
            or not raw_root
            or len(raw_root.encode("utf-8")) > _MAX_ROOT_BYTES
            or "\\" in raw_root
            or any(
                ord(character) < 32 or ord(character) == 127 for character in raw_root
            )
            or not os.path.isabs(raw_root)
            or raw_root.startswith("//")
            or os.path.normpath(raw_root) != raw_root
        ):
            raise ValueError("root")
        root = Path(raw_root)
        try:
            descriptor = _open_absolute_directory(root)
        except OSError:
            raise ValueError("root") from None
        else:
            os.close(descriptor)
        roots[config_key] = root
    return roots


def _read_bounded_regular_file(path: Path) -> bytes:
    descriptor = _open_absolute_file(path)
    try:
        before = os.fstat(descriptor)
        document = bytearray()
        while len(document) <= _MAX_CONFIG_BYTES:
            chunk = os.read(
                descriptor,
                min(8192, _MAX_CONFIG_BYTES + 1 - len(document)),
            )
            if not chunk:
                break
            document.extend(chunk)
        after = os.fstat(descriptor)
        if _file_identity(before) != _file_identity(after):
            raise OSError("configuration changed")
        return bytes(document)
    finally:
        with suppress(OSError):
            os.close(descriptor)


def _open_absolute_file(path: Path) -> int:
    components = _absolute_components(path)
    if len(components) < 2:
        raise OSError("configuration path is invalid")
    current = os.open("/", _directory_flags())
    leaf: int | None = None
    try:
        for component in components[1:-1]:
            descriptor = os.open(
                component,
                _directory_flags(),
                dir_fd=current,
            )
            os.close(current)
            current = descriptor
        leaf = os.open(
            components[-1],
            _file_flags(),
            dir_fd=current,
        )
        info = os.fstat(leaf)
        if not stat.S_ISREG(info.st_mode):
            raise OSError("configuration is not a regular file")
        result = leaf
        leaf = None
        return result
    finally:
        if leaf is not None:
            with suppress(OSError):
                os.close(leaf)
        with suppress(OSError):
            os.close(current)


def _open_absolute_directory(path: Path) -> int:
    components = _absolute_components(path)
    current = os.open("/", _directory_flags())
    try:
        for component in components[1:]:
            descriptor = os.open(
                component,
                _directory_flags(),
                dir_fd=current,
            )
            os.close(current)
            current = descriptor
        info = os.fstat(current)
        if not stat.S_ISDIR(info.st_mode):
            raise OSError("root is not a directory")
        return current
    except Exception:
        with suppress(OSError):
            os.close(current)
        raise


def _absolute_components(path: Path) -> tuple[str, ...]:
    raw_path = os.fspath(path)
    if (
        not isinstance(raw_path, str)
        or not os.path.isabs(raw_path)
        or raw_path.startswith("//")
        or os.path.normpath(raw_path) != raw_path
        or "\\" in raw_path
        or any(ord(character) < 32 or ord(character) == 127 for character in raw_path)
    ):
        raise OSError("path is invalid")
    components = Path(raw_path).parts
    if any(component in {"", ".", ".."} for component in components[1:]):
        raise OSError("path is invalid")
    return components


def _directory_flags() -> int:
    required = ("O_DIRECTORY", "O_NOFOLLOW", "O_CLOEXEC")
    if any(not hasattr(os, name) for name in required):
        raise OSError("safe directory flags unavailable")
    return os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW | os.O_CLOEXEC


def _file_flags() -> int:
    required = ("O_NOFOLLOW", "O_NONBLOCK", "O_CLOEXEC")
    if any(not hasattr(os, name) for name in required):
        raise OSError("safe file flags unavailable")
    return os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK | os.O_CLOEXEC


def _file_identity(info: os.stat_result) -> tuple[int, ...]:
    return (
        info.st_dev,
        info.st_ino,
        info.st_mode,
        info.st_nlink,
        info.st_size,
        info.st_mtime_ns,
        info.st_ctime_ns,
    )
