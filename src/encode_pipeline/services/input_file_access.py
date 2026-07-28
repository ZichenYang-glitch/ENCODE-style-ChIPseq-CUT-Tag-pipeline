"""Descriptor-owned observation of regular files beneath private pool roots."""

from __future__ import annotations

from contextlib import suppress
from dataclasses import dataclass
import hashlib
import os
from pathlib import Path
import stat

from encode_pipeline.platform.input_registry import (
    InputFileRevision,
    validate_input_file_relative_path,
)
from encode_pipeline.services.private_storage_pools import PrivateStoragePoolConfig


_READ_CHUNK_BYTES = 1024 * 1024
_PUBLIC_ERROR = "input file is missing, changed, or unsafe"

PathFingerprint = tuple[tuple[int, ...], ...]


class InputFileAccessError(RuntimeError):
    """Stable path-redacted local input observation failure."""

    def __init__(self) -> None:
        super().__init__(_PUBLIC_ERROR)


@dataclass(frozen=True)
class FileObservation:
    """Private registration evidence obtained from held descriptors."""

    relative_path: str
    size_bytes: int
    content_sha256: str
    path_fingerprint: PathFingerprint

    def __post_init__(self) -> None:
        validate_input_file_relative_path(self.relative_path)
        if (
            not isinstance(self.size_bytes, int)
            or isinstance(self.size_bytes, bool)
            or self.size_bytes < 0
        ):
            raise ValueError("size_bytes must be a non-negative integer")
        if (
            not isinstance(self.content_sha256, str)
            or len(self.content_sha256) != 64
            or any(
                character not in "0123456789abcdef" for character in self.content_sha256
            )
        ):
            raise ValueError("content_sha256 must be a lowercase SHA-256 digest")
        fingerprint = tuple(tuple(component) for component in self.path_fingerprint)
        if not fingerprint or any(
            not component
            or any(
                not isinstance(value, int) or isinstance(value, bool)
                for value in component
            )
            for component in fingerprint
        ):
            raise ValueError("path_fingerprint must contain descriptor identities")
        object.__setattr__(self, "path_fingerprint", fingerprint)


class InputFileAccess:
    """Open and reverify one configured local root without following symlinks."""

    def __init__(self, root: Path) -> None:
        self._root = root

    @classmethod
    def from_config(
        cls,
        config: PrivateStoragePoolConfig,
        config_key: str,
    ) -> InputFileAccess:
        """Compose access from a private opaque root mapping."""
        if not isinstance(config, PrivateStoragePoolConfig):
            raise InputFileAccessError()
        try:
            return cls(config.root_for(config_key))
        except Exception:
            raise InputFileAccessError() from None

    def observe(self, relative_path: str) -> FileObservation:
        """Hash one stable regular file and confirm its path identity by reopening."""
        descriptors: list[int] = []
        try:
            safe_path = validate_input_file_relative_path(relative_path)
            root_descriptor = _open_absolute_directory(self._root)
            descriptors.append(root_descriptor)
            current = root_descriptor
            fingerprint: list[tuple[int, ...]] = [
                _directory_identity(os.fstat(root_descriptor))
            ]
            components = safe_path.split("/")
            for component in components[:-1]:
                descriptor = os.open(
                    component,
                    _directory_flags(),
                    dir_fd=current,
                )
                descriptors.append(descriptor)
                current = descriptor
                info = os.fstat(descriptor)
                if not stat.S_ISDIR(info.st_mode):
                    raise OSError("unsafe directory")
                fingerprint.append(_directory_identity(info))

            file_descriptor = os.open(
                components[-1],
                _file_flags(),
                dir_fd=current,
            )
            descriptors.append(file_descriptor)
            before = os.fstat(file_descriptor)
            if not stat.S_ISREG(before.st_mode):
                raise OSError("unsafe file")
            before_identity = _file_identity(before)

            digest = hashlib.sha256()
            size_bytes = 0
            while True:
                chunk = os.read(file_descriptor, _READ_CHUNK_BYTES)
                if not chunk:
                    break
                size_bytes += len(chunk)
                digest.update(chunk)

            after = os.fstat(file_descriptor)
            after_identity = _file_identity(after)
            if before_identity != after_identity or size_bytes != after.st_size:
                raise OSError("file changed")
            fingerprint.append(after_identity)
            observed_fingerprint = tuple(fingerprint)
            if self._reopen_fingerprint(safe_path) != observed_fingerprint:
                raise OSError("path changed")
            return FileObservation(
                relative_path=safe_path,
                size_bytes=size_bytes,
                content_sha256=digest.hexdigest(),
                path_fingerprint=observed_fingerprint,
            )
        except (InputFileAccessError, OSError, TypeError, ValueError):
            raise InputFileAccessError() from None
        finally:
            for descriptor in reversed(descriptors):
                with suppress(OSError):
                    os.close(descriptor)

    def reverify(self, expected: InputFileRevision) -> FileObservation:
        """Reopen and compare the exact immutable revision before private use."""
        if not isinstance(expected, InputFileRevision):
            raise InputFileAccessError()
        observation = self.observe(expected.relative_path)
        if (
            observation.size_bytes != expected.size_bytes
            or observation.content_sha256 != expected.content_sha256
        ):
            raise InputFileAccessError()
        return observation

    def _reopen_fingerprint(self, relative_path: str) -> PathFingerprint:
        descriptors: list[int] = []
        try:
            root_descriptor = _open_absolute_directory(self._root)
            descriptors.append(root_descriptor)
            current = root_descriptor
            fingerprint: list[tuple[int, ...]] = [
                _directory_identity(os.fstat(root_descriptor))
            ]
            components = relative_path.split("/")
            for component in components[:-1]:
                descriptor = os.open(
                    component,
                    _directory_flags(),
                    dir_fd=current,
                )
                descriptors.append(descriptor)
                current = descriptor
                info = os.fstat(descriptor)
                if not stat.S_ISDIR(info.st_mode):
                    raise OSError("unsafe directory")
                fingerprint.append(_directory_identity(info))
            descriptor = os.open(
                components[-1],
                _file_flags(),
                dir_fd=current,
            )
            descriptors.append(descriptor)
            info = os.fstat(descriptor)
            if not stat.S_ISREG(info.st_mode):
                raise OSError("unsafe file")
            fingerprint.append(_file_identity(info))
            return tuple(fingerprint)
        finally:
            for descriptor in reversed(descriptors):
                with suppress(OSError):
                    os.close(descriptor)

    def __repr__(self) -> str:
        return f"{type(self).__name__}(private_root=<redacted>)"


def _open_absolute_directory(root: Path) -> int:
    raw_root = os.fspath(root)
    if (
        not isinstance(raw_root, str)
        or not os.path.isabs(raw_root)
        or raw_root.startswith("//")
        or os.path.normpath(raw_root) != raw_root
        or "\\" in raw_root
        or any(ord(character) < 32 or ord(character) == 127 for character in raw_root)
    ):
        raise OSError("unsafe root")
    components = Path(raw_root).parts
    if any(component in {"", ".", ".."} for component in components[1:]):
        raise OSError("unsafe root")
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
            raise OSError("unsafe root")
        return current
    except Exception:
        with suppress(OSError):
            os.close(current)
        raise


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


def _directory_identity(info: os.stat_result) -> tuple[int, int]:
    return (info.st_dev, info.st_ino)


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
