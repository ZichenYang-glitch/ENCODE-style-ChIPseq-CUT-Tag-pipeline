"""Small no-follow filesystem primitives used by deployment code."""

from __future__ import annotations

import hashlib
import os
from pathlib import Path
import stat

from encode_pipeline.deployment.errors import fail


READ_CHUNK_BYTES = 1024 * 1024


def require_directory(path: Path, *, code: str) -> os.stat_result:
    try:
        observed = path.lstat()
    except OSError:
        raise fail(code, "Deployment storage is unavailable.") from None
    if not stat.S_ISDIR(observed.st_mode) or stat.S_ISLNK(observed.st_mode):
        raise fail(code, "Deployment storage is invalid.")
    return observed


def create_directory(path: Path, *, mode: int = 0o755) -> None:
    """Create one absolute directory tree while refusing existing symlinks."""
    if not path.is_absolute():
        raise fail("DEPLOYMENT_PATH_INVALID", "Deployment path is invalid.")
    current = Path(path.anchor)
    for part in path.parts[1:]:
        current = current / part
        try:
            observed = current.lstat()
        except FileNotFoundError:
            try:
                current.mkdir(mode=mode)
            except OSError:
                raise fail(
                    "DEPLOYMENT_STORAGE_UNAVAILABLE",
                    "Deployment storage is unavailable.",
                ) from None
            continue
        except OSError:
            raise fail(
                "DEPLOYMENT_STORAGE_UNAVAILABLE",
                "Deployment storage is unavailable.",
            ) from None
        if not stat.S_ISDIR(observed.st_mode) or stat.S_ISLNK(observed.st_mode):
            raise fail(
                "DEPLOYMENT_PATH_UNSAFE", "Deployment storage boundary is unsafe."
            )


def read_regular_file(
    path: Path,
    *,
    max_bytes: int,
    code: str,
) -> tuple[bytes, os.stat_result]:
    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0)
    try:
        descriptor = os.open(path, flags)
    except OSError:
        raise fail(code, "Deployment file is unavailable.") from None
    try:
        before = os.fstat(descriptor)
        if (
            not stat.S_ISREG(before.st_mode)
            or before.st_nlink != 1
            or not 0 <= before.st_size <= max_bytes
        ):
            raise fail(code, "Deployment file is invalid.")
        chunks: list[bytes] = []
        remaining = max_bytes + 1
        while remaining > 0:
            chunk = os.read(descriptor, min(READ_CHUNK_BYTES, remaining))
            if not chunk:
                break
            chunks.append(chunk)
            remaining -= len(chunk)
        content = b"".join(chunks)
        after = os.fstat(descriptor)
        if (
            len(content) > max_bytes
            or _file_witness(before) != _file_witness(after)
            or len(content) != before.st_size
        ):
            raise fail(code, "Deployment file changed during verification.")
        return content, after
    finally:
        os.close(descriptor)


def hash_regular_file(
    path: Path,
    *,
    expected_size: int,
    code: str,
) -> tuple[str, os.stat_result]:
    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0)
    try:
        descriptor = os.open(path, flags)
    except OSError:
        raise fail(code, "Deployment file is unavailable.") from None
    try:
        before = os.fstat(descriptor)
        if (
            not stat.S_ISREG(before.st_mode)
            or before.st_nlink != 1
            or before.st_size != expected_size
        ):
            raise fail(code, "Deployment file is invalid.")
        digest = hashlib.sha256()
        observed_size = 0
        while True:
            chunk = os.read(descriptor, READ_CHUNK_BYTES)
            if not chunk:
                break
            observed_size += len(chunk)
            if observed_size > expected_size:
                raise fail(code, "Deployment file is invalid.")
            digest.update(chunk)
        after = os.fstat(descriptor)
        if observed_size != expected_size or _file_witness(before) != _file_witness(
            after
        ):
            raise fail(code, "Deployment file changed during verification.")
        return digest.hexdigest(), after
    finally:
        os.close(descriptor)


def _file_witness(value: os.stat_result) -> tuple[int, ...]:
    return (
        value.st_dev,
        value.st_ino,
        value.st_mode,
        value.st_nlink,
        value.st_uid,
        value.st_gid,
        value.st_size,
        value.st_mtime_ns,
        value.st_ctime_ns,
    )


def fsync_directory(path: Path) -> None:
    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_DIRECTORY", 0)
    try:
        descriptor = os.open(path, flags)
    except OSError:
        raise fail(
            "DEPLOYMENT_STORAGE_UNAVAILABLE", "Deployment storage is unavailable."
        ) from None
    try:
        os.fsync(descriptor)
    except OSError:
        raise fail(
            "DEPLOYMENT_STORAGE_UNAVAILABLE", "Deployment storage is unavailable."
        ) from None
    finally:
        os.close(descriptor)


__all__ = [
    "create_directory",
    "fsync_directory",
    "hash_regular_file",
    "read_regular_file",
    "require_directory",
]
