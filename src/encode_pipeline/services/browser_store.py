"""Descriptor-owned storage for browser display derivatives."""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass
import errno
import json
import os
from pathlib import Path
import stat
import tempfile


_DIRECTORY_FLAGS = os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW | os.O_CLOEXEC
_FILE_FLAGS = os.O_RDONLY | os.O_NOFOLLOW | os.O_CLOEXEC
_CONFLICT = frozenset(
    {errno.ENOENT, errno.ENOTDIR, errno.ELOOP, errno.ENXIO, errno.ENODEV}
)


class BrowserStoreError(Exception):
    """A derivative file is missing or no longer matches its summary."""

    def __init__(self, reason: str) -> None:
        super().__init__(reason)
        self.reason = reason


@dataclass
class StoredDerivative:
    """Open derivative bytes plus the public portion of its summary."""

    provenance: dict[str, object]
    size: int
    descriptor: int
    directory_descriptors: tuple[int, ...]
    closed: bool = False

    def close(self) -> None:
        if self.closed:
            return
        self.closed = True
        for descriptor in reversed((self.descriptor, *self.directory_descriptors)):
            try:
                os.close(descriptor)
            except OSError:
                pass


def derivative_directory(
    root: Path, run_id: str, artifact_id: str, revision: str
) -> Path:
    return root / "browser-derivatives" / run_id / artifact_id / revision


def write_derivative(
    root: Path,
    run_id: str,
    artifact_id: str,
    revision: str,
    *,
    data: bytes,
    index: bytes,
    provenance: dict[str, object],
) -> dict[str, object]:
    """Replace one revision's display files and return the stored summary."""

    directory = derivative_directory(root, run_id, artifact_id, revision)
    _require_real_ancestors(directory)
    directory.mkdir(parents=True, exist_ok=True)
    _require_real_ancestors(directory)
    _replace(directory / "track.gz", data)
    _replace(directory / "track.gz.tbi", index)
    summary = dict(provenance)
    summary.update(_file_stamp(directory / "track.gz"))
    summary["index_stamp"] = _file_stamp(directory / "track.gz.tbi")
    _replace(
        directory / "provenance.json",
        json.dumps(summary, sort_keys=True, separators=(",", ":")).encode("ascii"),
    )
    return _public_provenance(summary)


def open_derivative(
    root: Path,
    run_id: str,
    artifact_id: str,
    revision: str,
    *,
    index: bool,
) -> StoredDerivative:
    """Open the data or index file after its stored summary still matches."""

    directory = derivative_directory(root, run_id, artifact_id, revision)
    summary = _read_json(directory / "provenance.json")
    name = "track.gz.tbi" if index else "track.gz"
    stamp = summary["index_stamp"] if index else summary
    if not isinstance(stamp, dict):
        raise BrowserStoreError("derivative_invalid")
    descriptors, file_descriptor = _open_tree(directory / name)
    try:
        info = os.fstat(file_descriptor)
        if not stat.S_ISREG(info.st_mode) or not _stamp_matches(info, stamp):
            raise BrowserStoreError("derivative_changed")
        expected = summary["index_size"] if index else summary["data_size"]
        if info.st_size != expected:
            raise BrowserStoreError("derivative_changed")
    except BrowserStoreError:
        for descriptor in reversed((file_descriptor, *descriptors)):
            os.close(descriptor)
        raise
    return StoredDerivative(
        provenance=_public_provenance(summary),
        size=int(summary["index_size"] if index else summary["data_size"]),
        descriptor=file_descriptor,
        directory_descriptors=tuple(descriptors),
    )


def read_derivative_summary(
    root: Path,
    run_id: str,
    artifact_id: str,
    revision: str,
) -> dict[str, object] | None:
    path = derivative_directory(root, run_id, artifact_id, revision) / "provenance.json"
    try:
        info = path.lstat()
    except OSError:
        return None
    if not stat.S_ISREG(info.st_mode):
        return None
    try:
        return _public_provenance(_read_json(path))
    except (OSError, BrowserStoreError, json.JSONDecodeError):
        return None


def iter_descriptor(descriptor: int, start: int, length: int) -> Iterator[bytes]:
    """Yield one already-open file range."""

    os.lseek(descriptor, start, os.SEEK_SET)
    remaining = length
    while remaining:
        chunk = os.read(descriptor, min(64 * 1024, remaining))
        if not chunk:
            raise BrowserStoreError("derivative_changed")
        remaining -= len(chunk)
        yield chunk


def _require_real_ancestors(path: Path) -> None:
    current = Path("/")
    for part in path.parts[1:]:
        current /= part
        try:
            info = current.lstat()
        except OSError:
            continue
        if stat.S_ISLNK(info.st_mode):
            raise BrowserStoreError("derivative_invalid")


def _public_provenance(summary: dict[str, object]) -> dict[str, object]:
    hidden = {
        "device",
        "inode",
        "mtime_ns",
        "ctime_ns",
        "index_stamp",
    }
    return {key: value for key, value in summary.items() if key not in hidden}


def _file_stamp(path: Path) -> dict[str, int]:
    info = path.stat()
    return {
        "device": info.st_dev,
        "inode": info.st_ino,
        "mtime_ns": info.st_mtime_ns,
        "ctime_ns": info.st_ctime_ns,
    }


def _stamp_matches(info: os.stat_result, stamp: dict[str, object]) -> bool:
    try:
        return (
            info.st_dev == stamp["device"]
            and info.st_ino == stamp["inode"]
            and info.st_mtime_ns == stamp["mtime_ns"]
            and info.st_ctime_ns == stamp["ctime_ns"]
        )
    except (KeyError, TypeError):
        return False


def _replace(path: Path, payload: bytes) -> None:
    descriptor, name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        os.write(descriptor, payload)
        os.close(descriptor)
        descriptor = -1
        os.replace(name, path)
        name = ""
    finally:
        if descriptor >= 0:
            os.close(descriptor)
        if name:
            try:
                os.unlink(name)
            except OSError:
                pass


def _read_json(path: Path) -> dict[str, object]:
    try:
        value = json.loads(path.read_text(encoding="ascii"))
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        raise BrowserStoreError("derivative_invalid") from exc
    if not isinstance(value, dict):
        raise BrowserStoreError("derivative_invalid")
    return value


def _open_tree(path: Path) -> tuple[list[int], int]:
    if not path.is_absolute():
        raise BrowserStoreError("derivative_invalid")
    descriptors: list[int] = []
    try:
        current = os.open("/", _DIRECTORY_FLAGS)
        descriptors.append(current)
        parts = path.parts[1:]
        for index, part in enumerate(parts):
            if part in {"", ".", ".."}:
                raise BrowserStoreError("derivative_invalid")
            final = index == len(parts) - 1
            try:
                opened = os.open(
                    part,
                    _FILE_FLAGS if final else _DIRECTORY_FLAGS,
                    dir_fd=current,
                )
            except OSError as exc:
                if exc.errno in _CONFLICT:
                    raise BrowserStoreError("derivative_missing") from None
                raise BrowserStoreError("derivative_unavailable") from None
            descriptors.append(opened)
            current = opened
            info = os.fstat(opened)
            if final and not stat.S_ISREG(info.st_mode):
                raise BrowserStoreError("derivative_invalid")
            if not final and not stat.S_ISDIR(info.st_mode):
                raise BrowserStoreError("derivative_invalid")
    except BrowserStoreError:
        for descriptor in reversed(descriptors):
            os.close(descriptor)
        raise
    file_descriptor = descriptors.pop()
    return descriptors, file_descriptor
