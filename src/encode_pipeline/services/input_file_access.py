"""Descriptor-owned observation of regular files beneath private pool roots."""

from __future__ import annotations

from collections.abc import Iterator
from contextlib import contextmanager, suppress
import ctypes
from dataclasses import dataclass
import errno
import hashlib
import os
from pathlib import Path
import stat
import struct

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

    def verified_regular_path(self, expected: InputFileRevision) -> Path:
        """Return the registered leaf only while that same descriptor identity matches.

        The second observation uses the openat walk. It does not resolve the
        pathname, so a symlink planted after the first hash cannot choose a
        different regular file.
        """

        observation = self.reverify(expected)
        confirmed = self.observe(expected.relative_path)
        if (
            confirmed.relative_path != observation.relative_path
            or confirmed.size_bytes != observation.size_bytes
            or confirmed.content_sha256 != observation.content_sha256
            or confirmed.path_fingerprint != observation.path_fingerprint
        ):
            raise InputFileAccessError()
        descriptor = self._open_relative_leaf(confirmed.relative_path)
        try:
            if not _same_leaf(descriptor, confirmed.path_fingerprint[-1]):
                raise InputFileAccessError()
        finally:
            os.close(descriptor)
        leaf = self._root.joinpath(*confirmed.relative_path.split("/"))
        try:
            relative = leaf.relative_to(self._root)
        except ValueError:
            raise InputFileAccessError() from None
        if relative.as_posix() != confirmed.relative_path:
            raise InputFileAccessError()
        return leaf

    def bind_candidate(self, expected: InputFileRevision, candidate: str) -> int:
        """Open ``candidate`` only when it is still the reverified revision leaf.

        The caller owns the returned descriptor and must hold it until snapshot
        generation finishes. A later execution revalidation is separate.
        """

        if not isinstance(candidate, str) or not candidate.startswith("/"):
            raise InputFileAccessError()
        observation = self.reverify(expected)
        descriptor = _open_absolute_leaf(candidate)
        try:
            if not _same_leaf(descriptor, observation.path_fingerprint[-1]):
                raise InputFileAccessError()
            if os.fstat(descriptor).st_size != observation.size_bytes:
                raise InputFileAccessError()
            return os.dup(descriptor)
        except InputFileAccessError:
            raise
        except (OSError, TypeError, ValueError):
            raise InputFileAccessError() from None
        finally:
            os.close(descriptor)

    def candidate_still_matches(
        self,
        expected: InputFileRevision,
        candidate: str,
        held: os.stat_result,
    ) -> None:
        """Raise when the held leaf and the candidate path have diverged."""

        self.reverify(expected)
        descriptor = _open_absolute_leaf(candidate)
        try:
            current = os.fstat(descriptor)
        finally:
            os.close(descriptor)
        if (
            not stat.S_ISREG(held.st_mode)
            or not stat.S_ISREG(current.st_mode)
            or (held.st_dev, held.st_ino) != (current.st_dev, current.st_ino)
            or current.st_size != held.st_size
        ):
            raise InputFileAccessError()

    def _open_relative_leaf(self, relative_path: str) -> int:
        descriptors: list[int] = []
        try:
            safe_path = validate_input_file_relative_path(relative_path)
            root_descriptor = _open_absolute_directory(self._root)
            descriptors.append(root_descriptor)
            current = root_descriptor
            components = safe_path.split("/")
            for component in components[:-1]:
                descriptor = os.open(component, _directory_flags(), dir_fd=current)
                descriptors.append(descriptor)
                current = descriptor
                if not stat.S_ISDIR(os.fstat(descriptor).st_mode):
                    raise OSError("unsafe directory")
            leaf = os.open(components[-1], _file_flags(), dir_fd=current)
        except (OSError, TypeError, ValueError):
            raise InputFileAccessError() from None
        else:
            return leaf
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


def _same_leaf(descriptor: int, identity: tuple[int, ...]) -> bool:
    info = os.fstat(descriptor)
    return (
        stat.S_ISREG(info.st_mode)
        and info.st_dev == identity[0]
        and info.st_ino == identity[1]
    )


def _open_absolute_leaf(path: str) -> int:
    candidate = Path(path)
    if (
        not candidate.is_absolute()
        or not candidate.name
        or candidate.name in {".", ".."}
    ):
        raise InputFileAccessError()
    try:
        directory = os.open(candidate.parent, _directory_flags())
    except (OSError, TypeError, ValueError):
        raise InputFileAccessError() from None
    try:
        return os.open(candidate.name, _file_flags(), dir_fd=directory)
    except (OSError, TypeError, ValueError):
        raise InputFileAccessError() from None
    finally:
        with suppress(OSError):
            os.close(directory)


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


# Linux uapi inotify bits. IN_NONBLOCK and IN_CLOEXEC are the fcntl open flags.
_IN_MODIFY = 0x00000002
_IN_ATTRIB = 0x00000004
_IN_CLOSE_WRITE = 0x00000008
_IN_MOVED_FROM = 0x00000040
_IN_MOVED_TO = 0x00000080
_IN_CREATE = 0x00000100
_IN_DELETE = 0x00000200
_IN_Q_OVERFLOW = 0x00004000
_IN_IGNORED = 0x00008000
_IN_ONLYDIR = 0x01000000
_IN_NONBLOCK = getattr(os, "O_NONBLOCK", 0)
_IN_CLOEXEC = getattr(os, "O_CLOEXEC", 0)
_LEAF_CHANGE_MASK = (
    _IN_CREATE
    | _IN_DELETE
    | _IN_MODIFY
    | _IN_ATTRIB
    | _IN_CLOSE_WRITE
    | _IN_MOVED_FROM
    | _IN_MOVED_TO
)
_LEAF_WATCH_MASK = _LEAF_CHANGE_MASK | _IN_ONLYDIR
_INOTIFY_EVENT = struct.Struct("=iIII")
_INOTIFY_READ_BYTES = 4096
_INOTIFY_FUNCTIONS: tuple[object, object, object] | None = None


def _inotify_functions() -> tuple[object, object, object]:
    global _INOTIFY_FUNCTIONS
    if _INOTIFY_FUNCTIONS is None:
        libc = ctypes.CDLL(None, use_errno=True)
        init1 = ctypes.CFUNCTYPE(ctypes.c_int, ctypes.c_int, use_errno=True)(
            ("inotify_init1", libc)
        )
        add_watch = ctypes.CFUNCTYPE(
            ctypes.c_int,
            ctypes.c_int,
            ctypes.c_char_p,
            ctypes.c_uint32,
            use_errno=True,
        )(("inotify_add_watch", libc))
        read = ctypes.CFUNCTYPE(
            ctypes.c_ssize_t,
            ctypes.c_int,
            ctypes.c_void_p,
            ctypes.c_size_t,
            use_errno=True,
        )(("read", libc))
        _INOTIFY_FUNCTIONS = (init1, add_watch, read)
    return _INOTIFY_FUNCTIONS


def _inotify_fd_flags() -> int:
    if not _IN_NONBLOCK or not _IN_CLOEXEC:
        raise OSError("inotify flags unavailable")
    return _IN_NONBLOCK | _IN_CLOEXEC


def _single_leaf_name(leaf_name: str) -> bytes:
    if not isinstance(leaf_name, str) or "/" in leaf_name or "\x00" in leaf_name:
        raise ValueError("leaf_name must be one path component")
    if not leaf_name or leaf_name in {".", ".."}:
        raise ValueError("leaf_name must be one path component")
    encoded = os.fsencode(leaf_name)
    if not encoded or b"/" in encoded or b"\x00" in encoded:
        raise ValueError("leaf_name must be one path component")
    return encoded


class _LeafGenerationWatch:
    """One nonblocking inotify watch for a single directory entry."""

    __slots__ = ("_armed", "_changed", "_fd", "_leaf_name", "_read")

    def __init__(self) -> None:
        self._fd: int | None = None
        self._armed = False
        self._changed = False
        self._leaf_name = b""
        self._read = None

    def arm(self, parent: Path, leaf_name: str) -> None:
        try:
            self._arm(parent, leaf_name)
        except Exception:
            self._armed = False

    def changed(self) -> bool:
        if not self._armed:
            return True
        if self._fd is None:
            return self._changed
        try:
            self._drain()
        except Exception:
            self._changed = True
        return self._changed

    def close(self) -> None:
        descriptor = self._fd
        if descriptor is not None and self._armed and not self._changed:
            try:
                self._drain()
            except Exception:
                self._changed = True
        self._fd = None
        if descriptor is None:
            return
        with suppress(OSError):
            os.close(descriptor)

    def _arm(self, parent: Path, leaf_name: str) -> None:
        leaf_bytes = _single_leaf_name(leaf_name)
        init1, add_watch, read = _inotify_functions()
        descriptor = init1(_inotify_fd_flags())
        if descriptor < 0:
            raise OSError(ctypes.get_errno(), "inotify_init1")
        self._fd = descriptor
        self._read = read
        watched = add_watch(descriptor, os.fsencode(parent), _LEAF_WATCH_MASK)
        if watched < 0:
            raise OSError(ctypes.get_errno(), "inotify_add_watch")
        self._leaf_name = leaf_bytes
        self._armed = True

    def _drain(self) -> None:
        descriptor = self._fd
        read = self._read
        if descriptor is None or read is None:
            self._changed = True
            return
        buffer = ctypes.create_string_buffer(_INOTIFY_READ_BYTES)
        while True:
            count = read(descriptor, ctypes.addressof(buffer), len(buffer))
            if count < 0:
                error = ctypes.get_errno()
                if error == errno.EINTR:
                    continue
                if error in {errno.EAGAIN, errno.EWOULDBLOCK}:
                    return
                raise OSError(error, "inotify read")
            if count == 0:
                return
            self._consume(buffer.raw[:count])

    def _consume(self, payload: bytes) -> None:
        offset = 0
        total = len(payload)
        header = _INOTIFY_EVENT.size
        while offset < total:
            if total - offset < header:
                self._changed = True
                return
            _wd, mask, _cookie, name_len = _INOTIFY_EVENT.unpack_from(payload, offset)
            offset += header
            if name_len > total - offset:
                self._changed = True
                return
            raw_name = payload[offset : offset + name_len]
            offset += name_len
            if mask & _IN_Q_OVERFLOW:
                self._changed = True
                continue
            if mask & _IN_IGNORED and not (mask & _LEAF_CHANGE_MASK):
                continue
            name = raw_name.split(b"\x00", 1)[0]
            if name == self._leaf_name and mask & _LEAF_CHANGE_MASK:
                self._changed = True


@contextmanager
def watch_leaf_generation(
    parent: Path,
    leaf_name: str,
) -> Iterator[_LeafGenerationWatch]:
    """Yield an object with changed() -> bool.

    changed() is True when leaf_name was created, deleted, moved, written,
    or had its attributes changed since the watch started, or if the
    inotify queue overflowed. IN_IGNORED from removing our own watch is not
    a change. If inotify cannot be armed, changed() is True (fail closed).
    """

    watch = _LeafGenerationWatch()
    try:
        watch.arm(parent, leaf_name)
        yield watch
    finally:
        watch.close()
