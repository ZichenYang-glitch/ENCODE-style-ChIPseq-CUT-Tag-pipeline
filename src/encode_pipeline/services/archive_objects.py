"""Publish verified, immutable range-readable copies under an operator-owned root."""

from __future__ import annotations

from contextlib import contextmanager, suppress
import ctypes
import errno
import fcntl
import hashlib
import os
from pathlib import Path
import stat
from uuid import uuid4

from encode_pipeline.platform.archived_tracks import (
    ArchiveStoreError,
    CHUNK_BYTES,
    SealedMember,
)
from encode_pipeline.services.input_file_access import (
    _file_flags,
    _file_identity,
    _open_absolute_directory,
)

# Linux renameat2(2): publish a new directory entry without replacing one.
_RENAME_NOREPLACE = 1
_OBJECT_MODE = 0o440


def _storage_error(exc: OSError) -> ArchiveStoreError:
    if exc.errno in {errno.ELOOP, errno.ENOTDIR}:
        return ArchiveStoreError("archive_integrity")
    return ArchiveStoreError("archive_offline")


def _check_file(fd: int, member: SealedMember) -> None:
    if _file_identity(os.fstat(fd)) != member.identity:
        raise ArchiveStoreError("archive_integrity")


def _read_exact(fd: int, length: int, offset: int) -> bytes:
    """Accumulate POSIX short reads. Empty before ``length`` is real EOF."""
    parts = []
    remaining = length
    position = offset
    while remaining:
        try:
            block = os.pread(fd, remaining, position)
        except InterruptedError:
            continue
        if not block:
            raise ArchiveStoreError("archive_integrity")
        parts.append(block)
        remaining -= len(block)
        position += len(block)
    return b"".join(parts)


def _rename_noreplace(directory: int, source: str, destination: str) -> bool:
    """Publish ``source`` as ``destination`` in one call.

    Return False when that name is already taken. Do not replace a file,
    symlink, or directory.
    """
    libc = ctypes.CDLL(None, use_errno=True)
    try:
        renameat2 = libc.renameat2
    except AttributeError:
        raise OSError(errno.ENOSYS, "renameat2") from None
    renameat2.argtypes = (
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_uint,
    )
    renameat2.restype = ctypes.c_int
    ctypes.set_errno(0)
    result = renameat2(
        directory,
        os.fsencode(source),
        directory,
        os.fsencode(destination),
        _RENAME_NOREPLACE,
    )
    if result == 0:
        return True
    error = ctypes.get_errno()
    if error in {errno.EEXIST, errno.ENOTEMPTY, errno.EISDIR}:
        return False
    raise OSError(error or errno.EIO, "renameat2")


class ArchiveReader:
    """One pinned descriptor. Return slices of the exact chunks verified."""

    def __init__(self, fd: int, member: SealedMember):
        self._file = os.fdopen(fd, "rb", buffering=0)
        self._member = member

    def __enter__(self) -> ArchiveReader:
        return self

    def __exit__(self, *_exc) -> None:
        self.close()

    def close(self) -> None:
        self._file.close()

    def read(self, start: int, end: int) -> bytes:
        member = self._member
        if (
            type(start) is not int
            or type(end) is not int
            or not 0 <= start <= end <= member.size_bytes
        ):
            raise ArchiveStoreError("archive_invalid")
        if self._file.closed:
            raise ArchiveStoreError("archive_invalid")
        try:
            fd = self._file.fileno()
        except ValueError:
            raise ArchiveStoreError("archive_invalid") from None
        try:
            _check_file(fd, member)
            result = bytearray()
            offset = start
            while offset < end:
                number, within = divmod(offset, CHUNK_BYTES)
                chunk_start = number * CHUNK_BYTES
                expected = min(CHUNK_BYTES, member.size_bytes - chunk_start)
                chunk = _read_exact(fd, expected, chunk_start)
                if hashlib.sha256(chunk).hexdigest() != member.chunks[number]:
                    raise ArchiveStoreError("archive_integrity")
                length = min(len(chunk) - within, end - offset)
                if length <= 0:
                    raise ArchiveStoreError("archive_integrity")
                result.extend(chunk[within : within + length])
                offset += length
            _check_file(fd, member)
            return bytes(result)
        except OSError as exc:
            raise _storage_error(exc) from None


class ArchiveObjects:
    def __init__(self, root: Path | None):
        self._root = root

    @contextmanager
    def _directory(self):
        if self._root is None:
            raise ArchiveStoreError("archive_offline")
        try:
            fd = _open_absolute_directory(self._root)
            try:
                yield fd
            finally:
                os.close(fd)
        except OSError as exc:
            raise _storage_error(exc) from None

    def open_member(self, member: SealedMember) -> ArchiveReader:
        with self._directory() as directory:
            descriptor = os.open(member.sha256, _file_flags(), dir_fd=directory)
            try:
                _check_file(descriptor, member)
                return ArchiveReader(descriptor, member)
            except BaseException:
                os.close(descriptor)
                raise

    def seal_fd(self, fd: int, expected_sha256: str | None = None) -> SealedMember:
        """Copy a stable source from offset zero. Never replace an existing object."""
        with self._directory() as directory:
            temp = f".sealing-{uuid4().hex}"
            new_fd = os.open(
                temp,
                os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC,
                0o600,
                dir_fd=directory,
            )
            try:
                with os.fdopen(new_fd, "wb") as target:
                    size, digest, chunks = _copy(fd, target)
                    if expected_sha256 is not None and digest != expected_sha256:
                        raise ArchiveStoreError("archive_integrity")
                    target.flush()
                    os.fchmod(target.fileno(), _OBJECT_MODE)
                    os.fsync(target.fileno())
                # The digest name appears in one rename. A crash cannot leave
                # the staging name as a second hard link, and an existing
                # file, symlink, or directory stays in place for a later retry.
                fcntl.flock(directory, fcntl.LOCK_EX)
                try:
                    created = _rename_noreplace(directory, temp, digest)
                    if not created:
                        os.unlink(temp, dir_fd=directory)
                    os.fsync(directory)
                    installed = os.open(digest, _file_flags(), dir_fd=directory)
                finally:
                    fcntl.flock(directory, fcntl.LOCK_UN)
                try:
                    info = os.fstat(installed)
                    if (
                        not stat.S_ISREG(info.st_mode)
                        or info.st_size != size
                        or info.st_nlink != 1
                        or stat.S_IMODE(info.st_mode) != _OBJECT_MODE
                    ):
                        raise ArchiveStoreError("archive_integrity")
                    member = SealedMember(size, digest, chunks, _file_identity(info))
                    if not created:
                        os.lseek(installed, 0, os.SEEK_SET)
                        with os.fdopen(
                            installed, "rb", buffering=0, closefd=False
                        ) as handle:
                            observed = hashlib.file_digest(handle, "sha256").hexdigest()
                        if observed != digest:
                            raise ArchiveStoreError("archive_integrity")
                    _check_file(installed, member)
                    return member
                finally:
                    os.close(installed)
            finally:
                with suppress(FileNotFoundError):
                    os.unlink(temp, dir_fd=directory)


def _copy(fd: int, output) -> tuple[int, str, tuple[str, ...]]:
    before = os.fstat(fd)
    if not stat.S_ISREG(before.st_mode) or before.st_size <= 0:
        raise ArchiveStoreError("archive_invalid")
    digest = hashlib.sha256()
    chunks = []
    offset = 0
    while offset < before.st_size:
        length = min(CHUNK_BYTES, before.st_size - offset)
        block = _read_exact(fd, length, offset)
        output.write(block)
        digest.update(block)
        chunks.append(hashlib.sha256(block).hexdigest())
        offset += length
    if _file_identity(os.fstat(fd)) != _file_identity(before):
        raise ArchiveStoreError("archive_integrity")
    return offset, digest.hexdigest(), tuple(chunks)
