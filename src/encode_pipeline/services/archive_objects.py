"""Publish verified, immutable range-readable copies under an operator-owned root."""

from __future__ import annotations

from contextlib import contextmanager, suppress
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


def _storage_error(exc: OSError) -> ArchiveStoreError:
    if exc.errno in {errno.ELOOP, errno.ENOTDIR}:
        return ArchiveStoreError("archive_integrity")
    return ArchiveStoreError("archive_offline")


def _check_file(fd: int, member: SealedMember) -> None:
    if _file_identity(os.fstat(fd)) != member.identity:
        raise ArchiveStoreError("archive_integrity")


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
        if not 0 <= start <= end <= member.size_bytes:
            raise ArchiveStoreError("archive_invalid")
        fd = self._file.fileno()
        try:
            _check_file(fd, member)
            result = bytearray()
            offset = start
            while offset < end:
                number, within = divmod(offset, CHUNK_BYTES)
                chunk = os.pread(fd, CHUNK_BYTES, number * CHUNK_BYTES)
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
                    os.fchmod(target.fileno(), 0o440)
                    os.fsync(target.fileno())
                # Serialize final link/unlink and identity capture across publishers.
                # A hard link never replaces an existing digest. Removing our staging
                # link before recording ctime keeps identical retries stable.
                fcntl.flock(directory, fcntl.LOCK_EX)
                try:
                    try:
                        os.link(
                            temp, digest, src_dir_fd=directory, dst_dir_fd=directory
                        )
                        created = True
                    except FileExistsError:
                        created = False
                    os.unlink(temp, dir_fd=directory)
                    os.fsync(directory)
                    installed = os.open(digest, _file_flags(), dir_fd=directory)
                finally:
                    fcntl.flock(directory, fcntl.LOCK_UN)
                with os.fdopen(installed, "rb", buffering=0) as handle:
                    info = os.fstat(installed)
                    if (
                        not stat.S_ISREG(info.st_mode)
                        or info.st_size != size
                        or info.st_nlink != 1
                        or info.st_mode & 0o222
                    ):
                        raise ArchiveStoreError("archive_integrity")
                    member = SealedMember(size, digest, chunks, _file_identity(info))
                    if (
                        not created
                        and hashlib.file_digest(handle, "sha256").hexdigest() != digest
                    ):
                        raise ArchiveStoreError("archive_integrity")
                    _check_file(installed, member)
                    return member
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
        block = os.pread(fd, length, offset)
        if len(block) != length:
            raise ArchiveStoreError("archive_integrity")
        output.write(block)
        digest.update(block)
        chunks.append(hashlib.sha256(block).hexdigest())
        offset += length
    if _file_identity(os.fstat(fd)) != _file_identity(before):
        raise ArchiveStoreError("archive_integrity")
    return offset, digest.hexdigest(), tuple(chunks)
