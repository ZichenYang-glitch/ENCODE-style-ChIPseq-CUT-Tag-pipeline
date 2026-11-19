"""Copy one registered file into a response-private regular file."""

from __future__ import annotations

from contextlib import suppress
from hashlib import sha256
import os
from pathlib import Path
import stat
import tempfile

from encode_pipeline.platform.input_registry import validate_input_file_relative_path
from encode_pipeline.services.input_file_access import (
    _directory_flags,
    _file_flags,
    _open_absolute_directory,
)


_READ_CHUNK_BYTES = 1024 * 1024


class RegisteredResultBytesError(RuntimeError):
    """The registered bytes could not be sealed for this response."""


def seal_private_copy(
    *,
    root: Path,
    relative_path: str,
    expected_size: int,
    expected_sha256: str,
    response_root: Path,
) -> tuple[Path, Path]:
    """Hash the source into a 0600 file and return that file and its directory.

    The returned file is the only byte source for a successful response.
    Later changes to the source do not alter it. The caller deletes both
    paths when the response ends or when anything fails after this returns.
    """

    directory: Path | None = None
    copy_path: Path | None = None
    descriptors: list[int] = []
    copy_descriptor: int | None = None
    sealed = False
    try:
        safe_path = validate_input_file_relative_path(relative_path)
        directory = Path(
            tempfile.mkdtemp(prefix="helixweave-rreg-response-", dir=response_root)
        )
        copy_path = directory / "body"
        os.chmod(directory, 0o700)
        root_descriptor = _open_absolute_directory(root)
        descriptors.append(root_descriptor)
        current = root_descriptor
        components = safe_path.split("/")
        for component in components[:-1]:
            descriptor = os.open(component, _directory_flags(), dir_fd=current)
            descriptors.append(descriptor)
            current = descriptor
            if not stat.S_ISDIR(os.fstat(descriptor).st_mode):
                raise OSError("unsafe directory")
        source = os.open(components[-1], _file_flags(), dir_fd=current)
        descriptors.append(source)
        info = os.fstat(source)
        if not stat.S_ISREG(info.st_mode):
            raise OSError("unsafe file")
        copy_descriptor = os.open(
            copy_path,
            os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_CLOEXEC,
            0o600,
        )
        os.fchmod(copy_descriptor, 0o600)
        digest = sha256()
        size_bytes = 0
        while True:
            chunk = os.read(source, _READ_CHUNK_BYTES)
            if not chunk:
                break
            size_bytes += len(chunk)
            digest.update(chunk)
            _write_all(copy_descriptor, chunk)
        os.fsync(copy_descriptor)
        if size_bytes != expected_size or digest.hexdigest() != expected_sha256:
            raise RegisteredResultBytesError()
        sealed = True
        return copy_path, directory
    except RegisteredResultBytesError:
        raise
    except (OSError, TypeError, ValueError):
        raise RegisteredResultBytesError() from None
    finally:
        if copy_descriptor is not None:
            with suppress(OSError):
                os.close(copy_descriptor)
        for descriptor in reversed(descriptors):
            with suppress(OSError):
                os.close(descriptor)
        if not sealed and directory is not None and copy_path is not None:
            discard_private_copy(copy_path, directory)


def discard_private_copy(copy_path: Path, directory: Path) -> None:
    """Remove one response copy without reporting its location."""

    with suppress(OSError):
        os.unlink(copy_path)
    with suppress(OSError):
        os.rmdir(directory)


def _write_all(descriptor: int, chunk: bytes) -> None:
    view = memoryview(chunk)
    while view:
        written = os.write(descriptor, view)
        if written <= 0:
            raise OSError("short write")
        view = view[written:]
