"""Move one small artifact from a scratch root onto a durable root.

Copy, verify, and switch are separate durable stages. A crash resumes at the
last finished stage. Capture seals a byte string only after the scratch file
still satisfies the indexed download revision; the stored checksum is that
byte string. The run, artifact, generation, and revision stored for the
download link are not rewritten. After switch, that same link reads the
verified durable bytes. This is not a catalog migration, and it does not delete
the scratch file or choose a production retention disk.
"""

from __future__ import annotations

from dataclasses import dataclass
import errno
import json
import os
from hashlib import sha256
from pathlib import Path, PurePosixPath
import re
import stat

from encode_pipeline.platform.result_generations import (
    BOUNDED_ARTIFACT_CONTENT_REVISION_MAX_BYTES,
    build_artifact_path_identity,
    source_matches_artifact_revision,
    validate_artifact_generation,
    validate_artifact_path_identity,
    validate_artifact_revision,
)

SCHEME = "artifact-disk-relocation-v1"
MAX_RELOCATION_BYTES = BOUNDED_ARTIFACT_CONTENT_REVISION_MAX_BYTES
PUBLIC_BYTES_MISMATCH = "ARTIFACT_DOWNLOAD_RELOCATION_BYTES_MISMATCH"
PUBLIC_PERMISSION_MISMATCH = "ARTIFACT_DOWNLOAD_RELOCATION_PERMISSION_MISMATCH"
PUBLIC_UNAVAILABLE = "ARTIFACT_DOWNLOAD_RELOCATION_UNAVAILABLE"

_STAGES = frozenset({"captured", "copying", "copied", "verified", "switched"})
_JOURNAL_KEYS = frozenset(
    {
        "artifact_id",
        "content_sha256",
        "generation",
        "mode",
        "relative_path",
        "revision",
        "run_id",
        "scheme",
        "size_bytes",
        "stage",
    }
)
_HEX64 = re.compile(r"^[0-9a-f]{64}$")
_SAFE_OUTPUT_TYPE = re.compile(r"^[A-Za-z][A-Za-z0-9_.-]{0,127}$")
_SAFE_RUN_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,127}$")
_SAFE_ARTIFACT_ID = re.compile(r"^[A-Za-z][A-Za-z0-9_.-]{0,127}$")
_WINDOWS_ABSOLUTE = re.compile(r"^[A-Za-z]:")
_DIRECTORY_FLAGS = os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW | os.O_CLOEXEC
_FILE_READ_FLAGS = os.O_RDONLY | os.O_NOFOLLOW | os.O_CLOEXEC
_FILE_WRITE_FLAGS = os.O_WRONLY | os.O_CREAT | os.O_EXCL | os.O_NOFOLLOW | os.O_CLOEXEC


class ArtifactDiskRelocationError(RuntimeError):
    """Stage or byte failure that must not include a filesystem path."""

    def __init__(self, reason_code: str) -> None:
        super().__init__(reason_code)
        self.reason_code = reason_code


@dataclass(frozen=True)
class RelocationJournal:
    """One artifact's relocation stage and the bytes sealed before the move."""

    stage: str
    run_id: str
    artifact_id: str
    generation: str
    revision: str
    relative_path: str
    size_bytes: int
    content_sha256: str
    mode: int


@dataclass(frozen=True)
class DownloadRelocation:
    """How a saved download link should treat a durable copy."""

    kind: str
    content: bytes = b""
    reason_code: str | None = None


class ArtifactDiskRelocation:
    """Copy one regular file between two absolute roots without relabeling it."""

    def __init__(self, *, source_root: Path, durable_root: Path) -> None:
        self._source_root = _absolute_root(source_root, "source_root")
        self._durable_root = _absolute_root(durable_root, "durable_root")

    @property
    def source_root(self) -> Path:
        return self._source_root

    def capture(
        self,
        *,
        run_id: str,
        artifact_id: str,
        generation: str,
        revision: str,
        relative_path: str,
        output_type: str,
        path_identity: str,
        expected_size: int,
        content: bytes,
        mode: int,
    ) -> RelocationJournal:
        """Seal bytes only when this read still matches the indexed revision.

        ``content`` is the buffer the download revision check already accepted.
        The checksum is that buffer, and only when this read returns the same
        bytes and the same indexed revision. Copying later republishes the
        checksum, not a fresh unscoped read.
        """
        if (
            isinstance(expected_size, bool)
            or not isinstance(expected_size, int)
            or expected_size < 0
            or expected_size > MAX_RELOCATION_BYTES
            or isinstance(mode, bool)
            or not isinstance(mode, int)
            or not 0 <= mode <= 0o777
            or not isinstance(content, bytes)
            or not isinstance(output_type, str)
            or _SAFE_OUTPUT_TYPE.fullmatch(output_type) is None
        ):
            raise ArtifactDiskRelocationError("relocation_identity_invalid")
        try:
            checked_identity = validate_artifact_path_identity(path_identity)
        except (TypeError, ValueError) as exc:
            raise ArtifactDiskRelocationError("relocation_identity_invalid") from exc
        identity = _identity(
            run_id=run_id,
            artifact_id=artifact_id,
            generation=generation,
            revision=revision,
            relative_path=relative_path,
        )
        journal_path = self._journal_path(identity[0], identity[1])
        if _lstat(journal_path) is not None:
            raise ArtifactDiskRelocationError("relocation_journal_exists")
        observed, observed_mode, file_identity, parents = _read_indexed_file(
            _contained(self._source_root / identity[0], PurePosixPath(identity[4]))
        )
        if (
            observed != content
            or observed_mode != mode
            or len(observed) != expected_size
        ):
            raise ArtifactDiskRelocationError("relocation_source_changed")
        try:
            matches = source_matches_artifact_revision(
                output_type=output_type,
                relative_path=identity[4],
                descriptor_identity=":".join(str(value) for value in file_identity),
                content=observed,
                persisted_revision=identity[3],
            )
            observed_identity = build_artifact_path_identity(
                parent_identities=parents,
                file_identity=file_identity,
            )
        except (TypeError, ValueError) as exc:
            raise ArtifactDiskRelocationError("relocation_identity_invalid") from exc
        if not matches:
            raise ArtifactDiskRelocationError("relocation_indexed_revision_mismatch")
        if observed_identity != checked_identity:
            raise ArtifactDiskRelocationError("relocation_path_identity_mismatch")
        journal = RelocationJournal(
            stage="captured",
            run_id=identity[0],
            artifact_id=identity[1],
            generation=identity[2],
            revision=identity[3],
            relative_path=identity[4],
            size_bytes=len(observed),
            content_sha256=sha256(observed).hexdigest(),
            mode=observed_mode,
        )
        return self._write_journal(journal)

    def load(self, run_id: str, artifact_id: str) -> RelocationJournal:
        """Return the current journal for one artifact."""
        safe_run, safe_artifact = _ids(run_id, artifact_id)
        info = _lstat(self._journal_path(safe_run, safe_artifact))
        if info is None:
            raise ArtifactDiskRelocationError("relocation_journal_missing")
        return _parse_journal(self._journal_path(safe_run, safe_artifact), info)

    def payload_file(self, run_id: str, artifact_id: str) -> Path:
        """Return the durable regular-file path named by the journal."""
        return self._payload_path(self.load(run_id, artifact_id))

    def resume(self, *, run_id: str, artifact_id: str) -> RelocationJournal:
        """Advance exactly one unfinished stage, or keep an already switched seal."""
        journal = self.load(run_id, artifact_id)
        if journal.stage in {"captured", "copying"}:
            return self._copy(journal)
        if journal.stage == "copied":
            return self._verify(journal)
        if journal.stage == "verified":
            return self._switch(journal)
        if journal.stage == "switched":
            return journal
        raise ArtifactDiskRelocationError("relocation_stage_invalid")

    def copy(self, *, run_id: str, artifact_id: str) -> RelocationJournal:
        """Rewrite the durable payload from scratch and return to the copied stage."""
        return self._copy(self.load(run_id, artifact_id))

    def _copy(self, journal: RelocationJournal) -> RelocationJournal:
        if journal.stage == "switched":
            raise ArtifactDiskRelocationError("relocation_already_switched")
        self._write_journal(_replace_stage(journal, "copying"))
        content, mode = _read_tree(
            self._source_root / journal.run_id,
            PurePosixPath(journal.relative_path),
        )
        digest = sha256(content).hexdigest()
        if (
            len(content) != journal.size_bytes
            or digest != journal.content_sha256
            or mode != journal.mode
        ):
            raise ArtifactDiskRelocationError("relocation_source_changed")
        self._commit_payload(self._payload_path(journal), content, mode)
        return self._write_journal(_replace_stage(journal, "copied"))

    def _verify(self, journal: RelocationJournal) -> RelocationJournal:
        if journal.stage not in {"copied", "verified"}:
            raise ArtifactDiskRelocationError("relocation_stage_invalid")
        self._require_sealed(journal)
        return self._write_journal(_replace_stage(journal, "verified"))

    def _switch(self, journal: RelocationJournal) -> RelocationJournal:
        if journal.stage != "verified":
            raise ArtifactDiskRelocationError("relocation_stage_invalid")
        self._require_sealed(journal)
        return self._write_journal(_replace_stage(journal, "switched"))

    def _require_sealed(self, journal: RelocationJournal) -> None:
        content, mode = _read_tree(
            self._durable_root / "payload" / journal.run_id,
            PurePosixPath(journal.relative_path),
        )
        digest = sha256(content).hexdigest()
        if len(content) != journal.size_bytes or digest != journal.content_sha256:
            raise ArtifactDiskRelocationError("relocation_payload_mismatch")
        if mode != journal.mode:
            raise ArtifactDiskRelocationError("relocation_permission_mismatch")

    def _commit_payload(self, destination: Path, content: bytes, mode: int) -> None:
        parent = destination.parent
        parent.mkdir(parents=True, exist_ok=True)
        directory = os.open(parent, _DIRECTORY_FLAGS)
        temporary = destination.name + ".tmp"
        file_descriptor: int | None = None
        try:
            try:
                os.unlink(temporary, dir_fd=directory)
            except FileNotFoundError:
                pass
            file_descriptor = os.open(
                temporary,
                _FILE_WRITE_FLAGS,
                mode,
                dir_fd=directory,
            )
            view = memoryview(content)
            while view:
                written = os.write(file_descriptor, view)
                if written <= 0:
                    raise ArtifactDiskRelocationError("relocation_payload_mismatch")
                view = view[written:]
            os.fchmod(file_descriptor, mode)
            os.fsync(file_descriptor)
            os.rename(
                temporary,
                destination.name,
                src_dir_fd=directory,
                dst_dir_fd=directory,
            )
            os.fsync(directory)
        finally:
            if file_descriptor is not None:
                os.close(file_descriptor)
            os.close(directory)

    def _journal_path(self, run_id: str, artifact_id: str) -> Path:
        return self._durable_root / "journals" / run_id / f"{artifact_id}.json"

    def _payload_path(self, journal: RelocationJournal) -> Path:
        relative = PurePosixPath(journal.relative_path)
        return _contained(
            self._durable_root / "payload" / journal.run_id,
            relative,
        )

    def _write_journal(self, journal: RelocationJournal) -> RelocationJournal:
        path = self._journal_path(journal.run_id, journal.artifact_id)
        path.parent.mkdir(parents=True, exist_ok=True)
        payload = json.dumps(
            {
                "artifact_id": journal.artifact_id,
                "content_sha256": journal.content_sha256,
                "generation": journal.generation,
                "mode": journal.mode,
                "relative_path": journal.relative_path,
                "revision": journal.revision,
                "run_id": journal.run_id,
                "scheme": SCHEME,
                "size_bytes": journal.size_bytes,
                "stage": journal.stage,
            },
            allow_nan=False,
            separators=(",", ":"),
            sort_keys=True,
        ).encode("utf-8")
        directory = os.open(path.parent, _DIRECTORY_FLAGS)
        temporary = path.name + ".tmp"
        file_descriptor: int | None = None
        try:
            try:
                os.unlink(temporary, dir_fd=directory)
            except FileNotFoundError:
                pass
            file_descriptor = os.open(
                temporary,
                _FILE_WRITE_FLAGS,
                0o600,
                dir_fd=directory,
            )
            view = memoryview(payload)
            while view:
                written = os.write(file_descriptor, view)
                if written <= 0:
                    raise ArtifactDiskRelocationError("relocation_journal_invalid")
                view = view[written:]
            os.fsync(file_descriptor)
            os.rename(
                temporary,
                path.name,
                src_dir_fd=directory,
                dst_dir_fd=directory,
            )
            os.fsync(directory)
        finally:
            if file_descriptor is not None:
                os.close(file_descriptor)
            os.close(directory)
        stored = self.load(journal.run_id, journal.artifact_id)
        if stored != journal:
            raise ArtifactDiskRelocationError("relocation_journal_invalid")
        return stored


def lookup_switched_download(
    durable_root: Path,
    *,
    run_id: str,
    artifact_id: str,
    generation: str,
    revision: str,
    relative_path: str,
    size_bytes: int,
) -> DownloadRelocation:
    """Read a switched seal, or report that the scratch file is still authoritative."""
    try:
        root = _absolute_root(durable_root, "durable_root")
        safe_run, safe_artifact = _ids(run_id, artifact_id)
        journal_path = root / "journals" / safe_run / f"{safe_artifact}.json"
        info = _lstat(journal_path)
        if info is None:
            return DownloadRelocation(kind="absent")
        journal = _parse_journal(journal_path, info)
        if journal.stage != "switched":
            return DownloadRelocation(kind="absent")
        if (
            journal.run_id != run_id
            or journal.artifact_id != artifact_id
            or journal.generation != generation
            or journal.revision != revision
            or journal.relative_path != relative_path
            or journal.size_bytes != size_bytes
        ):
            return DownloadRelocation(kind="reject", reason_code=PUBLIC_UNAVAILABLE)
        content, mode = _read_tree(
            root / "payload" / journal.run_id,
            PurePosixPath(journal.relative_path),
        )
        digest = sha256(content).hexdigest()
        if len(content) != journal.size_bytes or digest != journal.content_sha256:
            return DownloadRelocation(
                kind="reject",
                reason_code=PUBLIC_BYTES_MISMATCH,
            )
        if mode != journal.mode:
            return DownloadRelocation(
                kind="reject",
                reason_code=PUBLIC_PERMISSION_MISMATCH,
            )
        return DownloadRelocation(kind="serve", content=content)
    except ArtifactDiskRelocationError as exc:
        return DownloadRelocation(
            kind="reject",
            reason_code=_public_reason(exc.reason_code),
        )
    except OSError:
        return DownloadRelocation(kind="reject", reason_code=PUBLIC_UNAVAILABLE)


def _public_reason(reason_code: str) -> str:
    if reason_code == "relocation_payload_mismatch":
        return PUBLIC_BYTES_MISMATCH
    if reason_code == "relocation_permission_mismatch":
        return PUBLIC_PERMISSION_MISMATCH
    return PUBLIC_UNAVAILABLE


def _replace_stage(journal: RelocationJournal, stage: str) -> RelocationJournal:
    return RelocationJournal(
        stage=stage,
        run_id=journal.run_id,
        artifact_id=journal.artifact_id,
        generation=journal.generation,
        revision=journal.revision,
        relative_path=journal.relative_path,
        size_bytes=journal.size_bytes,
        content_sha256=journal.content_sha256,
        mode=journal.mode,
    )


def _absolute_root(value: Path, name: str) -> Path:
    if (
        not isinstance(value, Path)
        or not value.is_absolute()
        or any(part in {"", ".", ".."} for part in value.parts[1:])
    ):
        raise ValueError(f"{name} must be an absolute Path")
    return value


def _ids(run_id: str, artifact_id: str) -> tuple[str, str]:
    if (
        not isinstance(run_id, str)
        or _SAFE_RUN_ID.fullmatch(run_id) is None
        or run_id in {".", ".."}
        or not isinstance(artifact_id, str)
        or _SAFE_ARTIFACT_ID.fullmatch(artifact_id) is None
    ):
        raise ArtifactDiskRelocationError("relocation_identity_invalid")
    return run_id, artifact_id


def _identity(
    *,
    run_id: str,
    artifact_id: str,
    generation: str,
    revision: str,
    relative_path: str,
) -> tuple[str, str, str, str, str]:
    safe_run, safe_artifact = _ids(run_id, artifact_id)
    try:
        validate_artifact_generation(generation)
        validate_artifact_revision(revision)
    except (TypeError, ValueError) as exc:
        raise ArtifactDiskRelocationError("relocation_identity_invalid") from exc
    return (
        safe_run,
        safe_artifact,
        generation,
        revision,
        _validate_relative_path(relative_path),
    )


def _validate_relative_path(value: object) -> str:
    if (
        not isinstance(value, str)
        or not value
        or len(value) > 2048
        or "\x00" in value
        or "\\" in value
        or value.startswith(("/", "~"))
        or _WINDOWS_ABSOLUTE.match(value) is not None
    ):
        raise ArtifactDiskRelocationError("relocation_identity_invalid")
    path = PurePosixPath(value)
    if (
        path.as_posix() != value
        or path.is_absolute()
        or len(path.parts) < 2
        or len(path.parts) > 64
        or path.parts[0] != "results"
        or any(component in {"", ".", ".."} for component in value.split("/"))
    ):
        raise ArtifactDiskRelocationError("relocation_identity_invalid")
    return value


def _contained(root: Path, relative: PurePosixPath) -> Path:
    candidate = root.joinpath(*relative.parts)
    try:
        candidate.relative_to(root)
    except ValueError as exc:
        raise ArtifactDiskRelocationError("relocation_identity_invalid") from exc
    return candidate


def _lstat(path: Path) -> os.stat_result | None:
    try:
        return os.lstat(path)
    except FileNotFoundError:
        return None


def _parse_journal(path: Path, info: os.stat_result) -> RelocationJournal:
    if stat.S_ISLNK(info.st_mode) or not stat.S_ISREG(info.st_mode):
        raise ArtifactDiskRelocationError("relocation_journal_invalid")
    if info.st_size > 8192:
        raise ArtifactDiskRelocationError("relocation_journal_invalid")
    descriptor = os.open(path, _FILE_READ_FLAGS)
    try:
        raw = _read_exact(descriptor, info.st_size)
    finally:
        os.close(descriptor)
    try:
        payload = json.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ArtifactDiskRelocationError("relocation_journal_invalid") from exc
    if not isinstance(payload, dict) or set(payload) != _JOURNAL_KEYS:
        raise ArtifactDiskRelocationError("relocation_journal_invalid")
    if payload["scheme"] != SCHEME or payload["stage"] not in _STAGES:
        raise ArtifactDiskRelocationError("relocation_journal_invalid")
    try:
        run_id, artifact_id, generation, revision, relative_path = _identity(
            run_id=payload["run_id"],
            artifact_id=payload["artifact_id"],
            generation=payload["generation"],
            revision=payload["revision"],
            relative_path=payload["relative_path"],
        )
    except ArtifactDiskRelocationError:
        raise
    expected_name = path.parent.name
    if run_id != expected_name or path.name != f"{artifact_id}.json":
        raise ArtifactDiskRelocationError("relocation_journal_invalid")
    size_bytes = payload["size_bytes"]
    mode = payload["mode"]
    digest = payload["content_sha256"]
    if (
        isinstance(size_bytes, bool)
        or not isinstance(size_bytes, int)
        or size_bytes < 0
        or size_bytes > MAX_RELOCATION_BYTES
        or isinstance(mode, bool)
        or not isinstance(mode, int)
        or mode < 0
        or mode > 0o777
        or not isinstance(digest, str)
        or _HEX64.fullmatch(digest) is None
    ):
        raise ArtifactDiskRelocationError("relocation_journal_invalid")
    return RelocationJournal(
        stage=payload["stage"],
        run_id=run_id,
        artifact_id=artifact_id,
        generation=generation,
        revision=revision,
        relative_path=relative_path,
        size_bytes=size_bytes,
        content_sha256=digest,
        mode=mode,
    )


def _descriptor_tuple(
    info: os.stat_result,
) -> tuple[int, int, int, int, int, int, int, int, int]:
    return (
        info.st_dev,
        info.st_ino,
        info.st_mode,
        info.st_nlink,
        info.st_uid,
        info.st_gid,
        info.st_size,
        info.st_mtime_ns,
        info.st_ctime_ns,
    )


def _read_indexed_file(
    absolute: Path,
) -> tuple[
    bytes,
    int,
    tuple[int, int, int, int, int, int, int, int, int],
    tuple[tuple[int, int, int], ...],
]:
    """Read one regular file from ``/`` and return bytes, mode, and identity."""
    if not absolute.is_absolute():
        raise ArtifactDiskRelocationError("relocation_source_unavailable")
    parts = absolute.parts[1:]
    if not parts:
        raise ArtifactDiskRelocationError("relocation_source_unavailable")
    descriptors: list[int] = []
    infos: list[os.stat_result] = []
    try:
        current = os.open("/", _DIRECTORY_FLAGS)
        descriptors.append(current)
        for index, component in enumerate(parts):
            final = index == len(parts) - 1
            flags = _FILE_READ_FLAGS | os.O_NONBLOCK if final else _DIRECTORY_FLAGS
            try:
                opened = os.open(component, flags, dir_fd=current)
            except OSError as exc:
                if exc.errno in {
                    errno.ENOENT,
                    errno.ELOOP,
                    errno.ENOTDIR,
                    errno.ENXIO,
                    errno.ENODEV,
                    errno.ESTALE,
                }:
                    raise ArtifactDiskRelocationError(
                        "relocation_source_unavailable"
                    ) from exc
                raise
            descriptors.append(opened)
            current = opened
            info = os.fstat(opened)
            if final:
                if not stat.S_ISREG(info.st_mode):
                    raise ArtifactDiskRelocationError("relocation_source_unavailable")
            elif not stat.S_ISDIR(info.st_mode):
                raise ArtifactDiskRelocationError("relocation_source_unavailable")
            infos.append(info)
        file_info = infos[-1]
        if file_info.st_size > MAX_RELOCATION_BYTES:
            raise ArtifactDiskRelocationError("relocation_too_large")
        file_identity = _descriptor_tuple(file_info)
        content = _read_exact(current, file_info.st_size)
        os.lseek(current, 0, os.SEEK_SET)
        if _descriptor_tuple(os.fstat(current)) != file_identity:
            raise ArtifactDiskRelocationError("relocation_source_changed")
        again = _read_exact(current, file_info.st_size)
        if again != content or _descriptor_tuple(os.fstat(current)) != file_identity:
            raise ArtifactDiskRelocationError("relocation_source_changed")
        for descriptor, info in zip(descriptors[1:-1], infos[:-1], strict=True):
            current_info = os.fstat(descriptor)
            if (current_info.st_dev, current_info.st_ino, current_info.st_mode) != (
                info.st_dev,
                info.st_ino,
                info.st_mode,
            ):
                raise ArtifactDiskRelocationError("relocation_source_changed")
        parents = tuple((info.st_dev, info.st_ino, info.st_mode) for info in infos[:-1])
        return content, stat.S_IMODE(file_info.st_mode), file_identity, parents
    except ArtifactDiskRelocationError:
        raise
    except OSError as exc:
        raise ArtifactDiskRelocationError("relocation_source_unavailable") from exc
    finally:
        for descriptor in reversed(descriptors):
            os.close(descriptor)


def _read_tree(root: Path, relative: PurePosixPath) -> tuple[bytes, int]:
    """Return the bytes and permission bits of one regular file, without following links."""
    if not root.is_absolute():
        raise ArtifactDiskRelocationError("relocation_source_unavailable")
    descriptors: list[int] = []
    try:
        current = os.open(root, _DIRECTORY_FLAGS)
        descriptors.append(current)
        info = os.fstat(current)
        if not stat.S_ISDIR(info.st_mode):
            raise ArtifactDiskRelocationError("relocation_source_unavailable")
        for index, component in enumerate(relative.parts):
            final = index == len(relative.parts) - 1
            try:
                opened = os.open(
                    component,
                    _FILE_READ_FLAGS if final else _DIRECTORY_FLAGS,
                    dir_fd=current,
                )
            except OSError as exc:
                if exc.errno in {errno.ENOENT, errno.ELOOP, errno.ENOTDIR}:
                    raise ArtifactDiskRelocationError(
                        "relocation_source_unavailable"
                    ) from exc
                raise
            descriptors.append(opened)
            current = opened
            info = os.fstat(opened)
            if not final and not stat.S_ISDIR(info.st_mode):
                raise ArtifactDiskRelocationError("relocation_source_unavailable")
            if final and not stat.S_ISREG(info.st_mode):
                raise ArtifactDiskRelocationError("relocation_source_unavailable")
        if info.st_size > MAX_RELOCATION_BYTES:
            raise ArtifactDiskRelocationError("relocation_too_large")
        first = _read_exact(current, info.st_size)
        os.lseek(current, 0, os.SEEK_SET)
        second_info = os.fstat(current)
        if stat.S_IMODE(second_info.st_mode) != stat.S_IMODE(info.st_mode):
            raise ArtifactDiskRelocationError("relocation_permission_mismatch")
        if second_info.st_size != info.st_size:
            raise ArtifactDiskRelocationError("relocation_payload_mismatch")
        second = _read_exact(current, second_info.st_size)
        if second != first:
            raise ArtifactDiskRelocationError("relocation_payload_mismatch")
        return first, stat.S_IMODE(info.st_mode)
    except ArtifactDiskRelocationError:
        raise
    except OSError as exc:
        raise ArtifactDiskRelocationError("relocation_source_unavailable") from exc
    finally:
        for descriptor in reversed(descriptors):
            os.close(descriptor)


def _read_exact(descriptor: int, size: int) -> bytes:
    chunks: list[bytes] = []
    remaining = size
    while remaining:
        chunk = os.read(descriptor, min(65536, remaining))
        if not chunk:
            break
        chunks.append(chunk)
        remaining -= len(chunk)
    content = b"".join(chunks)
    if len(content) != size:
        raise ArtifactDiskRelocationError("relocation_payload_mismatch")
    return content
