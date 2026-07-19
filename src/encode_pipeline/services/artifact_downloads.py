"""Descriptor-safe streaming of persisted run artifacts."""

from __future__ import annotations

from collections.abc import Iterator, Mapping
from dataclasses import dataclass
import errno
import os
from pathlib import Path, PurePosixPath
import re
import stat
from threading import Lock
from urllib.parse import quote

from encode_pipeline.platform.result_generations import (
    ARTIFACT_PATH_IDENTITY_METADATA_KEY,
    BOUNDED_ARTIFACT_CONTENT_REVISION_MAX_BYTES,
    build_artifact_content_revision,
    build_artifact_descriptor_revision,
    build_artifact_path_identity,
    validate_artifact_generation,
    validate_artifact_path_identity,
    validate_artifact_revision,
)
from encode_pipeline.platform.results import Issue, Result
from encode_pipeline.platform.runs import RunArtifactRef
from encode_pipeline.services.run_repositories import ResultGenerationChangedError
from encode_pipeline.services.runs import RunService


_CHUNK_SIZE = 64 * 1024
_MAX_ARTIFACT_BYTES = 2**63 - 1
_MAX_PATH_COMPONENTS = 64
_SAFE_RUN_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,127}$")
_SAFE_ARTIFACT_ID = re.compile(r"^[A-Za-z][A-Za-z0-9_.-]{0,127}$")
_SAFE_OUTPUT_TYPE = re.compile(r"^[A-Za-z][A-Za-z0-9_.-]{0,127}$")
_SAFE_MIME_TYPE = re.compile(
    r"^[A-Za-z0-9][A-Za-z0-9!#$&^_.+-]{0,126}/"
    r"[A-Za-z0-9][A-Za-z0-9!#$&^_.+-]{0,126}$"
)
_WINDOWS_ABSOLUTE = re.compile(r"^[A-Za-z]:")
_CONTROL_CHARACTER = re.compile(r"[\x00-\x1f\x7f]")
_CONFLICT_ERRNOS = frozenset(
    {
        errno.ENOENT,
        errno.ENOTDIR,
        errno.ELOOP,
        errno.ENXIO,
        errno.ENODEV,
        errno.ESTALE,
    }
)


class ArtifactDownloadStreamError(RuntimeError):
    """Public-safe signal that a prepared artifact stream became invalid."""

    def __init__(self) -> None:
        super().__init__("Artifact download could not continue safely.")


class _ArtifactSourceConflict(RuntimeError):
    def __init__(self, reason_code: str) -> None:
        super().__init__(reason_code)
        self.reason_code = reason_code


@dataclass(frozen=True)
class _Fingerprint:
    device: int
    inode: int
    mode: int
    links: int
    user: int
    group: int
    size: int
    modified_ns: int
    changed_ns: int


@dataclass(frozen=True)
class _OpenedNode:
    descriptor: int
    parent_index: int | None
    name: str
    fingerprint: _Fingerprint


class ArtifactDownloadPlan:
    """An already-open, bounded artifact stream with idempotent cleanup."""

    def __init__(
        self,
        *,
        nodes: tuple[_OpenedNode, ...] = (),
        content: bytes | None = None,
        size_bytes: int,
        media_type: str,
        filename: str,
        content_disposition: str,
    ) -> None:
        if (not nodes) == (content is None):
            raise ValueError("download plan requires exactly one byte source")
        self._nodes = nodes
        self._final_descriptor = nodes[-1].descriptor if nodes else None
        self._content = content
        self.size_bytes = size_bytes
        self.media_type = media_type
        self.filename = filename
        self.content_disposition = content_disposition
        self._lock = Lock()
        self._closed = False
        self._started = False

    @property
    def closed(self) -> bool:
        """Return whether every owned descriptor has been released."""
        with self._lock:
            return self._closed

    def iter_bytes(self) -> Iterator[bytes]:
        """Yield exactly the persisted byte count while identity stays stable."""
        with self._lock:
            if self._closed or self._started:
                raise ArtifactDownloadStreamError()
            self._started = True
        remaining = self.size_bytes
        try:
            if self._content is not None:
                for offset in range(0, len(self._content), _CHUNK_SIZE):
                    yield self._content[offset : offset + _CHUNK_SIZE]
                return
            assert self._final_descriptor is not None
            while remaining:
                _verify_descriptor_chain(self._nodes)
                chunk = os.read(self._final_descriptor, min(_CHUNK_SIZE, remaining))
                if not chunk:
                    raise ArtifactDownloadStreamError()
                _verify_descriptor_chain(self._nodes)
                remaining -= len(chunk)
                yield chunk
            _verify_descriptor_chain(self._nodes)
            if os.read(self._final_descriptor, 1):
                raise ArtifactDownloadStreamError()
            _verify_descriptor_chain(self._nodes)
        except ArtifactDownloadStreamError:
            raise
        except (OSError, TypeError, ValueError):
            raise ArtifactDownloadStreamError() from None
        finally:
            self.close()

    def close(self) -> None:
        """Release all descriptors once, in child-to-parent order."""
        with self._lock:
            if self._closed:
                return
            self._closed = True
            descriptors = tuple(node.descriptor for node in self._nodes)
        for descriptor in reversed(descriptors):
            try:
                os.close(descriptor)
            except OSError:
                pass


class _PreparedArtifactDownloadSource:
    """Keep the verified source open until the final persistence guard passes."""

    def __init__(
        self,
        *,
        nodes: tuple[_OpenedNode, ...],
        content: bytes | None,
        size_bytes: int,
        media_type: str,
        filename: str,
        content_disposition: str,
    ) -> None:
        if not nodes:
            raise ValueError("prepared source requires an open descriptor chain")
        self._nodes = nodes
        self._content = content
        self._size_bytes = size_bytes
        self._media_type = media_type
        self._filename = filename
        self._content_disposition = content_disposition

    def seal(self) -> ArtifactDownloadPlan:
        """Verify the source at the linearization point and transfer ownership."""
        if not self._nodes:
            raise ArtifactDownloadStreamError()
        _verify_descriptor_chain(self._nodes)
        try:
            if self._content is not None:
                plan = ArtifactDownloadPlan(
                    content=self._content,
                    size_bytes=self._size_bytes,
                    media_type=self._media_type,
                    filename=self._filename,
                    content_disposition=self._content_disposition,
                )
                self.close()
                return plan
            plan = ArtifactDownloadPlan(
                nodes=self._nodes,
                size_bytes=self._size_bytes,
                media_type=self._media_type,
                filename=self._filename,
                content_disposition=self._content_disposition,
            )
        except BaseException:
            self.close()
            raise
        self._nodes = ()
        return plan

    def close(self) -> None:
        """Release an unsealed source in child-to-parent order."""
        nodes, self._nodes = self._nodes, ()
        for descriptor in reversed(tuple(node.descriptor for node in nodes)):
            try:
                os.close(descriptor)
            except OSError:
                pass


class ArtifactDownloadService:
    """Resolve one persisted artifact into a descriptor-owned stream plan."""

    def __init__(self, *, run_service: RunService, workspace_root: Path) -> None:
        if not isinstance(run_service, RunService):
            raise ValueError("run_service must be a RunService")
        if (
            not isinstance(workspace_root, Path)
            or not workspace_root.is_absolute()
            or any(part in {"", ".", ".."} for part in workspace_root.parts[1:])
        ):
            raise ValueError("workspace_root must be an absolute Path")
        self._run_service = run_service
        self._workspace_root = workspace_root

    def prepare(
        self,
        run_id: str,
        artifact_id: str,
        *,
        expected_generation: str,
        expected_revision: str,
    ) -> Result[ArtifactDownloadPlan]:
        """Return a verified stream plan without exposing a workspace path."""
        try:
            validate_artifact_generation(expected_generation)
            validate_artifact_revision(expected_revision)
        except (TypeError, ValueError):
            return self._failure(
                "RUN_ARTIFACT_DOWNLOAD_CONFLICT",
                "Artifact content is no longer available as indexed.",
                path="generation",
                reason_code="ARTIFACT_DOWNLOAD_EXPECTATION_INVALID",
            )
        try:
            self._run_service.get_run(run_id)
        except KeyError:
            return self._failure(
                "RUN_NOT_FOUND",
                "Run was not found.",
                path="run_id",
            )
        except Exception:
            return self._unavailable("ARTIFACT_DOWNLOAD_REPOSITORY_UNAVAILABLE")

        try:
            _generation, artifact = self._run_service.get_artifact_at_generation(
                run_id,
                artifact_id,
                expected_generation=expected_generation,
            )
        except ResultGenerationChangedError:
            return self._generation_conflict()
        except KeyError:
            return self._failure(
                "RUN_ARTIFACT_NOT_FOUND",
                "Artifact was not found for this run.",
                path="artifact_id",
            )
        except (AttributeError, TypeError, ValueError):
            return self._data_invalid()
        except Exception:
            return self._unavailable("ARTIFACT_DOWNLOAD_REPOSITORY_UNAVAILABLE")

        if artifact.revision != expected_revision:
            return self._failure(
                "RUN_ARTIFACT_DOWNLOAD_CONFLICT",
                "Artifact content is no longer available as indexed.",
                path="revision",
                reason_code="ARTIFACT_DOWNLOAD_REVISION_CHANGED",
            )

        try:
            validated = self._validate_reference(
                run_id,
                artifact_id,
                artifact,
                expected_revision=expected_revision,
            )
        except (AttributeError, TypeError, ValueError):
            return self._data_invalid()

        (
            relative_path,
            expected_size,
            output_type,
            persisted_revision,
            path_identity,
            media_type,
            filename,
            disposition,
        ) = validated
        prepared: _PreparedArtifactDownloadSource | None = None
        try:
            prepared = self._open_plan(
                run_id=run_id,
                relative_path=relative_path,
                expected_size=expected_size,
                output_type=output_type,
                persisted_revision=persisted_revision,
                persisted_path_identity=path_identity,
                media_type=media_type,
                filename=filename,
                content_disposition=disposition,
            )
            try:
                _current_generation, current = (
                    self._run_service.get_artifact_at_generation(
                        run_id,
                        artifact_id,
                        expected_generation=expected_generation,
                    )
                )
                if current.revision != expected_revision:
                    raise ResultGenerationChangedError("artifact revision changed")
            except ResultGenerationChangedError:
                return self._generation_conflict()
            except KeyError:
                return self._generation_conflict()
            except (AttributeError, TypeError, ValueError):
                return self._data_invalid()
            except Exception:
                return self._unavailable("ARTIFACT_DOWNLOAD_REPOSITORY_UNAVAILABLE")
            plan = prepared.seal()
        except _ArtifactSourceConflict as exc:
            return self._failure(
                "RUN_ARTIFACT_DOWNLOAD_CONFLICT",
                "Artifact content is no longer available as indexed.",
                path="artifact_id",
                reason_code=exc.reason_code,
            )
        except ArtifactDownloadStreamError:
            return self._failure(
                "RUN_ARTIFACT_DOWNLOAD_CONFLICT",
                "Artifact content is no longer available as indexed.",
                path="artifact_id",
                reason_code="ARTIFACT_DOWNLOAD_SOURCE_CHANGED",
            )
        except Exception:
            return self._unavailable("ARTIFACT_DOWNLOAD_IO_UNAVAILABLE")
        finally:
            if prepared is not None:
                prepared.close()
        return Result.success(plan)

    @staticmethod
    def _validate_reference(
        run_id: str,
        artifact_id: str,
        artifact: RunArtifactRef,
        *,
        expected_revision: str,
    ) -> tuple[str, int, str, str, str, str, str, str]:
        if (
            not isinstance(run_id, str)
            or _SAFE_RUN_ID.fullmatch(run_id) is None
            or run_id in {".", ".."}
            or not isinstance(artifact_id, str)
            or _SAFE_ARTIFACT_ID.fullmatch(artifact_id) is None
            or not isinstance(artifact, RunArtifactRef)
            or artifact.run_id != run_id
            or artifact.artifact_id != artifact_id
            or artifact.artifact_type != "file"
        ):
            raise ValueError("artifact identity is invalid")
        expected_uri = f"run://runs/{quote(run_id, safe='')}/artifacts/{artifact_id}"
        if artifact.uri != expected_uri:
            raise ValueError("artifact URI is invalid")
        if not isinstance(artifact.metadata, Mapping):
            raise ValueError("artifact metadata is invalid")
        relative_path = _validate_relative_path(artifact.metadata.get("relative_path"))
        expected_size = _validate_size(artifact.metadata.get("size_bytes"))
        output_type = _validate_output_type(artifact.metadata.get("output_type"))
        path_identity_value = artifact.metadata.get(ARTIFACT_PATH_IDENTITY_METADATA_KEY)
        if not isinstance(path_identity_value, str):
            raise ValueError("artifact path identity is invalid")
        path_identity = validate_artifact_path_identity(path_identity_value)
        filename = _validate_filename(artifact.name, relative_path)
        media_type = artifact.mime_type or "application/octet-stream"
        if _SAFE_MIME_TYPE.fullmatch(media_type) is None:
            raise ValueError("artifact MIME type is invalid")
        disposition = _content_disposition(filename, artifact_id)
        return (
            relative_path,
            expected_size,
            output_type,
            artifact.revision,
            path_identity,
            media_type,
            filename,
            disposition,
        )

    def _open_plan(
        self,
        *,
        run_id: str,
        relative_path: str,
        expected_size: int,
        output_type: str,
        persisted_revision: str,
        persisted_path_identity: str,
        media_type: str,
        filename: str,
        content_disposition: str,
    ) -> _PreparedArtifactDownloadSource:
        required_flags = ("O_DIRECTORY", "O_NOFOLLOW", "O_CLOEXEC", "O_NONBLOCK")
        if any(not hasattr(os, name) for name in required_flags):
            raise OSError("safe descriptor flags are unavailable")
        directory_flags = os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW | os.O_CLOEXEC
        final_flags = os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK | os.O_CLOEXEC
        workspace = self._workspace_root / run_id
        components = (*workspace.parts[1:], *PurePosixPath(relative_path).parts)
        if not components or len(components) > _MAX_PATH_COMPONENTS * 2:
            raise ValueError("artifact path is too deep")

        nodes: list[_OpenedNode] = []
        owned_descriptors: list[int] = []
        try:
            root_descriptor = os.open("/", directory_flags)
            _own_descriptor(owned_descriptors, root_descriptor)
            nodes.append(
                _OpenedNode(
                    descriptor=root_descriptor,
                    parent_index=None,
                    name="/",
                    fingerprint=_fingerprint(os.fstat(root_descriptor)),
                )
            )
            current = root_descriptor
            for index, component in enumerate(components):
                final = index == len(components) - 1
                try:
                    descriptor = os.open(
                        component,
                        final_flags if final else directory_flags,
                        dir_fd=current,
                    )
                except OSError as exc:
                    if exc.errno in _CONFLICT_ERRNOS:
                        raise _ArtifactSourceConflict(
                            "ARTIFACT_DOWNLOAD_SOURCE_MISSING_OR_UNSAFE"
                        ) from None
                    raise
                _own_descriptor(owned_descriptors, descriptor)
                info = os.fstat(descriptor)
                node = _OpenedNode(
                    descriptor=descriptor,
                    parent_index=len(nodes) - 1,
                    name=component,
                    fingerprint=_fingerprint(info),
                )
                nodes.append(node)
                current = descriptor
                if not final and not stat.S_ISDIR(info.st_mode):
                    raise _ArtifactSourceConflict(
                        "ARTIFACT_DOWNLOAD_SOURCE_MISSING_OR_UNSAFE"
                    )
                if final:
                    if not stat.S_ISREG(info.st_mode):
                        raise _ArtifactSourceConflict(
                            "ARTIFACT_DOWNLOAD_SOURCE_TYPE_INVALID"
                        )
                    if info.st_size != expected_size:
                        raise _ArtifactSourceConflict("ARTIFACT_DOWNLOAD_SIZE_MISMATCH")
            stable_nodes = tuple(nodes)
            _verify_descriptor_chain(stable_nodes)
            file_info = os.fstat(stable_nodes[-1].descriptor)
            descriptor_revision = build_artifact_descriptor_revision(
                output_type=output_type,
                relative_path=relative_path,
                descriptor_identity=_descriptor_identity_text(file_info),
            )
            if descriptor_revision == persisted_revision:
                _require_path_identity(stable_nodes, persisted_path_identity)
                prepared = _PreparedArtifactDownloadSource(
                    nodes=stable_nodes,
                    content=None,
                    size_bytes=expected_size,
                    media_type=media_type,
                    filename=filename,
                    content_disposition=content_disposition,
                )
                owned_descriptors.clear()
                return prepared
            if expected_size > BOUNDED_ARTIFACT_CONTENT_REVISION_MAX_BYTES:
                raise _ArtifactSourceConflict("ARTIFACT_DOWNLOAD_REVISION_MISMATCH")
            content = _read_bounded_content(stable_nodes, expected_size)
            content_revision = build_artifact_content_revision(
                output_type=output_type,
                relative_path=relative_path,
                content=content,
            )
            if content_revision != persisted_revision:
                raise _ArtifactSourceConflict("ARTIFACT_DOWNLOAD_REVISION_MISMATCH")
            _require_path_identity(stable_nodes, persisted_path_identity)
            prepared = _PreparedArtifactDownloadSource(
                nodes=stable_nodes,
                content=content,
                size_bytes=expected_size,
                media_type=media_type,
                filename=filename,
                content_disposition=content_disposition,
            )
            owned_descriptors.clear()
            return prepared
        finally:
            for descriptor in reversed(owned_descriptors):
                try:
                    os.close(descriptor)
                except OSError:
                    pass

    @staticmethod
    def _data_invalid() -> Result[ArtifactDownloadPlan]:
        return ArtifactDownloadService._failure(
            "RUN_ARTIFACT_DOWNLOAD_DATA_INVALID",
            "Persisted artifact download metadata is unavailable.",
            path="artifact_id",
        )

    @staticmethod
    def _generation_conflict() -> Result[ArtifactDownloadPlan]:
        return ArtifactDownloadService._failure(
            "RUN_ARTIFACT_DOWNLOAD_CONFLICT",
            "Artifact content is no longer available as indexed.",
            path="generation",
            reason_code="ARTIFACT_DOWNLOAD_GENERATION_CHANGED",
        )

    @staticmethod
    def _unavailable(reason_code: str) -> Result[ArtifactDownloadPlan]:
        return ArtifactDownloadService._failure(
            "RUN_ARTIFACT_DOWNLOAD_UNAVAILABLE",
            "Artifact download is temporarily unavailable.",
            path="artifact_id",
            reason_code=reason_code,
        )

    @staticmethod
    def _failure(
        code: str,
        message: str,
        *,
        path: str,
        reason_code: str | None = None,
    ) -> Result[ArtifactDownloadPlan]:
        return Result.failure(
            [
                Issue(
                    code=code,
                    message=message,
                    severity="error",
                    path=path,
                    source="artifact_download_service",
                    technical_message=None,
                    hint=None,
                    context=(
                        {"reason_code": reason_code} if reason_code is not None else {}
                    ),
                )
            ]
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
        raise ValueError("artifact path is invalid")
    path = PurePosixPath(value)
    if (
        path.as_posix() != value
        or path.is_absolute()
        or len(path.parts) < 2
        or len(path.parts) > _MAX_PATH_COMPONENTS
        or path.parts[0] != "results"
        or any(component in {"", ".", ".."} for component in value.split("/"))
    ):
        raise ValueError("artifact path is invalid")
    return value


def _validate_size(value: object) -> int:
    if (
        not isinstance(value, int)
        or isinstance(value, bool)
        or value < 0
        or value > _MAX_ARTIFACT_BYTES
    ):
        raise ValueError("artifact size is invalid")
    return value


def _validate_output_type(value: object) -> str:
    if not isinstance(value, str) or _SAFE_OUTPUT_TYPE.fullmatch(value) is None:
        raise ValueError("artifact output type is invalid")
    return value


def _validate_filename(value: object, relative_path: str) -> str:
    if (
        not isinstance(value, str)
        or not value
        or len(value) > 255
        or len(value.encode("utf-8")) > 512
        or value in {".", ".."}
        or "/" in value
        or "\\" in value
        or _CONTROL_CHARACTER.search(value) is not None
        or value != PurePosixPath(relative_path).name
    ):
        raise ValueError("artifact filename is invalid")
    return value


def _content_disposition(filename: str, artifact_id: str) -> str:
    fallback = re.sub(r"[^A-Za-z0-9._-]", "_", filename)[:160]
    if fallback in {"", ".", ".."}:
        fallback = f"{artifact_id}.bin"
    encoded = quote(filename, safe="")
    return f"attachment; filename=\"{fallback}\"; filename*=UTF-8''{encoded}"


def _fingerprint(info: os.stat_result) -> _Fingerprint:
    is_directory = stat.S_ISDIR(info.st_mode)
    return _Fingerprint(
        device=info.st_dev,
        inode=info.st_ino,
        mode=info.st_mode,
        links=0 if is_directory else info.st_nlink,
        user=0 if is_directory else info.st_uid,
        group=0 if is_directory else info.st_gid,
        # Directory content changes do not change the opened path entry's
        # identity. Keeping their mutable size/mtime would make unrelated
        # sibling creation spuriously abort a valid artifact stream.
        size=0 if is_directory else info.st_size,
        modified_ns=0 if is_directory else info.st_mtime_ns,
        changed_ns=0 if is_directory else info.st_ctime_ns,
    )


def _descriptor_identity(
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


def _descriptor_identity_text(info: os.stat_result) -> str:
    return ":".join(str(value) for value in _descriptor_identity(info))


def _path_identity(nodes: tuple[_OpenedNode, ...]) -> str:
    if len(nodes) < 3:
        raise ArtifactDownloadStreamError()
    infos = tuple(os.fstat(node.descriptor) for node in nodes[1:])
    return build_artifact_path_identity(
        parent_identities=tuple(
            (info.st_dev, info.st_ino, info.st_mode) for info in infos[:-1]
        ),
        file_identity=_descriptor_identity(infos[-1]),
    )


def _require_path_identity(
    nodes: tuple[_OpenedNode, ...], persisted_path_identity: str
) -> None:
    if _path_identity(nodes) != persisted_path_identity:
        raise _ArtifactSourceConflict("ARTIFACT_DOWNLOAD_PATH_IDENTITY_MISMATCH")


def _read_bounded_content(nodes: tuple[_OpenedNode, ...], expected_size: int) -> bytes:
    descriptor = nodes[-1].descriptor
    first = _read_bounded_pass(nodes, expected_size)
    try:
        os.lseek(descriptor, 0, os.SEEK_SET)
    except OSError:
        raise ArtifactDownloadStreamError() from None
    _verify_descriptor_chain(nodes)
    second = _read_bounded_pass(nodes, expected_size)
    if second != first:
        raise _ArtifactSourceConflict("ARTIFACT_DOWNLOAD_SOURCE_CHANGED")
    return first


def _read_bounded_pass(nodes: tuple[_OpenedNode, ...], expected_size: int) -> bytes:
    descriptor = nodes[-1].descriptor
    chunks: list[bytes] = []
    remaining = expected_size + 1
    while remaining:
        _verify_descriptor_chain(nodes)
        chunk = os.read(descriptor, min(_CHUNK_SIZE, remaining))
        _verify_descriptor_chain(nodes)
        if not chunk:
            break
        chunks.append(chunk)
        remaining -= len(chunk)
    content = b"".join(chunks)
    _verify_descriptor_chain(nodes)
    if len(content) != expected_size:
        raise _ArtifactSourceConflict("ARTIFACT_DOWNLOAD_SIZE_MISMATCH")
    return content


def _own_descriptor(owned_descriptors: list[int], descriptor: int) -> None:
    """Take cleanup ownership immediately after a successful ``os.open``."""
    try:
        owned_descriptors.append(descriptor)
    except BaseException:
        try:
            os.close(descriptor)
        except OSError:
            pass
        raise


def _verify_descriptor_chain(nodes: tuple[_OpenedNode, ...]) -> None:
    try:
        for node in nodes:
            if _fingerprint(os.fstat(node.descriptor)) != node.fingerprint:
                raise ArtifactDownloadStreamError()
            if node.parent_index is not None:
                parent = nodes[node.parent_index]
                current_entry = os.stat(
                    node.name,
                    dir_fd=parent.descriptor,
                    follow_symlinks=False,
                )
                if _fingerprint(current_entry) != node.fingerprint:
                    raise ArtifactDownloadStreamError()
    except ArtifactDownloadStreamError:
        raise
    except (OSError, TypeError, ValueError):
        raise ArtifactDownloadStreamError() from None
