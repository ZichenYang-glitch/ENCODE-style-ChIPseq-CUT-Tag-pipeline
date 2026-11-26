"""Authorized reads of one sealed track version.

The immutable object store and the version table are supplied by the archive
owner. This service checks the live session, the project and experiment, and
the native or registered source, then copies one byte window from that version.
"""

from __future__ import annotations

from collections.abc import Callable, Iterator
from contextlib import ExitStack
from datetime import datetime, timezone
from pathlib import Path
import re
from tempfile import SpooledTemporaryFile
from typing import Any, BinaryIO

from encode_pipeline.adapters.browser.scaffolds import viewer_scaffold
from encode_pipeline.persistence.archived_tracks import ArchiveRepository
from encode_pipeline.persistence.models import (
    ExperimentRevisionRow,
    ExperimentRow,
    ProjectRow,
    RegisteredResultRevisionRow,
    RegisteredResultRow,
    RunExperimentBindingRow,
    RunProjectBindingRow,
)
from encode_pipeline.platform.archived_tracks import ArchiveStoreError, CHUNK_BYTES
from encode_pipeline.platform.authentication import (
    AuthenticatedPrincipal,
    UserRole,
)
from encode_pipeline.platform.result_generations import (
    validate_artifact_generation,
    validate_artifact_revision,
)
from encode_pipeline.platform.security_audit import (
    AuditAction,
    AuditActorKind,
    AuditOutcome,
    AuditResourceKind,
    SecurityAuditEvent,
    build_audit_resource,
    new_audit_event_id,
)
from encode_pipeline.services.archive_objects import ArchiveObjects, ArchiveReader
from encode_pipeline.services.browser_bytes import (
    BrowserByteError,
    ByteRange,
    parse_byte_range,
)


_PROJECT_ID = re.compile(r"^prj_[0-9a-f]{32}$")
_EXPERIMENT_ID = re.compile(r"^exp_[0-9a-f]{32}$")
_EXPERIMENT_REVISION_ID = re.compile(r"^expr_[0-9a-f]{32}$")
_RESULT_ID = re.compile(r"^rreg_[0-9a-f]{32}$")
_RESULT_REVISION_ID = re.compile(r"^rrev_[0-9a-f]{32}$")
_ARCHIVE_VERSION_ID = re.compile(r"^aver_[0-9a-f]{32}$")
_TOKEN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,127}$")
_ARTIFACT_ID = re.compile(r"^[A-Za-z][A-Za-z0-9_.-]{0,127}$")
_LABEL = re.compile(r"^[^\x00-\x1f\x7f\\/]{1,255}$")
_STORE_CODES = frozenset(
    {
        "archive_version_not_found",
        "archive_version_conflict",
        "archive_offline",
        "archive_integrity",
        "archive_invalid",
    }
)
_RUN_FIELDS = ("run_id", "artifact_id", "generation", "revision")
_RUN_INDEX_FIELDS = ("index_artifact_id", "index_generation", "index_revision")
_REGISTERED_FIELDS = ("result_id", "result_revision_id")
_REGISTERED_INDEX_FIELDS = ("index_result_id", "index_result_revision_id")
_TRACK_TYPES = frozenset({"bigwig", "bed", "longrange"})


class ArchivedTrackError(Exception):
    """A redacted archived-track failure. The message is safe to return."""

    def __init__(self, status_code: int, code: str, message: str) -> None:
        self.status_code = status_code
        self.code = code
        self.message = message
        super().__init__(message)


class ArchivedByteBody:
    """Headers and an optional body for one content or index response."""

    def __init__(
        self,
        status_code: int,
        headers: dict[str, str],
        chunks: Iterator[bytes] | None,
        close: Callable[[], None],
    ) -> None:
        self.status_code = status_code
        self.headers = headers
        self.chunks = chunks
        self.close = close


class SqlArchiveReadGate:
    """Project, experiment, and source checks shared with existing catalog reads."""

    def __init__(self, session_factory: Any) -> None:
        self._session_factory = session_factory

    def experiment_visible(self, project_id: str, experiment_id: str) -> bool:
        with self._session_factory() as session:
            if not _project_visible(session, project_id):
                return False
            experiment = session.get(ExperimentRow, experiment_id)
            return experiment is not None and experiment.project_id == project_id

    def allows(self, version: Any) -> bool:
        project_id = _text(getattr(version, "project_id", None))
        experiment_id = _text(getattr(version, "experiment_id", None))
        revision_id = _text(getattr(version, "experiment_revision_id", None))
        try:
            source = _public_source(version.source)
        except ArchivedTrackError:
            return False
        index_fields = (
            _RUN_INDEX_FIELDS if source["kind"] == "run" else _REGISTERED_INDEX_FIELDS
        )
        if (index_fields[0] in source) != (version.track_type == "longrange"):
            return False
        with self._session_factory() as session:
            if not _project_visible(session, project_id):
                return False
            experiment = session.get(ExperimentRow, experiment_id)
            if experiment is None or experiment.project_id != project_id:
                return False
            revision = session.get(ExperimentRevisionRow, revision_id)
            if (
                revision is None
                or revision.project_id != project_id
                or revision.experiment_id != experiment_id
            ):
                return False
            if source.get("kind") == "run":
                return _run_source_bound(session, project_id, revision, source)
            if source.get("kind") == "registered":
                return _registered_source_bound(session, project_id, source)
            return False


class ArchivedTrackService:
    """List and read one sealed version without substituting another copy."""

    def __init__(
        self,
        *,
        repository: Any,
        objects: Any,
        authentication: Any,
        access: SqlArchiveReadGate,
    ) -> None:
        self._repository = repository
        self._objects = objects
        self._authentication = authentication
        self._access = access

    def list_versions(
        self,
        principal: AuthenticatedPrincipal,
        project_id: str,
        experiment_id: str,
        *,
        limit: int,
        offset: int,
    ) -> dict[str, Any]:
        self._require_reader(principal)
        project_id = _require_id(project_id, _PROJECT_ID)
        experiment_id = _require_id(experiment_id, _EXPERIMENT_ID)
        limit = _page_limit(limit)
        offset = _page_offset(offset)
        self._require_store()
        if not self._experiment_visible(project_id, experiment_id):
            raise _not_found()
        # Any offset beyond SQLite's largest row count is an empty page.
        if offset > 2**63 - 1:
            return {"items": [], "next_offset": None}
        try:
            rows, next_offset = self._repository.list_versions(
                project_id,
                experiment_id,
                limit=limit,
                offset=offset,
            )
        except Exception as exc:
            raise _from_store(exc) from None
        if not isinstance(rows, (list, tuple)):
            raise _invalid()
        if next_offset is not None and (
            isinstance(next_offset, bool)
            or not isinstance(next_offset, int)
            or next_offset < 0
        ):
            raise _invalid()
        items = []
        for row in rows:
            self._require_scope(row, project_id, experiment_id)
            if not self._allows(row):
                continue
            items.append(_public_descriptor(row))
        return {"items": items, "next_offset": next_offset}

    def get_version(
        self,
        principal: AuthenticatedPrincipal,
        project_id: str,
        experiment_id: str,
        archive_version_id: str,
    ) -> dict[str, Any]:
        version = self._authorized_version(
            principal,
            project_id,
            experiment_id,
            archive_version_id,
        )
        return _public_descriptor(version)

    def open_bytes(
        self,
        principal: AuthenticatedPrincipal,
        project_id: str,
        experiment_id: str,
        archive_version_id: str,
        *,
        index: bool,
        range_header: str | None,
        head: bool,
    ) -> ArchivedByteBody:
        version = self._authorized_version(
            principal,
            project_id,
            experiment_id,
            archive_version_id,
        )
        member = _member(version, index=index)
        if self._objects is None:
            raise _offline()
        try:
            with ExitStack() as response_files:
                with ExitStack() as sources:
                    reader = sources.enter_context(self._objects.open_member(member))
                    companion = None
                    if "index" in version.members:
                        other = version.members["data" if index else "index"]
                        companion = sources.enter_context(
                            self._objects.open_member(other)
                        )
                    size = _member_size(member)
                    try:
                        selected = parse_byte_range(range_header, size)
                    except BrowserByteError:
                        return _range_error(size)
                    status, headers = _byte_headers(size, selected)
                    if head:
                        return ArchivedByteBody(status, headers, None, lambda: None)

                    # Keep small windows in memory. For a large range or full GET,
                    # spill only the requested bytes to a private temporary file.
                    body = response_files.enter_context(
                        SpooledTemporaryFile(max_size=CHUNK_BYTES)
                    )
                    _copy_window(reader, body, selected, size)
                    if companion is not None:
                        # The empty read rechecks the pinned companion's identity.
                        companion.read(0, 0)
                self._audit_registered_bytes(principal, version, index=index)
                chunks = iter(lambda: body.read(CHUNK_BYTES), b"")
                return ArchivedByteBody(
                    status, headers, chunks, response_files.pop_all().close
                )
        except Exception as exc:
            raise _from_store(exc) from None

    def _authorized_version(
        self,
        principal: AuthenticatedPrincipal,
        project_id: str,
        experiment_id: str,
        archive_version_id: str,
    ) -> Any:
        self._require_reader(principal)
        project_id = _require_id(project_id, _PROJECT_ID)
        experiment_id = _require_id(experiment_id, _EXPERIMENT_ID)
        archive_version_id = _require_id(archive_version_id, _ARCHIVE_VERSION_ID)
        self._require_store()
        if not self._experiment_visible(project_id, experiment_id):
            raise _not_found()
        try:
            version = self._repository.get_version(
                project_id,
                experiment_id,
                archive_version_id,
            )
        except Exception as exc:
            raise _from_store(exc) from None
        if version is None:
            raise _not_found()
        self._require_scope(version, project_id, experiment_id)
        if _text(getattr(version, "archive_version_id", None)) != archive_version_id:
            raise _not_found()
        if not self._allows(version):
            raise _not_found()
        return version

    def _audit_registered_bytes(
        self,
        principal: AuthenticatedPrincipal,
        version: Any,
        *,
        index: bool,
    ) -> None:
        """Record the same success audit as a registered-result download.

        Metadata, HEAD, and native run bytes are not downloads of a registered
        revision. A failed audit write does not start the response.
        """

        source = getattr(version, "source", None)
        if not isinstance(source, dict) or source.get("kind") != "registered":
            return
        result_id = source.get("index_result_id" if index else "result_id")
        revision_id = source.get(
            "index_result_revision_id" if index else "result_revision_id"
        )
        try:
            self._authentication.record_security_audit(
                SecurityAuditEvent(
                    event_id=new_audit_event_id(),
                    occurred_at=datetime.now(timezone.utc),
                    action=AuditAction.REGISTERED_RESULT_DOWNLOAD,
                    outcome=AuditOutcome.SUCCEEDED,
                    actor_kind=AuditActorKind.USER,
                    actor_user_id=principal.user_id,
                    resource=build_audit_resource(
                        AuditResourceKind.REGISTERED_RESULT,
                        getattr(version, "project_id", None),
                        result_id,
                        revision_id,
                    ),
                )
            )
        except Exception:
            raise _unavailable() from None

    def _require_reader(self, principal: AuthenticatedPrincipal) -> None:
        if self._authentication is None or not isinstance(
            principal, AuthenticatedPrincipal
        ):
            raise _unauthenticated()
        if principal.role not in {UserRole.ADMINISTRATOR, UserRole.MEMBER}:
            raise _unauthenticated()
        try:
            account = self._authentication.get_account_by_id(principal.user_id)
        except KeyError:
            raise _unauthenticated() from None
        except Exception:
            raise _unauthenticated() from None
        if (
            account is None
            or not getattr(account, "enabled", False)
            or getattr(account, "user_id", None) != principal.user_id
            or getattr(account, "role", None)
            not in {UserRole.ADMINISTRATOR, UserRole.MEMBER}
        ):
            raise _unauthenticated()

    def _require_store(self) -> None:
        if self._repository is None:
            raise _offline()

    def _allows(self, version: Any) -> bool:
        try:
            return self._access.allows(version)
        except Exception:
            raise _offline() from None

    def _experiment_visible(self, project_id: str, experiment_id: str) -> bool:
        try:
            return self._access.experiment_visible(project_id, experiment_id)
        except Exception:
            raise _offline() from None

    def _require_scope(self, version: Any, project_id: str, experiment_id: str) -> None:
        if _text(getattr(version, "project_id", None)) != project_id:
            raise _not_found()
        if _text(getattr(version, "experiment_id", None)) != experiment_id:
            raise _not_found()


def build_archived_track_service(
    *,
    session_factory: Any,
    authentication: Any,
    archive_root: Path | None,
) -> ArchivedTrackService:
    """Use the configured archive volume and the same catalog as experiment reads."""
    return ArchivedTrackService(
        repository=ArchiveRepository(session_factory),
        objects=ArchiveObjects(archive_root),
        authentication=authentication,
        access=SqlArchiveReadGate(session_factory),
    )


def _project_visible(session: Any, project_id: str) -> bool:
    project = session.get(ProjectRow, project_id)
    return (
        project is not None and project.archived_at is None and project.kind == "user"
    )


def _run_source_bound(
    session: Any,
    project_id: str,
    revision: ExperimentRevisionRow,
    source: dict[str, Any],
) -> bool:
    # The index pin names another artifact of this same run, not a second run.
    run_id = source.get("run_id")
    if not isinstance(run_id, str):
        return False
    binding = session.get(RunExperimentBindingRow, run_id)
    if (
        binding is None
        or binding.experiment_revision_id != revision.experiment_revision_id
        or binding.payload_digest != revision.payload_digest
    ):
        return False
    project_binding = session.get(RunProjectBindingRow, run_id)
    return project_binding is not None and project_binding.project_id == project_id


def _registered_source_bound(
    session: Any,
    project_id: str,
    source: dict[str, Any],
) -> bool:
    if not _revision_in_project(
        session,
        project_id,
        source.get("result_id"),
        source.get("result_revision_id"),
    ):
        return False
    if _has_any(source, _REGISTERED_INDEX_FIELDS):
        return _revision_in_project(
            session,
            project_id,
            source.get("index_result_id"),
            source.get("index_result_revision_id"),
        )
    return True


def _revision_in_project(
    session: Any,
    project_id: str,
    result_id: object,
    result_revision_id: object,
) -> bool:
    if not isinstance(result_id, str) or not isinstance(result_revision_id, str):
        return False
    result = session.get(RegisteredResultRow, result_id)
    if result is None or result.project_id != project_id:
        return False
    revision = session.get(RegisteredResultRevisionRow, result_revision_id)
    return revision is not None and revision.result_id == result_id


def _has_any(source: dict[str, Any], names: tuple[str, ...]) -> bool:
    return any(name in source and source[name] not in (None, "") for name in names)


def _public_descriptor(version: Any) -> dict[str, Any]:
    track_type = _text(getattr(version, "track_type", None))
    if track_type not in _TRACK_TYPES:
        raise _invalid()
    assembly = viewer_scaffold(_text(getattr(version, "assembly", None)))
    if assembly is None:
        raise _invalid()
    version_number = getattr(version, "version_number", None)
    if (
        isinstance(version_number, bool)
        or not isinstance(version_number, int)
        or version_number < 1
    ):
        raise _invalid()
    archive_version_id = _require_id(
        getattr(version, "archive_version_id", None),
        _ARCHIVE_VERSION_ID,
    )
    source = _public_source(getattr(version, "source", None))
    members = getattr(version, "members", None)
    return {
        "archive_version_id": archive_version_id,
        "project_id": _require_id(getattr(version, "project_id", None), _PROJECT_ID),
        "experiment_id": _require_id(
            getattr(version, "experiment_id", None),
            _EXPERIMENT_ID,
        ),
        "experiment_revision_id": _require_id(
            getattr(version, "experiment_revision_id", None),
            _EXPERIMENT_REVISION_ID,
        ),
        "track_key": _label(getattr(version, "track_key", None)),
        "version_number": version_number,
        "display_name": _label(getattr(version, "display_name", None)),
        "assembly": assembly,
        "track_type": track_type,
        "source": source,
        "has_index": _indexed(version, members),
    }


def _public_source(source: object) -> dict[str, str]:
    if not isinstance(source, dict):
        raise _invalid()
    kind = source.get("kind")
    if kind == "run":
        fields = _RUN_FIELDS
        index_fields = _RUN_INDEX_FIELDS
        patterns = {
            "run_id": _TOKEN,
            "artifact_id": _ARTIFACT_ID,
            "generation": _TOKEN,
            "revision": _TOKEN,
            "index_artifact_id": _ARTIFACT_ID,
            "index_generation": _TOKEN,
            "index_revision": _TOKEN,
        }
    elif kind == "registered":
        fields = _REGISTERED_FIELDS
        index_fields = _REGISTERED_INDEX_FIELDS
        patterns = {
            "result_id": _RESULT_ID,
            "result_revision_id": _RESULT_REVISION_ID,
            "index_result_id": _RESULT_ID,
            "index_result_revision_id": _RESULT_REVISION_ID,
        }
    else:
        raise _invalid()
    public = {"kind": kind}
    for name in fields:
        public[name] = _require_id(source.get(name), patterns[name])
    if _has_any(source, index_fields):
        for name in index_fields:
            public[name] = _require_id(source.get(name), patterns[name])
    if kind == "run":
        try:
            validate_artifact_generation(public["generation"])
            validate_artifact_revision(public["revision"])
            if "index_artifact_id" in public:
                validate_artifact_generation(public["index_generation"])
                validate_artifact_revision(public["index_revision"])
        except ValueError:
            raise _invalid() from None
    return public


def _indexed(version: Any, members: object) -> bool:
    """List rows omit members; detail rows carry the sealed data and index."""

    if not isinstance(members, dict):
        raise _invalid()
    if not members:
        declared = getattr(version, "has_index", None)
        if not isinstance(declared, bool):
            raise _invalid()
        return declared
    if "data" not in members:
        raise _invalid()
    indexed = "index" in members
    declared = getattr(version, "has_index", None)
    if isinstance(declared, bool) and declared != indexed:
        raise _integrity()
    return indexed


def _member(version: Any, *, index: bool) -> Any:
    members = getattr(version, "members", None)
    if not isinstance(members, dict) or "data" not in members:
        raise _invalid()
    has_index = _indexed(version, members)
    if index:
        if not has_index:
            raise _not_found()
        return members["index"]
    return members["data"]


def _member_size(member: Any) -> int:
    size = getattr(member, "size_bytes", None)
    if isinstance(size, bool) or not isinstance(size, int) or size < 0:
        raise _invalid()
    return size


def _copy_window(
    reader: ArchiveReader, body: BinaryIO, selected: ByteRange | None, size: int
) -> None:
    start, end = (0, size) if selected is None else (selected.start, selected.end + 1)
    while start < end:
        remaining_in_chunk = CHUNK_BYTES - start % CHUNK_BYTES
        stop = min(end, start + remaining_in_chunk)
        block = reader.read(start, stop)
        if not isinstance(block, (bytes, bytearray)) or len(block) != stop - start:
            raise _integrity()
        if body.write(block) != len(block):
            raise _unavailable()
        start = stop
    body.seek(0)


def _byte_headers(size: int, selected: ByteRange | None) -> tuple[int, dict[str, str]]:
    if selected is None:
        length = size
        status = 200
        headers = {
            "Accept-Ranges": "bytes",
            "Content-Length": str(length),
            "Cache-Control": "private, no-store",
            "X-Content-Type-Options": "nosniff",
        }
        return status, headers
    return 206, {
        "Accept-Ranges": "bytes",
        "Content-Length": str(selected.length),
        "Content-Range": f"bytes {selected.start}-{selected.end}/{size}",
        "Cache-Control": "private, no-store",
        "X-Content-Type-Options": "nosniff",
    }


def _range_error(size: int) -> ArchivedByteBody:
    return ArchivedByteBody(
        416,
        {
            "Accept-Ranges": "bytes",
            "Content-Range": f"bytes */{size}",
            "Content-Length": "0",
            "Cache-Control": "private, no-store",
            "X-Content-Type-Options": "nosniff",
        },
        None,
        lambda: None,
    )


def _from_store(exc: Exception) -> ArchivedTrackError:
    if isinstance(exc, ArchivedTrackError):
        return exc
    code = getattr(exc, "code", None)
    if isinstance(exc, ArchiveStoreError) and code in _STORE_CODES:
        return _mapped(str(code))
    if isinstance(exc, BrowserByteError):
        return _invalid()
    return _unavailable()


def _mapped(code: str) -> ArchivedTrackError:
    status, message = {
        "archive_version_not_found": (404, "The archived track was not found."),
        "archive_version_conflict": (409, "The archived track changed."),
        "archive_offline": (409, "The archived track is unavailable."),
        "archive_integrity": (409, "The archived track could not be verified."),
        "archive_invalid": (400, "The archived track request is invalid."),
    }[code]
    return ArchivedTrackError(status, code, message)


def _not_found() -> ArchivedTrackError:
    return _mapped("archive_version_not_found")


def _offline() -> ArchivedTrackError:
    return _mapped("archive_offline")


def _integrity() -> ArchivedTrackError:
    return _mapped("archive_integrity")


def _invalid() -> ArchivedTrackError:
    return _mapped("archive_invalid")


def _unavailable() -> ArchivedTrackError:
    return ArchivedTrackError(
        500,
        "INTERNAL_SERVER_ERROR",
        "The archived track is temporarily unavailable.",
    )


def _unauthenticated() -> ArchivedTrackError:
    return ArchivedTrackError(
        401,
        "AUTHENTICATION_REQUIRED",
        "Authentication is required.",
    )


def _require_id(value: object, pattern: re.Pattern[str]) -> str:
    if not isinstance(value, str) or pattern.fullmatch(value) is None:
        raise _invalid()
    return value


def _label(value: object) -> str:
    if not isinstance(value, str) or _LABEL.fullmatch(value) is None:
        raise _invalid()
    return value


def _text(value: object) -> str:
    return value if isinstance(value, str) else ""


def _page_limit(value: object) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or not 1 <= value <= 100:
        raise _invalid()
    return value


def _page_offset(value: object) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise _invalid()
    return value
