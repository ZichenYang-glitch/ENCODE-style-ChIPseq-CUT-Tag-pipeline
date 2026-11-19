"""Register external files and serve one exact revision from a private copy."""

from __future__ import annotations

from collections.abc import Callable, Mapping
from contextlib import suppress
from dataclasses import dataclass, field
from datetime import datetime, timezone
import base64
import json
import os
from pathlib import Path
import stat
import threading
from typing import Any
from uuid import uuid4

from encode_pipeline.persistence.registered_results import (
    RegisteredResultStoreError,
    RevisionDraft,
    SqlAlchemyRegisteredResultRepository,
)
from encode_pipeline.platform.authentication import AuthenticatedPrincipal
from encode_pipeline.platform.input_registry import validate_input_file_relative_path
from encode_pipeline.platform.registered_results import (
    ParentEdge,
    RegisteredResultPayloadError,
    RevisionRecord,
    build_revision_digest,
    metadata_from_canonical_json,
    parse_metadata,
    public_revision,
    request_sha256,
    validate_idempotency_key,
    validate_project_id,
    validate_result_id,
    validate_result_revision_id,
    validate_stable_key,
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
from encode_pipeline.services.input_file_access import (
    FileObservation,
    InputFileAccess,
    InputFileAccessError,
)
from encode_pipeline.services.private_storage_pools import (
    PrivateStoragePoolConfig,
    PrivateStoragePoolConfigError,
)
from encode_pipeline.services.registered_result_bytes import (
    RegisteredResultBytesError,
    discard_private_copy,
    seal_private_copy,
)


_MAX_PAGE = 50
_CURSOR_PREFIX = "rcur1."


def open_verified_copy(copy_path: Path) -> int:
    """Open one sealed response file for the caller to send and close."""

    return os.open(copy_path, os.O_RDONLY | os.O_NOFOLLOW | os.O_CLOEXEC)


class RegisteredResultAPIError(Exception):
    """Closed HTTP failure that carries no path, pool, or source text."""

    def __init__(self, status_code: int, code: str, message: str) -> None:
        super().__init__(code)
        self.status_code = status_code
        self.code = code
        self.message = message


class _ReleaseOnce:
    """Claim the single close of one response copy."""

    __slots__ = ("_done", "_lock")

    def __init__(self) -> None:
        self._lock = threading.Lock()
        self._done = False

    def claim(self) -> bool:
        with self._lock:
            if self._done:
                return False
            self._done = True
            return True


@dataclass(frozen=True)
class RegisteredResultDownload:
    """One verified response copy. The response owns the open descriptor."""

    descriptor: int
    copy_path: Path
    directory: Path
    size_bytes: int
    content_sha256: str
    filename: str
    _release: _ReleaseOnce = field(
        default_factory=_ReleaseOnce,
        init=False,
        repr=False,
        compare=False,
    )

    def release(self) -> None:
        """Close the descriptor and remove the private copy once."""

        if not self._release.claim():
            return
        with suppress(OSError):
            os.close(self.descriptor)
        discard_private_copy(self.copy_path, self.directory)


class RegisteredResultService:
    """Catalog writes and exact-revision reads for files that are not runs."""

    def __init__(
        self,
        repository: SqlAlchemyRegisteredResultRepository,
        *,
        input_registry: Any,
        reference_profiles: Any,
        authentication: Any,
        storage_pool_config_path: Path | None,
        response_root: Path,
    ) -> None:
        self._repository = repository
        self._input_registry = input_registry
        self._reference_profiles = reference_profiles
        self._authentication = authentication
        self._storage_pool_config_path = storage_pool_config_path
        self._response_root = response_root
        self.after_copy_verified: Callable[[Path], None] | None = None

    def list_current(
        self,
        principal: AuthenticatedPrincipal,
        project_id: str,
        *,
        cursor: str | None,
        limit: int,
        include_archived: bool,
    ) -> tuple[list[dict[str, Any]], str | None]:
        self._require_enabled(principal)
        project_id = self._project(project_id)
        limit = self._limit(limit)
        after = self._results_cursor(
            project_id,
            cursor,
            include_archived=include_archived,
        )
        try:
            rows = self._repository.list_current(
                project_id,
                include_archived=include_archived,
                after=after,
                limit=limit + 1,
            )
        except RegisteredResultStoreError as exc:
            raise self._from_store(exc) from None
        page = rows[:limit]
        next_cursor = None
        if len(rows) > limit and page:
            last = page[-1]
            next_cursor = _encode_cursor(
                {
                    "include_archived": include_archived,
                    "kind": "results",
                    "project_id": project_id,
                    "result_id": last.result_id,
                    "stable_key": last.stable_key,
                }
            )
        return [public_revision(row) for row in page], next_cursor

    def list_revisions(
        self,
        principal: AuthenticatedPrincipal,
        project_id: str,
        result_id: str,
        *,
        cursor: str | None,
        limit: int,
    ) -> tuple[list[dict[str, Any]], str | None]:
        self._require_enabled(principal)
        project_id = self._project(project_id)
        result_id = self._result(result_id)
        limit = self._limit(limit)
        after_number = self._revisions_cursor(project_id, result_id, cursor)
        try:
            rows = self._repository.list_revisions(
                project_id,
                result_id,
                after_revision_number=after_number,
                limit=limit + 1,
            )
        except RegisteredResultStoreError as exc:
            raise self._from_store(exc) from None
        page = rows[:limit]
        next_cursor = None
        if len(rows) > limit and page:
            last = page[-1]
            next_cursor = _encode_cursor(
                {
                    "kind": "revisions",
                    "project_id": project_id,
                    "result_id": result_id,
                    "revision_number": last.revision_number,
                }
            )
        return [public_revision(row) for row in page], next_cursor

    def get_revision(
        self,
        principal: AuthenticatedPrincipal,
        project_id: str,
        result_revision_id: str,
    ) -> dict[str, Any]:
        self._require_enabled(principal)
        record = self._visible_revision(project_id, result_revision_id)
        return public_revision(record)

    def prepare_download(
        self,
        principal: AuthenticatedPrincipal,
        project_id: str,
        result_revision_id: str,
    ) -> RegisteredResultDownload:
        self._require_enabled(principal)
        record = self._visible_revision(project_id, result_revision_id)
        metadata = self._stored_metadata(record)
        self._require_digest(record, metadata)
        self._ensure_response_root()
        try:
            copy_path, directory = seal_private_copy(
                root=self._pool_root(record.storage_pool_id),
                relative_path=record.relative_path,
                expected_size=record.size_bytes,
                expected_sha256=record.content_sha256,
                response_root=self._response_root,
            )
        except RegisteredResultBytesError as exc:
            raise self._bytes_unavailable() from exc
        try:
            if self.after_copy_verified is not None:
                self.after_copy_verified(copy_path)
        except Exception:
            discard_private_copy(copy_path, directory)
            raise
        descriptor: int | None = None
        try:
            descriptor = open_verified_copy(copy_path)
            info = os.fstat(descriptor)
            if not stat.S_ISREG(info.st_mode) or info.st_size != record.size_bytes:
                raise OSError
        except OSError:
            if descriptor is not None:
                with suppress(OSError):
                    os.close(descriptor)
            discard_private_copy(copy_path, directory)
            raise self._bytes_unavailable() from None
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
                        record.project_id,
                        record.result_id,
                        record.result_revision_id,
                    ),
                )
            )
        except Exception:
            with suppress(OSError):
                os.close(descriptor)
            discard_private_copy(copy_path, directory)
            raise
        return RegisteredResultDownload(
            descriptor=descriptor,
            copy_path=copy_path,
            directory=directory,
            size_bytes=record.size_bytes,
            content_sha256=record.content_sha256,
            filename=f"registered-result.{metadata.format}",
        )

    def register(
        self,
        principal: AuthenticatedPrincipal,
        project_id: str,
        *,
        idempotency_key: str,
        stable_key: str,
        metadata_payload: Mapping[str, Any],
        pool_relative_path: str,
    ) -> dict[str, Any]:
        self._require_enabled(principal)
        project_id = self._project(project_id)
        idempotency_key = self._idempotency(idempotency_key)
        metadata, relative_path, stable = self._registration_fields(
            metadata_payload,
            pool_relative_path,
            stable_key,
        )
        self._check_reference(metadata.reference)
        body = {
            "display_name": metadata.display_name,
            "format": metadata.format,
            "method": metadata.method,
            "origin": metadata.origin,
            "pool_relative_path": relative_path,
            "reference": metadata.reference,
            "stable_key": stable,
        }
        replay = self._matching_receipt(
            principal.user_id,
            idempotency_key,
            request_sha256("register", project_id, None, body),
        )
        if replay is not None:
            return public_revision(replay)
        binding = self._active_binding(project_id)
        observation = self._observe(binding.storage_pool_id, relative_path)
        parents = self._parent_edges(project_id, metadata)
        draft = self._draft(
            project_id=project_id,
            result_id=f"rreg_{uuid4().hex}",
            revision_number=1,
            metadata=metadata,
            observation=observation,
            storage_pool_id=binding.storage_pool_id,
            relative_path=relative_path,
            parents=parents,
        )
        try:
            record, _replayed = self._repository.register(
                actor_user_id=principal.user_id,
                idempotency_key=idempotency_key,
                request_sha256=request_sha256("register", project_id, None, body),
                project_id=project_id,
                stable_key=stable,
                draft=draft,
            )
        except RegisteredResultStoreError as exc:
            raise self._from_store(exc) from None
        return public_revision(record)

    def revise(
        self,
        principal: AuthenticatedPrincipal,
        project_id: str,
        result_id: str,
        *,
        idempotency_key: str,
        expected_revision_number: int,
        display_name: str | None,
        result_format: str | None,
        origin: Mapping[str, Any] | None,
        method: Mapping[str, Any] | None,
        reference: Mapping[str, Any] | None,
        pool_relative_path: str | None,
    ) -> dict[str, Any]:
        self._require_enabled(principal)
        project_id = self._project(project_id)
        result_id = self._result(result_id)
        idempotency_key = self._idempotency(idempotency_key)
        expected = self._revision_number(expected_revision_number)
        normalized = self._revise_request(
            display_name=display_name,
            result_format=result_format,
            origin=origin,
            method=method,
            reference=reference,
            pool_relative_path=pool_relative_path,
        )
        body = {"expected_revision_number": expected, **normalized}
        replay = self._matching_receipt(
            principal.user_id,
            idempotency_key,
            request_sha256("revise", project_id, result_id, body),
        )
        if replay is not None:
            return public_revision(replay)
        current = self._writable_head(project_id, result_id, expected)
        metadata = self._merged_metadata(current, normalized)
        self._check_reference(metadata.reference)
        if normalized["pool_relative_path"] is None:
            observation = self._observe(current.storage_pool_id, current.relative_path)
            if (
                observation.size_bytes != current.size_bytes
                or observation.content_sha256 != current.content_sha256
            ):
                raise self._bytes_unavailable()
            storage_pool_id = current.storage_pool_id
            relative_path = current.relative_path
        else:
            binding = self._active_binding(project_id)
            relative_path = normalized["pool_relative_path"]
            observation = self._observe(binding.storage_pool_id, relative_path)
            storage_pool_id = binding.storage_pool_id
        parents = self._parent_edges(project_id, metadata)
        draft = self._draft(
            project_id=project_id,
            result_id=result_id,
            revision_number=current.revision_number + 1,
            metadata=metadata,
            observation=observation,
            storage_pool_id=storage_pool_id,
            relative_path=relative_path,
            parents=parents,
        )
        try:
            record, _replayed = self._repository.revise(
                actor_user_id=principal.user_id,
                idempotency_key=idempotency_key,
                request_sha256=request_sha256("revise", project_id, result_id, body),
                project_id=project_id,
                expected_revision_number=expected,
                draft=draft,
            )
        except RegisteredResultStoreError as exc:
            raise self._from_store(exc) from None
        return public_revision(record)

    def archive(
        self,
        principal: AuthenticatedPrincipal,
        project_id: str,
        result_id: str,
        *,
        idempotency_key: str,
        expected_revision_number: int,
    ) -> dict[str, Any]:
        self._require_enabled(principal)
        project_id = self._project(project_id)
        result_id = self._result(result_id)
        idempotency_key = self._idempotency(idempotency_key)
        expected = self._revision_number(expected_revision_number)
        body = {"expected_revision_number": expected}
        request_hash = request_sha256("archive", project_id, result_id, body)
        replay = self._matching_receipt(
            principal.user_id,
            idempotency_key,
            request_hash,
        )
        if replay is not None:
            return self._archive_body(replay)
        archived_at = datetime.now(timezone.utc)
        try:
            record, _replayed = self._repository.archive(
                actor_user_id=principal.user_id,
                idempotency_key=idempotency_key,
                request_sha256=request_hash,
                project_id=project_id,
                result_id=result_id,
                expected_revision_number=expected,
                archived_at=archived_at,
            )
        except RegisteredResultStoreError as exc:
            raise self._from_store(exc) from None
        return self._archive_body(record)

    def _archive_body(self, record: RevisionRecord) -> dict[str, Any]:
        if record.result_archived_at is None:
            raise RegisteredResultAPIError(
                409,
                "catalog_revision_conflict",
                "The registered result changed.",
            )
        return {
            "archived_at": _zulu(record.result_archived_at),
            "result_id": record.result_id,
        }

    def _registration_fields(
        self,
        metadata_payload: Mapping[str, Any],
        pool_relative_path: str,
        stable_key: str,
    ):
        try:
            metadata = parse_metadata(metadata_payload)
            relative_path = validate_input_file_relative_path(pool_relative_path)
            stable = validate_stable_key(stable_key)
        except (RegisteredResultPayloadError, TypeError, ValueError):
            raise self._invalid() from None
        return metadata, relative_path, stable

    def _revise_request(
        self,
        *,
        display_name: str | None,
        result_format: str | None,
        origin: Mapping[str, Any] | None,
        method: Mapping[str, Any] | None,
        reference: Mapping[str, Any] | None,
        pool_relative_path: str | None,
    ) -> dict[str, Any]:
        try:
            return {
                "display_name": _normalized_part(
                    "display_name",
                    display_name,
                ),
                "format": _normalized_part("format", result_format),
                "method": _normalized_part("method", method),
                "origin": _normalized_part("origin", origin),
                "pool_relative_path": (
                    None
                    if pool_relative_path is None
                    else validate_input_file_relative_path(pool_relative_path)
                ),
                "reference": _normalized_part("reference", reference),
            }
        except (RegisteredResultPayloadError, TypeError, ValueError):
            raise self._invalid() from None

    def _merged_metadata(self, current: RevisionRecord, normalized: Mapping[str, Any]):
        current_metadata = self._stored_metadata(current)
        payload = {
            "display_name": _chosen(
                normalized["display_name"],
                current_metadata.display_name,
            ),
            "format": _chosen(normalized["format"], current_metadata.format),
            "method": _chosen(normalized["method"], current_metadata.method),
            "origin": _chosen(normalized["origin"], current_metadata.origin),
            "reference": _chosen(
                normalized["reference"],
                current_metadata.reference,
            ),
        }
        try:
            return parse_metadata(payload)
        except RegisteredResultPayloadError:
            raise self._invalid() from None

    def _draft(
        self,
        *,
        project_id: str,
        result_id: str,
        revision_number: int,
        metadata,
        observation: FileObservation,
        storage_pool_id: str,
        relative_path: str,
        parents: tuple[ParentEdge, ...],
    ) -> RevisionDraft:
        revision_id = f"rrev_{uuid4().hex}"
        canonical = metadata.canonical_json()
        return RevisionDraft(
            result_id=result_id,
            result_revision_id=revision_id,
            revision_number=revision_number,
            canonical_payload=canonical,
            revision_digest=build_revision_digest(
                project_id=project_id,
                result_id=result_id,
                result_revision_id=revision_id,
                revision_number=revision_number,
                canonical_payload=canonical,
                storage_pool_id=storage_pool_id,
                relative_path=relative_path,
                size_bytes=observation.size_bytes,
                content_sha256=observation.content_sha256,
                parents=parents,
            ),
            created_at=datetime.now(timezone.utc),
            storage_pool_id=storage_pool_id,
            relative_path=relative_path,
            size_bytes=observation.size_bytes,
            content_sha256=observation.content_sha256,
            parents=parents,
        )

    def _parent_edges(self, project_id: str, metadata) -> tuple[ParentEdge, ...]:
        edges: list[ParentEdge] = []
        for revision_id, role in metadata.parent_refs:
            parent = self._repository.get_revision(revision_id)
            if (
                parent is None
                or parent.project_id != project_id
                or parent.project_archived_at is not None
            ):
                raise RegisteredResultAPIError(
                    400,
                    "catalog_payload_invalid",
                    "Parent result revision was not found.",
                )
            edges.append(
                ParentEdge(
                    result_revision_id=revision_id,
                    role=role,
                    revision_digest=parent.revision_digest,
                )
            )
        return tuple(edges)

    def _check_reference(self, reference: Mapping[str, Any]) -> None:
        if reference.get("kind") != "reference_profile_revision":
            return
        revision_id = reference["reference_profile_revision_id"]
        try:
            stored = self._reference_profiles.get_revision(revision_id)
        except KeyError:
            raise self._invalid() from None
        if stored.public_identity_sha256 != reference["public_identity_sha256"]:
            raise self._invalid()

    def _writable_head(
        self,
        project_id: str,
        result_id: str,
        expected_revision_number: int,
    ) -> RevisionRecord:
        current = self._repository.get_head(result_id)
        if (
            current is None
            or current.project_id != project_id
            or current.project_archived_at is not None
        ):
            raise self._missing()
        if current.result_archived_at is not None:
            raise RegisteredResultAPIError(
                409,
                "catalog_already_archived",
                "The registered result is archived.",
            )
        if current.revision_number != expected_revision_number:
            raise RegisteredResultAPIError(
                409,
                "catalog_revision_conflict",
                "The registered result changed.",
            )
        return current

    def _visible_revision(
        self,
        project_id: str,
        result_revision_id: str,
    ) -> RevisionRecord:
        project_id = self._project(project_id)
        try:
            result_revision_id = validate_result_revision_id(result_revision_id)
        except RegisteredResultPayloadError:
            raise self._missing() from None
        record = self._repository.get_revision(result_revision_id)
        if (
            record is None
            or record.project_id != project_id
            or record.project_archived_at is not None
        ):
            raise self._missing()
        return record

    def _stored_metadata(self, record: RevisionRecord):
        try:
            return metadata_from_canonical_json(record.canonical_payload)
        except RegisteredResultPayloadError:
            raise self._bytes_unavailable() from None

    def _require_digest(self, record: RevisionRecord, metadata) -> None:
        parents = tuple(
            ParentEdge(
                result_revision_id=parent.result_revision_id,
                role=parent.role,
                revision_digest=parent.revision_digest,
            )
            for parent in record.parents
        )
        try:
            digest = build_revision_digest(
                project_id=record.project_id,
                result_id=record.result_id,
                result_revision_id=record.result_revision_id,
                revision_number=record.revision_number,
                canonical_payload=record.canonical_payload,
                storage_pool_id=record.storage_pool_id,
                relative_path=record.relative_path,
                size_bytes=record.size_bytes,
                content_sha256=record.content_sha256,
                parents=parents,
            )
        except RegisteredResultPayloadError:
            raise self._bytes_unavailable() from None
        if digest != record.revision_digest or len(parents) != len(
            metadata.parent_refs
        ):
            raise self._bytes_unavailable()

    def _observe(self, storage_pool_id: str, relative_path: str) -> FileObservation:
        try:
            access = InputFileAccess.from_config(
                self._pool_config(),
                self._pool_config_key(storage_pool_id),
            )
            return access.observe(relative_path)
        except (InputFileAccessError, PrivateStoragePoolConfigError, OSError):
            raise self._bytes_unavailable() from None

    def _pool_root(self, storage_pool_id: str) -> Path:
        try:
            return self._pool_config().root_for(self._pool_config_key(storage_pool_id))
        except (PrivateStoragePoolConfigError, KeyError, ValueError, OSError):
            raise self._bytes_unavailable() from None

    def _pool_config(self) -> PrivateStoragePoolConfig:
        if self._storage_pool_config_path is None:
            raise self._bytes_unavailable()
        try:
            return PrivateStoragePoolConfig.from_file(self._storage_pool_config_path)
        except (PrivateStoragePoolConfigError, OSError):
            raise self._bytes_unavailable() from None

    def _pool_config_key(self, storage_pool_id: str) -> str:
        try:
            pool = self._input_registry.get_storage_pool(storage_pool_id)
        except KeyError:
            raise self._bytes_unavailable() from None
        return pool.config_key

    def _active_binding(self, project_id: str):
        try:
            return self._input_registry.get_active_project_storage_pool_binding(
                project_id
            )
        except KeyError as exc:
            if str(exc).startswith("Unknown Project StoragePool binding"):
                raise self._invalid() from None
            raise self._missing() from None
        except ValueError:
            raise self._missing() from None

    def _matching_receipt(
        self,
        actor_user_id: str,
        idempotency_key: str,
        request_hash: str,
    ) -> RevisionRecord | None:
        try:
            return self._repository.lookup_receipt(
                actor_user_id,
                idempotency_key,
                request_hash,
            )
        except RegisteredResultStoreError as exc:
            raise self._from_store(exc) from None

    def _ensure_response_root(self) -> None:
        try:
            self._response_root.mkdir(parents=True, exist_ok=True)
            os.chmod(self._response_root, 0o700)
        except OSError:
            raise self._bytes_unavailable() from None

    def _require_enabled(self, principal: AuthenticatedPrincipal) -> None:
        try:
            account = self._authentication.get_account_by_id(principal.user_id)
        except KeyError:
            raise self._unauthenticated() from None
        if not account.enabled:
            raise self._unauthenticated()

    def _results_cursor(
        self,
        project_id: str,
        cursor: str | None,
        *,
        include_archived: bool,
    ) -> tuple[str, str] | None:
        if cursor is None:
            return None
        payload = _decode_cursor(cursor)
        if (
            payload.get("kind") != "results"
            or payload.get("project_id") != project_id
            or payload.get("include_archived") is not include_archived
            or set(payload)
            != {
                "include_archived",
                "kind",
                "project_id",
                "result_id",
                "stable_key",
            }
        ):
            raise RegisteredResultAPIError(
                400,
                "catalog_cursor_invalid",
                "The cursor is invalid.",
            )
        result_id = payload["result_id"]
        stable_key = payload["stable_key"]
        if not isinstance(result_id, str) or not isinstance(stable_key, str):
            raise RegisteredResultAPIError(
                400,
                "catalog_cursor_invalid",
                "The cursor is invalid.",
            )
        head = self._repository.get_head(result_id)
        visible = self._repository.current_anchor_matches(
            project_id,
            result_id,
            include_archived=include_archived,
        )
        if head is None or head.stable_key != stable_key or not visible:
            raise RegisteredResultAPIError(
                404,
                "catalog_cursor_not_found",
                "The cursor was not found.",
            )
        return stable_key, result_id

    def _revisions_cursor(
        self,
        project_id: str,
        result_id: str,
        cursor: str | None,
    ) -> int | None:
        if cursor is None:
            return None
        payload = _decode_cursor(cursor)
        if (
            payload.get("kind") != "revisions"
            or payload.get("project_id") != project_id
            or payload.get("result_id") != result_id
            or set(payload)
            != {
                "kind",
                "project_id",
                "result_id",
                "revision_number",
            }
        ):
            raise RegisteredResultAPIError(
                400,
                "catalog_cursor_invalid",
                "The cursor is invalid.",
            )
        number = payload["revision_number"]
        if not isinstance(number, int) or isinstance(number, bool) or number < 1:
            raise RegisteredResultAPIError(
                400,
                "catalog_cursor_invalid",
                "The cursor is invalid.",
            )
        if not self._repository.revision_anchor_matches(
            project_id,
            result_id,
            number,
        ):
            raise RegisteredResultAPIError(
                404,
                "catalog_cursor_not_found",
                "The cursor was not found.",
            )
        return number

    def _project(self, project_id: str) -> str:
        try:
            return validate_project_id(project_id)
        except RegisteredResultPayloadError:
            raise self._missing() from None

    def _result(self, result_id: str) -> str:
        try:
            return validate_result_id(result_id)
        except RegisteredResultPayloadError:
            raise self._missing() from None

    def _idempotency(self, value: str) -> str:
        try:
            return validate_idempotency_key(value)
        except RegisteredResultPayloadError:
            raise self._invalid() from None

    def _revision_number(self, value: int) -> int:
        if not isinstance(value, int) or isinstance(value, bool) or value < 1:
            raise self._invalid()
        return value

    def _limit(self, value: int) -> int:
        if (
            not isinstance(value, int)
            or isinstance(value, bool)
            or not 1 <= value <= _MAX_PAGE
        ):
            raise RegisteredResultAPIError(
                400,
                "catalog_filter_invalid",
                "The list filter is invalid.",
            )
        return value

    def _from_store(self, exc: RegisteredResultStoreError) -> RegisteredResultAPIError:
        return _api(exc.code)

    def _missing(self) -> RegisteredResultAPIError:
        return _api("not_found")

    def _invalid(self) -> RegisteredResultAPIError:
        return _api("payload")

    def _unauthenticated(self) -> RegisteredResultAPIError:
        return _api("unauthenticated")

    def _bytes_unavailable(self) -> RegisteredResultAPIError:
        return _api("unavailable")


def _api(code: str) -> RegisteredResultAPIError:
    status_code, error_code, message = _STORE_ERRORS.get(
        code, _STORE_ERRORS["unavailable"]
    )
    return RegisteredResultAPIError(status_code, error_code, message)


_STORE_ERRORS = {
    "not_found": (404, "catalog_not_found", "Registered result was not found."),
    "forbidden": (
        403,
        "ADMINISTRATOR_REQUIRED",
        "The operation requires an administrator.",
    ),
    "unauthenticated": (
        401,
        "AUTHENTICATION_REQUIRED",
        "Authentication is required.",
    ),
    "idempotency": (
        409,
        "catalog_idempotency_conflict",
        "Idempotency-Key was reused with a different request.",
    ),
    "revision": (
        409,
        "catalog_revision_conflict",
        "The registered result changed.",
    ),
    "archived": (
        409,
        "catalog_already_archived",
        "The registered result is archived.",
    ),
    "payload": (
        400,
        "catalog_payload_invalid",
        "Registered result metadata is invalid.",
    ),
    "unavailable": (
        409,
        "registered_result_bytes_unavailable",
        "Registered result bytes are unavailable.",
    ),
}


def _normalized_part(field: str, value: object) -> object:
    if value is None:
        return None
    payload = {
        "display_name": "placeholder",
        "format": "txt",
        "origin": {"kind": "external", "source": {"status": "unknown"}},
        "method": {"status": "unknown"},
        "reference": {"status": "unknown"},
    }
    payload[field] = value
    return getattr(parse_metadata(payload), field)


def _chosen(override: object, current: object) -> object:
    if override is None:
        return current
    return override


def _encode_cursor(payload: Mapping[str, Any]) -> str:
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    token = base64.urlsafe_b64encode(raw).decode("ascii").rstrip("=")
    return f"{_CURSOR_PREFIX}{token}"


def _decode_cursor(token: str) -> dict[str, Any]:
    if (
        not isinstance(token, str)
        or not token.startswith(_CURSOR_PREFIX)
        or len(token) > 4096
    ):
        raise RegisteredResultAPIError(
            400,
            "catalog_cursor_invalid",
            "The cursor is invalid.",
        )
    body = token[len(_CURSOR_PREFIX) :]
    padded = body + ("=" * (-len(body) % 4))
    try:
        decoded = base64.urlsafe_b64decode(padded.encode("ascii"))
        payload = json.loads(decoded)
    except (UnicodeError, ValueError, json.JSONDecodeError):
        raise RegisteredResultAPIError(
            400,
            "catalog_cursor_invalid",
            "The cursor is invalid.",
        ) from None
    if not isinstance(payload, dict):
        raise RegisteredResultAPIError(
            400,
            "catalog_cursor_invalid",
            "The cursor is invalid.",
        )
    return payload


def _zulu(value: datetime) -> str:
    if value.tzinfo is None or value.utcoffset() is None:
        value = value.replace(tzinfo=timezone.utc)
    return value.astimezone(timezone.utc).isoformat().replace("+00:00", "Z")
