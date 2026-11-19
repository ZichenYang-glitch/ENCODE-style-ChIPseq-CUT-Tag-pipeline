"""SQLite records for registered results. Byte paths never leave this module."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
import json

from sqlalchemy import select, tuple_
from sqlalchemy.exc import IntegrityError
from sqlalchemy.orm import Session, sessionmaker

from encode_pipeline.persistence.models import (
    ProjectRow,
    RegisteredResultParentRow,
    RegisteredResultRevisionHeadRow,
    RegisteredResultRevisionRow,
    RegisteredResultRow,
    RegisteredResultWriteReceiptRow,
    UserAccountRow,
)
from encode_pipeline.platform.registered_results import (
    REVISION_DIGEST_SCHEME,
    ParentEdge,
    ParentView,
    RegisteredResultPayloadError,
    RevisionRecord,
    metadata_from_canonical_json,
)


class RegisteredResultStoreError(RuntimeError):
    """Closed failure code for one registered-result write or read."""

    def __init__(self, code: str) -> None:
        super().__init__(code)
        self.code = code


@dataclass(frozen=True)
class RevisionDraft:
    """One new revision the service has already digested."""

    result_id: str
    result_revision_id: str
    revision_number: int
    canonical_payload: str
    revision_digest: str
    created_at: datetime
    storage_pool_id: str
    relative_path: str
    size_bytes: int
    content_sha256: str
    parents: tuple[ParentEdge, ...]


class SqlAlchemyRegisteredResultRepository:
    """Transactional catalog writes for results that are not platform runs."""

    def __init__(self, session_factory: sessionmaker[Session]) -> None:
        self._session_factory = session_factory

    def get_revision(self, result_revision_id: str) -> RevisionRecord | None:
        with self._session_factory() as session:
            return _load_revision(session, result_revision_id)

    def get_head(self, result_id: str) -> RevisionRecord | None:
        with self._session_factory() as session:
            head = session.get(RegisteredResultRevisionHeadRow, result_id)
            if head is None:
                return None
            return _load_revision(session, head.result_revision_id)

    def visible_project(self, project_id: str) -> bool:
        with self._session_factory() as session:
            return _visible_project(session, project_id) is not None

    def result_in_project(self, project_id: str, result_id: str) -> bool:
        with self._session_factory() as session:
            result = session.get(RegisteredResultRow, result_id)
            return result is not None and result.project_id == project_id

    def list_current(
        self,
        project_id: str,
        *,
        include_archived: bool,
        after: tuple[str, str] | None,
        limit: int,
    ) -> tuple[RevisionRecord, ...]:
        with self._session_factory() as session:
            if _visible_project(session, project_id) is None:
                raise RegisteredResultStoreError("not_found")
            statement = (
                select(RegisteredResultRow)
                .where(RegisteredResultRow.project_id == project_id)
                .order_by(
                    RegisteredResultRow.stable_key,
                    RegisteredResultRow.result_id,
                )
            )
            if not include_archived:
                statement = statement.where(RegisteredResultRow.archived_at.is_(None))
            if after is not None:
                statement = statement.where(
                    tuple_(
                        RegisteredResultRow.stable_key,
                        RegisteredResultRow.result_id,
                    )
                    > tuple_(after[0], after[1])
                )
            rows = session.scalars(statement.limit(limit)).all()
            records: list[RevisionRecord] = []
            for row in rows:
                head = session.get(RegisteredResultRevisionHeadRow, row.result_id)
                if head is None:
                    raise RegisteredResultStoreError("unavailable")
                record = _load_revision(session, head.result_revision_id)
                if record is None:
                    raise RegisteredResultStoreError("unavailable")
                records.append(record)
            return tuple(records)

    def list_revisions(
        self,
        project_id: str,
        result_id: str,
        *,
        after_revision_number: int | None,
        limit: int,
    ) -> tuple[RevisionRecord, ...]:
        with self._session_factory() as session:
            if _visible_project(session, project_id) is None:
                raise RegisteredResultStoreError("not_found")
            result = session.get(RegisteredResultRow, result_id)
            if result is None or result.project_id != project_id:
                raise RegisteredResultStoreError("not_found")
            statement = (
                select(RegisteredResultRevisionRow)
                .where(RegisteredResultRevisionRow.result_id == result_id)
                .order_by(RegisteredResultRevisionRow.revision_number)
            )
            if after_revision_number is not None:
                statement = statement.where(
                    RegisteredResultRevisionRow.revision_number > after_revision_number
                )
            rows = session.scalars(statement.limit(limit)).all()
            return tuple(
                _require_record(session, row.result_revision_id) for row in rows
            )

    def current_anchor_matches(
        self,
        project_id: str,
        result_id: str,
        *,
        include_archived: bool,
    ) -> bool:
        with self._session_factory() as session:
            result = session.get(RegisteredResultRow, result_id)
            if result is None or result.project_id != project_id:
                return False
            return include_archived or result.archived_at is None

    def lookup_receipt(
        self,
        actor_user_id: str,
        idempotency_key: str,
        request_sha256: str,
    ) -> RevisionRecord | None:
        """Return a prior write when the key and request still match."""

        with self._session_factory() as session:
            return _replay(
                session,
                actor_user_id,
                idempotency_key,
                request_sha256,
            )

    def revision_anchor_matches(
        self,
        project_id: str,
        result_id: str,
        revision_number: int,
    ) -> bool:
        with self._session_factory() as session:
            result = session.get(RegisteredResultRow, result_id)
            if result is None or result.project_id != project_id:
                return False
            row = session.scalar(
                select(RegisteredResultRevisionRow).where(
                    RegisteredResultRevisionRow.result_id == result_id,
                    RegisteredResultRevisionRow.revision_number == revision_number,
                )
            )
            return row is not None

    def register(
        self,
        *,
        actor_user_id: str,
        idempotency_key: str,
        request_sha256: str,
        project_id: str,
        stable_key: str,
        draft: RevisionDraft,
    ) -> tuple[RevisionRecord, bool]:
        try:
            with self._session_factory.begin() as session:
                _begin_write(session)
                _require_writer(session, actor_user_id)
                replay = _replay(
                    session,
                    actor_user_id,
                    idempotency_key,
                    request_sha256,
                )
                if replay is not None:
                    return replay, True
                if _visible_project(session, project_id) is None:
                    raise RegisteredResultStoreError("not_found")
                if draft.revision_number != 1:
                    raise RegisteredResultStoreError("payload")
                taken = session.scalar(
                    select(RegisteredResultRow.result_id).where(
                        RegisteredResultRow.project_id == project_id,
                        RegisteredResultRow.stable_key == stable_key,
                    )
                )
                if taken is not None:
                    raise RegisteredResultStoreError("payload")
                _require_parent_pins(session, project_id, draft.parents)
                session.add(
                    RegisteredResultRow(
                        result_id=draft.result_id,
                        project_id=project_id,
                        stable_key=stable_key,
                        created_at=draft.created_at,
                        archived_at=None,
                    )
                )
                session.flush()
                _insert_revision(session, draft)
                session.flush()
                session.add(
                    RegisteredResultRevisionHeadRow(
                        result_id=draft.result_id,
                        result_revision_id=draft.result_revision_id,
                        revision_number=draft.revision_number,
                    )
                )
                _insert_receipt(
                    session,
                    actor_user_id=actor_user_id,
                    idempotency_key=idempotency_key,
                    request_sha256=request_sha256,
                    draft=draft,
                )
                session.flush()
                return _require_record(session, draft.result_revision_id), False
        except IntegrityError as exc:
            raise RegisteredResultStoreError("payload") from exc

    def revise(
        self,
        *,
        actor_user_id: str,
        idempotency_key: str,
        request_sha256: str,
        project_id: str,
        expected_revision_number: int,
        draft: RevisionDraft,
    ) -> tuple[RevisionRecord, bool]:
        try:
            with self._session_factory.begin() as session:
                _begin_write(session)
                _require_writer(session, actor_user_id)
                replay = _replay(
                    session,
                    actor_user_id,
                    idempotency_key,
                    request_sha256,
                )
                if replay is not None:
                    return replay, True
                result = _writable_result(session, project_id, draft.result_id)
                head = session.get(
                    RegisteredResultRevisionHeadRow,
                    result.result_id,
                )
                if (
                    head is None
                    or head.revision_number != expected_revision_number
                    or draft.revision_number != head.revision_number + 1
                ):
                    raise RegisteredResultStoreError("revision")
                _require_parent_pins(session, project_id, draft.parents)
                _insert_revision(session, draft)
                session.flush()
                head.result_revision_id = draft.result_revision_id
                head.revision_number = draft.revision_number
                _insert_receipt(
                    session,
                    actor_user_id=actor_user_id,
                    idempotency_key=idempotency_key,
                    request_sha256=request_sha256,
                    draft=draft,
                )
                session.flush()
                return _require_record(session, draft.result_revision_id), False
        except IntegrityError as exc:
            raise RegisteredResultStoreError("payload") from exc

    def archive(
        self,
        *,
        actor_user_id: str,
        idempotency_key: str,
        request_sha256: str,
        project_id: str,
        result_id: str,
        expected_revision_number: int,
        archived_at: datetime,
    ) -> tuple[RevisionRecord, bool]:
        try:
            with self._session_factory.begin() as session:
                _begin_write(session)
                _require_writer(session, actor_user_id)
                replay = _replay(
                    session,
                    actor_user_id,
                    idempotency_key,
                    request_sha256,
                )
                if replay is not None:
                    return replay, True
                result = _writable_result(session, project_id, result_id)
                head = session.get(RegisteredResultRevisionHeadRow, result.result_id)
                if head is None or head.revision_number != expected_revision_number:
                    raise RegisteredResultStoreError("revision")
                result.archived_at = archived_at
                session.add(
                    RegisteredResultWriteReceiptRow(
                        actor_user_id=actor_user_id,
                        idempotency_key=idempotency_key,
                        request_sha256=request_sha256,
                        result_id=result.result_id,
                        result_revision_id=head.result_revision_id,
                        created_at=archived_at,
                    )
                )
                session.flush()
                return _require_record(session, head.result_revision_id), False
        except IntegrityError as exc:
            raise RegisteredResultStoreError("payload") from exc


def _begin_write(session: Session) -> None:
    if session.get_bind().dialect.name == "sqlite":
        session.connection().exec_driver_sql("BEGIN IMMEDIATE")


def _require_writer(session: Session, actor_user_id: str) -> None:
    account = session.get(UserAccountRow, actor_user_id)
    if account is None or account.status != "enabled":
        raise RegisteredResultStoreError("unauthenticated")
    if account.role != "administrator":
        raise RegisteredResultStoreError("forbidden")


def _visible_project(session: Session, project_id: str) -> ProjectRow | None:
    project = session.get(ProjectRow, project_id)
    if project is None or project.archived_at is not None or project.kind != "user":
        return None
    return project


def _writable_result(
    session: Session,
    project_id: str,
    result_id: str,
) -> RegisteredResultRow:
    if _visible_project(session, project_id) is None:
        raise RegisteredResultStoreError("not_found")
    result = session.get(RegisteredResultRow, result_id)
    if result is None or result.project_id != project_id:
        raise RegisteredResultStoreError("not_found")
    if result.archived_at is not None:
        raise RegisteredResultStoreError("archived")
    return result


def _replay(
    session: Session,
    actor_user_id: str,
    idempotency_key: str,
    request_sha256: str,
) -> RevisionRecord | None:
    receipt = session.get(
        RegisteredResultWriteReceiptRow,
        (actor_user_id, idempotency_key),
    )
    if receipt is None:
        return None
    if receipt.request_sha256 != request_sha256:
        raise RegisteredResultStoreError("idempotency")
    return _require_record(session, receipt.result_revision_id)


def _require_parent_pins(
    session: Session,
    project_id: str,
    parents: tuple[ParentEdge, ...],
) -> None:
    for parent in parents:
        revision = session.get(
            RegisteredResultRevisionRow,
            parent.result_revision_id,
        )
        if revision is None or revision.revision_digest != parent.revision_digest:
            raise RegisteredResultStoreError("payload")
        result = session.get(RegisteredResultRow, revision.result_id)
        if result is None or result.project_id != project_id:
            raise RegisteredResultStoreError("payload")


def _insert_revision(session: Session, draft: RevisionDraft) -> None:
    session.add(
        RegisteredResultRevisionRow(
            result_revision_id=draft.result_revision_id,
            result_id=draft.result_id,
            revision_number=draft.revision_number,
            canonical_payload=draft.canonical_payload,
            digest_scheme=REVISION_DIGEST_SCHEME,
            revision_digest=draft.revision_digest,
            created_at=draft.created_at,
            storage_pool_id=draft.storage_pool_id,
            relative_path=draft.relative_path,
            size_bytes=draft.size_bytes,
            content_sha256=draft.content_sha256,
        )
    )
    if not draft.parents:
        return
    session.flush()
    for ordinal, parent in enumerate(draft.parents):
        session.add(
            RegisteredResultParentRow(
                child_result_revision_id=draft.result_revision_id,
                ordinal=ordinal,
                parent_result_revision_id=parent.result_revision_id,
                parent_revision_digest=parent.revision_digest,
                role=parent.role,
            )
        )


def _insert_receipt(
    session: Session,
    *,
    actor_user_id: str,
    idempotency_key: str,
    request_sha256: str,
    draft: RevisionDraft,
) -> None:
    session.add(
        RegisteredResultWriteReceiptRow(
            actor_user_id=actor_user_id,
            idempotency_key=idempotency_key,
            request_sha256=request_sha256,
            result_id=draft.result_id,
            result_revision_id=draft.result_revision_id,
            created_at=draft.created_at,
        )
    )


def _require_record(session: Session, result_revision_id: str) -> RevisionRecord:
    record = _load_revision(session, result_revision_id)
    if record is None:
        raise RegisteredResultStoreError("unavailable")
    return record


def _load_revision(session: Session, result_revision_id: str) -> RevisionRecord | None:
    revision = session.get(RegisteredResultRevisionRow, result_revision_id)
    if revision is None:
        return None
    result = session.get(RegisteredResultRow, revision.result_id)
    if result is None:
        return None
    project = session.get(ProjectRow, result.project_id)
    head = session.get(RegisteredResultRevisionHeadRow, result.result_id)
    if project is None or head is None:
        return None
    parent_rows = session.scalars(
        select(RegisteredResultParentRow)
        .where(
            RegisteredResultParentRow.child_result_revision_id
            == revision.result_revision_id
        )
        .order_by(RegisteredResultParentRow.ordinal)
    ).all()
    parents: list[ParentView] = []
    for parent in parent_rows:
        parent_revision = session.get(
            RegisteredResultRevisionRow,
            parent.parent_result_revision_id,
        )
        if parent_revision is None:
            raise RegisteredResultStoreError("unavailable")
        parent_result = session.get(RegisteredResultRow, parent_revision.result_id)
        parent_head = session.get(
            RegisteredResultRevisionHeadRow,
            parent_revision.result_id,
        )
        if parent_result is None or parent_head is None:
            raise RegisteredResultStoreError("unavailable")
        try:
            metadata = metadata_from_canonical_json(parent_revision.canonical_payload)
        except RegisteredResultPayloadError as exc:
            raise RegisteredResultStoreError("unavailable") from exc
        status = (
            "current"
            if parent_head.result_revision_id == parent_revision.result_revision_id
            else "superseded"
        )
        parents.append(
            ParentView(
                result_id=parent_result.result_id,
                result_revision_id=parent_revision.result_revision_id,
                revision_number=parent_revision.revision_number,
                revision_status=status,
                current_result_revision_id=parent_head.result_revision_id,
                display_name=metadata.display_name,
                role=parent.role,
                revision_digest=parent.parent_revision_digest,
            )
        )
    try:
        json.loads(revision.canonical_payload)
    except json.JSONDecodeError as exc:
        raise RegisteredResultStoreError("unavailable") from exc
    return RevisionRecord(
        project_id=result.project_id,
        project_archived_at=project.archived_at,
        result_id=result.result_id,
        stable_key=result.stable_key,
        result_archived_at=result.archived_at,
        result_revision_id=revision.result_revision_id,
        revision_number=revision.revision_number,
        canonical_payload=revision.canonical_payload,
        revision_digest=revision.revision_digest,
        created_at=revision.created_at,
        storage_pool_id=revision.storage_pool_id,
        relative_path=revision.relative_path,
        size_bytes=revision.size_bytes,
        content_sha256=revision.content_sha256,
        head_result_revision_id=head.result_revision_id,
        parents=tuple(parents),
    )
