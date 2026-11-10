"""SQLite repository for authoring drafts, turns, and adoption receipts.

Write methods sample the clock only after ``BEGIN IMMEDIATE``. Conditional
updates decide completion, invalidation, and adoption races.
"""

from __future__ import annotations

import base64
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from secrets import token_hex, token_urlsafe
from typing import Any

from sqlalchemy import and_, delete, or_, select, update
from sqlalchemy.orm import Session, sessionmaker

from encode_pipeline.persistence.models import (
    AgentAuthoringSessionRow,
    AgentAuthoringTurnRow,
    AuthSessionRow,
    ProjectRow,
    UserAccountRow,
)
from encode_pipeline.platform.authoring import (
    ADMISSION_FAILED,
    AUTH_SESSION_ENDED,
    AUTHENTICATION_REQUIRED,
    AuthoringActor,
    AuthoringRejected,
    CAPABILITY_CHANGED,
    INVALID_REQUEST,
    CLOSED_BY_OWNER,
    CONNECTION_NOT_READY,
    CONNECTION_OWNED_BY_OTHER_LOGIN,
    DECISION_CONFLICT,
    DRAFT_EXPIRED,
    DRAFT_TTL,
    EXPIRED_RECORD_RETENTION,
    NOT_FOUND,
    PROPOSAL_EXPIRED,
    PROPOSAL_INVALIDATED,
    PROPOSAL_TTL,
    STALE_DRAFT,
    TURN_REQUEST_CONFLICT,
    TURN_RESULT_CONFLICT,
    draft_digest,
    empty_draft_content,
    sha256_text,
    utc_now,
)
from encode_pipeline.platform.data_registry import LEGACY_PROJECT_ID


Clock = Callable[[], datetime]
BeforeCommit = Callable[[], None]


@dataclass(frozen=True)
class AuthoringSessionRecord:
    sid: str
    owner_user_id: str
    project_id: str
    workflow_id: str
    schema_digest: str
    projection_digest: str | None
    content: dict[str, Any]
    draft_digest: str
    revision: int
    created_at: datetime
    content_updated_at: datetime
    draft_expires_at: datetime
    purge_after: datetime
    last_snapshot_id: str | None
    last_validated_revision: int | None
    validation_issue_codes: tuple[str, ...]
    validation_issues: tuple[dict[str, Any], ...]
    connection_state: str
    connection_reason_code: str | None
    originating_auth_session_ref: str | None
    instance_id: str | None
    generation: int
    bridge_owner_epoch: str | None
    connection_capability_digest: str | None


@dataclass(frozen=True)
class AuthoringTurnRecord:
    turn_id: str
    sid: str
    client_turn_id: str
    request_digest: str
    question: dict[str, Any]
    originating_auth_session_ref: str
    instance_id: str
    generation: int
    bridge_owner_epoch: str
    base_revision: int
    schema_digest: str
    projection_digest: str
    status: str
    failure_reason: str | None
    completed_at: datetime | None
    result_kind: str | None
    result_digest: str | None
    proposal_id: str | None
    proposal: dict[str, Any] | None
    proposal_digest: str | None
    proposal_expires_at: datetime | None
    decision: str | None
    invalidation_reason: str | None
    adopt_request_digest: str | None
    adoption_actor_user_id: str | None
    adoption_auth_session_ref: str | None
    adopted_at: datetime | None
    result_revision: int | None
    adoption_changed: bool | None
    adoption_content_digest: str | None
    adoption_draft_expires_at: datetime | None


@dataclass(frozen=True)
class AuthoringSessionPage:
    records: tuple[AuthoringSessionRecord, ...]
    next_cursor: str | None


@dataclass(frozen=True)
class ConnectionClaim:
    record: AuthoringSessionRecord
    reused: bool
    capability_token: str | None


@dataclass(frozen=True)
class TurnCompletion:
    kind: str
    result_digest: str
    proposal_id: str | None
    proposal: dict[str, Any] | None
    proposal_digest: str | None


@dataclass(frozen=True)
class AdoptionResult:
    record: AuthoringSessionRecord
    turn: AuthoringTurnRecord
    replayed: bool


class SqlAlchemyAuthoringRepository:
    """Conditional authoring store sharing the platform SQLite database."""

    def __init__(
        self,
        session_factory: sessionmaker[Session],
        *,
        clock: Clock | None = None,
        draft_ttl: timedelta = DRAFT_TTL,
        proposal_ttl: timedelta = PROPOSAL_TTL,
        retention: timedelta = EXPIRED_RECORD_RETENTION,
        before_commit: BeforeCommit | None = None,
        session_id_factory: Callable[[], str] | None = None,
        turn_id_factory: Callable[[], str] | None = None,
        proposal_id_factory: Callable[[], str] | None = None,
    ) -> None:
        self._session_factory = session_factory
        self._clock = clock
        self._draft_ttl = draft_ttl
        self._proposal_ttl = proposal_ttl
        self._retention = retention
        self._before_commit = before_commit
        self._session_ids = session_id_factory or (lambda: f"aas_{token_hex(16)}")
        self._turn_ids = turn_id_factory or (lambda: f"atr_{token_hex(16)}")
        self._proposal_ids = proposal_id_factory or (lambda: f"prop_{token_hex(16)}")

    def create_session(
        self,
        actor: AuthoringActor,
        *,
        project_id: str,
        workflow_id: str,
        schema_digest: str,
        projection_digest: str | None,
    ) -> AuthoringSessionRecord:
        content = empty_draft_content()
        with self._session_factory.begin() as session:
            _begin_write(session)
            now = self._now()
            self._require_actor(session, actor, now)
            if not _project_visible(session, project_id):
                raise AuthoringRejected(NOT_FOUND)
            expires = now + self._draft_ttl
            row = AgentAuthoringSessionRow(
                sid=self._session_ids(),
                owner_user_id=actor.user_id,
                project_id=project_id,
                workflow_id=workflow_id,
                schema_digest=schema_digest,
                projection_digest=projection_digest,
                draft_json=content,
                draft_digest=draft_digest(content),
                draft_revision=0,
                created_at=now,
                content_updated_at=now,
                draft_expires_at=expires,
                purge_after=expires + self._retention,
                connection_state="closed",
                generation=0,
            )
            session.add(row)
            session.flush()
            self._hook()
            return _session_record(row)

    def get_session(
        self,
        actor: AuthoringActor,
        sid: str,
    ) -> AuthoringSessionRecord:
        record: AuthoringSessionRecord | None = None
        with self._session_factory.begin() as session:
            _begin_write(session)
            now = self._now()
            self._require_actor(session, actor, now)
            row = session.get(AgentAuthoringSessionRow, sid)
            if (
                row is None
                or row.owner_user_id != actor.user_id
                or not _project_visible(session, row.project_id)
            ):
                record = None
            elif _as_utc(row.purge_after) <= now:
                _delete_session(session, sid)
                record = None
            else:
                record = _session_record(row)
            self._hook()
        if record is None:
            raise AuthoringRejected(NOT_FOUND)
        return record

    def list_sessions(
        self,
        actor: AuthoringActor,
        *,
        project_id: str | None,
        workflow_id: str | None,
        cursor_time: datetime | None,
        cursor_sid: str | None,
        limit: int,
    ) -> AuthoringSessionPage:
        with self._session_factory.begin() as session:
            _begin_write(session)
            now = self._now()
            self._require_actor(session, actor, now)
            self._purge_due(session, actor.user_id, now, limit=20)
            statement = (
                select(AgentAuthoringSessionRow)
                .join(
                    ProjectRow,
                    ProjectRow.project_id == AgentAuthoringSessionRow.project_id,
                )
                .where(
                    AgentAuthoringSessionRow.owner_user_id == actor.user_id,
                    AgentAuthoringSessionRow.purge_after > now,
                    ProjectRow.archived_at.is_(None),
                    ProjectRow.kind == "user",
                    ProjectRow.project_id != LEGACY_PROJECT_ID,
                )
            )
            if project_id is not None:
                statement = statement.where(
                    AgentAuthoringSessionRow.project_id == project_id
                )
            if workflow_id is not None:
                statement = statement.where(
                    AgentAuthoringSessionRow.workflow_id == workflow_id
                )
            if cursor_time is not None and cursor_sid is not None:
                statement = statement.where(
                    or_(
                        AgentAuthoringSessionRow.content_updated_at < cursor_time,
                        and_(
                            AgentAuthoringSessionRow.content_updated_at == cursor_time,
                            AgentAuthoringSessionRow.sid < cursor_sid,
                        ),
                    )
                )
            rows = session.scalars(
                statement.order_by(
                    AgentAuthoringSessionRow.content_updated_at.desc(),
                    AgentAuthoringSessionRow.sid.desc(),
                ).limit(limit + 1)
            ).all()
            self._hook()
            visible = tuple(_session_record(row) for row in rows[:limit])
            next_cursor = None
            if len(rows) > limit and visible:
                last = visible[-1]
                next_cursor = encode_authoring_cursor(
                    last.content_updated_at,
                    last.sid,
                )
            return AuthoringSessionPage(visible, next_cursor)

    def save_draft(
        self,
        actor: AuthoringActor,
        sid: str,
        *,
        expected_revision: int,
        content: Mapping[str, Any],
        schema_digest: str,
        projection_digest: str | None,
    ) -> tuple[AuthoringSessionRecord, bool]:
        normalized = dict(content)
        digest = draft_digest(normalized)
        with self._session_factory.begin() as session:
            _begin_write(session)
            now = self._now()
            row = self._owned(session, actor, sid, now, allow_expired=False)
            if row.draft_revision != expected_revision:
                raise AuthoringRejected(STALE_DRAFT)
            changed = row.draft_digest != digest
            schema_changed = (
                row.schema_digest != schema_digest
                or row.projection_digest != projection_digest
            )
            if changed:
                row.draft_json = normalized
                row.draft_digest = digest
                row.draft_revision += 1
                row.content_updated_at = now
                row.draft_expires_at = now + self._draft_ttl
                row.purge_after = row.draft_expires_at + self._retention
                _clear_validation(row)
                _invalidate_pending(session, sid, STALE_DRAFT)
            if schema_changed:
                row.schema_digest = schema_digest
                row.projection_digest = projection_digest
                if not changed:
                    _invalidate_pending(session, sid, PROPOSAL_INVALIDATED)
            session.flush()
            self._hook()
            return _session_record(row), changed

    def open_connection(
        self,
        actor: AuthoringActor,
        sid: str,
        *,
        owner_epoch: str,
        projection_digest: str,
    ) -> ConnectionClaim:
        with self._session_factory.begin() as session:
            _begin_write(session)
            now = self._now()
            row = self._owned(session, actor, sid, now, allow_expired=False)
            if row.projection_digest != projection_digest:
                raise AuthoringRejected(CAPABILITY_CHANGED)
            if row.connection_state in {"blocked", "ready"}:
                claim = self._reuse_or_fence(session, row, actor, owner_epoch, now)
                if claim is not None:
                    self._hook()
                    return claim
                session.refresh(row)
            token = token_urlsafe(32)
            row.connection_state = "blocked"
            row.generation += 1
            row.originating_auth_session_ref = actor.session_digest
            row.bridge_owner_epoch = owner_epoch
            row.instance_id = None
            row.connection_reason_code = None
            row.connection_capability_digest = sha256_text(token)
            session.flush()
            self._hook()
            return ConnectionClaim(_session_record(row), False, token)

    def mark_connection_ready(
        self,
        actor: AuthoringActor,
        sid: str,
        *,
        generation: int,
        instance_id: str,
        owner_epoch: str,
    ) -> AuthoringSessionRecord:
        with self._session_factory.begin() as session:
            _begin_write(session)
            now = self._now()
            row = self._owned(session, actor, sid, now, allow_expired=False)
            result = session.execute(
                update(AgentAuthoringSessionRow)
                .where(
                    AgentAuthoringSessionRow.sid == sid,
                    AgentAuthoringSessionRow.connection_state == "blocked",
                    AgentAuthoringSessionRow.generation == generation,
                    AgentAuthoringSessionRow.bridge_owner_epoch == owner_epoch,
                    AgentAuthoringSessionRow.originating_auth_session_ref
                    == actor.session_digest,
                )
                .values(connection_state="ready", instance_id=instance_id)
            )
            if result.rowcount != 1:
                raise AuthoringRejected(ADMISSION_FAILED)
            session.refresh(row)
            self._hook()
            return _session_record(row)

    def fail_connection(
        self,
        sid: str,
        *,
        generation: int,
        owner_epoch: str,
        reason_code: str,
    ) -> None:
        with self._session_factory.begin() as session:
            _begin_write(session)
            self._now()
            _fence_generation(
                session,
                sid,
                generation=generation,
                owner_epoch=owner_epoch,
                reason_code=reason_code,
            )
            self._hook()

    def close_connection(
        self,
        actor: AuthoringActor,
        sid: str,
    ) -> AuthoringSessionRecord:
        with self._session_factory.begin() as session:
            _begin_write(session)
            now = self._now()
            row = self._owned(session, actor, sid, now, allow_expired=True)
            if row.connection_state != "closed":
                _fence_generation(
                    session,
                    sid,
                    generation=row.generation,
                    owner_epoch=row.bridge_owner_epoch,
                    reason_code=CLOSED_BY_OWNER,
                )
                session.refresh(row)
            self._hook()
            return _session_record(row)

    def retire_foreign_epochs(self, owner_epoch: str) -> int:
        """Fence ready connections that do not belong to the current bridge owner."""

        with self._session_factory.begin() as session:
            _begin_write(session)
            self._now()
            rows = session.scalars(
                select(AgentAuthoringSessionRow).where(
                    AgentAuthoringSessionRow.connection_state.in_(("blocked", "ready")),
                    or_(
                        AgentAuthoringSessionRow.bridge_owner_epoch.is_(None),
                        AgentAuthoringSessionRow.bridge_owner_epoch != owner_epoch,
                    ),
                )
            ).all()
            retired = 0
            for row in rows:
                if _fence_generation(
                    session,
                    row.sid,
                    generation=row.generation,
                    owner_epoch=row.bridge_owner_epoch,
                    reason_code=AUTH_SESSION_ENDED,
                ):
                    retired += 1
            self._hook()
            return retired

    def retire_owner_epoch(self, owner_epoch: str) -> int:
        with self._session_factory.begin() as session:
            _begin_write(session)
            self._now()
            rows = session.scalars(
                select(AgentAuthoringSessionRow).where(
                    AgentAuthoringSessionRow.bridge_owner_epoch == owner_epoch,
                    AgentAuthoringSessionRow.connection_state.in_(("blocked", "ready")),
                )
            ).all()
            retired = 0
            for row in rows:
                if _fence_generation(
                    session,
                    row.sid,
                    generation=row.generation,
                    owner_epoch=owner_epoch,
                    reason_code=AUTH_SESSION_ENDED,
                ):
                    retired += 1
            self._hook()
            return retired

    def start_turn(
        self,
        actor: AuthoringActor,
        sid: str,
        *,
        client_turn_id: str,
        request_digest: str,
        question: Mapping[str, Any],
        expected_revision: int,
        schema_digest: str,
        projection_digest: str,
        owner_epoch: str,
    ) -> tuple[AuthoringTurnRecord, bool]:
        with self._session_factory.begin() as session:
            _begin_write(session)
            now = self._now()
            row = self._owned(session, actor, sid, now, allow_expired=False)
            self._require_ready(session, row, actor, owner_epoch, now)
            if row.draft_revision != expected_revision:
                raise AuthoringRejected(STALE_DRAFT)
            if (
                row.schema_digest != schema_digest
                or row.projection_digest != projection_digest
            ):
                raise AuthoringRejected(CAPABILITY_CHANGED)
            existing = session.scalar(
                select(AgentAuthoringTurnRow).where(
                    AgentAuthoringTurnRow.sid == sid,
                    AgentAuthoringTurnRow.client_turn_id == client_turn_id,
                )
            )
            if existing is not None:
                if existing.request_digest != request_digest:
                    raise AuthoringRejected(TURN_REQUEST_CONFLICT)
                self._hook()
                return _turn_record(existing), True
            if (
                session.scalar(
                    select(AgentAuthoringTurnRow.turn_id).where(
                        AgentAuthoringTurnRow.sid == sid,
                        AgentAuthoringTurnRow.status == "generating",
                    )
                )
                is not None
            ):
                raise AuthoringRejected(CONNECTION_NOT_READY)
            turn = AgentAuthoringTurnRow(
                turn_id=self._turn_ids(),
                sid=sid,
                client_turn_id=client_turn_id,
                request_digest=request_digest,
                question_json=dict(question),
                originating_auth_session_ref=actor.session_digest,
                instance_id=row.instance_id,
                generation=row.generation,
                bridge_owner_epoch=owner_epoch,
                base_revision=row.draft_revision,
                schema_digest=schema_digest,
                projection_digest=projection_digest,
                status="generating",
            )
            session.add(turn)
            session.flush()
            self._hook()
            return _turn_record(turn), False

    def complete_turn(
        self,
        actor: AuthoringActor,
        turn_id: str,
        completion: TurnCompletion,
    ) -> AuthoringTurnRecord:
        with self._session_factory.begin() as session:
            _begin_write(session)
            now = self._now()
            turn = session.get(AgentAuthoringTurnRow, turn_id)
            if turn is None:
                raise AuthoringRejected(NOT_FOUND)
            row = session.get(AgentAuthoringSessionRow, turn.sid)
            if row is None or row.owner_user_id != actor.user_id:
                raise AuthoringRejected(NOT_FOUND)
            authorized = self._callback_authorized(session, row, turn, actor, now)
            if turn.status != "generating":
                if not authorized:
                    raise AuthoringRejected(AUTHENTICATION_REQUIRED)
                if (
                    turn.status == "completed"
                    and turn.result_digest == completion.result_digest
                ):
                    self._hook()
                    return _turn_record(turn)
                raise AuthoringRejected(TURN_RESULT_CONFLICT)
            failure = self._completion_failure(session, row, turn, actor, now)
            if failure is not None:
                turn.status = "failed"
                turn.failure_reason = failure
                turn.completed_at = now
                session.flush()
                self._hook()
                return _turn_record(turn)
            expires = (
                now + self._proposal_ttl if completion.kind == "proposal" else None
            )
            decision = "pending" if completion.kind == "proposal" else None
            result = session.execute(
                update(AgentAuthoringTurnRow)
                .where(
                    AgentAuthoringTurnRow.turn_id == turn_id,
                    AgentAuthoringTurnRow.status == "generating",
                    AgentAuthoringTurnRow.generation == turn.generation,
                )
                .values(
                    status="completed",
                    completed_at=now,
                    result_kind=completion.kind,
                    result_digest=completion.result_digest,
                    proposal_id=completion.proposal_id,
                    proposal_json=completion.proposal,
                    proposal_digest=completion.proposal_digest,
                    proposal_expires_at=expires,
                    decision=decision,
                    failure_reason=None,
                )
            )
            if result.rowcount != 1:
                raise AuthoringRejected(TURN_RESULT_CONFLICT)
            session.refresh(turn)
            self._hook()
            return _turn_record(turn)

    def fail_turn(
        self,
        turn_id: str,
        *,
        generation: int,
        reason_code: str,
    ) -> None:
        with self._session_factory.begin() as session:
            _begin_write(session)
            now = self._now()
            session.execute(
                update(AgentAuthoringTurnRow)
                .where(
                    AgentAuthoringTurnRow.turn_id == turn_id,
                    AgentAuthoringTurnRow.status == "generating",
                    AgentAuthoringTurnRow.generation == generation,
                )
                .values(
                    status="failed",
                    failure_reason=reason_code,
                    completed_at=now,
                )
            )
            self._hook()

    def get_proposal(
        self,
        actor: AuthoringActor,
        sid: str,
        proposal_id: str,
    ) -> tuple[AuthoringSessionRecord, AuthoringTurnRecord]:
        with self._session_factory.begin() as session:
            _begin_write(session)
            now = self._now()
            row = self._owned(session, actor, sid, now, allow_expired=True)
            turn = session.scalar(
                select(AgentAuthoringTurnRow).where(
                    AgentAuthoringTurnRow.sid == sid,
                    AgentAuthoringTurnRow.proposal_id == proposal_id,
                )
            )
            if turn is None:
                raise AuthoringRejected(NOT_FOUND)
            self._hook()
            return _session_record(row), _turn_record(turn)

    def peek_session(self, sid: str) -> AuthoringSessionRecord | None:
        with self._session_factory.begin() as session:
            _begin_write(session)
            now = self._now()
            row = session.get(AgentAuthoringSessionRow, sid)
            if row is None:
                return None
            if _as_utc(row.purge_after) <= now:
                _delete_session(session, sid)
                self._hook()
                return None
            self._hook()
            return _session_record(row)

    def peek_turn(self, sid: str, turn_id: str) -> AuthoringTurnRecord | None:
        with self._session_factory.begin() as session:
            _begin_write(session)
            turn = session.get(AgentAuthoringTurnRow, turn_id)
            self._hook()
            if turn is None or turn.sid != sid:
                return None
            return _turn_record(turn)

    def get_turn(
        self,
        actor: AuthoringActor,
        sid: str,
        turn_id: str,
    ) -> tuple[AuthoringSessionRecord, AuthoringTurnRecord]:
        with self._session_factory.begin() as session:
            _begin_write(session)
            now = self._now()
            row = self._owned(session, actor, sid, now, allow_expired=True)
            turn = session.get(AgentAuthoringTurnRow, turn_id)
            if turn is None or turn.sid != sid:
                raise AuthoringRejected(NOT_FOUND)
            self._hook()
            return _session_record(row), _turn_record(turn)

    def adopt(
        self,
        actor: AuthoringActor,
        sid: str,
        proposal_id: str,
        *,
        expected_revision: int,
        proposal_digest_value: str,
        request_digest: str,
        content: Mapping[str, Any] | None,
        content_digest: str | None,
        changed: bool | None,
    ) -> AdoptionResult:
        with self._session_factory.begin() as session:
            _begin_write(session)
            now = self._now()
            self._require_actor(session, actor, now)
            row = session.get(AgentAuthoringSessionRow, sid)
            if (
                row is None
                or row.owner_user_id != actor.user_id
                or not _project_visible(session, row.project_id)
            ):
                raise AuthoringRejected(NOT_FOUND)
            if _as_utc(row.purge_after) <= now:
                _delete_session(session, sid)
                raise AuthoringRejected(NOT_FOUND)
            turn = session.scalar(
                select(AgentAuthoringTurnRow).where(
                    AgentAuthoringTurnRow.sid == sid,
                    AgentAuthoringTurnRow.proposal_id == proposal_id,
                )
            )
            if turn is None:
                raise AuthoringRejected(NOT_FOUND)
            if turn.decision == "adopted":
                if turn.adopt_request_digest != request_digest:
                    raise AuthoringRejected(
                        DECISION_CONFLICT,
                        detail={"proposal_state": "adopted"},
                    )
                self._hook()
                return AdoptionResult(_session_record(row), _turn_record(turn), True)
            if turn.decision == "rejected":
                raise AuthoringRejected(
                    DECISION_CONFLICT,
                    detail={"proposal_state": "rejected"},
                )
            if turn.decision != "pending":
                raise AuthoringRejected(PROPOSAL_INVALIDATED)
            if _as_utc(row.draft_expires_at) <= now:
                raise AuthoringRejected(DRAFT_EXPIRED)
            if (
                turn.proposal_expires_at is None
                or _as_utc(turn.proposal_expires_at) <= now
            ):
                raise AuthoringRejected(PROPOSAL_EXPIRED)
            if not self._proposal_binding_ok(session, row, turn, actor, now):
                raise AuthoringRejected(PROPOSAL_INVALIDATED)
            if row.draft_revision != turn.base_revision or (
                expected_revision != turn.base_revision
            ):
                raise AuthoringRejected(STALE_DRAFT)
            if turn.proposal_digest != proposal_digest_value:
                raise AuthoringRejected(DECISION_CONFLICT)
            if changed is None or content is None or content_digest is None:
                raise AuthoringRejected(INVALID_REQUEST)
            if changed:
                updated = session.execute(
                    update(AgentAuthoringSessionRow)
                    .where(
                        AgentAuthoringSessionRow.sid == sid,
                        AgentAuthoringSessionRow.draft_revision == turn.base_revision,
                    )
                    .values(
                        draft_json=dict(content),
                        draft_digest=content_digest,
                        draft_revision=turn.base_revision + 1,
                        content_updated_at=now,
                        draft_expires_at=now + self._draft_ttl,
                        purge_after=now + self._draft_ttl + self._retention,
                        last_snapshot_id=None,
                        last_validated_revision=None,
                        validation_issue_codes=None,
                        validation_issues=None,
                    )
                )
                if updated.rowcount != 1:
                    raise AuthoringRejected(STALE_DRAFT)
                _invalidate_pending(
                    session,
                    sid,
                    STALE_DRAFT,
                    except_turn_id=turn.turn_id,
                )
                result_revision = turn.base_revision + 1
                draft_expires = now + self._draft_ttl
            else:
                result_revision = row.draft_revision
                draft_expires = _as_utc(row.draft_expires_at)
            decided = session.execute(
                update(AgentAuthoringTurnRow)
                .where(
                    AgentAuthoringTurnRow.turn_id == turn.turn_id,
                    AgentAuthoringTurnRow.decision == "pending",
                )
                .values(
                    decision="adopted",
                    adopt_request_digest=request_digest,
                    adoption_actor_user_id=actor.user_id,
                    adoption_auth_session_ref=actor.session_digest,
                    adopted_at=now,
                    result_revision=result_revision,
                    adoption_changed=changed,
                    adoption_content_digest=content_digest,
                    adoption_draft_expires_at=draft_expires,
                    invalidation_reason=None,
                )
            )
            if decided.rowcount != 1:
                raise AuthoringRejected(DECISION_CONFLICT)
            session.refresh(row)
            session.refresh(turn)
            self._hook()
            return AdoptionResult(_session_record(row), _turn_record(turn), False)

    def reject(
        self,
        actor: AuthoringActor,
        sid: str,
        proposal_id: str,
        *,
        proposal_digest_value: str,
    ) -> AuthoringTurnRecord:
        with self._session_factory.begin() as session:
            _begin_write(session)
            now = self._now()
            row = self._owned(session, actor, sid, now, allow_expired=True)
            turn = session.scalar(
                select(AgentAuthoringTurnRow).where(
                    AgentAuthoringTurnRow.sid == sid,
                    AgentAuthoringTurnRow.proposal_id == proposal_id,
                )
            )
            if turn is None:
                raise AuthoringRejected(NOT_FOUND)
            if turn.proposal_digest != proposal_digest_value:
                raise AuthoringRejected(DECISION_CONFLICT)
            if turn.decision == "rejected":
                self._hook()
                return _turn_record(turn)
            if _as_utc(row.draft_expires_at) <= now:
                raise AuthoringRejected(DRAFT_EXPIRED)
            if turn.decision == "adopted":
                raise AuthoringRejected(
                    DECISION_CONFLICT,
                    detail={"proposal_state": "adopted"},
                )
            if turn.decision != "pending":
                raise AuthoringRejected(PROPOSAL_INVALIDATED)
            if (
                turn.proposal_expires_at is None
                or _as_utc(turn.proposal_expires_at) <= now
            ):
                raise AuthoringRejected(PROPOSAL_EXPIRED)
            del row
            decided = session.execute(
                update(AgentAuthoringTurnRow)
                .where(
                    AgentAuthoringTurnRow.turn_id == turn.turn_id,
                    AgentAuthoringTurnRow.decision == "pending",
                )
                .values(decision="rejected")
            )
            if decided.rowcount != 1:
                raise AuthoringRejected(DECISION_CONFLICT)
            session.refresh(turn)
            self._hook()
            return _turn_record(turn)

    def login_is_live(self, user_id: str, session_digest: str) -> bool:
        with self._session_factory.begin() as session:
            _begin_write(session)
            now = self._now()
            self._hook()
            return _auth_live(session, session_digest, user_id, now)

    def load_for_validation(
        self,
        actor: AuthoringActor,
        sid: str,
        *,
        expected_revision: int,
    ) -> AuthoringSessionRecord:
        with self._session_factory.begin() as session:
            _begin_write(session)
            now = self._now()
            row = self._owned(session, actor, sid, now, allow_expired=False)
            if row.draft_revision != expected_revision:
                raise AuthoringRejected(STALE_DRAFT)
            self._hook()
            return _session_record(row)

    def store_validation(
        self,
        actor: AuthoringActor,
        sid: str,
        *,
        expected_revision: int,
        schema_digest: str,
        snapshot_id: str | None,
        issue_codes: Sequence[str],
        issues: Sequence[Mapping[str, Any]],
    ) -> AuthoringSessionRecord:
        with self._session_factory.begin() as session:
            _begin_write(session)
            now = self._now()
            row = self._owned(session, actor, sid, now, allow_expired=False)
            if row.draft_revision != expected_revision:
                raise AuthoringRejected(STALE_DRAFT)
            if row.schema_digest != schema_digest:
                raise AuthoringRejected(CAPABILITY_CHANGED)
            row.last_snapshot_id = snapshot_id
            row.last_validated_revision = expected_revision
            row.validation_issue_codes = list(issue_codes)
            row.validation_issues = [dict(issue) for issue in issues]
            session.flush()
            self._hook()
            return _session_record(row)

    def _reuse_or_fence(
        self,
        session: Session,
        row: AgentAuthoringSessionRow,
        actor: AuthoringActor,
        owner_epoch: str,
        now: datetime,
    ) -> ConnectionClaim | None:
        occupying = row.originating_auth_session_ref
        same_login = occupying == actor.session_digest
        occupying_live = occupying is not None and _auth_live(
            session,
            occupying,
            row.owner_user_id,
            now,
        )
        same_epoch = row.bridge_owner_epoch == owner_epoch
        if same_login and occupying_live and same_epoch:
            return ConnectionClaim(_session_record(row), True, None)
        if not same_login and occupying_live:
            raise AuthoringRejected(CONNECTION_OWNED_BY_OTHER_LOGIN)
        _fence_generation(
            session,
            row.sid,
            generation=row.generation,
            owner_epoch=row.bridge_owner_epoch,
            reason_code=AUTH_SESSION_ENDED,
        )
        return None

    def _require_ready(
        self,
        session: Session,
        row: AgentAuthoringSessionRow,
        actor: AuthoringActor,
        owner_epoch: str,
        now: datetime,
    ) -> None:
        if row.connection_state != "ready" or row.instance_id is None:
            raise AuthoringRejected(CONNECTION_NOT_READY)
        if row.originating_auth_session_ref != actor.session_digest:
            if row.originating_auth_session_ref and _auth_live(
                session,
                row.originating_auth_session_ref,
                row.owner_user_id,
                now,
            ):
                raise AuthoringRejected(CONNECTION_OWNED_BY_OTHER_LOGIN)
            raise AuthoringRejected(CONNECTION_NOT_READY)
        if row.bridge_owner_epoch != owner_epoch:
            raise AuthoringRejected(CONNECTION_NOT_READY)
        if not _auth_live(session, actor.session_digest, actor.user_id, now):
            raise AuthoringRejected(AUTHENTICATION_REQUIRED)

    def _callback_authorized(
        self,
        session: Session,
        row: AgentAuthoringSessionRow,
        turn: AgentAuthoringTurnRow,
        actor: AuthoringActor,
        now: datetime,
    ) -> bool:
        return self._completion_failure(session, row, turn, actor, now) is None

    def _completion_failure(
        self,
        session: Session,
        row: AgentAuthoringSessionRow,
        turn: AgentAuthoringTurnRow,
        actor: AuthoringActor,
        now: datetime,
    ) -> str | None:
        if actor.session_digest != turn.originating_auth_session_ref:
            return AUTHENTICATION_REQUIRED
        if not _auth_live(
            session, turn.originating_auth_session_ref, actor.user_id, now
        ):
            return AUTH_SESSION_ENDED
        if _as_utc(row.draft_expires_at) <= now or _as_utc(row.purge_after) <= now:
            return DRAFT_EXPIRED
        if row.owner_user_id != actor.user_id or not _project_visible(
            session, row.project_id
        ):
            return NOT_FOUND
        if (
            row.connection_state != "ready"
            or row.generation != turn.generation
            or row.instance_id != turn.instance_id
            or row.bridge_owner_epoch != turn.bridge_owner_epoch
            or row.originating_auth_session_ref != turn.originating_auth_session_ref
        ):
            return PROPOSAL_INVALIDATED
        if row.draft_revision != turn.base_revision:
            return STALE_DRAFT
        if (
            row.schema_digest != turn.schema_digest
            or row.projection_digest != turn.projection_digest
        ):
            return CAPABILITY_CHANGED
        return None

    def _proposal_binding_ok(
        self,
        session: Session,
        row: AgentAuthoringSessionRow,
        turn: AgentAuthoringTurnRow,
        actor: AuthoringActor,
        now: datetime,
    ) -> bool:
        if actor.session_digest != turn.originating_auth_session_ref:
            return False
        if not _auth_live(session, actor.session_digest, actor.user_id, now):
            return False
        if (
            row.schema_digest != turn.schema_digest
            or row.projection_digest != turn.projection_digest
        ):
            return False
        return (
            row.connection_state == "ready"
            and row.generation == turn.generation
            and row.instance_id == turn.instance_id
            and row.bridge_owner_epoch == turn.bridge_owner_epoch
            and row.originating_auth_session_ref == turn.originating_auth_session_ref
        )

    def _owned(
        self,
        session: Session,
        actor: AuthoringActor,
        sid: str,
        now: datetime,
        *,
        allow_expired: bool,
    ) -> AgentAuthoringSessionRow:
        self._require_actor(session, actor, now)
        row = session.get(AgentAuthoringSessionRow, sid)
        if row is None or row.owner_user_id != actor.user_id:
            raise AuthoringRejected(NOT_FOUND)
        if not _project_visible(session, row.project_id):
            raise AuthoringRejected(NOT_FOUND)
        if _as_utc(row.purge_after) <= now:
            raise AuthoringRejected(NOT_FOUND)
        if not allow_expired and _as_utc(row.draft_expires_at) <= now:
            raise AuthoringRejected(DRAFT_EXPIRED)
        return row

    def _require_actor(
        self,
        session: Session,
        actor: AuthoringActor,
        now: datetime,
    ) -> None:
        if not _auth_live(session, actor.session_digest, actor.user_id, now):
            raise AuthoringRejected(AUTHENTICATION_REQUIRED)

    def _purge_due(
        self,
        session: Session,
        owner_user_id: str,
        now: datetime,
        *,
        limit: int,
    ) -> None:
        rows = session.scalars(
            select(AgentAuthoringSessionRow.sid)
            .where(
                AgentAuthoringSessionRow.owner_user_id == owner_user_id,
                AgentAuthoringSessionRow.purge_after <= now,
            )
            .limit(limit)
        ).all()
        for sid in rows:
            _delete_session(session, sid)

    def _now(self) -> datetime:
        return utc_now(self._clock)

    def _hook(self) -> None:
        if self._before_commit is not None:
            self._before_commit()


def fence_authoring_connections_for_auth_sessions(
    session: Session,
    session_digests: Sequence[str],
    revoked_at: datetime,
) -> None:
    """Close A1 connections inside the caller's authentication transaction."""

    digests = tuple(digest for digest in session_digests if digest)
    if not digests:
        return
    del revoked_at
    rows = session.scalars(
        select(AgentAuthoringSessionRow).where(
            AgentAuthoringSessionRow.originating_auth_session_ref.in_(digests),
            AgentAuthoringSessionRow.connection_state.in_(("blocked", "ready")),
        )
    ).all()
    for row in rows:
        _fence_generation(
            session,
            row.sid,
            generation=row.generation,
            owner_epoch=row.bridge_owner_epoch,
            reason_code=AUTH_SESSION_ENDED,
        )


def _fence_generation(
    session: Session,
    sid: str,
    *,
    generation: int,
    owner_epoch: str | None,
    reason_code: str,
) -> bool:
    result = session.execute(
        update(AgentAuthoringSessionRow)
        .where(
            AgentAuthoringSessionRow.sid == sid,
            AgentAuthoringSessionRow.generation == generation,
            AgentAuthoringSessionRow.connection_state.in_(("blocked", "ready")),
            AgentAuthoringSessionRow.bridge_owner_epoch == owner_epoch,
        )
        .values(
            connection_state="closed",
            connection_reason_code=reason_code,
            generation=generation + 1,
            connection_capability_digest=None,
        )
    )
    if result.rowcount != 1:
        return False
    session.execute(
        update(AgentAuthoringTurnRow)
        .where(
            AgentAuthoringTurnRow.sid == sid,
            AgentAuthoringTurnRow.generation == generation,
            AgentAuthoringTurnRow.status == "generating",
        )
        .values(status="failed", failure_reason=reason_code)
    )
    _invalidate_pending(session, sid, reason_code, generation=generation)
    return True


def _invalidate_pending(
    session: Session,
    sid: str,
    reason_code: str,
    *,
    generation: int | None = None,
    except_turn_id: str | None = None,
) -> None:
    statement = update(AgentAuthoringTurnRow).where(
        AgentAuthoringTurnRow.sid == sid,
        AgentAuthoringTurnRow.decision == "pending",
    )
    if generation is not None:
        statement = statement.where(AgentAuthoringTurnRow.generation == generation)
    if except_turn_id is not None:
        statement = statement.where(AgentAuthoringTurnRow.turn_id != except_turn_id)
    session.execute(
        statement.values(decision="invalidated", invalidation_reason=reason_code)
    )


def _delete_session(session: Session, sid: str) -> None:
    session.execute(
        delete(AgentAuthoringTurnRow).where(AgentAuthoringTurnRow.sid == sid)
    )
    session.execute(
        delete(AgentAuthoringSessionRow).where(AgentAuthoringSessionRow.sid == sid)
    )


def _clear_validation(row: AgentAuthoringSessionRow) -> None:
    row.last_snapshot_id = None
    row.last_validated_revision = None
    row.validation_issue_codes = None
    row.validation_issues = None


def _project_visible(session: Session, project_id: str) -> bool:
    row = session.get(ProjectRow, project_id)
    return (
        row is not None
        and row.kind == "user"
        and row.archived_at is None
        and row.project_id != LEGACY_PROJECT_ID
    )


def _auth_live(session: Session, digest: str, user_id: str, now: datetime) -> bool:
    row = session.get(AuthSessionRow, digest)
    if row is None or row.user_id != user_id or row.revoked_at is not None:
        return False
    if not (_as_utc(row.created_at) <= now < _as_utc(row.expires_at)):
        return False
    account = session.get(UserAccountRow, user_id)
    return account is not None and account.status == "enabled"


def _begin_write(session: Session) -> None:
    if session.get_bind().dialect.name == "sqlite":
        session.connection().exec_driver_sql("BEGIN IMMEDIATE")


def encode_authoring_cursor(updated_at: datetime, sid: str) -> str:
    raw = f"{_as_utc(updated_at).strftime('%Y-%m-%dT%H:%M:%SZ')}\n{sid}".encode()
    return base64.urlsafe_b64encode(raw).decode("ascii").rstrip("=")


def decode_authoring_cursor(value: str) -> tuple[datetime, str]:
    if not isinstance(value, str) or len(value) > 256:
        raise AuthoringRejected(INVALID_REQUEST)
    padded = value + ("=" * (-len(value) % 4))
    try:
        raw = base64.urlsafe_b64decode(padded.encode("ascii")).decode("utf-8")
        stamp, sid = raw.split("\n", 1)
        parsed = datetime.strptime(stamp, "%Y-%m-%dT%H:%M:%SZ").replace(
            tzinfo=timezone.utc
        )
    except (ValueError, UnicodeError):
        raise AuthoringRejected(INVALID_REQUEST) from None
    if len(sid) != 36 or not sid.startswith("aas_"):
        raise AuthoringRejected(INVALID_REQUEST)
    return parsed, sid


def _as_utc(value: datetime) -> datetime:
    if value.tzinfo is None or value.utcoffset() is None:
        value = value.replace(tzinfo=timezone.utc)
    return value.astimezone(timezone.utc).replace(microsecond=0)


def _session_record(row: AgentAuthoringSessionRow) -> AuthoringSessionRecord:
    codes = tuple(row.validation_issue_codes or ())
    issues = tuple(dict(issue) for issue in (row.validation_issues or ()))
    return AuthoringSessionRecord(
        sid=row.sid,
        owner_user_id=row.owner_user_id,
        project_id=row.project_id,
        workflow_id=row.workflow_id,
        schema_digest=row.schema_digest,
        projection_digest=row.projection_digest,
        content=dict(row.draft_json),
        draft_digest=row.draft_digest,
        revision=row.draft_revision,
        created_at=_as_utc(row.created_at),
        content_updated_at=_as_utc(row.content_updated_at),
        draft_expires_at=_as_utc(row.draft_expires_at),
        purge_after=_as_utc(row.purge_after),
        last_snapshot_id=row.last_snapshot_id,
        last_validated_revision=row.last_validated_revision,
        validation_issue_codes=tuple(str(code) for code in codes),
        validation_issues=issues,
        connection_state=row.connection_state,
        connection_reason_code=row.connection_reason_code,
        originating_auth_session_ref=row.originating_auth_session_ref,
        instance_id=row.instance_id,
        generation=row.generation,
        bridge_owner_epoch=row.bridge_owner_epoch,
        connection_capability_digest=row.connection_capability_digest,
    )


def _turn_record(row: AgentAuthoringTurnRow) -> AuthoringTurnRecord:
    proposal = None if row.proposal_json is None else dict(row.proposal_json)
    return AuthoringTurnRecord(
        turn_id=row.turn_id,
        sid=row.sid,
        client_turn_id=row.client_turn_id,
        request_digest=row.request_digest,
        question=dict(row.question_json),
        originating_auth_session_ref=row.originating_auth_session_ref,
        instance_id=row.instance_id,
        generation=row.generation,
        bridge_owner_epoch=row.bridge_owner_epoch,
        base_revision=row.base_revision,
        schema_digest=row.schema_digest,
        projection_digest=row.projection_digest,
        status=row.status,
        failure_reason=row.failure_reason,
        completed_at=None if row.completed_at is None else _as_utc(row.completed_at),
        result_kind=row.result_kind,
        result_digest=row.result_digest,
        proposal_id=row.proposal_id,
        proposal=proposal,
        proposal_digest=row.proposal_digest,
        proposal_expires_at=(
            None
            if row.proposal_expires_at is None
            else _as_utc(row.proposal_expires_at)
        ),
        decision=row.decision,
        invalidation_reason=row.invalidation_reason,
        adopt_request_digest=row.adopt_request_digest,
        adoption_actor_user_id=row.adoption_actor_user_id,
        adoption_auth_session_ref=row.adoption_auth_session_ref,
        adopted_at=None if row.adopted_at is None else _as_utc(row.adopted_at),
        result_revision=row.result_revision,
        adoption_changed=row.adoption_changed,
        adoption_content_digest=row.adoption_content_digest,
        adoption_draft_expires_at=(
            None
            if row.adoption_draft_expires_at is None
            else _as_utc(row.adoption_draft_expires_at)
        ),
    )
