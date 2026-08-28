"""SQLAlchemy persistence for LAN accounts, sessions, and security audit rows."""

from __future__ import annotations

from collections.abc import Callable
from datetime import datetime, timezone

from sqlalchemy import select
from sqlalchemy.exc import IntegrityError
from sqlalchemy.orm import Session, sessionmaker

from encode_pipeline.persistence.models import (
    AuthSessionRow,
    SecurityAuditEventRow,
    UserAccountRow,
)
from encode_pipeline.platform.authentication import (
    SessionRecord,
    SessionRevocationReason,
    UserAccount,
    UserRole,
    UserStatus,
    normalize_username,
)
from encode_pipeline.platform.security_audit import (
    AuditResource,
    AuditResourceKind,
    SecurityAuditEvent,
)
from encode_pipeline.services.authentication_repositories import (
    AuthenticationAccountConflictError,
)


class SqlAlchemyAuthenticationRepository:
    """Persist LAN authentication state without exposing ORM rows."""

    def __init__(self, session_factory: sessionmaker[Session]) -> None:
        if not isinstance(session_factory, sessionmaker):
            raise ValueError("session_factory must be a SQLAlchemy sessionmaker")
        self._session_factory = session_factory

    def create_account(
        self,
        account: UserAccount,
        *,
        audit: SecurityAuditEvent | None = None,
    ) -> UserAccount:
        _validate_account(account)
        _validate_audit(audit)
        try:
            with self._session_factory.begin() as session:
                _begin_write(session)
                session.add(_account_row(account))
                session.flush()
                if audit is not None:
                    insert_security_audit_event(session, audit)
        except IntegrityError as exc:
            raise AuthenticationAccountConflictError(
                "user account identity already exists"
            ) from exc
        return account

    def get_account_by_id(self, user_id: str) -> UserAccount:
        with self._session_factory() as session:
            row = session.get(UserAccountRow, user_id)
            if row is None:
                raise KeyError("unknown user account")
            return _account_from_row(row)

    def get_account_by_username(self, username: str) -> UserAccount:
        normalized = normalize_username(username)
        with self._session_factory() as session:
            row = session.scalar(
                select(UserAccountRow).where(UserAccountRow.username == normalized)
            )
            if row is None:
                raise KeyError("unknown user account")
            return _account_from_row(row)

    def list_accounts(self) -> tuple[UserAccount, ...]:
        with self._session_factory() as session:
            rows = session.scalars(
                select(UserAccountRow).order_by(UserAccountRow.username)
            ).all()
            return tuple(_account_from_row(row) for row in rows)

    def has_enabled_administrator(self) -> bool:
        with self._session_factory() as session:
            return (
                session.scalar(
                    select(UserAccountRow.user_id)
                    .where(
                        UserAccountRow.role == UserRole.ADMINISTRATOR.value,
                        UserAccountRow.status == UserStatus.ENABLED.value,
                    )
                    .limit(1)
                )
                is not None
            )

    def save_account(
        self,
        account: UserAccount,
        *,
        audit: SecurityAuditEvent | None = None,
        revoke_sessions_reason: SessionRevocationReason | None = None,
        revoked_at: datetime | None = None,
    ) -> int:
        _validate_account(account)
        _validate_audit(audit)
        if (revoke_sessions_reason is None) != (revoked_at is None):
            raise ValueError(
                "revoke_sessions_reason and revoked_at must be provided together"
            )
        reason = (
            _revocation_reason(revoke_sessions_reason)
            if revoke_sessions_reason is not None
            else None
        )
        revoked_time = _write_time(revoked_at) if revoked_at is not None else None
        try:
            with self._session_factory.begin() as session:
                _begin_write(session)
                row = session.get(UserAccountRow, account.user_id)
                if row is None:
                    raise KeyError("unknown user account")
                _apply_account_update(row, account)
                session.flush()
                revoked = 0
                if reason is not None:
                    active_rows = session.scalars(
                        select(AuthSessionRow).where(
                            AuthSessionRow.user_id == account.user_id,
                            AuthSessionRow.revoked_at.is_(None),
                        )
                    ).all()
                    for active in active_rows:
                        active.revoked_at = revoked_time
                        active.revocation_reason = reason.value
                    session.flush()
                    revoked = len(active_rows)
                if audit is not None:
                    insert_security_audit_event(session, audit)
        except IntegrityError as exc:
            raise AuthenticationAccountConflictError(
                "user account state conflicted"
            ) from exc
        return revoked

    def create_session(
        self,
        session: SessionRecord,
        *,
        updated_account: UserAccount | None = None,
        audit: SecurityAuditEvent | None = None,
    ) -> None:
        if not isinstance(session, SessionRecord):
            raise ValueError("session must be a SessionRecord")
        if updated_account is not None:
            _validate_account(updated_account)
        _validate_audit(audit)
        with self._session_factory.begin() as write_session:
            _begin_write(write_session)
            if updated_account is not None:
                row = write_session.get(UserAccountRow, updated_account.user_id)
                if row is None:
                    raise KeyError("unknown user account")
                _apply_account_update(row, updated_account)
                write_session.flush()
            write_session.add(_session_row(session))
            write_session.flush()
            if audit is not None:
                insert_security_audit_event(write_session, audit)

    def set_notification_email(
        self,
        user_id: str,
        notification_email: str | None,
        changed_at: datetime,
    ) -> UserAccount:
        return self._update_notification_account(
            user_id,
            lambda account: account.change_notification_email(
                notification_email,
                changed_at,
            ),
        )

    def set_terminal_email_enabled(
        self,
        user_id: str,
        enabled: bool,
        changed_at: datetime,
    ) -> UserAccount:
        return self._update_notification_account(
            user_id,
            lambda account: account.change_terminal_email_enabled(
                enabled,
                changed_at,
            ),
        )

    def _update_notification_account(
        self,
        user_id: str,
        update_account: Callable[[UserAccount], UserAccount],
    ) -> UserAccount:
        with self._session_factory.begin() as session:
            _begin_write(session)
            row = session.get(UserAccountRow, user_id)
            if row is None:
                raise KeyError("unknown user account")
            updated = update_account(_account_from_row(row))
            _apply_account_update(row, updated)
            session.flush()
            return updated

    def record_security_audit(self, event: SecurityAuditEvent) -> None:
        _validate_audit(event)
        with self._session_factory.begin() as session:
            _begin_write(session)
            insert_security_audit_event(session, event)

    def get_session(self, session_digest: str) -> SessionRecord:
        with self._session_factory() as session:
            row = session.get(AuthSessionRow, session_digest)
            if row is None:
                raise KeyError("unknown auth session")
            return _session_from_row(row)

    def revoke_session(
        self,
        session_digest: str,
        revoked_at: datetime,
        reason: SessionRevocationReason,
        *,
        audit: SecurityAuditEvent | None = None,
    ) -> bool:
        normalized_reason = _revocation_reason(reason)
        revoked_time = _write_time(revoked_at)
        _validate_audit(audit)
        with self._session_factory.begin() as session:
            _begin_write(session)
            row = session.get(AuthSessionRow, session_digest)
            if row is None:
                raise KeyError("unknown auth session")
            if row.revoked_at is not None:
                return False
            row.revoked_at = revoked_time
            row.revocation_reason = normalized_reason.value
            session.flush()
            if audit is not None:
                insert_security_audit_event(session, audit)
            return True

    def list_security_audit_events(
        self,
        *,
        limit: int = 100,
    ) -> tuple[SecurityAuditEvent, ...]:
        if (
            not isinstance(limit, int)
            or isinstance(limit, bool)
            or not 1 <= limit <= 1000
        ):
            raise ValueError("limit must be between 1 and 1000")
        with self._session_factory() as session:
            rows = session.scalars(
                select(SecurityAuditEventRow)
                .order_by(
                    SecurityAuditEventRow.occurred_at.desc(),
                    SecurityAuditEventRow.event_id.desc(),
                )
                .limit(limit)
            ).all()
            return tuple(_audit_from_row(row) for row in rows)


def insert_security_audit_event(
    session: Session,
    event: SecurityAuditEvent,
) -> None:
    """Append one closed audit row inside the caller's own transaction."""

    if not isinstance(event, SecurityAuditEvent):
        raise ValueError("event must be a SecurityAuditEvent")
    session.add(_audit_row(event))
    session.flush()


def _validate_account(account: UserAccount) -> None:
    if not isinstance(account, UserAccount):
        raise ValueError("account must be a UserAccount")


def _validate_audit(audit: SecurityAuditEvent | None) -> None:
    if audit is not None and not isinstance(audit, SecurityAuditEvent):
        raise ValueError("audit must be a SecurityAuditEvent")


def _revocation_reason(
    value: SessionRevocationReason | str,
) -> SessionRevocationReason:
    try:
        return (
            value
            if isinstance(value, SessionRevocationReason)
            else SessionRevocationReason(value)
        )
    except (TypeError, ValueError):
        raise ValueError("revocation reason is invalid") from None


def _write_time(value: datetime) -> datetime:
    if (
        not isinstance(value, datetime)
        or value.tzinfo is None
        or value.utcoffset() is None
    ):
        raise ValueError("revoked_at must be timezone-aware")
    return value.astimezone(timezone.utc)


def _apply_account_update(row: UserAccountRow, account: UserAccount) -> None:
    row.username = account.username
    row.role = account.role.value
    row.status = account.status.value
    row.password_hash = account.password_hash
    row.created_at = account.created_at
    row.updated_at = account.updated_at
    row.password_changed_at = account.password_changed_at
    row.notification_email = account.notification_email
    row.terminal_email_enabled = account.terminal_email_enabled


def _account_row(account: UserAccount) -> UserAccountRow:
    return UserAccountRow(
        user_id=account.user_id,
        username=account.username,
        role=account.role.value,
        status=account.status.value,
        password_hash=account.password_hash,
        created_at=account.created_at,
        updated_at=account.updated_at,
        password_changed_at=account.password_changed_at,
        notification_email=account.notification_email,
        terminal_email_enabled=account.terminal_email_enabled,
    )


def _account_from_row(row: UserAccountRow) -> UserAccount:
    return UserAccount(
        user_id=row.user_id,
        username=row.username,
        role=row.role,
        status=row.status,
        password_hash=row.password_hash,
        created_at=_as_utc(row.created_at),
        updated_at=_as_utc(row.updated_at),
        password_changed_at=_as_utc(row.password_changed_at),
        notification_email=row.notification_email,
        terminal_email_enabled=row.terminal_email_enabled,
    )


def _session_row(record: SessionRecord) -> AuthSessionRow:
    return AuthSessionRow(
        session_digest=record.session_digest,
        csrf_digest=record.csrf_digest,
        user_id=record.user_id,
        created_at=record.created_at,
        expires_at=record.expires_at,
        revoked_at=record.revoked_at,
        revocation_reason=(
            record.revocation_reason.value
            if record.revocation_reason is not None
            else None
        ),
    )


def _session_from_row(row: AuthSessionRow) -> SessionRecord:
    return SessionRecord(
        session_digest=row.session_digest,
        csrf_digest=row.csrf_digest,
        user_id=row.user_id,
        created_at=_as_utc(row.created_at),
        expires_at=_as_utc(row.expires_at),
        revoked_at=_as_utc(row.revoked_at) if row.revoked_at is not None else None,
        revocation_reason=row.revocation_reason,
    )


def _audit_row(event: SecurityAuditEvent) -> SecurityAuditEventRow:
    return SecurityAuditEventRow(
        event_id=event.event_id,
        occurred_at=event.occurred_at,
        action=event.action.value,
        outcome=event.outcome.value,
        actor_kind=event.actor_kind.value,
        actor_user_id=event.actor_user_id,
        resource_kind=(
            event.resource.kind.value if event.resource is not None else None
        ),
        resource_id=(
            event.resource.resource_id if event.resource is not None else None
        ),
        reason_code=(
            event.reason_code.value if event.reason_code is not None else None
        ),
    )


def _audit_from_row(row: SecurityAuditEventRow) -> SecurityAuditEvent:
    return SecurityAuditEvent(
        event_id=row.event_id,
        occurred_at=_as_utc(row.occurred_at),
        action=row.action,
        outcome=row.outcome,
        actor_kind=row.actor_kind,
        actor_user_id=row.actor_user_id,
        resource=(
            AuditResource._from_derived(
                AuditResourceKind(row.resource_kind),
                row.resource_id,
            )
            if row.resource_kind is not None
            else None
        ),
        reason_code=row.reason_code,
    )


def _as_utc(value: datetime) -> datetime:
    if value.tzinfo is None or value.utcoffset() is None:
        return value.replace(tzinfo=timezone.utc)
    return value.astimezone(timezone.utc)


def _begin_write(session: Session) -> None:
    if session.get_bind().dialect.name == "sqlite":
        session.connection().exec_driver_sql("BEGIN IMMEDIATE")
