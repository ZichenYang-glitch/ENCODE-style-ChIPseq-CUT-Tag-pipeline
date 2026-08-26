"""Persistence boundary for LAN accounts, sessions, and security audit rows."""

from __future__ import annotations

from datetime import datetime
from typing import Protocol

from encode_pipeline.platform.authentication import (
    SessionRecord,
    SessionRevocationReason,
    UserAccount,
)
from encode_pipeline.platform.security_audit import SecurityAuditEvent


class AuthenticationAccountConflictError(RuntimeError):
    """An account identity or unique login name invariant conflicted."""


class AuthenticationRepository(Protocol):
    """Atomic storage contract consumed by LAN authentication services."""

    def create_account(
        self,
        account: UserAccount,
        *,
        audit: SecurityAuditEvent | None = None,
    ) -> UserAccount: ...

    def get_account_by_id(self, user_id: str) -> UserAccount: ...

    def get_account_by_username(self, username: str) -> UserAccount: ...

    def list_accounts(self) -> tuple[UserAccount, ...]: ...

    def has_enabled_administrator(self) -> bool: ...

    def save_account(
        self,
        account: UserAccount,
        *,
        audit: SecurityAuditEvent | None = None,
        revoke_sessions_reason: SessionRevocationReason | None = None,
        revoked_at: datetime | None = None,
    ) -> int: ...

    def create_session(
        self,
        session: SessionRecord,
        *,
        updated_account: UserAccount | None = None,
        audit: SecurityAuditEvent | None = None,
    ) -> None: ...

    def record_security_audit(self, event: SecurityAuditEvent) -> None: ...

    def get_session(self, session_digest: str) -> SessionRecord: ...

    def revoke_session(
        self,
        session_digest: str,
        revoked_at: datetime,
        reason: SessionRevocationReason,
        *,
        audit: SecurityAuditEvent | None = None,
    ) -> bool: ...

    def list_security_audit_events(
        self,
        *,
        limit: int = 100,
    ) -> tuple[SecurityAuditEvent, ...]: ...
