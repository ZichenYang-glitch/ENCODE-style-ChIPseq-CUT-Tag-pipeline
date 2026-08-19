"""LAN login lifecycle and administrator account operations."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, replace
from datetime import datetime, timedelta, timezone

from encode_pipeline.platform.authentication import (
    AuthenticatedPrincipal,
    SessionRecord,
    SessionRevocationReason,
    UserAccount,
    UserRole,
    UserStatus,
    authenticated_principal_for_session,
)
from encode_pipeline.platform.security_audit import (
    AuditAction,
    AuditActorKind,
    AuditOutcome,
    AuditReasonCode,
    AuditResourceKind,
    SecurityAuditEvent,
    build_audit_resource,
    new_audit_event_id,
)
from encode_pipeline.services.authentication import (
    DEFAULT_SESSION_LIFETIME,
    BoundedLoginRateLimiter,
    PasswordManager,
    SessionSecrets,
    authenticate_login_candidate,
    csrf_request_is_valid,
    digest_session_token,
    new_session_record,
    new_session_secrets,
    new_user_id,
)
from encode_pipeline.services.authentication_repositories import (
    AuthenticationAccountConflictError,
    AuthenticationRepository,
)


AUTHENTICATION_ISSUE_MESSAGES = {
    "INVALID_CREDENTIALS": "The username or password is invalid.",
    "LOGIN_RATE_LIMITED": "The username or password is invalid.",
    "ADMINISTRATOR_REQUIRED": "The operation requires an administrator.",
    "OPERATION_CONFLICT": "The operation conflicts with the current account state.",
    "RESOURCE_NOT_FOUND": "The account does not exist.",
}


class AuthenticationError(RuntimeError):
    """Stable redacted authentication and account error."""

    def __init__(self, reason_code: str) -> None:
        self.reason_code = reason_code
        super().__init__(AUTHENTICATION_ISSUE_MESSAGES[reason_code])


@dataclass(frozen=True)
class AuthenticationActor:
    """The auditable caller of one account operation."""

    kind: AuditActorKind
    user_id: str | None = None

    @classmethod
    def for_principal(cls, principal: AuthenticatedPrincipal) -> AuthenticationActor:
        if not isinstance(principal, AuthenticatedPrincipal):
            raise ValueError("principal must be an AuthenticatedPrincipal")
        return cls(AuditActorKind.USER, principal.user_id)

    @classmethod
    def local_operator(cls) -> AuthenticationActor:
        return cls(AuditActorKind.LOCAL_OPERATOR, None)

    def __post_init__(self) -> None:
        if not isinstance(self.kind, AuditActorKind):
            raise ValueError("actor kind is invalid")
        if (self.kind is AuditActorKind.USER) != (self.user_id is not None):
            raise ValueError("only a user actor may carry actor_user_id")


@dataclass(frozen=True, repr=False)
class LoginSession:
    """One authenticated principal with its fresh browser credentials."""

    principal: AuthenticatedPrincipal
    secrets: SessionSecrets
    record: SessionRecord


class AuthenticationService:
    """Login, logout, and session resolution over the SQLite session authority."""

    def __init__(
        self,
        *,
        repository: AuthenticationRepository,
        password_manager: PasswordManager | None = None,
        rate_limiter: BoundedLoginRateLimiter | None = None,
        session_lifetime: timedelta = DEFAULT_SESSION_LIFETIME,
        audit_event_id_factory: Callable[[], str] = new_audit_event_id,
        now_factory: Callable[[], datetime] | None = None,
    ) -> None:
        self._repository = repository
        self._password_manager = password_manager or PasswordManager()
        self._rate_limiter = rate_limiter or BoundedLoginRateLimiter()
        self._session_lifetime = session_lifetime
        self._audit_event_id_factory = audit_event_id_factory
        self._now_factory = now_factory or (lambda: datetime.now(timezone.utc))

    def login(
        self,
        username: object,
        password: object,
        *,
        client_identity: str,
    ) -> LoginSession:
        """Verify one bounded login candidate and create its session atomically."""

        if not self._rate_limiter.allow_attempt(client_identity, username):
            self._record_login_failure(AuditReasonCode.LOGIN_RATE_LIMITED)
            raise AuthenticationError("LOGIN_RATE_LIMITED")
        try:
            account = self._repository.get_account_by_username(username)
        except (KeyError, ValueError):
            account = None
        authentication = authenticate_login_candidate(
            account,
            password,
            self._password_manager,
        )
        if not authentication.authenticated:
            self._record_login_failure(AuditReasonCode.INVALID_CREDENTIALS)
            raise AuthenticationError("INVALID_CREDENTIALS")
        principal = authentication.principal
        if principal is None or account is None:
            raise RuntimeError("authenticated login lost its principal")
        now = self._now()
        secrets = new_session_secrets()
        record = new_session_record(
            user_id=principal.user_id,
            secrets=secrets,
            created_at=now,
            lifetime=self._session_lifetime,
        )
        rehashed = (
            replace(
                account, password_hash=authentication.replacement_hash, updated_at=now
            )
            if authentication.replacement_hash is not None
            else None
        )
        self._repository.create_session(
            record,
            updated_account=rehashed,
            audit=self._event(
                AuditAction.LOGIN,
                AuditOutcome.SUCCEEDED,
                AuditActorKind.USER,
                actor_user_id=principal.user_id,
                occurred_at=now,
            ),
        )
        return LoginSession(principal=principal, secrets=secrets, record=record)

    def logout(self, session_token: object) -> bool:
        """Revoke one known session before cookies are cleared."""

        record = self._record_for_token(session_token)
        if record is None or record.revoked_at is not None:
            return False
        return self._repository.revoke_session(
            record.session_digest,
            self._now(),
            SessionRevocationReason.LOGOUT,
            audit=self._event(
                AuditAction.LOGOUT,
                AuditOutcome.SUCCEEDED,
                AuditActorKind.USER,
                actor_user_id=record.user_id,
                occurred_at=self._now(),
            ),
        )

    def setup_complete(self) -> bool:
        """Return whether any enabled administrator exists in the deployment."""

        return self._repository.has_enabled_administrator()

    def resolve_session(self, session_token: object) -> AuthenticatedPrincipal | None:
        """Resolve one browser session to its current enabled principal."""

        record = self._record_for_token(session_token)
        if record is None:
            return None
        try:
            account = self._repository.get_account_by_id(record.user_id)
        except KeyError:
            return None
        return authenticated_principal_for_session(record, account, self._now())

    def session_csrf_valid(
        self,
        session_token: object,
        cookie_token: object,
        header_token: object,
    ) -> bool:
        """Require the readable cookie, header, and durable digest to agree."""

        record = self._record_for_token(session_token)
        if record is None or not record.active_at(self._now()):
            return False
        return csrf_request_is_valid(cookie_token, header_token, record.csrf_digest)

    def _record_for_token(self, session_token: object) -> SessionRecord | None:
        try:
            digest = digest_session_token(session_token)
        except ValueError:
            return None
        try:
            return self._repository.get_session(digest)
        except KeyError:
            return None

    def _record_login_failure(self, reason: AuditReasonCode) -> None:
        self._repository.record_security_audit(
            self._event(
                AuditAction.LOGIN,
                AuditOutcome.FAILED,
                AuditActorKind.UNAUTHENTICATED,
                reason_code=reason,
                occurred_at=self._now(),
            )
        )

    def _event(
        self,
        action: AuditAction,
        outcome: AuditOutcome,
        actor_kind: AuditActorKind,
        *,
        actor_user_id: str | None = None,
        reason_code: AuditReasonCode | None = None,
        occurred_at: datetime,
    ) -> SecurityAuditEvent:
        return SecurityAuditEvent(
            event_id=self._audit_event_id_factory(),
            occurred_at=occurred_at,
            action=action,
            outcome=outcome,
            actor_kind=actor_kind,
            actor_user_id=actor_user_id,
            reason_code=reason_code,
        )

    def _now(self) -> datetime:
        value = self._now_factory()
        if (
            not isinstance(value, datetime)
            or value.tzinfo is None
            or value.utcoffset() is None
        ):
            raise ValueError("now_factory must return timezone-aware datetimes")
        return value.astimezone(timezone.utc)


class AccountAdministrationService:
    """Administrator and local-operator account lifecycle operations."""

    def __init__(
        self,
        *,
        repository: AuthenticationRepository,
        password_manager: PasswordManager | None = None,
        audit_event_id_factory: Callable[[], str] = new_audit_event_id,
        user_id_factory: Callable[[], str] = new_user_id,
        now_factory: Callable[[], datetime] | None = None,
    ) -> None:
        self._repository = repository
        self._password_manager = password_manager or PasswordManager()
        self._audit_event_id_factory = audit_event_id_factory
        self._user_id_factory = user_id_factory
        self._now_factory = now_factory or (lambda: datetime.now(timezone.utc))

    def list_accounts(self) -> tuple[UserAccount, ...]:
        """Return safe account projections for administrator surfaces."""

        return self._repository.list_accounts()

    def bootstrap_initial_administrator(
        self,
        username: object,
        password: object,
    ) -> UserAccount:
        """Create the unique first administrator; never an HTTP operation."""

        if any(
            account.role is UserRole.ADMINISTRATOR
            for account in self._repository.list_accounts()
        ):
            raise AuthenticationError("OPERATION_CONFLICT")
        return self._create_account(
            username,
            password,
            role=UserRole.ADMINISTRATOR,
            actor=AuthenticationActor.local_operator(),
        )

    def create_member(
        self,
        actor: AuthenticationActor,
        username: object,
        password: object,
    ) -> UserAccount:
        self._require_administrator(actor)
        return self._create_account(
            username,
            password,
            role=UserRole.MEMBER,
            actor=actor,
        )

    def set_account_status(
        self,
        actor: AuthenticationActor,
        user_id: str,
        *,
        enabled: bool,
    ) -> UserAccount:
        self._require_administrator(actor)
        if not isinstance(enabled, bool):
            raise ValueError("enabled must be boolean")
        account = self._account(user_id)
        now = self._now()
        if not enabled and account.enabled and account.role is UserRole.ADMINISTRATOR:
            others = any(
                other.role is UserRole.ADMINISTRATOR
                and other.enabled
                and other.user_id != account.user_id
                for other in self._repository.list_accounts()
            )
            if not others:
                raise AuthenticationError("OPERATION_CONFLICT")
        updated = account.enable(now) if enabled else account.disable(now)
        self._repository.save_account(
            updated,
            audit=self._account_event(
                AuditAction.ACCOUNT_ENABLE if enabled else AuditAction.ACCOUNT_DISABLE,
                actor,
                account.user_id,
                occurred_at=now,
            ),
            revoke_sessions_reason=(
                None if enabled else SessionRevocationReason.ACCOUNT_DISABLED
            ),
            revoked_at=None if enabled else now,
        )
        return updated

    def reset_password(
        self,
        actor: AuthenticationActor,
        user_id: str,
        new_password: object,
    ) -> UserAccount:
        self._require_administrator(actor)
        account = self._account(user_id)
        now = self._now()
        updated = account.change_password(
            self._password_manager.hash_password(new_password),
            now,
        )
        self._repository.save_account(
            updated,
            audit=self._account_event(
                AuditAction.ACCOUNT_PASSWORD_RESET,
                actor,
                account.user_id,
                occurred_at=now,
            ),
            revoke_sessions_reason=SessionRevocationReason.PASSWORD_RESET,
            revoked_at=now,
        )
        return updated

    def reset_password_for_username(
        self,
        username: object,
        new_password: object,
    ) -> UserAccount:
        """Reset one account by login name as a local operator operation."""

        try:
            account = self._repository.get_account_by_username(username)
        except KeyError:
            raise AuthenticationError("RESOURCE_NOT_FOUND") from None
        return self.reset_password(
            AuthenticationActor.local_operator(),
            account.user_id,
            new_password,
        )

    def revoke_sessions(
        self,
        actor: AuthenticationActor,
        user_id: str,
    ) -> int:
        self._require_administrator(actor)
        account = self._account(user_id)
        now = self._now()
        return self._repository.save_account(
            account,
            audit=self._account_event(
                AuditAction.ACCOUNT_SESSIONS_REVOKE,
                actor,
                account.user_id,
                occurred_at=now,
            ),
            revoke_sessions_reason=SessionRevocationReason.ALL_SESSIONS,
            revoked_at=now,
        )

    def _create_account(
        self,
        username: object,
        password: object,
        *,
        role: UserRole,
        actor: AuthenticationActor,
    ) -> UserAccount:
        now = self._now()
        account = UserAccount(
            user_id=self._user_id_factory(),
            username=username,
            role=role,
            status=UserStatus.ENABLED,
            password_hash=self._password_manager.hash_password(password),
            created_at=now,
            updated_at=now,
            password_changed_at=now,
        )
        try:
            return self._repository.create_account(
                account,
                audit=self._account_event(
                    AuditAction.ACCOUNT_CREATE,
                    actor,
                    account.user_id,
                    occurred_at=now,
                ),
            )
        except AuthenticationAccountConflictError:
            raise AuthenticationError("OPERATION_CONFLICT") from None

    def _account(self, user_id: str) -> UserAccount:
        try:
            return self._repository.get_account_by_id(user_id)
        except KeyError:
            raise AuthenticationError("RESOURCE_NOT_FOUND") from None

    def _require_administrator(self, actor: AuthenticationActor) -> None:
        if not isinstance(actor, AuthenticationActor):
            raise ValueError("actor must be an AuthenticationActor")
        if actor.kind is AuditActorKind.LOCAL_OPERATOR:
            return
        try:
            account = self._repository.get_account_by_id(actor.user_id or "")
        except KeyError:
            raise AuthenticationError("ADMINISTRATOR_REQUIRED") from None
        if not account.enabled or account.role is not UserRole.ADMINISTRATOR:
            raise AuthenticationError("ADMINISTRATOR_REQUIRED")

    def _account_event(
        self,
        action: AuditAction,
        actor: AuthenticationActor,
        user_id: str,
        *,
        occurred_at: datetime,
    ) -> SecurityAuditEvent:
        return SecurityAuditEvent(
            event_id=self._audit_event_id_factory(),
            occurred_at=occurred_at,
            action=action,
            outcome=AuditOutcome.SUCCEEDED,
            actor_kind=actor.kind,
            actor_user_id=actor.user_id,
            resource=build_audit_resource(AuditResourceKind.ACCOUNT, user_id),
        )

    def _now(self) -> datetime:
        value = self._now_factory()
        if (
            not isinstance(value, datetime)
            or value.tzinfo is None
            or value.utcoffset() is None
        ):
            raise ValueError("now_factory must return timezone-aware datetimes")
        return value.astimezone(timezone.utc)
