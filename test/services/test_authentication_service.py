"""LAN authentication and account administration service behavior."""

from __future__ import annotations

from datetime import datetime, timedelta, timezone

import pytest

from encode_pipeline.persistence import (
    SqlAlchemyAuthenticationRepository,
    create_database_engine,
    create_session_factory,
    upgrade_database,
)
from encode_pipeline.platform.authentication import (
    SessionRevocationReason,
    UserRole,
)
from encode_pipeline.platform.security_audit import (
    AuditAction,
    AuditActorKind,
    AuditOutcome,
    AuditReasonCode,
)
from encode_pipeline.services.authentication import (
    Argon2idParameters,
    BoundedLoginRateLimiter,
    LoginRateLimitPolicy,
    PasswordManager,
    digest_session_token,
)
from encode_pipeline.services.authentication_service import (
    AuthenticationActor,
    AuthenticationError,
    AccountAdministrationService,
    AuthenticationService,
)

NOW = datetime(2026, 8, 18, 12, 0, tzinfo=timezone.utc)
FAST_PROFILE = Argon2idParameters(memory_cost_kib=19456, time_cost=2, parallelism=1)
LEGACY_PROFILE = Argon2idParameters(
    memory_cost_kib=19456,
    time_cost=2,
    parallelism=1,
    hash_length=33,
)


class _Clock:
    def __init__(self) -> None:
        self.now = NOW

    def __call__(self) -> datetime:
        return self.now


@pytest.fixture
def clock() -> _Clock:
    return _Clock()


@pytest.fixture
def repository(tmp_path):
    database_url = f"sqlite:///{tmp_path / 'authentication.db'}"
    upgrade_database(database_url)
    engine = create_database_engine(database_url)
    yield SqlAlchemyAuthenticationRepository(create_session_factory(engine))
    engine.dispose()


@pytest.fixture
def password_manager() -> PasswordManager:
    return PasswordManager(FAST_PROFILE, legacy_parameters=(LEGACY_PROFILE,))


@pytest.fixture
def auth_service(repository, password_manager, clock) -> AuthenticationService:
    return AuthenticationService(
        repository=repository,
        password_manager=password_manager,
        now_factory=clock,
    )


@pytest.fixture
def admin_service(repository, password_manager, clock) -> AccountAdministrationService:
    return AccountAdministrationService(
        repository=repository,
        password_manager=password_manager,
        now_factory=clock,
    )


@pytest.fixture
def administrator(admin_service):
    return admin_service.bootstrap_initial_administrator(
        "root-admin", "correct horse battery staple"
    )


def _reject(callable_, *args, code: str, **kwargs) -> None:
    with pytest.raises(AuthenticationError) as captured:
        callable_(*args, **kwargs)
    assert captured.value.reason_code == code


def _audit_events(repository):
    return repository.list_security_audit_events(limit=100)


def test_bootstrap_creates_the_unique_administrator(
    repository, admin_service, administrator
) -> None:
    assert administrator.role is UserRole.ADMINISTRATOR
    assert repository.get_account_by_username("root-admin") == administrator

    _reject(
        admin_service.bootstrap_initial_administrator,
        "second-admin",
        "another long password",
        code="OPERATION_CONFLICT",
    )

    (event,) = _audit_events(repository)
    assert event.action is AuditAction.ACCOUNT_CREATE
    assert event.actor_kind is AuditActorKind.LOCAL_OPERATOR
    assert event.actor_user_id is None
    assert event.resource is not None


def test_login_creates_a_session_and_writes_the_audit(
    repository, auth_service, administrator
) -> None:
    login = auth_service.login(
        "root-admin", "correct horse battery staple", client_identity="127.0.0.1"
    )

    assert login.principal.user_id == administrator.user_id
    digest = digest_session_token(login.secrets.session_token)
    record = repository.get_session(digest)
    assert record.user_id == administrator.user_id
    assert record.active_at(NOW) is True

    logins = [
        event
        for event in _audit_events(repository)
        if event.action is AuditAction.LOGIN
    ]
    assert len(logins) == 1
    event = logins[0]
    assert event.outcome is AuditOutcome.SUCCEEDED
    assert event.actor_kind is AuditActorKind.USER
    assert event.actor_user_id == administrator.user_id


def test_login_failures_are_uniform_and_audited_without_identity(
    repository, auth_service, administrator
) -> None:
    for username, password in (
        ("root-admin", "wrong password entirely"),
        ("unknown-user", "correct horse battery staple"),
    ):
        _reject(
            auth_service.login,
            username,
            password,
            client_identity="127.0.0.1",
            code="INVALID_CREDENTIALS",
        )

    failures = [
        event
        for event in _audit_events(repository)
        if event.action is AuditAction.LOGIN and event.outcome is AuditOutcome.FAILED
    ]
    assert len(failures) == 2
    for event in failures:
        assert event.actor_kind is AuditActorKind.UNAUTHENTICATED
        assert event.actor_user_id is None
        assert event.resource is None
        assert event.reason_code is AuditReasonCode.INVALID_CREDENTIALS


def test_login_rejects_a_disabled_account_like_an_unknown_one(
    repository, auth_service, admin_service, administrator
) -> None:
    actor = AuthenticationActor.local_operator()
    member = admin_service.create_member(actor, "lab-member", "member password phrase")
    admin_service.set_account_status(actor, member.user_id, enabled=False)
    _reject(
        auth_service.login,
        "lab-member",
        "member password phrase",
        client_identity="127.0.0.1",
        code="INVALID_CREDENTIALS",
    )


def test_login_rate_limiter_closes_the_subject_window(
    repository, password_manager, clock, administrator
) -> None:
    limiter = BoundedLoginRateLimiter(
        LoginRateLimitPolicy(
            subject_attempts=2,
            client_attempts=10,
            window_seconds=300.0,
        ),
        clock=lambda: 1000.0,
    )
    service = AuthenticationService(
        repository=repository,
        password_manager=password_manager,
        rate_limiter=limiter,
        now_factory=clock,
    )
    for _attempt in range(2):
        _reject(
            service.login,
            "root-admin",
            "wrong password entirely",
            client_identity="127.0.0.1",
            code="INVALID_CREDENTIALS",
        )
    _reject(
        service.login,
        "root-admin",
        "correct horse battery staple",
        client_identity="127.0.0.1",
        code="LOGIN_RATE_LIMITED",
    )

    limited = [
        event
        for event in _audit_events(repository)
        if event.reason_code is AuditReasonCode.LOGIN_RATE_LIMITED
    ]
    assert len(limited) == 1
    assert limited[0].actor_kind is AuditActorKind.UNAUTHENTICATED


def test_login_upgrades_a_legacy_password_hash_atomically(
    repository, password_manager, clock
) -> None:
    legacy_manager = PasswordManager(LEGACY_PROFILE)
    admin_service = AccountAdministrationService(
        repository=repository,
        password_manager=legacy_manager,
        now_factory=clock,
    )
    account = admin_service.bootstrap_initial_administrator(
        "legacy-admin", "correct horse battery staple"
    )
    legacy_hash = account.password_hash

    service = AuthenticationService(
        repository=repository,
        password_manager=password_manager,
        now_factory=clock,
    )
    clock.now = NOW + timedelta(minutes=5)
    login = service.login(
        "legacy-admin", "correct horse battery staple", client_identity="127.0.0.1"
    )

    upgraded = repository.get_account_by_id(account.user_id)
    assert upgraded.password_hash != legacy_hash
    assert upgraded.password_changed_at == account.password_changed_at
    assert upgraded.updated_at == NOW + timedelta(minutes=5)
    assert (
        repository.get_session(
            digest_session_token(login.secrets.session_token)
        ).user_id
        == account.user_id
    )

    # The upgraded hash still authenticates through the current profile only.
    current_only = AuthenticationService(
        repository=repository,
        password_manager=PasswordManager(FAST_PROFILE),
        now_factory=clock,
    )
    assert (
        current_only.login(
            "legacy-admin",
            "correct horse battery staple",
            client_identity="127.0.0.2",
        ).principal.user_id
        == account.user_id
    )


def test_logout_revokes_one_session_once(repository, auth_service, administrator):
    login = auth_service.login(
        "root-admin", "correct horse battery staple", client_identity="127.0.0.1"
    )
    assert auth_service.logout(login.secrets.session_token) is True

    record = repository.get_session(digest_session_token(login.secrets.session_token))
    assert record.revocation_reason is SessionRevocationReason.LOGOUT
    assert auth_service.logout(login.secrets.session_token) is False
    assert auth_service.logout("A" * 43) is False
    assert auth_service.logout("not-a-token") is False

    logouts = [
        event
        for event in _audit_events(repository)
        if event.action is AuditAction.LOGOUT
    ]
    assert len(logouts) == 1
    assert logouts[0].actor_user_id == administrator.user_id


def test_resolve_session_enforces_expiry_revocation_and_disablement(
    repository, auth_service, admin_service, administrator, clock
) -> None:
    login = auth_service.login(
        "root-admin", "correct horse battery staple", client_identity="127.0.0.1"
    )
    token = login.secrets.session_token
    assert auth_service.resolve_session(token) is not None

    clock.now = NOW + timedelta(hours=9)
    assert auth_service.resolve_session(token) is None

    clock.now = NOW
    actor = AuthenticationActor.local_operator()
    member = admin_service.create_member(actor, "lab-member", "member password phrase")
    second = auth_service.login(
        "lab-member", "member password phrase", client_identity="127.0.0.1"
    )
    assert auth_service.resolve_session(second.secrets.session_token) is not None
    admin_service.set_account_status(actor, member.user_id, enabled=False)
    assert auth_service.resolve_session(second.secrets.session_token) is None
    assert auth_service.resolve_session("not-a-token") is None


def test_session_csrf_requires_all_three_values_to_agree(
    repository, auth_service, administrator, clock
) -> None:
    login = auth_service.login(
        "root-admin", "correct horse battery staple", client_identity="127.0.0.1"
    )
    session_token = login.secrets.session_token
    csrf_token = login.secrets.csrf_token

    assert auth_service.session_csrf_valid(session_token, csrf_token, csrf_token)
    assert not auth_service.session_csrf_valid(session_token, csrf_token, "B" * 43)
    assert not auth_service.session_csrf_valid(session_token, "B" * 43, csrf_token)
    assert not auth_service.session_csrf_valid("not-a-token", csrf_token, csrf_token)

    clock.now = NOW + timedelta(hours=9)
    assert not auth_service.session_csrf_valid(session_token, csrf_token, csrf_token)


def test_create_member_requires_an_administrator(
    repository, admin_service, administrator
) -> None:
    admin_actor = AuthenticationActor.for_principal(
        auth_principal(repository, administrator.user_id)
    )
    member = admin_service.create_member(
        admin_actor, "lab-member", "member password phrase"
    )
    assert member.role is UserRole.MEMBER

    member_actor = AuthenticationActor.for_principal(
        auth_principal(repository, member.user_id)
    )
    _reject(
        admin_service.create_member,
        member_actor,
        "another-member",
        "member password phrase",
        code="ADMINISTRATOR_REQUIRED",
    )

    created = [
        event
        for event in _audit_events(repository)
        if event.action is AuditAction.ACCOUNT_CREATE
    ]
    assert any(
        event.actor_kind is AuditActorKind.USER
        and event.actor_user_id == administrator.user_id
        for event in created
    )


def auth_principal(repository, user_id: str):
    from encode_pipeline.platform.authentication import AuthenticatedPrincipal

    return AuthenticatedPrincipal.from_account(repository.get_account_by_id(user_id))


def test_create_member_rejects_a_duplicate_username(
    repository, admin_service, administrator
) -> None:
    actor = AuthenticationActor.local_operator()
    admin_service.create_member(actor, "lab-member", "member password phrase")
    _reject(
        admin_service.create_member,
        actor,
        " Lab-Member ",
        "member password phrase",
        code="OPERATION_CONFLICT",
    )


def test_disable_revokes_sessions_and_the_last_administrator_is_protected(
    repository, auth_service, admin_service, administrator, clock
) -> None:
    actor = AuthenticationActor.local_operator()
    member = admin_service.create_member(actor, "lab-member", "member password phrase")
    login = auth_service.login(
        "lab-member", "member password phrase", client_identity="127.0.0.1"
    )
    _reject(
        admin_service.set_account_status,
        actor,
        administrator.user_id,
        enabled=False,
        code="OPERATION_CONFLICT",
    )

    updated = admin_service.set_account_status(actor, member.user_id, enabled=False)
    assert updated.enabled is False
    record = repository.get_session(digest_session_token(login.secrets.session_token))
    assert record.revocation_reason is SessionRevocationReason.ACCOUNT_DISABLED

    re_enabled = admin_service.set_account_status(actor, member.user_id, enabled=True)
    assert re_enabled.enabled is True

    actions = [
        event.action
        for event in _audit_events(repository)
        if event.action in {AuditAction.ACCOUNT_DISABLE, AuditAction.ACCOUNT_ENABLE}
    ]
    assert AuditAction.ACCOUNT_DISABLE in actions
    assert AuditAction.ACCOUNT_ENABLE in actions


def test_reset_password_rotates_the_hash_and_revokes_sessions(
    repository, auth_service, admin_service, administrator, clock
) -> None:
    actor = AuthenticationActor.local_operator()
    login = auth_service.login(
        "root-admin", "correct horse battery staple", client_identity="127.0.0.1"
    )
    clock.now = NOW + timedelta(minutes=10)
    updated = admin_service.reset_password(
        actor, administrator.user_id, "a brand new password phrase"
    )
    assert updated.password_changed_at == NOW + timedelta(minutes=10)

    record = repository.get_session(digest_session_token(login.secrets.session_token))
    assert record.revocation_reason is SessionRevocationReason.PASSWORD_RESET

    _reject(
        auth_service.login,
        "root-admin",
        "correct horse battery staple",
        client_identity="127.0.0.1",
        code="INVALID_CREDENTIALS",
    )
    assert (
        auth_service.login(
            "root-admin", "a brand new password phrase", client_identity="127.0.0.1"
        ).principal.user_id
        == administrator.user_id
    )

    resets = [
        event
        for event in _audit_events(repository)
        if event.action is AuditAction.ACCOUNT_PASSWORD_RESET
    ]
    assert len(resets) == 1
    assert resets[0].actor_kind is AuditActorKind.LOCAL_OPERATOR


def test_revoke_sessions_closes_every_active_session(
    repository, auth_service, admin_service, administrator
) -> None:
    actor = AuthenticationActor.local_operator()
    for client in ("127.0.0.1", "127.0.0.2"):
        auth_service.login(
            "root-admin", "correct horse battery staple", client_identity=client
        )
    revoked = admin_service.revoke_sessions(actor, administrator.user_id)
    assert revoked == 2
    assert admin_service.revoke_sessions(actor, administrator.user_id) == 0

    _reject(
        admin_service.revoke_sessions,
        actor,
        "usr_" + "9" * 32,
        code="RESOURCE_NOT_FOUND",
    )
