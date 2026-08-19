"""Disclosure-safe audit event contracts for authenticated actions."""

from __future__ import annotations

from datetime import datetime, timezone

import pytest

from encode_pipeline.platform.security_audit import (
    AuditAction,
    AuditActorKind,
    AuditOutcome,
    AuditReasonCode,
    AuditResource,
    AuditResourceKind,
    SecurityAuditEvent,
    build_audit_resource,
    new_audit_event_id,
)


NOW = datetime(2026, 8, 9, 12, 0, tzinfo=timezone.utc)
USER_ID = "usr_11111111111111111111111111111111"


def test_audit_action_vocabulary_covers_only_current_and_auth_owned_actions() -> None:
    assert {action.value for action in AuditAction} == {
        "auth.login",
        "auth.logout",
        "run.create",
        "run.start",
        "run.cancel",
        "artifact.download",
        "reference.register",
        "reference.enable",
        "reference.disable",
        "storage.register",
        "storage.archive",
        "account.create",
        "account.enable",
        "account.disable",
        "account.password_reset",
        "account.sessions_revoke",
        "run.fail",
        "run.requeue",
    }


def test_audit_actor_outcome_and_resource_vocabularies_are_closed() -> None:
    assert {kind.value for kind in AuditActorKind} == {
        "user",
        "local_operator",
        "unauthenticated",
    }
    assert {outcome.value for outcome in AuditOutcome} == {
        "succeeded",
        "denied",
        "failed",
    }
    assert {reason.value for reason in AuditReasonCode} == {
        "INVALID_CREDENTIALS",
        "LOGIN_RATE_LIMITED",
        "AUTHENTICATION_REQUIRED",
        "SESSION_INVALID",
        "CSRF_INVALID",
        "ADMINISTRATOR_REQUIRED",
        "ACCOUNT_DISABLED",
        "SETUP_REQUIRED",
        "OPERATION_CONFLICT",
        "RESOURCE_NOT_FOUND",
        "INTERNAL_FAILURE",
    }
    assert {kind.value for kind in AuditResourceKind} == {
        "account",
        "run",
        "artifact",
        "reference",
        "storage",
    }


def test_audit_event_has_only_stable_identity_and_result_fields() -> None:
    resource = build_audit_resource(AuditResourceKind.RUN, "run-123")
    event = SecurityAuditEvent(
        event_id="aevt_11111111111111111111111111111111",
        occurred_at=NOW,
        action="run.start",
        outcome="succeeded",
        actor_kind="user",
        actor_user_id=USER_ID,
        resource=resource,
        reason_code=None,
    )

    assert event.to_dict() == {
        "event_id": "aevt_11111111111111111111111111111111",
        "occurred_at": NOW,
        "action": "run.start",
        "outcome": "succeeded",
        "actor_kind": "user",
        "actor_user_id": USER_ID,
        "resource": {"kind": "run", "resource_id": resource.resource_id},
        "reason_code": None,
    }
    serialized = str(event.to_dict()).lower()
    for forbidden in (
        "password",
        "token",
        "authorization",
        "cookie",
        "exception",
        "/home/",
        "@example.com",
        "run-123",
    ):
        assert forbidden not in serialized


def test_failed_login_summary_has_constant_non_identifying_shape() -> None:
    event = SecurityAuditEvent(
        event_id="aevt_22222222222222222222222222222222",
        occurred_at=NOW,
        action=AuditAction.LOGIN,
        outcome=AuditOutcome.DENIED,
        actor_kind=AuditActorKind.UNAUTHENTICATED,
        reason_code="INVALID_CREDENTIALS",
    )
    assert event.actor_user_id is None
    assert event.resource is None
    assert event.to_dict()["reason_code"] == "INVALID_CREDENTIALS"


def test_audit_action_shape_prevents_identity_and_authority_spoofing() -> None:
    with pytest.raises(ValueError, match="non-identifying"):
        SecurityAuditEvent(
            event_id="aevt_66666666666666666666666666666666",
            occurred_at=NOW,
            action=AuditAction.LOGIN,
            outcome=AuditOutcome.DENIED,
            actor_kind=AuditActorKind.USER,
            actor_user_id=USER_ID,
            reason_code=AuditReasonCode.INVALID_CREDENTIALS,
        )
    with pytest.raises(ValueError, match="must not carry a resource"):
        SecurityAuditEvent(
            event_id="aevt_77777777777777777777777777777777",
            occurred_at=NOW,
            action=AuditAction.LOGIN,
            outcome=AuditOutcome.DENIED,
            actor_kind=AuditActorKind.UNAUTHENTICATED,
            resource=build_audit_resource(AuditResourceKind.ACCOUNT, USER_ID),
            reason_code=AuditReasonCode.INVALID_CREDENTIALS,
        )
    with pytest.raises(ValueError, match="not allowed"):
        _event(actor_kind=AuditActorKind.UNAUTHENTICATED)
    with pytest.raises(ValueError, match="resource kind"):
        SecurityAuditEvent(
            event_id="aevt_88888888888888888888888888888888",
            occurred_at=NOW,
            action=AuditAction.RUN_START,
            outcome=AuditOutcome.SUCCEEDED,
            actor_kind=AuditActorKind.USER,
            actor_user_id=USER_ID,
            resource=build_audit_resource(AuditResourceKind.ACCOUNT, USER_ID),
        )


def test_successful_login_requires_the_resolved_user_actor() -> None:
    with pytest.raises(ValueError, match="authenticated user"):
        SecurityAuditEvent(
            event_id="aevt_99999999999999999999999999999999",
            occurred_at=NOW,
            action=AuditAction.LOGIN,
            outcome=AuditOutcome.SUCCEEDED,
            actor_kind=AuditActorKind.UNAUTHENTICATED,
        )


@pytest.mark.parametrize(
    "resource_id",
    (
        "/private/runtime",
        "../run-1",
        "operator@example.com",
        "run\nsecret",
        "A" * 43,
        "a" * 64,
        "ares_run_" + "a" * 64,
        "correct-horse-battery-staple",
        "a" * 256,
    ),
)
def test_audit_resource_rejects_raw_identifiers_and_secret_shaped_values(
    resource_id: str,
) -> None:
    with pytest.raises(TypeError, match="bounded factory"):
        AuditResource(AuditResourceKind.RUN, resource_id)


def test_audit_resource_factory_is_stable_domain_separated_and_non_disclosing() -> None:
    run = build_audit_resource(AuditResourceKind.RUN, "run-123")
    same_run = build_audit_resource(AuditResourceKind.RUN, "run-123")
    artifact = build_audit_resource(
        AuditResourceKind.ARTIFACT,
        "run-123",
        "artifact-1",
    )

    assert run == same_run
    assert run != build_audit_resource(AuditResourceKind.RUN, "run-124")
    assert run.resource_id.startswith("ares_run_")
    assert len(run.resource_id) == len("ares_run_") + 64
    assert "run-123" not in run.resource_id
    assert artifact.resource_id.startswith("ares_artifact_")
    assert "run-123" not in artifact.resource_id
    assert "artifact-1" not in artifact.resource_id


def test_account_audit_resource_requires_an_opaque_user_id() -> None:
    resource = build_audit_resource("account", USER_ID)
    assert resource.resource_id.startswith("ares_account_")
    assert USER_ID not in resource.resource_id
    with pytest.raises(ValueError, match="opaque user"):
        build_audit_resource("account", "member")


@pytest.mark.parametrize(
    ("kind", "identities"),
    (
        (AuditResourceKind.RUN, ("../run-1",)),
        (AuditResourceKind.RUN, ("run-1", "extra")),
        (AuditResourceKind.ARTIFACT, ("run-1",)),
        (AuditResourceKind.ARTIFACT, ("run-1", "../artifact")),
        (AuditResourceKind.REFERENCE, ("refp_ABC",)),
        (AuditResourceKind.STORAGE, ("stgp_123",)),
    ),
)
def test_audit_resource_factory_rejects_invalid_or_incomplete_coordinates(
    kind: AuditResourceKind,
    identities: tuple[str, ...],
) -> None:
    with pytest.raises(ValueError, match="stable identity|wrong stable identity"):
        build_audit_resource(kind, *identities)


@pytest.mark.parametrize(
    "reason",
    (
        "invalid_credentials",
        "INVALID CREDENTIALS",
        "ValueError:/private/path",
        "TOKEN=secret",
        "A" * 65,
    ),
)
def test_audit_reason_is_a_bounded_code_not_raw_exception_text(reason: str) -> None:
    with pytest.raises(ValueError, match="reason_code"):
        SecurityAuditEvent(
            event_id="aevt_33333333333333333333333333333333",
            occurred_at=NOW,
            action=AuditAction.LOGIN,
            outcome=AuditOutcome.FAILED,
            actor_kind=AuditActorKind.UNAUTHENTICATED,
            reason_code=reason,
        )


def test_audit_event_rejects_arbitrary_payload_fields() -> None:
    with pytest.raises(TypeError):
        SecurityAuditEvent(
            event_id="aevt_44444444444444444444444444444444",
            occurred_at=NOW,
            action=AuditAction.LOGIN,
            outcome=AuditOutcome.DENIED,
            actor_kind=AuditActorKind.UNAUTHENTICATED,
            reason_code="INVALID_CREDENTIALS",
            password="must-not-fit-the-schema",  # type: ignore[call-arg]
        )


def test_audit_actor_identity_is_explicit_and_non_spoofable() -> None:
    with pytest.raises(ValueError, match="only a user actor"):
        _event(actor_kind=AuditActorKind.USER, actor_user_id=None)
    with pytest.raises(ValueError, match="only a user actor"):
        _event(actor_kind=AuditActorKind.LOCAL_OPERATOR, actor_user_id=USER_ID)

    operator_event = _event(actor_kind=AuditActorKind.LOCAL_OPERATOR)
    assert operator_event.actor_user_id is None


def test_unsuccessful_audit_event_requires_a_stable_reason_code() -> None:
    with pytest.raises(ValueError, match="require it"):
        _event(outcome=AuditOutcome.DENIED)


def _event(**overrides: object) -> SecurityAuditEvent:
    values: dict[str, object] = {
        "event_id": "aevt_55555555555555555555555555555555",
        "occurred_at": NOW,
        "action": AuditAction.ACCOUNT_CREATE,
        "outcome": AuditOutcome.SUCCEEDED,
        "actor_kind": AuditActorKind.LOCAL_OPERATOR,
        "resource": build_audit_resource(AuditResourceKind.ACCOUNT, USER_ID),
    }
    values.update(overrides)
    return SecurityAuditEvent(**values)  # type: ignore[arg-type]


def test_new_audit_event_id_is_random_and_path_safe() -> None:
    first = new_audit_event_id()
    second = new_audit_event_id()
    assert first != second
    assert first.startswith("aevt_")
    assert len(first) == 37
