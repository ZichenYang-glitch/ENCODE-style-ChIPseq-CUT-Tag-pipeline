"""Build closed administrator and operator security audit events."""

from __future__ import annotations

from collections.abc import Callable
from datetime import datetime

from encode_pipeline.platform.security_audit import (
    AuditAction,
    AuditOutcome,
    AuditResource,
    AuditResourceKind,
    SecurityAuditEvent,
    build_audit_resource,
    new_audit_event_id,
)
from encode_pipeline.services.authentication_service import AuthenticationActor

_REFERENCE_ACTIONS = frozenset(
    {
        AuditAction.REFERENCE_REGISTER,
        AuditAction.REFERENCE_ENABLE,
        AuditAction.REFERENCE_DISABLE,
    }
)
_STORAGE_ACTIONS = frozenset(
    {AuditAction.STORAGE_REGISTER, AuditAction.STORAGE_ARCHIVE}
)
_RECOVERY_ACTIONS = frozenset({AuditAction.RUN_FAIL, AuditAction.RUN_REQUEUE})


def _build(
    action: AuditAction,
    allowed: frozenset[AuditAction],
    actor: AuthenticationActor,
    resource: AuditResource,
    *,
    occurred_at: datetime,
    event_id_factory: Callable[[], str],
) -> SecurityAuditEvent:
    if action not in allowed:
        raise ValueError("action is not part of this audit family")
    if not isinstance(actor, AuthenticationActor):
        raise ValueError("actor must be an AuthenticationActor")
    if not callable(event_id_factory):
        raise ValueError("event_id_factory must be callable")
    return SecurityAuditEvent(
        event_id=event_id_factory(),
        occurred_at=occurred_at,
        action=action,
        outcome=AuditOutcome.SUCCEEDED,
        actor_kind=actor.kind,
        actor_user_id=actor.user_id,
        resource=resource,
    )


def build_reference_action_event(
    action: AuditAction,
    actor: AuthenticationActor,
    profile_id: str,
    *,
    occurred_at: datetime,
    event_id_factory: Callable[[], str] = new_audit_event_id,
) -> SecurityAuditEvent:
    """Return one reference-profile administration event."""

    return _build(
        action,
        _REFERENCE_ACTIONS,
        actor,
        build_audit_resource(AuditResourceKind.REFERENCE, profile_id),
        occurred_at=occurred_at,
        event_id_factory=event_id_factory,
    )


def build_storage_action_event(
    action: AuditAction,
    actor: AuthenticationActor,
    pool_id: str,
    *,
    occurred_at: datetime,
    event_id_factory: Callable[[], str] = new_audit_event_id,
) -> SecurityAuditEvent:
    """Return one storage-pool administration event."""

    return _build(
        action,
        _STORAGE_ACTIONS,
        actor,
        build_audit_resource(AuditResourceKind.STORAGE, pool_id),
        occurred_at=occurred_at,
        event_id_factory=event_id_factory,
    )


def build_recovery_action_event(
    action: AuditAction,
    actor: AuthenticationActor,
    run_id: str,
    *,
    occurred_at: datetime,
    event_id_factory: Callable[[], str] = new_audit_event_id,
) -> SecurityAuditEvent:
    """Return one guarded run-recovery administration event."""

    return _build(
        action,
        _RECOVERY_ACTIONS,
        actor,
        build_audit_resource(AuditResourceKind.RUN, run_id),
        occurred_at=occurred_at,
        event_id_factory=event_id_factory,
    )
