"""Build closed run-action security audit events at the service boundary."""

from __future__ import annotations

from collections.abc import Callable
from datetime import datetime

from encode_pipeline.platform.authentication import AuthenticatedPrincipal
from encode_pipeline.platform.security_audit import (
    AuditAction,
    AuditActorKind,
    AuditOutcome,
    AuditResourceKind,
    SecurityAuditEvent,
    build_audit_resource,
    new_audit_event_id,
)

_RUN_ACTIONS = frozenset(
    {AuditAction.RUN_CREATE, AuditAction.RUN_START, AuditAction.RUN_CANCEL}
)


def build_run_action_event(
    action: AuditAction,
    principal: AuthenticatedPrincipal,
    run_id: str,
    *,
    occurred_at: datetime,
    event_id_factory: Callable[[], str] = new_audit_event_id,
) -> SecurityAuditEvent:
    """Return one authenticated run-action event bound to its run target."""

    if action not in _RUN_ACTIONS:
        raise ValueError("action must be a run create, start, or cancel")
    if not isinstance(principal, AuthenticatedPrincipal):
        raise ValueError("principal must be an AuthenticatedPrincipal")
    if not callable(event_id_factory):
        raise ValueError("event_id_factory must be callable")
    return SecurityAuditEvent(
        event_id=event_id_factory(),
        occurred_at=occurred_at,
        action=action,
        outcome=AuditOutcome.SUCCEEDED,
        actor_kind=AuditActorKind.USER,
        actor_user_id=principal.user_id,
        resource=build_audit_resource(AuditResourceKind.RUN, run_id),
    )


def build_artifact_download_event(
    principal: AuthenticatedPrincipal,
    run_id: str,
    artifact_id: str,
    *,
    occurred_at: datetime,
    event_id_factory: Callable[[], str] = new_audit_event_id,
) -> SecurityAuditEvent:
    """Return one authenticated artifact-download event."""

    if not isinstance(principal, AuthenticatedPrincipal):
        raise ValueError("principal must be an AuthenticatedPrincipal")
    if not callable(event_id_factory):
        raise ValueError("event_id_factory must be callable")
    return SecurityAuditEvent(
        event_id=event_id_factory(),
        occurred_at=occurred_at,
        action=AuditAction.ARTIFACT_DOWNLOAD,
        outcome=AuditOutcome.SUCCEEDED,
        actor_kind=AuditActorKind.USER,
        actor_user_id=principal.user_id,
        resource=build_audit_resource(
            AuditResourceKind.ARTIFACT,
            run_id,
            artifact_id,
        ),
    )
