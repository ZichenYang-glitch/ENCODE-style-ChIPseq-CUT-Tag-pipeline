"""Run-action security audit construction and service forwarding."""

from __future__ import annotations

from datetime import datetime, timezone

import pytest

from encode_pipeline.platform.authentication import (
    AuthenticatedPrincipal,
    UserRole,
)
from encode_pipeline.platform.adapters import WorkflowInputs
from encode_pipeline.platform.runs import RunStatus
from encode_pipeline.platform.security_audit import (
    AuditAction,
    AuditActorKind,
    AuditOutcome,
    AuditResourceKind,
)
from encode_pipeline.services.defaults import create_default_workflow_registry
from encode_pipeline.services.run_repositories import InMemoryRunRepository
from encode_pipeline.services.run_security_audit import (
    build_artifact_download_event,
    build_run_action_event,
)
from encode_pipeline.services.runs import RunService

NOW = datetime(2026, 8, 19, 12, 0, tzinfo=timezone.utc)
PRINCIPAL = AuthenticatedPrincipal(
    user_id="usr_" + "1" * 32,
    username="root-admin",
    role=UserRole.ADMINISTRATOR,
)
WORKFLOW_ID = "encode-style-chipseq-cuttag-atac-mnase"


@pytest.fixture
def service() -> tuple[RunService, InMemoryRunRepository]:
    repository = InMemoryRunRepository()
    return (
        RunService(
            registry=create_default_workflow_registry(),
            id_factory=lambda: "run-audit-1",
            repository=repository,
        ),
        repository,
    )


def test_build_run_action_event_binds_actor_and_opaque_run_target() -> None:
    event = build_run_action_event(
        AuditAction.RUN_START,
        PRINCIPAL,
        "run-1",
        occurred_at=NOW,
    )
    assert event.action is AuditAction.RUN_START
    assert event.outcome is AuditOutcome.SUCCEEDED
    assert event.actor_kind is AuditActorKind.USER
    assert event.actor_user_id == PRINCIPAL.user_id
    assert event.resource is not None
    assert event.resource.kind is AuditResourceKind.RUN
    assert event.resource.resource_id.startswith("ares_run_")
    assert "run-1" not in event.resource.resource_id

    with pytest.raises(ValueError):
        build_run_action_event(
            AuditAction.LOGIN,
            PRINCIPAL,
            "run-1",
            occurred_at=NOW,
        )


def test_build_artifact_download_event_uses_the_two_part_target() -> None:
    event = build_artifact_download_event(
        PRINCIPAL,
        "run-1",
        "artifact-1",
        occurred_at=NOW,
    )
    assert event.action is AuditAction.ARTIFACT_DOWNLOAD
    assert event.resource is not None
    assert event.resource.kind is AuditResourceKind.ARTIFACT
    assert event.resource.resource_id.startswith("ares_artifact_")


def test_create_run_records_the_audit_in_memory(service) -> None:
    run_service, repository = service
    record = run_service.create_run(
        WORKFLOW_ID,
        WorkflowInputs(config={}, samples=None, options={}),
        security_audit_actor=PRINCIPAL,
    )
    (event,) = repository._security_audits
    assert event.action is AuditAction.RUN_CREATE
    assert event.actor_user_id == PRINCIPAL.user_id

    unauthenticated = RunService(
        registry=create_default_workflow_registry(),
        id_factory=lambda: "run-audit-2",
        repository=repository,
    )
    unauthenticated.create_run(
        WORKFLOW_ID,
        WorkflowInputs(config={"key": "other"}, samples=None, options={}),
    )
    assert repository._security_audits == [event]


def test_transition_and_cancel_forward_the_closed_event(service) -> None:
    run_service, repository = service
    record = run_service.create_run(
        WORKFLOW_ID,
        WorkflowInputs(config={}, samples=None, options={}),
    )
    started = build_run_action_event(
        AuditAction.RUN_START,
        PRINCIPAL,
        record.run_id,
        occurred_at=NOW,
    )
    run_service.transition_run(
        record.run_id,
        RunStatus.VALIDATING,
        security_audit=started,
    )
    cancelled = build_run_action_event(
        AuditAction.RUN_CANCEL,
        PRINCIPAL,
        record.run_id,
        occurred_at=NOW,
    )
    run_service.cancel_run(record.run_id, reason="test", security_audit=cancelled)

    assert repository._security_audits == [started, cancelled]
    assert repository.get_run(record.run_id).status is RunStatus.CANCELLED
