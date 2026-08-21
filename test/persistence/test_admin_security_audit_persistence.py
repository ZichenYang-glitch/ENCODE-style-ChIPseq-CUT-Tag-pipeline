"""Atomic reference, storage, and recovery security-audit persistence."""

from __future__ import annotations

from datetime import datetime, timezone

import pytest

from encode_pipeline.persistence import (
    SqlAlchemyAuthenticationRepository,
    SqlAlchemyInputRegistryRepository,
    SqlAlchemyReferenceProfileRepository,
    SqlAlchemyRunRepository,
    create_database_engine,
    create_session_factory,
    upgrade_database,
)
from encode_pipeline.platform.execution import RunExecutionAssignment
from encode_pipeline.platform.input_registry import StoragePool
from encode_pipeline.platform.reference_profiles import (
    ReferenceProfile,
    ReferenceProfileRevision,
    ReferenceProfileWorkflowBinding,
    build_reference_profile_revision,
)
from encode_pipeline.platform.runs import RunRecord, RunStatus
from encode_pipeline.platform.security_audit import (
    AuditAction,
    AuditActorKind,
    AuditOutcome,
    AuditResourceKind,
    SecurityAuditEvent,
    build_audit_resource,
    new_audit_event_id,
)
from encode_pipeline.services.runs import RunEventDraft

NOW = datetime(2026, 8, 19, 12, 0, tzinfo=timezone.utc)
LATER = datetime(2026, 8, 19, 13, 0, tzinfo=timezone.utc)
PROFILE_ID = "refp_" + "1" * 32
REVISION_ID = "refpr_" + "2" * 32
REVISION_ID_2 = "refpr_" + "3" * 32
POOL_ID = "stgp_" + "4" * 32
RUN_ID = "run-recovery-audit-1"
JOB_ID = "run-execution-" + "5" * 64


@pytest.fixture
def repositories(tmp_path):
    database_url = f"sqlite:///{tmp_path / 'platform.db'}"
    upgrade_database(database_url)
    engine = create_database_engine(database_url)
    session_factory = create_session_factory(engine)
    yield {
        "reference": SqlAlchemyReferenceProfileRepository(session_factory),
        "input": SqlAlchemyInputRegistryRepository(session_factory),
        "run": SqlAlchemyRunRepository(session_factory),
        "auth": SqlAlchemyAuthenticationRepository(session_factory),
    }
    engine.dispose()


def _event(
    action: AuditAction, kind: AuditResourceKind, *ids: str
) -> SecurityAuditEvent:
    return SecurityAuditEvent(
        event_id=new_audit_event_id(),
        occurred_at=NOW,
        action=action,
        outcome=AuditOutcome.SUCCEEDED,
        actor_kind=AuditActorKind.LOCAL_OPERATOR,
        resource=build_audit_resource(kind, *ids),
    )


def _profile() -> ReferenceProfile:
    return ReferenceProfile(
        profile_id=PROFILE_ID,
        safe_key="grch38-audit",
        created_at=NOW,
    )


def _revision(revision_id: str, number: int, created_at) -> ReferenceProfileRevision:
    return build_reference_profile_revision(
        revision_id=revision_id,
        profile_id=PROFILE_ID,
        revision_number=number,
        display_name="GRCh38 audit",
        organism="Homo sapiens",
        assembly="GRCh38",
        config_key="grch38-audit-private",
        workflow_bindings=(
            ReferenceProfileWorkflowBinding(
                workflow_id="encode-style-chipseq-cuttag-atac-mnase",
                contract_version="test-contract-v1",
                identity_sha256="a" * 64,
            ),
        ),
        created_at=created_at,
    )


def _running_run_record() -> tuple[RunRecord, RunExecutionAssignment]:
    record = RunRecord(
        run_id=RUN_ID,
        workflow_id="encode-style-chipseq-cuttag-atac-mnase",
        inputs={"config": {}, "samples": None, "options": {}},
        status=RunStatus.QUEUED,
        created_at=NOW,
        updated_at=NOW,
        started_at=NOW,
        ended_at=None,
        current_stage="execution",
        cancellation_reason=None,
        error=None,
        tags={},
    )
    assignment = RunExecutionAssignment(
        run_id=RUN_ID,
        job_id=JOB_ID,
        backend="rq",
        queue_name="default",
        created_at=NOW,
        dispatched_at=NOW,
        claimed_at=None,
        cancellation_requested_at=None,
        cancellation_reason=None,
        cancellation_acknowledged_at=None,
        managed_container_scope=None,
        managed_container_endpoint_identity=None,
        requeue_requested_at=None,
        requeue_confirmed_at=None,
    )
    return record, assignment


def test_reference_create_enable_disable_write_operator_audits(repositories) -> None:
    reference = repositories["reference"]
    auth = repositories["auth"]
    created = _event(
        AuditAction.REFERENCE_REGISTER, AuditResourceKind.REFERENCE, PROFILE_ID
    )
    reference.create_profile(
        _profile(), _revision(REVISION_ID, 1, NOW), security_audit=created
    )

    enabled = _event(
        AuditAction.REFERENCE_ENABLE, AuditResourceKind.REFERENCE, PROFILE_ID
    )
    reference.set_enabled_revision(PROFILE_ID, REVISION_ID, security_audit=enabled)

    appended = _event(
        AuditAction.REFERENCE_REGISTER, AuditResourceKind.REFERENCE, PROFILE_ID
    )
    reference.append_revision(
        _revision(REVISION_ID_2, 2, LATER),
        expected_previous_revision_number=1,
        security_audit=appended,
    )

    disabled = _event(
        AuditAction.REFERENCE_DISABLE, AuditResourceKind.REFERENCE, PROFILE_ID
    )
    reference.set_enabled_revision(PROFILE_ID, None, security_audit=disabled)

    actions = sorted(event.action.value for event in auth.list_security_audit_events())
    assert actions == [
        AuditAction.REFERENCE_DISABLE.value,
        AuditAction.REFERENCE_ENABLE.value,
        AuditAction.REFERENCE_REGISTER.value,
        AuditAction.REFERENCE_REGISTER.value,
    ]


def test_reference_append_rolls_back_revision_with_audit_failure(repositories) -> None:
    reference = repositories["reference"]
    auth = repositories["auth"]
    reference.create_profile(_profile(), _revision(REVISION_ID, 1, NOW))
    duplicate = _event(
        AuditAction.REFERENCE_REGISTER, AuditResourceKind.REFERENCE, PROFILE_ID
    )
    auth.record_security_audit(duplicate)

    from encode_pipeline.services.reference_profile_repositories import (
        ReferenceProfileConflictError,
    )

    with pytest.raises(ReferenceProfileConflictError):
        reference.append_revision(
            _revision(REVISION_ID_2, 2, LATER),
            expected_previous_revision_number=1,
            security_audit=duplicate,
        )
    assert [
        revision.revision_id for revision in reference.list_revisions(PROFILE_ID)
    ] == [REVISION_ID]
    assert auth.list_security_audit_events() == (duplicate,)


def test_storage_register_writes_operator_audit(repositories) -> None:
    input_repository = repositories["input"]
    auth = repositories["auth"]
    event = _event(AuditAction.STORAGE_REGISTER, AuditResourceKind.STORAGE, POOL_ID)
    input_repository.create_storage_pool(
        StoragePool(
            storage_pool_id=POOL_ID,
            config_key="audit-ingress",
            display_name="Audit ingress",
            created_at=NOW,
        ),
        security_audit=event,
    )
    (persisted,) = auth.list_security_audit_events()
    assert persisted == event


def _stage_running_run(repositories) -> RunExecutionAssignment:
    run_repository = repositories["run"]
    record, assignment = _running_run_record()
    run_repository.create_run(
        record,
        RunEventDraft(
            event_type="status_changed",
            message="Run created.",
            status=RunStatus.CREATED,
            context={"previous_status": None, "new_status": RunStatus.CREATED.value},
        ),
    )
    run_repository.ensure_execution_assignment(
        assignment,
        expected_status=RunStatus.QUEUED,
    )
    return assignment


def test_recovery_requeue_and_fail_write_operator_audits(repositories) -> None:
    auth = repositories["auth"]
    run_repository = repositories["run"]
    assignment = _stage_running_run(repositories)

    requeue = _event(AuditAction.RUN_REQUEUE, AuditResourceKind.RUN, RUN_ID)
    preparation = run_repository.prepare_execution_requeue(
        RUN_ID,
        expected_status=RunStatus.QUEUED,
        expected_assignment=assignment,
        requested_at=LATER,
        event=RunEventDraft(
            event_type="requeue_requested",
            message="Requeue requested.",
            status=RunStatus.RUNNING,
        ),
        security_audit=requeue,
    )
    assert preparation.created is True

    fail = _event(AuditAction.RUN_FAIL, AuditResourceKind.RUN, RUN_ID)
    failed_record = _running_run_record()[0]
    from dataclasses import replace

    from encode_pipeline.platform.results import Issue

    failed_record = replace(
        failed_record,
        status=RunStatus.FAILED,
        updated_at=LATER,
        ended_at=LATER,
        error=Issue(
            code="RUN_FAILED_BY_ADMIN_RECOVERY",
            message="Run was failed by an administrator.",
            severity="error",
        ),
    )
    changed = run_repository.fail_run_by_recovery(
        failed_record,
        expected_status=RunStatus.QUEUED,
        expected_assignment=preparation.assignment,
        event=RunEventDraft(
            event_type="status_changed",
            message="Run failed by recovery.",
            status=RunStatus.FAILED,
        ),
        security_audit=fail,
    )
    assert changed is True

    actions = sorted(event.action.value for event in auth.list_security_audit_events())
    assert actions == [AuditAction.RUN_FAIL.value, AuditAction.RUN_REQUEUE.value]
