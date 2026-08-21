"""Atomic run-mutation plus security-audit persistence behavior."""

from __future__ import annotations

from datetime import datetime, timezone

import pytest

from encode_pipeline.persistence import (
    SqlAlchemyAuthenticationRepository,
    SqlAlchemyRunRepository,
    create_database_engine,
    create_session_factory,
    upgrade_database,
)
from encode_pipeline.platform.adapters import WorkflowInputs
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
from encode_pipeline.services.run_repositories import ConcurrentRunUpdateError
from encode_pipeline.services.runs import RunEventDraft

NOW = datetime(2026, 8, 19, 12, 0, tzinfo=timezone.utc)
RUN_ID = "run-audit-1"
USER_ID = "usr_" + "1" * 32


@pytest.fixture
def repositories(tmp_path):
    database_url = f"sqlite:///{tmp_path / 'platform.db'}"
    upgrade_database(database_url)
    engine = create_database_engine(database_url)
    session_factory = create_session_factory(engine)
    yield (
        SqlAlchemyRunRepository(session_factory),
        SqlAlchemyAuthenticationRepository(session_factory),
    )
    engine.dispose()


def _record(run_id: str = RUN_ID, status: RunStatus = RunStatus.CREATED) -> RunRecord:
    return RunRecord(
        run_id=run_id,
        workflow_id="encode-style-chipseq-cuttag-atac-mnase",
        inputs=WorkflowInputs(config={}, samples=None, options={}).to_dict(),
        status=status,
        created_at=NOW,
        updated_at=NOW,
        started_at=None,
        ended_at=None,
        current_stage=None,
        cancellation_reason=None,
        error=None,
        tags={},
    )


def _event(action: AuditAction = AuditAction.RUN_CREATE) -> SecurityAuditEvent:
    return SecurityAuditEvent(
        event_id=new_audit_event_id(),
        occurred_at=NOW,
        action=action,
        outcome=AuditOutcome.SUCCEEDED,
        actor_kind=AuditActorKind.USER,
        actor_user_id=USER_ID,
        resource=build_audit_resource(AuditResourceKind.RUN, RUN_ID),
    )


def _draft(status: RunStatus = RunStatus.CREATED) -> RunEventDraft:
    return RunEventDraft(
        event_type="status_changed",
        message="Run created.",
        status=status,
        context={"previous_status": None, "new_status": status.value},
    )


def test_create_run_commits_the_audit_event_in_one_transaction(repositories) -> None:
    run_repository, auth_repository = repositories
    event = _event()
    run_repository.create_run(_record(), _draft(), security_audit=event)

    assert run_repository.get_run(RUN_ID).status is RunStatus.CREATED
    (persisted,) = auth_repository.list_security_audit_events()
    assert persisted == event


def test_create_run_rolls_back_the_run_when_the_audit_insert_fails(
    repositories,
) -> None:
    run_repository, auth_repository = repositories
    duplicate = _event()
    auth_repository.record_security_audit(duplicate)

    with pytest.raises(ValueError):
        run_repository.create_run(_record(), _draft(), security_audit=duplicate)

    with pytest.raises(KeyError):
        run_repository.get_run(RUN_ID)
    (persisted,) = auth_repository.list_security_audit_events()
    assert persisted == duplicate


def test_update_run_commits_status_and_audit_atomically(repositories) -> None:
    run_repository, auth_repository = repositories
    run_repository.create_run(_record(), _draft())
    cancelled = _record(status=RunStatus.CANCELLED)
    event = _event(AuditAction.RUN_CANCEL)
    run_repository.update_run(
        cancelled,
        expected_status=RunStatus.CREATED,
        event=_draft(RunStatus.CANCELLED),
        security_audit=event,
    )
    assert run_repository.get_run(RUN_ID).status is RunStatus.CANCELLED
    (persisted,) = auth_repository.list_security_audit_events()
    assert persisted == event

    replay = _event(AuditAction.RUN_START)
    with pytest.raises(ConcurrentRunUpdateError):
        run_repository.update_run(
            _record(status=RunStatus.QUEUED),
            expected_status=RunStatus.CREATED,
            event=_draft(RunStatus.QUEUED),
            security_audit=replay,
        )
    assert auth_repository.list_security_audit_events() == (event,)


def test_security_audit_must_be_a_closed_event(repositories) -> None:
    run_repository, _auth_repository = repositories
    with pytest.raises(ValueError):
        run_repository.create_run(
            _record("run-audit-2"),
            _draft(),
            security_audit={"action": "run.create"},  # type: ignore[arg-type]
        )
