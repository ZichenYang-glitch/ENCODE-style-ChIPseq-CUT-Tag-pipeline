"""Focused orchestration tests for best-effort terminal run email."""

from __future__ import annotations

from datetime import datetime, timezone
from decimal import Decimal
from types import SimpleNamespace

import pytest

from encode_pipeline.platform.authentication import UserRole, UserStatus
from encode_pipeline.platform.notifications import SmtpTerminalEmailSettings
from encode_pipeline.platform.runs import RunQcMetric, RunRecord, RunStatus
from encode_pipeline.services.terminal_notifications import TerminalNotificationService


NOW = datetime(2026, 8, 27, 12, 0, tzinfo=timezone.utc)


def _settings() -> SmtpTerminalEmailSettings:
    return SmtpTerminalEmailSettings(
        admin_recipients=("admin@example.test",),
        sender="helixweave@example.test",
        application_base_url="https://lab.example.test",
        smtp_host="smtp.example.test",
        smtp_port=587,
        tls_mode="starttls",
        smtp_username="smtp-user",
        smtp_password="smtp-password",
    )


def _record(status: RunStatus = RunStatus.SUCCEEDED) -> RunRecord:
    return RunRecord(
        run_id="run-1",
        workflow_id="workflow-neutral-adapter",
        inputs={"private": "/workspace/private.fastq.gz"},
        status=status,
        created_at=NOW,
        updated_at=NOW,
        started_at=NOW,
        ended_at=NOW,
        current_stage="execution",
        cancellation_reason="private cancellation text",
        error=None,
        tags={"private": "payload"},
    )


def _metric() -> RunQcMetric:
    return RunQcMetric(
        metric_id="qcmetric-" + "1" * 64,
        run_id="run-1",
        metric_key="alignment.mapping_rate",
        display_name="Mapping rate",
        value=Decimal("0.987"),
        unit="fraction",
        scope="sample",
        sample_id="S1",
        experiment_id=None,
        assay=None,
        qc_flag="pass",
        source_artifact_id="artifact-1",
        produced_at=NOW,
    )


class _RunRepository:
    def __init__(self, record: RunRecord) -> None:
        self.record = record
        self.requester_user_id = "usr_11111111111111111111111111111111"
        self.events = []
        self.metrics = (_metric(),)
        self.total_metrics = 14
        self.qc_read_raises = False

    def get_run(self, run_id: str) -> RunRecord:
        assert run_id == "run-1"
        return self.record

    def get_run_requester_user_id(self, run_id: str) -> str | None:
        assert run_id == "run-1"
        return self.requester_user_id

    def get_result_state(self, run_id: str):
        assert run_id == "run-1"
        if self.qc_read_raises:
            raise RuntimeError("/workspace/secret indexing failure")
        return SimpleNamespace(
            qc_outcome="succeeded",
            qc_generation="qcgen-" + "2" * 64,
        )

    def list_qc_metrics_for_terminal_notification(
        self,
        run_id: str,
        *,
        expected_generation: str,
        limit: int,
    ):
        assert run_id == "run-1"
        assert expected_generation == "qcgen-" + "2" * 64
        assert limit == 12
        return self.total_metrics, self.metrics

    def add_event(self, run_id: str, event):
        assert run_id == "run-1"
        self.events.append(event)


class _AuthenticationRepository:
    def __init__(self) -> None:
        self.account = SimpleNamespace(
            role=UserRole.MEMBER,
            status=UserStatus.ENABLED,
            terminal_email_enabled=True,
            notification_email="member@example.test",
        )

    def get_account_by_id(self, user_id: str):
        assert user_id == "usr_11111111111111111111111111111111"
        return self.account


class _Transport:
    def __init__(self, *, raises: bool = False) -> None:
        self.raises = raises
        self.deliveries = []

    def send(self, message, envelope_recipients) -> None:
        if self.raises:
            raise RuntimeError("smtp reply contains secret credentials")
        self.deliveries.append((message, envelope_recipients))


def _service(repository, authentication, transport) -> TerminalNotificationService:
    return TerminalNotificationService(
        settings=_settings(),
        run_repository=repository,
        authentication_repository=authentication,
        transport=transport,
    )


@pytest.mark.parametrize("delivery_fails", [False, True])
@pytest.mark.parametrize("recording_fails", [False, True])
def test_outcome_diagnostic_preserves_delivery_and_recording_semantics(
    delivery_fails, recording_fails, tmp_path, monkeypatch, assert_failure_diagnostics
):
    private = str(tmp_path / "PRIVATE_EVENT_FAILURE")
    monkeypatch.setenv("HELIX_TEST_PRIVATE", "PRIVATE_ENV_VALUE")
    repository = _RunRepository(_record(RunStatus.FAILED))
    original_record = repository.record
    authentication = _AuthenticationRepository()
    authentication.account.notification_email = "PRIVATE_RECIPIENT@example.test"
    calls = []

    class Transport:
        def send(self, message, envelope_recipients):
            calls.append(("send", message.get_content(), envelope_recipients))
            if delivery_fails:
                raise RuntimeError(private + " PRIVATE_EXCEPTION_TEXT")

    add_event = repository.add_event

    def record(run_id, event):
        calls.append(("record", event.event_type))
        if recording_fails:
            raise RuntimeError(private + " PRIVATE_EXCEPTION_TEXT")
        return add_event(run_id, event)

    monkeypatch.setattr(repository, "add_event", record)
    assert (
        _service(repository, authentication, Transport()).notify_terminal_run(
            "run-1", RunStatus.FAILED
        )
        is None
    )
    outcome = "terminal_email_failed" if delivery_fails else "terminal_email_sent"
    assert [call[0] for call in calls] == ["send", "record"]
    assert calls[1] == ("record", outcome)
    assert repository.record == original_record
    assert [event.event_type for event in repository.events] == (
        [] if recording_fails else [outcome]
    )
    assert_failure_diagnostics(
        [
            (
                "terminal_notifications",
                "record_outcome",
                "TERMINAL_EMAIL_EVENT_RECORD_FAILED",
            )
        ]
        if recording_fails
        else [],
        private=(
            private,
            "PRIVATE_EXCEPTION_TEXT",
            "PRIVATE_ENV_VALUE",
            "PRIVATE_RECIPIENT",
            calls[0][1],
            "/workspace/private.fastq.gz",
            "smtp-password",
        ),
    )


def test_success_email_uses_private_envelope_and_current_bounded_qc_projection() -> (
    None
):
    repository = _RunRepository(_record())
    transport = _Transport()

    _service(repository, _AuthenticationRepository(), transport).notify_terminal_run(
        "run-1",
        RunStatus.SUCCEEDED,
        include_qc=True,
    )

    assert len(transport.deliveries) == 1
    message, recipients = transport.deliveries[0]
    assert recipients == ("admin@example.test", "member@example.test")
    assert message["To"] == "undisclosed-recipients:;"
    body = message.get_content()
    assert "Mapping rate (sample S1): 0.987 fraction [PASS]" in body
    assert "… and 13 more metrics" in body
    assert "https://lab.example.test/runs/run-1" in body
    for private_text in (
        "/workspace/private.fastq.gz",
        "private cancellation text",
        "payload",
        "smtp-password",
    ):
        assert private_text not in body
    (event,) = repository.events
    assert event.event_type == "terminal_email_sent"
    assert event.context == {
        "reason_code": "TERMINAL_EMAIL_SENT",
        "recipient_count": 2,
        "metric_count": 1,
    }


def test_missing_or_failed_qc_still_sends_success_without_private_failure_text() -> (
    None
):
    repository = _RunRepository(_record())
    repository.qc_read_raises = True
    transport = _Transport()

    _service(repository, _AuthenticationRepository(), transport).notify_terminal_run(
        "run-1",
        RunStatus.SUCCEEDED,
        include_qc=True,
    )

    message, _recipients = transport.deliveries[0]
    assert "No persisted QC metrics are available." in message.get_content()
    assert "/workspace/secret" not in message.get_content()
    assert repository.events[0].context == {
        "reason_code": "TERMINAL_EMAIL_SENT_WITHOUT_QC",
        "recipient_count": 2,
        "metric_count": 0,
    }


def test_opted_out_member_is_excluded_and_delivery_failure_is_redacted() -> None:
    repository = _RunRepository(_record(RunStatus.FAILED))
    authentication = _AuthenticationRepository()
    authentication.account.terminal_email_enabled = False
    transport = _Transport(raises=True)

    _service(repository, authentication, transport).notify_terminal_run(
        "run-1",
        RunStatus.FAILED,
    )

    (event,) = repository.events
    assert event.event_type == "terminal_email_failed"
    assert event.context == {
        "reason_code": "TERMINAL_EMAIL_DELIVERY_FAILED",
        "recipient_count": 1,
        "metric_count": 0,
    }
    assert "secret" not in repr(event)
    assert "example.test" not in repr(event)


def test_stale_terminal_callback_is_skipped_without_delivery() -> None:
    repository = _RunRepository(_record(RunStatus.CANCELLED))
    transport = _Transport()

    _service(repository, _AuthenticationRepository(), transport).notify_terminal_run(
        "run-1",
        RunStatus.FAILED,
    )

    assert transport.deliveries == []
    assert repository.events[0].context == {
        "reason_code": "TERMINAL_EMAIL_STATE_CHANGED",
        "recipient_count": 0,
        "metric_count": 0,
    }


def test_member_cannot_push_the_private_envelope_over_the_total_limit() -> None:
    repository = _RunRepository(_record(RunStatus.FAILED))
    transport = _Transport()
    settings = SmtpTerminalEmailSettings(
        admin_recipients=tuple(f"admin-{index}@example.test" for index in range(32)),
        sender="helixweave@example.test",
        application_base_url="https://lab.example.test",
        smtp_host="smtp.example.test",
        smtp_port=587,
        tls_mode="starttls",
        smtp_username="smtp-user",
        smtp_password="smtp-password",
    )
    service = TerminalNotificationService(
        settings=settings,
        run_repository=repository,
        authentication_repository=_AuthenticationRepository(),
        transport=transport,
    )

    service.notify_terminal_run("run-1", RunStatus.FAILED)

    assert transport.deliveries == []
    assert repository.events[0].context == {
        "reason_code": "TERMINAL_EMAIL_RECIPIENT_LIMIT",
        "recipient_count": 33,
        "metric_count": 0,
    }
