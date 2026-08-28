"""Best-effort terminal run email orchestration over durable public-safe state."""

from __future__ import annotations

from collections.abc import Callable, Mapping
from email.message import EmailMessage
from urllib.parse import quote

from encode_pipeline.platform.authentication import UserRole, UserStatus
from encode_pipeline.platform.notifications import (
    DisabledTerminalEmailSettings,
    DisabledTerminalRunNotifier,
    MAX_NOTIFICATION_RECIPIENTS,
    SmtpTerminalEmailSettings,
    TerminalRunNotifier,
    load_terminal_email_settings,
)
from encode_pipeline.platform.runs import RunQcMetric, RunRecord, RunStatus
from encode_pipeline.services.authentication_repositories import (
    AuthenticationRepository,
)
from encode_pipeline.services.run_repositories import RunEventDraft, RunRepository
from encode_pipeline.services.smtp_transport import (
    EmailTransport,
    SmtpEmailTransport,
    UNDISCLOSED_RECIPIENTS_HEADER,
)


MAX_RENDERED_QC_METRICS = 12
MAX_SAFE_TEXT_LENGTH = 160


class TerminalNotificationService:
    """Attempt one terminal email without ever changing run or result state."""

    def __init__(
        self,
        *,
        settings: SmtpTerminalEmailSettings,
        run_repository: RunRepository,
        authentication_repository: AuthenticationRepository,
        transport: EmailTransport,
    ) -> None:
        if not isinstance(settings, SmtpTerminalEmailSettings):
            raise ValueError("settings must enable terminal email")
        if run_repository is None:
            raise ValueError("run_repository is required")
        if authentication_repository is None:
            raise ValueError("authentication_repository is required")
        if transport is None or not callable(getattr(transport, "send", None)):
            raise ValueError("transport is required")
        self._settings = settings
        self._run_repository = run_repository
        self._authentication_repository = authentication_repository
        self._transport = transport

    def notify_terminal_run(
        self,
        run_id: str,
        status: RunStatus,
        *,
        include_qc: bool = False,
    ) -> None:
        """Make one isolated delivery attempt for an already-committed terminal run."""

        recipient_count = 0
        metric_count = 0
        try:
            record = self._run_repository.get_run(run_id)
            if (
                not isinstance(status, RunStatus)
                or not status.is_terminal
                or record.status is not status
            ):
                self._record_outcome(
                    run_id,
                    status if isinstance(status, RunStatus) else None,
                    event_type="terminal_email_skipped",
                    reason_code="TERMINAL_EMAIL_STATE_CHANGED",
                    recipient_count=0,
                    metric_count=0,
                )
                return
            recipients = self._recipients(run_id)
            recipient_count = len(recipients)
            if recipient_count > MAX_NOTIFICATION_RECIPIENTS:
                self._record_outcome(
                    run_id,
                    status,
                    event_type="terminal_email_skipped",
                    reason_code="TERMINAL_EMAIL_RECIPIENT_LIMIT",
                    recipient_count=recipient_count,
                    metric_count=0,
                )
                return
            if not recipients:
                self._record_outcome(
                    run_id,
                    status,
                    event_type="terminal_email_skipped",
                    reason_code="TERMINAL_EMAIL_NO_RECIPIENTS",
                    recipient_count=0,
                    metric_count=0,
                )
                return
            total_metrics, metrics = self._qc_projection(
                run_id,
                record,
                include_qc=include_qc,
            )
            metric_count = len(metrics)
            message = self._message(
                record, total_metrics=total_metrics, metrics=metrics
            )
            self._transport.send(message, recipients)
        except Exception:
            self._record_outcome(
                run_id,
                status if isinstance(status, RunStatus) else None,
                event_type="terminal_email_failed",
                reason_code="TERMINAL_EMAIL_DELIVERY_FAILED",
                recipient_count=recipient_count,
                metric_count=metric_count,
            )
            return
        self._record_outcome(
            run_id,
            status,
            event_type="terminal_email_sent",
            reason_code=(
                "TERMINAL_EMAIL_SENT"
                if status is not RunStatus.SUCCEEDED or metrics
                else "TERMINAL_EMAIL_SENT_WITHOUT_QC"
            ),
            recipient_count=recipient_count,
            metric_count=metric_count,
        )

    def _recipients(self, run_id: str) -> tuple[str, ...]:
        recipients = list(self._settings.admin_recipients)
        requester = self._run_repository.get_run_requester_user_id(run_id)
        if requester is not None:
            try:
                account = self._authentication_repository.get_account_by_id(requester)
            except Exception:
                account = None
            if (
                account is not None
                and account.role is UserRole.MEMBER
                and account.status is UserStatus.ENABLED
                and account.terminal_email_enabled
                and account.notification_email is not None
            ):
                recipients.append(account.notification_email)
        return tuple(dict.fromkeys(address.casefold() for address in recipients))

    def _qc_projection(
        self,
        run_id: str,
        record: RunRecord,
        *,
        include_qc: bool,
    ) -> tuple[int | None, tuple[RunQcMetric, ...]]:
        if record.status is not RunStatus.SUCCEEDED or not include_qc:
            return None, ()
        try:
            result_state = self._run_repository.get_result_state(run_id)
            generation = result_state.qc_generation
            if result_state.qc_outcome != "succeeded" or generation is None:
                return None, ()
            total, metrics = (
                self._run_repository.list_qc_metrics_for_terminal_notification(
                    run_id,
                    expected_generation=generation,
                    limit=MAX_RENDERED_QC_METRICS,
                )
            )
            return total, metrics
        except Exception:
            return None, ()

    def _message(
        self,
        record: RunRecord,
        *,
        total_metrics: int | None,
        metrics: tuple[RunQcMetric, ...],
    ) -> EmailMessage:
        status = record.status.value.upper()
        run_id = _safe_text(record.run_id)
        workflow_id = _safe_text(record.workflow_id)
        run_link = (
            f"{self._settings.application_base_url}/runs/"
            f"{quote(record.run_id, safe='')}"
        )
        ended_at = (
            "unknown"
            if record.ended_at is None
            else record.ended_at.isoformat(timespec="seconds")
        )
        lines = [
            f"Run status: {status}",
            f"Run ID: {run_id}",
            f"Workflow: {workflow_id}",
            f"Finished: {ended_at}",
            f"View run: {run_link}",
        ]
        if record.status is RunStatus.SUCCEEDED:
            lines.extend(["", "Persisted QC summary:"])
            if total_metrics is None or not metrics:
                lines.append("No persisted QC metrics are available.")
            else:
                lines.extend(_metric_line(metric) for metric in metrics)
                if total_metrics > len(metrics):
                    lines.append(f"… and {total_metrics - len(metrics)} more metrics")

        message = EmailMessage()
        message["From"] = self._settings.sender
        message["To"] = UNDISCLOSED_RECIPIENTS_HEADER
        message["Subject"] = f"HelixWeave run {status}: {run_id}"
        message.set_content("\n".join(lines) + "\n", subtype="plain", charset="utf-8")
        return message

    def _record_outcome(
        self,
        run_id: str,
        status: RunStatus | None,
        *,
        event_type: str,
        reason_code: str,
        recipient_count: int,
        metric_count: int,
    ) -> None:
        try:
            self._run_repository.add_event(
                run_id,
                RunEventDraft(
                    event_type=event_type,
                    message=_EVENT_MESSAGES[event_type],
                    status=status,
                    stage="notification",
                    context={
                        "reason_code": reason_code,
                        "recipient_count": recipient_count,
                        "metric_count": metric_count,
                    },
                ),
            )
        except Exception:
            return


def compose_terminal_run_notifier(
    *,
    environ: Mapping[str, str],
    run_repository: RunRepository,
    authentication_repository: AuthenticationRepository,
    transport_factory: Callable[[SmtpTerminalEmailSettings], EmailTransport] = (
        SmtpEmailTransport
    ),
) -> TerminalRunNotifier:
    """Return an explicit disabled or fully composed notifier; never optional None."""

    settings = load_terminal_email_settings(environ)
    if isinstance(settings, DisabledTerminalEmailSettings):
        return DisabledTerminalRunNotifier()
    if not isinstance(settings, SmtpTerminalEmailSettings):
        raise ValueError("terminal email settings are invalid")
    transport = transport_factory(settings)
    return TerminalNotificationService(
        settings=settings,
        run_repository=run_repository,
        authentication_repository=authentication_repository,
        transport=transport,
    )


_EVENT_MESSAGES = {
    "terminal_email_sent": "Terminal run email sent.",
    "terminal_email_failed": "Terminal run email delivery failed.",
    "terminal_email_skipped": "Terminal run email skipped.",
}


def _safe_text(value: object) -> str:
    if not isinstance(value, str):
        return "unavailable"
    normalized = " ".join(value.split())
    normalized = "".join(
        character
        for character in normalized
        if ord(character) >= 32 and ord(character) != 127
    )
    return normalized[:MAX_SAFE_TEXT_LENGTH] or "unavailable"


def _metric_line(metric: RunQcMetric) -> str:
    flag = "" if metric.qc_flag is None else f" [{metric.qc_flag.upper()}]"
    coordinate = metric.experiment_id or metric.sample_id
    scope = (
        metric.scope
        if coordinate is None
        else f"{metric.scope} {_safe_text(coordinate)}"
    )
    unit = "" if not metric.unit else f" {_safe_text(metric.unit)}"
    return (
        f"- {_safe_text(metric.display_name)} ({_safe_text(scope)}): "
        f"{metric.value}{unit}{flag}"
    )
