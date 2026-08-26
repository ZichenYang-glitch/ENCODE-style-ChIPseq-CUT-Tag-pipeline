"""Small standard-library SMTP transport for terminal run notifications."""

from __future__ import annotations

import smtplib
import ssl
from email.message import EmailMessage
from typing import Protocol

from encode_pipeline.platform.notifications import (
    MAX_NOTIFICATION_RECIPIENTS,
    SmtpTerminalEmailSettings,
    SmtpTlsMode,
    normalize_notification_email,
)


UNDISCLOSED_RECIPIENTS_HEADER = "undisclosed-recipients:;"
MAX_ENVELOPE_RECIPIENTS = MAX_NOTIFICATION_RECIPIENTS


class EmailTransport(Protocol):
    """Deliver an already-rendered message to private envelope recipients."""

    def send(
        self,
        message: EmailMessage,
        envelope_recipients: tuple[str, ...],
    ) -> None: ...


class SmtpDeliveryError(RuntimeError):
    """A stable public-safe transport failure without SMTP response details."""


class SmtpEmailTransport:
    """Deliver plain-text messages through one validated SMTP endpoint."""

    def __init__(self, settings: SmtpTerminalEmailSettings) -> None:
        if not isinstance(settings, SmtpTerminalEmailSettings):
            raise ValueError("settings must enable SMTP terminal email")
        self._settings = settings

    def send(
        self,
        message: EmailMessage,
        envelope_recipients: tuple[str, ...],
    ) -> None:
        recipients = _private_envelope_recipients(envelope_recipients)
        _validate_message(message, sender=self._settings.sender)
        try:
            if self._settings.tls_mode is SmtpTlsMode.IMPLICIT_TLS:
                context = ssl.create_default_context()
                client = smtplib.SMTP_SSL(
                    self._settings.smtp_host,
                    self._settings.smtp_port,
                    timeout=self._settings.timeout_seconds,
                    context=context,
                )
            else:
                client = smtplib.SMTP(
                    self._settings.smtp_host,
                    self._settings.smtp_port,
                    timeout=self._settings.timeout_seconds,
                )

            with client as smtp:
                if self._settings.tls_mode is SmtpTlsMode.STARTTLS:
                    context = ssl.create_default_context()
                    smtp.ehlo()
                    smtp.starttls(context=context)
                    smtp.ehlo()
                if self._settings.smtp_username is not None:
                    smtp.login(
                        self._settings.smtp_username,
                        self._settings.smtp_password,
                    )
                refused = smtp.send_message(
                    message,
                    from_addr=self._settings.sender,
                    to_addrs=list(recipients),
                )
                if refused:
                    raise SmtpDeliveryError("SMTP delivery failed")
        except SmtpDeliveryError:
            raise
        except (OSError, smtplib.SMTPException):
            raise SmtpDeliveryError("SMTP delivery failed") from None


def _private_envelope_recipients(value: object) -> tuple[str, ...]:
    if not isinstance(value, tuple) or not value:
        raise ValueError("envelope_recipients must be a non-empty tuple")
    recipients: list[str] = []
    seen: set[str] = set()
    for candidate in value:
        recipient = normalize_notification_email(candidate)
        if recipient not in seen:
            recipients.append(recipient)
            seen.add(recipient)
    if len(recipients) > MAX_ENVELOPE_RECIPIENTS:
        raise ValueError(
            f"envelope_recipients must contain at most {MAX_ENVELOPE_RECIPIENTS}"
        )
    return tuple(recipients)


def _validate_message(message: object, *, sender: str) -> None:
    if not isinstance(message, EmailMessage):
        raise ValueError("message must be an EmailMessage")
    if message.get("From") != sender:
        raise ValueError("message From must match the configured sender")
    if message.get("To") != UNDISCLOSED_RECIPIENTS_HEADER:
        raise ValueError("message must not expose envelope recipients")
    if message.get("Cc") is not None or message.get("Bcc") is not None:
        raise ValueError("message must not expose envelope recipients")
    if message.defects:
        raise ValueError("message contains invalid headers")
    if message.is_multipart() or message.get_content_type() != "text/plain":
        raise ValueError("message must be plain text without attachments")
