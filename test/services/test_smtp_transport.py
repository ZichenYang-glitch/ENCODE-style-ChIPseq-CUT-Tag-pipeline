"""Tests for the private-envelope standard-library SMTP transport."""

from __future__ import annotations

import smtplib
from email.message import EmailMessage
from unittest.mock import ANY

import pytest

import encode_pipeline.services.smtp_transport as smtp_transport_module
from encode_pipeline.platform.notifications import SmtpTerminalEmailSettings
from encode_pipeline.services.smtp_transport import (
    UNDISCLOSED_RECIPIENTS_HEADER,
    SmtpDeliveryError,
    SmtpEmailTransport,
)


def _settings(**overrides: object) -> SmtpTerminalEmailSettings:
    values = {
        "admin_recipients": ("admin@example.test",),
        "sender": "helixweave@example.test",
        "application_base_url": "https://lab.example.test",
        "smtp_host": "smtp.example.test",
        "smtp_port": 587,
        "tls_mode": "starttls",
        "smtp_username": "smtp-user",
        "smtp_password": "smtp-password",
        "timeout_seconds": 4.0,
    }
    values.update(overrides)
    return SmtpTerminalEmailSettings(**values)


def _message() -> EmailMessage:
    message = EmailMessage()
    message["From"] = "helixweave@example.test"
    message["To"] = UNDISCLOSED_RECIPIENTS_HEADER
    message["Subject"] = "HelixWeave run succeeded"
    message.set_content("Run run_123 succeeded.\n")
    return message


class _FakeSmtp:
    def __init__(self, *arguments: object, **keywords: object) -> None:
        self.arguments = arguments
        self.keywords = keywords
        self.calls: list[tuple[object, ...]] = []
        self.refused: dict[str, tuple[int, bytes]] = {}

    def __enter__(self):
        self.calls.append(("enter",))
        return self

    def __exit__(self, *_arguments: object) -> None:
        self.calls.append(("exit",))

    def ehlo(self) -> None:
        self.calls.append(("ehlo",))

    def starttls(self, *, context: object) -> None:
        self.calls.append(("starttls", context))

    def login(self, username: str, password: str) -> None:
        self.calls.append(("login", username, password))

    def send_message(
        self,
        message: EmailMessage,
        *,
        from_addr: str,
        to_addrs: list[str],
    ) -> dict[str, tuple[int, bytes]]:
        self.calls.append(("send_message", message, from_addr, to_addrs))
        return self.refused


def test_starttls_transport_uses_private_deduplicated_envelope(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    clients: list[_FakeSmtp] = []

    def create_client(*args: object, **kwargs: object) -> _FakeSmtp:
        client = _FakeSmtp(*args, **kwargs)
        clients.append(client)
        return client

    tls_context = object()
    monkeypatch.setattr(smtp_transport_module.smtplib, "SMTP", create_client)
    monkeypatch.setattr(
        smtp_transport_module.ssl,
        "create_default_context",
        lambda: tls_context,
    )

    SmtpEmailTransport(_settings()).send(
        _message(),
        ("ADMIN@example.test", "member@example.test", "admin@example.test"),
    )

    assert len(clients) == 1
    client = clients[0]
    assert client.arguments == ("smtp.example.test", 587)
    assert client.keywords == {"timeout": 4.0}
    assert client.calls == [
        ("enter",),
        ("ehlo",),
        ("starttls", tls_context),
        ("ehlo",),
        ("login", "smtp-user", "smtp-password"),
        (
            "send_message",
            ANY,
            "helixweave@example.test",
            ["admin@example.test", "member@example.test"],
        ),
        ("exit",),
    ]


def test_implicit_tls_uses_verified_ssl_constructor(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    clients: list[_FakeSmtp] = []

    def create_client(*args: object, **kwargs: object) -> _FakeSmtp:
        client = _FakeSmtp(*args, **kwargs)
        clients.append(client)
        return client

    tls_context = object()
    monkeypatch.setattr(smtp_transport_module.smtplib, "SMTP_SSL", create_client)
    monkeypatch.setattr(
        smtp_transport_module.ssl,
        "create_default_context",
        lambda: tls_context,
    )

    SmtpEmailTransport(_settings(smtp_port=465, tls_mode="implicit_tls")).send(
        _message(), ("admin@example.test",)
    )

    assert clients[0].arguments == ("smtp.example.test", 465)
    assert clients[0].keywords == {"timeout": 4.0, "context": tls_context}
    assert not any(call[0] == "starttls" for call in clients[0].calls)


def test_local_plaintext_transport_uses_no_tls_or_auth(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    clients: list[_FakeSmtp] = []

    def create_client(*args: object, **kwargs: object) -> _FakeSmtp:
        client = _FakeSmtp(*args, **kwargs)
        clients.append(client)
        return client

    monkeypatch.setattr(smtp_transport_module.smtplib, "SMTP", create_client)

    SmtpEmailTransport(
        _settings(
            smtp_host="localhost",
            smtp_port=25,
            tls_mode="local_plaintext",
            smtp_username=None,
            smtp_password=None,
        )
    ).send(_message(), ("admin@example.test",))

    assert not any(call[0] in {"starttls", "login"} for call in clients[0].calls)


@pytest.mark.parametrize(
    "unsafe_kind",
    (
        "visible_to",
        "visible_cc",
        "attachment",
    ),
)
def test_transport_rejects_recipient_headers_and_attachments(
    unsafe_kind: str,
) -> None:
    message = _message()
    if unsafe_kind == "visible_to":
        del message["To"]
        message["To"] = "admin@example.test"
    elif unsafe_kind == "visible_cc":
        message["Cc"] = "member@example.test"
    else:
        message.add_attachment(
            b"private",
            maintype="application",
            subtype="octet-stream",
            filename="result.txt",
        )

    with pytest.raises(ValueError):
        SmtpEmailTransport(_settings()).send(message, ("admin@example.test",))


def test_transport_propagates_smtp_failure_without_output(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    def fail_connect(*_args: object, **_kwargs: object) -> None:
        raise smtplib.SMTPConnectError(421, b"private server response")

    monkeypatch.setattr(smtp_transport_module.smtplib, "SMTP", fail_connect)

    with pytest.raises(SmtpDeliveryError, match="SMTP delivery failed") as caught:
        SmtpEmailTransport(_settings()).send(
            _message(),
            ("admin@example.test",),
        )

    assert "private server response" not in str(caught.value)
    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err == ""


def test_transport_treats_partial_recipient_refusal_as_failure(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    client = _FakeSmtp()
    client.refused = {"member@example.test": (550, b"private reply")}
    monkeypatch.setattr(
        smtp_transport_module.smtplib,
        "SMTP",
        lambda *_args, **_kwargs: client,
    )
    monkeypatch.setattr(
        smtp_transport_module.ssl,
        "create_default_context",
        lambda: object(),
    )

    with pytest.raises(SmtpDeliveryError, match="SMTP delivery failed") as caught:
        SmtpEmailTransport(_settings()).send(
            _message(),
            ("admin@example.test", "member@example.test"),
        )

    assert "private reply" not in str(caught.value)
