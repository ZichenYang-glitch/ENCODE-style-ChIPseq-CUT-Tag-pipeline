"""Tests for workflow-neutral terminal email contracts and settings."""

from __future__ import annotations

import pytest

from encode_pipeline.platform.notifications import (
    DisabledTerminalEmailSettings,
    DisabledTerminalRunNotifier,
    MAX_NOTIFICATION_RECIPIENTS,
    SMTP_HOST_ENV,
    SMTP_PASSWORD_ENV,
    SMTP_TLS_MODE_ENV,
    SMTP_USERNAME_ENV,
    TERMINAL_EMAIL_ADMIN_RECIPIENTS_ENV,
    TERMINAL_EMAIL_APPLICATION_BASE_URL_ENV,
    TERMINAL_EMAIL_ENABLED_ENV,
    TERMINAL_EMAIL_FROM_ENV,
    SmtpTerminalEmailSettings,
    SmtpTlsMode,
    load_terminal_email_settings,
    normalize_application_base_url,
    normalize_notification_email,
)
from encode_pipeline.platform.runs import RunStatus


def _enabled_environment(**overrides: str) -> dict[str, str]:
    environment = {
        TERMINAL_EMAIL_ENABLED_ENV: "true",
        TERMINAL_EMAIL_ADMIN_RECIPIENTS_ENV: "admin@example.test",
        TERMINAL_EMAIL_FROM_ENV: "helixweave@example.test",
        TERMINAL_EMAIL_APPLICATION_BASE_URL_ENV: "https://lab.example.test:8443",
        SMTP_HOST_ENV: "smtp.example.test",
        SMTP_TLS_MODE_ENV: "starttls",
        SMTP_USERNAME_ENV: "smtp-user",
        SMTP_PASSWORD_ENV: "smtp-password",
    }
    environment.update(overrides)
    return environment


def test_disabled_notifier_is_an_explicit_noop() -> None:
    notifier = DisabledTerminalRunNotifier()

    assert (
        notifier.notify_terminal_run(
            "run_ignored",
            RunStatus.SUCCEEDED,
            include_qc=True,
        )
        is None
    )


def test_settings_default_to_explicit_disabled_state() -> None:
    assert load_terminal_email_settings({}) == DisabledTerminalEmailSettings()


def test_enabled_settings_are_complete_normalized_and_redacted() -> None:
    environment = _enabled_environment(
        **{
            TERMINAL_EMAIL_ADMIN_RECIPIENTS_ENV: (
                "ADMIN@example.test,member@example.test,admin@example.test"
            )
        }
    )

    settings = load_terminal_email_settings(environment)

    assert isinstance(settings, SmtpTerminalEmailSettings)
    assert settings.admin_recipients == (
        "admin@example.test",
        "member@example.test",
    )
    assert settings.sender == "helixweave@example.test"
    assert settings.application_base_url == "https://lab.example.test:8443"
    assert settings.smtp_port == 587
    assert settings.tls_mode is SmtpTlsMode.STARTTLS
    rendered = repr(settings)
    for private_value in (
        "admin@example.test",
        "helixweave@example.test",
        "lab.example.test",
        "smtp.example.test",
        "smtp-user",
        "smtp-password",
    ):
        assert private_value not in rendered


@pytest.mark.parametrize(
    "environment",
    (
        {TERMINAL_EMAIL_ENABLED_ENV: "maybe"},
        {TERMINAL_EMAIL_ENABLED_ENV: "true"},
        _enabled_environment(**{SMTP_USERNAME_ENV: ""}),
        _enabled_environment(**{SMTP_PASSWORD_ENV: ""}),
        _enabled_environment(**{SMTP_TLS_MODE_ENV: "plaintext"}),
        _enabled_environment(
            **{
                SMTP_TLS_MODE_ENV: "local_plaintext",
                SMTP_HOST_ENV: "smtp.example.test",
            }
        ),
    ),
)
def test_enabled_settings_fail_closed_for_incomplete_or_invalid_config(
    environment: dict[str, str],
) -> None:
    with pytest.raises(ValueError):
        load_terminal_email_settings(environment)


def test_local_plaintext_smtp_is_limited_to_loopback() -> None:
    settings = load_terminal_email_settings(
        _enabled_environment(
            **{
                SMTP_TLS_MODE_ENV: "local_plaintext",
                SMTP_HOST_ENV: "127.0.0.1",
                SMTP_USERNAME_ENV: "",
                SMTP_PASSWORD_ENV: "",
            }
        )
    )

    assert isinstance(settings, SmtpTerminalEmailSettings)
    assert settings.smtp_port == 25
    assert settings.tls_mode is SmtpTlsMode.LOCAL_PLAINTEXT


def test_local_plaintext_smtp_rejects_credentials() -> None:
    with pytest.raises(ValueError, match="does not accept credentials"):
        load_terminal_email_settings(
            _enabled_environment(
                **{
                    SMTP_TLS_MODE_ENV: "local_plaintext",
                    SMTP_HOST_ENV: "localhost",
                }
            )
        )


@pytest.mark.parametrize(
    "value",
    (
        "Display Name <admin@example.test>",
        "admin@example.test\r\nBcc: stolen@example.test",
        "admin@@example.test",
        ".admin@example.test",
        "admin..user@example.test",
        "admin@example..test",
        "ad min@example.test",
        "δοκιμή@example.test",
    ),
)
def test_notification_email_rejects_identity_and_header_ambiguity(value: str) -> None:
    with pytest.raises(ValueError):
        normalize_notification_email(value)


def test_admin_recipient_limit_is_applied_after_deduplication() -> None:
    recipients = tuple(f"admin-{index}@example.test" for index in range(33))

    with pytest.raises(ValueError, match=str(MAX_NOTIFICATION_RECIPIENTS)):
        SmtpTerminalEmailSettings(
            admin_recipients=recipients,
            sender="sender@example.test",
            application_base_url="https://lab.example.test",
            smtp_host="smtp.example.test",
            smtp_port=587,
            tls_mode="starttls",
        )


@pytest.mark.parametrize(
    "value",
    (
        "http://lab.example.test",
        "https://user:password@lab.example.test",
        "https://lab.example.test/runs",
        "https://lab.example.test?private=value",
        "javascript:alert(1)",
    ),
)
def test_application_base_url_rejects_non_origin_or_insecure_remote_values(
    value: str,
) -> None:
    with pytest.raises(ValueError):
        normalize_application_base_url(value)


def test_application_base_url_allows_explicit_loopback_http() -> None:
    assert normalize_application_base_url("http://localhost:5173/") == (
        "http://localhost:5173"
    )
