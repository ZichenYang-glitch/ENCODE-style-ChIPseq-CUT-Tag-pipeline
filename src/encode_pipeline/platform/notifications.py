"""Workflow-neutral terminal email notification contracts and settings."""

from __future__ import annotations

import ipaddress
import math
import os
import re
from collections.abc import Mapping
from dataclasses import dataclass, field
from enum import Enum
from typing import Protocol
from urllib.parse import urlsplit, urlunsplit

from encode_pipeline.platform.runs import RunStatus


TERMINAL_EMAIL_ENABLED_ENV = "HELIXWEAVE_TERMINAL_EMAIL_ENABLED"
TERMINAL_EMAIL_ADMIN_RECIPIENTS_ENV = "HELIXWEAVE_TERMINAL_EMAIL_ADMIN_RECIPIENTS"
TERMINAL_EMAIL_FROM_ENV = "HELIXWEAVE_TERMINAL_EMAIL_FROM"
TERMINAL_EMAIL_APPLICATION_BASE_URL_ENV = (
    "HELIXWEAVE_TERMINAL_EMAIL_APPLICATION_BASE_URL"
)
SMTP_HOST_ENV = "HELIXWEAVE_SMTP_HOST"
SMTP_PORT_ENV = "HELIXWEAVE_SMTP_PORT"
SMTP_TLS_MODE_ENV = "HELIXWEAVE_SMTP_TLS_MODE"
SMTP_USERNAME_ENV = "HELIXWEAVE_SMTP_USERNAME"
SMTP_PASSWORD_ENV = "HELIXWEAVE_SMTP_PASSWORD"
SMTP_TIMEOUT_SECONDS_ENV = "HELIXWEAVE_SMTP_TIMEOUT_SECONDS"
TERMINAL_EMAIL_ENV_NAMES = frozenset(
    {
        TERMINAL_EMAIL_ENABLED_ENV,
        TERMINAL_EMAIL_ADMIN_RECIPIENTS_ENV,
        TERMINAL_EMAIL_FROM_ENV,
        TERMINAL_EMAIL_APPLICATION_BASE_URL_ENV,
        SMTP_HOST_ENV,
        SMTP_PORT_ENV,
        SMTP_TLS_MODE_ENV,
        SMTP_USERNAME_ENV,
        SMTP_PASSWORD_ENV,
        SMTP_TIMEOUT_SECONDS_ENV,
    }
)

MAX_NOTIFICATION_EMAIL_LENGTH = 254
MAX_NOTIFICATION_RECIPIENTS = 32
MAX_SMTP_TIMEOUT_SECONDS = 10.0
DEFAULT_SMTP_TIMEOUT_SECONDS = 5.0

_EMAIL_LOCAL_PART = re.compile(r"^[A-Za-z0-9!#$%&'*+/=?^_`{|}~.-]+$")
_DNS_LABEL = re.compile(r"^[A-Za-z0-9](?:[A-Za-z0-9-]{0,61}[A-Za-z0-9])?$")


class TerminalRunNotifier(Protocol):
    """Observe one terminal run after the relevant durable commit."""

    def notify_terminal_run(
        self,
        run_id: str,
        status: RunStatus,
        *,
        include_qc: bool = False,
    ) -> None: ...


class DisabledTerminalRunNotifier:
    """Explicit no-op used when operator-owned terminal email is disabled."""

    def notify_terminal_run(
        self,
        run_id: str,
        status: RunStatus,
        *,
        include_qc: bool = False,
    ) -> None:
        del run_id, status, include_qc


class SmtpTlsMode(str, Enum):
    """The complete set of supported SMTP connection modes."""

    STARTTLS = "starttls"
    IMPLICIT_TLS = "implicit_tls"
    LOCAL_PLAINTEXT = "local_plaintext"


@dataclass(frozen=True)
class DisabledTerminalEmailSettings:
    """Explicit configuration state for a disabled notifier."""

    enabled: bool = field(default=False, init=False)


@dataclass(frozen=True)
class SmtpTerminalEmailSettings:
    """Validated operator-owned SMTP configuration for API and worker processes."""

    admin_recipients: tuple[str, ...] = field(repr=False)
    sender: str = field(repr=False)
    application_base_url: str = field(repr=False)
    smtp_host: str = field(repr=False)
    smtp_port: int
    tls_mode: SmtpTlsMode | str
    smtp_username: str | None = field(default=None, repr=False)
    smtp_password: str | None = field(default=None, repr=False)
    timeout_seconds: float = DEFAULT_SMTP_TIMEOUT_SECONDS
    enabled: bool = field(default=True, init=False)

    def __post_init__(self) -> None:
        recipients = _normalize_recipients(self.admin_recipients)
        sender = normalize_notification_email(self.sender)
        application_base_url = normalize_application_base_url(self.application_base_url)
        smtp_host = _smtp_host(self.smtp_host)
        smtp_port = _smtp_port(self.smtp_port)
        tls_mode = _tls_mode(self.tls_mode)
        username = _optional_secret(self.smtp_username, "smtp_username")
        password = _optional_secret(self.smtp_password, "smtp_password")
        if (username is None) != (password is None):
            raise ValueError(
                "smtp_username and smtp_password must be configured together"
            )
        timeout_seconds = _smtp_timeout(self.timeout_seconds)
        if tls_mode is SmtpTlsMode.LOCAL_PLAINTEXT and not _is_loopback(smtp_host):
            raise ValueError("local_plaintext SMTP requires a loopback host")
        if tls_mode is SmtpTlsMode.LOCAL_PLAINTEXT and username is not None:
            raise ValueError("local_plaintext SMTP does not accept credentials")

        object.__setattr__(self, "admin_recipients", recipients)
        object.__setattr__(self, "sender", sender)
        object.__setattr__(self, "application_base_url", application_base_url)
        object.__setattr__(self, "smtp_host", smtp_host)
        object.__setattr__(self, "smtp_port", smtp_port)
        object.__setattr__(self, "tls_mode", tls_mode)
        object.__setattr__(self, "smtp_username", username)
        object.__setattr__(self, "smtp_password", password)
        object.__setattr__(self, "timeout_seconds", timeout_seconds)


TerminalEmailSettings = DisabledTerminalEmailSettings | SmtpTerminalEmailSettings


def load_terminal_email_settings(
    environ: Mapping[str, str] | None = None,
) -> TerminalEmailSettings:
    """Load the closed operator-owned notification environment contract."""

    source = os.environ if environ is None else environ
    enabled = _environment_bool(source.get(TERMINAL_EMAIL_ENABLED_ENV, "false"))
    if not enabled:
        return DisabledTerminalEmailSettings()

    tls_mode = _required_environment(source, SMTP_TLS_MODE_ENV)
    normalized_mode = _tls_mode(tls_mode)
    default_port = 465 if normalized_mode is SmtpTlsMode.IMPLICIT_TLS else 587
    if normalized_mode is SmtpTlsMode.LOCAL_PLAINTEXT:
        default_port = 25

    recipients = tuple(
        part.strip()
        for part in _required_environment(
            source,
            TERMINAL_EMAIL_ADMIN_RECIPIENTS_ENV,
        ).split(",")
    )
    return SmtpTerminalEmailSettings(
        admin_recipients=recipients,
        sender=_required_environment(source, TERMINAL_EMAIL_FROM_ENV),
        application_base_url=_required_environment(
            source,
            TERMINAL_EMAIL_APPLICATION_BASE_URL_ENV,
        ),
        smtp_host=_required_environment(source, SMTP_HOST_ENV),
        smtp_port=_environment_int(
            source.get(SMTP_PORT_ENV, str(default_port)),
            "smtp_port",
        ),
        tls_mode=normalized_mode,
        smtp_username=source.get(SMTP_USERNAME_ENV),
        smtp_password=source.get(SMTP_PASSWORD_ENV),
        timeout_seconds=_environment_float(
            source.get(
                SMTP_TIMEOUT_SECONDS_ENV,
                str(DEFAULT_SMTP_TIMEOUT_SECONDS),
            ),
            "smtp_timeout_seconds",
        ),
    )


def normalize_notification_email(value: object) -> str:
    """Return one bounded ASCII mailbox used only as notification contact data."""

    if not isinstance(value, str):
        raise ValueError("notification email must be a string")
    if value != value.strip() or not value.isascii():
        raise ValueError("notification email must be canonical ASCII")
    if not 3 <= len(value) <= MAX_NOTIFICATION_EMAIL_LENGTH:
        raise ValueError("notification email has an invalid length")
    if any(character.isspace() or ord(character) < 32 for character in value):
        raise ValueError("notification email contains unsupported characters")
    if value.count("@") != 1:
        raise ValueError("notification email must be one mailbox")
    local_part, domain = value.rsplit("@", 1)
    if (
        not local_part
        or len(local_part) > 64
        or _EMAIL_LOCAL_PART.fullmatch(local_part) is None
        or local_part.startswith(".")
        or local_part.endswith(".")
        or ".." in local_part
        or not _valid_mail_domain(domain)
    ):
        raise ValueError("notification email must be one mailbox")
    return f"{local_part.casefold()}@{domain.casefold()}"


def normalize_application_base_url(value: object) -> str:
    """Validate the operator-provided origin used to construct safe run links."""

    if not isinstance(value, str) or not value or value != value.strip():
        raise ValueError("application_base_url must be an absolute origin")
    if not value.isascii() or any(ord(character) < 32 for character in value):
        raise ValueError("application_base_url must be an absolute origin")
    try:
        parsed = urlsplit(value)
        port = parsed.port
    except ValueError:
        raise ValueError("application_base_url must be an absolute origin") from None
    if (
        parsed.scheme not in {"http", "https"}
        or not parsed.hostname
        or parsed.username is not None
        or parsed.password is not None
        or parsed.path not in {"", "/"}
        or parsed.query
        or parsed.fragment
        or not _valid_network_host(parsed.hostname)
        or port is not None
        and not 1 <= port <= 65_535
    ):
        raise ValueError("application_base_url must be an absolute origin")
    if parsed.scheme == "http" and not _is_loopback(parsed.hostname):
        raise ValueError("HTTP application_base_url requires a loopback host")
    return urlunsplit((parsed.scheme, parsed.netloc, "", "", ""))


def _normalize_recipients(value: object) -> tuple[str, ...]:
    if not isinstance(value, tuple) or not value:
        raise ValueError("admin_recipients must be a non-empty tuple")
    normalized: list[str] = []
    seen: set[str] = set()
    for candidate in value:
        address = normalize_notification_email(candidate)
        if address not in seen:
            normalized.append(address)
            seen.add(address)
    if len(normalized) > MAX_NOTIFICATION_RECIPIENTS:
        raise ValueError(
            f"admin_recipients must contain at most {MAX_NOTIFICATION_RECIPIENTS}"
        )
    return tuple(normalized)


def _valid_mail_domain(value: str) -> bool:
    if not value or len(value) > 253 or value.startswith(".") or value.endswith("."):
        return False
    labels = value.split(".")
    return all(_DNS_LABEL.fullmatch(label) is not None for label in labels)


def _environment_bool(value: object) -> bool:
    if not isinstance(value, str):
        raise ValueError("terminal_email_enabled must be true or false")
    normalized = value.strip().casefold()
    if normalized in {"1", "true", "yes", "on"}:
        return True
    if normalized in {"0", "false", "no", "off"}:
        return False
    raise ValueError("terminal_email_enabled must be true or false")


def _required_environment(environ: Mapping[str, str], name: str) -> str:
    value = environ.get(name)
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be configured when terminal email is enabled")
    if value != value.strip() or any(character in value for character in "\x00\r\n"):
        raise ValueError(f"{name} must be canonical single-line text")
    return value


def _environment_int(value: object, name: str) -> int:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be an integer")
    try:
        normalized = int(value)
    except (TypeError, ValueError):
        raise ValueError(f"{name} must be an integer") from None
    return normalized


def _environment_float(value: object, name: str) -> float:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a finite number")
    try:
        normalized = float(value)
    except (TypeError, ValueError):
        raise ValueError(f"{name} must be a finite number") from None
    if not math.isfinite(normalized):
        raise ValueError(f"{name} must be a finite number")
    return normalized


def _smtp_port(value: object) -> int:
    if (
        isinstance(value, bool)
        or not isinstance(value, int)
        or not 1 <= value <= 65_535
    ):
        raise ValueError("smtp_port must be between 1 and 65535")
    return value


def _smtp_timeout(value: object) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError("timeout_seconds must be a positive finite number")
    normalized = float(value)
    if (
        not math.isfinite(normalized)
        or normalized <= 0
        or normalized > MAX_SMTP_TIMEOUT_SECONDS
    ):
        raise ValueError(
            f"timeout_seconds must be no greater than {MAX_SMTP_TIMEOUT_SECONDS:g}"
        )
    return normalized


def _tls_mode(value: object) -> SmtpTlsMode:
    try:
        return value if isinstance(value, SmtpTlsMode) else SmtpTlsMode(value)
    except (TypeError, ValueError):
        raise ValueError("smtp_tls_mode is invalid") from None


def _smtp_host(value: object) -> str:
    if (
        not isinstance(value, str)
        or not value
        or value != value.strip()
        or not value.isascii()
        or len(value) > 253
        or any(character.isspace() or ord(character) < 32 for character in value)
        or any(character in value for character in "/\\@")
        or not _valid_network_host(value)
    ):
        raise ValueError("smtp_host must be a bounded hostname or IP address")
    return value.casefold()


def _valid_network_host(value: str) -> bool:
    normalized = value.rstrip(".")
    if not normalized:
        return False
    try:
        ipaddress.ip_address(normalized)
        return True
    except ValueError:
        return all(
            _DNS_LABEL.fullmatch(label) is not None for label in normalized.split(".")
        )


def _optional_secret(value: object, name: str) -> str | None:
    if value is None or value == "":
        return None
    if (
        not isinstance(value, str)
        or len(value) > 4096
        or any(character in value for character in "\x00\r\n")
    ):
        raise ValueError(f"{name} must be bounded single-line text")
    return value


def _is_loopback(host: str) -> bool:
    normalized = host.rstrip(".").casefold()
    if normalized == "localhost":
        return True
    try:
        return ipaddress.ip_address(normalized).is_loopback
    except ValueError:
        return False
