"""Deployment boundary tests for terminal email operator configuration."""

from __future__ import annotations

import runpy
from pathlib import Path

import pytest

import encode_pipeline.services.process_runner as process_runner_module


TEMPLATES = (
    Path(__file__).parents[2] / "src" / "encode_pipeline" / "deployment" / "templates"
)
NOTIFICATION_ENVIRONMENT = {
    "HELIXWEAVE_TERMINAL_EMAIL_ENABLED": "true",
    "HELIXWEAVE_TERMINAL_EMAIL_ADMIN_RECIPIENTS": "admin@example.test",
    "HELIXWEAVE_TERMINAL_EMAIL_FROM": "helixweave@example.test",
    "HELIXWEAVE_TERMINAL_EMAIL_APPLICATION_BASE_URL": "https://lab.example.test",
    "HELIXWEAVE_SMTP_HOST": "smtp.example.test",
    "HELIXWEAVE_SMTP_PORT": "587",
    "HELIXWEAVE_SMTP_TLS_MODE": "starttls",
    "HELIXWEAVE_SMTP_USERNAME": "smtp-user",
    "HELIXWEAVE_SMTP_PASSWORD": "smtp-password",
    "HELIXWEAVE_SMTP_TIMEOUT_SECONDS": "5",
}


def test_api_and_worker_units_load_but_cannot_read_notification_environment() -> None:
    for unit_name in ("helixweave-api.service.in", "helixweave-worker.service.in"):
        unit = (TEMPLATES / unit_name).read_text()
        assert "EnvironmentFile=-/etc/helixweave/notifications.env" in unit
        assert "InaccessiblePaths=/etc/helixweave/notifications.env" in unit


def test_tmpfiles_owns_notifications_environment_as_root_only() -> None:
    content = (TEMPLATES / "helixweave.tmpfiles.conf").read_text()

    assert "f /etc/helixweave/notifications.env 0600 root root -" in content
    assert "z /etc/helixweave/notifications.env 0600 root root -" in content


def test_candidate_launcher_forwards_notifications_only_to_api_and_worker(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    namespace = runpy.run_path(str(TEMPLATES / "helixweave-service"))
    with monkeypatch.context() as scoped:
        for name, value in NOTIFICATION_ENVIRONMENT.items():
            scoped.setenv(name, value)

        for command in ("api", "worker"):
            selected = namespace["_environment"](command)
            assert all(
                selected[name] == value
                for name, value in NOTIFICATION_ENVIRONMENT.items()
            )

        for command in (
            "db-prepare",
            "operator-action",
            "encode-runtime-prepare",
            "bulk-runtime-prepare",
        ):
            selected = namespace["_environment"](command)
            assert set(selected).isdisjoint(NOTIFICATION_ENVIRONMENT)


def test_candidate_launcher_rejects_multiline_notification_values(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    namespace = runpy.run_path(str(TEMPLATES / "helixweave-service"))
    with monkeypatch.context() as scoped:
        scoped.setenv("HELIXWEAVE_SMTP_PASSWORD", "private\nvalue")
        with pytest.raises(ValueError):
            namespace["_environment"]("worker")


def test_scientific_process_does_not_inherit_notification_environment(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    for name, value in NOTIFICATION_ENVIRONMENT.items():
        monkeypatch.setenv(name, value)

    result = process_runner_module._subprocess_environment({})

    assert result.is_success is True
    assert set(result.value).isdisjoint(NOTIFICATION_ENVIRONMENT)
