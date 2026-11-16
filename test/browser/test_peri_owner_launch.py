"""The browser runtime admits one fixed Peri owner and does not choose a fallback."""

from __future__ import annotations

import os
from pathlib import Path

import pytest

import platform_runtime


def test_absent_peri_binary_does_not_start_an_owner(tmp_path: Path) -> None:
    environment: dict[str, str] = {}

    platform_runtime.start_admitted_peri_owner(tmp_path, environment)

    assert "HELIXWEAVE_PERI_OWNER_SOCKET" not in environment
    assert not (tmp_path / "peri-owner.pid").exists()


def test_long_owner_socket_path_is_rejected_before_launch(monkeypatch) -> None:
    monkeypatch.setattr(
        platform_runtime.subprocess,
        "Popen",
        lambda *_args, **_kwargs: pytest.fail("owner launch must not start"),
    )
    root = Path("/tmp") / ("peri-owner-path-" + "x" * 120)

    with pytest.raises(RuntimeError, match="socket path"):
        platform_runtime.start_admitted_peri_owner(
            root, {"HELIXWEAVE_PERI_BINARY": "peri"}
        )


def test_rejected_peri_binary_does_not_publish_a_socket(
    tmp_path: Path, monkeypatch
) -> None:
    class Exited:
        pid = 4321

        def poll(self) -> int:
            return 1

    monkeypatch.setattr(
        platform_runtime.subprocess, "Popen", lambda *_args, **_kwargs: Exited()
    )
    environment = {"HELIXWEAVE_PERI_BINARY": str(tmp_path / "not-peri")}

    with pytest.raises(RuntimeError, match="not admitted"):
        platform_runtime.start_admitted_peri_owner(tmp_path, environment)

    assert "HELIXWEAVE_PERI_OWNER_SOCKET" not in environment
    assert not (tmp_path / "peri-owner.pid").exists()


def test_admitted_owner_records_the_socket_it_opened(
    tmp_path: Path, monkeypatch
) -> None:
    class Running:
        pid = 5151

        def poll(self) -> None:
            return None

    def popen(_argv: list[str], **kwargs: object) -> Running:
        env = kwargs["env"]
        assert isinstance(env, dict)
        socket_path = Path(str(env["HELIXWEAVE_PERI_OWNER_SOCKET"]))
        socket_path.write_text("", encoding="utf-8")
        return Running()

    monkeypatch.setattr(platform_runtime.subprocess, "Popen", popen)
    environment = {"HELIXWEAVE_PERI_BINARY": str(tmp_path / "peri")}

    platform_runtime.start_admitted_peri_owner(tmp_path, environment)

    assert environment["HELIXWEAVE_PERI_OWNER_SOCKET"].endswith("/peri-owner.sock")
    assert (tmp_path / "peri-owner.pid").read_text(encoding="utf-8") == "5151"


def test_owner_cleanup_does_not_signal_an_unrelated_pid(tmp_path: Path) -> None:
    (tmp_path / "peri-owner.pid").write_text(str(os.getpid()), encoding="utf-8")

    platform_runtime.stop_peri_owner(tmp_path)

    os.kill(os.getpid(), 0)
