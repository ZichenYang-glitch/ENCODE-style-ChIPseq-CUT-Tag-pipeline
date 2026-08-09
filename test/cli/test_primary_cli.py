"""Primary HelixWeave CLI compatibility contracts."""

from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys

import pytest

from encode_pipeline.cli import app
from encode_pipeline.cli import local_platform


REPO_ROOT = Path(__file__).resolve().parents[2]


def test_primary_cli_delegates_existing_arguments_to_local_platform(
    monkeypatch,
) -> None:
    observed: list[list[str]] = []
    monkeypatch.setattr(
        local_platform,
        "main",
        lambda arguments: observed.append(list(arguments)) or 17,
    )

    assert app.main(["--doctor"]) == 17
    assert observed == [["--doctor"]]


def test_primary_cli_dispatches_only_the_admin_namespace(monkeypatch) -> None:
    observed: list[list[str]] = []
    monkeypatch.setattr(
        app.admin,
        "main",
        lambda arguments: observed.append(list(arguments)) or 19,
    )

    assert app.main(["admin", "--database-url", "sqlite:////tmp/platform.db"]) == 19
    assert observed == [["--database-url", "sqlite:////tmp/platform.db"]]


def test_primary_cli_rejects_unknown_commands(capsys) -> None:
    with pytest.raises(SystemExit, match="2"):
        app.main(["unknown"])

    assert "unrecognized arguments: unknown" in capsys.readouterr().err


def test_primary_cli_help_covers_doctor_and_demo(capsys) -> None:
    with pytest.raises(SystemExit, match="0"):
        app.main(["--help"])

    output = capsys.readouterr().out
    assert "HelixWeave" in output
    assert "--doctor" in output
    assert "--input-authoring-demo" in output


def test_primary_cli_exposes_the_admin_run_recovery_commands(capsys) -> None:
    with pytest.raises(SystemExit, match="0"):
        app.main(
            [
                "admin",
                "--database-url",
                "sqlite:////tmp/platform.db",
                "run",
                "--help",
            ]
        )

    output = capsys.readouterr().out
    assert "diagnose" in output
    assert "fail" in output
    assert "requeue" in output


def test_module_and_compatibility_script_share_exact_help(tmp_path: Path) -> None:
    environment = {
        **os.environ,
        "PYTHONPATH": str(REPO_ROOT / "src"),
        "PYTHONNOUSERSITE": "1",
        "PYTHONDONTWRITEBYTECODE": "1",
    }
    module = subprocess.run(
        [sys.executable, "-m", "encode_pipeline", "--help"],
        cwd=tmp_path,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )
    script = subprocess.run(
        [
            sys.executable,
            "-I",
            "-S",
            str(REPO_ROOT / "scripts/checkout_bootstrap.py"),
            "--repository-root",
            str(REPO_ROOT),
            "local-platform",
            "--help",
        ],
        cwd=tmp_path,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )

    assert module.returncode == 0, module.stderr
    assert script.returncode == 0, script.stderr
    assert module.stdout == script.stdout
    assert module.stdout.startswith("usage: helixweave ")
