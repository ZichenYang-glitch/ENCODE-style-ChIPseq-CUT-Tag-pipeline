"""Drift-gate tests for the frontend generated client."""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[2]
FRONTEND_DIR = REPO_ROOT / "frontend"
GENERATED_DIR = FRONTEND_DIR / "src" / "api" / "generated"


@pytest.fixture
def empty_generated_index():
    """Yield while guaranteeing no stray untracked files exist under generated/."""
    _assert_no_untracked_generated_files()
    try:
        yield
    finally:
        _assert_no_untracked_generated_files()


def _assert_no_untracked_generated_files() -> None:
    result = subprocess.run(
        [
            "git",
            "ls-files",
            "--others",
            "--exclude-standard",
            "--",
            str(GENERATED_DIR),
        ],
        check=True,
        capture_output=True,
        text=True,
        cwd=REPO_ROOT,
    )
    assert result.stdout.strip() == "", "stray untracked generated files detected"


def test_untracked_generated_file_fails_drift_gate(empty_generated_index):
    """A newly generated but untracked file under src/api/generated/ must fail the drift gate.

    This mirrors the CI check and ensures the gate cannot pass with stale/orphaned
    generated artifacts.
    """
    probe = GENERATED_DIR / "probe_untracked_drift.ts"
    assert not probe.exists()

    try:
        probe.write_text("// intentionally untracked drift probe\n")

        ls_result = subprocess.run(
            [
                "git",
                "ls-files",
                "--others",
                "--exclude-standard",
                "--",
                str(GENERATED_DIR),
            ],
            check=True,
            capture_output=True,
            text=True,
            cwd=REPO_ROOT,
        )
        assert "probe_untracked_drift.ts" in ls_result.stdout

        gate_result = subprocess.run(
            "git diff --exit-code -- openapi.json src/api/generated/ && "
            'test -z "$(git ls-files --others --exclude-standard -- src/api/generated/)"',
            shell=True,
            cwd=FRONTEND_DIR,
            capture_output=True,
            text=True,
        )
        assert gate_result.returncode != 0, (
            "drift gate must fail when an untracked generated file exists"
        )
    finally:
        probe.unlink(missing_ok=True)
