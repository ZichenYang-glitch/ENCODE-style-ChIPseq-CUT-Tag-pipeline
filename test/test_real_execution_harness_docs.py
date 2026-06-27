"""Guard that the real-execution harness policy doc stays in sync.

This test does not run any real execution, containers, or heavy compute.
It only checks that the policy document exists and references the concrete
entry points and legacy classifications it claims to describe.
"""

import sys
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parent.parent
DOCS_PATH = REPO_ROOT / "docs" / "development" / "real-execution-harness.md"

sys.path.insert(0, str(REPO_ROOT / "test"))
from test_stage_shim import LEGACY_STAGE_CLASSIFICATION


REQUIRED_MENTIONS = (
    "test/test_stage8b_tiny_execution.py",
    "workflow_dispatch",
    "scripts/smoke_container_runner.sh",
    "profiles/default",
    "profiles/hpc",
    "docs/container-usage.md",
)


def test_real_execution_harness_doc_exists():
    assert DOCS_PATH.is_file(), f"Missing harness policy doc: {DOCS_PATH}"


def test_real_execution_harness_doc_mentions_required_entry_points():
    text = DOCS_PATH.read_text(encoding="utf-8")
    missing = [m for m in REQUIRED_MENTIONS if m not in text]
    assert not missing, f"Doc missing required references: {missing}"


def test_real_execution_only_scripts_are_mentioned():
    """Every real-execution-only legacy script must appear in the policy doc."""
    text = DOCS_PATH.read_text(encoding="utf-8")
    real_exec_scripts = [
        name
        for name, meta in LEGACY_STAGE_CLASSIFICATION.items()
        if meta.get("category") == "real-execution-only"
    ]
    missing = [s for s in real_exec_scripts if s not in text]
    assert not missing, f"Doc must mention each real-execution-only script: {missing}"


def test_manual_integration_count_is_documented():
    """The doc must state the current number of manual-integration scripts."""
    text = DOCS_PATH.read_text(encoding="utf-8")
    manual_integration = [
        name
        for name, meta in LEGACY_STAGE_CLASSIFICATION.items()
        if meta.get("category") == "manual-integration"
    ]
    expected_count = len(manual_integration)
    expected_phrase = f"manual-integration scripts ({expected_count})"
    assert expected_phrase in text, (
        f"Doc must contain the phrase {expected_phrase!r} so the count stays in sync"
    )
