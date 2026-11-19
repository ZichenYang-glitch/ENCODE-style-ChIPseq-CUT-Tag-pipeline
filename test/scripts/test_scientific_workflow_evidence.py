"""Tests for the scientific-workflow evidence ledger checker."""

from __future__ import annotations

import importlib.util
from pathlib import Path

_SCRIPT = (
    Path(__file__).resolve().parents[2]
    / "scripts/check_scientific_workflow_evidence.py"
)
_SPEC = importlib.util.spec_from_file_location(
    "check_scientific_workflow_evidence",
    _SCRIPT,
)
checker = importlib.util.module_from_spec(_SPEC)
assert _SPEC.loader is not None
_SPEC.loader.exec_module(checker)

LEDGER = checker.ROOT / "docs/development/scientific-workflow-evidence-ledger.md"
QUALIFICATION = (
    checker.ROOT
    / "src/encode_pipeline/contracts/nfcore_rnaseq"
    / "default-execution-qualification-1.1.0.json"
)


def test_repository_ledger_passes():
    failures, _summary = checker.evaluate(checker.ROOT, None)
    assert failures == []
    assert checker.main() == 0


def test_missing_current_manifest_citation_fails():
    ledger = LEDGER.read_text(encoding="utf-8")
    broken = ledger.replace("current-manifest-file: ", "dropped-manifest-file: ", 1)
    failures, _summary = checker.evaluate(checker.ROOT, broken)
    assert "ledger is missing current-manifest-file" in failures
    assert checker.main(broken) == 1


def test_historical_identity_rejects_the_current_file_hash():
    ledger = LEDGER.read_text(encoding="utf-8")
    historical = checker.HISTORICAL_IDENTITY["historical-qualification-file"]
    current = checker.digest(QUALIFICATION)
    assert historical != current
    old_line = f"historical-qualification-file: {historical}"
    new_line = f"historical-qualification-file: {current}"
    broken = ledger.replace(old_line, new_line, 1)
    broken = broken.replace(checker.NOT_RERUN, checker.RERUN_CLAIM, 1)
    failures, _summary = checker.evaluate(checker.ROOT, broken)
    assert (
        "ledger historical-qualification-file does not match "
        "the historical evidence identity"
    ) in failures
    assert (
        "ledger claims historical evidence was rerun on the current software identity"
    ) in failures
    assert checker.main(broken) == 1
