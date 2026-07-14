"""Guard that the real-execution harness policy doc stays in sync.

This test does not run any real execution, containers, or heavy compute.
It only checks that the policy document references its concrete entry points.
"""

from pathlib import Path

import pytest
import yaml

from _tool_resolver import require_external_tools


REPO_ROOT = Path(__file__).resolve().parent.parent
DOCS_PATH = REPO_ROOT / "docs" / "development" / "real-execution-harness.md"
CI_PATH = REPO_ROOT / ".github" / "workflows" / "ci.yml"
CHIPSEQ_SPEC_PATH = REPO_ROOT / "workflow" / "envs" / "chipseq.yml"
CHIPSEQ_LOCK_PATH = REPO_ROOT / "workflow" / "envs" / "chipseq.lock"


REQUIRED_MENTIONS = (
    "test/real_execution/test_scientific_tiny_preprocessing.py",
    "test/real_execution/test_pseudoreplicate_splitting.py",
    "test/real_execution/test_cuttag_fragment_size.py",
    "test/real_execution",
    "HELIXWEAVE_REQUIRE_REAL_EXECUTION",
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


def test_ci_runs_the_complete_real_execution_tier_fail_closed():
    workflow = yaml.safe_load(CI_PATH.read_text(encoding="utf-8"))
    steps = workflow["jobs"]["real-execution"]["steps"]
    execution_step = next(step for step in steps if "run" in step)

    assert execution_step["env"]["HELIXWEAVE_REQUIRE_REAL_EXECUTION"] == "1"
    assert "python3 -m pytest -m real_execution" in execution_step["run"]
    assert "test/real_execution -v" in execution_step["run"]
    assert "test_scientific_tiny_preprocessing.py" not in execution_step["run"]


def test_chipseq_exact_environment_contains_pytest_and_its_dependencies():
    spec = yaml.safe_load(CHIPSEQ_SPEC_PATH.read_text(encoding="utf-8"))
    assert "pytest >=8,<10" in spec["dependencies"]

    lock = CHIPSEQ_LOCK_PATH.read_text(encoding="utf-8")
    for package in ("iniconfig", "pluggy", "pytest"):
        assert f"/noarch/{package}-" in lock


def test_missing_real_tools_skip_locally_but_fail_a_required_gate(monkeypatch):
    monkeypatch.delenv("HELIXWEAVE_REQUIRE_REAL_EXECUTION", raising=False)
    with pytest.raises(pytest.skip.Exception, match="samtools"):
        require_external_tools(["samtools"])

    monkeypatch.setenv("HELIXWEAVE_REQUIRE_REAL_EXECUTION", "1")
    with pytest.raises(pytest.fail.Exception, match="samtools"):
        require_external_tools(["samtools"])
