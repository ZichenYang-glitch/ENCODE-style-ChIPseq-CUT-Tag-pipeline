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
PYPROJECT_PATH = REPO_ROOT / "pyproject.toml"


REQUIRED_MENTIONS = (
    "platform-real-execution",
    "platform_real_execution",
    "test/workers/test_redis_process_integration.py",
    "test/workers/test_tiny_execution_e2e.py",
    "test/workers/test_cancellation_e2e.py",
    "test/real_execution/test_scientific_tiny_preprocessing.py",
    "test/real_execution/test_pseudoreplicate_splitting.py",
    "test/real_execution/test_cuttag_fragment_size.py",
    "test/real_execution",
    "HELIXWEAVE_REQUIRE_REAL_EXECUTION",
    "workflow_dispatch",
    "nightly",
    "published",
    "container-smoke",
    "scripts/smoke_container_runner.sh",
    "profiles/default",
    "profiles/hpc",
    "docs/container-usage.md",
)


def _run_text(job):
    return "\n".join(str(step["run"]) for step in job["steps"] if "run" in step)


def _step_with(job, fragment):
    return next(step for step in job["steps"] if fragment in str(step.get("run", "")))


def test_real_execution_harness_doc_exists():
    assert DOCS_PATH.is_file(), f"Missing harness policy doc: {DOCS_PATH}"


def test_real_execution_harness_doc_mentions_required_entry_points():
    text = DOCS_PATH.read_text(encoding="utf-8")
    missing = [m for m in REQUIRED_MENTIONS if m not in text]
    assert not missing, f"Doc missing required references: {missing}"


def test_ci_defines_complete_real_execution_tiers():
    workflow = yaml.safe_load(CI_PATH.read_text(encoding="utf-8"))
    jobs = workflow["jobs"]
    assert {
        "platform-real-execution",
        "real-execution",
        "container-smoke",
    } <= jobs.keys()

    platform = jobs["platform-real-execution"]
    platform_text = _run_text(platform)
    assert "pytest" in platform_text
    assert "platform_real_execution" in platform_text
    assert "test/workers" in platform_text

    scientific = jobs["real-execution"]
    scientific_text = _run_text(scientific)
    assert "pytest" in scientific_text
    assert "real_execution" in scientific_text
    assert "test/real_execution" in scientific_text
    assert "test_scientific_tiny_preprocessing.py" not in scientific_text

    container_text = _run_text(jobs["container-smoke"])
    assert "containers/Dockerfile.runner" in container_text
    assert "scripts/smoke_container_runner.sh" in container_text

    for job in (platform, scientific, jobs["container-smoke"]):
        condition = str(job["if"])
        assert "workflow_dispatch" in condition
        assert "schedule" in condition
        assert "release" in condition
        assert "pull_request" not in condition


def test_real_pytest_jobs_provide_prerequisites_and_junit():
    workflow = yaml.safe_load(CI_PATH.read_text(encoding="utf-8"))
    jobs = workflow["jobs"]

    platform_job = jobs["platform-real-execution"]
    platform_step = _step_with(platform_job, "platform_real_execution")
    platform_env = {**platform_job.get("env", {}), **platform_step.get("env", {})}
    assert platform_env["ENCODE_PIPELINE_TEST_REDIS_URL"].startswith(
        "redis://127.0.0.1:6379/"
    )
    assert "--junitxml" in platform_step["run"]

    scientific_job = jobs["real-execution"]
    scientific_step = _step_with(scientific_job, "real_execution")
    scientific_env = {
        **scientific_job.get("env", {}),
        **scientific_step.get("env", {}),
    }
    assert scientific_env["HELIXWEAVE_REQUIRE_REAL_EXECUTION"] == "1"
    assert "--junitxml" in scientific_step["run"]


def test_pytest_marker_contract_matches_documented_tiers():
    config = PYPROJECT_PATH.read_text(encoding="utf-8")

    for marker in ("full_main", "platform_real_execution", "real_execution"):
        assert f'"{marker}:' in config
    assert '"--strict-markers"' in config
    assert "not real_execution and not platform_real_execution" in config
    assert "xfail_strict = true" in config


def test_chipseq_exact_environment_contains_pytest_and_its_dependencies():
    spec = yaml.safe_load(CHIPSEQ_SPEC_PATH.read_text(encoding="utf-8"))
    assert "pytest >=8,<10" in spec["dependencies"]

    lock = CHIPSEQ_LOCK_PATH.read_text(encoding="utf-8")
    for package in ("iniconfig", "pluggy", "pytest"):
        assert f"/noarch/{package}-" in lock


def test_missing_real_tools_skip_locally_but_fail_the_ci_gate(monkeypatch):
    monkeypatch.delenv("HELIXWEAVE_REQUIRE_REAL_EXECUTION", raising=False)
    with pytest.raises(pytest.skip.Exception, match="samtools"):
        require_external_tools(["samtools"])

    monkeypatch.setenv("HELIXWEAVE_REQUIRE_REAL_EXECUTION", "1")
    with pytest.raises(pytest.fail.Exception, match="samtools"):
        require_external_tools(["samtools"])
