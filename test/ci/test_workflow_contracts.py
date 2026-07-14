"""Durable contracts for the CI tier and coverage-ratchet topology."""

from pathlib import Path
import re

import yaml


REPO_ROOT = Path(__file__).resolve().parents[2]
WORKFLOW_DIR = REPO_ROOT / ".github" / "workflows"
EVENT_SHA = "${{ github.sha }}"
PR_HEAD_SHA = "${{ github.event.pull_request.head.sha }}"
PR_BASE_SHA = "${{ github.event.pull_request.base.sha }}"


def _load(name: str) -> dict:
    return yaml.safe_load((WORKFLOW_DIR / name).read_text(encoding="utf-8"))


def _runs(job: dict) -> str:
    return "\n".join(str(step.get("run", "")) for step in job["steps"])


def _checkout_steps(workflow: dict):
    for job in workflow["jobs"].values():
        for step in job["steps"]:
            if str(step.get("uses", "")).startswith("actions/checkout@"):
                yield step


def _checkout_jobs(workflow: dict):
    for job_id, job in workflow["jobs"].items():
        if any(
            str(step.get("uses", "")).startswith("actions/checkout@")
            for step in job["steps"]
        ):
            yield job_id, job


def test_all_workflows_parse_and_required_jobs_keep_stable_ids():
    workflows = {
        path.name: yaml.safe_load(path.read_text(encoding="utf-8"))
        for path in sorted(WORKFLOW_DIR.glob("*.yml"))
    }
    assert workflows
    assert {
        "fast-checks",
        "coverage",
        "frontend",
        "browser-e2e",
        "platform-real-execution",
        "real-execution",
        "container-smoke",
    } <= _load("ci.yml")["jobs"].keys()
    assert "lint" in _load("lint.yml")["jobs"]
    assert "lock-check" in _load("lock-check.yml")["jobs"]


def test_ci_event_matrix_and_concurrency_are_tier_aware():
    workflow = _load("ci.yml")
    assert set(workflow["on"]) == {
        "pull_request",
        "push",
        "workflow_dispatch",
        "schedule",
        "release",
    }
    assert workflow["on"]["push"]["branches"] == ["main"]
    assert workflow["on"]["release"]["types"] == ["published"]
    assert workflow["on"]["schedule"]
    assert "github.event_name" in workflow["concurrency"]["group"]


def test_pull_request_jobs_checkout_the_merge_result_and_report_pr_identity():
    for name in ("ci.yml", "lint.yml", "lock-check.yml"):
        workflow = _load(name)
        assert "pull_request" in workflow["on"]
        jobs = list(_checkout_jobs(workflow))
        assert jobs, name
        for job_id, job in jobs:
            checkout = next(_checkout_steps({"jobs": {job_id: job}}))
            assert checkout["with"]["fetch-depth"] == 0
            # On pull_request, github.sha is GitHub's synthetic merge commit.
            assert checkout["with"]["ref"] == EVENT_SHA

            identity = next(
                step
                for step in job["steps"]
                if str(step.get("name", "")).startswith("Confirm exact checkout")
            )
            assert identity["env"]["EXPECTED_SHA"] == EVENT_SHA
            assert identity["env"]["PR_HEAD_SHA"] == PR_HEAD_SHA
            assert identity["env"]["PR_BASE_SHA"] == PR_BASE_SHA
            assert "tested-sha=" in identity["run"]
            assert "pr-head-sha=" in identity["run"]
            assert "pr-base-sha=" in identity["run"]
            assert "Tested SHA:" in identity["run"]
            assert "PR head SHA:" in identity["run"]
            assert "PR base SHA:" in identity["run"]


def test_non_pr_jobs_checkout_the_exact_event_sha_without_a_head_fallback():
    for name in ("ci.yml", "lint.yml", "lock-check.yml"):
        workflow = _load(name)
        assert "workflow_dispatch" in workflow["on"]
        assert "push" in workflow["on"]
        for checkout in _checkout_steps(workflow):
            # On push, dispatch, schedule, and release, github.sha is the exact
            # event commit; no pull-request-head fallback may override it.
            assert checkout["with"]["ref"] == EVENT_SHA
            assert "pull_request.head.sha" not in checkout["with"]["ref"]


def test_fast_checks_is_the_only_deterministic_pytest_coverage_producer():
    jobs = _load("ci.yml")["jobs"]
    producer = _runs(jobs["fast-checks"])
    consumer = _runs(jobs["coverage"])

    assert producer.count("python3 -m pytest test") == 1
    assert (
        "not full_main and not platform_real_execution and not real_execution"
        in producer
    )
    assert "not platform_real_execution and not real_execution" in producer
    assert "--cov-fail-under=0" in producer
    assert "--junitxml=pytest-report.xml" in producer
    assert "check_junit_outcomes.py pytest-report.xml" in producer
    assert "pytest" not in consumer
    assert jobs["coverage"]["needs"] == "fast-checks"
    assert "always()" in jobs["coverage"]["if"]


def test_coverage_artifact_and_ratchets_are_stable_and_nonduplicative():
    jobs = _load("ci.yml")["jobs"]
    producer_steps = jobs["fast-checks"]["steps"]
    upload = next(
        step
        for step in producer_steps
        if step.get("uses") == "actions/upload-artifact@v4"
    )
    settings = upload["with"]
    assert settings["name"] == "python-coverage-${{ github.run_id }}"
    assert "run_attempt" not in settings["name"]
    assert settings["include-hidden-files"] is True
    assert settings["overwrite"] is True
    assert settings["if-no-files-found"] == "error"
    for artifact in (".coverage", "coverage.xml", "coverage.json", "pytest-report.xml"):
        assert artifact in settings["path"]

    coverage_runs = _runs(jobs["coverage"])
    assert "--fail-under=80" in coverage_runs
    assert "coverage report --fail-under=82" in coverage_runs
    for floor in ("88.45", "87.28", "89.06", "82.37"):
        assert f"--fail-under={floor}" in coverage_runs

    config = (REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8")
    coverage_report = re.search(
        r"(?ms)^\[tool\.coverage\.report\]\n(.*?)(?=^\[|\Z)", config
    )
    assert coverage_report
    assert re.search(r"(?m)^fail_under\s*=\s*82\s*$", coverage_report.group(1))


def test_all_pytest_tiers_enforce_zero_skip_junit_outcomes():
    jobs = _load("ci.yml")["jobs"]
    for job_id, report in (
        ("fast-checks", "pytest-report.xml"),
        ("platform-real-execution", "platform-real-report.xml"),
        ("real-execution", "scientific-real-report.xml"),
    ):
        runs = _runs(jobs[job_id])
        assert f"--junitxml={report}" in runs
        assert f"check_junit_outcomes.py {report}" in runs

    assert "xfail_strict = true" in (REPO_ROOT / "pyproject.toml").read_text(
        encoding="utf-8"
    )


def test_real_execution_jobs_are_non_pr_conditional_and_fail_closed():
    jobs = _load("ci.yml")["jobs"]
    for job_id in ("platform-real-execution", "real-execution", "container-smoke"):
        condition = jobs[job_id]["if"]
        assert "workflow_dispatch" in condition
        assert "schedule" in condition
        assert "release" in condition
        assert "pull_request" not in condition
        assert "push" not in condition

    platform_runs = _runs(jobs["platform-real-execution"])
    scientific_runs = _runs(jobs["real-execution"])
    container_runs = _runs(jobs["container-smoke"])
    assert "-m platform_real_execution test/workers" in platform_runs
    assert "-m real_execution test/real_execution" in scientific_runs
    assert "HELIXWEAVE_REQUIRE_REAL_EXECUTION" in str(jobs["real-execution"])
    assert "docker build" in container_runs
    assert "scripts/smoke_container_runner.sh" in container_runs


def test_lint_and_lock_workflows_cover_the_maintained_contracts():
    lint = _load("lint.yml")
    lock = _load("lock-check.yml")
    lint_runs = _runs(lint["jobs"]["lint"])
    lock_runs = _runs(lock["jobs"]["lock-check"])

    maintained = "src scripts test containers workflow/lib"
    assert f"ruff check {maintained}" in lint_runs
    assert f"ruff format --check {maintained}" in lint_runs
    assert "snakefmt --check workflow/" in lint_runs
    assert "test/check_snakemake_lint.py" in lint_runs

    assert set(lock["on"]) == {"pull_request", "push", "workflow_dispatch"}
    assert "github.event_name" in lock["concurrency"]["group"]
    assert "github.event.pull_request.base.sha" in str(lock)
    assert "github.event.before" in str(lock)
    assert "git diff --name-status -z --find-renames" in lock_runs
    assert "':(glob)workflow/envs/*.yml'" in lock_runs
    assert "require_deleted_lock" in lock_runs
    assert "Push base is empty or all-zero" in lock_runs
    assert "Pull-request base SHA is unexpectedly empty" in lock_runs
