"""Durable contracts for the CI tier and coverage-ratchet topology."""

import copy
import os
from pathlib import Path
import re
import socket
import subprocess
import sys
import tempfile

from coverage import Coverage
import pytest
import yaml


REPO_ROOT = Path(__file__).resolve().parents[2]
WORKFLOW_DIR = REPO_ROOT / ".github" / "workflows"
EVENT_SHA = "${{ github.sha }}"
PR_HEAD_SHA = "${{ github.event.pull_request.head.sha }}"
PR_BASE_SHA = "${{ github.event.pull_request.base.sha }}"
REDIS_SECRET_EXPR = "${{ secrets.ENCODE_PIPELINE_TEST_REDIS_URL }}"
REDIS_VAR_EXPR = "${{ vars.ENCODE_PIPELINE_TEST_REDIS_URL }}"
REDIS_SECRET_STEPS = (
    "Require runner-owned runtime coordinates",
    "Run rapid, STAR Salmon, and lifecycle acceptance",
    "Run public registry-to-results Bulk RNA-seq product Gate",
)
HARNESS_DOC = REPO_ROOT / "docs" / "development" / "real-execution-harness.md"


def _load(name: str) -> dict:
    return yaml.safe_load((WORKFLOW_DIR / name).read_text(encoding="utf-8"))


def _runs(job: dict) -> str:
    return "\n".join(str(step.get("run", "")) for step in job["steps"])


def _checkout_steps(workflow: dict):
    for job in workflow["jobs"].values():
        for step in job["steps"]:
            if str(step.get("uses", "")).startswith("actions/checkout@"):
                yield step


@pytest.mark.parametrize(
    ("uid", "gid", "expected", "succeeds"),
    (
        ("", "", "", True),
        (
            "0",
            "0",
            "ENCODE_PIPELINE_BULK_CONTAINER_UID=0\nENCODE_PIPELINE_BULK_CONTAINER_GID=0\n",
            True,
        ),
        (
            "1234",
            "5678",
            "ENCODE_PIPELINE_BULK_CONTAINER_UID=1234\nENCODE_PIPELINE_BULK_CONTAINER_GID=5678\n",
            True,
        ),
        ("0", "", "ENCODE_PIPELINE_BULK_CONTAINER_UID=0\n", True),
        ("0", "1\nINJECTED=1", "", False),
    ),
)
def test_protected_bulk_step_exports_validated_optional_container_identity(
    tmp_path: Path, uid: str, gid: str, expected: str, succeeds: bool
):
    steps = _load("ci.yml")["jobs"]["bulk-rnaseq-real-execution"]["steps"]
    step = next(
        s for s in steps if s.get("name") == "Bind optional Bulk container identity"
    )
    assert step["env"] == {
        "CONFIGURED_BULK_CONTAINER_UID": "${{ vars.ENCODE_PIPELINE_BULK_CONTAINER_UID }}",
        "CONFIGURED_BULK_CONTAINER_GID": "${{ vars.ENCODE_PIPELINE_BULK_CONTAINER_GID }}",
    }
    names = [s.get("name") for s in steps]
    assert names.index("Install platform package") < names.index(step["name"])
    assert names.index(step["name"]) < names.index(
        "Run rapid, STAR Salmon, and lifecycle acceptance"
    )
    github_env = tmp_path / "github-env"
    github_env.write_text("", encoding="utf-8")
    environment = {
        **os.environ,
        "PATH": f"{Path(sys.executable).parent}{os.pathsep}{os.defpath}",
        "PYTHONPATH": str(REPO_ROOT / "src"),
        "GITHUB_ENV": str(github_env),
        "CONFIGURED_BULK_CONTAINER_UID": uid,
        "CONFIGURED_BULK_CONTAINER_GID": gid,
        "ENCODE_PIPELINE_BULK_CONTAINER_UID": "42",
        "ENCODE_PIPELINE_BULK_CONTAINER_GID": "43",
    }
    result = subprocess.run(
        ["bash", "-euo", "pipefail", "-c", step["run"]],
        cwd=REPO_ROOT,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )
    assert (result.returncode == 0) is succeeds, result.stderr
    assert github_env.read_text(encoding="utf-8") == expected
    if not succeeds:
        assert "ENCODE_PIPELINE_BULK_CONTAINER_GID" in result.stderr
        assert "INJECTED" not in result.stderr


def _runner_coordinate_step() -> dict:
    steps = _load("ci.yml")["jobs"]["bulk-rnaseq-real-execution"]["steps"]
    return next(
        step
        for step in steps
        if step.get("name") == "Require runner-owned runtime coordinates"
    )


def _run_coordinate_step(
    root: Path, redis_url: str
) -> subprocess.CompletedProcess[str]:
    runtime = root / "runtime"
    runtime.mkdir()
    fixture = root / "fixture.json"
    fixture.write_text("{}\n", encoding="utf-8")
    docker = root / "docker"
    docker.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
    docker.chmod(0o755)
    runner_temp = root / "runner-temp"
    runner_temp.mkdir()
    socket_dir = Path(tempfile.mkdtemp(prefix="hwc", dir="/tmp"))
    endpoint = socket.socket(socket.AF_UNIX)
    endpoint.bind(str(socket_dir / "d.sock"))
    try:
        environment = {
            "PATH": os.defpath,
            "HELIXWEAVE_BULK_RNASEQ_RUNTIME_ROOT": str(runtime),
            "HELIXWEAVE_BULK_RNASEQ_FIXTURE_MANIFEST": str(fixture),
            "ENCODE_PIPELINE_TEST_REDIS_URL": redis_url,
            "ENCODE_PIPELINE_MANAGED_DOCKER_EXECUTABLE": str(docker),
            "ENCODE_PIPELINE_MANAGED_DOCKER_SOCKET": str(socket_dir / "d.sock"),
            "RUNNER_TEMP": str(runner_temp),
        }
        return subprocess.run(
            ["bash", "-euo", "pipefail", "-c", _runner_coordinate_step()["run"]],
            cwd=REPO_ROOT,
            env=environment,
            capture_output=True,
            text=True,
            check=False,
        )
    finally:
        endpoint.close()
        (socket_dir / "d.sock").unlink(missing_ok=True)
        socket_dir.rmdir()


def _redis_env_refs(job: dict) -> list[tuple[str, str, str]]:
    refs = []
    for step in job["steps"]:
        for key, value in (step.get("env") or {}).items():
            if key in {
                "ENCODE_PIPELINE_TEST_REDIS_URL",
                "ENCODE_PIPELINE_E2E_REDIS_URL",
            }:
                refs.append((step.get("name"), key, str(value)))
    return refs


def assert_bulk_redis_secret_contract(job: dict) -> None:
    """Redis is one environment secret on three steps, never a variable."""

    assert set(job.get("env", {})) == {"PYTHONDONTWRITEBYTECODE"}
    assert REDIS_VAR_EXPR not in str(job)
    assert "secrets.ENCODE_PIPELINE_TEST_REDIS_URL ||" not in str(job)
    assert str(job).count(REDIS_SECRET_EXPR) == 4
    assert _redis_env_refs(job) == [
        (REDIS_SECRET_STEPS[0], "ENCODE_PIPELINE_TEST_REDIS_URL", REDIS_SECRET_EXPR),
        (REDIS_SECRET_STEPS[1], "ENCODE_PIPELINE_TEST_REDIS_URL", REDIS_SECRET_EXPR),
        (REDIS_SECRET_STEPS[2], "ENCODE_PIPELINE_TEST_REDIS_URL", REDIS_SECRET_EXPR),
        (REDIS_SECRET_STEPS[2], "ENCODE_PIPELINE_E2E_REDIS_URL", REDIS_SECRET_EXPR),
    ]
    precheck = next(
        step
        for step in job["steps"]
        if step.get("name") == "Verify Redis before scientific execution"
    )
    assert "ENCODE_PIPELINE_TEST_REDIS_URL" not in (precheck.get("env") or {})


def _expand_redis(value: str, marker: str) -> str:
    return str(value).replace(REDIS_SECRET_EXPR, marker).replace(REDIS_VAR_EXPR, marker)


def _mask_registered(text: str, secrets: list[str]) -> str:
    masked = text
    for secret in secrets:
        if secret:
            masked = masked.replace(secret, "***")
    return masked


def runner_env_headers(job: dict, marker: str) -> dict[str, str]:
    """Model the runner printing step env before the script starts.

    Secret values are registered before the first header. A configuration
    variable is not, so its first header still contains the value. ``add-mask``
    inside a script only affects headers printed after that script.
    """

    document = str(job)
    registered = (
        [marker]
        if marker and REDIS_SECRET_EXPR in document and REDIS_VAR_EXPR not in document
        else []
    )
    base = {
        key: _expand_redis(value, marker)
        for key, value in (job.get("env") or {}).items()
    }
    headers: dict[str, str] = {}
    for step in job["steps"]:
        shown = dict(base)
        for key, value in (step.get("env") or {}).items():
            shown[key] = _expand_redis(value, marker)
        script = _expand_redis(str(step.get("run") or ""), marker)
        lines = [_mask_registered(line, registered) for line in script.splitlines()]
        if shown:
            lines.append("env:")
            lines.extend(
                f"  {key}: {_mask_registered(value, registered)}"
                for key, value in shown.items()
            )
        name = step.get("name")
        if name:
            headers[name] = "\n".join(lines)
        if "::add-mask::" in str(step.get("run") or "") and marker not in registered:
            registered.append(marker)
    return headers


def _f9217e9_vars_shape(job: dict) -> dict:
    """The Redis wiring at f9217e9: four variable references, no secret."""

    rewritten = copy.deepcopy(job)
    for step in rewritten["steps"]:
        for key, value in (step.get("env") or {}).items():
            if value == REDIS_SECRET_EXPR:
                step["env"][key] = REDIS_VAR_EXPR
    return rewritten


def test_runner_coordinate_masks_redis_and_files_it_for_precheck(tmp_path: Path):
    marker = "unix:///helixweave-private-redis-marker/%2e/redis.sock"
    result = _run_coordinate_step(tmp_path, marker)
    visible = "\n".join(
        line
        for line in f"{result.stdout}\n{result.stderr}".splitlines()
        if not line.startswith("::")
    )
    stored_path = tmp_path / "runner-temp" / "helixweave-gate-redis-url"
    stored = stored_path.read_text(encoding="utf-8") if stored_path.is_file() else ""
    mask_ok = result.stdout.splitlines() == [
        "::add-mask::" + marker.replace("%", "%25")
    ]
    assert result.returncode == 0
    assert not visible.strip()
    assert "helixweave-private-redis-marker" not in visible
    assert stored == marker
    assert mask_ok
    assert stored_path.stat().st_mode & 0o777 == 0o600


def test_runner_coordinate_rejects_multiline_redis_without_echoing_it(tmp_path: Path):
    result = _run_coordinate_step(
        tmp_path,
        "unix:///helixweave-private-redis-marker/redis.sock\ninjected",
    )
    blob = f"{result.stdout}\n{result.stderr}"
    assert result.returncode != 0
    assert "helixweave-private-redis-marker" not in blob
    assert "injected" not in blob
    assert "must be a single line" in blob
    assert not (tmp_path / "runner-temp" / "helixweave-gate-redis-url").exists()


def test_runner_coordinate_missing_secret_fails_closed(tmp_path: Path):
    result = _run_coordinate_step(tmp_path, "")
    blob = f"{result.stdout}\n{result.stderr}"
    assert result.returncode != 0
    assert (
        "ENCODE_PIPELINE_TEST_REDIS_URL must be configured by the protected environment"
    ) in blob
    assert "redis://" not in blob
    assert "unix:" not in blob
    assert not (tmp_path / "runner-temp" / "helixweave-gate-redis-url").exists()


def test_secret_redis_is_masked_before_the_first_env_header():
    job = _load("ci.yml")["jobs"]["bulk-rnaseq-real-execution"]
    assert_bulk_redis_secret_contract(job)
    marker = "unix:///helixweave-private-redis-marker/redis.sock"
    headers = runner_env_headers(job, marker)
    for name in REDIS_SECRET_STEPS:
        header = headers[name]
        assert marker not in header
        assert "ENCODE_PIPELINE_TEST_REDIS_URL: ***" in header
    assert "ENCODE_PIPELINE_E2E_REDIS_URL: ***" in headers[REDIS_SECRET_STEPS[2]]
    precheck = headers["Verify Redis before scientific execution"]
    assert "ENCODE_PIPELINE_TEST_REDIS_URL" not in precheck
    document = HARNESS_DOC.read_text(encoding="utf-8")
    assert "environment secret\n`ENCODE_PIPELINE_TEST_REDIS_URL`" in document
    assert "Create\nthe secret before any protected run" in document
    assert "authorized operator" in document


def test_f9217e9_vars_shape_fails_the_secret_contract_and_leaks():
    marker = "unix:///helixweave-private-redis-marker/redis.sock"
    current = _load("ci.yml")["jobs"]["bulk-rnaseq-real-execution"]
    old = _f9217e9_vars_shape(current)
    with pytest.raises(AssertionError):
        assert_bulk_redis_secret_contract(old)
    headers = runner_env_headers(old, marker)
    assert marker in headers[REDIS_SECRET_STEPS[0]]
    assert marker not in headers[REDIS_SECRET_STEPS[1]]
    assert marker not in headers[REDIS_SECRET_STEPS[2]]
    assert str(old).count(REDIS_VAR_EXPR) == 4
    assert REDIS_SECRET_EXPR not in str(old)


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
        "fast-checks-shard",
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
    assert workflow["on"]["push"]["tags"] == ["v*"]
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


def test_deterministic_suite_is_two_disjoint_shards_with_one_canonical_consumer():
    """One logical suite: two exhaustive shards, one combined coverage artifact."""
    jobs = _load("ci.yml")["jobs"]
    shard_job = jobs["fast-checks-shard"]
    aggregate = jobs["fast-checks"]
    consumer = _runs(jobs["coverage"])

    assert shard_job["strategy"]["matrix"]["shard"] == [1, 2]
    assert shard_job["strategy"]["fail-fast"] is False
    assert shard_job["timeout-minutes"] == 45
    producer = _runs(shard_job)
    assert "Soft budget: 1200s (informational)" in producer
    assert (
        producer.count(
            "python3 -I -S scripts/checkout_bootstrap.py --repository-root . pytest"
        )
        == 1
    )
    assert "split_deterministic_shards.py" in producer
    assert "--root test --shard ${{ matrix.shard }} --of 2" in producer
    assert (
        "not platform_real_execution and not real_execution "
        "and not bulk_rnaseq_real_execution" in producer
    )
    assert "not full_main" not in producer
    shard_run = next(
        step["run"]
        for step in shard_job["steps"]
        if step.get("name") == "Run deterministic Python shard"
    )
    assert "EVENT_NAME" not in shard_run
    assert "--cov-fail-under=0" in producer
    assert "--cov-report=" in producer
    assert "--cov-report=xml" not in producer
    assert "--cov-report=json" not in producer
    assert "--cov-report=term" not in producer
    assert "--junitxml=pytest-report-shard-${{ matrix.shard }}.xml" in producer
    assert (
        "check_junit_outcomes.py" in producer
        and "pytest-report-shard-${{ matrix.shard }}.xml" in producer
    )
    shard_envs = [
        step.get("env", {})
        for step in shard_job["steps"]
        if step.get("name") == "Run deterministic Python shard"
    ]
    assert shard_envs[0]["COVERAGE_FILE"] == ".coverage.shard-${{ matrix.shard }}"

    assert aggregate["needs"] == "fast-checks-shard"
    assert "always()" in aggregate["if"]
    assert aggregate["timeout-minutes"] == 10
    aggregate_runs = _runs(aggregate)
    assert "checkout_bootstrap.py --repository-root . pytest" not in aggregate_runs
    assert "coverage combine" in aggregate_runs
    assert "python-coverage-shard-*-${{ github.run_id }}" in str(aggregate)
    for evidence in (
        ".coverage.shard-1",
        ".coverage.shard-2",
        "pytest-report-shard-1.xml",
        "pytest-report-shard-2.xml",
    ):
        assert evidence in aggregate_runs
    assert (
        "check_junit_outcomes.py" in aggregate_runs
        and "pytest-report-shard-1.xml pytest-report-shard-2.xml" in aggregate_runs
    )
    assert "pytest" not in consumer
    assert jobs["coverage"]["needs"] == "fast-checks"
    assert "always()" in jobs["coverage"]["if"]


def test_fast_checks_proves_python_provenance_before_importing_product():
    shard_steps = _load("ci.yml")["jobs"]["fast-checks-shard"]["steps"]
    shard_provenance_index = next(
        index
        for index, step in enumerate(shard_steps)
        if step.get("name") == "Verify current-checkout Python provenance"
    )
    shard_pytest_index = next(
        index
        for index, step in enumerate(shard_steps)
        if step.get("name") == "Run deterministic Python shard"
    )
    assert shard_provenance_index < shard_pytest_index
    for step_name in (
        "Verify current-checkout Python provenance",
        "Run deterministic Python shard",
    ):
        run = str(
            next(step for step in shard_steps if step.get("name") == step_name)["run"]
        )
        assert "python3 -I -S scripts/checkout_bootstrap.py" in run
        assert "--repository-root ." in run

    aggregate_steps = _load("ci.yml")["jobs"]["fast-checks"]["steps"]
    aggregate_provenance_index = next(
        index
        for index, step in enumerate(aggregate_steps)
        if step.get("name") == "Verify current-checkout Python provenance"
    )
    validation_index = next(
        index
        for index, step in enumerate(aggregate_steps)
        if step.get("name") == "Validate default config"
    )
    local_platform_index = next(
        index
        for index, step in enumerate(aggregate_steps)
        if step.get("name") == "Verify local platform CLI through clean bootstrap"
    )
    assert aggregate_provenance_index < validation_index
    assert aggregate_provenance_index < local_platform_index
    provenance = str(aggregate_steps[aggregate_provenance_index]["run"])
    validation = str(aggregate_steps[validation_index]["run"])
    local_platform = str(aggregate_steps[local_platform_index]["run"])
    for command in (provenance, validation, local_platform):
        assert "python3 -I -S scripts/checkout_bootstrap.py" in command
        assert "--repository-root ." in command
    assert " verify-checkout" in provenance
    assert " validate" in validation
    assert " local-platform --help" in local_platform


def test_protected_pytest_tiers_use_only_the_checkout_bootstrap():
    jobs = _load("ci.yml")["jobs"]
    protected_jobs = (
        "fast-checks-shard",
        "platform-real-execution",
        "real-execution",
        "bulk-rnaseq-real-execution",
    )
    canonical = "python3 -I -S scripts/checkout_bootstrap.py --repository-root . pytest"
    bare_pytest = (
        r"(?m)(?:^|[;&|]\s*|\n\s*)"
        r"(?:python(?:\d+(?:\.\d+)*)?\s+-m\s+pytest|pytest)\b"
    )

    for job_name in protected_jobs:
        suite_steps = [
            str(step.get("run", ""))
            for step in jobs[job_name]["steps"]
            if "--junitxml=" in str(step.get("run", ""))
        ]
        assert len(suite_steps) == 1, job_name
        suite = suite_steps[0]
        assert suite.count(canonical) == 1, job_name
        assert re.search(bare_pytest, suite) is None, job_name
        assert "-p no:cacheprovider" not in suite, job_name


def test_scientific_real_separates_the_test_runner_from_locked_toolchains():
    steps = _load("ci.yml")["jobs"]["real-execution"]["steps"]
    python_setup = next(
        step for step in steps if step.get("name") == "Setup Python test environment"
    )
    scientific_setup = next(
        step for step in steps if step.get("name") == "Setup scientific toolchain"
    )
    multiqc_setup = next(
        step for step in steps if step.get("name") == "Setup MultiQC toolchain"
    )
    install = next(
        step
        for step in steps
        if step.get("name") == "Install HelixWeave scientific test package"
    )
    suite = next(
        step for step in steps if step.get("name") == "Scientific real-execution suite"
    )

    assert python_setup["with"]["environment-file"] == "workflow/envs/ci-fast.lock"
    assert scientific_setup["with"]["environment-file"] == "workflow/envs/chipseq.lock"
    assert scientific_setup["with"]["environment-name"] == "chipseq"
    assert multiqc_setup["with"]["environment-file"] == "workflow/envs/multiqc.lock"
    assert multiqc_setup["with"]["environment-name"] == "chipseq-multiqc"
    for setup in (scientific_setup, multiqc_setup):
        assert setup["with"]["init-shell"] == "none"
        assert setup["with"]["generate-run-shell"] is False
    assert (
        steps.index(python_setup)
        < steps.index(install)
        < steps.index(scientific_setup)
        < steps.index(multiqc_setup)
        < steps.index(suite)
    )
    assert "micromamba run -n ci-fast" in suite["run"]
    for variable, environment, executable in (
        ("SAMTOOLS", "chipseq", "samtools"),
        ("FASTQC", "chipseq", "fastqc"),
        ("TRIM_GALORE", "chipseq", "trim_galore"),
        ("MULTIQC", "chipseq-multiqc", "multiqc"),
    ):
        assert (
            f'{variable}="$MAMBA_ROOT_PREFIX/envs/{environment}/bin/{executable}"'
            in suite["run"]
        )
    conda_setup = next(
        s for s in steps if s.get("name") == "Setup scientific conda CLI"
    )
    assert conda_setup["with"]["environment-file"] == "workflow/envs/ci-conda.lock"
    assert conda_setup["with"]["environment-name"] == "scientific-conda"
    assert steps.index(conda_setup) < steps.index(suite)
    assert 'SNAKEMAKE="$RUNNER_TEMP/scientific-harness-bin/snakemake"' in suite["run"]
    assert "--workflow-profile none --profile none" in suite["run"]
    assert "envs/scientific-conda/bin:" in suite["run"]
    assert "conda info --json" in suite["run"]
    assert "--maxfail=1" in suite["run"]
    assert _load("ci.yml")["jobs"]["real-execution"]["timeout-minutes"] == 45


def test_documented_python_timing_budgets_match_the_workflow():
    harness = (REPO_ROOT / "docs" / "development" / "harness.md").read_text(
        encoding="utf-8"
    )
    assert "| Deterministic Python shards | 20 min | 45 min |" in harness
    assert "| `fast-checks` aggregate | 1 min | 10 min |" in harness


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
    for artifact in (
        ".coverage",
        "coverage.xml",
        "coverage.json",
        "pytest-report-shard-1.xml",
        "pytest-report-shard-2.xml",
    ):
        assert artifact in settings["path"]

    steps = {step.get("name"): step for step in jobs["coverage"]["steps"]}
    ratchets = steps["Enforce complete-suite coverage ratchets"]
    changed_lines = steps["Enforce changed-line coverage"]
    assert ratchets["if"] == "needs.fast-checks.result == 'success'"
    assert "set -euo pipefail" in ratchets["run"]
    assert 'exit "$result"' in ratchets["run"]
    for floor in ("83", "88.45", "87.28", "89.06", "82.37"):
        assert f"--fail-under={floor} || result=1" in ratchets["run"]
    assert steps["Setup micromamba"]["id"] == "coverage-environment"
    assert "!cancelled()" in changed_lines["if"]
    assert "steps.coverage-environment.outcome == 'success'" in changed_lines["if"]
    assert "github.event_name == 'pull_request'" in changed_lines["if"]
    assert "diff-cover coverage.xml" in changed_lines["run"]
    assert "github.event.pull_request.base.sha" in changed_lines["run"]
    assert "--fail-under=80" in changed_lines["run"]

    config = (REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8")
    coverage_report = re.search(
        r"(?ms)^\[tool\.coverage\.report\]\n(.*?)(?=^\[|\Z)", config
    )
    assert coverage_report
    assert re.search(r"(?m)^fail_under\s*=\s*83\s*$", coverage_report.group(1))


def test_migration_snapshot_coverage_is_attributed_without_omitting_sources():
    jobs = _load("ci.yml")["jobs"]
    producer_runs = _runs(jobs["fast-checks-shard"])
    coverage_runs = _runs(jobs["coverage"])

    assert (
        'coverage_migration_root="$RUNNER_TEMP/'
        "helixweave-coverage-migrations-$GITHUB_RUN_ID-$GITHUB_RUN_ATTEMPT"
        '-shard-${{ matrix.shard }}"' in producer_runs
    )
    assert 'mkdir -m 700 -- "$coverage_migration_root"' in producer_runs
    assert (
        'export HELIXWEAVE_COVERAGE_MIGRATION_ROOT="$coverage_migration_root"'
        in producer_runs
    )
    assert "trap cleanup_coverage_migration_root EXIT" in producer_runs
    for source_root in (
        "src/encode_pipeline",
        "scripts",
        "workflow/lib",
        "containers",
    ):
        assert f"--cov={source_root}" in producer_runs
    assert '--cov="$coverage_migration_root"' in producer_runs

    config_path = REPO_ROOT / "pyproject.toml"
    config_text = config_path.read_text(encoding="utf-8")
    assert (
        '"${HELIXWEAVE_COVERAGE_MIGRATION_ROOT-'
        "/__helixweave_coverage_migration_root_not_configured__}/"
        'helixweave-migration-snapshot-*/alembic"' in config_text
    )
    configured_root = os.environ.get(
        "HELIXWEAVE_COVERAGE_MIGRATION_ROOT",
        "/__helixweave_coverage_migration_root_not_configured__",
    )
    config = Coverage(config_file=str(config_path))
    config.load()
    assert config.config.paths["migration_execution"] == [
        "src/encode_pipeline/persistence/alembic",
        f"{configured_root}/helixweave-migration-snapshot-*/alembic",
    ]
    assert config.config.run_omit == []
    assert config.config.report_omit is None
    assert "--omit" not in coverage_runs


def test_all_pytest_tiers_enforce_zero_skip_junit_outcomes():
    jobs = _load("ci.yml")["jobs"]
    shard_runs = _runs(jobs["fast-checks-shard"])
    shard_report = "pytest-report-shard-${{ matrix.shard }}.xml"
    assert f"--junitxml={shard_report}" in shard_runs
    assert "check_junit_outcomes.py" in shard_runs
    assert shard_report in shard_runs

    aggregate_runs = _runs(jobs["fast-checks"])
    assert (
        "check_junit_outcomes.py" in aggregate_runs
        and "pytest-report-shard-1.xml pytest-report-shard-2.xml" in aggregate_runs
    )

    for job_id, report in (
        ("platform-real-execution", "platform-real-report.xml"),
        ("real-execution", "scientific-real-report.xml"),
    ):
        runs = _runs(jobs[job_id])
        assert f"--junitxml={report}" in runs
        assert f"check_junit_outcomes.py {report}" in runs

    bulk_runs = _runs(jobs["bulk-rnaseq-real-execution"])
    assert "bulk-rnaseq-real-report.xml" in bulk_runs
    assert "check_junit_outcomes.py" in bulk_runs

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

    bulk_job = jobs["bulk-rnaseq-real-execution"]
    bulk_condition = bulk_job["if"]
    assert "workflow_dispatch" in bulk_condition
    assert "inputs.bulk_rnaseq_real_execution" in bulk_condition
    assert "pull_request" not in bulk_condition
    assert bulk_job["runs-on"] == [
        "self-hosted",
        "linux",
        "x64",
        "helixweave-bulk-rnaseq",
    ]
    bulk_runs = _runs(bulk_job)
    assert "stage_bulk_rnaseq_runtime_assets.py" in bulk_runs
    assert "--phase verify" in bulk_runs
    assert "-m bulk_rnaseq_real_execution test/bulk_rnaseq_real_execution" in bulk_runs
    assert "HELIXWEAVE_REQUIRE_BULK_RNASEQ_REAL_EXECUTION" in str(bulk_job)
    redis_precheck = next(
        step
        for step in bulk_job["steps"]
        if step.get("name") == "Verify Redis before scientific execution"
    )
    scientific = next(
        step
        for step in bulk_job["steps"]
        if step.get("name") == "Run rapid, STAR Salmon, and lifecycle acceptance"
    )
    names = [step.get("name") for step in bulk_job["steps"]]
    assert names.index("Install platform package") < names.index(redis_precheck["name"])
    assert names.index(
        "Verify pinned runtime admission without mutation or network"
    ) < names.index(redis_precheck["name"])
    assert names.index(redis_precheck["name"]) < names.index(scientific["name"])
    assert redis_precheck["shell"] == "micromamba-shell {0}"
    assert "continue-on-error" not in redis_precheck
    # The runner prints this map in the step log header. Redis stays off it;
    # the coordinate step files the secret for this precheck.
    assert set(redis_precheck["env"]) == set(scientific["env"]) - {
        "ENCODE_PIPELINE_TEST_REDIS_URL"
    }
    assert "ENCODE_PIPELINE_TEST_REDIS_URL" not in redis_precheck["env"]
    assert "${{ vars.ENCODE_PIPELINE_TEST_REDIS_URL }}" not in str(redis_precheck)
    precheck_run = redis_precheck["run"]
    assert "require_gate_settings()" in precheck_run
    assert "sys.tracebacklimit = 0" in precheck_run
    assert "raise SystemExit(1) from None" in precheck_run
    assert "echo" not in precheck_run
    assert "print(" not in precheck_run
    assert "ENCODE_PIPELINE_TEST_REDIS_URL" not in precheck_run
    assert "redis://" not in precheck_run
    assert "unix:" not in precheck_run
    coordinate = next(
        step
        for step in bulk_job["steps"]
        if step.get("name") == "Require runner-owned runtime coordinates"
    )
    assert names.index(coordinate["name"]) < names.index(redis_precheck["name"])
    coordinate_run = coordinate["run"]
    assert "printf '::add-mask::%s\\n' \"$redis_mask\"" in coordinate_run
    assert '>"$RUNNER_TEMP/helixweave-gate-redis-url"' in coordinate_run.replace(
        " ",
        "",
    )
    assert "GITHUB_ENV" not in coordinate_run
    assert "${{ vars.ENCODE_PIPELINE_TEST_REDIS_URL }}" not in coordinate_run


def test_bulk_product_gate_is_exact_head_path_free_and_coordinate_scoped():
    workflow = _load("ci.yml")
    dispatch_inputs = workflow["on"]["workflow_dispatch"]["inputs"]
    assert "bulk_rnaseq_expected_sha" in dispatch_inputs

    job = workflow["jobs"]["bulk-rnaseq-real-execution"]
    assert set(job.get("env", {})) == {"PYTHONDONTWRITEBYTECODE"}
    steps = job["steps"]
    checkout = next(
        step
        for step in steps
        if str(step.get("uses", "")).startswith("actions/checkout@")
    )
    assert checkout["uses"] == (
        "actions/checkout@fbc6f3992d24b796d5a048ff273f7fcc4a7b6c09"
    )
    assert checkout["with"]["persist-credentials"] is False
    assert checkout["with"]["token"] == "local-mirror-only"

    confirm = next(
        step
        for step in steps
        if step.get("name") == "Confirm exact checkout and clean source tree"
    )
    assert confirm["env"]["EXPECTED_REVIEWED_SHA"] == (
        "${{ inputs.bulk_rnaseq_expected_sha }}"
    )
    assert "^[0-9a-f]{40}$" in confirm["run"]
    assert 'test "$actual_sha" = "$EXPECTED_REVIEWED_SHA"' in confirm["run"]

    admitted_toolchain = next(
        step
        for step in steps
        if step.get("name") == "Require pre-admitted protected toolchain"
    )
    assert admitted_toolchain["env"] == {
        "EXPECTED_MICROMAMBA_VERSION": "2.8.1",
        "EXPECTED_MICROMAMBA_SHA256": (
            "9689782d863c05a1bf5d2d371ba527104e7a4eb4310c1637d8653b751aed9c82"
        ),
    }
    assert "micromamba --version" in admitted_toolchain["run"]
    assert "sha256sum" in admitted_toolchain["run"]

    micromamba = next(step for step in steps if step.get("name") == "Setup micromamba")
    assert micromamba["uses"] == (
        "mamba-org/setup-micromamba@d7c9bd84e824b79d2af72a2d4196c7f4300d3476"
    )
    assert micromamba["with"]["download-micromamba"] is False
    assert micromamba["with"]["cache-downloads"] is False
    assert micromamba["with"]["cache-environment"] is False
    assert micromamba["with"]["create-args"] == "--offline"
    assert "cache-downloads-key" not in micromamba["with"]
    assert "cache-environment-key" not in micromamba["with"]
    assert "micromamba-version" not in micromamba["with"]
    assert "micromamba-url" not in micromamba["with"]
    assert "micromamba-binary-path" not in micromamba["with"]

    node = next(
        step
        for step in steps
        if step.get("name") == "Setup Node.js for the protected product Gate"
    )
    admitted_node = next(
        step
        for step in steps
        if step.get("name") == "Require pre-admitted Node.js toolchain"
    )
    assert admitted_node["env"] == {
        "EXPECTED_NODE_VERSION": "v22.23.1",
        "EXPECTED_NPM_VERSION": "10.9.8",
    }
    assert "RUNNER_TOOL_CACHE" in admitted_node["run"]
    assert "x64.complete" in admitted_node["run"]
    assert 'bin/node" --version' in admitted_node["run"]
    assert 'bin/npm" --version' in admitted_node["run"]
    assert node["uses"] == (
        "actions/setup-node@49933ea5288caeca8642d1e84afbd3f7d6820020"
    )
    assert node["with"] == {
        "node-version": "22.23.1",
        "check-latest": False,
    }

    product = next(
        step
        for step in steps
        if step.get("name")
        == "Run public registry-to-results Bulk RNA-seq product Gate"
    )
    assert product["env"]["HELIXWEAVE_REQUIRE_BULK_RNASEQ_PRODUCT_GATE"] == "1"
    assert product["env"]["HELIXWEAVE_REQUIRE_BULK_RNASEQ_REAL_EXECUTION"] == "1"
    assert "bulk-rnaseq-unavailable.spec.ts" in product["run"]
    assert "continue-on-error" not in product

    protected_coordinates = {
        "HELIXWEAVE_BULK_RNASEQ_RUNTIME_ROOT",
        "HELIXWEAVE_BULK_RNASEQ_FIXTURE_MANIFEST",
        "ENCODE_PIPELINE_TEST_REDIS_URL",
        "ENCODE_PIPELINE_MANAGED_DOCKER_EXECUTABLE",
        "ENCODE_PIPELINE_MANAGED_DOCKER_SOCKET",
    }
    assert protected_coordinates <= set(product["env"])
    assert all(name not in job["env"] for name in protected_coordinates)

    dependencies = next(
        step
        for step in steps
        if step.get("name") == "Install protected product Gate dependencies"
    )
    assert dependencies["env"] == {
        "DO_NOT_TRACK": "1",
        "NPM_CONFIG_AUDIT": "false",
        "NPM_CONFIG_CACHE": "${{ runner.tool_cache }}/helixweave-npm-cache",
        "NPM_CONFIG_FUND": "false",
        "NPM_CONFIG_OFFLINE": "true",
        "PLAYWRIGHT_BROWSERS_PATH": (
            "${{ runner.tool_cache }}/helixweave-playwright/1.61.1"
        ),
        "SCARF_ANALYTICS": "false",
        "SCARF_NO_ANALYTICS": "true",
    }
    assert "unshare --user --map-root-user --net" in dependencies["run"]
    assert "npm ci --offline --no-audit --no-fund" in dependencies["run"]
    assert (
        "./node_modules/.bin/playwright install-deps chromium --dry-run"
        in (dependencies["run"])
    )
    assert "./node_modules/.bin/playwright install chromium" in dependencies["run"]
    assert "chromium.executablePath()" in dependencies["run"]
    assert "npx" not in dependencies["run"]
    assert "playwright install-deps chromium --dry-run" in dependencies["run"]
    assert "playwright install chromium" in dependencies["run"]
    assert dependencies["run"].index(
        "playwright install-deps chromium --dry-run"
    ) < dependencies["run"].index("playwright install chromium")
    assert "--with-deps" not in dependencies["run"]
    assert "continue-on-error" not in dependencies
    assert "sudo" not in str(job)
    assert product["env"]["PLAYWRIGHT_BROWSERS_PATH"] == (
        "${{ runner.tool_cache }}/helixweave-playwright/1.61.1"
    )

    audit = next(
        step
        for step in steps
        if step.get("name") == "Require path-free product evidence"
    )
    assert audit["if"] == "always()"
    assert "bulk-rnaseq-product-evidence.json" in audit["env"]["PRODUCT_EVIDENCE"]
    assert "set -euo pipefail" in audit["run"]
    assert "json.tool" in audit["run"]
    assert "continue-on-error" not in audit

    container_audit = next(
        step for step in steps if step.get("name") == "Require no managed containers"
    )
    assert container_audit["if"] == "always()"
    assert "set -euo pipefail" in container_audit["run"]
    assert "ps -aq" in container_audit["run"]
    assert "--filter" not in container_audit["run"]
    assert "continue-on-error" not in container_audit

    upload = next(
        step
        for step in steps
        if step.get("name") == "Upload path-free bulk RNA-seq evidence"
    )
    assert upload["if"] == "always()"
    assert upload["uses"] == (
        "actions/upload-artifact@ea165f8d65b6e75b540449e92b4886f43607fa02"
    )
    assert "**/evidence/**" in upload["with"]["path"]
    assert upload["with"]["if-no-files-found"] == "error"
    assert steps.index(audit) < steps.index(container_audit) < steps.index(upload)

    product_spec = (
        REPO_ROOT / "frontend" / "e2e" / "bulk-rnaseq-unavailable.spec.ts"
    ).read_text(encoding="utf-8")
    assert "test.skip" not in product_spec
    assert "test.fixme" not in product_spec


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


# The chipseq environment stays free of the HelixWeave distribution. Rule
# scripts in that environment import the checkout by inserting src/. The
# driver Snakefile rejects a foreign encode_pipeline before rules run.
HELIXWEAVE_DISTRIBUTION_NAMES = ("helixweave", "encode-pipeline", "encode_pipeline")


def _dependency_names(dependencies) -> set:
    names = set()
    for entry in dependencies:
        if isinstance(entry, str):
            names.add(re.split(r"[=<>\s]", entry, maxsplit=1)[0].lower())
        elif isinstance(entry, dict):
            for pip_entry in entry.get("pip", []):
                names.add(re.split(r"[=<>\[;\s]", pip_entry, maxsplit=1)[0].lower())
    return names


def _is_helixweave_distribution(name: str) -> bool:
    return any(
        name == dist or name.startswith(dist + "-") or name.startswith(dist + "_")
        for dist in HELIXWEAVE_DISTRIBUTION_NAMES
    )


def test_browser_e2e_installs_os_deps_without_the_azure_apt_mirror() -> None:
    job = _load("ci.yml")["jobs"]["browser-e2e"]
    step = next(
        item
        for item in job["steps"]
        if item.get("name") == "Install Playwright operating-system dependencies"
    )
    run = step["run"]
    assert "set -euo pipefail" in run
    assert "azure.archive.ubuntu.com" not in run
    assert "https://archive.ubuntu.com/ubuntu/" in run
    assert "https://security.ubuntu.com/ubuntu/" in run
    assert run.strip().endswith("npx playwright install-deps chromium")
    assert "--dry-run" not in run
    assert "continue-on-error" not in step
    bulk = _load("ci.yml")["jobs"]["bulk-rnaseq-real-execution"]
    assert bulk["runs-on"] == ["self-hosted", "linux", "x64", "helixweave-bulk-rnaseq"]
    assert "inputs.bulk_rnaseq_real_execution" in str(bulk["if"])
    assert "apt-mirrors.txt" not in _runs(bulk)


def test_chipseq_env_yml_excludes_helixweave() -> None:
    env = yaml.safe_load(
        (REPO_ROOT / "workflow" / "envs" / "chipseq.yml").read_text(encoding="utf-8")
    )
    names = _dependency_names(env.get("dependencies", []))
    offenders = [name for name in names if _is_helixweave_distribution(name)]
    assert not offenders, f"chipseq.yml must not install HelixWeave: {offenders}"


def test_chipseq_env_lock_excludes_helixweave() -> None:
    lock_path = REPO_ROOT / "workflow" / "envs" / "chipseq.lock"
    offenders = []
    for raw in lock_path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#") or line.startswith("@"):
            continue
        token = line.split(" @ ", 1)[0] if " @ " in line else line.rsplit("/", 1)[-1]
        token = token.split("#", 1)[0].lower()  # strip conda-md5 fragment
        if _is_helixweave_distribution(token):
            offenders.append(raw)
    assert not offenders, f"chipseq.lock must not pin HelixWeave: {offenders}"
