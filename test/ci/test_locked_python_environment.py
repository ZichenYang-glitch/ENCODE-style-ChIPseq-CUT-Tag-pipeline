"""Contracts for the reproducible Python CI environment."""

from pathlib import Path
import re

from coverage import Coverage
import yaml


REPO_ROOT = Path(__file__).resolve().parents[2]
ENVIRONMENT_YAML = REPO_ROOT / "workflow" / "envs" / "ci-fast.yml"
ENVIRONMENT_LOCK = ENVIRONMENT_YAML.with_suffix(".lock")
LINT_ENVIRONMENT_YAML = REPO_ROOT / "workflow" / "envs" / "ci-lint.yml"
LINT_ENVIRONMENT_LOCK = LINT_ENVIRONMENT_YAML.with_suffix(".lock")
WORKFLOW_DIR = REPO_ROOT / ".github" / "workflows"

REQUIRED_DIRECT_DEPENDENCIES = {
    "alembic",
    "coverage",
    "diff-cover",
    "fakeredis",
    "fastapi",
    "httpx",
    "jsonschema",
    "lupa",
    "mypy",
    "packaging",
    "pip",
    "pytest",
    "pytest-cov",
    "python-build",
    "pyyaml",
    "redis-py",
    "rq",
    "ruff",
    "setuptools",
    "snakemake-minimal",
    "sqlalchemy",
    "uvicorn",
    "wheel",
}
REQUIRED_LINT_DEPENDENCIES = {
    "packaging",
    "pyyaml",
    "ruff",
    "snakefmt",
    "snakemake-minimal",
}


def _dependency_name(specification: str) -> str:
    return re.split(r"[ <>=!]", specification, maxsplit=1)[0]


def test_ci_environment_declares_build_runtime_test_and_gate_dependencies():
    environment = yaml.safe_load(ENVIRONMENT_YAML.read_text(encoding="utf-8"))
    specifications = environment["dependencies"]
    names = {_dependency_name(specification) for specification in specifications}

    assert REQUIRED_DIRECT_DEPENDENCIES <= names
    assert "coverage >=7.10.6,<8" in specifications
    assert "diff-cover >=9,<10" in specifications
    assert "pytest-cov >=7,<8" in specifications
    assert "python-build >=1.2,<2" in specifications
    assert "redis-py >=7,<9" in specifications
    assert "rq >=2.10,<2.11" in specifications
    assert "ruff ==0.15.21" in specifications
    assert "packaging <26" in specifications


def test_lint_environment_locks_formatters_without_an_editable_install():
    environment = yaml.safe_load(LINT_ENVIRONMENT_YAML.read_text(encoding="utf-8"))
    specifications = environment["dependencies"]
    names = {_dependency_name(specification) for specification in specifications}

    assert REQUIRED_LINT_DEPENDENCIES <= names
    assert "ruff ==0.15.21" in specifications
    assert "snakefmt ==2.0.3" in specifications
    lint_workflow = (WORKFLOW_DIR / "lint.yml").read_text(encoding="utf-8")
    assert "workflow/envs/ci-lint.lock" in lint_workflow
    assert "pip install" not in lint_workflow


def test_explicit_lock_contains_every_declared_gate_and_runtime_package():
    lock = ENVIRONMENT_LOCK.read_text(encoding="utf-8")

    assert re.search(r"^# input_hash: [0-9a-f]{64}$", lock, re.MULTILINE)
    assert "\n@EXPLICIT\n" in lock
    for package in REQUIRED_DIRECT_DEPENDENCIES:
        assert re.search(rf"/{re.escape(package)}-[0-9]", lock), package

    lint_lock = LINT_ENVIRONMENT_LOCK.read_text(encoding="utf-8")
    assert re.search(r"^# input_hash: [0-9a-f]{64}$", lint_lock, re.MULTILINE)
    assert "\n@EXPLICIT\n" in lint_lock
    for package in REQUIRED_LINT_DEPENDENCIES:
        assert re.search(rf"/{re.escape(package)}-[0-9]", lint_lock), package


def test_editable_ci_installs_are_offline_and_do_not_resolve_dependencies():
    workflows = {
        path.name: yaml.safe_load(path.read_text(encoding="utf-8"))
        for path in WORKFLOW_DIR.glob("*.yml")
    }
    install_steps = []
    for workflow in workflows.values():
        for job in workflow["jobs"].values():
            for step in job["steps"]:
                command = str(step.get("run", ""))
                if "pip install" in command:
                    install_steps.append(command)

    assert install_steps
    for command in install_steps:
        assert "python3 -m pip install" in command
        assert "--no-index" in command
        assert "--no-deps" in command
        assert "--no-build-isolation" in command
        assert "python3 -m pip check" in command


def test_coverage_source_includes_all_authored_production_python_roots():
    config = Coverage(config_file=str(REPO_ROOT / "pyproject.toml"))
    config.load()

    assert set(config.config.source) == {
        "src/encode_pipeline",
        "scripts",
        "workflow/lib",
        "containers",
    }
