"""Shared pytest fixtures for the test suite."""

import os

import pytest


def pytest_addoption(parser):
    """Add custom command-line options."""
    parser.addoption(
        "--update-snapshots",
        action="store_true",
        default=False,
        help="Update DAG snapshot fixtures from current dry-run output",
    )


@pytest.fixture(scope="session")
def repo_root():
    """Return the repository root directory."""
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


@pytest.fixture(scope="session")
def snakefile(repo_root):
    """Return the path to the workflow Snakefile."""
    return os.path.join(repo_root, "workflow", "Snakefile")


@pytest.fixture(scope="session")
def profiles_dir(repo_root):
    """Return the directory containing test profiles."""
    return os.path.join(repo_root, "test", "profiles")


@pytest.fixture(scope="session")
def snapshots_dir(repo_root):
    """Return the directory for DAG snapshot fixtures."""
    return os.path.join(repo_root, "test", "fixtures", "dag_snapshots")


@pytest.fixture(scope="session")
def smoke_profiles():
    """Return the list of smoke-test profile names."""
    return [
        "chipseq_se_noctrl",
        "chipseq_pe_noctrl",
        "chipseq_pe_ctrlsample",
        "cuttag_pe_noctrl",
        "cuttag_pe_seacr",
        "chipseq_idr_dryrun",
        "chipseq_pe_external_ctrlbam",
        "mnase_pe_noctrl",
    ]
