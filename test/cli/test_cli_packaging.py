"""Tests for CLI package boundaries and repo-root handling."""

import ast
import importlib
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import pytest

from encode_pipeline.cli import dag as dag_cli

REPO_ROOT = Path(__file__).resolve().parents[2]


# ---------------------------------------------------------------------------
# Package boundary: no CLI module should import scripts.* or use sys.path.insert
# ---------------------------------------------------------------------------

CLI_MODULES = [
    "encode_pipeline.cli.validate",
    "encode_pipeline.cli.manifest",
    "encode_pipeline.cli.dag",
    "encode_pipeline.cli._logging",
]


def _read_module_source(module_name):
    mod = importlib.import_module(module_name)
    return Path(mod.__file__).read_text(encoding="utf-8")


@pytest.mark.parametrize("module_name", CLI_MODULES)
def test_cli_module_does_not_import_scripts(module_name):
    source = _read_module_source(module_name)
    tree = ast.parse(source)
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                assert not alias.name.startswith("scripts"), (
                    f"{module_name} imports {alias.name}"
                )
        if isinstance(node, ast.ImportFrom):
            assert node.module is None or not node.module.startswith("scripts"), (
                f"{module_name} imports from {node.module}"
            )


def test_manifest_cli_has_no_sys_path_insert():
    source = _read_module_source("encode_pipeline.cli.manifest")
    assert "sys.path.insert" not in source


# ---------------------------------------------------------------------------
# encode-dag repo-root resolution
# ---------------------------------------------------------------------------


def _dag_cli_environment(**updates):
    snakemake = dag_cli._find_snakemake()
    assert snakemake is not None
    environment = {
        **os.environ,
        "PYTHONDONTWRITEBYTECODE": "1",
        "SNAKEMAKE": snakemake,
    }
    environment.update(updates)
    return environment


def test_dag_snakemake_resolution_has_explicit_precedence_and_fallback(
    tmp_path, monkeypatch
):
    monkeypatch.setenv("SNAKEMAKE", "/configured/snakemake")
    monkeypatch.setattr(
        dag_cli.shutil,
        "which",
        lambda _executable: pytest.fail("PATH lookup must follow SNAKEMAKE"),
    )
    assert dag_cli._find_snakemake() == "/configured/snakemake"

    monkeypatch.delenv("SNAKEMAKE")
    monkeypatch.setattr(dag_cli.shutil, "which", lambda _executable: "/path/snakemake")
    assert dag_cli._find_snakemake() == "/path/snakemake"

    home = tmp_path / "home"
    candidate = home / "miniconda3" / "envs" / "chipseq" / "bin" / "snakemake"
    candidate.parent.mkdir(parents=True)
    candidate.touch()
    candidate.chmod(0o755)
    monkeypatch.setattr(dag_cli.shutil, "which", lambda _executable: None)
    monkeypatch.setattr(dag_cli.Path, "home", classmethod(lambda _cls: home))

    assert dag_cli._find_snakemake() == str(candidate)

    candidate.unlink()
    assert dag_cli._find_snakemake() is None


def test_dag_cli_resolves_explicit_repo_root():
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "encode_pipeline.cli.dag",
            "diff",
            "--profile",
            "chipseq_se_noctrl",
            "--repo-root",
            str(REPO_ROOT),
        ],
        cwd=tempfile.gettempdir(),
        capture_output=True,
        text=True,
        env=_dag_cli_environment(),
    )
    assert result.returncode == 0, result.stderr
    assert "No differences for chipseq_se_noctrl" in result.stdout


def test_dag_cli_fails_without_repo_root_outside_repo():
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "encode_pipeline.cli.dag",
            "diff",
            "--profile",
            "chipseq_se_noctrl",
        ],
        cwd=tempfile.gettempdir(),
        capture_output=True,
        text=True,
        env=_dag_cli_environment(),
    )
    assert result.returncode == 2
    assert "Could not locate encode-pipeline repository root" in result.stderr


def test_dag_cli_resolves_env_repo_root():
    env = _dag_cli_environment(ENCODE_PIPELINE_REPO_ROOT=str(REPO_ROOT))
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "encode_pipeline.cli.dag",
            "diff",
            "--profile",
            "chipseq_se_noctrl",
        ],
        cwd=tempfile.gettempdir(),
        capture_output=True,
        text=True,
        env=env,
    )
    assert result.returncode == 0, result.stderr
    assert "No differences for chipseq_se_noctrl" in result.stdout


# ---------------------------------------------------------------------------
# encode-validate structured logging
# ---------------------------------------------------------------------------


def test_validate_cli_emits_single_structured_log_line():
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "encode_pipeline.cli.validate",
            "--config",
            str(REPO_ROOT / "config" / "config.yaml"),
        ],
        capture_output=True,
        text=True,
        env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
    )
    assert result.returncode == 0, result.stderr
    log_lines = [
        line
        for line in result.stderr.splitlines()
        if line.startswith("[") and "encode-pipeline" in line
    ]
    assert len(log_lines) == 1, f"Expected 1 log line, got: {log_lines}"
    assert "config_hash=" in log_lines[0]
    assert "config_hash=na" not in log_lines[0]


# ---------------------------------------------------------------------------
# encode-manifest package boundary and output parity
# ---------------------------------------------------------------------------


def test_manifest_cli_matches_legacy_script(tmp_path):
    legacy = tmp_path / "legacy.tsv"
    new_cli = tmp_path / "new.tsv"

    subprocess.run(
        [
            sys.executable,
            str(REPO_ROOT / "scripts" / "make_manifest.py"),
            "--config",
            str(REPO_ROOT / "config" / "config.yaml"),
            "--output",
            str(legacy),
        ],
        check=True,
        capture_output=True,
        text=True,
        env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
    )
    subprocess.run(
        [
            sys.executable,
            "-m",
            "encode_pipeline.cli.manifest",
            "--config",
            str(REPO_ROOT / "config" / "config.yaml"),
            "--output",
            str(new_cli),
        ],
        check=True,
        capture_output=True,
        text=True,
        env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"},
    )

    assert legacy.read_text() == new_cli.read_text()
