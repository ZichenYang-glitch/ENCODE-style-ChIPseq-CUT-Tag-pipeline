"""Tests for the sentinel-owned results-visibility demo project."""

from __future__ import annotations

import csv
import io
from pathlib import Path
import subprocess
import sys

import pytest

from encode_pipeline.adapters.encode_qc import _QC_HEADER
from scripts.results_visibility_fixture import (
    FIXTURE_SENTINEL_CONTENT,
    FIXTURE_SENTINEL_NAME,
    prepare_results_visibility_fixture,
)


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]


def test_fixture_contains_deterministic_project_and_four_inputs(tmp_path):
    project_root = tmp_path / "runtime" / "results-visibility-project"

    inputs = prepare_results_visibility_fixture(project_root)

    configs = (
        inputs.results_config,
        inputs.cancel_config,
        inputs.empty_config,
        inputs.malformed_config,
    )
    assert [config["threads"] for config in configs] == [1, 2, 3, 4]
    assert {config["outdir"] for config in configs} == {"results"}
    assert all("samples" not in config for config in configs)
    assert all(
        config["replicate_analysis"] == {"enabled": False}
        and config["chipseq_idr"] == {"enabled": False}
        for config in configs
    )
    assert all("stage4b" not in config and "stage5" not in config for config in configs)
    assert inputs.samples_path.is_absolute()
    assert inputs.samples_path == project_root / "samples.tsv"
    with inputs.samples_path.open(encoding="utf-8", newline="") as handle:
        samples = list(csv.DictReader(handle, delimiter="\t"))
    assert len(samples) == 1
    assert samples[0]["sample"] == "C1"
    assert samples[0]["role"] == "treatment"
    assert samples[0]["target"] == "CTCF"

    required = (
        "pyproject.toml",
        "src/encode_pipeline/__init__.py",
        "docs/architecture/artifact-inventory.yaml",
        "profiles/default/config.yaml",
        "workflow/Snakefile",
        "scripts/results_visibility_task.py",
        "samples.tsv",
    )
    assert all((project_root / relative).is_file() for relative in required)
    assert (project_root / FIXTURE_SENTINEL_NAME).read_text(encoding="utf-8") == (
        FIXTURE_SENTINEL_CONTENT
    )

    reader = csv.DictReader(
        io.StringIO(inputs.expected_qc_summary, newline=""),
        delimiter="\t",
    )
    rows = list(reader)
    assert tuple(reader.fieldnames or ()) == _QC_HEADER
    assert len(rows) == 1
    assert rows[0] == {
        **{column: "NA" for column in _QC_HEADER},
        "sample": "C1",
        "assay": "chipseq",
        "total_reads": "1000",
        "frip": "0.125",
        "peak_count": "50",
        "percent_duplication": "0.2",
        "estimated_library_size": "9007199254740993",
        "nrf": "0.8",
        "pbc1": "0.75",
        "pbc2": "3.0",
    }


def test_fixture_refuses_unowned_existing_directory_without_deleting_it(tmp_path):
    project_root = tmp_path / "runtime" / "results-visibility-project"
    project_root.mkdir(parents=True)
    unrelated = project_root / "keep.txt"
    unrelated.write_text("keep\n", encoding="utf-8")

    with pytest.raises(ValueError, match="not owned"):
        prepare_results_visibility_fixture(project_root)

    assert unrelated.read_text(encoding="utf-8") == "keep\n"


def test_fixture_replaces_only_owned_project_and_preserves_runtime_state(tmp_path):
    runtime_root = tmp_path / "runtime"
    project_root = runtime_root / "results-visibility-project"
    first = prepare_results_visibility_fixture(project_root)
    old_project_file = project_root / "old.txt"
    old_project_file.write_text("replace me\n", encoding="utf-8")
    database = runtime_root / "platform.db"
    database.write_bytes(b"canonical-sqlite-marker")

    second = prepare_results_visibility_fixture(project_root)

    assert not old_project_file.exists()
    assert database.read_bytes() == b"canonical-sqlite-marker"
    assert second.expected_qc_summary == first.expected_qc_summary
    assert second.results_config == first.results_config
    assert "samples" not in second.results_config


@pytest.mark.parametrize("unsafe", [Path("/"), Path("/tmp"), REPOSITORY_ROOT])
def test_fixture_rejects_dangerous_project_roots(unsafe):
    with pytest.raises(ValueError, match="dedicated"):
        prepare_results_visibility_fixture(unsafe)


@pytest.mark.parametrize(
    ("mode", "qc_expected", "valid_qc"),
    [
        ("results", True, True),
        ("empty", False, False),
        ("malformed", True, False),
    ],
)
def test_fixture_task_writes_truthful_mode_outputs(
    tmp_path,
    mode,
    qc_expected,
    valid_qc,
):
    project_root = tmp_path / "runtime" / "results-visibility-project"
    inputs = prepare_results_visibility_fixture(project_root)
    workspace = tmp_path / "workspace" / mode
    workspace.mkdir(parents=True)
    complete = Path("result/complete.txt")
    manifest = Path("results/multiqc/result_manifest.tsv")
    qc_summary = Path("results/C1/01_qc/C1.qc_summary.tsv")

    completed = subprocess.run(
        [
            sys.executable,
            str(project_root / "scripts/results_visibility_task.py"),
            mode,
            str(complete),
            str(manifest),
            str(qc_summary),
        ],
        cwd=workspace,
        capture_output=True,
        text=True,
        check=False,
        timeout=10,
    )

    assert completed.returncode == 0, completed.stderr
    assert (workspace / complete).read_text(encoding="utf-8") == "success\n"
    manifest_text = (workspace / manifest).read_text(encoding="utf-8")
    assert "result_manifest\tpresent" in manifest_text
    assert (workspace / qc_summary).exists() is qc_expected
    assert ("qc_summary\tpresent" in manifest_text) is qc_expected
    if valid_qc:
        assert (workspace / qc_summary).read_text(encoding="utf-8") == (
            inputs.expected_qc_summary
        )
    elif qc_expected:
        assert (workspace / qc_summary).read_text(encoding="utf-8") != (
            inputs.expected_qc_summary
        )
