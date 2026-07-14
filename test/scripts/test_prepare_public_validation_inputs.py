"""Offline CLI contracts for the public validation input inventory."""

import csv
import io
import json
import os
from pathlib import Path
import subprocess
import sys


REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "scripts" / "prepare_public_validation_inputs.py"
QUEUES = [
    "tf_chip_cebpb",
    "broad_histone_h3k27me3",
    "atac_keratinocyte",
    "cuttag_h3k27me3",
]
ACCESSIONS = [
    "ENCSR000DYI",
    "ENCSR000AKB",
    "ENCSR254KDA",
    "GSE145187",
]


def _run_inventory(tmp_path: Path, *arguments: str) -> str:
    environment = os.environ.copy()
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    result = subprocess.run(
        [sys.executable, str(SCRIPT), *arguments],
        cwd=tmp_path,
        env=environment,
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stderr
    assert result.stderr == ""
    assert list(tmp_path.iterdir()) == []
    return result.stdout


def test_default_tsv_inventory_is_deterministic_and_ordered(tmp_path):
    first = _run_inventory(tmp_path)
    second = _run_inventory(tmp_path)

    assert first == second
    rows = list(csv.DictReader(io.StringIO(first), delimiter="\t"))
    assert list(rows[0]) == [
        "queue",
        "assay",
        "target",
        "accession",
        "source_url",
        "genome",
        "layout",
        "peak_mode",
        "n_treatment_bioreps",
        "n_control_samples",
        "stage5_idr",
        "seacr_enabled",
        "notes",
    ]
    assert [row["queue"] for row in rows] == QUEUES
    assert [row["accession"] for row in rows] == ACCESSIONS
    assert [row["stage5_idr"] for row in rows] == [
        "true",
        "false",
        "false",
        "false",
    ]


def test_json_inventory_preserves_typed_dataset_contract(tmp_path):
    datasets = json.loads(_run_inventory(tmp_path, "--json"))

    assert [dataset["queue"] for dataset in datasets] == QUEUES
    assert [dataset["accession"] for dataset in datasets] == ACCESSIONS
    assert datasets[0]["has_control"] is True
    assert datasets[0]["stage5_idr"] is True
    assert datasets[2]["has_control"] is False
    assert datasets[3]["seacr_enabled"] is True


def test_dry_run_describes_every_dataset_without_downloading(tmp_path):
    output = _run_inventory(tmp_path, "--dry-run")

    assert "[dry-run] 4 dataset(s) defined" in output
    assert [
        line.strip().split(":", maxsplit=1)[0]
        for line in output.splitlines()
        if line.startswith("  ")
    ] == QUEUES
    assert output.rstrip().endswith("[dry-run] No downloads performed")


def test_metadata_checklist_is_offline_and_covers_every_accession(tmp_path):
    output = _run_inventory(tmp_path, "--check-metadata")

    assert "Public Validation Dataset Metadata Checklist" in output
    assert "=== Acceptance Criteria ===" in output
    assert "=== Manual Verification Steps ===" in output
    for queue, accession in zip(QUEUES, ACCESSIONS, strict=True):
        assert f"--- {queue} ({accession}) ---" in output
    assert output.rstrip().endswith("No network requests were made by this script.")


def test_report_stub_listing_preserves_the_public_cli_contract(tmp_path):
    output = _run_inventory(tmp_path, "--report-stubs")

    for queue, accession in zip(QUEUES, ACCESSIONS, strict=True):
        assert queue in output
        assert accession in output
        assert f"docs/release-checks/public-data-runs/{queue}.md" in output
    assert (
        "Template: docs/release-checks/public-data-execution-report-template.md"
        in output
    )
    assert output.rstrip().endswith("No downloads performed.")
