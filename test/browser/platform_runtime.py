#!/usr/bin/env python3
"""Prepare a controlled tiny workflow and exec the real local platform stack."""

from __future__ import annotations

import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import sys
from uuid import uuid4

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from scripts.results_visibility_fixture import (  # noqa: E402
    ResultsVisibilityInputs,
    prepare_results_visibility_fixture,
)
from bulk_product_runtime import (  # noqa: E402
    prepare_bulk_product_browser_runtime,
)


OWNERSHIP_SENTINEL = ".encode-platform-playwright-owned"
PRODUCT_BULK_GATE_ENV = "HELIXWEAVE_REQUIRE_BULK_RNASEQ_PRODUCT_GATE"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare_bulk_authoring_fixture(runtime_root: Path) -> dict[str, object]:
    """Create tiny synthetic authoring inputs that are never used for execution."""
    fixture_root = runtime_root / "bulk-authoring"
    reads_root = fixture_root / "reads"
    reference_root = fixture_root / "reference"
    rrna_root = fixture_root / "rrna"
    sortmerna_index_root = rrna_root / "sortmerna-index"
    for directory in (reads_root, reference_root, sortmerna_index_root):
        directory.mkdir(parents=True, exist_ok=True)

    read_names = (
        "pairedA_1.lib1.L001.fastq.gz",
        "pairedA_2.lib1.L001.fastq.gz",
        "pairedA_1.lib1.L002.fastq.gz",
        "pairedA_2.lib1.L002.fastq.gz",
        "singleB.lib1.L001.fastq.gz",
    )
    synthetic_fastq = b"@synthetic-read\nACGTACGT\n+\nIIIIIIII\n"
    for read_name in read_names:
        (reads_root / read_name).write_bytes(gzip.compress(synthetic_fastq, mtime=0))

    fasta_path = reference_root / "synthetic.fa"
    fasta_path.write_text(">chrSynthetic\nACGTACGTACGT\n", encoding="utf-8")
    gtf_path = reference_root / "synthetic.gtf"
    gtf_path.write_text(
        'chrSynthetic\tHelixWeave\texon\t1\t12\t.\t+\t.\tgene_id "GENE1"; '
        'transcript_id "TX1";\n',
        encoding="utf-8",
    )
    rrna_fasta_path = rrna_root / "synthetic-rrna.fa"
    rrna_fasta_path.write_text(">rrnaSynthetic\nACGTACGT\n", encoding="utf-8")
    rrna_manifest_path = rrna_root / "database-manifest.json"
    rrna_manifest_path.write_text(
        json.dumps(
            {
                "database": rrna_fasta_path.name,
                "sha256": _sha256(rrna_fasta_path),
            },
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    sortmerna_manifest_path = sortmerna_index_root / "identity.json"
    sortmerna_manifest_path.write_text(
        json.dumps(
            {
                "database_manifest_sha256": _sha256(rrna_manifest_path),
                "fixture": "authoring-only",
            },
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    samples_path = fixture_root / "samples.tsv"
    header = "sample\tlibrary\tlane\tlayout\tfastq_1\tfastq_2\tstrandedness\tplatform"
    rows = (
        (
            "pairedA",
            "lib1",
            "L001",
            "PE",
            reads_root / read_names[0],
            reads_root / read_names[1],
            "reverse",
            "ILLUMINA",
        ),
        (
            "pairedA",
            "lib1",
            "L002",
            "PE",
            reads_root / read_names[2],
            reads_root / read_names[3],
            "reverse",
            "ILLUMINA",
        ),
        (
            "singleB",
            "lib1",
            "L001",
            "SE",
            reads_root / read_names[4],
            "",
            "auto",
            "ILLUMINA",
        ),
    )
    samples_path.write_text(
        "\n".join([header, *("\t".join(map(str, row)) for row in rows)]) + "\n",
        encoding="utf-8",
    )

    config = {
        "standard": {
            "analysis": {"alignment": "star", "quantification": "salmon"},
            "reference": {
                "reference_id": "HelixWeave-synthetic",
                "fasta": str(fasta_path),
                "fasta_sha256": _sha256(fasta_path),
                "gtf": str(gtf_path),
                "gtf_sha256": _sha256(gtf_path),
                "annotation_style": "gencode",
            },
            "trimming": {"enabled": True, "tool": "fastp"},
            "ribosomal_rna_removal": {
                "enabled": True,
                "tool": "sortmerna",
                "save_filtered_reads": True,
                "database_manifest": {
                    "path": str(rrna_manifest_path),
                    "identity_sha256": _sha256(rrna_manifest_path),
                },
                "sortmerna_index": {
                    "path": str(sortmerna_index_root),
                    "identity_sha256": _sha256(sortmerna_manifest_path),
                },
            },
            "umi": {
                "enabled": True,
                "mode": "read_name",
                "deduplication_tool": "umitools",
                "read_name_separator": ":",
                "grouping_method": "directional",
                "emit_dedup_stats": False,
                "primary_alignments_only": False,
            },
            "outputs": {"trimmed_reads": True, "umi_intermediates": True},
        },
        "advanced": {"min_trimmed_reads": 1},
    }
    return {
        "bulkWorkflowId": "bulk-rnaseq",
        "bulkSamplesPath": str(samples_path),
        "bulkConfig": config,
        "bulkOptions": {},
        "bulkExpectedExecution": "not_configured",
        "bulkRequiredArtifactOutputTypes": [],
        "bulkRequiredQcMetricKeys": [],
    }


def prepare_bulk_browser_fixture(
    runtime_root: Path,
    environ: dict[str, str] | None = None,
) -> tuple[dict[str, object], dict[str, str]]:
    """Select ordinary authoring or the explicitly protected product fixture."""
    source = dict(os.environ if environ is None else environ)
    requested = source.get(PRODUCT_BULK_GATE_ENV)
    if requested is None:
        return prepare_bulk_authoring_fixture(runtime_root), {}
    if requested != "1":
        raise ValueError(f"{PRODUCT_BULK_GATE_ENV} must be exactly 1 when configured")
    prepared = prepare_bulk_product_browser_runtime(runtime_root, source)
    return dict(prepared.manifest_fields), dict(prepared.deployment_environment)


def write_manifest(
    runtime_root: Path,
    queue_name: str,
    inputs: ResultsVisibilityInputs,
    bulk_authoring: dict[str, object],
) -> None:
    """Write the bounded browser-visible fixture contract."""
    manifest = {
        "workflowId": "encode-style-chipseq-cuttag-atac-mnase",
        "samplesPath": str(inputs.samples_path),
        "resultsConfig": inputs.results_config,
        "cancelConfig": inputs.cancel_config,
        "emptyConfig": inputs.empty_config,
        "malformedConfig": inputs.malformed_config,
        "expectedQcSummary": inputs.expected_qc_summary,
        "runtimeRoot": str(runtime_root),
        "workspaceRoot": str(runtime_root / "workspaces"),
        "markerRoot": str(runtime_root / "tmp"),
        "queueName": queue_name,
    }
    manifest.update(bulk_authoring)
    (runtime_root / "runtime.json").write_text(
        json.dumps(manifest, sort_keys=True), encoding="utf-8"
    )


def prepare_owned_runtime_root(runtime_value: str, runtime_owner: str) -> Path:
    """Reset only the invocation-owned temporary runtime directory."""
    runtime_root = Path(runtime_value).resolve()
    if runtime_root in {Path("/"), Path("/tmp"), Path.home(), REPOSITORY_ROOT}:
        raise ValueError("E2E runtime root must be a dedicated temporary directory")
    sentinel = runtime_root / OWNERSHIP_SENTINEL
    if not sentinel.is_file() or sentinel.read_text(encoding="utf-8") != runtime_owner:
        raise ValueError("E2E runtime root is not owned by this Playwright invocation")
    shutil.rmtree(runtime_root)
    runtime_root.mkdir(parents=True)
    (runtime_root / OWNERSHIP_SENTINEL).write_text(runtime_owner, encoding="utf-8")
    return runtime_root


def main() -> None:
    runtime_value = os.environ.get("ENCODE_PIPELINE_E2E_ROOT")
    runtime_owner = os.environ.get("ENCODE_PIPELINE_E2E_OWNER")
    if not runtime_value or not runtime_owner:
        raise ValueError("Playwright runtime root and owner token are required")
    runtime_root = prepare_owned_runtime_root(runtime_value, runtime_owner)
    project_root = runtime_root / "project"
    inputs = prepare_results_visibility_fixture(project_root)
    bulk_authoring, deployment_environment = prepare_bulk_browser_fixture(runtime_root)
    queue_name = f"encode-pipeline-browser-{uuid4().hex}"
    write_manifest(runtime_root, queue_name, inputs, bulk_authoring)
    redis_url = os.environ.get(
        "ENCODE_PIPELINE_E2E_REDIS_URL", "redis://127.0.0.1:6380/0"
    )
    argv = [
        sys.executable,
        str(REPOSITORY_ROOT / "scripts" / "run_local_platform.py"),
        "--project-root",
        str(project_root),
        "--frontend-root",
        str(REPOSITORY_ROOT / "frontend"),
        "--runtime-root",
        str(runtime_root),
        "--redis-url",
        redis_url,
        "--queue-name",
        queue_name,
        "--api-port",
        "8010",
        "--frontend-port",
        "4173",
        "--readiness-timeout",
        "300" if deployment_environment else "120",
    ]
    child_environment = os.environ.copy()
    child_environment.update(deployment_environment)
    os.execvpe(sys.executable, argv, child_environment)


if __name__ == "__main__":
    main()
