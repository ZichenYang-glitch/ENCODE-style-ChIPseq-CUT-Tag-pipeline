"""Production-rule evidence for general PE analysis-template counts.

Each case seeds a fresh run identity, writes real synthetic inputs with
samtools, and then executes the production ``analysis_template_metrics`` and
``analysis_policy`` rules.  No direct writer call stands in for a production
rule.  The cases prove that ChIP-seq, ATAC-seq and CUT&Tag with fragment QC
disabled all record analysis template counts, while SE is explicitly not
applicable.
"""

from __future__ import annotations

import csv
from decimal import Decimal
import json
import os
import shutil
import subprocess
from pathlib import Path

import pytest

from _tool_resolver import require_external_tools, resolve_tool
from conftest import seed_run_identity

from encode_pipeline.adapters.encode import EncodeStyleWorkflowAdapter
from encode_pipeline.platform.adapters import (
    QcSourceArtifact,
    QcSourceDocument,
    WorkflowInputs,
)
from encode_pipeline.qc.contract import (
    ANALYSIS_POLICY_V2_INPUT_STAGE_KEYS,
    ANALYSIS_TEMPLATE_METRICS_COLUMNS,
    MITO_QC_V2_COLUMNS,
)

pytestmark = pytest.mark.real_execution

REPO = Path(__file__).resolve().parents[2]
SNAKEFILE = REPO / "workflow" / "Snakefile"
SAMPLE = "S1"

SAMPLE_COLUMNS = (
    "sample",
    "fastq_1",
    "fastq_2",
    "layout",
    "assay",
    "target",
    "peak_mode",
    "genome",
    "bowtie2_index",
    "experiment",
    "biological_replicate",
)

ASSAY_TARGET = {
    "chipseq": ("H3K27ac", "narrow"),
    "atac": ("ATAC", "narrow"),
    "cuttag": ("H3K4me3", "narrow"),
}

SAM_HEADER = [
    "@HD\tVN:1.6\tSO:coordinate",
    "@SQ\tSN:chr1\tLN:100000",
    "@SQ\tSN:chr2\tLN:100000",
]

# Four valid templates (one with legal negative TLEN first) and four invalid
# templates covering each invalid reason exactly once.
ANALYSIS_RECORDS = [
    "v1\t99\tchr1\t100\t60\t50M\t=\t250\t150\t*\t*",
    "v1\t147\tchr1\t250\t60\t50M\t=\t100\t-150\t*\t*",
    "v2\t99\tchr1\t400\t60\t50M\t=\t650\t250\t*\t*",
    "v2\t147\tchr1\t650\t60\t50M\t=\t400\t-250\t*\t*",
    "v3\t99\tchr1\t800\t60\t50M\t=\t1200\t400\t*\t*",
    "v3\t147\tchr1\t1200\t60\t50M\t=\t800\t-400\t*\t*",
    "n1\t147\tchr1\t3800\t60\t50M\t=\t4000\t-200\t*\t*",
    "n1\t99\tchr1\t4000\t60\t50M\t=\t3800\t200\t*\t*",
    "z1\t99\tchr1\t1500\t60\t50M\t=\t1600\t0\t*\t*",
    "z1\t147\tchr1\t1600\t60\t50M\t=\t1500\t0\t*\t*",
    "o1\t99\tchr1\t1800\t60\t50M\t=\t1950\t150\t*\t*",
    "c1\t99\tchr1\t2100\t60\t50M\tchr2\t2300\t200\t*\t*",
    "c1\t147\tchr2\t2300\t60\t50M\tchr1\t2100\t-200\t*\t*",
    "i1\t99\tchr1\t2500\t60\t50M\t=\t2700\t200\t*\t*",
    "i1\t147\tchr1\t2700\t60\t50M\t=\t2500\t-250\t*\t*",
]


@pytest.fixture(scope="module")
def samtools_path() -> str:
    candidate = resolve_tool("samtools", "SAMTOOLS")
    resolved = shutil.which(candidate)
    if not resolved and os.path.isabs(candidate) and os.access(candidate, os.X_OK):
        resolved = candidate
    if not resolved:
        require_external_tools(["samtools"], "analysis-template real execution")
    return resolved


@pytest.fixture(scope="module")
def snakemake_executable() -> str:
    candidate = resolve_tool("snakemake", "SNAKEMAKE")
    resolved = shutil.which(candidate)
    if not resolved and os.path.isabs(candidate) and os.access(candidate, os.X_OK):
        resolved = candidate
    if not resolved:
        require_external_tools(["snakemake"], "analysis-template real execution")
    return resolved


def _run(argv, cwd: Path):
    result = subprocess.run(
        [str(value) for value in argv],
        cwd=cwd,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, (
        f"Command failed: {argv!r}\n"
        f"stdout:\n{result.stdout[-3000:]}\n"
        f"stderr:\n{result.stderr[-3000:]}"
    )
    return result


def _write_tsv(path: Path, columns, row) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerow({column: row.get(column, "NA") for column in columns})


def _write_project(
    tmp_path: Path,
    *,
    assay: str,
    layout: str,
    cuttag_fragment_size: bool,
    profile: str,
) -> Path:
    target, peak_mode = ASSAY_TARGET[assay]
    samples = tmp_path / "config" / "samples.tsv"
    samples.parent.mkdir(parents=True)
    with samples.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=SAMPLE_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerow(
            {
                "sample": SAMPLE,
                "fastq_1": "R1.fq.gz",
                "fastq_2": "" if layout == "SE" else "R2.fq.gz",
                "layout": layout,
                "assay": assay,
                "target": target,
                "peak_mode": peak_mode,
                "genome": "hs",
                "bowtie2_index": "idx",
                "experiment": "E1",
                "biological_replicate": "1",
            }
        )
    config = {
        "samples": "config/samples.tsv",
        "outdir": "results",
        "threads": 1,
        "use_control": False,
        "multiqc": False,
        "replicate_analysis": False,
        "analysis_profile": profile,
        "config_schema_version": "2.0.0",
        "blacklist_mode": "off",
        "genome_resources": {"hs": {"effective_genome_size": "hs"}},
        "qc": {
            "frip": False,
            "signal_tracks": False,
            "summary": False,
            "cuttag_fragment_size": cuttag_fragment_size,
            "cross_correlation": False,
        },
    }
    if profile == "updated":
        config["duplicate_backend"] = "picard"
        config["control_mode"] = "matched"
    config_path = tmp_path / "config" / "config.yaml"
    config_path.write_text(json.dumps(config, indent=2), encoding="utf-8")
    return config_path


def _write_stage_metrics(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(
            {
                "total_alignments": 6,
                "kept_alignments": 2,
                "dropped_alignments": 4,
                "dropped_templates": 2,
                "dropped_alignments_by_reason": {},
                "dropped_templates_by_reason": {},
            }
        ),
        encoding="utf-8",
    )


def _write_mito(path: Path) -> None:
    _write_tsv(
        path,
        MITO_QC_V2_COLUMNS,
        {
            "schema_version": "2",
            "sample": SAMPLE,
            "organelle_reads": "2",
            "total_mapped_reads": "8",
            "fraction": "0.250000",
            "mode": "report",
            "status": "ok",
            "input_stage": "pre_duplicate_mapq_bam",
            "matched_contigs": "chrM",
            "organelle_chroms_explicit": "false",
        },
    )


def _write_analysis_bam(path: Path, samtools_path: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    sam = path.with_suffix(".sam")
    sam.write_text("\n".join(SAM_HEADER + ANALYSIS_RECORDS) + "\n", encoding="utf-8")
    _run([samtools_path, "view", "-bS", "-o", path, sam], cwd=path.parent)
    sam.unlink()


def _run_policy(
    tmp_path: Path,
    snakemake_executable: str,
    config_path: Path,
    *,
    run_fragment_rule: bool,
):
    rules = ["analysis_template_metrics", "analysis_policy"]
    if run_fragment_rule:
        rules.append("cuttag_fragment_size")
    target = f"results/{SAMPLE}/05_qc/{SAMPLE}.analysis_policy.json"
    return _run(
        [
            snakemake_executable,
            "-s",
            SNAKEFILE,
            "--workflow-profile",
            "none",
            "--directory",
            tmp_path,
            "--configfile",
            config_path,
            "--cores",
            "1",
            "--printshellcmds",
            "--allowed-rules",
            *rules,
            "--",
            target,
        ],
        cwd=tmp_path,
    )


def _read_row(path: Path) -> dict:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        assert list(reader.fieldnames or []) == list(ANALYSIS_TEMPLATE_METRICS_COLUMNS)
        rows = list(reader)
    assert len(rows) == 1
    return rows[0]


def _policy_source(policy_path: Path, assay: str) -> QcSourceDocument:
    return QcSourceDocument(
        source=QcSourceArtifact(
            artifact_id="analysis-policy",
            output_type="analysis_policy",
            relative_path=f"results/{SAMPLE}/05_qc/{SAMPLE}.analysis_policy.json",
            metadata={"scope": "sample", "sample_id": SAMPLE, "assay": assay},
        ),
        content=policy_path.read_bytes(),
    )


@pytest.mark.parametrize("profile", ["legacy", "updated"])
@pytest.mark.parametrize(
    ("assay", "cuttag_fragment_size"),
    [
        ("chipseq", False),
        ("atac", False),
        ("cuttag", False),
        ("cuttag", True),
    ],
)
def test_pe_assays_record_general_analysis_template_counts(
    tmp_path,
    snakemake_executable,
    samtools_path,
    assay,
    cuttag_fragment_size,
    profile,
):
    config_path = _write_project(
        tmp_path,
        assay=assay,
        layout="PE",
        cuttag_fragment_size=cuttag_fragment_size,
        profile=profile,
    )
    seed_run_identity(tmp_path, snakemake_executable, SNAKEFILE, config_path)
    base = tmp_path / "results" / SAMPLE
    _write_stage_metrics(
        base / "05_qc" / f"{SAMPLE}.blacklist_dropped.stage_metrics.json"
    )
    _write_mito(base / "05_qc" / f"{SAMPLE}.mito_fraction.tsv")
    _write_analysis_bam(base / "02_align" / f"{SAMPLE}.analysis.bam", samtools_path)

    _run_policy(
        tmp_path,
        snakemake_executable,
        config_path,
        run_fragment_rule=cuttag_fragment_size,
    )

    component_path = base / "05_qc" / f"{SAMPLE}.analysis_template_metrics.tsv"
    component = _read_row(component_path)
    assert component["schema_version"] == "2"
    assert component["sample"] == SAMPLE
    assert component["layout"] == "PE"
    assert component["analysis_status"] == "ok"
    assert component["input_stage"] == "analysis_bam"
    assert component["analysis_total_templates"] == "8"
    assert component["analysis_valid_templates"] == "4"
    assert component["analysis_invalid_templates"] == "4"
    assert component["analysis_tlen_zero_templates"] == "1"
    assert component["analysis_orphan_templates"] == "1"
    assert component["analysis_cross_contig_templates"] == "1"
    assert component["analysis_inconsistent_templates"] == "1"
    assert component["alignments_total"] == "15"
    assert component["upstream_dropped_templates"] == "2"
    assert component["upstream_dropped_alignments"] == "4"
    assert int(component["analysis_valid_templates"]) + int(
        component["analysis_invalid_templates"]
    ) == int(component["analysis_total_templates"])

    policy_path = base / "05_qc" / f"{SAMPLE}.analysis_policy.json"
    policy = json.loads(policy_path.read_text(encoding="utf-8"))
    assert policy["schema_version"] == 2
    assert policy["profile"] == profile
    assert policy["analysis_template_status"] == "ok"
    assert set(policy["input_stages"]) == set(ANALYSIS_POLICY_V2_INPUT_STAGE_KEYS)
    assert policy["input_stages"]["analysis_template_metrics"] == "analysis_bam"
    assert policy["stage_counts"]["analysis_valid_templates"] == 4
    assert policy["stage_counts"]["analysis_invalid_templates"] == 4
    assert policy["stage_counts"]["upstream_dropped_templates"] == 2
    assert policy["stage_counts"]["upstream_dropped_alignments"] == 4

    if cuttag_fragment_size:
        assert (base / "01_qc" / f"{SAMPLE}.cuttag_fragment_size.tsv").is_file()

    result = EncodeStyleWorkflowAdapter().extract_qc_metrics(
        WorkflowInputs(config={}), (_policy_source(policy_path, assay),)
    )
    assert result.is_success
    metrics = {item.metric_key: item.value for item in result.value}
    assert metrics["analysis.stage.valid_templates"] == Decimal("4")
    assert metrics["analysis.stage.invalid_templates"] == Decimal("4")


@pytest.mark.parametrize("profile", ["legacy", "updated"])
def test_single_end_is_not_applicable_without_fabricated_zero(
    tmp_path,
    snakemake_executable,
    samtools_path,
    profile,
):
    config_path = _write_project(
        tmp_path,
        assay="chipseq",
        layout="SE",
        cuttag_fragment_size=False,
        profile=profile,
    )
    seed_run_identity(tmp_path, snakemake_executable, SNAKEFILE, config_path)
    base = tmp_path / "results" / SAMPLE
    _write_stage_metrics(
        base / "05_qc" / f"{SAMPLE}.blacklist_dropped.stage_metrics.json"
    )
    _write_mito(base / "05_qc" / f"{SAMPLE}.mito_fraction.tsv")
    _write_analysis_bam(base / "02_align" / f"{SAMPLE}.analysis.bam", samtools_path)

    _run_policy(
        tmp_path,
        snakemake_executable,
        config_path,
        run_fragment_rule=False,
    )

    component = _read_row(base / "05_qc" / f"{SAMPLE}.analysis_template_metrics.tsv")
    assert component["layout"] == "SE"
    assert component["analysis_status"] == "not_applicable"
    for column in (
        "analysis_total_templates",
        "analysis_valid_templates",
        "analysis_invalid_templates",
        "analysis_tlen_zero_templates",
        "analysis_orphan_templates",
        "analysis_cross_contig_templates",
        "analysis_inconsistent_templates",
        "alignments_total",
    ):
        assert component[column] == "NA"

    policy_path = base / "05_qc" / f"{SAMPLE}.analysis_policy.json"
    policy = json.loads(policy_path.read_text(encoding="utf-8"))
    assert policy["analysis_template_status"] == "not_applicable"
    assert policy["input_stages"]["analysis_template_metrics"] == "not_applicable"
    assert policy["stage_counts"]["analysis_valid_templates"] is None
    assert policy["stage_counts"]["analysis_invalid_templates"] is None

    result = EncodeStyleWorkflowAdapter().extract_qc_metrics(
        WorkflowInputs(config={}), (_policy_source(policy_path, "chipseq"),)
    )
    assert result.is_success
    metrics = {item.metric_key: item.value for item in result.value}
    assert "analysis.stage.valid_templates" not in metrics
    assert "analysis.stage.invalid_templates" not in metrics
