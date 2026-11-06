"""Real-tool contract for the opt-in fragment-occupancy TSS profile."""

import gzip
import json
import os
import shutil
import subprocess
from pathlib import Path

import pytest
from _tool_resolver import require_external_tools, resolve_tool

pytestmark = pytest.mark.real_execution

REPO_ROOT = Path(__file__).resolve().parents[2]
SNAKEFILE = REPO_ROOT / "workflow" / "Snakefile"
HEADER = [
    "@HD\tVN:1.6\tSO:coordinate",
    "@SQ\tSN:chr1\tLN:10000",
]
SEQ = "A" * 50
QUAL = "I" * 50


def _tool(name, env_var):
    candidate = resolve_tool(name, env_var)
    resolved = shutil.which(candidate)
    if not resolved and os.path.isabs(candidate) and os.access(candidate, os.X_OK):
        resolved = candidate
    if not resolved:
        require_external_tools([name], "TSS fragment real execution")
    return resolved


@pytest.fixture(scope="module")
def samtools_path():
    return _tool("samtools", "SAMTOOLS")


@pytest.fixture(scope="module")
def snakemake_path():
    return _tool("snakemake", "SNAKEMAKE")


def _run(command, *, cwd=None, env=None):
    result = subprocess.run(command, capture_output=True, text=True, cwd=cwd, env=env)
    assert result.returncode == 0, (
        f"Command failed: {command!r}\n"
        f"stdout:\n{result.stdout[-2000:]}\n"
        f"stderr:\n{result.stderr[-2000:]}"
    )
    return result


def _pair(name, start, mate_start, tlen):
    return [
        "\t".join(
            [
                name,
                "99",
                "chr1",
                str(start),
                "60",
                "50M",
                "=",
                str(mate_start),
                str(tlen),
                SEQ,
                QUAL,
            ]
        ),
        "\t".join(
            [
                name,
                "147",
                "chr1",
                str(mate_start),
                "60",
                "50M",
                "=",
                str(start),
                str(-tlen),
                SEQ,
                QUAL,
            ]
        ),
    ]


def test_fragment_tss_profile_runs_end_to_end(tmp_path, samtools_path, snakemake_path):
    gtf = tmp_path / "annotation.gtf"
    gtf.write_text(
        'chr1\ttest\ttranscript\t1\t100\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n',
        encoding="utf-8",
    )
    (tmp_path / "samples.tsv").write_text(
        "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\tpeak_mode\tgenome\t"
        "bowtie2_index\texperiment\tcondition\treplicate\tbiological_replicate\t"
        "technical_replicate\trole\tcontrol_sample\tcontrol_bam\n"
        "T1\tR1.fq\tR2.fq\tPE\tcuttag\tH3K4me3\tnarrow\ths\tbt2_idx\t"
        "exp1\tH3K4me3\t1\t1\t1\ttreatment\t\t\n",
        encoding="utf-8",
    )
    (tmp_path / "config.yaml").write_text(
        "\n".join(
            [
                "samples: samples.tsv",
                "outdir: results",
                "threads: 1",
                "use_control: false",
                "multiqc: false",
                "replicate_analysis: false",
                "extend_reads: auto",
                "binsize: 1",
                "genome_resources:",
                "  hs:",
                '    effective_genome_size: "hs"',
                f'    gtf: "{gtf}"',
                "config_schema_version: 2.0.0",
                "analysis_profile: updated",
                "fragment_cpm: auto",
                "tss_fragment_profile: true",
                "organelle_mode: auto",
                "blacklist_mode: off",
                "duplicate_backend: picard",
                "duplicate_analysis: exclude",
                "control_mode: matched",
                "reference_mode: legacy",
                "external_control_policy: warn",
                "qc:",
                "  summary: false",
                "  frip: false",
                "  signal_tracks: false",
                "  tss_enrichment: false",
                "",
            ]
        ),
        encoding="utf-8",
    )
    _run(
        [
            snakemake_path,
            "-s",
            str(SNAKEFILE),
            "--configfile",
            "config.yaml",
            "--workflow-profile",
            "none",
            "-j",
            "1",
            "results/.helixweave/analysis_identity.json",
        ],
        cwd=tmp_path,
    )
    sam = tmp_path / "input.sam"
    raw = tmp_path / "raw.bam"
    bam = tmp_path / "results/T1/02_align/T1.analysis.bam"
    bam.parent.mkdir(parents=True)
    sam.write_text(
        "\n".join(HEADER + _pair("frag1", 100, 300, 250)) + "\n",
        encoding="utf-8",
    )
    _run([samtools_path, "view", "-b", "-o", str(raw), str(sam)])
    _run([samtools_path, "sort", "-o", str(bam), str(raw)])
    _run([samtools_path, "index", str(bam)])
    for path in (bam, Path(str(bam) + ".bai")):
        os.utime(path, (1, 1))
    _run(
        [
            snakemake_path,
            "-s",
            str(SNAKEFILE),
            "--configfile",
            "config.yaml",
            "--workflow-profile",
            "none",
            "-j",
            "1",
            "results/T1/05_qc/tss/T1.fragment.tss_profile.tsv",
        ],
        cwd=tmp_path,
    )
    for relative in (
        "05_qc/tss/T1.fragment.tss_matrix.gz",
        "05_qc/tss/T1.fragment.tss_profile.tsv",
        "05_qc/tss/T1.fragment.tss_profile.pdf",
    ):
        path = tmp_path / "results/T1" / relative
        assert path.exists() and path.stat().st_size > 0

    matrix = tmp_path / "results/T1/05_qc/tss/T1.fragment.tss_matrix.gz"
    with gzip.open(matrix, "rt", encoding="utf-8") as handle:
        header = json.loads(handle.readline().lstrip("@"))
        data = handle.readline().rstrip("\n").split("\t")
    assert header["upstream"] == [3000]
    assert header["downstream"] == [3000]
    assert header["bin size"] == [10]
    assert header["sample_labels"] == ["T1.fragment.CPM"]
    assert header["sample_boundaries"] == [0, 600]
    values = [float(value) for value in data[6:]]
    assert len(values) == 600
    assert any(value > 0 for value in values)
