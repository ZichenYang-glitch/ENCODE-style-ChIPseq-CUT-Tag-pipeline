"""Real Snakemake execution for the analysis_filtered BAM rule."""

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
    "@SQ\tSN:chrM\tLN:100",
]
SEQ = "A" * 50
QUAL = "I" * 50


def _tool(name, env_var):
    candidate = resolve_tool(name, env_var)
    resolved = shutil.which(candidate)
    if not resolved and os.path.isabs(candidate) and os.access(candidate, os.X_OK):
        resolved = candidate
    if not resolved:
        require_external_tools([name], "analysis_filtered rule execution")
    return resolved


@pytest.fixture(scope="module")
def snakemake_path():
    return _tool("snakemake", "SNAKEMAKE")


@pytest.fixture(scope="module")
def samtools_path():
    return _tool("samtools", "SAMTOOLS")


def _pair(name, chrom, start, mate_start, tlen, extra_flag=0):
    return [
        "\t".join(
            [
                name,
                str(99 + extra_flag),
                chrom,
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
                str(147 + extra_flag),
                chrom,
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


def _run(command, *, cwd=None):
    result = subprocess.run(command, capture_output=True, text=True, cwd=cwd)
    assert result.returncode == 0, (
        f"Command failed: {command!r}\n"
        f"stdout:\n{result.stdout[-2000:]}\n"
        f"stderr:\n{result.stderr[-2000:]}"
    )
    return result


def _write_bam(tmp_path, samtools_path):
    sam = tmp_path / "input.sam"
    raw = tmp_path / "raw.bam"
    bam = tmp_path / "results" / "T1" / "02_align" / "T1.final.bam"
    bam.parent.mkdir(parents=True)
    records = []
    records += _pair("normal", "chr1", 500, 600, 150)
    records += _pair("duplicate1", "chr1", 700, 800, 150, extra_flag=1024)
    records += _pair("duplicate2", "chr1", 700, 800, 150, extra_flag=1024)
    records += _pair("blacklist", "chr1", 5, 100, 145)
    records += _pair("mito", "chrM", 10, 100, 140)
    sam.write_text("\n".join(HEADER + records) + "\n", encoding="utf-8")
    _run([samtools_path, "view", "-b", "-o", str(raw), str(sam)])
    _run([samtools_path, "sort", "-o", str(bam), str(raw)])
    _run([samtools_path, "index", str(bam)])
    return bam


def test_prepare_analysis_bam_applies_all_filters(
    tmp_path, snakemake_path, samtools_path
):
    blacklist = tmp_path / "blacklist.bed"
    blacklist.write_text("chr1\t0\t10\n", encoding="utf-8")
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
                "genome_resources:",
                "  hs:",
                '    effective_genome_size: "hs"',
                f'    blacklist: "{blacklist}"',
                '    organelle_chroms: ["chrM"]',
                "config_schema_version: 2.0.0",
                "analysis_profile: updated",
                "blacklist_mode: analysis",
                "organelle_mode: auto",
                "duplicate_backend: picard",
                "duplicate_analysis: exclude",
                "control_mode: matched",
                "reference_mode: legacy",
                "external_control_policy: warn",
                'fragment_cpm: "off"',
                "tss_fragment_profile: false",
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
    _write_bam(tmp_path, samtools_path)
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
            "results/analysis_filtered/T1/02_align/T1.analysis.bam",
        ],
        cwd=tmp_path,
    )

    output_bam = tmp_path / "results/analysis_filtered/T1/02_align/T1.analysis.bam"
    records = _run([samtools_path, "view", str(output_bam)]).stdout.strip().splitlines()
    assert len(records) == 2
    assert {line.split("\t")[0] for line in records} == {"normal"}

    dropped = (
        tmp_path / "results/analysis_filtered/T1/05_qc/T1.blacklist_dropped.tsv"
    ).read_text(encoding="utf-8")
    fields = dropped.strip().splitlines()[1].split("\t")
    assert fields[0] == "T1"
    assert fields[1] == "PE"
    assert fields[2] == "1"
    assert (
        tmp_path / "results/analysis_filtered/T1/02_align/T1.analysis.bam.bai"
    ).exists()


def test_peak_counts_preserves_caller_and_filtered_counts(
    tmp_path, snakemake_path, samtools_path
):
    blacklist = tmp_path / "blacklist.bed"
    blacklist.write_text("chr1\t0\t10\n", encoding="utf-8")
    (tmp_path / "samples.tsv").write_text(
        "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\tpeak_mode\tgenome\t"
        "bowtie2_index\texperiment\tcondition\treplicate\tbiological_replicate\t"
        "technical_replicate\trole\tcontrol_sample\tcontrol_bam\n"
        "T1\tR1.fq\tR2.fq\tPE\tchipseq\tH3K27ac\tnarrow\ths\tbt2_idx\t"
        "exp1\tH3K27ac\t1\t1\t1\ttreatment\t\t\n",
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
                "genome_resources:",
                "  hs:",
                '    effective_genome_size: "hs"',
                f'    blacklist: "{blacklist}"',
                "config_schema_version: 2.0.0",
                "analysis_profile: updated",
                "blacklist_mode: analysis",
                "organelle_mode: auto",
                "duplicate_backend: picard",
                "duplicate_analysis: exclude",
                "control_mode: matched",
                "reference_mode: legacy",
                "external_control_policy: warn",
                'fragment_cpm: "off"',
                "tss_fragment_profile: false",
                "qc:",
                "  summary: false",
                "  frip: false",
                "  signal_tracks: false",
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
    peaks_dir = tmp_path / "results/analysis_filtered/T1/04_peaks/T1"
    filtered_dir = (
        tmp_path / "results/analysis_filtered/T1/04_peaks/T1_blacklist_filtered"
    )
    peaks_dir.mkdir(parents=True)
    filtered_dir.mkdir(parents=True)
    (peaks_dir / "T1_peaks.narrowPeak").write_text(
        "chr1\t0\t100\tbad\t100\t.\t1\t2\t3\t50\n"
        "chr1\t500\t600\tgood\t100\t.\t1\t2\t3\t50\n",
        encoding="utf-8",
    )
    (filtered_dir / "T1_peaks.blacklist_filtered.narrowPeak").write_text(
        "chr1\t500\t600\tgood\t100\t.\t1\t2\t3\t50\n",
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
            "results/analysis_filtered/T1/01_qc/T1.peak_counts.tsv",
        ],
        cwd=tmp_path,
    )
    rows = (
        (tmp_path / "results/analysis_filtered/T1/01_qc/T1.peak_counts.tsv")
        .read_text(encoding="utf-8")
        .strip()
        .splitlines()
    )
    fields = rows[1].split("\t")
    assert fields[0] == "T1"
    assert fields[1] == "narrow"
    # ``peaks`` always keeps the true caller count; the blacklist-filtered
    # count is reported separately for the delivered peak set.
    assert fields[2] == "2"
    assert fields[3] == "1"
