"""Real-tool contract for Picard duplicate marking.

Picard MarkDuplicates marks duplicates and keeps every alignment in
final.bam. A PATH that has samtools but not Picard fails closed.
"""

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
    "@RG\tID:rg1\tSM:S1\tLB:lib1\tPL:ILLUMINA",
]
SEQ = "A" * 50
QUAL = "I" * 50


def _tool(name, env_var):
    candidate = resolve_tool(name, env_var)
    resolved = shutil.which(candidate)
    if not resolved and os.path.isabs(candidate) and os.access(candidate, os.X_OK):
        resolved = candidate
    if not resolved:
        require_external_tools([name], "legacy duplicate backend real execution")
    return resolved


@pytest.fixture(scope="module")
def snakemake_path():
    return _tool("snakemake", "SNAKEMAKE")


@pytest.fixture(scope="module")
def samtools_path():
    return _tool("samtools", "SAMTOOLS")


@pytest.fixture(scope="module")
def picard_path():
    return _tool("picard", "PICARD")


def _record(name, flag, pos, mate_pos, tlen):
    return "\t".join(
        [
            name,
            str(flag),
            "chr1",
            str(pos),
            "60",
            "50M",
            "=",
            str(mate_pos),
            str(tlen),
            SEQ,
            QUAL,
            "RG:Z:rg1",
        ]
    )


def _make_duplicate_bam(tmp_path, samtools_path):
    sam = tmp_path / "input.sam"
    raw = tmp_path / "raw.bam"
    bam = tmp_path / "results" / "T1" / "02_align" / "T1.mapq30.bam"
    bam.parent.mkdir(parents=True)
    sam.write_text(
        "\n".join(
            HEADER
            + [
                _record("r1", 99, 100, 200, 150),
                _record("r1", 147, 200, 100, -150),
                _record("r2", 99, 100, 200, 150),
                _record("r2", 147, 200, 100, -150),
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    _run([samtools_path, "view", "-b", "-o", str(raw), str(sam)])
    _run([samtools_path, "sort", "-o", str(bam), str(raw)])
    _run([samtools_path, "index", str(bam)])
    return bam


def _run(command, *, env=None, cwd=None, check=True):
    result = subprocess.run(command, capture_output=True, text=True, env=env, cwd=cwd)
    if check:
        assert result.returncode == 0, (
            f"Command failed: {command!r}\n"
            f"stdout:\n{result.stdout[-1000:]}\n"
            f"stderr:\n{result.stderr[-1000:]}"
        )
    return result


def _write_profile(tmp_path):
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
                "mapq: 30",
                'duplicate_analysis: "exclude"',
                "trim: true",
                "use_control: false",
                "multiqc: false",
                "blacklist_mode: off",
                "genome_resources:",
                "  hs:",
                '    effective_genome_size: "hs"',
                '    chrom_sizes: ""',
                '    blacklist: ""',
                '    gtf: ""',
                '    reference_fasta: ""',
                "",
            ]
        ),
        encoding="utf-8",
    )


def _run_marker(tmp_path, snakemake_path):
    return _run(
        [
            snakemake_path,
            "-s",
            str(SNAKEFILE),
            "--profile",
            "none",
            "--workflow-profile",
            "none",
            "--configfile",
            str(tmp_path / "config.yaml"),
            "-j",
            "1",
            "--allowed-rules",
            "analysis_identity_marker",
            "--",
            "results/.helixweave/analysis_identity.json",
        ],
        cwd=tmp_path,
    )


def _run_duplicate_rule(tmp_path, snakemake_path, *, path, check=True):
    env = os.environ.copy()
    env["PATH"] = f"{path}:/usr/bin:/bin"
    return _run(
        [
            snakemake_path,
            "-s",
            str(SNAKEFILE),
            "--profile",
            "none",
            "--workflow-profile",
            "none",
            "--configfile",
            str(tmp_path / "config.yaml"),
            "-j",
            "1",
            "results/T1/02_align/T1.final.bam",
        ],
        env=env,
        cwd=tmp_path,
        check=check,
    )


def test_picard_marks_duplicates_without_removing_them(
    tmp_path, snakemake_path, samtools_path, picard_path
):
    _write_profile(tmp_path)
    marker = _run_marker(tmp_path, snakemake_path)
    assert marker.returncode == 0, marker.stderr[-1000:]
    _make_duplicate_bam(tmp_path, samtools_path)

    bin_dir = tmp_path / "picard-bin"
    bin_dir.mkdir()
    (bin_dir / "samtools").symlink_to(samtools_path)
    (bin_dir / "picard").symlink_to(picard_path)
    conda = shutil.which("conda")
    if conda:
        (bin_dir / "conda").symlink_to(conda)

    result = _run_duplicate_rule(tmp_path, snakemake_path, path=str(bin_dir))
    assert result.returncode == 0, result.stderr[-1000:]

    metrics = (tmp_path / "results/T1/01_qc/T1.dup_metrics.txt").read_text(
        encoding="utf-8"
    )
    assert "picard.sam.DuplicationMetrics" in metrics
    assert "PERCENT_DUPLICATION" in metrics
    assert "samtools markdup" not in metrics
    final_bam = tmp_path / "results/T1/02_align/T1.final.bam"
    count = int(_run([samtools_path, "view", "-c", str(final_bam)]).stdout.strip())
    duplicate_count = int(
        _run(
            [samtools_path, "view", "-c", "-f", "0x400", str(final_bam)]
        ).stdout.strip()
    )
    assert count == 4
    assert duplicate_count > 0


def test_missing_picard_fails_without_samtools_fallback(
    tmp_path, snakemake_path, samtools_path
):
    _write_profile(tmp_path)
    marker = _run_marker(tmp_path, snakemake_path)
    assert marker.returncode == 0, marker.stderr[-1000:]
    _make_duplicate_bam(tmp_path, samtools_path)

    bin_dir = tmp_path / "samtools-bin"
    bin_dir.mkdir()
    (bin_dir / "samtools").symlink_to(samtools_path)
    conda = shutil.which("conda")
    if conda:
        (bin_dir / "conda").symlink_to(conda)

    result = _run_duplicate_rule(
        tmp_path, snakemake_path, path=str(bin_dir), check=False
    )
    output = result.stdout + result.stderr
    assert result.returncode != 0, output
    assert "Picard MarkDuplicates is required" in output
    assert "samtools markdup used" not in output
    assert not (tmp_path / "results/T1/02_align/T1.final.bam").exists()
