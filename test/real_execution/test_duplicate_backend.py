"""Real-tool contract for legacy duplicate backend dispatch.

The legacy profile must keep the original Picard-available dispatch and the
samtools-only fallback. The updated profile uses Picard mark-only and is
covered by test_analysis_chain_b1.py.
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


def _run(command, *, env=None, cwd=None):
    result = subprocess.run(command, capture_output=True, text=True, env=env, cwd=cwd)
    assert result.returncode == 0, (
        f"Command failed: {command!r}\n"
        f"stdout:\n{result.stdout[-1000:]}\n"
        f"stderr:\n{result.stderr[-1000:]}"
    )
    return result


def _write_profile(tmp_path, *, remove_dup):
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
                f'remove_dup: "{remove_dup}"',
                "trim: true",
                "use_control: false",
                "multiqc: false",
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


def _run_duplicate_rule(tmp_path, snakemake_path, *, path):
    env = os.environ.copy()
    env["PATH"] = f"{path}:/usr/bin:/bin"
    return _run(
        [
            snakemake_path,
            "-s",
            str(SNAKEFILE),
            "--configfile",
            str(tmp_path / "config.yaml"),
            "-j",
            "1",
            "results/T1/02_align/T1.final.bam",
        ],
        env=env,
        cwd=tmp_path,
    )


@pytest.mark.parametrize(
    ("remove_dup", "with_picard"),
    [
        ("yes", True),
        ("yes", False),
        ("no", True),
        ("no", False),
    ],
)
def test_legacy_duplicate_dispatch(
    tmp_path, snakemake_path, samtools_path, picard_path, remove_dup, with_picard
):
    _make_duplicate_bam(tmp_path, samtools_path)
    _write_profile(tmp_path, remove_dup=remove_dup)

    if with_picard:
        tool_path = str(Path(picard_path).parent)
    else:
        bin_dir = tmp_path / "samtools-bin"
        bin_dir.mkdir()
        (bin_dir / "samtools").symlink_to(samtools_path)
        tool_path = str(bin_dir)

    _run_duplicate_rule(tmp_path, snakemake_path, path=tool_path)

    metrics = (tmp_path / "results/T1/01_qc/T1.dup_metrics.txt").read_text(
        encoding="utf-8"
    )
    final_bam = tmp_path / "results/T1/02_align/T1.final.bam"
    if with_picard:
        assert "picard.sam.DuplicationMetrics" in metrics
        assert "PERCENT_DUPLICATION" in metrics
    else:
        if remove_dup == "yes":
            assert "samtools markdup used (picard not found)." in metrics
        else:
            assert "picard not found; duplicate metrics unavailable." in metrics

    count = int(_run([samtools_path, "view", "-c", str(final_bam)]).stdout.strip())
    if remove_dup == "yes":
        assert count == 2
    else:
        assert count == 4
