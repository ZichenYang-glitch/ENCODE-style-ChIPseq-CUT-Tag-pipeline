"""Real-tool regression for MACS3 extra_args migration boundaries (SCI-2).

The validator normalizes legacy --keep-dup tokens out of
tool_parameters.macs3.extra_args. This test proves against the real macs3
executable and the production macs3_callpeak rule that:

- a quoted value with spaces (e.g. --tempdir '/path/temp with spaces')
  survives normalization and still works when the production rule
  interpolates extra_args into its shell command;
- an abbreviated --keep-d 1 conflict is rejected at workflow startup
  instead of silently overriding the canonical --keep-dup all.

Usage:
    SNAKEMAKE=/path/to/snakemake MACS3=/path/to/macs3 \
    SAMTOOLS=/path/to/samtools \
        python3 -m pytest -m real_execution \
        test/real_execution/test_macs3_extra_args_migration.py -v
"""

import os
import shlex
import shutil
import subprocess
from pathlib import Path

import pytest
from _tool_resolver import require_external_tools, resolve_tool

from encode_pipeline.config.validate import validate_config

pytestmark = pytest.mark.real_execution

REPO_ROOT = Path(__file__).resolve().parents[2]
SNAKEFILE = REPO_ROOT / "workflow" / "Snakefile"

CHROM_LEN = 1_000_000
GENOME_SIZE = "10000000"
READ_LEN = 50
FRAGMENT_LEN = 200
N_FRAGMENTS = 60
CLUSTER_START = 100_000

SAMPLES_TSV = (
    "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\tpeak_mode\tgenome\t"
    "bowtie2_index\texperiment\tcondition\treplicate\tbiological_replicate\t"
    "technical_replicate\trole\tcontrol_sample\tcontrol_bam\n"
    "T1\tR1.fq\tR2.fq\tPE\tchipseq\tH3K27ac\tnarrow\ths\tbt2_idx\t"
    "exp1\tH3K27ac\t1\t1\t1\ttreatment\t\t\n"
)


def _tool(name, env_var):
    candidate = resolve_tool(name, env_var)
    resolved = shutil.which(candidate)
    if not resolved and os.path.isabs(candidate) and os.access(candidate, os.X_OK):
        resolved = candidate
    if not resolved:
        require_external_tools([name], "MACS3 extra_args migration real execution")
    return resolved


@pytest.fixture(scope="module")
def snakemake_path():
    return _tool("snakemake", "SNAKEMAKE")


@pytest.fixture(scope="module")
def macs3_path():
    return _tool("macs3", "MACS3")


@pytest.fixture(scope="module")
def samtools_path():
    return _tool("samtools", "SAMTOOLS")


def _run(command, *, env=None, cwd=None, check=True):
    result = subprocess.run(command, capture_output=True, text=True, env=env, cwd=cwd)
    if check:
        assert result.returncode == 0, (
            f"Command failed: {command!r}\n"
            f"stdout:\n{result.stdout[-1000:]}\n"
            f"stderr:\n{result.stderr[-1000:]}"
        )
    return result


def _make_pe_bam(bam_path, samtools_path):
    """Coordinate-sorted PE BAM with one cluster of proper pairs."""
    lines = [
        "@HD\tVN:1.6\tSO:coordinate",
        f"@SQ\tSN:chr1\tLN:{CHROM_LEN}",
    ]
    seq = "A" * READ_LEN
    qual = "I" * READ_LEN
    for i in range(N_FRAGMENTS):
        pos1 = CLUSTER_START + i * 10
        pos2 = pos1 + FRAGMENT_LEN - READ_LEN
        lines.append(
            f"frag{i}\t99\tchr1\t{pos1}\t60\t{READ_LEN}M\t=\t{pos2}\t{FRAGMENT_LEN}\t{seq}\t{qual}"
        )
        lines.append(
            f"frag{i}\t147\tchr1\t{pos2}\t60\t{READ_LEN}M\t=\t{pos1}\t{-FRAGMENT_LEN}\t{seq}\t{qual}"
        )
    sam_path = bam_path.with_suffix(".sam")
    sam_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    raw_path = bam_path.with_suffix(".raw.bam")
    _run([samtools_path, "view", "-b", "-o", str(raw_path), str(sam_path)])
    _run([samtools_path, "sort", "-o", str(bam_path), str(raw_path)])
    return bam_path


def _normalized_macs3_extra_args(tmp_path, extra):
    """Run the real config validator the Snakefile uses, return extra_args."""
    samples_path = tmp_path / "samples.tsv"
    samples_path.write_text(SAMPLES_TSV, encoding="utf-8")
    config = {
        "samples": str(samples_path),
        "use_control": False,
        "tool_parameters": {"macs3": {"extra_args": extra}},
    }
    return validate_config(config)["tool_parameters"]["macs3"]["extra_args"]


def test_quoted_tempdir_survives_normalization_into_shell_command(
    tmp_path, macs3_path, samtools_path
):
    """Normalized extra_args keeps shell argument boundaries for real macs3.

    Mirrors the production rule: extra_args is interpolated into the shell
    command line verbatim, exactly like {params.extra} in macs3_callpeak.
    """
    tempdir = tmp_path / "temp dir with spaces"
    tempdir.mkdir()
    extra = _normalized_macs3_extra_args(
        tmp_path, f"--tempdir '{tempdir}' --keep-dup all"
    )
    # The synonymous --keep-dup all is migrated out; the quoted path keeps
    # its protecting quotes so the shell still sees one argument.
    assert "--keep-dup" not in shlex.split(extra)
    assert shlex.split(extra) == ["--tempdir", str(tempdir)]

    bam = _make_pe_bam(tmp_path / "mini.bam", samtools_path)
    outdir = tmp_path / "peaks"
    outdir.mkdir()
    command = (
        f"{shlex.quote(macs3_path)} callpeak -t {shlex.quote(str(bam))} "
        f"-f BAMPE -g {GENOME_SIZE} -q 0.05 -B --nomodel "
        f"--keep-dup all -n test --outdir {shlex.quote(str(outdir))} {extra}"
    )
    result = subprocess.run(
        command, shell=True, capture_output=True, text=True, executable="/bin/bash"
    )
    assert result.returncode == 0, (
        f"macs3 via shell failed: {command!r}\nstderr:\n{result.stderr[-1000:]}"
    )
    assert (outdir / "test_peaks.narrowPeak").exists()


def _write_project(tmp_path, extra_args):
    (tmp_path / "samples.tsv").write_text(SAMPLES_TSV, encoding="utf-8")
    (tmp_path / "config.yaml").write_text(
        "\n".join(
            [
                "samples: samples.tsv",
                "outdir: results",
                "threads: 1",
                "mapq: 30",
                'remove_dup: "no"',
                "trim: true",
                "use_control: false",
                "multiqc: false",
                "tool_parameters:",
                "  macs3:",
                f'    extra_args: "{extra_args}"',
                "genome_resources:",
                "  hs:",
                '    effective_genome_size: "10000000"',
                '    chrom_sizes: ""',
                '    blacklist: ""',
                '    gtf: ""',
                '    reference_fasta: ""',
                "",
            ]
        ),
        encoding="utf-8",
    )


def _snakemake_env(*tool_paths):
    env = os.environ.copy()
    dirs = [str(Path(tool).parent) for tool in tool_paths]
    # Rules carry conda: directives; snakemake probes the conda frontend at
    # DAG construction even without --use-conda, so keep it reachable.
    conda = shutil.which("conda")
    if conda:
        dirs.append(str(Path(conda).parent))
    env["PATH"] = ":".join([*dirs, "/usr/bin", "/bin"])
    return env


def _snakemake_run(snakemake_path, tmp_path, target, env):
    return _run(
        [
            snakemake_path,
            "-s",
            str(SNAKEFILE),
            "--configfile",
            str(tmp_path / "config.yaml"),
            "-j",
            "1",
            target,
        ],
        env=env,
        cwd=tmp_path,
        check=False,
    )


def test_production_macs3_rule_accepts_quoted_tempdir(
    tmp_path, snakemake_path, macs3_path, samtools_path
):
    """The production macs3_callpeak rule runs real macs3 with a quoted
    space-containing --tempdir from extra_args (legacy profile, staged
    final.bam)."""
    tempdir = tmp_path / "temp dir with spaces"
    tempdir.mkdir()
    _write_project(tmp_path, f"--nomodel --extsize 200 --tempdir '{tempdir}'")
    env = _snakemake_env(macs3_path, samtools_path)

    marker = "results/.helixweave/analysis_identity.json"
    result = _snakemake_run(snakemake_path, tmp_path, marker, env)
    assert result.returncode == 0, (
        f"marker run failed:\n{result.stdout[-800:]}\n{result.stderr[-800:]}"
    )

    # Stage the analysis input after the run identity exists so the SCI-1
    # outdir guard sees a marked directory.
    final_bam = tmp_path / "results" / "T1" / "02_align" / "T1.final.bam"
    final_bam.parent.mkdir(parents=True, exist_ok=True)
    _make_pe_bam(final_bam, samtools_path)
    _run([samtools_path, "index", str(final_bam)])

    result = _snakemake_run(snakemake_path, tmp_path, "results/T1/04_peaks/T1", env)
    assert result.returncode == 0, (
        f"macs3_callpeak failed:\n{result.stdout[-800:]}\n{result.stderr[-800:]}"
    )
    assert (
        tmp_path / "results" / "T1" / "04_peaks" / "T1" / "T1_peaks.narrowPeak"
    ).exists()


def test_production_run_rejects_abbreviated_keep_dup_conflict(
    tmp_path, snakemake_path, macs3_path, samtools_path
):
    """--keep-d 1 is an argparse abbreviation of --keep-dup and conflicts
    with the canonical --keep-dup all; startup validation must refuse it."""
    _write_project(tmp_path, "--keep-d 1")
    env = _snakemake_env(macs3_path, samtools_path)
    result = _snakemake_run(
        snakemake_path,
        tmp_path,
        "results/.helixweave/analysis_identity.json",
        env,
    )
    assert result.returncode != 0
    output = result.stdout + result.stderr
    assert "conflicts with the BAM-layer" in output
