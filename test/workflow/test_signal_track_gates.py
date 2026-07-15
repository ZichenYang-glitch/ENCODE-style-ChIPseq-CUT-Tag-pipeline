#!/usr/bin/env python3
"""Protect sample and pooled FE/ppois signal-track DAG gates.

Verifies that BigWig targets are correctly gated on:
- chrom_sizes empty → no BigWig targets
- chrom_sizes configured → BigWig targets appear
- chrom_sizes non-empty but file missing → validation fails
- qc.signal_tracks: false → no bedGraph or BigWig targets
- Pooled experiment with chrom_sizes → pooled BigWig targets appear
- Pooled experiment without chrom_sizes → no pooled BigWig targets

"""

import os
import subprocess
import sys

from _tool_resolver import resolve_tool


REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SNAKEFILE = os.path.join(REPO_ROOT, "workflow", "Snakefile")
VALIDATOR = os.path.join(REPO_ROOT, "scripts", "validate_samples.py")
SNAKEMAKE = resolve_tool("snakemake", "SNAKEMAKE")


def _env(workdir):
    e = os.environ.copy()
    e["PYTHONDONTWRITEBYTECODE"] = "1"
    e["XDG_CACHE_HOME"] = os.path.join(workdir, ".cache")
    return e


def _write_samples(workdir, with_experiment=False, with_control=False):
    """Write a minimal 2-sample samples.tsv with placeholder FASTQ paths."""
    lines = [
        "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\tpeak_mode\tgenome\tbowtie2_index",
    ]
    if with_experiment:
        lines[0] += "\texperiment\tbiological_replicate"
        lines.append(
            "S1\tS1_R1.fq\tS1_R2.fq\tPE\tchipseq\tCTCF\tnarrow\ths\t/opt/reference/bt2_idx\tEXP1\t1"
        )
        lines.append(
            "S2\tS2_R1.fq\tS2_R2.fq\tPE\tchipseq\tCTCF\tnarrow\ths\t/opt/reference/bt2_idx\tEXP1\t2"
        )
    elif with_control:
        lines.append(
            "S1\tS1_R1.fq\tS1_R2.fq\tPE\tchipseq\tCTCF\tnarrow\ths\t/opt/reference/bt2_idx"
        )
        lines[0] += "\trole"
        lines[1] += "\ttreatment"
        lines.append(
            "CTRL1\tCTRL1_R1.fq\tCTRL1_R2.fq\tPE\tchipseq\tCTCF\tnarrow\ths\t/opt/reference/bt2_idx"
        )
        lines[2] += "\tcontrol"
    else:
        lines.append(
            "S1\tS1_R1.fq\tS1_R2.fq\tPE\tchipseq\tCTCF\tnarrow\ths\t/opt/reference/bt2_idx"
        )

    path = os.path.join(workdir, "samples.tsv")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    return path


def _write_config(
    workdir,
    chrom_sizes_value='""',
    signal_tracks=True,
    stage4b="true",
    with_experiment=False,
    with_control=False,
):
    """Write a config.yaml that supports the test case."""
    signal_tracks_str = "true" if signal_tracks else "false"

    content = f"""samples: "{os.path.join(workdir, "samples.tsv")}"
outdir: "results"
threads: 1
mapq: 30
binsize: 10
remove_dup: "auto"
trim: true
extend_reads: "auto"
use_control: false
multiqc: false
stage4b: {stage4b}
stage5: false
qc:
  signal_tracks: {signal_tracks_str}
  blacklist_filter: false
  frip: false
  library_complexity: false
  nrf_pbc: false
  summary: false
  cuttag_fragment_size: false
genome_resources:
  hs:
    effective_genome_size: "hs"
    chrom_sizes: {chrom_sizes_value}
    blacklist: ""
"""

    path = os.path.join(workdir, "config.yaml")
    with open(path, "w") as f:
        f.write(content)
    return path


def _create_placeholder_fastqs(workdir):
    """Create empty placeholder FASTQs so Snakemake can resolve input paths."""
    for fname in [
        "S1_R1.fq",
        "S1_R2.fq",
        "S2_R1.fq",
        "S2_R2.fq",
        "CTRL1_R1.fq",
        "CTRL1_R2.fq",
    ]:
        fpath = os.path.join(workdir, fname)
        if not os.path.exists(fpath):
            with open(fpath, "w"):
                pass


def _run_snakemake_dryrun(workdir):
    """Run Snakemake dry-run. Returns (returncode, stdout+stderr)."""
    result = subprocess.run(
        [
            SNAKEMAKE,
            "-s",
            SNAKEFILE,
            "--configfile",
            os.path.join(workdir, "config.yaml"),
            "-n",
        ],
        cwd=workdir,
        capture_output=True,
        text=True,
        env=_env(workdir),
    )
    return result.returncode, result.stdout + result.stderr


def _run_validation(workdir):
    """Run validate_samples.py. Returns (returncode, stderr)."""
    result = subprocess.run(
        [sys.executable, VALIDATOR, "--config", os.path.join(workdir, "config.yaml")],
        cwd=workdir,
        capture_output=True,
        text=True,
        env=_env(workdir),
    )
    return result.returncode, result.stderr


def test_no_chrom_sizes_no_bigwig(tmp_path):
    workdir = str(tmp_path)
    _write_samples(workdir)
    _write_config(workdir, chrom_sizes_value='""', signal_tracks=True)
    _create_placeholder_fastqs(workdir)

    rc, output = _run_snakemake_dryrun(workdir)
    assert rc == 0, output
    assert "signal_track_fe_bw" not in output
    assert "signal_track_ppois_bw" not in output
    assert "signal_track_fe" in output
    assert "signal_track_ppois" in output


def test_chrom_sizes_configured_bigwig(tmp_path):
    workdir = str(tmp_path)
    chrom_sizes_path = tmp_path / "chrom.sizes"
    chrom_sizes_path.write_text("chr1\t1000000\n", encoding="utf-8")
    _write_samples(workdir)
    _write_config(
        workdir,
        chrom_sizes_value=f'"{chrom_sizes_path}"',
        signal_tracks=True,
    )
    _create_placeholder_fastqs(workdir)

    rc, output = _run_snakemake_dryrun(workdir)
    assert rc == 0, output
    assert "signal_track_fe_bw" in output
    assert "signal_track_ppois_bw" in output


def test_chrom_sizes_missing_file(tmp_path):
    workdir = str(tmp_path)
    _write_samples(workdir)
    _write_config(
        workdir,
        chrom_sizes_value='"/nonexistent/chrom.sizes"',
        signal_tracks=True,
    )

    rc, stderr = _run_validation(workdir)
    assert rc != 0
    assert "file not found" in stderr.lower()


def test_signal_tracks_false(tmp_path):
    workdir = str(tmp_path)
    chrom_sizes_path = tmp_path / "chrom.sizes"
    chrom_sizes_path.write_text("chr1\t1000000\n", encoding="utf-8")
    _write_samples(workdir)
    _write_config(
        workdir,
        chrom_sizes_value=f'"{chrom_sizes_path}"',
        signal_tracks=False,
    )
    _create_placeholder_fastqs(workdir)

    rc, output = _run_snakemake_dryrun(workdir)
    assert rc == 0, output
    assert "signal_track_fe" not in output
    assert "signal_track_ppois" not in output
    assert "signal_track_fe_bw" not in output
    assert "signal_track_ppois_bw" not in output


def test_pooled_chrom_sizes_pooled_bigwig(tmp_path):
    workdir = str(tmp_path)
    chrom_sizes_path = tmp_path / "chrom.sizes"
    chrom_sizes_path.write_text("chr1\t1000000\n", encoding="utf-8")
    _write_samples(workdir, with_experiment=True)
    _write_config(
        workdir,
        chrom_sizes_value=f'"{chrom_sizes_path}"',
        signal_tracks=True,
        stage4b="true",
        with_experiment=True,
    )
    _create_placeholder_fastqs(workdir)

    rc, output = _run_snakemake_dryrun(workdir)
    assert rc == 0, output
    assert "pooled_signal_track_fe_bw" in output
    assert "pooled_signal_track_ppois_bw" in output


def test_pooled_no_chrom_sizes_no_pooled_bigwig(tmp_path):
    workdir = str(tmp_path)
    _write_samples(workdir, with_experiment=True)
    _write_config(
        workdir,
        chrom_sizes_value='""',
        signal_tracks=True,
        stage4b="true",
        with_experiment=True,
    )
    _create_placeholder_fastqs(workdir)

    rc, output = _run_snakemake_dryrun(workdir)
    assert rc == 0, output
    assert "pooled_signal_track_fe_bw" not in output
    assert "pooled_signal_track_ppois_bw" not in output
    assert "pooled_signal_track_fe" in output
    assert "pooled_signal_track_ppois" in output
