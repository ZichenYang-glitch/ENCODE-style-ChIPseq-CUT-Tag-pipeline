#!/usr/bin/env python3
"""Stage 22 stress tests for FE/ppois BigWig DAG gating.

Verifies that BigWig targets are correctly gated on:
- chrom_sizes empty → no BigWig targets
- chrom_sizes configured → BigWig targets appear
- chrom_sizes non-empty but file missing → validation fails
- qc.signal_tracks: false → no bedGraph or BigWig targets
- Pooled experiment with chrom_sizes → pooled BigWig targets appear
- Pooled experiment without chrom_sizes → no pooled BigWig targets

Usage:
    SNAKEMAKE=/path/to/snakemake python3 test/test_stage22_bigwig_stress.py
"""

import subprocess
import tempfile
import shutil
import os
import sys
from _tool_resolver import resolve_tool


REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SNAKEFILE = os.path.join(REPO_ROOT, "workflow", "Snakefile")
VALIDATOR = os.path.join(REPO_ROOT, "scripts", "validate_samples.py")
SNAKEMAKE = resolve_tool("snakemake", "SNAKEMAKE")

PASSED = 0
TOTAL = 0


def _env():
    e = os.environ.copy()
    e["PYTHONDONTWRITEBYTECODE"] = "1"
    return e


def _mk_workdir():
    return tempfile.mkdtemp(prefix="smoke_stage22_", dir="/tmp")


def _write_samples(workdir, with_experiment=False, with_control=False):
    """Write a minimal 2-sample samples.tsv with placeholder FASTQ paths."""
    lines = [
        "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\tpeak_mode\tgenome\tbowtie2_index",
    ]
    if with_experiment:
        lines[0] += "\texperiment\tbiological_replicate"
        lines.append("S1\tS1_R1.fq\tS1_R2.fq\tPE\tchipseq\tCTCF\tnarrow\ths\t/opt/reference/bt2_idx\tEXP1\t1")
        lines.append("S2\tS2_R1.fq\tS2_R2.fq\tPE\tchipseq\tCTCF\tnarrow\ths\t/opt/reference/bt2_idx\tEXP1\t2")
    elif with_control:
        lines.append("S1\tS1_R1.fq\tS1_R2.fq\tPE\tchipseq\tCTCF\tnarrow\ths\t/opt/reference/bt2_idx")
        lines[0] += "\trole"
        lines[1] += "\ttreatment"
        lines.append("CTRL1\tCTRL1_R1.fq\tCTRL1_R2.fq\tPE\tchipseq\tCTCF\tnarrow\ths\t/opt/reference/bt2_idx")
        lines[2] += "\tcontrol"
    else:
        lines.append("S1\tS1_R1.fq\tS1_R2.fq\tPE\tchipseq\tCTCF\tnarrow\ths\t/opt/reference/bt2_idx")

    path = os.path.join(workdir, "samples.tsv")
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")
    return path


def _write_config(workdir, chrom_sizes_value='""', signal_tracks=True, stage4b="true",
                   with_experiment=False, with_control=False):
    """Write a config.yaml that supports the test case."""
    signal_tracks_str = "true" if signal_tracks else "false"

    content = f"""samples: "{os.path.join(workdir, 'samples.tsv')}"
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
    for fname in ["S1_R1.fq", "S1_R2.fq", "S2_R1.fq", "S2_R2.fq",
                   "CTRL1_R1.fq", "CTRL1_R2.fq"]:
        fpath = os.path.join(workdir, fname)
        if not os.path.exists(fpath):
            with open(fpath, "w"):
                pass


def _run_snakemake_dryrun(workdir):
    """Run Snakemake dry-run. Returns (returncode, stdout+stderr)."""
    result = subprocess.run(
        [SNAKEMAKE, "-s", SNAKEFILE, "--configfile",
         os.path.join(workdir, "config.yaml"), "-n"],
        cwd=workdir,
        capture_output=True,
        text=True,
        env=_env(),
    )
    return result.returncode, result.stdout + result.stderr


def _run_validation(workdir):
    """Run validate_samples.py. Returns (returncode, stderr)."""
    result = subprocess.run(
        [sys.executable, VALIDATOR, "--config",
         os.path.join(workdir, "config.yaml")],
        cwd=workdir,
        capture_output=True,
        text=True,
        env=_env(),
    )
    return result.returncode, result.stderr


def _record(name, passed):
    global PASSED, TOTAL
    TOTAL += 1
    if passed:
        PASSED += 1
        print("PASS: %s" % name)
    else:
        print("FAIL: %s" % name)


# ---------------------------------------------------------------------------
# Test 1: chrom_sizes empty → no BigWig targets in DAG
# ---------------------------------------------------------------------------

def test_no_chrom_sizes_no_bigwig():
    workdir = _mk_workdir()
    try:
        _write_samples(workdir)
        _write_config(workdir, chrom_sizes_value='""', signal_tracks=True)
        _create_placeholder_fastqs(workdir)

        rc, output = _run_snakemake_dryrun(workdir)
        if rc != 0:
            _record("1-chrom_sizes_empty_no_bw_targets", False)
            return

        # BigWig rules should NOT appear
        has_fe_bw = "signal_track_fe_bw" in output
        has_ppois_bw = "signal_track_ppois_bw" in output
        # bedGraph rules SHOULD still appear
        has_fe_bdg = "signal_track_fe" in output
        has_ppois_bdg = "signal_track_ppois" in output

        passed = (not has_fe_bw and not has_ppois_bw
                  and has_fe_bdg and has_ppois_bdg)
        if not passed:
            print("  has_fe_bw=%s has_ppois_bw=%s has_fe_bdg=%s has_ppois_bdg=%s"
                  % (has_fe_bw, has_ppois_bw, has_fe_bdg, has_ppois_bdg))
        _record("1-chrom_sizes_empty_no_bw_targets", passed)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


# ---------------------------------------------------------------------------
# Test 2: chrom_sizes configured → BigWig targets appear
# ---------------------------------------------------------------------------

def test_chrom_sizes_configured_bigwig():
    workdir = _mk_workdir()
    try:
        # Use a valid file as chrom_sizes (config.yaml itself exists)
        chrom_sizes_path = os.path.join(workdir, "chrom.sizes")
        with open(chrom_sizes_path, "w") as f:
            f.write("chr1\t1000000\n")

        _write_samples(workdir)
        _write_config(workdir, chrom_sizes_value=f'"{chrom_sizes_path}"',
                      signal_tracks=True)
        _create_placeholder_fastqs(workdir)

        rc, output = _run_snakemake_dryrun(workdir)
        if rc != 0:
            _record("2-chrom_sizes_configured_bw_targets", False)
            return

        has_fe_bw = "signal_track_fe_bw" in output
        has_ppois_bw = "signal_track_ppois_bw" in output

        passed = has_fe_bw and has_ppois_bw
        if not passed:
            print("  has_fe_bw=%s has_ppois_bw=%s" % (has_fe_bw, has_ppois_bw))
        _record("2-chrom_sizes_configured_bw_targets", passed)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


# ---------------------------------------------------------------------------
# Test 3: chrom_sizes non-empty but file missing → validation fails
# ---------------------------------------------------------------------------

def test_chrom_sizes_missing_file():
    workdir = _mk_workdir()
    try:
        _write_samples(workdir)
        _write_config(workdir, chrom_sizes_value='"/nonexistent/chrom.sizes"',
                      signal_tracks=True)

        rc, stderr = _run_validation(workdir)

        passed = (rc != 0 and "file not found" in stderr.lower())
        if not passed:
            print("  rc=%s stderr=%s" % (rc, stderr[-200:]))
        _record("3-chrom_sizes_missing_file_validation_fails", passed)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


# ---------------------------------------------------------------------------
# Test 4: qc.signal_tracks: false → no bedGraph or BigWig targets
# ---------------------------------------------------------------------------

def test_signal_tracks_false():
    workdir = _mk_workdir()
    try:
        chrom_sizes_path = os.path.join(workdir, "chrom.sizes")
        with open(chrom_sizes_path, "w") as f:
            f.write("chr1\t1000000\n")

        _write_samples(workdir)
        _write_config(workdir, chrom_sizes_value=f'"{chrom_sizes_path}"',
                      signal_tracks=False)
        _create_placeholder_fastqs(workdir)

        rc, output = _run_snakemake_dryrun(workdir)
        if rc != 0:
            _record("4-signal_tracks_false_no_targets", False)
            return

        has_fe_bdg = "signal_track_fe" in output
        has_ppois_bdg = "signal_track_ppois" in output
        has_fe_bw = "signal_track_fe_bw" in output
        has_ppois_bw = "signal_track_ppois_bw" in output

        passed = (not has_fe_bdg and not has_ppois_bdg
                  and not has_fe_bw and not has_ppois_bw)
        if not passed:
            print("  has_fe_bdg=%s has_ppois_bdg=%s has_fe_bw=%s has_ppois_bw=%s"
                  % (has_fe_bdg, has_ppois_bdg, has_fe_bw, has_ppois_bw))
        _record("4-signal_tracks_false_no_targets", passed)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


# ---------------------------------------------------------------------------
# Test 5: Pooled experiment with chrom_sizes → pooled BigWig targets appear
# ---------------------------------------------------------------------------

def test_pooled_chrom_sizes_pooled_bigwig():
    workdir = _mk_workdir()
    try:
        chrom_sizes_path = os.path.join(workdir, "chrom.sizes")
        with open(chrom_sizes_path, "w") as f:
            f.write("chr1\t1000000\n")

        _write_samples(workdir, with_experiment=True)
        _write_config(workdir, chrom_sizes_value=f'"{chrom_sizes_path}"',
                      signal_tracks=True, stage4b="true", with_experiment=True)
        _create_placeholder_fastqs(workdir)

        rc, output = _run_snakemake_dryrun(workdir)
        if rc != 0:
            _record("5-pooled_chrom_sizes_pooled_bw", False)
            return

        has_pooled_fe_bw = "pooled_signal_track_fe_bw" in output
        has_pooled_ppois_bw = "pooled_signal_track_ppois_bw" in output

        passed = has_pooled_fe_bw and has_pooled_ppois_bw
        if not passed:
            print("  has_pooled_fe_bw=%s has_pooled_ppois_bw=%s"
                  % (has_pooled_fe_bw, has_pooled_ppois_bw))
        _record("5-pooled_chrom_sizes_pooled_bw", passed)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


# ---------------------------------------------------------------------------
# Test 6: Pooled experiment without chrom_sizes → no pooled BigWig targets
# ---------------------------------------------------------------------------

def test_pooled_no_chrom_sizes_no_pooled_bigwig():
    workdir = _mk_workdir()
    try:
        _write_samples(workdir, with_experiment=True)
        _write_config(workdir, chrom_sizes_value='""', signal_tracks=True,
                      stage4b="true", with_experiment=True)
        _create_placeholder_fastqs(workdir)

        rc, output = _run_snakemake_dryrun(workdir)
        if rc != 0:
            _record("6-pooled_no_chrom_sizes_no_pooled_bw", False)
            return

        has_pooled_fe_bw = "pooled_signal_track_fe_bw" in output
        has_pooled_ppois_bw = "pooled_signal_track_ppois_bw" in output
        # bedGraph pooled tracks SHOULD appear
        has_pooled_fe_bdg = "pooled_signal_track_fe" in output
        has_pooled_ppois_bdg = "pooled_signal_track_ppois" in output

        passed = (not has_pooled_fe_bw and not has_pooled_ppois_bw
                  and has_pooled_fe_bdg and has_pooled_ppois_bdg)
        if not passed:
            print("  has_pooled_fe_bw=%s has_pooled_ppois_bw=%s"
                  " has_pooled_fe_bdg=%s has_pooled_ppois_bdg=%s"
                  % (has_pooled_fe_bw, has_pooled_ppois_bw,
                     has_pooled_fe_bdg, has_pooled_ppois_bdg))
        _record("6-pooled_no_chrom_sizes_no_pooled_bw", passed)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    print("Starting Stage 22 BigWig Stress Tests\n")

    test_no_chrom_sizes_no_bigwig()
    test_chrom_sizes_configured_bigwig()
    test_chrom_sizes_missing_file()
    test_signal_tracks_false()
    test_pooled_chrom_sizes_pooled_bigwig()
    test_pooled_no_chrom_sizes_no_pooled_bigwig()

    print("\nSummary: %d/%d tests passed." % (PASSED, TOTAL))

    # Clean __pycache__
    pycache = os.path.join(REPO_ROOT, "scripts", "__pycache__")
    if os.path.isdir(pycache):
        shutil.rmtree(pycache)

    if PASSED < TOTAL:
        sys.exit(1)


if __name__ == "__main__":
    main()
