#!/usr/bin/env python3
"""Stage 30 stress tests — strict input validation."""

import os
import shutil
import subprocess
import sys
import tempfile

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
VALIDATOR = os.path.join(REPO_ROOT, "scripts", "validate_samples.py")

PASSED = 0
TOTAL = 0


def _record(name, passed):
    global PASSED, TOTAL
    TOTAL += 1
    if passed:
        PASSED += 1
        print("PASS: %s" % name)
    else:
        print("FAIL: %s" % name)


def _write_files(workdir, config_extra="", samples_lines=None):
    config_path = os.path.join(workdir, "config.yaml")
    samples_path = os.path.join(workdir, "samples.tsv")

    with open(config_path, "w") as f:
        f.write(f'samples: "{samples_path}"\n')
        f.write("use_control: false\n")
        f.write(config_extra)

    if samples_lines is None:
        samples_lines = [
            "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\t"
            "peak_mode\tgenome\tbowtie2_index",
            "S1\t{workdir}/R1.fq\t{workdir}/R2.fq\tPE\tchipseq\tCTCF\t"
            "narrow\ths\t{workdir}/idx",
        ]
    # Replace {workdir} placeholder
    resolved = [line.format(workdir=workdir) for line in samples_lines]
    with open(samples_path, "w") as f:
        f.write("\n".join(resolved) + "\n")

    return config_path


def _run_validator(config_path, strict=False):
    cmd = [sys.executable, VALIDATOR, "--config", config_path]
    if strict:
        cmd.append("--strict-inputs")
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.returncode, result.stdout + result.stderr


# ---------------------------------------------------------------------------
# Test 1: --strict-inputs with missing fastq_1 → fails
# ---------------------------------------------------------------------------

def test_strict_missing_fastq1():
    workdir = tempfile.mkdtemp(prefix="s30_t1_", dir="/tmp")
    try:
        fq1 = os.path.join(workdir, "R1.fq")
        fq2 = os.path.join(workdir, "R2.fq")
        # Create fq2 but not fq1
        with open(fq2, "w") as f:
            f.write("")

        config_path = _write_files(workdir)
        rc, output = _run_validator(config_path, strict=True)

        passed = (rc != 0
                  and "fastq_1" in output.lower()
                  and "not found" in output.lower())
        if not passed:
            print("  rc=%s output=%s" % (rc, output[-200:]))
        _record("1-strict_missing_fastq1", passed)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


# ---------------------------------------------------------------------------
# Test 2: --strict-inputs PE missing fastq_2 → fails
# ---------------------------------------------------------------------------

def test_strict_missing_fastq2():
    workdir = tempfile.mkdtemp(prefix="s30_t2_", dir="/tmp")
    try:
        fq1 = os.path.join(workdir, "R1.fq")
        with open(fq1, "w") as f:
            f.write("")

        config_path = _write_files(workdir)
        rc, output = _run_validator(config_path, strict=True)

        passed = (rc != 0
                  and "fastq_2" in output.lower()
                  and "not found" in output.lower())
        if not passed:
            print("  rc=%s output=%s" % (rc, output[-200:]))
        _record("2-strict_missing_fastq2", passed)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


# ---------------------------------------------------------------------------
# Test 3: --strict-inputs with missing Bowtie2 index → fails
# ---------------------------------------------------------------------------

def test_strict_missing_bt2():
    workdir = tempfile.mkdtemp(prefix="s30_t3_", dir="/tmp")
    try:
        fq1 = os.path.join(workdir, "R1.fq")
        fq2 = os.path.join(workdir, "R2.fq")
        for f in (fq1, fq2):
            with open(f, "w") as fh:
                fh.write("")

        config_path = _write_files(workdir)
        rc, output = _run_validator(config_path, strict=True)

        passed = (rc != 0
                  and "bowtie2 index" in output.lower()
                  and ".1.bt2" in output)
        if not passed:
            print("  rc=%s output=%s" % (rc, output[-300:]))
        _record("3-strict_missing_bt2", passed)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


# ---------------------------------------------------------------------------
# Test 4: --strict-inputs with complete .bt2 set → passes
# ---------------------------------------------------------------------------

def test_strict_bt2_set_passes():
    workdir = tempfile.mkdtemp(prefix="s30_t4_", dir="/tmp")
    try:
        fq1 = os.path.join(workdir, "R1.fq")
        fq2 = os.path.join(workdir, "R2.fq")
        for f in (fq1, fq2):
            with open(f, "w") as fh:
                fh.write("")

        idx_prefix = os.path.join(workdir, "idx")
        for suffix in (".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2",
                       ".rev.1.bt2", ".rev.2.bt2"):
            with open(idx_prefix + suffix, "w") as fh:
                fh.write("")

        config_path = _write_files(workdir)
        rc, output = _run_validator(config_path, strict=True)

        passed = (rc == 0 and "OK" in output)
        if not passed:
            print("  rc=%s output=%s" % (rc, output[-200:]))
        _record("4-strict_bt2_set_passes", passed)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


# ---------------------------------------------------------------------------
# Test 5: --strict-inputs with complete .bt2l set → passes
# ---------------------------------------------------------------------------

def test_strict_bt2l_set_passes():
    workdir = tempfile.mkdtemp(prefix="s30_t5_", dir="/tmp")
    try:
        fq1 = os.path.join(workdir, "R1.fq")
        fq2 = os.path.join(workdir, "R2.fq")
        for f in (fq1, fq2):
            with open(f, "w") as fh:
                fh.write("")

        idx_prefix = os.path.join(workdir, "idx")
        for suffix in (".1.bt2l", ".2.bt2l", ".3.bt2l", ".4.bt2l",
                       ".rev.1.bt2l", ".rev.2.bt2l"):
            with open(idx_prefix + suffix, "w") as fh:
                fh.write("")

        config_path = _write_files(workdir)
        rc, output = _run_validator(config_path, strict=True)

        passed = (rc == 0 and "OK" in output)
        if not passed:
            print("  rc=%s output=%s" % (rc, output[-200:]))
        _record("5-strict_bt2l_set_passes", passed)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


# ---------------------------------------------------------------------------
# Test 6: Non-strict placeholder paths → passes (smoke profile compatibility)
# ---------------------------------------------------------------------------

def test_nonstrict_placeholder_paths():
    workdir = tempfile.mkdtemp(prefix="s30_t6_", dir="/tmp")
    try:
        config_path = _write_files(workdir)
        rc, output = _run_validator(config_path, strict=False)

        passed = rc == 0 and "OK" in output
        if not passed:
            print("  rc=%s output=%s" % (rc, output[-200:]))
        _record("6-nonstrict_placeholder_paths", passed)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


# ---------------------------------------------------------------------------
# Test 7: .fq.gz existing file → passes in strict mode
# ---------------------------------------------------------------------------

def test_strict_fq_gz_passes():
    workdir = tempfile.mkdtemp(prefix="s30_t7_", dir="/tmp")
    try:
        fq1 = os.path.join(workdir, "R1.fq.gz")
        fq2 = os.path.join(workdir, "R2.fq.gz")
        for f in (fq1, fq2):
            with open(f, "w") as fh:
                fh.write("")

        idx_prefix = os.path.join(workdir, "idx")
        for suffix in (".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2",
                       ".rev.1.bt2", ".rev.2.bt2"):
            with open(idx_prefix + suffix, "w") as fh:
                fh.write("")

        config_path = _write_files(
            workdir,
            samples_lines=[
                "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\t"
                "peak_mode\tgenome\tbowtie2_index",
                "S1\t{workdir}/R1.fq.gz\t{workdir}/R2.fq.gz\tPE\t"
                "chipseq\tCTCF\tnarrow\ths\t{workdir}/idx",
            ])
        rc, output = _run_validator(config_path, strict=True)

        passed = (rc == 0 and "OK" in output)
        if not passed:
            print("  rc=%s output=%s" % (rc, output[-200:]))
        _record("7-strict_fq_gz_passes", passed)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


# ---------------------------------------------------------------------------
# Test 8: MACS3 fragment fallback emits warning (code-level smoke)
# ---------------------------------------------------------------------------
# The MACS3 fallback warning is in chipseq.smk. We verify the code
# path exists by checking the file contains the expected warning text.

def test_macs3_fallback_warning():
    chipseq_path = os.path.join(REPO_ROOT, "workflow", "rules", "chipseq.smk")
    with open(chipseq_path) as fh:
        content = fh.read()

    passed = ("MACS3 fragment size unavailable" in content
              and "extendReads 200" in content
              and "sys.stderr" in content)
    _record("8-macs3_fallback_warning", passed)


# ---------------------------------------------------------------------------
def main():
    print("Starting Stage 30 Strict Input Validation Tests\n")

    test_strict_missing_fastq1()
    test_strict_missing_fastq2()
    test_strict_missing_bt2()
    test_strict_bt2_set_passes()
    test_strict_bt2l_set_passes()
    test_nonstrict_placeholder_paths()
    test_strict_fq_gz_passes()
    test_macs3_fallback_warning()

    print("\nSummary: %d/%d tests passed." % (PASSED, TOTAL))

    pycache = os.path.join(REPO_ROOT, "scripts", "__pycache__")
    if os.path.isdir(pycache):
        shutil.rmtree(pycache)

    if PASSED < TOTAL:
        sys.exit(1)


if __name__ == "__main__":
    main()
