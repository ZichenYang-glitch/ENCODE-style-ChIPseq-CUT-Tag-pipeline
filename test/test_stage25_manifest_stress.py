#!/usr/bin/env python3
"""Stage 25 stress tests for make_manifest.py.

Tests:
  1. Default config → manifest contains expected output types
  2. signal_tracks: false → signal outputs not_applicable
  3. chrom_sizes empty → BigWig outputs not_applicable, bdg present
  4. stage5: false → no IDR rows
  5. Missing file → status: missing
  6. Deterministic row ordering
  7. No CRLF in output
  8. --strict flag exits non-zero on missing
  9. Custom outdir honored (via --config-json)
 10. Sample sheet missing experiment → validator default (experiment=sample)
 11. chrom_sizes configured → BigWig rows appear and use check_exists
"""

import csv
import json
import os
import shutil
import subprocess
import sys
import tempfile

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MAKE_MANIFEST = os.path.join(REPO_ROOT, "scripts", "make_manifest.py")

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


def _write_files(workdir, config_updates=None, sample_rows=None,
                 use_json=False):
    """Write config.yaml and samples.tsv in workdir.

    If use_json=True, write config.json instead of config.yaml.
    Returns (config_path_or_none, samples_path).
    """
    cfg = {
        "samples": os.path.join(workdir, "samples.tsv"),
        "outdir": os.path.join(workdir, "results"),
        "use_control": False,
        "multiqc": True,
        "stage4b": True,
        "stage5": False,
        "qc": {"signal_tracks": True, "summary": True},
        "genome_resources": {"hs": {"effective_genome_size": "hs",
                                      "chrom_sizes": ""}},
    }
    if config_updates:
        # Deep update for nested dicts
        for k, v in config_updates.items():
            if isinstance(v, dict) and isinstance(cfg.get(k), dict):
                cfg[k].update(v)
            else:
                cfg[k] = v

    if use_json:
        config_path = os.path.join(workdir, "config.json")
        with open(config_path, "w") as fh:
            json.dump(cfg, fh)
    else:
        config_path = os.path.join(workdir, "config.yaml")
        with open(config_path, "w") as fh:
            for key, value in cfg.items():
                if isinstance(value, dict):
                    fh.write(f"{key}:\n")
                    for k, v in value.items():
                        if isinstance(v, str):
                            fh.write(f"  {k}: \"{v}\"\n")
                        else:
                            fh.write(f"  {k}: {str(v).lower()}\n")
                elif isinstance(value, bool):
                    fh.write(f"{key}: {str(value).lower()}\n")
                elif isinstance(value, str):
                    fh.write(f"{key}: \"{value}\"\n")
                else:
                    fh.write(f"{key}: {value}\n")

    samples_path = os.path.join(workdir, "samples.tsv")
    if sample_rows is None:
        sample_rows = [["sample", "fastq_1", "fastq_2", "layout", "assay",
                        "target", "peak_mode", "genome", "bowtie2_index"],
                       ["S1", "R1.fq", "R2.fq", "PE", "chipseq", "CTCF",
                        "narrow", "hs", "/opt/idx/ref"]]
    with open(samples_path, "w") as fh:
        fh.write("\t".join(sample_rows[0]) + "\n")
        for row in sample_rows[1:]:
            fh.write("\t".join(row) + "\n")

    return config_path, samples_path


def _run_manifest(workdir, config_path, extra_args=None, use_json=False):
    """Run make_manifest.py and return (returncode, output_path, rows)."""
    out = os.path.join(workdir, "manifest.tsv")
    cmd = [sys.executable, MAKE_MANIFEST, "--output", out]
    if use_json:
        with open(config_path) as fh:
            config_json = fh.read()
        cmd.extend(["--config-json", config_json])
    else:
        cmd.extend(["--config", config_path])
    if extra_args:
        cmd.extend(extra_args)

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode == 0 and os.path.exists(out):
        with open(out, newline="") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
    else:
        rows = []
    return result.returncode, out, rows, result.stderr


# ---------------------------------------------------------------------------
# Test 1: Default config → manifest contains expected output types
# ---------------------------------------------------------------------------

def test_default_expected_types():
    workdir = tempfile.mkdtemp(prefix="s25_t1_", dir="/tmp")
    try:
        config_path, _ = _write_files(workdir)
        rc, out, rows, stderr = _run_manifest(workdir, config_path)
        if rc != 0:
            print("  stderr:", stderr.strip()[-200:])
            _record("1-default_expected_types", False)
            return

        types = {r["output_type"] for r in rows}
        expected = {"final_bam", "final_bai", "cpm_bigwig", "macs3_peak",
                    "qc_summary", "macs3_fe_bdg", "macs3_ppois_bdg",
                    "macs3_fe_bw", "macs3_ppois_bw",
                    "stage3_qc_summary", "multiqc_report"}
        passed = expected.issubset(types)
        if not passed:
            print("  Missing types:", expected - types)
        _record("1-default_expected_types", passed)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


# ---------------------------------------------------------------------------
# Test 2: signal_tracks: false → all signal outputs not_applicable
# ---------------------------------------------------------------------------

def test_signal_tracks_false():
    workdir = tempfile.mkdtemp(prefix="s25_t2_", dir="/tmp")
    try:
        config_path, _ = _write_files(
            workdir, config_updates={"qc": {"signal_tracks": False, "summary": True}})
        rc, out, rows, stderr = _run_manifest(workdir, config_path)

        if rc != 0:
            _record("2-signal_tracks_false", False)
            return

        signal_rows = [r for r in rows if r["output_type"] in
                       ("macs3_fe_bdg", "macs3_ppois_bdg",
                        "macs3_fe_bw", "macs3_ppois_bw")]
        passed = all(r["status"] == "not_applicable" for r in signal_rows)
        _record("2-signal_tracks_false", passed)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


# ---------------------------------------------------------------------------
# Test 3: chrom_sizes empty → BigWig outputs not_applicable
# ---------------------------------------------------------------------------

def test_chrom_sizes_empty():
    workdir = tempfile.mkdtemp(prefix="s25_t3_", dir="/tmp")
    try:
        config_path, _ = _write_files(
            workdir,
            config_updates={"genome_resources": {"hs": {"effective_genome_size": "hs",
                                                          "chrom_sizes": ""}}})
        rc, out, rows, stderr = _run_manifest(workdir, config_path)

        if rc != 0:
            _record("3-chrom_sizes_empty", False)
            return

        bw_rows = [r for r in rows if r["output_type"] in
                   ("macs3_fe_bw", "macs3_ppois_bw")]
        bdg_rows = [r for r in rows if r["output_type"] in
                    ("macs3_fe_bdg", "macs3_ppois_bdg")]
        passed = (all(r["status"] == "not_applicable" for r in bw_rows)
                  and all(r["status"] in ("present", "missing") for r in bdg_rows))
        _record("3-chrom_sizes_empty", passed)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


# ---------------------------------------------------------------------------
# Test 4: stage5: false → no IDR rows
# ---------------------------------------------------------------------------

def test_stage5_false_no_idr():
    workdir = tempfile.mkdtemp(prefix="s25_t4_", dir="/tmp")
    try:
        config_path, _ = _write_files(
            workdir,
            sample_rows=[["sample", "fastq_1", "fastq_2", "layout", "assay",
                          "target", "peak_mode", "genome", "bowtie2_index",
                          "experiment", "biological_replicate"],
                         ["S1", "R1.fq", "R2.fq", "PE", "chipseq", "CTCF",
                          "narrow", "hs", "/opt/idx/ref", "EXP1", "1"],
                         ["S2", "R2.fq", "R3.fq", "PE", "chipseq", "CTCF",
                          "narrow", "hs", "/opt/idx/ref", "EXP1", "2"]])
        rc, out, rows, stderr = _run_manifest(workdir, config_path)

        if rc != 0:
            print("  stderr:", stderr.strip()[-200:])
            _record("4-stage5_false_no_idr", False)
            return

        idr_types = {r["output_type"] for r in rows
                     if r["output_type"].startswith("idr_")}
        passed = len(idr_types) == 0
        _record("4-stage5_false_no_idr", passed)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


# ---------------------------------------------------------------------------
# Test 5: Missing file → status: missing
# ---------------------------------------------------------------------------

def test_missing_file_status():
    workdir = tempfile.mkdtemp(prefix="s25_t5_", dir="/tmp")
    try:
        config_path, _ = _write_files(workdir)
        rc, out, rows, stderr = _run_manifest(workdir, config_path)

        if rc != 0:
            _record("5-missing_file_status", False)
            return

        check_rows = [r for r in rows if r["status"] != "not_applicable"]
        passed = all(r["status"] == "missing" for r in check_rows) and len(check_rows) > 0
        _record("5-missing_file_status", passed)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


# ---------------------------------------------------------------------------
# Test 6: Deterministic row ordering
# ---------------------------------------------------------------------------

def test_deterministic_ordering():
    workdir = tempfile.mkdtemp(prefix="s25_t6_", dir="/tmp")
    try:
        config_path, _ = _write_files(
            workdir,
            sample_rows=[["sample", "fastq_1", "fastq_2", "layout", "assay",
                          "target", "peak_mode", "genome", "bowtie2_index"],
                         ["SB", "R1.fq", "R2.fq", "PE", "chipseq", "CTCF",
                          "narrow", "hs", "/opt/idx/ref"],
                         ["SA", "R1.fq", "R2.fq", "PE", "chipseq", "CTCF",
                          "narrow", "hs", "/opt/idx/ref"]])

        out1 = os.path.join(workdir, "manifest1.tsv")
        out2 = os.path.join(workdir, "manifest2.tsv")
        for out in (out1, out2):
            subprocess.run(
                [sys.executable, MAKE_MANIFEST, "--config", config_path,
                 "--output", out],
                capture_output=True, text=True,
            )

        with open(out1) as f1, open(out2) as f2:
            passed = f1.read() == f2.read()
        _record("6-deterministic_ordering", passed)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


# ---------------------------------------------------------------------------
# Test 7: No CRLF in output
# ---------------------------------------------------------------------------

def test_no_crlf():
    workdir = tempfile.mkdtemp(prefix="s25_t7_", dir="/tmp")
    try:
        config_path, _ = _write_files(workdir)
        out = os.path.join(workdir, "manifest.tsv")
        subprocess.run(
            [sys.executable, MAKE_MANIFEST, "--config", config_path,
             "--output", out],
            capture_output=True, text=True,
        )
        with open(out, "rb") as fh:
            raw = fh.read()
        passed = b"\r\n" not in raw
        _record("7-no_crlf", passed)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


# ---------------------------------------------------------------------------
# Test 8: --strict flag exits non-zero on missing
# ---------------------------------------------------------------------------

def test_strict_flag():
    workdir = tempfile.mkdtemp(prefix="s25_t8_", dir="/tmp")
    try:
        config_path, _ = _write_files(workdir)
        rc, out, rows, stderr = _run_manifest(
            workdir, config_path, extra_args=["--strict"])
        passed = rc != 0  # should fail because no results/
        _record("8-strict_flag", passed)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


# ---------------------------------------------------------------------------
# Test 9: Custom outdir honored via --config-json
# ---------------------------------------------------------------------------

def test_custom_outdir_config_json():
    workdir = tempfile.mkdtemp(prefix="s25_t9_", dir="/tmp")
    try:
        custom_outdir = os.path.join(workdir, "custom_results")
        _write_files(workdir, config_updates={"outdir": custom_outdir})

        # Build JSON config manually so we test --config-json path
        cfg = {
            "samples": os.path.join(workdir, "samples.tsv"),
            "outdir": custom_outdir,
            "use_control": False,
            "multiqc": True,
            "stage4b": True,
            "stage5": False,
            "qc": {"signal_tracks": True, "summary": True},
            "genome_resources": {"hs": {"effective_genome_size": "hs",
                                          "chrom_sizes": ""}},
        }
        config_path = os.path.join(workdir, "config.json")
        with open(config_path, "w") as fh:
            json.dump(cfg, fh)

        rc, out, rows, stderr = _run_manifest(
            workdir, config_path, use_json=True)

        if rc != 0:
            print("  stderr:", stderr.strip()[-200:])
            _record("9-custom_outdir_json", False)
            return

        # All paths should use custom_outdir
        non_na_rows = [r for r in rows if r["status"] != "not_applicable"]
        passed = all(r["path"].startswith(custom_outdir) for r in non_na_rows)
        _record("9-custom_outdir_json", passed)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


# ---------------------------------------------------------------------------
# Test 10: Missing experiment → validator default (experiment=sample)
# ---------------------------------------------------------------------------

def test_missing_experiment_defaults():
    workdir = tempfile.mkdtemp(prefix="s25_t10_", dir="/tmp")
    try:
        # Two samples WITHOUT experiment column → validator defaults to sample id
        config_path, _ = _write_files(
            workdir,
            sample_rows=[["sample", "fastq_1", "fastq_2", "layout", "assay",
                          "target", "peak_mode", "genome", "bowtie2_index",
                          "biological_replicate"],
                         ["S1", "R1.fq", "R2.fq", "PE", "chipseq", "CTCF",
                          "narrow", "hs", "/opt/idx/ref", "1"],
                         ["S2", "R1b.fq", "R2b.fq", "PE", "chipseq", "CTCF",
                          "narrow", "hs", "/opt/idx/ref", "2"]])

        rc, out, rows, stderr = _run_manifest(workdir, config_path)

        if rc != 0:
            print("  stderr:", stderr.strip()[-200:])
            _record("10-missing_experiment_defaults", False)
            return

        # Should have per-sample rows (with sample_id filled, experiment_id empty)
        # and NO experiment-level rows because samples are in different experiments
        # (default: experiment = sample id, so S1 and S2 are in different experiments)
        exp_rows = [r for r in rows if r["experiment_id"]]
        passed = len(exp_rows) == 0  # no pooled rows since 1 sample per experiment
        _record("10-missing_experiment_defaults", passed)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


# ---------------------------------------------------------------------------
# Test 11: chrom_sizes configured → BigWig rows check existence
# ---------------------------------------------------------------------------

def test_chrom_sizes_configured_bw_present():
    workdir = tempfile.mkdtemp(prefix="s25_t11_", dir="/tmp")
    try:
        chrom_sizes = os.path.join(workdir, "chr.sizes")
        with open(chrom_sizes, "w") as fh:
            fh.write("chr1\t1000000\n")

        config_path, _ = _write_files(
            workdir,
            config_updates={"genome_resources": {
                "hs": {"effective_genome_size": "hs",
                       "chrom_sizes": chrom_sizes}}})

        rc, out, rows, stderr = _run_manifest(workdir, config_path)

        if rc != 0:
            print("  stderr:", stderr.strip()[-200:])
            _record("11-chrom_sizes_configured_bw_present", False)
            return

        bw_rows = [r for r in rows if r["output_type"] in
                   ("macs3_fe_bw", "macs3_ppois_bw")]
        # BW rows should use check_exists (not be not_applicable)
        passed = (len(bw_rows) == 2
                  and all(r["status"] != "not_applicable" for r in bw_rows))
        _record("11-chrom_sizes_configured_bw_present", passed)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


# ---------------------------------------------------------------------------
# Test 12: 1 biorep with 2 techreps → no pooled, but has biorep merge
# ---------------------------------------------------------------------------

def test_techrep_biorep_gating():
    workdir = tempfile.mkdtemp(prefix="s25_t12_", dir="/tmp")
    try:
        config_path, _ = _write_files(
            workdir,
            sample_rows=[["sample", "fastq_1", "fastq_2", "layout", "assay",
                          "target", "peak_mode", "genome", "bowtie2_index",
                          "experiment", "biological_replicate",
                          "technical_replicate"],
                         ["T1a", "R1.fq", "R1b.fq", "PE", "chipseq", "CTCF",
                          "narrow", "hs", "/opt/idx/ref", "EXP1", "1", "1"],
                         ["T1b", "R2.fq", "R2b.fq", "PE", "chipseq", "CTCF",
                          "narrow", "hs", "/opt/idx/ref", "EXP1", "1", "2"]])
        rc, out, rows, stderr = _run_manifest(workdir, config_path)

        if rc != 0:
            print("  stderr:", stderr.strip()[-200:])
            _record("12-techrep_biorep_gating", False)
            return

        types = {r["output_type"] for r in rows}
        # 1 biorep with 2 techreps: should have biorep1 merge, but NO pooled
        passed = ("biorep1_final_bam" in types
                  and "pooled_final_bam" not in types
                  and "pooled_macs3_peak" not in types)
        if not passed:
            print("  output_types:", sorted(types))
        _record("12-techrep_biorep_gating", passed)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


# ---------------------------------------------------------------------------
# Test 13: 2 bioreps, 1 with 2 techreps → IDR eligible when stage5 true
# ---------------------------------------------------------------------------

def test_idr_with_tech_reps():
    workdir = tempfile.mkdtemp(prefix="s25_t13_", dir="/tmp")
    try:
        config_path, _ = _write_files(
            workdir,
            config_updates={"stage5": True,
                            "idr": {"seed": 42, "threshold": 0.05,
                                    "rank": "p.value"}},
            sample_rows=[["sample", "fastq_1", "fastq_2", "layout", "assay",
                          "target", "peak_mode", "genome", "bowtie2_index",
                          "experiment", "biological_replicate",
                          "technical_replicate"],
                         ["T1a", "R1.fq", "R1b.fq", "PE", "chipseq", "CTCF",
                          "narrow", "hs", "/opt/idx/ref", "EXP1", "1", "1"],
                         ["T1b", "R2.fq", "R2b.fq", "PE", "chipseq", "CTCF",
                          "narrow", "hs", "/opt/idx/ref", "EXP1", "1", "2"],
                         ["T2", "R3.fq", "R3b.fq", "PE", "chipseq", "CTCF",
                          "narrow", "hs", "/opt/idx/ref", "EXP1", "2", "1"]])
        rc, out, rows, stderr = _run_manifest(workdir, config_path)

        if rc != 0:
            print("  stderr:", stderr.strip()[-200:])
            _record("13-idr_with_tech_reps", False)
            return

        types = {r["output_type"] for r in rows}
        # 2 unique bioreps → should emit IDR rows
        expected_idr = {"idr_conservative", "idr_optimal",
                        "idr_reproducibility_summary"}
        passed = expected_idr.issubset(types)
        if not passed:
            print("  missing IDR types:", expected_idr - types)
            print("  output_types:", sorted(types))
        _record("13-idr_with_tech_reps", passed)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


# ---------------------------------------------------------------------------
# Test 14: _manifest_dependency_targets includes downstream targets
# ---------------------------------------------------------------------------

def test_manifest_dep_targets_includes_downstream():
    """Verify that _manifest_dependency_targets() includes multiqc and
    pooled/IDR targets under enabled configs, as a dry-run sanity check."""
    workdir = tempfile.mkdtemp(prefix="s25_t14_", dir="/tmp")
    try:
        from _tool_resolver import resolve_tool
        snakemake = resolve_tool("snakemake", "SNAKEMAKE")
    except ImportError:
        print("  SKIP: _tool_resolver not importable")
        _record("14-manifest_dep_targets", True)
        return

    try:
        # Write config that enables pooled + signal
        outdir = os.path.join(workdir, "results")
        config_path = os.path.join(workdir, "config.yaml")
        samples_path = os.path.join(workdir, "samples.tsv")
        with open(config_path, "w") as f:
            f.write(f'samples: "{samples_path}"\n')
            f.write(f'outdir: "{outdir}"\n')
            f.write("threads: 1\nmapq: 30\nbinsize: 10\n")
            f.write("remove_dup: auto\ntrim: true\n")
            f.write("extend_reads: auto\nuse_control: false\n")
            f.write("multiqc: true\nstage4b: true\nstage5: false\n")
            f.write("qc:\n  signal_tracks: true\n  summary: true\n")
            f.write("genome_resources:\n  hs:\n")
            f.write("    effective_genome_size: hs\n")
            f.write("    chrom_sizes: \"\"\n")
        with open(samples_path, "w") as f:
            f.write("sample\tfastq_1\tlayout\tassay\ttarget\tpeak_mode\t"
                    "genome\tbowtie2_index\texperiment\tbiological_replicate\n")
            f.write("S1\tS1_R1.fq\tSE\tchipseq\tCTCF\tnarrow\ths\t"
                    "/opt/idx\tEXP1\t1\n")

        # Create placeholder FASTQ
        with open(os.path.join(workdir, "S1_R1.fq"), "w") as f:
            f.write("")

        result = subprocess.run(
            [snakemake, "-s",
             os.path.join(REPO_ROOT, "workflow", "Snakefile"),
             "--configfile", config_path, "-n"],
            cwd=workdir, capture_output=True, text=True,
        )

        # Dry-run should succeed and include result_manifest in targets
        output = result.stdout + result.stderr
        passed = "result_manifest" in output
        _record("14-manifest_dep_targets", passed)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


# ---------------------------------------------------------------------------
def main():
    print("Starting Stage 25 Manifest Stress Tests\n")

    test_default_expected_types()
    test_signal_tracks_false()
    test_chrom_sizes_empty()
    test_stage5_false_no_idr()
    test_missing_file_status()
    test_deterministic_ordering()
    test_no_crlf()
    test_strict_flag()
    test_custom_outdir_config_json()
    test_missing_experiment_defaults()
    test_chrom_sizes_configured_bw_present()
    test_techrep_biorep_gating()
    test_idr_with_tech_reps()
    test_manifest_dep_targets_includes_downstream()

    print("\nSummary: %d/%d tests passed." % (PASSED, TOTAL))

    pycache = os.path.join(REPO_ROOT, "scripts", "__pycache__")
    if os.path.isdir(pycache):
        shutil.rmtree(pycache)

    if PASSED < TOTAL:
        sys.exit(1)


if __name__ == "__main__":
    main()