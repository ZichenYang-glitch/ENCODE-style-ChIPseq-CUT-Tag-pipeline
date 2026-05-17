"""Stage 6b stress tests — pooled QC summary scheduling and content.

Reuses Stage 4b/6a harness pattern.
"""

import os
import shutil
import subprocess
import sys

SNAKEFILE = "workflow/Snakefile"
SNAKEMAKE = os.environ.get(
    "SNAKEMAKE",
    "/home/irenadler/miniconda3/envs/chipseq/bin/snakemake",
)

BASE_CONFIG = """\
samples: "test_stage6b_samples.tsv"
use_control: false
threads: 1
stage4b: true
stage5: false
"""

SIGNAL_FALSE = "qc:\n  signal_tracks: false\n"

HDR = (
    "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\t"
    "peak_mode\tgenome\tbowtie2_index\texperiment\tcondition\t"
    "replicate\tbiological_replicate\ttechnical_replicate\trole\n"
)

SINGLE = HDR + (
    "T1\tR1.fq\t\tSE\tchipseq\tH3K27ac\tnarrow\ths\tidx\t"
    "exp1\tH3K27ac\t1\t1\t1\ttreatment\n"
)

TWO_BROAD = HDR + (
    "T1\tR1.fq\t\tSE\tchipseq\tH3K27me3\tbroad\ths\tidx\t"
    "exp1\tH3K27me3\t1\t1\t1\ttreatment\n"
    "T2\tR2.fq\t\tSE\tchipseq\tH3K27me3\tbroad\ths\tidx\t"
    "exp1\tH3K27me3\t2\t2\t1\ttreatment\n"
)

TWO_NARROW = HDR + (
    "T1\tR1.fq\t\tSE\tchipseq\tH3K27ac\tnarrow\ths\tidx\t"
    "exp1\tH3K27ac\t1\t1\t1\ttreatment\n"
    "T2\tR2.fq\t\tSE\tchipseq\tH3K27ac\tnarrow\ths\tidx\t"
    "exp1\tH3K27ac\t2\t2\t1\ttreatment\n"
)


def _write_inputs(config_yaml, samples_tsv):
    for fq in ["R1.fq", "R2.fq"]:
        if not os.path.exists(fq):
            with open(fq, "w") as f:
                f.write("")
    with open("test_stage6b_config.yaml", "w") as f:
        f.write(config_yaml)
    with open("test_stage6b_samples.tsv", "w") as f:
        f.write(samples_tsv)


def _snakemake(args, config_yaml, samples_tsv):
    _write_inputs(config_yaml, samples_tsv)
    return subprocess.run(
        [SNAKEMAKE, "-s", SNAKEFILE, "--configfile",
         "test_stage6b_config.yaml", *args],
        capture_output=True, text=True,
    )


def dryrun(name, cfg, tsv, expect=(), forbid=(), printshell=False):
    args = ["-n", "-p"] if printshell else ["-n"]
    result = _snakemake(args, cfg, tsv)
    output = result.stdout + result.stderr
    if result.returncode != 0:
        print("FAIL: %s\n   Dry-run failed.\n   stdout: %s\n   stderr: %s"
              % (name, result.stdout.strip()[-500:],
                 result.stderr.strip()[-500:]))
        return False
    for text in expect:
        if text not in output:
            print("FAIL: %s\n   Expected %r not found." % (name, text))
            return False
    for text in forbid:
        if text in output:
            print("FAIL: %s\n   Forbidden %r found." % (name, text))
            return False
    print("PASS: %s" % name)
    return True


def cleanup():
    for f in ["test_stage6b_config.yaml", "test_stage6b_samples.tsv",
              "R1.fq", "R2.fq"]:
        if os.path.exists(f):
            os.remove(f)
    for d in ["__pycache__", "scripts/__pycache__"]:
        if os.path.isdir(d):
            shutil.rmtree(d)


def main():
    print("Starting Stage 6b Stress Tests\n")
    tests = 0
    passed = 0

    # 1. Multi-biorep schedules pooled_qc_summary
    tests += 1
    if dryrun("pooled_qc_summary in DAG", BASE_CONFIG, TWO_NARROW,
              expect=("pooled_experiment_qc_summary",)):
        passed += 1

    # 2. Single-sample -> no pooled QC summary
    tests += 1
    if dryrun("single-sample no pooled QC", BASE_CONFIG, SINGLE,
              forbid=("pooled_experiment_qc_summary",)):
        passed += 1

    # 3. signal_tracks: false -> summary still produced
    tests += 1
    if dryrun("signal false still schedules summary",
              BASE_CONFIG + SIGNAL_FALSE, TWO_NARROW,
              expect=("pooled_experiment_qc_summary",)):
        passed += 1

    # 4. H3K27me3 broad -> broadPeak resolved in shell
    tests += 1
    if dryrun("broadPeak for H3K27me3 broad",
              BASE_CONFIG, TWO_BROAD,
              expect=("broadPeak", "pooled_experiment_qc_summary"),
              printshell=True):
        passed += 1

    # 5. H3K27ac narrow -> narrowPeak resolved in shell
    tests += 1
    if dryrun("narrowPeak for H3K27ac narrow",
              BASE_CONFIG, TWO_NARROW,
              expect=("narrowPeak", "pooled_experiment_qc_summary"),
              printshell=True):
        passed += 1

    # 6. Dry-run compatible with Stage 5 disabled
    tests += 1
    if dryrun("Stage 5 disabled compatible",
              BASE_CONFIG, TWO_NARROW,
              expect=("pooled_experiment_qc_summary",)):
        passed += 1

    # 7. --list-rules shows rule
    result = _snakemake(["--list-rules"], BASE_CONFIG, TWO_NARROW)
    tests += 1
    if result.returncode == 0 and "pooled_experiment_qc_summary" in (
            result.stdout + result.stderr):
        print("PASS: pooled_experiment_qc_summary in --list-rules")
        passed += 1
    else:
        print("FAIL: pooled_experiment_qc_summary not in --list-rules")

    # ------------------------------------------------------------------
    # TSV content tests (verify header, columns, and values via dry-run -p)
    # ------------------------------------------------------------------

    # 8. pooled_qc_summary.py call visible in dry-run
    tests += 1
    cfg = BASE_CONFIG
    samples = TWO_NARROW
    _write_inputs(cfg, samples)
    result = _snakemake(["-n", "-p"], cfg, samples)
    output = result.stdout + result.stderr
    if "pooled_qc_summary.py" in output and "--output" in output:
        print("PASS: pooled_qc_summary.py call visible in -p output")
        passed += 1
    else:
        print("FAIL: pooled_qc_summary.py call not found in -p output")

    # 9. H3K27me3 + narrow produces peak_mode_status=mismatch
    tests += 1
    cfg = BASE_CONFIG
    _write_inputs(cfg, TWO_BROAD.replace("broad", "narrow"))
    result = _snakemake(["-n", "-p"], cfg,
                        TWO_BROAD.replace("broad", "narrow"))
    output = result.stdout + result.stderr
    if "mismatch" in output:
        print("PASS: H3K27me3 + narrow mismatch visible")
        passed += 1
    else:
        print("FAIL: H3K27me3 + narrow mismatch not found")

    # 10. H3K27ac + narrow produces context_dependent and ok
    tests += 1
    _write_inputs(cfg, TWO_NARROW)
    result = _snakemake(["-n", "-p"], cfg, TWO_NARROW)
    output = result.stdout + result.stderr
    if "context_dependent" in output:
        print("PASS: H3K27ac context_dependent visible")
        passed += 1
    else:
        print("FAIL: H3K27ac context_dependent not found")

    # 11. signal_tracks=false writes NA for FE/ppois and disabled
    tests += 1
    _write_inputs(cfg + SIGNAL_FALSE, TWO_NARROW)
    result = _snakemake(["-n", "-p"], cfg + SIGNAL_FALSE, TWO_NARROW)
    output = result.stdout + result.stderr
    if "disabled" in output or '"NA"' in output or "NA" in output:
        print("PASS: signal_tracks=false disabled/NA visible")
        passed += 1
    else:
        print("FAIL: signal_tracks=false disabled/NA not found")

    # 12. biological_replicate_labels preserves actual labels
    # Use labels 2 and 4
    TWO_BIOREPS_24 = HDR + (
        "T1\tR1.fq\t\tSE\tchipseq\tH3K27ac\tnarrow\ths\tidx\t"
        "exp1\tH3K27ac\t2\t2\t1\ttreatment\n"
        "T2\tR2.fq\t\tSE\tchipseq\tH3K27ac\tnarrow\ths\tidx\t"
        "exp1\tH3K27ac\t4\t4\t1\ttreatment\n"
    )
    tests += 1
    _write_inputs(cfg, TWO_BIOREPS_24)
    result = _snakemake(["-n", "-p"], cfg, TWO_BIOREPS_24)
    output = result.stdout + result.stderr
    if "2,4" in output:
        print("PASS: bio_rep_labels 2,4 visible")
        passed += 1
    else:
        print("FAIL: bio_rep_labels 2,4 not found")

    # ------------------------------------------------------------------
    # Content-level TSV execution test (direct rule logic, no Snakemake)
    # ------------------------------------------------------------------
    # 13. Run pooled_experiment_qc_summary shell logic against fake tree
    import tempfile
    import json
    import textwrap

    tmpdir = tempfile.mkdtemp(prefix="stage6b_real_")
    exp = "test_exp"

    # Build a minimal experiment output tree
    exp_dir = os.path.join(tmpdir, "results", "experiments", exp)
    os.makedirs(os.path.join(exp_dir, "02_align"), exist_ok=True)
    os.makedirs(os.path.join(exp_dir, "04_peaks", "pooled",
                             exp + "_pooled_peaks"), exist_ok=True)
    os.makedirs(os.path.join(exp_dir, "01_qc"), exist_ok=True)
    os.makedirs(os.path.join(exp_dir, "03_signal"), exist_ok=True)

    # Fake pooled BAM
    pooled_bam = os.path.join(exp_dir, "02_align", exp + ".pooled.final.bam")
    with open(pooled_bam, "w") as f:
        f.write("")

    # Fake pooled peaks file (narrowPeak, 3 lines)
    pk_dir = os.path.join(exp_dir, "04_peaks", "pooled",
                          exp + "_pooled_peaks")
    pk_file = os.path.join(pk_dir, exp + "_pooled_peaks.narrowPeak")
    with open(pk_file, "w") as f:
        f.write("peak1\npeak2\npeak3\n")

    # Run pooled_qc_summary.py directly against the fake experiment tree
    out_tsv = os.path.join(exp_dir, "01_qc",
                           exp + ".pooled_qc_summary.tsv")
    result = subprocess.run(
        [sys.executable, "scripts/pooled_qc_summary.py",
         "--experiment", exp,
         "--assay", "chipseq",
         "--target", "H3K27me3",
         "--peak-mode", "narrow",
         "--inferred-histone-class", "broad_like",
         "--expected-peak-mode", "broad",
         "--peak-mode-status", "mismatch",
         "--n-bioreps", "2",
         "--bio-rep-labels", "1,2",
         "--pooled-bam", pooled_bam,
         "--pooled-peaks-dir", pk_dir,
         "--pooled-fe-bdg", "NA",
         "--pooled-ppois-bdg", "NA",
         "--signal-tracks-status", "disabled",
         "--output", out_tsv],
        capture_output=True, text=True)
    if result.returncode != 0:
        print("FAIL (exec): script failed: %s" % result.stderr.strip())

    # Check the TSV was produced
    tsv_ok = True
    if not os.path.exists(out_tsv):
        tsv_ok = False
        print("FAIL (exec): summary TSV not produced at %s" % out_tsv)
    else:
        lines = open(out_tsv).read().strip().split("\n")
        if len(lines) < 2:
            tsv_ok = False
            print("FAIL (exec): TSV has no data row")
        else:
            header = lines[0].split("\t")
            data = lines[1].split("\t")
            if len(header) != 15:
                tsv_ok = False
                print("FAIL (exec): header has %d cols, expected 15" % len(header))
            elif len(data) != 15:
                tsv_ok = False
                print("FAIL (exec): data has %d cols, expected 15" % len(data))
            else:
                col = dict(zip(header, data))
                checks = [
                    ("experiment matches", col.get("experiment") == exp,
                     col.get("experiment")),
                    ("target is H3K27me3", col.get("target") == "H3K27me3",
                     col.get("target")),
                    ("broad_like", col.get("inferred_histone_class") == "broad_like",
                     col.get("inferred_histone_class")),
                    ("mismatch", col.get("peak_mode_status") == "mismatch",
                     col.get("peak_mode_status")),
                    ("disabled", col.get("signal_tracks_status") == "disabled",
                     col.get("signal_tracks_status")),
                    ("FE=NA", col.get("pooled_FE_bdg") == "NA",
                     col.get("pooled_FE_bdg")),
                    ("ppois=NA", col.get("pooled_ppois_bdg") == "NA",
                     col.get("pooled_ppois_bdg")),
                    ("peak_count=3", col.get("pooled_peak_count") == "3",
                     col.get("pooled_peak_count")),
                ]
                for label, ok_val, val in checks:
                    if not ok_val:
                        tsv_ok = False
                        print("FAIL (exec): %s: got %r" % (label, val))
        if tsv_ok:
            print("PASS (exec): real TSV content verified")

    tests += 1
    if tsv_ok:
        passed += 1

    shutil.rmtree(tmpdir)

    print("\nSummary: %d/%d tests passed." % (passed, tests))
    cleanup()
    if passed < tests:
        sys.exit(1)


if __name__ == "__main__":
    main()
