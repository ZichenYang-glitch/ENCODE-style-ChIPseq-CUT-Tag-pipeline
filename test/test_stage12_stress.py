#!/usr/bin/env python3
"""Stage 12 stress tests — QC hardening with new modules.

Tests validation acceptance/rejection and DAG scheduling for:
- cross_correlation (phantompeakqualtools)
- preseq_complexity (preseq lc_extrap)
- picard_metrics (Picard CollectMultipleMetrics)

Uses real validate_samples.py calls and Snakemake dry-runs,
following the style of existing stage stress tests.
"""

import os
import sys
import shutil
import tempfile
import subprocess

# --- Tool resolution --------------------------------------------------------
_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_REPO, "test"))
from _tool_resolver import resolve_tool

SNAKEMAKE = resolve_tool("snakemake", "SNAKEMAKE")

# --- Helpers ----------------------------------------------------------------

def _write_yaml(path, data):
    """Write a minimal YAML dict to *path*."""
    with open(path, "w") as fh:
        for key, value in data.items():
            if isinstance(value, bool):
                value = str(value).lower()
            elif isinstance(value, dict):
                _write_nested_yaml(fh, key, value)
                continue
            fh.write(f"{key}: {value}\n")


def _write_nested_yaml(fh, key, data, indent=0):
    """Recursively write a nested YAML mapping."""
    prefix = "  " * indent
    fh.write(f"{prefix}{key}:\n")
    for k, v in data.items():
        if isinstance(v, dict):
            _write_nested_yaml(fh, k, v, indent + 1)
        elif isinstance(v, bool):
            fh.write(f"{prefix}  {k}: {str(v).lower()}\n")
        elif isinstance(v, str):
            if v:
                fh.write(f"{prefix}  {k}: \"{v}\"\n")
            else:
                fh.write(f"{prefix}  {k}: \"\"\n")
        else:
            fh.write(f"{prefix}  {k}: {v}\n")


def _write_samples_tsv(path, rows):
    """Write a sample sheet TSV with header line."""
    with open(path, "w") as fh:
        fh.write("\t".join(rows[0].keys()) + "\n")
        for row in rows:
            fh.write("\t".join(str(v) for v in row.values()) + "\n")


def _make_config(qc_block=None, genome_resources=None):
    """Return a baseline config dict with optional qc and genome_resources."""
    cfg = {
        "samples": "config/samples.tsv",
        "outdir": "results",
        "threads": 4,
        "mapq": 30,
        "binsize": 10,
        "remove_dup": "auto",
        "trim": "true",
        "extend_reads": "auto",
        "use_control": "false",
        "multiqc": "false",
        "stage4b": "false",
        "stage5": "false",
    }
    if qc_block is not None:
        cfg["qc"] = qc_block
    if genome_resources is not None:
        cfg["genome_resources"] = genome_resources
    return cfg


def _chipseq_sample(sid="s1"):
    """Return a minimal ChIP-seq treatment sample row."""
    return {
        "sample": sid,
        "fastq_1": f"/tmp/fake_{sid}_1.fq.gz",
        "fastq_2": "",
        "layout": "SE",
        "assay": "chipseq",
        "target": "CTCF",
        "peak_mode": "narrow",
        "genome": "hg38",
        "bowtie2_index": "/tmp/fake_hg38",
    }


def _run_validate(config_path):
    """Run validate_samples.py on *config_path*. Returns (rc, stdout, stderr)."""
    script = os.path.join(_REPO, "scripts", "validate_samples.py")
    p = subprocess.run(
        [sys.executable, script, "--config", config_path],
        cwd=_REPO,
        capture_output=True,
        text=True,
    )
    return p.returncode, p.stdout, p.stderr


def _run_snakemake_dryrun(config_path, samples_path, extra_args=None, quiet=True):
    """Run snakemake --dry-run. Returns (rc, stdout, stderr).

    When *quiet* is True, uses --quiet (shows warnings but no job list).
    When False, shows full job list for DAG inspection.
    When *extra_args* is provided, --quiet is forced off to avoid
    Snakemake CLI argument collision (--quiet expects mode not target).
    """
    if extra_args:
        quiet = False
    cmd = [
        SNAKEMAKE,
        "-s", os.path.join(_REPO, "workflow", "Snakefile"),
        "--configfile", config_path,
        "--dry-run",
        "--use-conda",
    ]
    if quiet:
        cmd.append("--quiet")
    if extra_args:
        cmd.extend(extra_args)
    p = subprocess.run(
        cmd,
        cwd=_REPO,
        capture_output=True,
        text=True,
    )
    return p.returncode, p.stdout, p.stderr


# --- Validation tests -------------------------------------------------------

def test_missing_qc_block_defaults_new_keys_false():
    """Missing qc block → new keys default false, validation passes."""
    td = tempfile.mkdtemp(prefix="s12_test1_")
    try:
        samples_tsv = os.path.join(td, "samples.tsv")
        _write_samples_tsv(samples_tsv, [_chipseq_sample()])
        config_yml = os.path.join(td, "config.yaml")
        config = _make_config()  # no qc block
        config["samples"] = samples_tsv
        _write_yaml(config_yml, config)

        rc, stdout, stderr = _run_validate(config_yml)
        assert rc == 0, f"Expected pass, got rc={rc}\nstderr: {stderr}"
        assert "OK:" in stdout, f"Expected OK in stdout, got: {stdout}"
    finally:
        shutil.rmtree(td, ignore_errors=True)


def test_all_new_keys_true_validates():
    """All 3 new qc keys set to true → validation passes."""
    td = tempfile.mkdtemp(prefix="s12_test2_")
    try:
        # Create a real reference FASTA file so path validation passes
        fake_ref = os.path.join(td, "fake_ref.fa")
        with open(fake_ref, "w") as fh:
            fh.write(">chr1\nACGT\n")

        samples_tsv = os.path.join(td, "samples.tsv")
        _write_samples_tsv(samples_tsv, [_chipseq_sample()])
        config_yml = os.path.join(td, "config.yaml")
        config = _make_config(qc_block={
            "cross_correlation": "true",
            "preseq_complexity": "true",
            "picard_metrics": "true",
        })
        config["samples"] = samples_tsv
        config["genome_resources"] = {
            "hg38": {
                "effective_genome_size": 2913022398,
                "reference_fasta": fake_ref,
            },
        }
        _write_yaml(config_yml, config)

        rc, stdout, stderr = _run_validate(config_yml)
        assert rc == 0, f"Expected pass, got rc={rc}\nstderr: {stderr}"
        assert "OK:" in stdout
    finally:
        shutil.rmtree(td, ignore_errors=True)


def test_all_new_keys_false_validates():
    """All 3 new qc keys set to false → validation passes."""
    td = tempfile.mkdtemp(prefix="s12_test3_")
    try:
        samples_tsv = os.path.join(td, "samples.tsv")
        _write_samples_tsv(samples_tsv, [_chipseq_sample()])
        config_yml = os.path.join(td, "config.yaml")
        config = _make_config(qc_block={
            "cross_correlation": "false",
            "preseq_complexity": "false",
            "picard_metrics": "false",
        })
        config["samples"] = samples_tsv
        _write_yaml(config_yml, config)

        rc, stdout, stderr = _run_validate(config_yml)
        assert rc == 0, f"Expected pass, got rc={rc}\nstderr: {stderr}"
        assert "OK:" in stdout
    finally:
        shutil.rmtree(td, ignore_errors=True)


def test_invalid_boolean_cross_correlation_rejected():
    """Invalid boolean value for cross_correlation → rejected."""
    td = tempfile.mkdtemp(prefix="s12_test4_")
    try:
        samples_tsv = os.path.join(td, "samples.tsv")
        _write_samples_tsv(samples_tsv, [_chipseq_sample()])
        config_yml = os.path.join(td, "config.yaml")
        config = _make_config(qc_block={"cross_correlation": "maybe"})
        config["samples"] = samples_tsv
        _write_yaml(config_yml, config)

        rc, stdout, stderr = _run_validate(config_yml)
        assert rc != 0, f"Expected failure, got rc={rc}"
        assert "cross_correlation" in stderr, (
            f"Expected 'cross_correlation' in error, got: {stderr}"
        )
    finally:
        shutil.rmtree(td, ignore_errors=True)


def test_picard_metrics_true_missing_ref_fails():
    """picard_metrics: true + empty reference_fasta → validation error."""
    td = tempfile.mkdtemp(prefix="s12_test5_")
    try:
        samples_tsv = os.path.join(td, "samples.tsv")
        _write_samples_tsv(samples_tsv, [_chipseq_sample()])
        config_yml = os.path.join(td, "config.yaml")
        config = _make_config(
            qc_block={"picard_metrics": "true"},
            genome_resources={
                "hg38": {
                    "effective_genome_size": 2913022398,
                    "reference_fasta": "",
                },
            },
        )
        config["samples"] = samples_tsv
        _write_yaml(config_yml, config)

        rc, stdout, stderr = _run_validate(config_yml)
        assert rc != 0, f"Expected failure, got rc={rc}"
        assert "reference_fasta" in stderr.lower(), (
            f"Expected 'reference_fasta' in error, got: {stderr}"
        )
    finally:
        shutil.rmtree(td, ignore_errors=True)


def test_picard_metrics_false_missing_ref_passes():
    """picard_metrics: false + empty reference_fasta → validation passes."""
    td = tempfile.mkdtemp(prefix="s12_test6_")
    try:
        samples_tsv = os.path.join(td, "samples.tsv")
        _write_samples_tsv(samples_tsv, [_chipseq_sample()])
        config_yml = os.path.join(td, "config.yaml")
        config = _make_config(
            qc_block={"picard_metrics": "false"},
            genome_resources={
                "hg38": {
                    "effective_genome_size": 2913022398,
                    "reference_fasta": "",
                },
            },
        )
        config["samples"] = samples_tsv
        _write_yaml(config_yml, config)

        rc, stdout, stderr = _run_validate(config_yml)
        assert rc == 0, f"Expected pass, got rc={rc}\nstderr: {stderr}"
        assert "OK:" in stdout
    finally:
        shutil.rmtree(td, ignore_errors=True)


# --- DAG scheduling tests ---------------------------------------------------

_DISABLED_QC = {
    "blacklist_filter": "false",
    "frip": "false",
    "library_complexity": "false",
    "nrf_pbc": "false",
    "signal_tracks": "false",
    "summary": "false",
    "cuttag_fragment_size": "false",
}


def _create_dryrun_env():
    """Create a temp dir with samples.tsv, config.yaml, and placeholder files.
    Returns (td, config_yml, samples_tsv). Caller must cleanup with rmtree.
    """
    td = tempfile.mkdtemp(prefix="s12_dag_")
    samples_tsv = os.path.join(td, "samples.tsv")
    s1 = _chipseq_sample("s1")
    _write_samples_tsv(samples_tsv, [s1])

    # Create placeholder file paths so validate_samples passes
    for fpath in [s1["fastq_1"], s1["bowtie2_index"] + ".1.bt2"]:
        d = os.path.dirname(fpath)
        if d and not os.path.isdir(d):
            os.makedirs(d, exist_ok=True)
        if not os.path.exists(fpath):
            with open(fpath, "w") as fh:
                fh.write("")

    return td


def test_cross_correlation_true_in_dag():
    """cross_correlation: true → cross_correlation rule in DAG job list."""
    td = _create_dryrun_env()
    try:
        samples_tsv = os.path.join(td, "samples.tsv")
        config_yml = os.path.join(td, "config.yaml")
        qc = dict(_DISABLED_QC)
        qc["cross_correlation"] = "true"
        config = _make_config(qc_block=qc)
        config["samples"] = samples_tsv
        config["outdir"] = os.path.join(td, "results")
        _write_yaml(config_yml, config)

        rc, stdout, stderr = _run_snakemake_dryrun(
            config_yml, samples_tsv, quiet=False,
        )
        assert rc == 0, f"Dry-run failed: rc={rc}\nstderr: {stderr}"

        output = stdout + stderr
        assert "cross_correlation" in output, (
            f"Expected cross_correlation rule in dry-run jobs list"
        )
    finally:
        shutil.rmtree(td, ignore_errors=True)


def test_preseq_complexity_true_in_dag():
    """preseq_complexity: true → preseq_complexity rule in DAG job list."""
    td = _create_dryrun_env()
    try:
        samples_tsv = os.path.join(td, "samples.tsv")
        config_yml = os.path.join(td, "config.yaml")
        qc = dict(_DISABLED_QC)
        qc["preseq_complexity"] = "true"
        config = _make_config(qc_block=qc)
        config["samples"] = samples_tsv
        config["outdir"] = os.path.join(td, "results")
        _write_yaml(config_yml, config)

        rc, stdout, stderr = _run_snakemake_dryrun(
            config_yml, samples_tsv, quiet=False,
        )
        assert rc == 0, f"Dry-run failed: rc={rc}\nstderr: {stderr}"

        output = stdout + stderr
        assert "preseq_complexity" in output, (
            f"Expected preseq_complexity rule in dry-run jobs list"
        )
    finally:
        shutil.rmtree(td, ignore_errors=True)


def test_picard_metrics_true_in_dag():
    """picard_metrics: true + reference_fasta → picard rule in DAG job list."""
    td = _create_dryrun_env()
    try:
        samples_tsv = os.path.join(td, "samples.tsv")
        config_yml = os.path.join(td, "config.yaml")
        fake_ref = os.path.join(td, "fake_ref.fa")
        with open(fake_ref, "w") as fh:
            fh.write(">chr1\nACGT\n")
        qc = dict(_DISABLED_QC)
        qc["picard_metrics"] = "true"
        config = _make_config(
            qc_block=qc,
            genome_resources={
                "hg38": {
                    "effective_genome_size": 2913022398,
                    "reference_fasta": fake_ref,
                },
            },
        )
        config["samples"] = samples_tsv
        config["outdir"] = os.path.join(td, "results")
        _write_yaml(config_yml, config)

        rc, stdout, stderr = _run_snakemake_dryrun(
            config_yml, samples_tsv, quiet=False,
        )
        assert rc == 0, f"Dry-run failed: rc={rc}\nstderr: {stderr}"

        output = stdout + stderr
        assert "picard_collect_multiple_metrics" in output, (
            f"Expected picard_collect_multiple_metrics rule in dry-run jobs list"
        )
    finally:
        shutil.rmtree(td, ignore_errors=True)


def test_all_three_enabled_all_in_dag():
    """All 3 modules enabled → all rules in DAG, dry-run succeeds."""
    td = _create_dryrun_env()
    try:
        samples_tsv = os.path.join(td, "samples.tsv")
        config_yml = os.path.join(td, "config.yaml")
        fake_ref = os.path.join(td, "fake_ref.fa")
        with open(fake_ref, "w") as fh:
            fh.write(">chr1\nACGT\n")
        qc = dict(_DISABLED_QC)
        qc.update({
            "cross_correlation": "true",
            "preseq_complexity": "true",
            "picard_metrics": "true",
        })
        config = _make_config(
            qc_block=qc,
            genome_resources={
                "hg38": {
                    "effective_genome_size": 2913022398,
                    "reference_fasta": fake_ref,
                },
            },
        )
        config["samples"] = samples_tsv
        config["outdir"] = os.path.join(td, "results")
        _write_yaml(config_yml, config)

        rc, stdout, stderr = _run_snakemake_dryrun(
            config_yml, samples_tsv, quiet=False,
        )
        assert rc == 0, f"Dry-run failed: rc={rc}\nstderr: {stderr}"

        output = stdout + stderr
        assert "cross_correlation" in output
        assert "preseq_complexity" in output
        assert "picard_collect_multiple_metrics" in output
    finally:
        shutil.rmtree(td, ignore_errors=True)


def test_all_disabled_dag_unchanged():
    """All new modules disabled → no stage12 rules in DAG."""
    td = _create_dryrun_env()
    try:
        samples_tsv = os.path.join(td, "samples.tsv")
        config_yml = os.path.join(td, "config.yaml")
        qc = dict(_DISABLED_QC)
        qc.update({
            "cross_correlation": "false",
            "preseq_complexity": "false",
            "picard_metrics": "false",
        })
        config = _make_config(qc_block=qc)
        config["samples"] = samples_tsv
        config["outdir"] = os.path.join(td, "results")
        _write_yaml(config_yml, config)

        rc, stdout, stderr = _run_snakemake_dryrun(
            config_yml, samples_tsv, quiet=False,
        )
        assert rc == 0, f"Dry-run failed: rc={rc}\nstderr: {stderr}"

        output = stdout + stderr
        assert "cross_correlation" not in output, (
            f"Should not see cross_correlation in DAG when disabled"
        )
        assert "preseq_complexity" not in output, (
            f"Should not see preseq_complexity in DAG when disabled"
        )
        assert "picard_collect_multiple_metrics" not in output, (
            f"Should not see picard_collect_multiple_metrics in DAG when disabled"
        )
    finally:
        shutil.rmtree(td, ignore_errors=True)


def test_explicit_cross_correlation_target_resolves():
    """Explicit concrete target for cross_correlation resolves successfully."""
    td = _create_dryrun_env()
    try:
        samples_tsv = os.path.join(td, "samples.tsv")
        config_yml = os.path.join(td, "config.yaml")
        qc = dict(_DISABLED_QC)
        qc["cross_correlation"] = "true"
        config = _make_config(qc_block=qc)
        config["outdir"] = os.path.join(td, "results")
        config["samples"] = samples_tsv
        _write_yaml(config_yml, config)

        # Explicit concrete output path (not wildcard rule name)
        target = os.path.join(
            td, "results", "s1", "05_qc", "cross_correlation", "s1.cc.qc",
        )
        rc, stdout, stderr = _run_snakemake_dryrun(
            config_yml, samples_tsv,
            extra_args=[target],
        )
        assert rc == 0, (
            f"Explicit target failed: rc={rc}\nstderr: {stderr}"
        )
    finally:
        shutil.rmtree(td, ignore_errors=True)


# --- Stage 14: cross_correlation_summary DAG tests -----------------------

def test_cross_corr_summary_in_dag():
    """cross_correlation: true → cross_correlation_summary.tsv in DAG."""
    td = _create_dryrun_env()
    try:
        samples_tsv = os.path.join(td, "samples.tsv")
        config_yml = os.path.join(td, "config.yaml")
        qc = dict(_DISABLED_QC)
        qc["cross_correlation"] = "true"
        config = _make_config(qc_block=qc)
        config["samples"] = samples_tsv
        config["outdir"] = os.path.join(td, "results")
        _write_yaml(config_yml, config)

        rc, stdout, stderr = _run_snakemake_dryrun(
            config_yml, samples_tsv, quiet=False,
        )
        assert rc == 0, f"Dry-run failed: rc={rc}\nstderr: {stderr}"

        output = stdout + stderr
        assert "cross_correlation_summary" in output, (
            "Expected cross_correlation_summary rule in DAG when enabled"
        )
    finally:
        shutil.rmtree(td, ignore_errors=True)


def test_cross_corr_summary_not_in_dag_when_disabled():
    """cross_correlation: false → cross_correlation_summary NOT in DAG."""
    td = _create_dryrun_env()
    try:
        samples_tsv = os.path.join(td, "samples.tsv")
        config_yml = os.path.join(td, "config.yaml")
        qc = dict(_DISABLED_QC)
        qc["cross_correlation"] = "false"
        config = _make_config(qc_block=qc)
        config["samples"] = samples_tsv
        config["outdir"] = os.path.join(td, "results")
        _write_yaml(config_yml, config)

        rc, stdout, stderr = _run_snakemake_dryrun(
            config_yml, samples_tsv, quiet=False,
        )
        assert rc == 0, f"Dry-run failed: rc={rc}\nstderr: {stderr}"

        output = stdout + stderr
        assert "cross_correlation_summary" not in output, (
            "Should not see cross_correlation_summary in DAG when disabled"
        )
    finally:
        shutil.rmtree(td, ignore_errors=True)


def test_cross_corr_summary_explicit_target():
    """Explicit target cross_correlation_summary.tsv resolves."""
    td = _create_dryrun_env()
    try:
        samples_tsv = os.path.join(td, "samples.tsv")
        config_yml = os.path.join(td, "config.yaml")
        qc = dict(_DISABLED_QC)
        qc["cross_correlation"] = "true"
        config = _make_config(qc_block=qc)
        config["outdir"] = os.path.join(td, "results")
        config["samples"] = samples_tsv
        _write_yaml(config_yml, config)

        target = os.path.join(td, "results", "multiqc",
                              "cross_correlation_summary.tsv")
        rc, stdout, stderr = _run_snakemake_dryrun(
            config_yml, samples_tsv,
            extra_args=[target],
        )
        assert rc == 0, (
            f"Explicit target failed: rc={rc}\nstderr: {stderr}"
        )
    finally:
        shutil.rmtree(td, ignore_errors=True)


# --- Stage 15: MultiQC custom content tests ------------------------------

def test_multiqc_config_declares_cross_corr_custom_content():
    """workflow/multiqc_config.yaml declares cross-correlation custom data."""
    config_path = os.path.join(_REPO, "workflow", "multiqc_config.yaml")
    with open(config_path) as fh:
        text = fh.read()

    assert "custom_data:" in text
    assert "cross_correlation_qc:" in text
    assert "mnase_qc:" in text
    assert "*mnase_qc_summary.tsv" in text
    assert "plot_type: \"table\"" in text
    assert "cross_correlation_summary.tsv" in text


def test_multiqc_rule_uses_custom_config_and_summary_search_path():
    """multiqc rule passes the custom config and summary TSV to MultiQC."""
    rule_path = os.path.join(_REPO, "workflow", "rules", "report.smk")
    with open(rule_path) as fh:
        text = fh.read()

    assert "multiqc_config.yaml" in text
    assert "--config {params.multiqc_config:q}" in text
    assert "--filename multiqc_report.html" in text
    assert "--force" in text
    assert "params.search_paths:q" in text
    assert "cross_correlation_summary.tsv" in text


def test_multiqc_report_target_with_cross_corr_resolves():
    """MultiQC report target resolves with cross-correlation custom content."""
    td = _create_dryrun_env()
    try:
        samples_tsv = os.path.join(td, "samples.tsv")
        config_yml = os.path.join(td, "config.yaml")
        qc = dict(_DISABLED_QC)
        qc["cross_correlation"] = "true"
        config = _make_config(qc_block=qc)
        config["multiqc"] = "true"
        config["outdir"] = os.path.join(td, "results")
        config["samples"] = samples_tsv
        _write_yaml(config_yml, config)

        target = os.path.join(td, "results", "multiqc", "multiqc_report.html")
        rc, stdout, stderr = _run_snakemake_dryrun(
            config_yml, samples_tsv,
            extra_args=[target],
        )
        assert rc == 0, (
            f"MultiQC report target failed: rc={rc}\nstderr: {stderr}"
        )
        output = stdout + stderr
        assert "cross_correlation_summary" in output
        assert "multiqc" in output
    finally:
        shutil.rmtree(td, ignore_errors=True)


# --- Runner ---------------------------------------------------------------

if __name__ == "__main__":
    import traceback

    tests = [
        # Validation tests
        ("missing_qc_defaults_false", test_missing_qc_block_defaults_new_keys_false),
        ("all_new_true_validates", test_all_new_keys_true_validates),
        ("all_new_false_validates", test_all_new_keys_false_validates),
        ("invalid_bool_rejected", test_invalid_boolean_cross_correlation_rejected),
        ("picard_true_missing_ref_fails", test_picard_metrics_true_missing_ref_fails),
        ("picard_false_missing_ref_passes", test_picard_metrics_false_missing_ref_passes),
        # DAG tests
        ("cross_corr_true_in_dag", test_cross_correlation_true_in_dag),
        ("preseq_true_in_dag", test_preseq_complexity_true_in_dag),
        ("picard_true_in_dag", test_picard_metrics_true_in_dag),
        ("all_three_in_dag", test_all_three_enabled_all_in_dag),
        ("all_disabled_no_dag", test_all_disabled_dag_unchanged),
        ("explicit_target_resolves", test_explicit_cross_correlation_target_resolves),
        # Stage 14: cross_correlation_summary DAG tests
        ("cross_corr_summary_in_dag", test_cross_corr_summary_in_dag),
        ("cross_corr_summary_disabled_not_in_dag", test_cross_corr_summary_not_in_dag_when_disabled),
        ("cross_corr_summary_explicit_target", test_cross_corr_summary_explicit_target),
        # Stage 15: MultiQC custom content tests
        ("multiqc_config_cross_corr_custom_content", test_multiqc_config_declares_cross_corr_custom_content),
        ("multiqc_rule_uses_custom_config", test_multiqc_rule_uses_custom_config_and_summary_search_path),
        ("multiqc_report_target_with_cross_corr", test_multiqc_report_target_with_cross_corr_resolves),
    ]

    passed = 0
    failed = 0
    for name, fn in tests:
        try:
            fn()
            print(f"PASS: {name}")
            passed += 1
        except AssertionError as exc:
            print(f"FAIL: {name} — {exc}")
            failed += 1
        except Exception:
            print(f"FAIL: {name} — unexpected error:")
            traceback.print_exc()
            failed += 1

    print(f"\n{passed} passed, {failed} failed, {passed + failed} total")
    sys.exit(0 if failed == 0 else 1)
