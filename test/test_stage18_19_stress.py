#!/usr/bin/env python3
"""Stage 18/19 stress tests — TSS profiles and ATAC baseline support."""

import os
import shutil
import sys
import tempfile

_REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(_REPO, "test"))

from test_stage12_stress import (
    _DISABLED_QC,
    _make_config,
    _run_snakemake_dryrun,
    _run_validate,
    _write_samples_tsv,
    _write_yaml,
)


def _sample(td, sid="s1", assay="chipseq", peak_mode="narrow"):
    """Return a minimal PE treatment sample row with temp input paths."""
    return {
        "sample": sid,
        "fastq_1": os.path.join(td, f"{sid}_R1.fq.gz"),
        "fastq_2": os.path.join(td, f"{sid}_R2.fq.gz"),
        "layout": "PE",
        "assay": assay,
        "target": "ATAC" if assay == "atac" else "CTCF",
        "peak_mode": peak_mode,
        "genome": "hg38",
        "bowtie2_index": os.path.join(td, "bt2", "hg38"),
    }


def _touch_inputs(row):
    """Create placeholder FASTQs and one Bowtie2 index shard."""
    os.makedirs(os.path.dirname(row["bowtie2_index"]), exist_ok=True)
    for path in (
        row["fastq_1"],
        row["fastq_2"],
        row["bowtie2_index"] + ".1.bt2",
    ):
        with open(path, "w") as fh:
            fh.write("")


def _write_gtf(path):
    with open(path, "w") as fh:
        fh.write(
            'chr1\tsrc\ttranscript\t10\t50\t.\t+\t.\t'
            'gene_id "g1"; transcript_id "tx1";\n'
        )


def _base_genome_resources(gtf=""):
    return {
        "hg38": {
            "effective_genome_size": 2913022398,
            "gtf": gtf,
        },
    }


def test_tss_true_missing_gtf_fails_validation():
    td = tempfile.mkdtemp(prefix="s18_validate_")
    try:
        row = _sample(td)
        samples_tsv = os.path.join(td, "samples.tsv")
        config_yml = os.path.join(td, "config.yaml")
        _write_samples_tsv(samples_tsv, [row])
        config = _make_config(
            qc_block={"tss_enrichment": "true"},
            genome_resources=_base_genome_resources(gtf=""),
        )
        config["samples"] = samples_tsv
        _write_yaml(config_yml, config)

        rc, stdout, stderr = _run_validate(config_yml)
        assert rc != 0, "Expected tss_enrichment without GTF to fail"
        assert "gtf" in stderr.lower(), stderr
    finally:
        shutil.rmtree(td, ignore_errors=True)


def test_tss_false_missing_gtf_passes_validation():
    td = tempfile.mkdtemp(prefix="s18_validate_")
    try:
        row = _sample(td)
        samples_tsv = os.path.join(td, "samples.tsv")
        config_yml = os.path.join(td, "config.yaml")
        _write_samples_tsv(samples_tsv, [row])
        config = _make_config(
            qc_block={"tss_enrichment": "false"},
            genome_resources=_base_genome_resources(gtf=""),
        )
        config["samples"] = samples_tsv
        _write_yaml(config_yml, config)

        rc, stdout, stderr = _run_validate(config_yml)
        assert rc == 0, f"Expected pass, got stderr: {stderr}"
        assert "OK:" in stdout
    finally:
        shutil.rmtree(td, ignore_errors=True)


def test_tss_true_schedules_tss_rules():
    td = tempfile.mkdtemp(prefix="s18_dag_")
    try:
        row = _sample(td)
        _touch_inputs(row)
        gtf = os.path.join(td, "annotation.gtf")
        _write_gtf(gtf)
        samples_tsv = os.path.join(td, "samples.tsv")
        config_yml = os.path.join(td, "config.yaml")
        _write_samples_tsv(samples_tsv, [row])
        qc = dict(_DISABLED_QC)
        qc["tss_enrichment"] = "true"
        config = _make_config(
            qc_block=qc,
            genome_resources=_base_genome_resources(gtf=gtf),
        )
        config["samples"] = samples_tsv
        config["outdir"] = os.path.join(td, "results")
        _write_yaml(config_yml, config)

        rc, stdout, stderr = _run_snakemake_dryrun(config_yml, samples_tsv, quiet=False)
        assert rc == 0, f"Dry-run failed: {stderr}"
        output = stdout + stderr
        assert "tss_bed_from_gtf" in output
        assert "tss_enrichment_profile" in output
    finally:
        shutil.rmtree(td, ignore_errors=True)


def test_tss_false_does_not_schedule_tss_rules():
    td = tempfile.mkdtemp(prefix="s18_dag_")
    try:
        row = _sample(td)
        _touch_inputs(row)
        samples_tsv = os.path.join(td, "samples.tsv")
        config_yml = os.path.join(td, "config.yaml")
        _write_samples_tsv(samples_tsv, [row])
        qc = dict(_DISABLED_QC)
        qc["tss_enrichment"] = "false"
        config = _make_config(qc_block=qc)
        config["samples"] = samples_tsv
        config["outdir"] = os.path.join(td, "results")
        _write_yaml(config_yml, config)

        rc, stdout, stderr = _run_snakemake_dryrun(config_yml, samples_tsv, quiet=False)
        assert rc == 0, f"Dry-run failed: {stderr}"
        output = stdout + stderr
        assert "tss_bed_from_gtf" not in output
        assert "tss_enrichment_profile" not in output
    finally:
        shutil.rmtree(td, ignore_errors=True)


def test_atac_narrow_validates():
    td = tempfile.mkdtemp(prefix="s19_validate_")
    try:
        row = _sample(td, assay="atac", peak_mode="narrow")
        samples_tsv = os.path.join(td, "samples.tsv")
        config_yml = os.path.join(td, "config.yaml")
        _write_samples_tsv(samples_tsv, [row])
        config = _make_config()
        config["samples"] = samples_tsv
        _write_yaml(config_yml, config)

        rc, stdout, stderr = _run_validate(config_yml)
        assert rc == 0, f"Expected ATAC narrow to pass: {stderr}"
        assert "OK:" in stdout
    finally:
        shutil.rmtree(td, ignore_errors=True)


def test_atac_broad_rejected():
    td = tempfile.mkdtemp(prefix="s19_validate_")
    try:
        row = _sample(td, assay="atac", peak_mode="broad")
        samples_tsv = os.path.join(td, "samples.tsv")
        config_yml = os.path.join(td, "config.yaml")
        _write_samples_tsv(samples_tsv, [row])
        config = _make_config()
        config["samples"] = samples_tsv
        _write_yaml(config_yml, config)

        rc, stdout, stderr = _run_validate(config_yml)
        assert rc != 0, "Expected ATAC broad to fail"
        assert "atac" in stderr.lower()
        assert "narrow" in stderr.lower()
    finally:
        shutil.rmtree(td, ignore_errors=True)


def test_atac_dryrun_uses_shared_pipeline_without_cuttag_sidecars():
    td = tempfile.mkdtemp(prefix="s19_dag_")
    try:
        row = _sample(td, assay="atac", peak_mode="narrow")
        _touch_inputs(row)
        samples_tsv = os.path.join(td, "samples.tsv")
        config_yml = os.path.join(td, "config.yaml")
        _write_samples_tsv(samples_tsv, [row])
        qc = dict(_DISABLED_QC)
        config = _make_config(qc_block=qc)
        config["samples"] = samples_tsv
        config["outdir"] = os.path.join(td, "results")
        _write_yaml(config_yml, config)

        rc, stdout, stderr = _run_snakemake_dryrun(config_yml, samples_tsv, quiet=False)
        assert rc == 0, f"Dry-run failed: {stderr}"
        output = stdout + stderr
        assert "macs3_callpeak" in output
        assert "cuttag_fragment_size" not in output
        assert "seacr_call" not in output
    finally:
        shutil.rmtree(td, ignore_errors=True)


def test_atac_macs3_shell_contains_tn5_policy():
    td = tempfile.mkdtemp(prefix="s19_dag_")
    try:
        row = _sample(td, assay="atac", peak_mode="narrow")
        _touch_inputs(row)
        samples_tsv = os.path.join(td, "samples.tsv")
        config_yml = os.path.join(td, "config.yaml")
        _write_samples_tsv(samples_tsv, [row])
        qc = dict(_DISABLED_QC)
        config = _make_config(qc_block=qc)
        config["samples"] = samples_tsv
        config["outdir"] = os.path.join(td, "results")
        _write_yaml(config_yml, config)

        target = os.path.join(td, "results", "s1", "04_peaks", "s1")
        rc, stdout, stderr = _run_snakemake_dryrun(
            config_yml,
            samples_tsv,
            extra_args=["--printshellcmds", target],
        )
        assert rc == 0, f"Dry-run failed: {stderr}"
        output = stdout + stderr
        assert "--nomodel --shift -100 --extsize 200" in output
    finally:
        shutil.rmtree(td, ignore_errors=True)


def run_all():
    tests = [
        test_tss_true_missing_gtf_fails_validation,
        test_tss_false_missing_gtf_passes_validation,
        test_tss_true_schedules_tss_rules,
        test_tss_false_does_not_schedule_tss_rules,
        test_atac_narrow_validates,
        test_atac_broad_rejected,
        test_atac_dryrun_uses_shared_pipeline_without_cuttag_sidecars,
        test_atac_macs3_shell_contains_tn5_policy,
    ]
    for test in tests:
        test()
    print(f"OK: {len(tests)}/{len(tests)} Stage 18/19 stress tests passed")


if __name__ == "__main__":
    run_all()
