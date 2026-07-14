"""Snakemake dry-run contracts for replicate merging and experiment pooling."""

from __future__ import annotations

import pytest


SAMPLE_COLUMNS = [
    "sample",
    "fastq_1",
    "fastq_2",
    "layout",
    "assay",
    "target",
    "peak_mode",
    "genome",
    "bowtie2_index",
    "experiment",
    "condition",
    "replicate",
    "biological_replicate",
    "technical_replicate",
    "role",
]


def _replicate_rows(*replicate_pairs):
    """Build treatment rows from ``(biorep, techrep)`` pairs."""

    return [
        {
            "sample": f"T{index}",
            "layout": "SE",
            "assay": "chipseq",
            "target": "H3K27ac",
            "peak_mode": "narrow",
            "experiment": "exp1",
            "condition": "H3K27ac",
            "replicate": str(biorep),
            "biological_replicate": str(biorep),
            "technical_replicate": str(techrep),
            "role": "treatment",
        }
        for index, (biorep, techrep) in enumerate(replicate_pairs, start=1)
    ]


def _config(tmp_path, *, stage4b=True, signal_tracks=True):
    return {
        "outdir": str(tmp_path / "results"),
        "use_control": False,
        "threads": 1,
        "stage4b": stage4b,
        "stage5": False,
        "qc": {"signal_tracks": signal_tracks},
    }


@pytest.fixture
def replicate_dry_run(tmp_path, tmp_config, run_snakemake, monkeypatch):
    """Run one isolated replicate DAG without executing workflow tools."""

    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / ".cache"))
    monkeypatch.setenv("PYTHONDONTWRITEBYTECODE", "1")

    def _run(sample_rows, config):
        rows = ["\t".join(SAMPLE_COLUMNS)]
        for sample in sample_rows:
            fastq = tmp_path / f"{sample['sample']}.fastq"
            fastq.write_text("", encoding="utf-8")
            row = {
                **sample,
                "fastq_1": str(fastq),
                "fastq_2": "",
                "genome": "hs",
                "bowtie2_index": "idx",
            }
            rows.append("\t".join(row[column] for column in SAMPLE_COLUMNS))

        _, config_path, _ = tmp_config(
            config=config,
            samples="\n".join(rows) + "\n",
        )
        result = run_snakemake(config_path, quiet=False)
        output = result.stdout + result.stderr
        assert result.returncode == 0, (
            "Replicate/pooling DAG dry-run failed:\n"
            f"stdout:\n{result.stdout}\n"
            f"stderr:\n{result.stderr}"
        )
        return output

    return _run


def test_disabled_replicate_processing_schedules_no_merge_or_pool(
    replicate_dry_run, tmp_path
):
    output = replicate_dry_run(
        _replicate_rows((1, 1), (2, 1)),
        _config(tmp_path, stage4b=False),
    )

    for rule in (
        "merge_biorep_bam",
        "pool_treatment_bam",
        "macs3_pooled_peaks",
        "pooled_experiment_qc_summary",
    ):
        assert rule not in output


def test_technical_replicates_merge_within_one_biorep_without_pooling(
    replicate_dry_run, tmp_path
):
    output = replicate_dry_run(
        _replicate_rows((1, 1), (1, 2)),
        _config(tmp_path),
    )

    assert "merge_biorep_bam" in output
    assert "experiments/exp1/02_align/biorep1.final.bam" in output
    assert "pool_treatment_bam" not in output
    assert "experiments/exp1/02_align/exp1.pooled.final.bam" not in output
    assert "macs3_pooled_peaks" not in output
    assert "pooled_experiment_qc_summary" not in output


def test_multiple_bioreps_schedule_pooled_bam_peaks_and_qc_summary(
    replicate_dry_run, tmp_path
):
    output = replicate_dry_run(
        _replicate_rows((1, 1), (2, 1)),
        _config(tmp_path),
    )

    assert "merge_biorep_bam" in output
    assert "pool_treatment_bam" in output
    assert "experiments/exp1/02_align/exp1.pooled.final.bam" in output
    assert "macs3_pooled_peaks" in output
    assert "pooled_experiment_qc_summary" in output
    assert "experiments/exp1/01_qc/exp1.pooled_qc_summary.tsv" in output


def test_single_sample_schedules_no_pooled_qc_summary(replicate_dry_run, tmp_path):
    output = replicate_dry_run(
        _replicate_rows((1, 1)),
        _config(tmp_path),
    )

    assert "pool_treatment_bam" not in output
    assert "pooled_experiment_qc_summary" not in output


def test_pooled_qc_summary_does_not_require_signal_tracks(replicate_dry_run, tmp_path):
    output = replicate_dry_run(
        _replicate_rows((1, 1), (2, 1)),
        _config(tmp_path, signal_tracks=False),
    )

    assert "pool_treatment_bam" in output
    assert "pooled_experiment_qc_summary" in output
