"""Snakemake DAG contracts for pooled experiment signal tracks."""

from __future__ import annotations

import pytest


SAMPLE_COLUMNS = (
    "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\tpeak_mode\tgenome\t"
    "bowtie2_index\texperiment\tbiological_replicate\n"
)


def _samples(tmp_path, *, peak_mode, count):
    rows = [SAMPLE_COLUMNS.rstrip("\n")]
    for index in range(1, count + 1):
        read = tmp_path / f"reads_{index}.fastq"
        read.write_text("", encoding="utf-8")
        target = "H3K27me3" if peak_mode == "broad" else "H3K27ac"
        rows.append(
            f"T{index}\t{read}\t\tSE\tchipseq\t{target}\t{peak_mode}\ths\t"
            f"idx\texp1\t{index}"
        )
    return "\n".join(rows) + "\n"


def _dry_run(
    tmp_path,
    tmp_config,
    run_snakemake,
    *,
    peak_mode="narrow",
    count=2,
    signal_tracks=True,
    printshellcmds=False,
):
    config = {
        "outdir": str(tmp_path / "results"),
        "use_control": False,
        "threads": 1,
        "stage4b": True,
        "stage5": False,
        "qc": {"signal_tracks": signal_tracks},
    }
    _, config_path, _ = tmp_config(
        config=config,
        samples=_samples(tmp_path, peak_mode=peak_mode, count=count),
    )
    extra_args = ["--printshellcmds"] if printshellcmds else None
    result = run_snakemake(config_path, extra_args=extra_args, quiet=False)
    assert result.returncode == 0, result.stdout + result.stderr
    return result.stdout + result.stderr


@pytest.mark.parametrize(
    ("count", "signal_tracks"),
    [(1, True), (2, False)],
    ids=["single-biorep", "signal-tracks-disabled"],
)
def test_pooled_signal_rules_require_multiple_bioreps_and_enabled_gate(
    tmp_path, tmp_config, run_snakemake, count, signal_tracks
):
    output = _dry_run(
        tmp_path,
        tmp_config,
        run_snakemake,
        count=count,
        signal_tracks=signal_tracks,
    )

    assert "pooled_signal_track_fe" not in output
    assert "pooled_signal_track_ppois" not in output


@pytest.mark.parametrize("peak_mode", ["narrow", "broad"])
def test_pooled_signal_rules_schedule_mode_specific_outputs(
    tmp_path, tmp_config, run_snakemake, peak_mode
):
    output = _dry_run(
        tmp_path,
        tmp_config,
        run_snakemake,
        peak_mode=peak_mode,
        printshellcmds=True,
    )

    assert "pooled_signal_track_fe" in output
    assert "pooled_signal_track_ppois" in output
    assert "exp1.pooled.FE.bdg" in output
    assert "exp1.pooled.ppois.bdg" in output
    assert "exp1_pooled" in output
    assert "-B" in output
    if peak_mode == "broad":
        assert "broadPeak" in output
