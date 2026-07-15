"""Snakemake-to-producer contracts for pooled experiment QC metadata."""

from __future__ import annotations

from pathlib import Path

import pytest


SAMPLE_HEADER = (
    "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\tpeak_mode\tgenome\t"
    "bowtie2_index\texperiment\tbiological_replicate\n"
)


def _sample_sheet(
    tmp_path: Path,
    *,
    target: str,
    peak_mode: str,
    biological_replicates: tuple[int, int],
) -> str:
    rows = [SAMPLE_HEADER.rstrip("\n")]
    for index, biological_replicate in enumerate(biological_replicates, start=1):
        read = tmp_path / f"reads_{index}.fastq"
        read.write_text("", encoding="utf-8")
        rows.append(
            f"T{index}\t{read}\t\tSE\tchipseq\t{target}\t{peak_mode}\ths\tidx"
            f"\texp1\t{biological_replicate}"
        )
    return "\n".join(rows) + "\n"


def _dry_run(
    tmp_path,
    tmp_config,
    run_snakemake,
    *,
    target,
    peak_mode,
    biological_replicates=(1, 2),
    signal_tracks=True,
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
        samples=_sample_sheet(
            tmp_path,
            target=target,
            peak_mode=peak_mode,
            biological_replicates=biological_replicates,
        ),
    )
    result = run_snakemake(
        config_path,
        extra_args=["--printshellcmds"],
        quiet=False,
    )
    output = result.stdout + result.stderr
    assert result.returncode == 0, output[-5000:]
    assert "pooled_experiment_qc_summary" in output
    assert "scripts/pooled_qc_summary.py" in output
    return output


@pytest.mark.parametrize(
    (
        "target",
        "peak_mode",
        "biological_replicates",
        "signal_tracks",
        "expected_arguments",
    ),
    [
        pytest.param(
            "H3K27me3",
            "broad",
            (1, 2),
            True,
            (
                "--peak-mode broad",
                "--inferred-histone-class broad_like",
                "--expected-peak-mode broad",
                "--peak-mode-status ok",
                "--bio-rep-labels 1,2",
            ),
            id="broad-histone-compatible",
        ),
        pytest.param(
            "H3K27me3",
            "narrow",
            (2, 4),
            True,
            (
                "--peak-mode narrow",
                "--inferred-histone-class broad_like",
                "--expected-peak-mode broad",
                "--peak-mode-status mismatch",
                "--bio-rep-labels 2,4",
            ),
            id="narrow-histone-mismatch-noncontiguous-replicates",
        ),
        pytest.param(
            "H3K27ac",
            "narrow",
            (2, 4),
            False,
            (
                "--peak-mode narrow",
                "--inferred-histone-class context_dependent",
                "--expected-peak-mode broad_or_narrow",
                "--peak-mode-status ok",
                "--bio-rep-labels 2,4",
                "--pooled-fe-bdg NA",
                "--pooled-ppois-bdg NA",
            ),
            id="signal-disabled-uses-na",
        ),
    ],
)
def test_pooled_qc_rule_passes_validated_metadata_to_the_producer(
    tmp_path,
    tmp_config,
    run_snakemake,
    target,
    peak_mode,
    biological_replicates,
    signal_tracks,
    expected_arguments,
):
    output = _dry_run(
        tmp_path,
        tmp_config,
        run_snakemake,
        target=target,
        peak_mode=peak_mode,
        biological_replicates=biological_replicates,
        signal_tracks=signal_tracks,
    )

    assert "--n-bioreps 2" in output
    assert '--signal-tracks-status "$SIG_STATUS"' in output
    for argument in expected_arguments:
        assert argument in output
    if signal_tracks:
        assert 'if [[ "True" == "True" ]]' in output
        assert "--pooled-fe-bdg NA" not in output
        assert "--pooled-ppois-bdg NA" not in output
    else:
        assert 'if [[ "False" == "True" ]]' in output
        assert 'SIG_STATUS="disabled"' in output
