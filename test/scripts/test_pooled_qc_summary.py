"""Producer contract for pooled experiment QC summaries."""

import csv
from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import pytest  # noqa: E402

from scripts import pooled_qc_summary  # noqa: E402


EXPECTED_COLUMNS = [
    "experiment",
    "assay",
    "target",
    "inferred_histone_class",
    "expected_peak_mode",
    "configured_peak_mode",
    "peak_mode_status",
    "biological_replicates",
    "biological_replicate_labels",
    "pooled_bam",
    "pooled_peaks",
    "pooled_peak_count",
    "pooled_FE_bdg",
    "pooled_ppois_bdg",
    "signal_tracks_status",
]


@pytest.mark.parametrize(
    ("peak_mode", "suffix"),
    [("narrow", "narrowPeak"), ("broad", "broadPeak")],
)
def test_count_peaks_uses_mode_specific_file(tmp_path, peak_mode, suffix):
    peak_file = tmp_path / f"exp_pooled_peaks.{suffix}"
    peak_file.write_text("peak1\npeak2\npeak3\n", encoding="utf-8")

    assert pooled_qc_summary.count_peaks(tmp_path, peak_mode, "exp") == 3


def test_count_peaks_fails_closed_when_peak_file_is_missing(tmp_path, capsys):
    with pytest.raises(SystemExit, match="1"):
        pooled_qc_summary.count_peaks(tmp_path, "narrow", "missing")

    assert "pooled peak file not found" in capsys.readouterr().err


def test_count_peaks_supports_direct_file_path(tmp_path):
    peak_file = tmp_path / "custom_filtered_peaks.narrowPeak"
    peak_file.write_text("peak1\npeak2\n", encoding="utf-8")
    assert pooled_qc_summary.count_peaks(peak_file, "narrow", "exp") == 2


def test_main_writes_exact_fifteen_column_summary(tmp_path, monkeypatch):
    peak_dir = tmp_path / "peaks"
    peak_dir.mkdir()
    (peak_dir / "exp_pooled_peaks.narrowPeak").write_text(
        "peak1\npeak2\npeak3\n",
        encoding="utf-8",
    )
    output = tmp_path / "nested" / "pooled.tsv"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "pooled_qc_summary.py",
            "--experiment",
            "exp",
            "--assay",
            "chipseq",
            "--target",
            "H3K27me3",
            "--peak-mode",
            "narrow",
            "--inferred-histone-class",
            "broad_like",
            "--expected-peak-mode",
            "broad",
            "--peak-mode-status",
            "mismatch",
            "--n-bioreps",
            "2",
            "--bio-rep-labels",
            "2,4",
            "--pooled-bam",
            "results/exp.pooled.bam",
            "--pooled-peaks-dir",
            str(peak_dir),
            "--output",
            str(output),
        ],
    )

    pooled_qc_summary.main()

    with output.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert list(rows[0]) == EXPECTED_COLUMNS
    assert rows == [
        {
            "experiment": "exp",
            "assay": "chipseq",
            "target": "H3K27me3",
            "inferred_histone_class": "broad_like",
            "expected_peak_mode": "broad",
            "configured_peak_mode": "narrow",
            "peak_mode_status": "mismatch",
            "biological_replicates": "2",
            "biological_replicate_labels": "2,4",
            "pooled_bam": "results/exp.pooled.bam",
            "pooled_peaks": str(peak_dir),
            "pooled_peak_count": "3",
            "pooled_FE_bdg": "NA",
            "pooled_ppois_bdg": "NA",
            "signal_tracks_status": "disabled",
        }
    ]
