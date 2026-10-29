"""Behavior tests for the CUT&Tag retain/exclude comparison."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

_SCRIPT = (
    Path(__file__).resolve().parents[2] / "scripts/compare_cuttag_duplicate_arms.py"
)
_SPEC = importlib.util.spec_from_file_location("compare_cuttag_duplicate_arms", _SCRIPT)
cmp = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(cmp)


def _qc(policy: str, frip: str = "0.100000") -> dict[str, str]:
    return {
        "sample": "T2",
        "duplicate_analysis": policy,
        "resolved_max_insert": "700",
        "blacklist_mode": "analysis",
        "organelle_mode": "exclude",
        "layout": "PE",
        "assay": "cuttag",
        "peak_mode": "narrow",
        "genome": "hs",
        "analysis_valid_fragments": "10",
        "analysis_input_stage": "analysis_bam",
        "percent_duplication": "0.200000",
        "nrf": "0.900000",
        "pbc1": "0.800000",
        "pbc2": "2.000000",
        "frip": frip,
        "total_reads": "100",
        "reads_in_peaks": "10",
    }


def test_merge_includes_overlap_and_book_ended_intervals():
    merged = cmp.merge_intervals(
        [("chr1", 0, 10), ("chr1", 10, 15), ("chr1", 20, 25), ("chr2", 0, 5)]
    )
    assert merged == [("chr1", 0, 15), ("chr1", 20, 25), ("chr2", 0, 5)]


def test_identical_peaks_have_unit_jaccard_and_full_overlap():
    peaks = [("chr1", 0, 100), ("chr1", 200, 250)]
    result = cmp.bp_jaccard(peaks, list(peaks))
    assert result["jaccard"] == 1
    assert result["intersection_bp"] == 150
    assert cmp.overlap_fraction(peaks, peaks) == 1


def test_disjoint_peaks_have_zero_jaccard_and_zero_overlap():
    retain = [("chr1", 0, 10)]
    exclude = [("chr1", 20, 30)]
    result = cmp.bp_jaccard(retain, exclude)
    assert result["jaccard"] == 0
    assert result["union_bp"] == 20
    assert cmp.overlap_fraction(retain, exclude) == 0
    assert cmp.overlap_fraction(exclude, retain) == 0


def test_partial_overlap_jaccard_uses_merged_bases():
    retain = [("chr1", 0, 100), ("chr1", 90, 120)]
    exclude = [("chr1", 100, 200)]
    result = cmp.bp_jaccard(retain, exclude)
    assert result["retain_merged_bp"] == 120
    assert result["intersection_bp"] == 20
    assert result["union_bp"] == 200
    assert result["jaccard"] == pytest.approx(0.1)
    assert cmp.overlap_fraction(retain, exclude) == 0.5


def test_empty_query_overlap_is_undefined_not_zero():
    assert cmp.overlap_fraction([], [("chr1", 0, 10)]) is None
    assert cmp.bp_jaccard([], [])["jaccard"] is None


def test_frip_effects_use_fixed_epsilon_and_percentage_points():
    effects = cmp.frip_effects(0.04, 0.05)
    assert effects["absolute_difference_percentage_points"] == pytest.approx(1.0)
    assert effects["relative_change"] == pytest.approx(0.25)
    assert effects["log2_ratio"] == pytest.approx(
        cmp._log2((0.05 + 1e-6) / (0.04 + 1e-6))
    )


def test_frip_rejects_values_outside_unit_interval():
    with pytest.raises(ValueError):
        cmp.frip_effects(1.2, 0.1)


def test_density_uses_common_universe_and_remaining_reference_bases():
    density = cmp.density_enrichment(20, 100, 10, 110)
    assert density["reads_per_peak_base"] == 2
    assert density["reads_per_offpeak_base"] == pytest.approx(0.8)
    assert density["enrichment"] == pytest.approx(2.5)
    assert density["offpeak_bp"] == 100


def test_density_is_undefined_when_the_universe_covers_the_reference():
    density = cmp.density_enrichment(10, 10, 5, 5)
    assert density["reads_per_offpeak_base"] is None
    assert density["enrichment"] is None


def test_configuration_guard_rejects_a_drifting_fixed_field():
    exclude = _qc("exclude")
    exclude["resolved_max_insert"] = "500"
    with pytest.raises(ValueError, match="resolved_max_insert"):
        cmp.assert_fixed_configuration(_qc("retain"), exclude)


def test_record_marks_own_peak_frip_as_not_cross_arm_comparable():
    record = cmp.comparison_record(
        _qc("retain"),
        _qc("exclude", frip="0.200000"),
        [("chr1", 0, 50)],
        [("chr1", 25, 75)],
    )
    assert record["own_peak_frip"]["retain"]["comparable_across_arms"] is False
    assert record["common_universe_frip"] == {"status": "not_measured"}
    assert record["replicate_consistency"] == "not_applicable_single_pair"
    assert record["peak_count"]["role"] == "descriptive_only"
    assert record["merged_bp_jaccard"]["jaccard"] == pytest.approx(25 / 75)


def test_cli_writes_common_universe_effects(tmp_path: Path, monkeypatch):
    def write_qc(path: Path, policy: str):
        header = "\t".join(_qc(policy))
        values = "\t".join(_qc(policy).values())
        path.write_text(f"{header}\n{values}\n", encoding="utf-8")

    retain_qc = tmp_path / "retain.qc.tsv"
    exclude_qc = tmp_path / "exclude.qc.tsv"
    write_qc(retain_qc, "retain")
    write_qc(exclude_qc, "exclude")
    retain_peaks = tmp_path / "retain.narrowPeak"
    exclude_peaks = tmp_path / "exclude.narrowPeak"
    retain_peaks.write_text("chr1\t0\t100\t.\t0\t.\t1\t1\t1\t50\n", encoding="utf-8")
    exclude_peaks.write_text("chr1\t50\t150\t.\t0\t.\t1\t1\t1\t50\n", encoding="utf-8")
    sizes = tmp_path / "chrom.sizes"
    sizes.write_text("chr1\t1000\n", encoding="utf-8")

    def fake_frip(script, sample, bam, peaks, output):
        reads = "40" if "retain" in output.name else "20"
        output.write_text(
            "sample\ttotal_reads\treads_in_peaks\tfrip\tbam\tpeaks\n"
            f"{sample}\t100\t{reads}\t{int(reads) / 100:.6f}\t{bam}\t{peaks}\n",
            encoding="utf-8",
        )
        return cmp.load_qc_summary(output, sample)

    monkeypatch.setattr(cmp, "count_common_frip", fake_frip)
    output = tmp_path / "comparison.json"
    assert (
        cmp.main(
            [
                "--sample",
                "T2",
                "--retain-qc",
                str(retain_qc),
                "--exclude-qc",
                str(exclude_qc),
                "--retain-peaks",
                str(retain_peaks),
                "--exclude-peaks",
                str(exclude_peaks),
                "--chrom-sizes",
                str(sizes),
                "--retain-bam",
                str(tmp_path / "retain.bam"),
                "--exclude-bam",
                str(tmp_path / "exclude.bam"),
                "--output",
                str(output),
            ]
        )
        == 0
    )
    record = json.loads(output.read_text(encoding="utf-8"))
    assert record["common_universe_frip"]["status"] == "measured"
    assert record["common_universe_bp"] == 150
    effects = record["common_universe_frip"]["effects_exclude_minus_retain"]
    assert effects["absolute_difference_percentage_points"] == pytest.approx(-20)
    assert record["signal_background"]["retain"]["enrichment"] == pytest.approx(
        (40 / 150) / (60 / 850)
    )
