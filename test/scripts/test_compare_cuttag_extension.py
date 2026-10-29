"""Behavior tests for CUT&Tag replicate and max_insert metrics."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path


_SCRIPT = Path(__file__).resolve().parents[2] / "scripts/compare_cuttag_extension.py"
_SPEC = importlib.util.spec_from_file_location("compare_cuttag_extension", _SCRIPT)
ext = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(ext)


def _pair(qname: str, r1_pos: int, r2_pos: int, mapq: int, length: int) -> str:
    return (
        f"{qname}\t99\tchr1\t{r1_pos}\t{mapq}\t25M\t=\t{r2_pos}\t{length}\t*\t*\n"
        f"{qname}\t147\tchr1\t{r2_pos}\t{mapq}\t25M\t=\t{r1_pos}\t{-length}\t*\t*\n"
    )


def test_parser_keeps_primary_proper_pairs_and_reports_singletons():
    text = _pair("keep", 10, 40, 40, 55) + (
        "alone\t99\tchr1\t10\t40\t25M\t=\t40\t55\t*\t*\n"
        "skip\t355\tchr1\t10\t40\t25M\t=\t40\t55\t*\t*\n"
    )
    templates, incomplete = ext.templates_from_sam_text(text)
    assert set(templates) == {"keep"}
    assert templates["keep"].length == 55
    assert templates["keep"].r1.pos == 10
    assert templates["keep"].r2.pos == 40
    assert incomplete == ["alone"]


def test_duplicate_primary_is_rejected():
    text = _pair("q", 10, 40, 40, 55) + ("q\t99\tchr1\t12\t40\t25M\t=\t40\t55\t*\t*\n")
    try:
        ext.templates_from_sam_text(text)
    except ValueError as error:
        assert "duplicate primary" in str(error)
    else:
        raise AssertionError("duplicate primary was accepted")


def test_max_insert_reports_admission_loss_and_short_moves(tmp_path: Path):
    shorter = tmp_path / "500.sam"
    longer = tmp_path / "700.sam"
    shorter.write_text(
        _pair("stable", 10, 40, 40, 55)
        + _pair("moved", 10, 80, 40, 90)
        + _pair("lowmapq", 10, 40, 10, 55)
        + _pair("lost", 20, 60, 40, 70),
        encoding="utf-8",
    )
    longer.write_text(
        _pair("stable", 10, 40, 40, 55)
        + _pair("moved", 15, 80, 40, 90)
        + _pair("lowmapq", 10, 40, 10, 55)
        + _pair("long", 100, 700, 8, 650),
        encoding="utf-8",
    )
    output = tmp_path / "max.json"
    assert (
        ext.main(
            [
                "max-insert",
                "--sam-500",
                str(shorter),
                "--sam-700",
                str(longer),
                "--output",
                str(output),
            ]
        )
        == 0
    )
    payload = json.loads(output.read_text(encoding="utf-8"))
    admitted = payload["insert_contrast"]["newly_admitted_templates"]
    assert admitted["total"] == 1
    assert admitted["by_length_bin"]["501-700"] == 1
    assert admitted["by_mapq_bin"]["<30"] == 1
    assert payload["qnames_only_in_500"] == 1
    stability = payload["insert_contrast"]["short_template_stability"]
    assert stability["denominator"] == 3
    assert stability["changed"] == 1
    assert stability["fraction"] == 1 / 3
    assert payload["short_template_trigger"]["triggered"] is True


def test_exact_signal_thresholds_do_not_trigger():
    at_frip_limit = ext.signal_regression_triggered(0.45, 0.5)
    assert abs(at_frip_limit["relative_frip_change"] + 0.1) < 1e-12
    assert at_frip_limit["relative_offpeak_change"] > 0.05
    assert at_frip_limit["triggered"] is False
    beyond_both = ext.signal_regression_triggered(0.15, 0.2)
    assert beyond_both["relative_frip_change"] < -0.10
    assert beyond_both["relative_offpeak_change"] > 0.05
    assert beyond_both["triggered"] is True


def test_zero_reference_frip_is_not_comparable():
    result = ext.signal_regression_triggered(0.0, 0.0)
    assert result["status"] == "not_comparable"
    assert result["triggered"] is False


def test_replicate_jaccard_reuses_merged_base_definition():
    left = [("chr1", 0, 100), ("chr1", 90, 120)]
    right = [("chr1", 100, 200)]
    result = ext.peak_jaccard(left, right)
    assert result["left_merged_bp"] == 120
    assert result["intersection_bp"] == 20
    assert result["jaccard"] == 0.1
    assert result["left_peaks_touched"] == 0.5


def test_coverage_correlation_ignores_shared_empty_windows():
    left = {
        "a": ext.Template("a", ext.Mate("chr1", 1, 40), ext.Mate("chr1", 50, 40), 80),
        "b": ext.Template("b", ext.Mate("chr1", 1, 40), ext.Mate("chr1", 40, 40), 70),
        "c": ext.Template(
            "c", ext.Mate("chr1", 100_001, 40), ext.Mate("chr1", 100_050, 40), 80
        ),
    }
    right = {
        "a": ext.Template("a", ext.Mate("chr1", 1, 40), ext.Mate("chr1", 50, 40), 80),
        "b": ext.Template("b", ext.Mate("chr1", 10, 40), ext.Mate("chr1", 40, 40), 70),
        "d": ext.Template(
            "d", ext.Mate("chr1", 200_001, 40), ext.Mate("chr1", 200_040, 40), 80
        ),
    }
    report = ext.length_binned_coverage(left, right)
    bin_report = report["bins"]["1-100"]
    assert bin_report["empty_windows_included"] is False
    assert bin_report["quality_evidence"] is False
    assert bin_report["occupied_windows"] == 3
    assert bin_report["pearson"] is not None


def test_max_insert_pair_rejects_a_second_variable():
    def row(max_insert: str, duplicate: str = "retain") -> dict[str, str]:
        return {
            "duplicate_analysis": duplicate,
            "resolved_max_insert": max_insert,
            "blacklist_mode": "analysis",
            "organelle_mode": "exclude",
            "layout": "PE",
            "assay": "cuttag",
            "peak_mode": "narrow",
            "genome": "hs",
        }

    ext.assert_single_variable_max_insert(row("500"), row("700"))
    try:
        ext.assert_single_variable_max_insert(row("500", "exclude"), row("700"))
    except ValueError as error:
        assert "duplicate_analysis" in str(error)
    else:
        raise AssertionError("a second changed variable was accepted")


def test_freeze_stays_closed_without_confirmation_or_replicates():
    closed = ext.evaluate_freeze(
        hard_failure=False,
        signal_triggered_on_every_independent_replicate=False,
        independent_replicate_count=1,
        confirmation_evaluated=False,
        short_template_triggered=False,
    )
    assert closed["freeze_default"] is False
    assert "confirmation_not_evaluated" in closed["blocking_reasons"]
    assert "fewer_than_two_independent_replicates" in closed["blocking_reasons"]


def test_replicate_qc_accepts_different_sample_names(tmp_path: Path):
    header = (
        "sample\tduplicate_analysis\tresolved_max_insert\tblacklist_mode\t"
        "organelle_mode\tlayout\tassay\tpeak_mode\tgenome\n"
    )
    left_qc = tmp_path / "left.tsv"
    right_qc = tmp_path / "right.tsv"
    row = "retain\t700\tanalysis\texclude\tPE\tcuttag\tnarrow\ths\n"
    left_qc.write_text(header + "T1\t" + row, encoding="utf-8")
    right_qc.write_text(header + "T2\t" + row, encoding="utf-8")
    left_peaks = tmp_path / "left.narrowPeak"
    right_peaks = tmp_path / "right.narrowPeak"
    left_peaks.write_text("chr1\t0\t10\n", encoding="utf-8")
    right_peaks.write_text("chr1\t0\t10\n", encoding="utf-8")
    output = tmp_path / "replicates.json"
    assert (
        ext.main(
            [
                "replicates",
                "--left-sample",
                "T1",
                "--right-sample",
                "T2",
                "--left-label",
                "rep1",
                "--right-label",
                "rep2",
                "--left-qc",
                str(left_qc),
                "--right-qc",
                str(right_qc),
                "--left-peaks",
                str(left_peaks),
                "--right-peaks",
                str(right_peaks),
                "--output",
                str(output),
            ]
        )
        == 0
    )
    payload = json.loads(output.read_text(encoding="utf-8"))
    assert payload["peak_jaccard"]["jaccard"] == 1
    assert payload["biological_independence"] == "not_established_by_this_script"


def test_cli_rejects_a_single_bam_side(tmp_path: Path):
    try:
        ext.main(
            [
                "max-insert",
                "--sam-500",
                str(tmp_path / "missing.sam"),
                "--output",
                str(tmp_path / "out.json"),
            ]
        )
    except SystemExit as error:
        assert error.code != 0
    else:
        raise AssertionError("one-sided max_insert comparison was accepted")
