"""Behavior tests for the CUT&Tag max_insert 500 vs 700 template contrast."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest


_SCRIPT = (
    Path(__file__).resolve().parents[2]
    / "scripts"
    / "compare_cuttag_insert_stability.py"
)
_SPEC = importlib.util.spec_from_file_location(
    "compare_cuttag_insert_stability", _SCRIPT
)
cmp = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(cmp)

_COLUMNS = (
    "qname",
    "r1_chrom",
    "r1_pos",
    "r2_chrom",
    "r2_pos",
    "length",
    "mapq",
)


def _template(
    qname: str,
    *,
    r1_chrom: str = "chr1",
    r1_pos: int = 100,
    r2_chrom: str = "chr1",
    r2_pos: int = 250,
    length: int = 150,
    mapq: int = 40,
) -> dict[str, object]:
    return {
        "qname": qname,
        "r1_chrom": r1_chrom,
        "r1_pos": r1_pos,
        "r2_chrom": r2_chrom,
        "r2_pos": r2_pos,
        "length": length,
        "mapq": mapq,
    }


def _write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    lines = ["\t".join(_COLUMNS)]
    for row in rows:
        lines.append("\t".join(str(row[column]) for column in _COLUMNS))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def test_newly_admitted_templates_are_counted_by_length_and_mapq_bins():
    admitted = [
        _template("a", length=0, mapq=29),
        _template("b", length=150, mapq=0),
        _template("c", length=151, mapq=30),
        _template("d", length=300, mapq=29),
        _template("e", length=301, mapq=60),
        _template("f", length=500, mapq=29),
        _template("g", length=501, mapq=30),
        _template("h", length=700, mapq=10),
        _template("i", length=701, mapq=30),
        _template("j", length=1000, mapq=60),
    ]
    counts = cmp.classify_new_templates([], admitted)
    assert counts["total"] == 10
    assert counts["by_length_bin"] == {
        "<=150": 2,
        "151-300": 2,
        "301-500": 2,
        "501-700": 2,
        ">700": 2,
    }
    assert counts["by_mapq_bin"] == {"<30": 5, ">=30": 5}
    assert counts["by_length_and_mapq"]["<=150"] == {"<30": 2, ">=30": 0}
    assert counts["by_length_and_mapq"]["151-300"] == {"<30": 1, ">=30": 1}
    assert counts["by_length_and_mapq"]["301-500"] == {"<30": 1, ">=30": 1}
    assert counts["by_length_and_mapq"]["501-700"] == {"<30": 1, ">=30": 1}
    assert counts["by_length_and_mapq"][">700"] == {"<30": 0, ">=30": 2}
    assert sum(counts["by_length_bin"].values()) == counts["total"]
    assert sum(counts["by_mapq_bin"].values()) == counts["total"]


def test_templates_present_in_the_500_arm_are_not_newly_admitted():
    shared_500 = _template("shared", length=100, mapq=40)
    shared_700 = _template("shared", length=800, mapq=10, r1_pos=5)
    only_500 = _template("only500", length=80, mapq=60)
    only_700 = _template("only700", length=620, mapq=55)
    counts = cmp.classify_new_templates(
        [shared_500, only_500],
        [shared_700, only_700],
    )
    assert counts["total"] == 1
    assert counts["by_length_bin"]["501-700"] == 1
    assert counts["by_length_bin"][">700"] == 0
    assert counts["by_length_bin"]["<=150"] == 0
    assert counts["by_mapq_bin"] == {"<30": 0, ">=30": 1}
    assert counts["by_length_and_mapq"]["501-700"][">=30"] == 1


def test_empty_arms_report_zero_admissions_and_every_bin():
    counts = cmp.classify_new_templates([], [])
    assert counts["total"] == 0
    assert set(counts["by_length_bin"]) == {
        "<=150",
        "151-300",
        "301-500",
        "501-700",
        ">700",
    }
    assert all(value == 0 for value in counts["by_length_bin"].values())
    assert counts["by_mapq_bin"] == {"<30": 0, ">=30": 0}


def test_identical_templates_contribute_no_change():
    record = _template("kept", length=250, mapq=42, r1_pos=10, r2_pos=260)
    same = cmp.short_template_stability([record], [dict(record)])
    assert same == {"denominator": 1, "changed": 0, "fraction": 0.0}
    drifted = _template("kept", length=480, mapq=60, r1_pos=10, r2_pos=260)
    placement = cmp.short_template_stability([record], [drifted])
    assert placement == {"denominator": 1, "changed": 0, "fraction": 0.0}


def test_r1_or_r2_coordinate_change_counts_as_changed():
    r1_500 = _template("r1", r1_pos=100, r2_pos=200)
    r1_700 = _template("r1", r1_pos=110, r2_pos=200)
    r1 = cmp.short_template_stability([r1_500], [r1_700])
    assert r1 == {"denominator": 1, "changed": 1, "fraction": 1.0}
    r2_500 = _template("r2", r1_pos=100, r2_pos=200)
    r2_700 = _template("r2", r1_pos=100, r2_pos=240)
    r2 = cmp.short_template_stability([r2_500], [r2_700])
    assert r2 == {"denominator": 1, "changed": 1, "fraction": 1.0}


def test_changed_mate_contig_counts_when_positions_are_equal():
    reference = _template(
        "mates",
        r1_chrom="chr1",
        r2_chrom="chr1",
        r1_pos=100,
        r2_pos=200,
    )
    observed = _template(
        "mates",
        r1_chrom="chr1",
        r2_chrom="chr2",
        r1_pos=100,
        r2_pos=200,
    )
    result = cmp.short_template_stability([reference], [observed])
    assert result == {"denominator": 1, "changed": 1, "fraction": 1.0}
    swapped = _template(
        "mates",
        r1_chrom="chr2",
        r2_chrom="chr1",
        r1_pos=100,
        r2_pos=100,
    )
    reference_equal_positions = _template(
        "mates",
        r1_chrom="chr1",
        r2_chrom="chr2",
        r1_pos=100,
        r2_pos=100,
    )
    swap = cmp.short_template_stability([reference_equal_positions], [swapped])
    assert swap["changed"] == 1


def test_low_mapq_or_long_templates_are_excluded_from_stability():
    stable = _template("stable", mapq=40, length=200)
    low = _template("low", mapq=29, length=100, r1_pos=100)
    long = _template("long", mapq=60, length=501, r1_pos=100)
    at_mapq = _template("at-mapq", mapq=30, length=100, r1_pos=10)
    at_length = _template("at-length", mapq=40, length=500, r1_pos=10)
    arm500 = [stable, low, long, at_mapq, at_length]
    arm700 = [
        _template("stable", mapq=40, length=200),
        _template("low", mapq=29, length=100, r1_pos=999),
        _template("long", mapq=60, length=501, r1_pos=999),
        _template("at-mapq", mapq=30, length=100, r1_pos=11),
        _template("at-length", mapq=40, length=500, r1_pos=11),
    ]
    result = cmp.short_template_stability(arm500, arm700)
    assert result == {"denominator": 3, "changed": 2, "fraction": pytest.approx(2 / 3)}


def test_stability_fraction_is_none_when_the_denominator_is_zero():
    low = _template("low", mapq=10, length=100)
    long = _template("long", mapq=40, length=600)
    result = cmp.short_template_stability([low, long], [low, long])
    assert result == {"denominator": 0, "changed": 0, "fraction": None}
    assert cmp.short_template_stability([], []) == {
        "denominator": 0,
        "changed": 0,
        "fraction": None,
    }


def test_700_only_template_is_admitted_and_not_a_stability_failure():
    stable = _template("stable", length=200, mapq=40)
    admitted = _template("new", length=650, mapq=12, r1_pos=9, r2_pos=659)
    counts = cmp.classify_new_templates([stable], [stable, admitted])
    assert counts["total"] == 1
    assert counts["by_length_bin"]["501-700"] == 1
    assert counts["by_mapq_bin"]["<30"] == 1
    assert counts["by_length_and_mapq"]["501-700"]["<30"] == 1
    stability = cmp.short_template_stability([stable], [stable, admitted])
    assert stability == {"denominator": 1, "changed": 0, "fraction": 0.0}
    only_new = cmp.short_template_stability([], [admitted])
    assert only_new == {"denominator": 0, "changed": 0, "fraction": None}


def test_qualifying_template_missing_from_700_is_not_a_coordinate_change():
    kept = _template("kept")
    lost = _template("lost", r1_pos=50, r2_pos=80)
    result = cmp.short_template_stability([kept, lost], [kept])
    assert result == {"denominator": 2, "changed": 0, "fraction": 0.0}


def test_duplicate_qname_in_one_arm_is_rejected():
    duplicate = [_template("t"), _template("t", r1_pos=1)]
    with pytest.raises(ValueError, match="duplicate template qname: t"):
        cmp.classify_new_templates(duplicate, [])
    with pytest.raises(ValueError, match="duplicate template qname: t"):
        cmp.short_template_stability([], duplicate)


def test_cli_summarizes_two_template_tsvs(tmp_path: Path):
    arm500 = [
        _template("same", length=180, mapq=40, r1_pos=100, r2_pos=280),
        _template("moved", length=200, mapq=50, r1_pos=100, r2_pos=300),
    ]
    arm700 = [
        _template("same", length=180, mapq=40, r1_pos=100, r2_pos=280),
        _template("moved", length=200, mapq=50, r1_pos=130, r2_pos=300),
        _template("new", length=720, mapq=20, r1_pos=1, r2_pos=721),
    ]
    arm500_path = tmp_path / "arm500.tsv"
    arm700_path = tmp_path / "arm700.tsv"
    _write_tsv(arm500_path, arm500)
    _write_tsv(arm700_path, arm700)
    output = tmp_path / "evidence" / "insert-stability.json"
    assert not output.parent.exists()
    assert (
        cmp.main(
            [
                "--arm-500",
                str(arm500_path),
                "--arm-700",
                str(arm700_path),
                "--output",
                str(output),
            ]
        )
        == 0
    )
    payload = json.loads(output.read_text(encoding="utf-8"))
    assert payload["reference_arm"] == "max_insert_500"
    assert payload["candidate_arm"] == "max_insert_700"
    assert payload["newly_admitted_templates"] == cmp.classify_new_templates(
        arm500, arm700
    )
    assert payload["short_template_stability"] == cmp.short_template_stability(
        arm500, arm700
    )
    assert payload["newly_admitted_templates"]["total"] == 1
    assert payload["newly_admitted_templates"]["by_length_bin"][">700"] == 1
    assert payload["newly_admitted_templates"]["by_mapq_bin"]["<30"] == 1
    assert payload["short_template_stability"]["denominator"] == 2
    assert payload["short_template_stability"]["changed"] == 1
    assert payload["short_template_stability"]["fraction"] == pytest.approx(0.5)


def test_template_tsv_rejects_an_unexpected_header(tmp_path: Path):
    path = tmp_path / "templates.tsv"
    path.write_text("qname\tlength\nread1\t100\n", encoding="utf-8")
    with pytest.raises(ValueError, match="unexpected TSV header"):
        cmp.load_templates(path)
