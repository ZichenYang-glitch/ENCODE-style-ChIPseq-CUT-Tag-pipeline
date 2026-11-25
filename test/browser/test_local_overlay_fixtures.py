"""Check the fixture contract in ci-fast and run native format IO in the scientific tier.

For the two real_execution tests, set WASHU_FIXTURE_PYTHON to the Python
from chipseq.lock. The test runner itself only needs the ci-fast environment.
"""

from __future__ import annotations

import gzip
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

from encode_pipeline.services.browser_derivatives import (
    MAX_SOURCE_BYTES,
    MAX_SOURCE_RECORDS,
)
from encode_pipeline.services.browser_bgzf import IndexedText, query_indexed_bed


_GENERATOR = Path(__file__).with_name("local_overlay_fixtures.py")
_CONSENSUS = Path(__file__).parents[1] / "fixtures/browser/local-overlay-consensus.json"
_EGFR_ROWS = [
    "chr7\t55101000\t55102000\tchr7:55106000-55107000,1",
    "chr7\t55106000\t55107000\tchr7:55101000-55102000,1",
]
_TP53_ROWS = [
    "chr17\t7672000\t7673000\tchr17:7677000-7678000,1",
    "chr17\t7677000\t7678000\tchr17:7672000-7673000,1",
]


def _generate(dest: Path, *, filler_pairs: int | None = None) -> dict:
    python = os.environ.get("WASHU_FIXTURE_PYTHON", sys.executable)
    command = [python, str(_GENERATOR), "--dest", str(dest)]
    if filler_pairs is not None:
        command.extend(["--filler-pairs", str(filler_pairs)])
    subprocess.run(command, check=True)
    return json.loads((dest / "observation.json").read_text(encoding="utf-8"))


def _inside(gene: dict, start: int, end: int) -> bool:
    return (
        gene["zero_based_half_open"][0] <= start
        and end <= gene["zero_based_half_open"][1]
    )


def test_chr7_matches_the_drawn_fixture_and_tp53_stays_inside_the_gene() -> None:
    spec = json.loads(_CONSENSUS.read_text(encoding="utf-8"))
    assert spec["filenames"]["longrange"] == "pet.longrange.gz"
    assert spec["filenames"]["longrange_index"] == "pet.longrange.gz.tbi"
    assert spec["filenames"]["signal"] == "signal.bw"
    assert spec["filenames"]["peaks"] == "peak.bed"
    assert spec["longrange_pairs"]["egfr"] == [
        {"a": [55101000, 55102000], "b": [55106000, 55107000]}
    ]
    assert spec["peaks"][0] == {
        "chrom": "chr7",
        "start": 55103200,
        "end": 55104100,
        "name": "Peak-one",
    }
    assert spec["signals"][0] == {
        "chrom": "chr7",
        "start": 55104400,
        "end": 55105900,
        "value": 4.5,
    }
    assert spec["signals"][1]["value"] == 2.25
    assert spec["peaks"][1]["name"] == "Peak-tp53"
    egfr = spec["genes"]["EGFR"]
    tp53 = spec["genes"]["TP53"]
    for key, gene in (("egfr", egfr), ("tp53", tp53)):
        pair = spec["longrange_pairs"][key][0]
        assert _inside(gene, pair["a"][0], pair["a"][1])
        assert _inside(gene, pair["b"][0], pair["b"][1])
    for item in spec["signals"] + spec["peaks"]:
        gene = egfr if item["chrom"] == "chr7" else tp53
        assert _inside(gene, item["start"], item["end"])


@pytest.mark.real_execution
def test_small_build_is_deterministic_and_rejects_a_foreign_index(
    tmp_path: Path,
) -> None:
    dest = tmp_path / "a"
    first = _generate(dest, filler_pairs=8)
    second = _generate(tmp_path / "b", filler_pairs=8)
    assert first["data_sha256"] == second["data_sha256"]
    assert first["index_sha256"] == second["index_sha256"]
    names = first["filenames"]
    sample = IndexedText(
        bgzf=(dest / names["longrange"]).read_bytes(),
        tabix=(dest / names["longrange_index"]).read_bytes(),
        record_count=first["record_count"],
    )
    for chrom, start, end, expected in (
        ("chr7", 55101000, 55107000, _EGFR_ROWS),
        ("chr17", 7672000, 7678000, _TP53_ROWS),
    ):
        rows, truncated = query_indexed_bed(sample, chrom, start, end, limit=10)
        assert list(rows) == expected
        assert not truncated
    assert first["window"]["hits"] == _EGFR_ROWS
    assert first["tp53_window"]["hits"] == _TP53_ROWS
    assert first["signal_intervals"] == [
        ["chr7", 55104400, 55105900, 4.5],
        ["chr17", 7675400, 7676900, 2.25],
    ]
    peaks = (tmp_path / "a" / names["peaks"]).read_text(encoding="ascii")
    assert peaks == (
        "chr7\t55103200\t55104100\tPeak-one\nchr17\t7674200\t7675100\tPeak-tp53\n"
    )
    wrong = dest / names["other_index"]
    assert gzip.decompress(wrong.read_bytes()).startswith(b"TBI\x01")
    other = IndexedText(
        bgzf=(dest / names["other_longrange"]).read_bytes(),
        tabix=wrong.read_bytes(),
        record_count=1,
    )
    rows, truncated = query_indexed_bed(other, "chr1", 10, 20, limit=10)
    assert list(rows) == ["chr1\t10\t20\tchr1:30-40,1"]
    assert not truncated
    assert first["wrong_index"]["raised"] is True
    assert first["wrong_index"]["hits"] == []


@pytest.mark.real_execution
def test_default_fixture_exceeds_both_online_caps(tmp_path: Path) -> None:
    observation = _generate(tmp_path)
    assert observation["online_source_byte_limit"] == MAX_SOURCE_BYTES
    assert observation["online_source_record_limit"] == MAX_SOURCE_RECORDS
    assert observation["filler_pairs"] == 1_010_000
    assert observation["record_count"] == 2_020_004
    assert observation["record_count"] > MAX_SOURCE_RECORDS
    assert observation["data_file_bytes"] > MAX_SOURCE_BYTES
    assert observation["index_file_bytes"] != observation["data_file_bytes"]
    assert observation["window"]["hits"] == _EGFR_ROWS
    assert observation["tp53_window"]["hits"] == _TP53_ROWS
    assert observation["window"]["elapsed_ms"] >= 0
    assert observation["wrong_index"]["raised"] is True
    saved = json.loads((tmp_path / "observation.json").read_text(encoding="utf-8"))
    assert saved["record_count"] == observation["record_count"]
    assert saved["data_file_bytes"] == observation["data_file_bytes"]
    assert saved["index_file_bytes"] == observation["index_file_bytes"]
