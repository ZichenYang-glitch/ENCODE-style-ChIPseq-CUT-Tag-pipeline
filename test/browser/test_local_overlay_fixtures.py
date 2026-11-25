"""Indexed synthetic tracks, checked by a real tabix fetch."""

from __future__ import annotations

import importlib.util
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
from local_overlay_fixtures import (
    DEFAULT_FILLER_PAIRS,
    ONLINE_SOURCE_BYTES,
    ONLINE_SOURCE_RECORDS,
    load_consensus,
    pair_lines,
)


def _inside(gene: dict, start: int, end: int) -> bool:
    return (
        gene["zero_based_half_open"][0] <= start
        and end <= gene["zero_based_half_open"][1]
    )


def test_chr7_matches_the_drawn_fixture_and_tp53_stays_inside_the_gene() -> None:
    spec = load_consensus()
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
    assert ONLINE_SOURCE_BYTES == MAX_SOURCE_BYTES
    assert ONLINE_SOURCE_RECORDS == MAX_SOURCE_RECORDS


@pytest.fixture(scope="module")
def fixture_python() -> str:
    # The browser data writers live in workflow/envs/chipseq.lock; ci-fast
    # can still check the fixture coordinates and online converter contract.
    interpreter = os.environ.get("WASHU_FIXTURE_PYTHON")
    if interpreter:
        return interpreter
    if all(importlib.util.find_spec(name) for name in ("pysam", "pyBigWig")):
        return sys.executable
    pytest.skip("set WASHU_FIXTURE_PYTHON to the locked chipseq Python")


def _build(dest: Path, pairs: int, interpreter: str) -> dict:
    subprocess.run(
        [
            interpreter,
            str(Path(__file__).with_name("local_overlay_fixtures.py")),
            "--dest",
            str(dest),
            "--filler-pairs",
            str(pairs),
        ],
        check=True,
    )
    return json.loads((dest / "observation.json").read_text(encoding="utf-8"))


def test_small_build_is_deterministic_and_rejects_a_foreign_index(
    tmp_path: Path,
    fixture_python: str,
) -> None:
    spec = load_consensus()
    first = _build(tmp_path / "a", 8, fixture_python)
    second = _build(tmp_path / "b", 8, fixture_python)
    assert first["data_sha256"] == second["data_sha256"]
    assert first["index_sha256"] == second["index_sha256"]
    names = spec["filenames"]
    expected = pair_lines("chr7", spec["longrange_pairs"]["egfr"][0])
    assert first["window"]["hits"] == expected
    assert first["tp53_window"]["hits"] == pair_lines(
        "chr17", spec["longrange_pairs"]["tp53"][0]
    )
    assert first["signal_intervals"] == [
        ["chr7", 55104400, 55105900, 4.5],
        ["chr17", 7675400, 7676900, 2.25],
    ]
    peaks = (tmp_path / "a" / names["peaks"]).read_text(encoding="ascii")
    assert peaks == (
        "chr7\t55103200\t55104100\tPeak-one\nchr17\t7674200\t7675100\tPeak-tp53\n"
    )
    assert first["wrong_index"]["raised"] is True
    assert first["wrong_index"]["hits"] == []


def test_default_fixture_exceeds_both_online_caps(
    tmp_path: Path, fixture_python: str
) -> None:
    spec = load_consensus()
    observation = _build(tmp_path, DEFAULT_FILLER_PAIRS, fixture_python)
    assert observation["filler_pairs"] == 1_010_000
    assert observation["record_count"] == 2_020_004
    assert observation["record_count"] > MAX_SOURCE_RECORDS
    assert observation["data_file_bytes"] > MAX_SOURCE_BYTES
    assert observation["index_file_bytes"] != observation["data_file_bytes"]
    assert observation["window"]["hits"] == pair_lines(
        "chr7", spec["longrange_pairs"]["egfr"][0]
    )
    assert observation["tp53_window"]["hits"] == pair_lines(
        "chr17", spec["longrange_pairs"]["tp53"][0]
    )
    assert observation["window"]["elapsed_ms"] >= 0
    assert observation["wrong_index"]["raised"] is True
    saved = json.loads((tmp_path / "observation.json").read_text(encoding="utf-8"))
    assert saved["record_count"] == observation["record_count"]
    assert saved["data_file_bytes"] == observation["data_file_bytes"]
    assert saved["index_file_bytes"] == observation["index_file_bytes"]
