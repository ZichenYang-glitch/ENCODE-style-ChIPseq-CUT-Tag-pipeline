"""Indexed synthetic tracks, checked by a real tabix fetch."""

from __future__ import annotations

import importlib.util
import json
import os
from pathlib import Path
import re
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
    interpreter = os.environ.get("WASHU_FIXTURE_PYTHON")
    if interpreter:
        return interpreter
    if all(importlib.util.find_spec(name) for name in ("pysam", "pyBigWig")):
        return sys.executable
    pytest.fail(
        "set WASHU_FIXTURE_PYTHON to the Python from workflow/envs/deeptools.lock"
    )


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


def test_mm10_consensus_places_fixed_values_inside_egfr_and_trp53() -> None:
    human = load_consensus()
    spec = load_consensus("mm10")
    assert spec["assembly"] == "mm10"
    assert spec["annotation"]["release"] == "M25"
    assert spec["annotation"]["assembly"] == "GRCm38"
    assert spec["annotation"]["assembly_accession"] == "GCA_000001635.2"
    assert spec["chrom_sizes"] == {"chr11": 122082543}
    assert spec["filenames"] == human["filenames"]
    egfr = spec["genes"]["Egfr"]
    trp53 = spec["genes"]["Trp53"]
    assert egfr == {
        "gene_id": "ENSMUSG00000020122.16",
        "chrom": "chr11",
        "gtf_start": 16752203,
        "gtf_end": 16918158,
        "strand": "+",
        "mgi_id": "MGI:95294",
        "hub_region": "chr11:16752203-16918158",
        "zero_based_half_open": [16752202, 16918158],
    }
    assert trp53 == {
        "gene_id": "ENSMUSG00000059552.13",
        "chrom": "chr11",
        "gtf_start": 69580359,
        "gtf_end": 69591873,
        "strand": "+",
        "mgi_id": "MGI:98834",
        "hub_region": "chr11:69580359-69591873",
        "zero_based_half_open": [69580358, 69591873],
    }
    for gene in (egfr, trp53):
        assert gene["zero_based_half_open"] == [gene["gtf_start"] - 1, gene["gtf_end"]]
        assert (
            gene["hub_region"]
            == f"{gene['chrom']}:{gene['gtf_start']}-{gene['gtf_end']}"
        )
    assert spec["longrange_pairs"]["egfr"] == [
        {"a": [16780000, 16781000], "b": [16820000, 16821000]}
    ]
    assert spec["longrange_pairs"]["trp53"] == [
        {"a": [69582000, 69583000], "b": [69588000, 69589000]}
    ]
    assert spec["signals"] == [
        {"chrom": "chr11", "start": 16790000, "end": 16791500, "value": 7.5},
        {"chrom": "chr11", "start": 69584000, "end": 69585500, "value": 3.25},
    ]
    assert spec["peaks"] == [
        {"chrom": "chr11", "start": 16785000, "end": 16786000, "name": "Peak-egfr"},
        {"chrom": "chr11", "start": 69586000, "end": 69587000, "name": "Peak-trp53"},
    ]
    assert spec["signals"][0]["value"] != human["signals"][0]["value"]
    # 90 kb still places both anchors. A 55 kb window around the same
    # anchors leaves the longrange track without a placement.
    assert spec["coordinate_probe"] == {
        "gene": "Egfr",
        "region": "chr11:16760001-16850000",
    }
    probe_start, probe_end = (
        int(part)
        for part in spec["coordinate_probe"]["region"].split(":")[1].split("-")
    )
    # Typed text is 1-based and closed. GenomeHub keeps the end and subtracts
    # one from the start, which still covers the 0-based Egfr anchors.
    assert probe_start - 1 <= spec["longrange_pairs"]["egfr"][0]["a"][0]
    assert spec["longrange_pairs"]["egfr"][0]["b"][1] <= probe_end
    assert probe_start != spec["longrange_pairs"]["egfr"][0]["a"][0]
    for key, gene in (("egfr", egfr), ("trp53", trp53)):
        pair = spec["longrange_pairs"][key][0]
        assert _inside(gene, pair["a"][0], pair["a"][1])
        assert _inside(gene, pair["b"][0], pair["b"][1])
    for item in spec["signals"] + spec["peaks"]:
        owners = [
            gene for gene in (egfr, trp53) if _inside(gene, item["start"], item["end"])
        ]
        assert len(owners) == 1


def _build_mouse(dest: Path, interpreter: str, pairs: int | None = None) -> dict:
    command = [
        interpreter,
        str(Path(__file__).with_name("local_overlay_fixtures.py")),
        "--dest",
        str(dest),
        "--assembly",
        "mm10",
    ]
    if pairs is not None:
        command.extend(["--filler-pairs", str(pairs)])
    subprocess.run(command, check=True)
    return json.loads((dest / "observation.json").read_text(encoding="utf-8"))


def test_mm10_build_reads_back_both_genes(
    tmp_path: Path,
    fixture_python: str,
) -> None:
    spec = load_consensus("mm10")
    first = _build_mouse(tmp_path / "a", fixture_python)
    second = _build_mouse(tmp_path / "b", fixture_python)
    assert first["assembly"] == "mm10"
    assert first["filler_pairs"] == 0
    assert first["record_count"] == 4
    assert first["feature_records"] == 4
    assert first["exceeds_source_bytes"] is False
    assert first["exceeds_source_records"] is False
    assert first["data_sha256"] == second["data_sha256"]
    assert first["index_sha256"] == second["index_sha256"]
    assert first["window"]["hits"] == pair_lines(
        "chr11", spec["longrange_pairs"]["egfr"][0]
    )
    assert first["trp53_window"]["hits"] == pair_lines(
        "chr11", spec["longrange_pairs"]["trp53"][0]
    )
    assert first["signal_intervals"] == [
        ["chr11", 16790000, 16791500, 7.5],
        ["chr11", 69584000, 69585500, 3.25],
    ]
    peaks = (tmp_path / "a" / spec["filenames"]["peaks"]).read_text(encoding="ascii")
    assert peaks == (
        "chr11\t16785000\t16786000\tPeak-egfr\nchr11\t69586000\t69587000\tPeak-trp53\n"
    )
    assert first["wrong_index"]["raised"] is True
    assert first["wrong_index"]["hits"] == []


def test_pinned_mouse_genes_match_released_gencode_m25_records() -> None:
    fixture = (
        Path(__file__).resolve().parents[1]
        / "fixtures/browser/gencode-m25-fixture-genes.gtf"
    )
    spec = load_consensus("mm10")
    seen = set()
    for line in fixture.read_text(encoding="utf-8").splitlines():
        if line.startswith("#"):
            continue
        chrom, _, feature, start, end, _, strand, _, attributes = line.split("\t")
        assert feature == "gene"
        attrs = dict(re.findall(r'([A-Za-z0-9_]+) "([^"]*)"', attributes))
        name = attrs["gene_name"]
        assert name not in seen
        seen.add(name)
        gene = spec["genes"][name]
        assert (chrom, int(start), int(end), strand) == (
            gene["chrom"],
            gene["gtf_start"],
            gene["gtf_end"],
            gene["strand"],
        )
        assert attrs["gene_id"] == gene["gene_id"]
        assert attrs["mgi_id"] == gene["mgi_id"]
        assert attrs["gene_type"] == "protein_coding"
    assert seen == set(spec["genes"])

    provenance = (
        Path(__file__).resolve().parents[2]
        / "frontend/src/viewer/search/mm10-gene-index.provenance.json"
    )
    scaffold = json.loads(provenance.read_text(encoding="utf-8"))["viewer_scaffold"]
    lengths = dict(zip(scaffold["sequences"], scaffold["lengths"], strict=True))
    assert lengths["chr11"] == spec["chrom_sizes"]["chr11"]
