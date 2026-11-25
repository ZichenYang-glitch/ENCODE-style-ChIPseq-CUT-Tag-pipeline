"""Windowed reads of the synthetic local tracks, counted from real file IO."""

from __future__ import annotations

import hashlib
from pathlib import Path

import pytest

from encode_pipeline.services.browser_bgzf import (
    BrowserIndexError,
    IndexedText,
    query_indexed_bed,
)
from encode_pipeline.services.browser_derivatives import (
    MAX_SOURCE_BYTES,
    MAX_SOURCE_RECORDS,
    DerivativeLimits,
    DisplayDerivativeError,
    build_display_derivative,
)
from local_overlay_fixtures import (
    DEFAULT_FILLER_RECORDS,
    build_local_overlay_fixtures,
    load_consensus,
    measure_window,
    pair_lines,
    read_signal_intervals,
)


def _inside(gene: dict, start: int, end: int) -> bool:
    return gene["gtf_start"] - 1 <= start and end <= gene["gtf_end"]


def test_planted_features_sit_inside_the_gencode_genes() -> None:
    spec = load_consensus()
    egfr = spec["genes"]["EGFR"]
    tp53 = spec["genes"]["TP53"]
    for pair in spec["longrange_pairs"]["egfr"]:
        assert _inside(egfr, pair["a"][0], pair["a"][1])
        assert _inside(egfr, pair["b"][0], pair["b"][1])
    for pair in spec["longrange_pairs"]["tp53"]:
        assert _inside(tp53, pair["a"][0], pair["a"][1])
        assert _inside(tp53, pair["b"][0], pair["b"][1])
    window = spec["windows"]["egfr"]
    first = spec["longrange_pairs"]["egfr"][0]
    later = spec["longrange_pairs"]["egfr"][1]
    assert window["start"] <= first["a"][0] and window["end"] >= first["b"][1]
    assert later["a"][0] >= window["end"]


def test_online_converter_stops_above_its_record_cap() -> None:
    source = (
        b"chr7\t10\t20\tchr7\t30\t40\tread\t60\t+\t+\n"
        b"chr7\t50\t60\tchr7\t70\t80\tread\t60\t+\t+\n"
    )
    with pytest.raises(DisplayDerivativeError) as caught:
        build_display_derivative(
            kind="pet_records",
            source=source,
            source_artifact_id="artifact-pet",
            source_revision="rev",
            source_generation="gen",
            source_output_type="hitrac_bedpe_all",
            limits=DerivativeLimits(source_records=1),
        )
    assert caught.value.reason == "source_too_large"
    assert MAX_SOURCE_BYTES == 32 * 1024 * 1024
    assert MAX_SOURCE_RECORDS == 2_000_000


def test_small_window_matches_the_repository_reader_and_is_deterministic(
    tmp_path: Path,
) -> None:
    spec = load_consensus()
    first = build_local_overlay_fixtures(tmp_path / "a", filler_records=25)
    second = build_local_overlay_fixtures(tmp_path / "b", filler_records=25)
    assert first["data_sha256"] == second["data_sha256"]
    assert first["index_sha256"] == second["index_sha256"]
    names = spec["filenames"]
    data = (tmp_path / "a" / names["longrange"]).read_bytes()
    index = (tmp_path / "a" / names["longrange_index"]).read_bytes()
    window = spec["windows"]["egfr"]
    expected = pair_lines("chr7", spec["longrange_pairs"]["egfr"][0])
    hits, truncated = query_indexed_bed(
        IndexedText(bgzf=data, tabix=index, record_count=first["record_count"]),
        window["chrom"],
        window["start"],
        window["end"],
        limit=20,
    )
    assert truncated is False
    assert hits == expected
    assert tuple(first["window"]["hits"]) == expected
    assert first["window"]["data_bytes_read"] > 0
    assert first["window"]["data_bytes_read"] < first["data_file_bytes"]
    assert first["miss"]["data_bytes_read"] == 0
    assert first["miss"]["hits"] == []
    later = pair_lines("chr7", spec["longrange_pairs"]["egfr"][1])[0]
    assert later not in first["window"]["hits"]
    intervals = read_signal_intervals(tmp_path / "a" / names["signal"])
    assert (intervals[0][0], intervals[0][1], intervals[0][2], intervals[0][3]) == (
        "chr7",
        55019040,
        55019080,
        4.5,
    )
    assert "egfr-peak\t80" in (tmp_path / "a" / names["peaks"]).read_text(
        encoding="ascii"
    )
    try:
        wrong = measure_window(
            tmp_path / "a" / names["longrange_wrong"],
            tmp_path / "a" / names["longrange_wrong_index"],
            "chr7",
            window["start"],
            window["end"],
        )
    except BrowserIndexError:
        wrong = None
    if wrong is not None:
        assert expected[0] not in wrong.hits


def test_oversized_prebuilt_index_reads_one_window(tmp_path: Path) -> None:
    spec = load_consensus()
    observation = build_local_overlay_fixtures(
        tmp_path, filler_records=DEFAULT_FILLER_RECORDS
    )
    assert observation["record_count"] > MAX_SOURCE_RECORDS
    assert observation["data_file_bytes"] > MAX_SOURCE_BYTES
    assert observation["anchor_records"] == 12
    window = observation["window"]
    assert tuple(window["hits"]) == pair_lines(
        "chr7", spec["longrange_pairs"]["egfr"][0]
    )
    assert window["data_bytes_read"] > 0
    assert window["data_bytes_read"] * 1000 < observation["data_file_bytes"]
    assert window["index_bytes_read"] == observation["index_file_bytes"]
    assert window["index_bytes_read"] != window["data_bytes_read"]
    assert window["elapsed_ms"] > 0
    assert observation["miss"]["data_bytes_read"] == 0
    wide = measure_window(
        tmp_path / spec["filenames"]["longrange"],
        tmp_path / spec["filenames"]["longrange_index"],
        "chr7",
        spec["genes"]["EGFR"]["gtf_start"] - 1,
        spec["genes"]["EGFR"]["gtf_end"],
    )
    assert len(wide.hits) == 6
    assert wide.data_bytes_read * 100 < observation["data_file_bytes"]
    digest = hashlib.sha256(
        (tmp_path / spec["filenames"]["longrange"]).read_bytes()
    ).hexdigest()
    assert digest == observation["data_sha256"]


@pytest.mark.skipif(
    not Path(
        "/home/yangzichen/reference/genomes/GRCh38/GRCh38.annotation.gtf"
    ).is_file(),
    reason="local GENCODE v49 GTF is not on this machine",
)
def test_consensus_matches_the_local_gencode_gtf() -> None:
    path = Path("/home/yangzichen/reference/genomes/GRCh38/GRCh38.annotation.gtf")
    spec = load_consensus()
    wanted = {
        "EGFR": spec["genes"]["EGFR"],
        "TP53": spec["genes"]["TP53"],
        "DDX11L16": spec["genes"]["ambiguous"]["loci"][0],
    }
    found: dict[str, tuple[str, int, int]] = {}
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if "\tgene\t" not in line:
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            name = ""
            for field in parts[8].split(";"):
                field = field.strip()
                if field.startswith('gene_name "'):
                    name = field[len('gene_name "') : -1]
            if name in wanted and name not in found:
                found[name] = (parts[0], int(parts[3]), int(parts[4]))
            if len(found) == len(wanted):
                break
    for name, gene in wanted.items():
        assert found[name] == (gene["chrom"], gene["gtf_start"], gene["gtf_end"])
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    assert digest.hexdigest() == spec["annotation"]["local_sha256"]


def test_revision_bump_is_confined_to_the_changed_run(tmp_path: Path) -> None:
    from washu_local_overlay_server import build_app, revise_run_b

    app, info = build_app(tmp_path / "platform.db", tmp_path / "workspaces")
    before = info["generation_b"]
    after = revise_run_b(app)
    assert after != before
    assert info["run_a"] != info["run_b"] != info["run_mm10"]
    assert info["signal_a"] != info["signal_mm10"]
