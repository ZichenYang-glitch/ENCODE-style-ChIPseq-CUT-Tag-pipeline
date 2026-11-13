"""Display-derivative coordinates, scores, and indexed window queries."""

from __future__ import annotations

import gzip

import pytest

from encode_pipeline.services.browser_bgzf import IndexedText, query_indexed_bed
from encode_pipeline.services.browser_derivatives import (
    MAX_SOURCE_BYTES,
    DerivativeLimits,
    DisplayDerivativeError,
    build_display_derivative,
)


REVISION = "artifactrev-" + "a" * 64
GENERATION = "artifactgen-" + "b" * 64


def _build(
    kind: str,
    text: str,
    *,
    output_type: str = "hitrac_bedpe_all",
    limits: DerivativeLimits | None = None,
):
    return build_display_derivative(
        kind=kind,
        source=text.encode("ascii"),
        source_artifact_id="artifact-pet",
        source_revision=REVISION,
        source_generation=GENERATION,
        source_output_type=output_type,
        limits=limits,
    )


def test_pet_derivative_keeps_pairs_without_copying_mapq() -> None:
    built = _build(
        "pet_records",
        "\n".join(
            (
                "chr1\t1200\t1400\tchr1\t1600\t1800\tread\t211\t+\t+",
                "chr1\t3000\t3100\tchr1\t3000\t3100\tread\t7\t+\t+",
                "chr1\t5000\t5200\tchr2\t100\t180\tread\t40\t+\t-",
            )
        )
        + "\n",
    )

    decoded = built.data
    assert b"211" not in decoded
    assert b"\t40\t" not in decoded
    assert built.provenance["pair_count"] == 3
    assert built.provenance["anchor_row_count"] == 5
    assert built.provenance["score_meaning"] == "display_flag"
    assert built.provenance["coordinate_convention"] == "bed-0-based-half-open"
    rows, truncated = query_indexed_bed(
        IndexedText(built.data, built.index, 5),
        "chr1",
        1000,
        2000,
        limit=10,
    )
    assert truncated is False
    assert rows == (
        "chr1\t1200\t1400\tchr1:1600-1800,1",
        "chr1\t1600\t1800\tchr1:1200-1400,1",
    )
    assert sum(row.endswith(",1") for row in rows) == 2


def test_window_without_records_is_empty_not_an_error() -> None:
    built = _build(
        "pet_records",
        "chr1\t1200\t1400\tchr1\t1600\t1800\tread\t211\t+\t+\n",
    )
    rows, truncated = query_indexed_bed(
        IndexedText(built.data, built.index, 2),
        "chr1",
        8000,
        9000,
        limit=10,
    )
    assert rows == ()
    assert truncated is False


def test_peak_score_is_not_used_for_pet_records() -> None:
    built = _build(
        "peak_intervals",
        "chr1\t1500\t1700\tpeak-a\t80\t+\n",
        output_type="macs3_peak",
    )
    rows, _truncated = query_indexed_bed(
        IndexedText(built.data, built.index, 1),
        "chr1",
        1400,
        1800,
        limit=5,
    )
    assert rows == ("chr1\t1500\t1700\tpeak-a\t80\t+",)


def test_source_over_the_display_cap_is_rejected() -> None:
    with pytest.raises(DisplayDerivativeError) as caught:
        build_display_derivative(
            kind="pet_records",
            source=b"x" * (MAX_SOURCE_BYTES + 1),
            source_artifact_id="artifact-pet",
            source_revision=REVISION,
            source_generation=GENERATION,
            source_output_type="hitrac_bedpe_all",
        )
    assert caught.value.reason == "source_too_large"


def test_gzip_expansion_stops_at_the_decompressed_cap() -> None:
    line = b"chr1\t1000\t1100\tchr1\t2000\t2100\tread\t40\t+\t-\n"
    raw = line * 20_000
    source = gzip.compress(raw, compresslevel=9)
    assert len(source) < 64 * 1024
    assert len(raw) > 8192
    with pytest.raises(DisplayDerivativeError) as caught:
        build_display_derivative(
            kind="pet_records",
            source=source,
            source_artifact_id="artifact-pet",
            source_revision=REVISION,
            source_generation=GENERATION,
            source_output_type="hitrac_bedpe_all",
            limits=DerivativeLimits(decompressed_bytes=8192),
        )
    assert caught.value.reason == "source_too_large"


def test_line_over_the_cap_is_rejected_before_a_derivative_exists() -> None:
    with pytest.raises(DisplayDerivativeError) as caught:
        _build(
            "pet_records",
            "chr1\t1000\t1100\tchr1\t2000\t2100\tread\t40\t+\t-\n",
            limits=DerivativeLimits(line_bytes=32),
        )
    assert caught.value.reason == "source_too_large"


def test_display_rows_include_both_pet_anchors() -> None:
    with pytest.raises(DisplayDerivativeError) as caught:
        _build(
            "pet_records",
            "chr1\t1200\t1400\tchr1\t1600\t1800\tread\t211\t+\t+\n",
            limits=DerivativeLimits(output_rows=1),
        )
    assert caught.value.reason == "source_too_large"
    built = _build(
        "pet_records",
        "chr1\t1200\t1400\tchr1\t1600\t1800\tread\t211\t+\t+\n",
        limits=DerivativeLimits(output_rows=2, source_records=1),
    )
    assert built.provenance["pair_count"] == 1
    assert built.provenance["anchor_row_count"] == 2
    assert b"211" not in built.data


def test_source_record_cap_stops_before_sorting_extra_pairs() -> None:
    with pytest.raises(DisplayDerivativeError) as caught:
        _build(
            "pet_records",
            "\n".join(
                (
                    "chr1\t1200\t1400\tchr1\t1600\t1800\tread\t10\t+\t+",
                    "chr1\t3000\t3100\tchr1\t3200\t3300\tread\t11\t+\t+",
                )
            ),
            limits=DerivativeLimits(source_records=1, output_rows=10),
        )
    assert caught.value.reason == "source_too_large"
