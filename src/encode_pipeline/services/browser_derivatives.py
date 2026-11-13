"""Build traceable WashU display files without replacing source artifacts."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import re
import zlib


CONVERTER_ID = "washu-display-derivative-1"
MAX_SOURCE_BYTES = 32 * 1024 * 1024
MAX_DECOMPRESSED_BYTES = 256 * 1024 * 1024
MAX_LINE_BYTES = 64 * 1024
MAX_SOURCE_RECORDS = 2_000_000
MAX_OUTPUT_ROWS = 2_000_000
DISPLAY_FLAG = "1"
_CHROM = re.compile(r"^[A-Za-z0-9_.:-]{1,128}$")
_INTEGER = re.compile(r"(?:0|[1-9][0-9]*)$")
_GZIP_WBITS = 16 + zlib.MAX_WBITS


class DisplayDerivativeError(ValueError):
    """The source cannot become a display derivative."""

    def __init__(self, reason: str) -> None:
        super().__init__(reason)
        self.reason = reason


@dataclass(frozen=True)
class DerivativeLimits:
    """Caps applied while inflating and parsing, before any derivative is returned."""

    source_bytes: int = MAX_SOURCE_BYTES
    decompressed_bytes: int = MAX_DECOMPRESSED_BYTES
    line_bytes: int = MAX_LINE_BYTES
    source_records: int = MAX_SOURCE_RECORDS
    output_rows: int = MAX_OUTPUT_ROWS


@dataclass(frozen=True)
class DisplayDerivative:
    """Indexed display bytes plus the source summary stored beside them."""

    data: bytes
    index: bytes
    provenance: dict[str, object]


def build_display_derivative(
    *,
    kind: str,
    source: bytes,
    source_artifact_id: str,
    source_revision: str,
    source_generation: str,
    source_output_type: str,
    limits: DerivativeLimits | None = None,
) -> DisplayDerivative:
    """Convert PET or interval text into a tabix-indexed display file."""

    active = limits or DerivativeLimits()
    if len(source) > active.source_bytes:
        raise DisplayDerivativeError("source_too_large")
    if kind in {"pet_records", "paired_intervals"}:
        rows, pair_count = _pet_rows(source, active)
        score_meaning = "display_flag"
        detail = (
            "Two longrange rows are anchors for one source pair. Their "
            "scores are the flag 1 and must not be added. MAPQ and any "
            "source score column were not copied."
        )
    elif kind in {"peak_intervals", "annotation"}:
        rows = _interval_rows(
            source, copy_score=kind == "peak_intervals", limits=active
        )
        pair_count = len(rows)
        score_meaning = (
            "source_interval_score" if kind == "peak_intervals" else "display_flag"
        )
        detail = "Indexed intervals use BED 0-based half-open coordinates."
    else:
        raise DisplayDerivativeError("derivative_kind_unsupported")
    if not rows:
        raise DisplayDerivativeError("source_has_no_records")
    rows.sort(key=_sort_key)
    from encode_pipeline.services.browser_bgzf import write_indexed_bed

    indexed = write_indexed_bed(tuple(rows))
    data_sha = hashlib.sha256(indexed.bgzf).hexdigest()
    index_sha = hashlib.sha256(indexed.tabix).hexdigest()
    first = rows[0].split("\t")
    suggested = _suggested_region(first[0], int(first[1]), int(first[2]))
    return DisplayDerivative(
        data=indexed.bgzf,
        index=indexed.tabix,
        provenance={
            "converter_id": CONVERTER_ID,
            "kind": kind,
            "source_artifact_id": source_artifact_id,
            "source_revision": source_revision,
            "source_generation": source_generation,
            "source_output_type": source_output_type,
            "coordinate_convention": "bed-0-based-half-open",
            "score_meaning": score_meaning,
            "score_policy": detail,
            "pair_count": pair_count,
            "anchor_row_count": indexed.record_count,
            "data_sha256": data_sha,
            "index_sha256": index_sha,
            "data_size": len(indexed.bgzf),
            "index_size": len(indexed.tabix),
            "suggested_region": suggested,
        },
    )


def _gzip_chunks(source: bytes, limits: DerivativeLimits):
    decompressor = zlib.decompressobj(_GZIP_WBITS)
    pending = source
    produced = 0
    while True:
        budget = limits.decompressed_bytes - produced
        try:
            chunk = decompressor.decompress(pending, budget + 1)
        except zlib.error as exc:
            raise DisplayDerivativeError("source_unreadable") from exc
        pending = decompressor.unconsumed_tail
        if len(chunk) > budget:
            raise DisplayDerivativeError("source_too_large")
        produced += len(chunk)
        if chunk:
            yield chunk
        if decompressor.eof:
            if pending:
                raise DisplayDerivativeError("source_unreadable")
            return
        if pending:
            continue
        try:
            flushed = decompressor.flush(budget + 1)
        except zlib.error as exc:
            raise DisplayDerivativeError("source_unreadable") from exc
        if len(flushed) > budget:
            raise DisplayDerivativeError("source_too_large")
        if flushed:
            yield flushed
        if not decompressor.eof:
            raise DisplayDerivativeError("source_unreadable")
        return


def _iter_text_lines(source: bytes, limits: DerivativeLimits):
    if source.startswith(b"\x1f\x8b"):
        chunks = _gzip_chunks(source, limits)
    else:
        if len(source) > limits.decompressed_bytes:
            raise DisplayDerivativeError("source_too_large")
        chunks = (source,)
    pending = bytearray()
    for chunk in chunks:
        pending.extend(chunk)
        while True:
            newline = pending.find(b"\n")
            if newline < 0:
                break
            if newline > limits.line_bytes:
                raise DisplayDerivativeError("source_too_large")
            raw = bytes(pending[:newline])
            del pending[: newline + 1]
            yield _text_line(raw)
        if len(pending) > limits.line_bytes:
            raise DisplayDerivativeError("source_too_large")
    if pending:
        yield _text_line(bytes(pending))


def _text_line(raw: bytes) -> str:
    if raw.endswith(b"\r"):
        raw = raw[:-1]
    if b"\x00" in raw or b"\r" in raw:
        raise DisplayDerivativeError("source_unreadable")
    try:
        return raw.decode("ascii")
    except UnicodeError as exc:
        raise DisplayDerivativeError("source_unreadable") from exc


def _reserve_record(count: int, limits: DerivativeLimits) -> int:
    if count >= limits.source_records:
        raise DisplayDerivativeError("source_too_large")
    return count + 1


def _append_rows(
    rows: list[str], additions: tuple[str, ...], limits: DerivativeLimits
) -> None:
    if len(rows) + len(additions) > limits.output_rows:
        raise DisplayDerivativeError("source_too_large")
    rows.extend(additions)


def _pet_rows(source: bytes, limits: DerivativeLimits) -> tuple[list[str], int]:
    rows: list[str] = []
    pairs = 0
    for raw in _iter_text_lines(source, limits):
        if not raw or raw.startswith("#"):
            continue
        columns = raw.split("\t")
        if len(columns) < 6:
            raise DisplayDerivativeError("source_unreadable")
        left = _interval(columns[0], columns[1], columns[2])
        right = _interval(columns[3], columns[4], columns[5])
        if len(columns) >= 10:
            _reject_mapq_as_score(columns[7])
        pairs = _reserve_record(pairs, limits)
        anchors = [_anchor(left, right)]
        if left != right:
            anchors.append(_anchor(right, left))
        _append_rows(rows, tuple(anchors), limits)
    return rows, pairs


def _interval_rows(
    source: bytes, *, copy_score: bool, limits: DerivativeLimits
) -> list[str]:
    rows: list[str] = []
    records = 0
    for raw in _iter_text_lines(source, limits):
        if not raw or raw.startswith(("#", "track", "browser")):
            continue
        columns = raw.split("\t")
        if len(columns) < 3:
            raise DisplayDerivativeError("source_unreadable")
        chrom, start, end = _interval(columns[0], columns[1], columns[2])
        name = columns[3] if len(columns) > 3 and columns[3] else "."
        if any(char in name for char in "\t\r\n"):
            raise DisplayDerivativeError("source_unreadable")
        score = "."
        if copy_score and len(columns) > 4 and _INTEGER.fullmatch(columns[4]):
            score = columns[4]
        strand = (
            columns[5] if len(columns) > 5 and columns[5] in {"+", "-", "."} else "."
        )
        records = _reserve_record(records, limits)
        _append_rows(
            rows, (f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}",), limits
        )
    return rows


def _interval(chrom: str, start_text: str, end_text: str) -> tuple[str, int, int]:
    if _CHROM.fullmatch(chrom) is None:
        raise DisplayDerivativeError("source_unreadable")
    if _INTEGER.fullmatch(start_text) is None or _INTEGER.fullmatch(end_text) is None:
        raise DisplayDerivativeError("source_unreadable")
    start = int(start_text)
    end = int(end_text)
    if end <= start:
        raise DisplayDerivativeError("source_unreadable")
    return chrom, start, end


def _reject_mapq_as_score(value: str) -> None:
    if _INTEGER.fullmatch(value) is None or int(value) > 255:
        raise DisplayDerivativeError("source_unreadable")


def _anchor(
    origin: tuple[str, int, int],
    mate: tuple[str, int, int],
) -> str:
    chrom, start, end = origin
    mate_chrom, mate_start, mate_end = mate
    target = f"{mate_chrom}:{mate_start}-{mate_end},{DISPLAY_FLAG}"
    return f"{chrom}\t{start}\t{end}\t{target}"


def _sort_key(line: str) -> tuple[str, int, int]:
    chrom, start, end, *_rest = line.split("\t")
    return chrom, int(start), int(end)


def _suggested_region(chrom: str, start: int, end: int) -> str:
    flank = max(1000, end - start)
    view_start = max(0, start - flank)
    view_end = end + flank
    return f"{chrom}:{view_start}-{view_end}"
