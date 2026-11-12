"""BGZF and tabix for browser display derivatives.

The writer follows the BGZF EOF marker and TBI layout used by htslib. It is
limited to coordinate-sorted BED-like text so PET and interval tracks can be
queried by genomic window.
"""

from __future__ import annotations

from dataclasses import dataclass
import struct
import zlib


_BGZF_EOF = bytes(
    (
        0x1F,
        0x8B,
        0x08,
        0x04,
        0x00,
        0x00,
        0x00,
        0x00,
        0x00,
        0xFF,
        0x06,
        0x00,
        0x42,
        0x43,
        0x02,
        0x00,
        0x1B,
        0x00,
        0x03,
        0x00,
        0x00,
        0x00,
        0x00,
        0x00,
        0x00,
        0x00,
        0x00,
        0x00,
    )
)
_MAX_UNCOMPRESSED = 65280
_UNSET_OFFSET = (1 << 64) - 1
_LINEAR_SHIFT = 14


class BrowserIndexError(ValueError):
    """A display file could not be indexed."""


@dataclass(frozen=True)
class IndexedText:
    """BGZF payload, TBI payload, and the records that were indexed."""

    bgzf: bytes
    tabix: bytes
    record_count: int


@dataclass(frozen=True)
class _LineSpan:
    chrom: str
    start: int
    end: int
    virtual_start: int
    virtual_end: int
    text: bytes


def write_indexed_bed(lines: tuple[str, ...]) -> IndexedText:
    """BGZIP and tabix-index already sorted BED-like lines."""

    encoded = tuple(_encode_line(line) for line in lines)
    bgzf, spans = _pack_bgzf(encoded)
    tabix = _bgzf_compress(_write_tabix(spans))
    return IndexedText(bgzf=bgzf + _BGZF_EOF, tabix=tabix, record_count=len(spans))


def query_indexed_bed(
    indexed: IndexedText,
    chrom: str,
    start: int,
    end: int,
    *,
    limit: int,
) -> tuple[tuple[str, ...], bool]:
    """Return overlapping lines and whether the limit stopped the scan."""

    if limit < 1:
        raise BrowserIndexError("query limit is invalid")
    if start < 0 or end <= start:
        raise BrowserIndexError("query interval is invalid")
    matches: list[str] = []
    truncated = False
    for span in _overlapping_spans(indexed, chrom, start, end):
        if len(matches) >= limit:
            truncated = True
            break
        if span.start < end and span.end > start:
            matches.append(span.text.decode("ascii").rstrip("\n"))
    return tuple(matches), truncated


def _encode_line(line: str) -> _LineSpan:
    if not isinstance(line, str) or "\n" in line or "\r" in line:
        raise BrowserIndexError("indexed line is invalid")
    parts = line.split("\t")
    if len(parts) < 3:
        raise BrowserIndexError("indexed line is invalid")
    chrom = parts[0]
    if not chrom.isascii() or not chrom or any(char.isspace() for char in chrom):
        raise BrowserIndexError("indexed contig is invalid")
    try:
        start = int(parts[1])
        end = int(parts[2])
    except ValueError as exc:
        raise BrowserIndexError("indexed coordinate is invalid") from exc
    if start < 0 or end <= start:
        raise BrowserIndexError("indexed coordinate is invalid")
    return _LineSpan(chrom, start, end, 0, 0, (line + "\n").encode("ascii"))


def _pack_bgzf(
    lines: tuple[_LineSpan, ...],
) -> tuple[bytes, tuple[_LineSpan, ...]]:
    blocks: list[bytes] = []
    spans: list[_LineSpan] = []
    buffer = bytearray()
    block_offset = 0

    def flush() -> None:
        nonlocal buffer, block_offset
        if not buffer:
            return
        block = _bgzf_block(bytes(buffer))
        blocks.append(block)
        block_offset += len(block)
        buffer = bytearray()

    for line in lines:
        if len(line.text) > _MAX_UNCOMPRESSED:
            raise BrowserIndexError("indexed line is too long")
        if buffer and len(buffer) + len(line.text) > _MAX_UNCOMPRESSED:
            flush()
        virtual_start = (block_offset << 16) | len(buffer)
        buffer.extend(line.text)
        virtual_end = (block_offset << 16) | len(buffer)
        spans.append(
            _LineSpan(
                line.chrom,
                line.start,
                line.end,
                virtual_start,
                virtual_end,
                line.text,
            )
        )
    flush()
    payload = b"".join(blocks)
    if spans:
        eof_virtual = len(payload) << 16
        last = spans[-1]
        spans[-1] = _LineSpan(
            last.chrom,
            last.start,
            last.end,
            last.virtual_start,
            eof_virtual,
            last.text,
        )
    return payload, tuple(spans)


def _bgzf_compress(data: bytes) -> bytes:
    pieces = [
        _bgzf_block(data[offset : offset + _MAX_UNCOMPRESSED])
        for offset in range(0, max(len(data), 1), _MAX_UNCOMPRESSED)
    ]
    if not data:
        pieces = [_bgzf_block(b"")]
    return b"".join(pieces) + _BGZF_EOF


def _bgzf_block(data: bytes) -> bytes:
    compressor = zlib.compressobj(6, zlib.DEFLATED, -15)
    compressed = compressor.compress(data) + compressor.flush()
    total = len(compressed) + 26
    if total > 65536:
        raise BrowserIndexError("compressed block is too large")
    header = struct.pack(
        "<BBBBIBBH",
        0x1F,
        0x8B,
        8,
        4,
        0,
        0,
        255,
        6,
    )
    extra = struct.pack("<BBHH", 66, 67, 2, total - 1)
    trailer = struct.pack("<II", zlib.crc32(data) & 0xFFFFFFFF, len(data))
    return header + extra + compressed + trailer


def _write_tabix(spans: tuple[_LineSpan, ...]) -> bytes:
    grouped: dict[str, list[_LineSpan]] = {}
    order: list[str] = []
    for span in spans:
        if span.chrom not in grouped:
            order.append(span.chrom)
            grouped[span.chrom] = []
        grouped[span.chrom].append(span)
    names = b"".join(chrom.encode("ascii") + b"\0" for chrom in order)
    pieces = [
        b"TBI\x01",
        struct.pack("<i", len(order)),
        struct.pack("<6i", 0x10000, 1, 2, 3, ord("#"), 0),
        struct.pack("<i", len(names)),
        names,
    ]
    for chrom in order:
        pieces.append(_reference_index(grouped[chrom]))
    return b"".join(pieces)


def _reference_index(spans: list[_LineSpan]) -> bytes:
    bins: dict[int, list[tuple[int, int]]] = {}
    max_end = 0
    for span in spans:
        max_end = max(max_end, span.end)
        bins.setdefault(_reg2bin(span.start, span.end), []).append(
            (span.virtual_start, span.virtual_end)
        )
    n_intv = (max_end >> _LINEAR_SHIFT) + 1
    linear = [_UNSET_OFFSET] * n_intv
    for span in spans:
        last = max(span.start, span.end - 1) >> _LINEAR_SHIFT
        for index in range(span.start >> _LINEAR_SHIFT, last + 1):
            if span.virtual_start < linear[index]:
                linear[index] = span.virtual_start
    body = [struct.pack("<I", len(bins))]
    for bin_id in sorted(bins):
        chunks = sorted(bins[bin_id])
        body.append(struct.pack("<II", bin_id, len(chunks)))
        for start, end in chunks:
            body.append(struct.pack("<QQ", start, end))
    body.append(struct.pack("<I", n_intv))
    for value in linear:
        stored = 0 if value == _UNSET_OFFSET else value
        body.append(struct.pack("<Q", stored))
    return b"".join(body)


def _reg2bin(start: int, end: int) -> int:
    end -= 1
    if start >> 14 == end >> 14:
        return ((1 << 15) - 1) // 7 + (start >> 14)
    if start >> 17 == end >> 17:
        return ((1 << 12) - 1) // 7 + (start >> 17)
    if start >> 20 == end >> 20:
        return ((1 << 9) - 1) // 7 + (start >> 20)
    if start >> 23 == end >> 23:
        return ((1 << 6) - 1) // 7 + (start >> 23)
    if start >> 26 == end >> 26:
        return ((1 << 3) - 1) // 7 + (start >> 26)
    return 0


def _overlapping_spans(
    indexed: IndexedText,
    chrom: str,
    start: int,
    end: int,
) -> tuple[_LineSpan, ...]:
    references = _read_tabix(indexed.tabix)
    wanted = references.get(chrom)
    if wanted is None:
        return ()
    bin_ids = set(_reg2bins(start, end))
    chunks: list[tuple[int, int]] = []
    for bin_id, bin_chunks in wanted:
        if bin_id in bin_ids:
            chunks.extend(bin_chunks)
    found: list[_LineSpan] = []
    for virtual_start, virtual_end in sorted(chunks):
        raw = _read_virtual(indexed.bgzf, virtual_start, virtual_end)
        for line in raw.splitlines():
            parsed = _encode_line(line.decode("ascii"))
            if parsed.chrom == chrom and parsed.start < end and parsed.end > start:
                found.append(parsed)
    return tuple(found)


def _reg2bins(start: int, end: int) -> tuple[int, ...]:
    end -= 1
    bins = [0]
    for shift, base in ((26, 1), (23, 9), (20, 73), (17, 585), (14, 4681)):
        for key in range(start >> shift, (end >> shift) + 1):
            bins.append(base + key)
    return tuple(bins)


def _read_tabix(
    payload: bytes,
) -> dict[str, tuple[tuple[int, tuple[tuple[int, int], ...]], ...]]:
    payload = _decompress_bgzf(payload)
    if not payload.startswith(b"TBI\x01"):
        raise BrowserIndexError("tabix magic is invalid")
    offset = 4
    (n_ref,) = struct.unpack_from("<i", payload, offset)
    offset += 4 + 24
    (name_length,) = struct.unpack_from("<i", payload, offset)
    offset += 4
    names = payload[offset : offset + name_length].split(b"\0")[:-1]
    offset += name_length
    if len(names) != n_ref:
        raise BrowserIndexError("tabix names are invalid")
    references: list[tuple[tuple[int, tuple[tuple[int, int], ...]], ...]] = []
    for _reference in range(n_ref):
        (n_bin,) = struct.unpack_from("<I", payload, offset)
        offset += 4
        bins: list[tuple[int, tuple[tuple[int, int], ...]]] = []
        for _bin in range(n_bin):
            bin_id, n_chunk = struct.unpack_from("<II", payload, offset)
            offset += 8
            chunks = tuple(
                struct.unpack_from("<QQ", payload, offset + index * 16)
                for index in range(n_chunk)
            )
            offset += n_chunk * 16
            bins.append((bin_id, chunks))
        (n_intv,) = struct.unpack_from("<I", payload, offset)
        offset += 4 + n_intv * 8
        references.append(tuple(bins))
    return {
        name.decode("ascii"): reference for name, reference in zip(names, references)
    }


def _read_virtual(payload: bytes, start: int, end: int) -> bytes:
    pieces: list[bytes] = []
    cursor = start
    while cursor < end:
        block_offset = cursor >> 16
        uncompressed_offset = cursor & 0xFFFF
        raw, block_size = _decompress_block(payload, block_offset)
        if (end >> 16) == block_offset:
            pieces.append(raw[uncompressed_offset : end & 0xFFFF])
            break
        pieces.append(raw[uncompressed_offset:])
        cursor = (block_offset + block_size) << 16
    return b"".join(pieces)


def _decompress_bgzf(payload: bytes) -> bytes:
    if payload.startswith(b"TBI\x01"):
        return payload
    offset = 0
    parts: list[bytes] = []
    while offset < len(payload):
        raw, size = _decompress_block(payload, offset)
        if not raw:
            break
        parts.append(raw)
        offset += size
    return b"".join(parts)


def _decompress_block(payload: bytes, offset: int) -> tuple[bytes, int]:
    if payload[offset : offset + 2] != b"\x1f\x8b":
        raise BrowserIndexError("bgzf block is invalid")
    extra_length = struct.unpack_from("<H", payload, offset + 10)[0]
    block_size = struct.unpack_from("<H", payload, offset + 16)[0] + 1
    block = payload[offset : offset + block_size]
    compressed = block[12 + extra_length : -8]
    raw = zlib.decompress(compressed, wbits=-15)
    crc, isize = struct.unpack_from("<II", block, block_size - 8)
    if (zlib.crc32(raw) & 0xFFFFFFFF) != crc or (len(raw) & 0xFFFFFFFF) != isize:
        raise BrowserIndexError("bgzf block is invalid")
    return raw, block_size
