"""Synthetic local-track fixtures and a real windowed byte count.

The longrange file is a prebuilt BGZF+tabix track. Filler records push it past
the online converter caps (32 MiB source or 2,000,000 records). A window query
reads the tabix file and only the BGZF blocks that overlap that window.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path
import struct
import sys
import time
import zlib

_ROOT = Path(__file__).resolve().parents[2]
_SRC = _ROOT / "src"
if str(_SRC) not in sys.path:
    sys.path.insert(0, str(_SRC))

from encode_pipeline.services.browser_bgzf import (  # noqa: E402
    BrowserIndexError,
    _BGZF_EOF,
    _LINEAR_SHIFT,
    _UNSET_OFFSET,
    _bgzf_block,
    _bgzf_compress,
    _decompress_block,
    _read_tabix,
    _reg2bin,
    _reg2bins,
)
from encode_pipeline.services.browser_derivatives import (  # noqa: E402
    MAX_SOURCE_BYTES,
    MAX_SOURCE_RECORDS,
)

CONSENSUS_PATH = (
    _ROOT / "test" / "fixtures" / "browser" / "local-overlay-consensus.json"
)
DEFAULT_FILLER_RECORDS = 2_000_001
_BLOCK_LIMIT = 16384
_BIGWIG_MAGIC = 0x888FFC26
_BPT_MAGIC = 0x78CA8C91
_CIR_MAGIC = 0x2468ACE0


def load_consensus() -> dict:
    return json.loads(CONSENSUS_PATH.read_text(encoding="utf-8"))


def anchor_line(
    chrom: str, start: int, end: int, mate_start: int, mate_end: int
) -> str:
    return f"{chrom}\t{start}\t{end}\t{chrom}:{mate_start}-{mate_end},1"


def pair_lines(chrom: str, pair: dict) -> tuple[str, str]:
    a_start, a_end = pair["a"]
    b_start, b_end = pair["b"]
    return (
        anchor_line(chrom, a_start, a_end, b_start, b_end),
        anchor_line(chrom, b_start, b_end, a_start, a_end),
    )


def peak_line(peak: dict) -> str:
    return (
        f"{peak['chrom']}\t{peak['start']}\t{peak['end']}\t"
        f"{peak['name']}\t{peak['score']}\t+"
    )


def region_string(window: dict) -> str:
    return f"{window['chrom']}:{window['start']}-{window['end']}"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _filler_line(index: int) -> bytes:
    start = 1_000_000 + index * 20
    end = start + 10
    mate = end
    mate_end = mate + 10
    mixed = (index + 1) * 0x9E3779B1
    token = "".join(f"{(mixed * (step + 3)) & 0xFFFFFFFF:08x}" for step in range(8))
    return (f"chr1\t{start}\t{end}\tchr1:{mate}-{mate_end},1\t{token}").encode("ascii")


class _StreamedLongrange:
    """Write coordinate-sorted BED text as BGZF and a tabix index."""

    def __init__(self, data_path: Path) -> None:
        self._handle = data_path.open("wb")
        self.block_offset = 0
        self.buffer = bytearray()
        self.count = 0
        self._order: list[str] = []
        self._bins: dict[str, dict[int, list[list[int]]]] = {}
        self._linear: dict[str, list[int]] = {}
        self._last_chrom = ""
        self._last_start = -1
        self._last_chunk: list[int] | None = None

    def add_text(self, line: bytes) -> None:
        if b"\n" in line or b"\r" in line or not line.isascii():
            raise BrowserIndexError("indexed line is invalid")
        parts = line.split(b"\t")
        if len(parts) < 3:
            raise BrowserIndexError("indexed line is invalid")
        chrom = parts[0].decode("ascii")
        start = int(parts[1])
        end = int(parts[2])
        if start < 0 or end <= start:
            raise BrowserIndexError("indexed coordinate is invalid")
        if chrom != self._last_chrom:
            if chrom in self._bins:
                raise BrowserIndexError("indexed contig is out of order")
            self._order.append(chrom)
            self._bins[chrom] = {}
            self._linear[chrom] = []
            self._last_chrom = chrom
            self._last_start = -1
        elif start < self._last_start:
            raise BrowserIndexError("indexed lines are not sorted")
        text = line + b"\n"
        if len(text) > _BLOCK_LIMIT:
            raise BrowserIndexError("indexed line is too long")
        if self.buffer and len(self.buffer) + len(text) > _BLOCK_LIMIT:
            self._flush()
        virtual_start = (self.block_offset << 16) | len(self.buffer)
        self.buffer.extend(text)
        virtual_end = (self.block_offset << 16) | len(self.buffer)
        self._note(chrom, start, end, virtual_start, virtual_end)
        self._last_start = start
        self.count += 1

    def flush(self) -> None:
        self._flush()

    def finish(self) -> int:
        self._flush()
        if self._last_chunk is not None:
            self._last_chunk[1] = self.block_offset << 16
        self._handle.write(_BGZF_EOF)
        self._handle.close()
        return self.count

    def tabix_bytes(self) -> bytes:
        names = b"".join(chrom.encode("ascii") + b"\0" for chrom in self._order)
        pieces = [
            b"TBI\x01",
            struct.pack("<i", len(self._order)),
            struct.pack("<6i", 0x10000, 1, 2, 3, ord("#"), 0),
            struct.pack("<i", len(names)),
            names,
        ]
        for chrom in self._order:
            pieces.append(self._reference_index(chrom))
        return _bgzf_compress(b"".join(pieces))

    def _flush(self) -> None:
        if not self.buffer:
            return
        block = _bgzf_block(bytes(self.buffer))
        self._handle.write(block)
        self.block_offset += len(block)
        self.buffer = bytearray()

    def _note(
        self,
        chrom: str,
        start: int,
        end: int,
        virtual_start: int,
        virtual_end: int,
    ) -> None:
        needed = (max(start, end - 1) >> _LINEAR_SHIFT) + 1
        linear = self._linear[chrom]
        if len(linear) < needed:
            linear.extend([_UNSET_OFFSET] * (needed - len(linear)))
        last = max(start, end - 1) >> _LINEAR_SHIFT
        for index in range(start >> _LINEAR_SHIFT, last + 1):
            if virtual_start < linear[index]:
                linear[index] = virtual_start
        bin_id = _reg2bin(start, end)
        chunks = self._bins[chrom].setdefault(bin_id, [])
        if chunks and chunks[-1][1] == virtual_start:
            chunks[-1][1] = virtual_end
            self._last_chunk = chunks[-1]
        else:
            chunks.append([virtual_start, virtual_end])
            self._last_chunk = chunks[-1]

    def _reference_index(self, chrom: str) -> bytes:
        bins = self._bins[chrom]
        linear = self._linear[chrom]
        body = [struct.pack("<I", len(bins))]
        for bin_id in sorted(bins):
            chunks = bins[bin_id]
            body.append(struct.pack("<II", bin_id, len(chunks)))
            for start, end in chunks:
                body.append(struct.pack("<QQ", start, end))
        body.append(struct.pack("<I", len(linear)))
        for value in linear:
            stored = 0 if value == _UNSET_OFFSET else value
            body.append(struct.pack("<Q", stored))
        return b"".join(body)


class _CountedReader:
    def __init__(self, path: Path) -> None:
        self._handle = path.open("rb")
        self.bytes_read = 0

    def read_all(self) -> bytes:
        data = self._handle.read()
        self.bytes_read += len(data)
        return data

    def read_at(self, offset: int, size: int) -> bytes:
        self._handle.seek(offset)
        data = self._handle.read(size)
        self.bytes_read += len(data)
        return data

    def close(self) -> None:
        self._handle.close()


def _read_block(
    reader: _CountedReader, cache: dict[int, tuple[bytes, int]], offset: int
):
    cached = cache.get(offset)
    if cached is not None:
        return cached
    header = reader.read_at(offset, 18)
    if len(header) < 18 or header[:2] != b"\x1f\x8b":
        raise BrowserIndexError("bgzf block is invalid")
    block_size = struct.unpack_from("<H", header, 16)[0] + 1
    rest = reader.read_at(offset + 18, block_size - 18)
    raw, size = _decompress_block(header + rest, 0)
    if size != block_size:
        raise BrowserIndexError("bgzf block size mismatch")
    cache[offset] = (raw, size)
    return cache[offset]


def _read_virtual_file(
    reader: _CountedReader,
    cache: dict[int, tuple[bytes, int]],
    start: int,
    end: int,
) -> bytes:
    pieces: list[bytes] = []
    cursor = start
    while cursor < end:
        block_offset = cursor >> 16
        uncompressed_offset = cursor & 0xFFFF
        raw, block_size = _read_block(reader, cache, block_offset)
        if (end >> 16) == block_offset:
            pieces.append(raw[uncompressed_offset : end & 0xFFFF])
            break
        pieces.append(raw[uncompressed_offset:])
        cursor = (block_offset + block_size) << 16
    return b"".join(pieces)


@dataclass(frozen=True)
class WindowRead:
    region: str
    hits: tuple[str, ...]
    data_bytes_read: int
    index_bytes_read: int
    elapsed_ms: float


def measure_window(
    data_path: Path, index_path: Path, chrom: str, start: int, end: int
) -> WindowRead:
    """Query one window and count bytes actually read from each file."""

    started = time.perf_counter()
    index_reader = _CountedReader(index_path)
    try:
        references = _read_tabix(index_reader.read_all())
        index_bytes = index_reader.bytes_read
    finally:
        index_reader.close()
    hits: list[str] = []
    data_bytes = 0
    wanted = references.get(chrom)
    if wanted is not None:
        bin_ids = set(_reg2bins(start, end))
        chunks: list[tuple[int, int]] = []
        for bin_id, bin_chunks in wanted:
            if bin_id in bin_ids:
                chunks.extend(bin_chunks)
        reader = _CountedReader(data_path)
        cache: dict[int, tuple[bytes, int]] = {}
        seen: set[str] = set()
        try:
            for virtual_start, virtual_end in sorted(chunks):
                raw = _read_virtual_file(reader, cache, virtual_start, virtual_end)
                for encoded in raw.splitlines():
                    text = encoded.decode("ascii")
                    columns = text.split("\t")
                    if len(columns) < 3 or columns[0] != chrom:
                        continue
                    if (
                        int(columns[1]) < end
                        and int(columns[2]) > start
                        and text not in seen
                    ):
                        seen.add(text)
                        hits.append(text)
        finally:
            data_bytes = reader.bytes_read
            reader.close()
    elapsed_ms = (time.perf_counter() - started) * 1000
    return WindowRead(
        region=f"{chrom}:{start}-{end}",
        hits=tuple(hits),
        data_bytes_read=data_bytes,
        index_bytes_read=index_bytes,
        elapsed_ms=elapsed_ms,
    )


def write_signal_bigwig(
    path: Path,
    chrom_sizes: list[tuple[str, int]],
    intervals: list[tuple[str, int, int, float]],
) -> None:
    """Write an uncompressed-section bigWig that pyBigWig and WashU can read."""

    chrom_ids = {name: index for index, (name, _size) in enumerate(chrom_sizes)}
    grouped: dict[str, list[tuple[int, int, float]]] = {
        name: [] for name, _size in chrom_sizes
    }
    for chrom, start, end, value in intervals:
        grouped[chrom].append((start, end, value))
    sections: list[tuple[int, int, int, bytes, bytes]] = []
    covered = 0
    minimum = None
    maximum = None
    total = 0.0
    squares = 0.0
    for chrom, rows in grouped.items():
        if not rows:
            continue
        rows.sort()
        raw = _bedgraph_section(chrom_ids[chrom], rows)
        compressed = zlib.compress(raw)
        section_start = rows[0][0]
        section_end = rows[-1][1]
        sections.append((chrom_ids[chrom], section_start, section_end, raw, compressed))
        for start, end, value in rows:
            bases = end - start
            covered += bases
            total += value * bases
            squares += value * value * bases
            minimum = value if minimum is None else min(minimum, value)
            maximum = value if maximum is None else max(maximum, value)
    if not sections or minimum is None or maximum is None:
        raise BrowserIndexError("signal bigWig has no intervals")
    summary = struct.pack("<Qdddd", covered, minimum, maximum, total, squares)
    key_size = max(len(name) for name, _size in chrom_sizes)
    chrom_tree = _chrom_tree(chrom_sizes, key_size)
    prefix = struct.pack("<II", len(sections), 0)
    data = bytearray(prefix)
    index_items: list[tuple[int, int, int, int, int]] = []
    data_offset = 64 + len(summary) + len(chrom_tree)
    cursor = data_offset + len(prefix)
    for chrom_id, start, end, _raw, compressed in sections:
        index_items.append((chrom_id, start, end, cursor, len(compressed)))
        data.extend(compressed)
        cursor += len(compressed)
    index = _cir_tree(index_items, len(data))
    header = bytearray(64)
    struct.pack_into("<IHH", header, 0, _BIGWIG_MAGIC, 4, 0)
    struct.pack_into("<Q", header, 8, 64 + len(summary))
    struct.pack_into("<Q", header, 16, data_offset)
    struct.pack_into("<Q", header, 24, data_offset + len(data))
    struct.pack_into("<QQ", header, 36, 0, 64)
    struct.pack_into("<I", header, 52, 32768)
    path.write_bytes(
        bytes(header)
        + summary
        + chrom_tree
        + bytes(data)
        + index
        + struct.pack("<I", _BIGWIG_MAGIC)
    )


def read_signal_intervals(path: Path) -> list[tuple[str, int, int, float]]:
    """Return bedGraph intervals from a bigWig written by this module."""

    data = path.read_bytes()
    magic, version, zoom = struct.unpack_from("<IHH", data, 0)
    if magic != _BIGWIG_MAGIC or version != 4 or zoom != 0:
        raise BrowserIndexError("signal bigWig header is invalid")
    chrom_offset, data_offset, index_offset = struct.unpack_from("<QQQ", data, 8)
    names = _read_chrom_names(data, chrom_offset)
    (
        _magic,
        _block,
        item_count,
        start_chrom,
        _start,
        end_chrom,
        _end,
        _file_size,
        _slot,
        _reserved,
    ) = struct.unpack_from("<IIQIIIIQII", data, index_offset)
    if _magic != _CIR_MAGIC:
        raise BrowserIndexError("signal bigWig index is invalid")
    cursor = index_offset + 48
    is_leaf, _reserved_byte, count = struct.unpack_from("<BBH", data, cursor)
    if is_leaf != 1 or count != item_count:
        raise BrowserIndexError("signal bigWig index is invalid")
    cursor += 4
    found: list[tuple[str, int, int, float]] = []
    for _item in range(count):
        chrom_id, _start, _end_chrom, _item_end, offset, size = struct.unpack_from(
            "<IIIIQQ", data, cursor
        )
        cursor += 32
        raw = zlib.decompress(data[offset : offset + size])
        found.extend(_decode_section(names[chrom_id], raw))
    if end_chrom >= len(names) or start_chrom < 0:
        raise BrowserIndexError("signal bigWig index is invalid")
    if data_offset > len(data):
        raise BrowserIndexError("signal bigWig data is invalid")
    return found


def _bedgraph_section(chrom_id: int, rows: list[tuple[int, int, float]]) -> bytes:
    header = struct.pack("<IIIII", chrom_id, rows[0][0], rows[-1][1], 0, 0)
    header += struct.pack("<BBH", 1, 0, len(rows))
    body = b"".join(
        struct.pack("<IIf", start, end, value) for start, end, value in rows
    )
    return header + body


def _decode_section(chrom: str, raw: bytes) -> list[tuple[str, int, int, float]]:
    _chrom, _start, _end, _step, _reserved, kind, _pad, count = struct.unpack_from(
        "<IIIII BBH", raw, 0
    )
    if kind != 1:
        raise BrowserIndexError("signal section is not bedGraph")
    cursor = 24
    rows: list[tuple[str, int, int, float]] = []
    for _index in range(count):
        start, end, value = struct.unpack_from("<IIf", raw, cursor)
        rows.append((chrom, start, end, value))
        cursor += 12
    return rows


def _chrom_tree(chrom_sizes: list[tuple[str, int]], key_size: int) -> bytes:
    header = struct.pack(
        "<IIIIQQ", _BPT_MAGIC, len(chrom_sizes), key_size, 8, len(chrom_sizes), 0
    )
    leaf = struct.pack("<BBH", 1, 0, len(chrom_sizes))
    items = bytearray()
    for index, (name, size) in enumerate(chrom_sizes):
        items.extend(name.encode("ascii").ljust(key_size, b"\0"))
        items.extend(struct.pack("<II", index, size))
    return header + leaf + bytes(items)


def _read_chrom_names(data: bytes, offset: int) -> list[str]:
    magic, _block, key_size, val_size, count, _reserved = struct.unpack_from(
        "<IIIIQQ", data, offset
    )
    if magic != _BPT_MAGIC or val_size != 8:
        raise BrowserIndexError("signal bigWig chromosomes are invalid")
    cursor = offset + 32
    is_leaf, _pad, leaf_count = struct.unpack_from("<BBH", data, cursor)
    if is_leaf != 1 or leaf_count != count:
        raise BrowserIndexError("signal bigWig chromosomes are invalid")
    cursor += 4
    names: list[str] = []
    for _index in range(count):
        names.append(
            data[cursor : cursor + key_size].split(b"\0", 1)[0].decode("ascii")
        )
        cursor += key_size + val_size
    return names


def _cir_tree(items: list[tuple[int, int, int, int, int]], data_size: int) -> bytes:
    header = struct.pack(
        "<IIQIIIIQII",
        _CIR_MAGIC,
        64,
        len(items),
        items[0][0],
        items[0][1],
        items[-1][0],
        items[-1][2],
        data_size,
        1,
        0,
    )
    leaf = struct.pack("<BBH", 1, 0, len(items))
    body = bytearray()
    for chrom_id, start, end, offset, size in items:
        body.extend(
            struct.pack("<IIIIQQ", chrom_id, start, chrom_id, end, offset, size)
        )
    return header + leaf + bytes(body)


def _feature_groups(consensus: dict) -> list[list[bytes]]:
    groups: list[list[bytes]] = []
    for key, chrom in (("egfr", "chr7"), ("tp53", "chr17")):
        for pair in consensus["longrange_pairs"][key]:
            groups.append([line.encode("ascii") for line in pair_lines(chrom, pair)])
    return groups


def _link_or_copy(source: Path, dest: Path) -> None:
    if dest.exists():
        dest.unlink()
    try:
        os.link(source, dest)
    except OSError:
        dest.write_bytes(source.read_bytes())


def build_local_overlay_fixtures(
    dest: Path,
    *,
    filler_records: int = DEFAULT_FILLER_RECORDS,
    consensus: dict | None = None,
) -> dict:
    """Write the synthetic files and the measured window read."""

    started = time.perf_counter()
    spec = consensus if consensus is not None else load_consensus()
    dest.mkdir(parents=True, exist_ok=True)
    names = spec["filenames"]
    data_path = dest / names["longrange"]
    index_path = dest / names["longrange_index"]
    stream = _StreamedLongrange(data_path)
    for index in range(filler_records):
        stream.add_text(_filler_line(index))
    # Each planted pair gets its own BGZF block so a window cannot
    # accidentally decompress a pair that sits elsewhere on the chromosome.
    for group in _feature_groups(spec):
        stream.flush()
        for line in group:
            stream.add_text(line)
    record_count = stream.finish()
    index_path.write_bytes(stream.tabix_bytes())
    _link_or_copy(data_path, dest / names["longrange_missing_index"])
    _link_or_copy(data_path, dest / names["longrange_wrong"])
    wrong = _StreamedLongrange(dest / "wrong-source.longrange.gz")
    wrong.add_text(b"chr1\t100\t200\tchr1:300-400,1")
    wrong.finish()
    (dest / names["longrange_wrong_index"]).write_bytes(wrong.tabix_bytes())
    (dest / "wrong-source.longrange.gz").unlink()
    chrom_sizes = list(spec["chrom_sizes_hg38"].items())
    signal_intervals = [
        (item["chrom"], item["start"], item["end"], float(item["value"]))
        for item in spec["signals"]
        if item["chrom"] in spec["chrom_sizes_hg38"]
    ]
    # Header chromosomes follow the order WashU expects for the chromosomes we use.
    used = []
    for name, size in chrom_sizes:
        if any(item[0] == name for item in signal_intervals):
            used.append((name, size))
    write_signal_bigwig(dest / names["signal"], used, signal_intervals)
    peak_text = "".join(peak_line(peak) + "\n" for peak in spec["peaks"])
    (dest / names["peaks"]).write_text(peak_text, encoding="ascii")
    window = spec["windows"]["egfr"]
    measured = measure_window(
        data_path, index_path, window["chrom"], window["start"], window["end"]
    )
    missed = measure_window(data_path, index_path, "chr7", 1000, 2000)
    observation = {
        "filler_records": filler_records,
        "record_count": record_count,
        "anchor_records": record_count - filler_records,
        "data_file_bytes": data_path.stat().st_size,
        "index_file_bytes": index_path.stat().st_size,
        "data_sha256": _sha256(data_path),
        "index_sha256": _sha256(index_path),
        "signal_file_bytes": (dest / names["signal"]).stat().st_size,
        "peaks_file_bytes": (dest / names["peaks"]).stat().st_size,
        "online_source_byte_limit": MAX_SOURCE_BYTES,
        "online_source_record_limit": MAX_SOURCE_RECORDS,
        "exceeds_source_bytes": data_path.stat().st_size > MAX_SOURCE_BYTES,
        "exceeds_source_records": record_count > MAX_SOURCE_RECORDS,
        "window": {
            "region": measured.region,
            "data_bytes_read": measured.data_bytes_read,
            "index_bytes_read": measured.index_bytes_read,
            "elapsed_ms": round(measured.elapsed_ms, 3),
            "hits": list(measured.hits),
        },
        "miss": {
            "region": missed.region,
            "data_bytes_read": missed.data_bytes_read,
            "index_bytes_read": missed.index_bytes_read,
            "hits": list(missed.hits),
        },
        "build_seconds": round(time.perf_counter() - started, 3),
    }
    (dest / "observation.json").write_text(
        json.dumps(observation, indent=2) + "\n", encoding="utf-8"
    )
    return observation


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dest", type=Path, required=True)
    parser.add_argument("--filler-records", type=int, default=DEFAULT_FILLER_RECORDS)
    args = parser.parse_args()
    observation = build_local_overlay_fixtures(
        args.dest, filler_records=args.filler_records
    )
    window = observation["window"]
    if not window["hits"] or window["data_bytes_read"] <= 0:
        raise SystemExit("window query did not read the planted anchors")
    if window["data_bytes_read"] >= observation["data_file_bytes"]:
        raise SystemExit("window query read the whole longrange file")
    if (
        not observation["exceeds_source_records"]
        and not observation["exceeds_source_bytes"]
    ):
        raise SystemExit("synthetic longrange does not exceed the online converter cap")
    print(
        "LOCAL_OVERLAY_FIXTURE "
        + json.dumps(
            {
                "records": observation["record_count"],
                "data_file_bytes": observation["data_file_bytes"],
                "index_file_bytes": observation["index_file_bytes"],
                "window_data_bytes_read": window["data_bytes_read"],
                "window_index_bytes_read": window["index_bytes_read"],
                "window_elapsed_ms": window["elapsed_ms"],
                "window_hits": len(window["hits"]),
                "miss_data_bytes_read": observation["miss"]["data_bytes_read"],
            },
            separators=(",", ":"),
        )
    )


if __name__ == "__main__":
    main()
