"""Check that a sealed tabix index locates records in the selected WashU track."""

from __future__ import annotations

from functools import lru_cache
import gzip
import io
import math
import re
import struct
import zlib

from encode_pipeline.platform.archived_tracks import ArchiveStoreError, CHUNK_BYTES
from encode_pipeline.services.archive_objects import ArchiveReader
from encode_pipeline.services.browser_bgzf import _reg2bin


_BLOCK = 65536
_MATE = re.compile(
    rb"[\w.]+:(\d+)-(\d+),([+-]?(?:\d+\.?\d*|\.\d+)(?:e[+-]?\d+)?)", re.I
)


class _MemberStream(io.RawIOBase):
    def __init__(self, member: ArchiveReader, size: int):
        self.member, self.size, self.position = member, size, 0

    def readable(self):
        return True

    def readinto(self, target):
        end = min(self.position + len(target), self.size)
        chunk = self.member.read(self.position, end)
        self.position = end
        target[: len(chunk)] = chunk
        return len(chunk)


def _block_size(header: bytes) -> int:
    if header[:4] != b"\x1f\x8b\x08\x04":
        raise ValueError("expected BGZF")
    end = 12 + struct.unpack_from("<H", header, 10)[0]
    pos = 12
    while pos + 4 <= end:
        tag, length = struct.unpack_from("<2sH", header, pos)
        pos += 4
        if pos + length > end:
            break
        if tag == b"BC" and length == 2:
            return struct.unpack_from("<H", header, pos)[0] + 1
        pos += length
    raise ValueError("missing BGZF block size")


def verify_longrange_index(
    data: ArchiveReader, data_size: int, index: ArchiveReader, index_size: int
) -> None:
    """Follow indexed virtual offsets on verified archive bytes before publication."""

    @lru_cache(maxsize=8)
    def block(offset):
        if not 0 <= offset < data_size:
            raise ValueError("offset outside data")
        compressed = data.read(offset, min(offset + _BLOCK, data_size))
        size = _block_size(compressed)
        decoder = zlib.decompressobj(31)
        raw = decoder.decompress(compressed[:size], _BLOCK + 1)
        if (
            not decoder.eof
            or decoder.unused_data
            or len(raw) > _BLOCK
            or size > len(compressed)
        ):
            raise ValueError("invalid BGZF block")
        return raw, offset + size

    def check_chunk(name, bin_id, start, end, zero_based):
        if start >= end or end > data_size << 16:
            raise ValueError("invalid chunk")
        stop, within_end = divmod(end, _BLOCK)
        if stop < data_size:
            tail, after = block(stop)
            if within_end and (
                within_end > len(tail)
                or (
                    tail[within_end - 1] != 10
                    and (
                        within_end != len(tail)
                        or (after < data_size and block(after)[0])
                    )
                )
            ):
                raise ValueError("chunk ends inside a record")
        offset, within = divmod(start, _BLOCK)
        raw, offset = block(offset)
        if within >= len(raw) or (within and raw[within - 1] != 10):
            raise ValueError("chunk starts outside a record")
        record = raw[within:]
        while b"\n" not in record and len(record) <= _BLOCK and offset < data_size:
            raw, offset = block(offset)
            record += raw
        fields = record.split(b"\n", 1)[0].rstrip(b"\r").split(b"\t")
        if len(fields) < 4 or fields[0] != name:
            raise ValueError("wrong reference")
        left, right = int(fields[1]) - int(not zero_based), int(fields[2])
        mate = _MATE.fullmatch(fields[3])
        if (
            not 0 <= left < right <= 1 << 29
            or not mate
            or int(mate[1]) >= int(mate[2])
            or not math.isfinite(float(mate[3]))
        ):
            raise ValueError("invalid longrange record")
        feature = _reg2bin(left, right)
        while feature > bin_id:
            feature = (feature - 1) >> 3
        if feature != bin_id:
            raise ValueError("wrong coordinate bin")

    try:
        with io.BufferedReader(
            _MemberStream(index, index_size), CHUNK_BYTES
        ) as compressed:
            _block_size(compressed.peek(_BLOCK))
            with gzip.GzipFile(fileobj=compressed) as stream:
                _check_tabix(stream, data_size, check_chunk)
    except (EOFError, OSError, ValueError, struct.error, zlib.error, OverflowError):
        raise ArchiveStoreError("archive_invalid") from None


def _check_tabix(stream, data_size, check_chunk):
    magic, refs, flags, seq, beg, end, meta, skip, name_size = struct.unpack(
        "<4s8i", stream.read(36)
    )
    if (
        magic != b"TBI\x01"
        or refs <= 0
        or (seq, beg, end) != (1, 2, 3)
        or flags not in (0, 0x10000)
        or skip < 0
        or not 0 <= meta <= 255
    ):
        raise ValueError("invalid BED tabix header")
    names = bytearray()
    while len(names) < name_size:
        part = stream.read(min(_BLOCK, name_size - len(names)))
        if not part:
            raise ValueError("short reference names")
        names.extend(part)
    reference_names = bytes(names).split(b"\0")
    if (
        len(reference_names) != refs + 1
        or reference_names[-1]
        or not all(reference_names[:-1])
    ):
        raise ValueError("invalid references")
    chunks = 0
    for name in reference_names[:-1]:
        (bins,) = struct.unpack("<I", stream.read(4))
        for _ in range(bins):
            bin_id, count = struct.unpack("<II", stream.read(8))
            if bin_id > 37450:
                raise ValueError("invalid bin")
            for _ in range(count):
                first, last = struct.unpack("<QQ", stream.read(16))
                if (
                    bin_id != 37450
                ):  # htslib reserves this bin for per-reference statistics.
                    check_chunk(name, bin_id, first, last, flags == 0x10000)
                    chunks += 1
        (intervals,) = struct.unpack("<I", stream.read(4))
        for _ in range(intervals):
            (offset,) = struct.unpack("<Q", stream.read(8))
            if offset >> 16 >= data_size:
                raise ValueError("linear offset outside data")
    if not chunks or len(stream.read(9)) not in (0, 8):
        raise ValueError("incomplete index")
