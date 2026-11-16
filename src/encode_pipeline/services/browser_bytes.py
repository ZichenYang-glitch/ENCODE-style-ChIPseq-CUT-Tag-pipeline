"""Authorized byte ranges over an already verified artifact descriptor."""

from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass
import os
import re


_RANGE = re.compile(r"^bytes=(\d*)-(\d*)$")
_CHUNK = 64 * 1024


class BrowserByteError(Exception):
    """A range could not be satisfied. The caller still closes the plan."""

    def __init__(self, reason: str) -> None:
        super().__init__(reason)
        self.reason = reason


@dataclass(frozen=True)
class ByteRange:
    """Inclusive HTTP byte interval."""

    start: int
    end: int

    @property
    def length(self) -> int:
        return self.end - self.start + 1


def parse_byte_range(header: str | None, size: int) -> ByteRange | None:
    """Return one range, None for a full read, or raise for 416 cases."""

    if size < 0:
        raise BrowserByteError("invalid_range")
    if header is None:
        return None
    if "\n" in header or "\r" in header or "," in header:
        raise BrowserByteError("multiple_ranges" if "," in header else "invalid_range")
    match = _RANGE.fullmatch(header.strip())
    if match is None or size == 0:
        raise BrowserByteError("invalid_range")
    start_text, end_text = match.group(1), match.group(2)
    if start_text == "" and end_text == "":
        raise BrowserByteError("invalid_range")
    if start_text == "":
        suffix = int(end_text)
        if suffix <= 0:
            raise BrowserByteError("unsatisfiable_range")
        if suffix >= size:
            return ByteRange(0, size - 1)
        return ByteRange(size - suffix, size - 1)
    start = int(start_text)
    end = size - 1 if end_text == "" else int(end_text)
    if start >= size or start > end:
        raise BrowserByteError("unsatisfiable_range")
    return ByteRange(start, min(end, size - 1))


def iter_plan_bytes(plan: object, byte_range: ByteRange | None) -> Iterator[bytes]:
    """Yield a verified plan's full body or one range, then close it."""

    try:
        content = getattr(plan, "_content", None)
        size = int(getattr(plan, "size_bytes"))
        if byte_range is None:
            start, length = 0, size
        else:
            start, length = byte_range.start, byte_range.length
        if content is not None:
            view = content[start : start + length]
            if len(view) != length:
                raise BrowserByteError("short_read")
            for offset in range(0, len(view), _CHUNK):
                yield view[offset : offset + _CHUNK]
            return
        descriptor = getattr(plan, "_final_descriptor", None)
        nodes = getattr(plan, "_nodes", ())
        if not isinstance(descriptor, int):
            raise BrowserByteError("short_read")
        _verify(nodes)
        os.lseek(descriptor, start, os.SEEK_SET)
        remaining = length
        while remaining:
            _verify(nodes)
            chunk = os.read(descriptor, min(_CHUNK, remaining))
            if not chunk:
                raise BrowserByteError("short_read")
            _verify(nodes)
            remaining -= len(chunk)
            yield chunk
        _verify(nodes)
    finally:
        close = getattr(plan, "close", None)
        if callable(close):
            close()


def _verify(nodes: object) -> None:
    if not nodes:
        return
    from encode_pipeline.services.artifact_downloads import (
        _verify_descriptor_chain,
    )

    _verify_descriptor_chain(nodes)
