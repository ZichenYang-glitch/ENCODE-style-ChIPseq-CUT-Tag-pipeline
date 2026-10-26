"""Complete-template parsing shared by fragment-size and MNase QC writers."""

from __future__ import annotations

from dataclasses import dataclass
import os
import subprocess


SAMTOOLS = os.environ.get("SAMTOOLS", "samtools")

_FPAIRED = 0x1
_FPROPER_PAIR = 0x2
_FUNMAP = 0x4
_FMUNMAP = 0x8
_FFIRST = 0x40
_FSECOND = 0x80
_FSECONDARY = 0x100
_FSUPPLEMENTARY = 0x800


@dataclass(frozen=True)
class SamRecord:
    qname: str
    flag: int
    rname: str
    pos: int
    cigar: str
    rnext: str
    pnext: int
    tlen: int


@dataclass(frozen=True)
class Template:
    qname: str
    records: tuple[SamRecord, ...]


@dataclass(frozen=True)
class TemplateClassification:
    valid: bool
    length: int | None
    reason: str | None


def cigar_ref_length(cigar: str) -> int:
    """Return the reference-consuming length of a CIGAR string."""
    length = 0
    digits = ""
    for char in cigar:
        if char.isdigit():
            digits += char
            continue
        if not digits:
            continue
        count = int(digits)
        digits = ""
        if char in ("M", "D", "N", "=", "X"):
            length += count
    return length


def _parse_record(line: str) -> SamRecord | None:
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 9:
        return None
    try:
        return SamRecord(
            qname=fields[0],
            flag=int(fields[1]),
            rname=fields[2],
            pos=int(fields[3]),
            cigar=fields[5],
            rnext=fields[6],
            pnext=int(fields[7]),
            tlen=int(fields[8]),
        )
    except (ValueError, IndexError):
        return None


class TemplateReader:
    """Stream templates from one BAM and retain the alignment count."""

    def __init__(self, bam: str, *, samtools: str = SAMTOOLS):
        self.bam = bam
        self.samtools = samtools
        self.alignments_total = 0
        self._process: subprocess.Popen | None = None

    def __enter__(self) -> "TemplateReader":
        self._process = subprocess.Popen(
            [self.samtools, "view", self.bam],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        return self

    def __iter__(self):
        process = self._process
        if process is None or process.stdout is None:
            raise RuntimeError("template reader is not open")
        pending: dict[str, tuple[SamRecord, ...]] = {}
        last_completed = ""
        for line in process.stdout:
            if not line or line.startswith("@"):
                continue
            self.alignments_total += 1
            record = _parse_record(line)
            if record is None:
                raise ValueError("BAM alignment record is malformed")
            if record.flag & (_FSECONDARY | _FSUPPLEMENTARY):
                continue
            if record.qname == last_completed:
                raise ValueError("template has more than two primary records")
            existing = pending.get(record.qname)
            if existing is None:
                pending[record.qname] = (record,)
            elif len(existing) == 1:
                del pending[record.qname]
                last_completed = record.qname
                yield Template(record.qname, (existing[0], record))
            else:  # pragma: no cover - defensive invariant
                raise ValueError("template has more than two primary records")
        for qname, records in pending.items():
            yield Template(qname, records)

    def __exit__(self, exc_type, exc, traceback) -> bool:
        process = self._process
        self._process = None
        if process is None:
            return False
        if process.stdout is not None:
            process.stdout.close()
        stderr = process.stderr.read() if process.stderr is not None else ""
        returncode = process.wait()
        if exc_type is not None:
            return False
        if returncode != 0:
            raise RuntimeError(f"samtools view failed on {self.bam}: {stderr.strip()}")
        return False


def classify_pe_template(template: Template) -> TemplateClassification:
    """Classify one complete PE template using the analysis-BAM contract.

    Valid templates keep legal negative TLEN and use ``abs(TLEN)`` exactly
    once.  Zero TLEN, orphan/missing mates, cross-contig pairs and mutually
    inconsistent mates are invalid and never produce a biological length.
    """
    records = template.records
    if len(records) != 2:
        return TemplateClassification(False, None, "orphan")
    first, second = records
    if (first.flag | second.flag) & (_FSECONDARY | _FSUPPLEMENTARY):
        return TemplateClassification(False, None, "inconsistent")
    if not (first.flag & _FPAIRED) or not (second.flag & _FPAIRED):
        return TemplateClassification(False, None, "inconsistent")
    if (first.flag | second.flag) & (_FUNMAP | _FMUNMAP):
        return TemplateClassification(False, None, "inconsistent")
    if not (first.flag & _FPROPER_PAIR) or not (second.flag & _FPROPER_PAIR):
        return TemplateClassification(False, None, "inconsistent")
    if not (
        (first.flag & _FFIRST and second.flag & _FSECOND)
        or (first.flag & _FSECOND and second.flag & _FFIRST)
    ):
        return TemplateClassification(False, None, "inconsistent")
    if first.rname != second.rname:
        return TemplateClassification(False, None, "cross_contig")
    if first.rnext not in ("=", first.rname) or second.rnext not in (
        "=",
        second.rname,
    ):
        return TemplateClassification(False, None, "cross_contig")
    if first.pos <= 0 or second.pos <= 0:
        return TemplateClassification(False, None, "inconsistent")
    if cigar_ref_length(first.cigar) <= 0 or cigar_ref_length(second.cigar) <= 0:
        return TemplateClassification(False, None, "inconsistent")
    if first.pos != second.pnext or second.pos != first.pnext:
        return TemplateClassification(False, None, "inconsistent")
    if first.tlen == 0 or second.tlen == 0:
        return TemplateClassification(False, None, "tlen_zero")
    if first.tlen != -second.tlen:
        return TemplateClassification(False, None, "inconsistent")
    length = abs(first.tlen)
    if length <= 0:  # pragma: no cover - guarded by nonzero reciprocal TLEN
        return TemplateClassification(False, None, "tlen_zero")
    return TemplateClassification(True, length, None)
