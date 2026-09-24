"""Check the BAM origins of emitted PETs without changing scientific outputs.

One invocation covers one sample. Each emitted BEDPE row needs two distinct BAM
records with its exact qname, one unambiguous read1 and one unambiguous read2.
All matching alignments are eligible, including secondary/supplementary records;
proper-pair, mapping quality and duplicate flags are not additional filters.

This is an existence check, not reconstruction of the converter's history. A
matching record may support multiple emitted rows, which are each checked and
never deduplicated. It does not establish unique multiplicity or that every BAM
record was emitted. Ambiguous-role records cannot supply support, but unrelated
orphans/ambiguous records alone do not reject an otherwise supported result.
"""

from collections.abc import Iterable
import os
from pathlib import Path
import re
import sqlite3

_CIGAR_TOKEN = re.compile(r"([1-9][0-9]*)([MIDNSHP=X])")
_REFERENCE_OPS = frozenset("MDN=X")
_BATCH_SIZE = 256


class PairError(ValueError):
    """A controlled reason code; raw data and paths are never in its message."""

    def __init__(self, reason_code: str) -> None:
        self.reason_code = reason_code
        super().__init__(reason_code)


def _token(value: str) -> bool:
    return bool(value) and not any(character.isspace() for character in value)


def _reference_span(cigar: str) -> int:
    end = 0
    span = 0
    for match in _CIGAR_TOKEN.finditer(cigar):
        if match.start() != end:
            raise PairError("pair_sam_invalid")
        end = match.end()
        if match.group(2) in _REFERENCE_OPS:
            span += int(match.group(1))
    if end != len(cigar) or span == 0:
        raise PairError("pair_sam_invalid")
    return span


def _alignment(line: str) -> tuple[str, int, str, int, int, str] | None:
    if line.startswith("@"):
        return None
    fields = line.rstrip("\r\n").split("\t")
    if len(fields) < 11 or not _token(fields[0]):
        raise PairError("pair_sam_invalid")
    try:
        flag = int(fields[1])
        position = int(fields[3])
    except ValueError:
        raise PairError("pair_sam_invalid") from None
    if not 0 <= flag <= 65535 or position < 0:
        raise PairError("pair_sam_invalid")
    if flag & 0x4:
        return None
    if position == 0 or fields[2] == "*" or not _token(fields[2]):
        raise PairError("pair_sam_invalid")
    span = _reference_span(fields[5])
    role = flag & (0x40 | 0x80)
    if role not in (0x40, 0x80):
        return None
    start = position - 1
    return fields[0], role, fields[2], start, start + span, "-" if flag & 0x10 else "+"


def _pet(
    fields: tuple[str, ...],
) -> tuple[str, tuple[str, int, int, str], tuple[str, int, int, str]]:
    if len(fields) != 10 or not _token(fields[6]):
        raise PairError("pair_bedpe_invalid")
    ends = []
    for offset, strand in ((0, fields[8]), (3, fields[9])):
        try:
            start, end = int(fields[offset + 1]), int(fields[offset + 2])
        except ValueError:
            raise PairError("pair_bedpe_invalid") from None
        contig = fields[offset]
        if (
            not _token(contig)
            or contig == "*"
            or start < 0
            or end <= start
            or strand not in ("+", "-")
        ):
            raise PairError("pair_bedpe_invalid")
        ends.append((contig, start, end, strand))
    return fields[6], ends[0], ends[1]


def verify_pairs(
    sam_lines: Iterable[str],
    bedpe_rows: Iterable[tuple[str, ...]],
    database: Path,
) -> None:
    """Verify one sample using a new, private on-disk association database.

    SAM ordering is irrelevant: in particular natural qname sort is never
    treated as lexicographic order. SQLite's bounded page cache plus small insert
    batches avoid keeping a sample, or a high-multiplicity qname, in Python RAM.
    The database remains as private evidence on both success and failure; its
    connection and journal are closed before returning or raising.
    """
    connection = None
    try:
        descriptor = os.open(database, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
        os.close(descriptor)
        connection = sqlite3.connect(database)
        connection.execute("PRAGMA cache_size = -2048")
        connection.execute("PRAGMA temp_store = FILE")
        connection.execute(
            "CREATE TABLE alignments (qname TEXT NOT NULL, role INTEGER NOT NULL, "
            "contig TEXT NOT NULL, start INTEGER NOT NULL, end INTEGER NOT NULL, "
            "strand TEXT NOT NULL)"
        )
        connection.execute(
            "CREATE INDEX support ON alignments(qname, role, contig, start, end, strand)"
        )
        pending = []
        for line in sam_lines:
            alignment = _alignment(line)
            if alignment is not None:
                pending.append(alignment)
            if len(pending) == _BATCH_SIZE:
                connection.executemany(
                    "INSERT INTO alignments VALUES (?, ?, ?, ?, ?, ?)", pending
                )
                connection.commit()
                pending.clear()
        if pending:
            connection.executemany(
                "INSERT INTO alignments VALUES (?, ?, ?, ?, ?, ?)", pending
            )
        connection.commit()

        def supports(
            qname: str, role: int, endpoint: tuple[str, int, int, str]
        ) -> bool:
            # The lookup considers every matching record, not a chosen first
            # alignment. LIMIT 1 only stops once existence has been established.
            return (
                connection.execute(
                    "SELECT 1 FROM alignments WHERE "
                    "qname = ? AND role = ? AND contig = ? AND start = ? "
                    "AND end = ? AND strand = ? LIMIT 1",
                    (qname, role, *endpoint),
                ).fetchone()
                is not None
            )

        for row in bedpe_rows:
            qname, first, second = _pet(row)
            if not (
                supports(qname, 0x40, first) and supports(qname, 0x80, second)
            ) and not (supports(qname, 0x40, second) and supports(qname, 0x80, first)):
                raise PairError("pair_support_missing")
    except (OSError, sqlite3.Error, OverflowError):
        raise PairError("pair_index_unavailable") from None
    finally:
        if connection is not None:
            connection.close()
