"""Private completion checks for pinned tracPre2 outputs, without rewriting them.

The on-disk sets below reproduce the keys in tracPre2.getUniqueBedpe and
cLoops2.ds.PET/qc.evaBedpe. They verify recorded results; they do not generate
replacement scientific output. Exact tuples avoid accepting hash collisions.
Pair provenance and child-process completion are separate required checks.
"""

from __future__ import annotations

from collections.abc import Iterator, Mapping
from dataclasses import dataclass
import csv
from decimal import Decimal, InvalidOperation
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import re
import sqlite3
import stat
import uuid


class OutputRejected(ValueError):
    """Controlled diagnostic; never includes file names or exception strings."""

    def __init__(
        self, reason_code: str, sample: str | None = None, collection: str | None = None
    ) -> None:
        super().__init__(reason_code)
        self.reason_code = reason_code
        self.sample = sample
        self.collection = collection


def summary_columns(mapq: int) -> list[str]:
    """The fixed upstream's final fifteen metrics, including its dynamic MAPQ."""
    if type(mapq) is not int or not 0 <= mapq <= 255:
        raise OutputRejected("invalid_mapq")
    return [
        "total raw sequences",
        "after linker removing sequences",
        "mapping ratio",
        f"total mapped PETs (mapq>={mapq})",
        "total mapped PETs redundancy",
        "total mapped PETs intra-chromosomal ratio",
        "total mapped PETs close ratio (distance<=1kb)",
        "total mapped PETs middle ratio (1kb<distance<=10kb)",
        "total mapped PETs distal ratio (10kb<distance)",
        "unique noBg PETs",
        "Yield",
        "unique noBg mapped PETs intra-chromosomal ratio",
        "unique noBg mapped PETs close ratio (distance<=1kb)",
        "unique noBg mapped PETs middle ratio (1kb<distance<=10kb)",
        "unique noBg mapped PETs distal ratio (10kb<distance)",
    ]


def _require(condition: bool, code: str) -> None:
    if not condition:
        raise OutputRejected(code)


def _regular(path: Path) -> None:
    _require(stat.S_ISREG(path.lstat().st_mode), "output_not_regular")


def _directory(path: Path) -> None:
    _require(stat.S_ISDIR(path.lstat().st_mode), "output_not_directory")


def _rid(value: str) -> tuple[int, int, int]:
    _require(
        re.fullmatch(r"[1-9][0-9]*_(?:-1|0|[1-9][0-9]*)_(?:-1|0|[1-9][0-9]*)", value)
        is not None,
        "invalid_internal_read_id",
    )
    a, b, c = (int(part) for part in value.split("_"))
    return a, b, c


def iter_bedpe(path: Path, contigs: Mapping[str, int]) -> Iterator[tuple[str, ...]]:
    """Read strict ten-column original BEDPE with bounded per-record memory."""
    try:
        _regular(path)
        with gzip.open(path, "rt", encoding="ascii", newline="") as stream:
            for line in stream:
                _require(line.endswith("\n"), "bedpe_truncated_line")
                row = tuple(line.removesuffix("\n").split("\t"))
                _require(len(row) == 10, "bedpe_columns")
                for chrom, start, end, strand in ((0, 1, 2, 8), (3, 4, 5, 9)):
                    _require(row[chrom] in contigs, "bedpe_unknown_contig")
                    _require(
                        re.fullmatch(r"(?:0|[1-9][0-9]*)", row[start]) is not None
                        and re.fullmatch(r"(?:0|[1-9][0-9]*)", row[end]) is not None,
                        "bedpe_coordinate",
                    )
                    _require(
                        0 <= int(row[start]) < int(row[end]) <= contigs[row[chrom]],
                        "bedpe_reference_bounds",
                    )
                    _require(row[strand] in {"+", "-"}, "bedpe_strand")
                _require(
                    re.fullmatch(r"(?:0|[1-9][0-9]*)", row[7]) is not None
                    and 0 <= int(row[7]) <= 255,
                    "bedpe_mapq",
                )
                _rid(row[6])
                yield row
    except (OSError, EOFError, UnicodeError, ValueError) as exc:
        if isinstance(exc, OutputRejected):
            raise
        raise OutputRejected("bedpe_unreadable") from None


def _trimmed_count(root: Path, sample: str, raw: int, db: sqlite3.Connection) -> int:
    """Original cutLinker writes exactly four lines and identical internal IDs."""
    paths = [root / f"{sample}_R{mate}.fastq.gz" for mate in (1, 2)]
    for path in paths:
        _regular(path)
    count = 0
    previous = 0
    with (
        gzip.open(paths[0], "rt", encoding="ascii", newline="") as r1,
        gzip.open(paths[1], "rt", encoding="ascii", newline="") as r2,
    ):
        while True:
            names = [r1.readline(), r2.readline()]
            if names == ["", ""]:
                break
            _require(all(names), "trimmed_mate_count")
            _require(names[0] == names[1], "trimmed_mate_id")
            _require(
                names[0].startswith("@") and names[0].endswith("\n"),
                "trimmed_fastq_header",
            )
            name = names[0][1:-1]
            read_number, linker1, linker2 = _rid(name)
            _require(previous < read_number <= raw, "trimmed_read_order")
            previous = read_number
            for stream, linker in ((r1, linker1), (r2, linker2)):
                seq, plus, quality = [stream.readline() for _ in range(3)]
                _require(
                    seq.endswith("\n") and quality.endswith("\n") and plus == "+\n",
                    "trimmed_fastq_format",
                )
                seq, quality = seq[:-1], quality[:-1]
                _require(
                    len(seq) >= 10
                    and len(seq) == len(quality)
                    and re.fullmatch(r"[A-Za-z]+", seq) is not None
                    and all(33 <= ord(char) <= 126 for char in quality),
                    "trimmed_fastq_format",
                )
                _require(linker == -1 or linker == len(seq), "trimmed_linker_length")
            db.execute("INSERT INTO trimmed VALUES (?)", (name,))
            count += 1
    return count


def _mapping(root: Path, samples: set[str]) -> dict[str, tuple[int, float]]:
    """Read only the fixed upstream FLAG_A tool report blocks, without guessing."""
    result: dict[str, tuple[int, float]] = {}
    logs = list(root.glob("*_tracPre2.py.log"))
    _require(len(logs) == 1, "mapping_log_count")
    _regular(logs[0])
    active: str | None = None
    count: int | None = None
    ratio: float | None = None
    with logs[0].open(encoding="utf-8") as stream:
        for line in stream:
            if "FLAG_A:" in line:
                _require(active is None, "mapping_log_nested")
                active = line.split("FLAG_A:", 1)[1].strip()
                _require(
                    active in samples and active not in result, "mapping_log_sample"
                )
                count, ratio = None, None
            elif line.strip() == "FLAG_A":
                _require(
                    active is not None and count is not None and ratio is not None,
                    "mapping_log_incomplete",
                )
                result[active] = (count, ratio)
                active = None
            elif active is not None:
                header = re.fullmatch(r"([0-9]+) reads; of these:\n", line)
                rate = re.fullmatch(
                    r"([0-9]+(?:\.[0-9]+)?)% overall alignment rate\n", line
                )
                if header:
                    _require(count is None, "mapping_log_duplicate")
                    count = int(header.group(1))
                if rate:
                    _require(ratio is None, "mapping_log_duplicate")
                    ratio = float(rate.group(1)) / 100
                    _require(0 <= ratio <= 1, "mapping_ratio_range")
    _require(active is None and set(result) == samples, "mapping_log_incomplete")
    return result


@dataclass
class _Stats:
    total: int = 0
    unique: int = 0
    cis: int = 0
    close: int = 0
    middle: int = 0
    far: int = 0

    def add(
        self, row: tuple[str, ...], db: sqlite3.Connection, collection: str
    ) -> None:
        self.total += 1
        a = (row[0], int(row[1]), int(row[2]), row[8])
        b = (row[3], int(row[4]), int(row[5]), row[9])
        cis = a[0] == b[0]
        if (cis and a[1] + a[2] > b[1] + b[2]) or (not cis and a[0] > b[0]):
            a, b = b, a
        key = json.dumps((a, b), separators=(",", ":"))
        inserted = db.execute(
            "INSERT OR IGNORE INTO qc_keys VALUES (?, ?)", (collection, key)
        ).rowcount
        if not inserted:
            return
        self.unique += 1
        if cis:
            self.cis += 1
            # PET truncates each nonnegative endpoint midpoint before subtraction.
            distance = abs(int((a[1] + a[2]) / 2) - int((b[1] + b[2]) / 2))
            self.close += distance <= 1000
            self.middle += 1000 < distance <= 10000
            self.far += 10000 < distance

    def ratios(self) -> list[float]:
        return [
            self.cis / self.unique,
            self.close / self.cis,
            self.middle / self.cis,
            self.far / self.cis,
        ]


def _no_bg_eligible(row: tuple[str, ...]) -> bool:
    _, first, second = _rid(row[6])
    if first + second != -2:
        return True
    if row[0] != row[3]:
        return False
    # getUniqueBedpe uses untruncated midpoints (unlike the QC distance above).
    distance = abs((int(row[1]) + int(row[2])) / 2 - (int(row[4]) + int(row[5])) / 2)
    return distance >= 1000


def _policy(stats: _Stats, sample: str, collection: str) -> None:
    if stats.total == 0:
        raise OutputRejected("empty_pet_set", sample, collection)
    if stats.cis == 0:
        raise OutputRejected("no_cis_denominator", sample, collection)


def _artifact(path: Path, root: Path) -> dict[str, str | int]:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return {
        "path": path.relative_to(root).as_posix(),
        "size_bytes": path.stat().st_size,
        "sha256": digest.hexdigest(),
    }


def _sample(
    root: Path,
    sample: str,
    raw: int,
    mapq: int,
    contigs: Mapping[str, int],
    mapping: tuple[int, float],
    metrics: list[float],
    db: sqlite3.Connection,
) -> dict[str, object]:
    directory = root / sample
    _directory(directory)
    _regular(directory / f"{sample}.bam")
    trim = _trimmed_count(directory, sample, raw, db)
    _require(mapping[0] == trim, "mapping_trim_count")
    paths = [directory / f"{sample}_{suffix}.bedpe.gz" for suffix in ("all", "unique")]
    all_stats, no_bg_stats = _Stats(), _Stats()
    for row in iter_bedpe(paths[0], contigs):
        _require(int(row[7]) >= mapq, "bedpe_below_mapq")
        _require(
            db.execute("SELECT 1 FROM trimmed WHERE name=?", (row[6],)).fetchone()
            is not None,
            "bedpe_unknown_trimmed_read",
        )
        all_stats.add(row, db, "all")
        if _no_bg_eligible(row):
            # Background removal keeps the first row for the raw first six fields;
            # it ignores strand, unlike the separately normalized QC identity.
            coordinates = "\t".join(row[:6])
            db.execute(
                "INSERT OR IGNORE INTO expected_no_bg VALUES (?, ?, 0)",
                (coordinates, "\t".join(row)),
            )
    _policy(all_stats, sample, "all")
    for row in iter_bedpe(paths[1], contigs):
        changed = db.execute(
            "UPDATE expected_no_bg SET seen=1 WHERE row=? AND seen=0", ("\t".join(row),)
        ).rowcount
        _require(changed == 1, "noBg_content_mismatch")
        no_bg_stats.add(row, db, "noBg")
    _policy(no_bg_stats, sample, "noBg")
    _require(
        db.execute("SELECT 1 FROM expected_no_bg WHERE seen=0 LIMIT 1").fetchone()
        is None,
        "noBg_missing_pet",
    )
    # Number of BAM alignments/PETs need not be bounded by input pairs when a
    # future tool emits multiple alignments. Pair provenance is checked separately.
    expected = [
        raw,
        trim,
        mapping[1],
        all_stats.total,
        1 - all_stats.unique / all_stats.total,
        *all_stats.ratios(),
        no_bg_stats.total,
        no_bg_stats.total / raw,
        *no_bg_stats.ratios(),
    ]
    _require(
        all(metrics[index] == expected[index] for index in (0, 1, 3, 9)),
        "summary_count_ratio_mismatch",
    )
    _require(
        all(
            math.isclose(a, b, rel_tol=0, abs_tol=1e-12)
            for a, b in zip(metrics, expected)
        ),
        "summary_count_ratio_mismatch",
    )
    return {
        "raw": raw,
        "trimmed": trim,
        "all": all_stats.total,
        "noBg": no_bg_stats.total,
        "qc_all_unique": all_stats.unique,
        "qc_all_cis": all_stats.cis,
        "metrics": metrics,
        "artifacts": [_artifact(path, root) for path in paths],
    }


def verify_outputs(
    output_dir: Path,
    raw_counts: Mapping[str, int],
    contigs: Mapping[str, int],
    mapq: int,
    work_dir: Path,
) -> dict[str, object]:
    """Verify all samples or reject the attempt; preserve every input/output byte.

    SQLite scratch stays in the caller's private work_dir, including on failure.
    This check alone never writes a completion marker or authorizes publication.
    """
    columns = summary_columns(mapq)
    _require(bool(raw_counts), "missing_samples")
    _require(
        all(re.fullmatch(r"s[0-9]{6}", sample) for sample in raw_counts), "sample_token"
    )
    _require(
        all(type(count) is int and count >= 0 for count in raw_counts.values()),
        "raw_count_invalid",
    )
    sample_context: str | None = None
    try:
        _directory(output_dir)
        _directory(work_dir)
        directories = {path.name for path in output_dir.iterdir() if path.is_dir()}
        _require(directories == set(raw_counts), "output_samples")
        summary = output_dir / "tracPre_summary.txt"
        _regular(summary)
        table: dict[str, list[float]] = {}
        with summary.open(encoding="ascii", newline="") as stream:
            rows = csv.reader(stream, delimiter="\t")
            _require(next(rows, None) == [""] + columns, "final_summary_header")
            for row in rows:
                _require(len(row) == 16, "summary_width")
                _require(
                    row[0] in raw_counts and row[0] not in table, "summary_samples"
                )
                values = [float(value) for value in row[1:]]
                _require(
                    all(math.isfinite(value) for value in values), "summary_nonfinite"
                )
                for index in (0, 1, 3, 9):
                    count = Decimal(row[index + 1])
                    _require(
                        count == count.to_integral_value() and count >= 0,
                        "summary_count",
                    )
                    values[index] = int(count)
                table[row[0]] = values
        _require(set(table) == set(raw_counts), "summary_samples")
        mapping = _mapping(output_dir, set(raw_counts))
        results = {}
        for sample, raw in raw_counts.items():
            sample_context = sample
            path = work_dir / f"output-{sample}-{uuid.uuid4().hex}.sqlite"
            descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
            os.close(descriptor)
            db = sqlite3.connect(path)
            try:
                db.execute("PRAGMA cache_size=-4096")
                db.execute("PRAGMA temp_store=FILE")
                db.executescript(
                    "CREATE TABLE trimmed (name TEXT PRIMARY KEY);"
                    "CREATE TABLE qc_keys (collection TEXT, key TEXT, "
                    "PRIMARY KEY(collection,key));"
                    "CREATE TABLE expected_no_bg (coordinates TEXT PRIMARY KEY, "
                    "row TEXT UNIQUE, seen INTEGER);"
                )
                results[sample] = _sample(
                    output_dir,
                    sample,
                    raw,
                    mapq,
                    contigs,
                    mapping[sample],
                    table[sample],
                    db,
                )
                db.commit()
            finally:
                db.close()
        return {"samples": results, "summary": _artifact(summary, output_dir)}
    except OutputRejected as exc:
        if exc.sample is None:
            exc.sample = sample_context
        raise
    except (
        OSError,
        EOFError,
        UnicodeError,
        ValueError,
        InvalidOperation,
        sqlite3.Error,
        csv.Error,
    ):
        raise OutputRejected("output_unreadable", sample_context) from None
