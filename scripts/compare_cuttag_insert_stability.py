#!/usr/bin/env python3
"""Compare CUT&Tag paired-end templates at max_insert 500 and 700.

Template identity is the read query name. A record has ``qname``,
``r1_chrom``, ``r1_pos``, ``r2_chrom``, ``r2_pos``, ``length``, and
``mapq``. ``length`` is the absolute template length. ``mapq`` is already
the minimum of the two mates. This script does not read BAM files, does
not shell out to samtools, and does not compute peak Jaccard or FRiP.

Newly admitted templates are query names present in the 700 arm and
absent from the 500 arm, counted by insert-length bin and MAPQ bin.
Short-template stability is the fraction of 500-arm templates with
MAPQ >= 30 and length <= 500 whose R1 position, R2 position, or ordered
mate contigs differ in the 700 arm. Counts are reported only. This script
does not decide a default.
"""

from __future__ import annotations

import argparse
import json
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import TypedDict


TEMPLATE_COLUMNS = (
    "qname",
    "r1_chrom",
    "r1_pos",
    "r2_chrom",
    "r2_pos",
    "length",
    "mapq",
)
_INT_FIELDS = ("r1_pos", "r2_pos", "length", "mapq")
_LENGTH_BIN_UPPER = (
    (150, "<=150"),
    (300, "151-300"),
    (500, "301-500"),
    (700, "501-700"),
)
LENGTH_BINS = tuple(name for _upper, name in _LENGTH_BIN_UPPER) + (">700",)
MAPQ_BINS = ("<30", ">=30")
HIGH_CONFIDENCE_MAPQ = 30
SHORT_TEMPLATE_MAX_LENGTH = 500


class Template(TypedDict):
    qname: str
    r1_chrom: str
    r1_pos: int
    r2_chrom: str
    r2_pos: int
    length: int
    mapq: int


class NewTemplateCounts(TypedDict):
    total: int
    by_length_bin: dict[str, int]
    by_mapq_bin: dict[str, int]
    by_length_and_mapq: dict[str, dict[str, int]]


class ShortTemplateStability(TypedDict):
    denominator: int
    changed: int
    fraction: float | None


def _require_text(value: object, field: str) -> str:
    if not isinstance(value, str) or value == "":
        raise ValueError(f"{field} must be a non-empty string")
    return value


def _require_int(value: object, qname: str, field: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"{qname} {field} must be an int")
    return value


def _validate_template(record: Mapping[str, object]) -> Template:
    missing = [key for key in TEMPLATE_COLUMNS if key not in record]
    if missing:
        raise ValueError("template record missing " + ", ".join(missing))
    qname = _require_text(record["qname"], "qname")
    r1_chrom = _require_text(record["r1_chrom"], "r1_chrom")
    r2_chrom = _require_text(record["r2_chrom"], "r2_chrom")
    values = {field: _require_int(record[field], qname, field) for field in _INT_FIELDS}
    if values["length"] < 0:
        raise ValueError(f"{qname} length must be non-negative")
    if values["mapq"] < 0:
        raise ValueError(f"{qname} mapq must be non-negative")
    return {
        "qname": qname,
        "r1_chrom": r1_chrom,
        "r1_pos": values["r1_pos"],
        "r2_chrom": r2_chrom,
        "r2_pos": values["r2_pos"],
        "length": values["length"],
        "mapq": values["mapq"],
    }


def _index_templates(
    arm: Sequence[Mapping[str, object]],
) -> dict[str, Template]:
    indexed: dict[str, Template] = {}
    for record in arm:
        template = _validate_template(record)
        qname = template["qname"]
        if qname in indexed:
            raise ValueError(f"duplicate template qname: {qname}")
        indexed[qname] = template
    return indexed


def length_bin(length: int) -> str:
    """Return the insert-length bin for an absolute template length."""
    if length < 0:
        raise ValueError("template length must be non-negative")
    for upper, name in _LENGTH_BIN_UPPER:
        if length <= upper:
            return name
    return ">700"


def mapq_bin(mapq: int) -> str:
    """Return ``<30`` or ``>=30`` for a template MAPQ."""
    if mapq < 0:
        raise ValueError("template MAPQ must be non-negative")
    if mapq < HIGH_CONFIDENCE_MAPQ:
        return "<30"
    return ">=30"


def _empty_counts() -> NewTemplateCounts:
    by_length = {name: 0 for name in LENGTH_BINS}
    by_mapq = {name: 0 for name in MAPQ_BINS}
    joint = {name: {mapq_name: 0 for mapq_name in MAPQ_BINS} for name in LENGTH_BINS}
    return {
        "total": 0,
        "by_length_bin": by_length,
        "by_mapq_bin": by_mapq,
        "by_length_and_mapq": joint,
    }


def classify_new_templates(
    arm500: Sequence[Mapping[str, object]],
    arm700: Sequence[Mapping[str, object]],
) -> NewTemplateCounts:
    """Count 700-arm templates whose query name is absent from the 500 arm.

    Length bins are ``<=150``, ``151-300``, ``301-500``, ``501-700``, and
    ``>700``. MAPQ bins are ``<30`` and ``>=30``. Empty bins are returned
    as zeros. ``by_length_and_mapq`` is the joint count. A query name
    present in both arms is not newly admitted, even if its coordinates
    differ.
    """
    present_in_500 = set(_index_templates(arm500))
    incoming = _index_templates(arm700)
    counts = _empty_counts()
    for qname, template in incoming.items():
        if qname in present_in_500:
            continue
        length_name = length_bin(template["length"])
        mapq_name = mapq_bin(template["mapq"])
        counts["by_length_bin"][length_name] += 1
        counts["by_mapq_bin"][mapq_name] += 1
        counts["by_length_and_mapq"][length_name][mapq_name] += 1
        counts["total"] += 1
    return counts


def _is_high_confidence_short(template: Template) -> bool:
    return (
        template["mapq"] >= HIGH_CONFIDENCE_MAPQ
        and template["length"] <= SHORT_TEMPLATE_MAX_LENGTH
    )


def _assignment_changed(reference: Template, observed: Template) -> bool:
    """Return True when mate contigs or coordinates differ.

    Mate assignment is the ordered ``(r1_chrom, r2_chrom)`` pair. A change
    of either contig counts even when both positions are numerically equal.
    Coordinate changes are position changes while those contigs are kept.
    """
    same_mates = (
        reference["r1_chrom"] == observed["r1_chrom"]
        and reference["r2_chrom"] == observed["r2_chrom"]
    )
    if not same_mates:
        return True
    return (
        reference["r1_pos"] != observed["r1_pos"]
        or reference["r2_pos"] != observed["r2_pos"]
    )


def short_template_stability(
    arm500: Sequence[Mapping[str, object]],
    arm700: Sequence[Mapping[str, object]],
) -> ShortTemplateStability:
    """Return how often high-confidence short templates move.

    The denominator is every 500-arm template with MAPQ >= 30 and length
    <= 500. MAPQ < 30 or length > 500 is excluded. A denominator template
    counts as changed when the 700-arm record has a different R1 position,
    R2 position, or ordered mate-contig pair. Identical placement
    contributes no change. A template present only in the 700 arm is not
    in this denominator and is not a stability failure. A qualifying
    template with no 700-arm record stays in the denominator and does not
    increment ``changed``; loss is a separate count-conservation check.
    ``fraction`` is ``None`` when the denominator is 0.
    """
    reference = _index_templates(arm500)
    observed = _index_templates(arm700)
    denominator = 0
    changed = 0
    for qname, template in reference.items():
        if not _is_high_confidence_short(template):
            continue
        denominator += 1
        match = observed.get(qname)
        if match is not None and _assignment_changed(template, match):
            changed += 1
    fraction = None if denominator == 0 else changed / denominator
    return {"denominator": denominator, "changed": changed, "fraction": fraction}


def _parse_int(text: str, qname: str, field: str) -> int:
    try:
        return int(text, 10)
    except ValueError as error:
        raise ValueError(f"{qname} {field} is not an integer") from error


def load_templates(path: Path) -> list[Template]:
    """Load template rows from a TSV with ``TEMPLATE_COLUMNS``."""
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines:
        raise ValueError(f"empty TSV: {path}")
    header = tuple(lines[0].split("\t"))
    if header != TEMPLATE_COLUMNS:
        expected = "\t".join(TEMPLATE_COLUMNS)
        raise ValueError(f"unexpected TSV header in {path}; expected {expected}")
    records: list[Template] = []
    for line_number, line in enumerate(lines[1:], start=2):
        if not line.strip():
            continue
        cells = line.split("\t")
        if len(cells) != len(TEMPLATE_COLUMNS):
            raise ValueError(f"column count mismatch in {path}:{line_number}")
        raw = dict(zip(TEMPLATE_COLUMNS, cells, strict=True))
        qname = raw["qname"]
        label = qname or f"line {line_number}"
        parsed: dict[str, object] = {
            "qname": qname,
            "r1_chrom": raw["r1_chrom"],
            "r2_chrom": raw["r2_chrom"],
        }
        for field in _INT_FIELDS:
            parsed[field] = _parse_int(raw[field], label, field)
        records.append(_validate_template(parsed))
    return records


def comparison_record(
    arm500: Sequence[Mapping[str, object]],
    arm700: Sequence[Mapping[str, object]],
) -> dict[str, object]:
    """Bundle the two max_insert contrasts without choosing a default."""
    return {
        "scope": "single max_insert contrast; not a default decision",
        "reference_arm": "max_insert_500",
        "candidate_arm": "max_insert_700",
        "template_identity": "qname",
        "length_definition": "absolute template length",
        "mapq_definition": "minimum of the two mates",
        "length_bins": list(LENGTH_BINS),
        "mapq_bins": list(MAPQ_BINS),
        "newly_admitted_templates": classify_new_templates(arm500, arm700),
        "short_template_stability": short_template_stability(arm500, arm700),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--arm-500", type=Path, required=True)
    parser.add_argument("--arm-700", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    record = comparison_record(
        load_templates(args.arm_500),
        load_templates(args.arm_700),
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(record, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
