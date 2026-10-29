#!/usr/bin/env python3
"""Compare one CUT&Tag retain arm with one exclude arm.

The formulas are the SCI-7 comparison design applied to a single
treatment/control pair. Peak count is descriptive. This script does not
decide a default.

FRiP uses analysis-BAM read records, the same denominator as
``scripts/calc_frip.py`` (``samtools view -c``). The numerator is reads
overlapping the common peak universe. Off-peak fraction is ``1 - FRiP``
on that universe. Per-base densities divide those read counts by the
common-universe base pairs and by the remaining reference bases in
``chrom.sizes``. Enrichment is signal density / background density.
"""

from __future__ import annotations

import argparse
import json
import statistics
import subprocess
import sys
from pathlib import Path


FRIP_EPSILON = 1e-6
RETAIN = "retain"
EXCLUDE = "exclude"
FIXED_QC_FIELDS = (
    "resolved_max_insert",
    "blacklist_mode",
    "organelle_mode",
    "layout",
    "assay",
    "peak_mode",
    "genome",
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines:
        raise ValueError(f"empty TSV: {path}")
    header = lines[0].split("\t")
    rows = []
    for line in lines[1:]:
        if not line.strip():
            continue
        cells = line.split("\t")
        if len(cells) != len(header):
            raise ValueError(f"column count mismatch in {path}")
        rows.append(dict(zip(header, cells, strict=True)))
    return rows


def load_qc_summary(path: Path, sample: str) -> dict[str, str]:
    rows = [row for row in read_tsv(path) if row.get("sample") == sample]
    if len(rows) != 1:
        raise ValueError(f"{path} must contain exactly one row for {sample}")
    return rows[0]


def parse_narrow_peaks(path: Path) -> list[tuple[str, int, int]]:
    peaks = []
    for line_number, line in enumerate(
        path.read_text(encoding="utf-8").splitlines(), 1
    ):
        if not line.strip() or line.startswith("#") or line.startswith("track"):
            continue
        cells = line.split("\t")
        if len(cells) < 3:
            raise ValueError(f"{path}:{line_number} is not a narrowPeak interval")
        chrom, start, end = cells[0], int(cells[1]), int(cells[2])
        if end < start:
            raise ValueError(f"{path}:{line_number} has end < start")
        peaks.append((chrom, start, end))
    return peaks


def merge_intervals(
    peaks: list[tuple[str, int, int]],
) -> list[tuple[str, int, int]]:
    """Merge overlaps and book-ended intervals, matching bedtools merge."""
    merged: list[tuple[str, int, int]] = []
    for chrom, start, end in sorted(peaks):
        if merged and merged[-1][0] == chrom and start <= merged[-1][2]:
            prev_chrom, prev_start, prev_end = merged[-1]
            merged[-1] = (prev_chrom, prev_start, max(prev_end, end))
        else:
            merged.append((chrom, start, end))
    return merged


def interval_bp(peaks: list[tuple[str, int, int]]) -> int:
    return sum(end - start for _chrom, start, end in peaks)


def intersection_bp(
    left: list[tuple[str, int, int]],
    right: list[tuple[str, int, int]],
) -> int:
    by_chrom: dict[str, list[tuple[int, int]]] = {}
    for chrom, start, end in right:
        by_chrom.setdefault(chrom, []).append((start, end))
    total = 0
    for chrom, start, end in left:
        for other_start, other_end in by_chrom.get(chrom, []):
            overlap = min(end, other_end) - max(start, other_start)
            if overlap > 0:
                total += overlap
    return total


def overlaps_any(
    peak: tuple[str, int, int],
    others: list[tuple[str, int, int]],
) -> bool:
    chrom, start, end = peak
    for other_chrom, other_start, other_end in others:
        if other_chrom != chrom:
            continue
        if min(end, other_end) - max(start, other_start) > 0:
            return True
    return False


def overlap_fraction(
    query: list[tuple[str, int, int]],
    target: list[tuple[str, int, int]],
) -> float | None:
    if not query:
        return None
    hits = sum(1 for peak in query if overlaps_any(peak, target))
    return hits / len(query)


def peak_widths(peaks: list[tuple[str, int, int]]) -> dict[str, float | None]:
    widths = [end - start for _chrom, start, end in peaks]
    if not widths:
        return {"n": 0, "mean_bp": None, "median_bp": None}
    return {
        "n": len(widths),
        "mean_bp": statistics.fmean(widths),
        "median_bp": statistics.median(widths),
    }


def bp_jaccard(
    left: list[tuple[str, int, int]],
    right: list[tuple[str, int, int]],
) -> dict[str, float | int | None]:
    left_merged = merge_intervals(left)
    right_merged = merge_intervals(right)
    left_bp = interval_bp(left_merged)
    right_bp = interval_bp(right_merged)
    shared = intersection_bp(left_merged, right_merged)
    union = left_bp + right_bp - shared
    return {
        "retain_merged_bp": left_bp,
        "exclude_merged_bp": right_bp,
        "intersection_bp": shared,
        "union_bp": union,
        "jaccard": None if union == 0 else shared / union,
    }


def frip_effects(reference: float, candidate: float) -> dict[str, float]:
    if not 0 <= reference <= 1 or not 0 <= candidate <= 1:
        raise ValueError("FRiP values must be fractions in [0, 1]")
    return {
        "absolute_difference_percentage_points": (candidate - reference) * 100,
        "log2_ratio": _log2((candidate + FRIP_EPSILON) / (reference + FRIP_EPSILON)),
        "relative_change": None
        if reference == 0
        else (candidate - reference) / reference,
    }


def _log2(value: float) -> float:
    import math

    return math.log2(value)


def offpeak_fraction(frip: float) -> float:
    if not 0 <= frip <= 1:
        raise ValueError("FRiP values must be fractions in [0, 1]")
    return 1.0 - frip


def density_enrichment(
    reads_in_peaks: int,
    total_reads: int,
    peak_bp: int,
    genome_bp: int,
) -> dict[str, float | None]:
    if total_reads < reads_in_peaks or reads_in_peaks < 0 or total_reads < 0:
        raise ValueError("read counts are inconsistent")
    if peak_bp < 0 or genome_bp < 0:
        raise ValueError("base counts must be non-negative")
    offpeak_reads = total_reads - reads_in_peaks
    offpeak_bp = genome_bp - peak_bp
    signal = None if peak_bp == 0 else reads_in_peaks / peak_bp
    background = None if offpeak_bp <= 0 else offpeak_reads / offpeak_bp
    enrichment = None
    if signal is not None and background not in (None, 0):
        enrichment = signal / background
    return {
        "reads_per_peak_base": signal,
        "reads_per_offpeak_base": background,
        "enrichment": enrichment,
        "offpeak_bp": offpeak_bp,
    }


def assert_fixed_configuration(retain: dict[str, str], exclude: dict[str, str]) -> None:
    if retain.get("duplicate_analysis") != RETAIN:
        raise ValueError("retain arm duplicate_analysis is not retain")
    if exclude.get("duplicate_analysis") != EXCLUDE:
        raise ValueError("exclude arm duplicate_analysis is not exclude")
    for field in FIXED_QC_FIELDS:
        if retain.get(field) != exclude.get(field):
            raise ValueError(
                f"fixed field {field} differs: {retain.get(field)} vs {exclude.get(field)}"
            )


def load_chrom_sizes(path: Path) -> int:
    total = 0
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        _chrom, length = line.split()[:2]
        total += int(length)
    if total <= 0:
        raise ValueError(f"chrom.sizes has no positive bases: {path}")
    return total


def write_bed(path: Path, peaks: list[tuple[str, int, int]]) -> None:
    path.write_text(
        "".join(f"{chrom}\t{start}\t{end}\n" for chrom, start, end in peaks),
        encoding="utf-8",
    )


def count_common_frip(
    script: Path, sample: str, bam: Path, peaks: Path, output: Path
) -> dict[str, str]:
    completed = subprocess.run(
        [
            sys.executable,
            str(script),
            "--sample",
            sample,
            "--bam",
            str(bam),
            "--peaks",
            str(peaks),
            "--output",
            str(output),
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        raise RuntimeError(completed.stderr.strip() or "calc_frip.py failed")
    return load_qc_summary(output, sample)


def comparison_record(
    retain_qc: dict[str, str],
    exclude_qc: dict[str, str],
    retain_peaks: list[tuple[str, int, int]],
    exclude_peaks: list[tuple[str, int, int]],
    common_frip: dict[str, dict[str, str]] | None = None,
    genome_bp: int | None = None,
) -> dict[str, object]:
    assert_fixed_configuration(retain_qc, exclude_qc)
    jaccard = bp_jaccard(retain_peaks, exclude_peaks)
    common = merge_intervals(retain_peaks + exclude_peaks)
    record: dict[str, object] = {
        "scope": "single treatment/control pair; not biological acceptance",
        "reference_arm": RETAIN,
        "candidate_arm": EXCLUDE,
        "replicate_consistency": "not_applicable_single_pair",
        "peak_count": {
            "retain": len(retain_peaks),
            "exclude": len(exclude_peaks),
            "role": "descriptive_only",
        },
        "peak_width_bp": {
            "retain": peak_widths(retain_peaks),
            "exclude": peak_widths(exclude_peaks),
        },
        "overlap_fraction_at_least_one_bp": {
            "retain_peaks_hit_by_exclude": overlap_fraction(
                retain_peaks, exclude_peaks
            ),
            "exclude_peaks_hit_by_retain": overlap_fraction(
                exclude_peaks, retain_peaks
            ),
        },
        "merged_bp_jaccard": jaccard,
        "common_universe_bp": interval_bp(common),
        "analysis_templates": {
            arm: {
                "analysis_valid_fragments": qc.get("analysis_valid_fragments"),
                "stage": qc.get("analysis_input_stage"),
                "duplicate_fraction_picard": qc.get("percent_duplication"),
                "duplicate_metric_stage": "pre_duplicate_markdup_bam",
            }
            for arm, qc in ((RETAIN, retain_qc), (EXCLUDE, exclude_qc))
        },
        "library_complexity": {
            arm: {
                "nrf": qc.get("nrf"),
                "pbc1": qc.get("pbc1"),
                "pbc2": qc.get("pbc2"),
                "stage": "pre_duplicate_mapq_bam",
                "denominator": "MAPQ-filtered fragments before duplicate removal",
            }
            for arm, qc in ((RETAIN, retain_qc), (EXCLUDE, exclude_qc))
        },
        "own_peak_frip": {
            arm: {
                "frip": qc.get("frip"),
                "total_reads": qc.get("total_reads"),
                "reads_in_peaks": qc.get("reads_in_peaks"),
                "denominator": "analysis BAM read records from calc_frip.py",
                "peaks": "arm-specific blacklist-filtered narrowPeak",
                "comparable_across_arms": False,
            }
            for arm, qc in ((RETAIN, retain_qc), (EXCLUDE, exclude_qc))
        },
    }
    if common_frip is None:
        record["common_universe_frip"] = {"status": "not_measured"}
        record["signal_background"] = {"status": "not_measured"}
        return record
    retain_frip = float(common_frip[RETAIN]["frip"])
    exclude_frip = float(common_frip[EXCLUDE]["frip"])
    effects = frip_effects(retain_frip, exclude_frip)
    retain_off = offpeak_fraction(retain_frip)
    exclude_off = offpeak_fraction(exclude_frip)
    record["common_universe_frip"] = {
        "status": "measured",
        "normalization": "analysis BAM read records / common merged narrowPeak universe",
        "epsilon": FRIP_EPSILON,
        "retain": common_frip[RETAIN],
        "exclude": common_frip[EXCLUDE],
        "effects_exclude_minus_retain": effects,
        "offpeak_fraction": {RETAIN: retain_off, EXCLUDE: exclude_off},
        "offpeak_background_ratio_exclude_over_retain": (
            None if retain_off == 0 else exclude_off / retain_off
        ),
    }
    if genome_bp is None:
        record["signal_background"] = {"status": "not_measured_no_chrom_sizes"}
        return record
    peak_bp = int(record["common_universe_bp"])
    record["signal_background"] = {
        "status": "measured",
        "normalization": "analysis-BAM read records per reference base",
        "genome_bp": genome_bp,
        "retain": density_enrichment(
            int(common_frip[RETAIN]["reads_in_peaks"]),
            int(common_frip[RETAIN]["total_reads"]),
            peak_bp,
            genome_bp,
        ),
        "exclude": density_enrichment(
            int(common_frip[EXCLUDE]["reads_in_peaks"]),
            int(common_frip[EXCLUDE]["total_reads"]),
            peak_bp,
            genome_bp,
        ),
    }
    return record


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--retain-qc", type=Path, required=True)
    parser.add_argument("--exclude-qc", type=Path, required=True)
    parser.add_argument("--retain-peaks", type=Path, required=True)
    parser.add_argument("--exclude-peaks", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--chrom-sizes", type=Path)
    parser.add_argument("--retain-bam", type=Path)
    parser.add_argument("--exclude-bam", type=Path)
    parser.add_argument(
        "--frip-script",
        type=Path,
        default=Path(__file__).with_name("calc_frip.py"),
    )
    args = parser.parse_args(argv)
    retain_qc = load_qc_summary(args.retain_qc, args.sample)
    exclude_qc = load_qc_summary(args.exclude_qc, args.sample)
    retain_peaks = parse_narrow_peaks(args.retain_peaks)
    exclude_peaks = parse_narrow_peaks(args.exclude_peaks)
    common_frip = None
    if args.retain_bam or args.exclude_bam:
        if not (args.retain_bam and args.exclude_bam):
            raise SystemExit("both analysis BAMs are required for common-universe FRiP")
        universe = args.output.with_suffix(".common.bed")
        write_bed(universe, merge_intervals(retain_peaks + exclude_peaks))
        common_frip = {
            RETAIN: count_common_frip(
                args.frip_script,
                args.sample,
                args.retain_bam,
                universe,
                args.output.with_suffix(".retain.frip.tsv"),
            ),
            EXCLUDE: count_common_frip(
                args.frip_script,
                args.sample,
                args.exclude_bam,
                universe,
                args.output.with_suffix(".exclude.frip.tsv"),
            ),
        }
    genome_bp = load_chrom_sizes(args.chrom_sizes) if args.chrom_sizes else None
    record = comparison_record(
        retain_qc,
        exclude_qc,
        retain_peaks,
        exclude_peaks,
        common_frip,
        genome_bp,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
