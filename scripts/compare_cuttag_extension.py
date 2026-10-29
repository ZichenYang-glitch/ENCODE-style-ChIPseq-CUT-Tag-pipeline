#!/usr/bin/env python3
"""Replicate and max_insert metrics for the CUT&Tag extension.

Interval merge, base Jaccard and FRiP effects come from
``compare_cuttag_duplicate_arms.py``. This script does not freeze a default.
``auto`` is not an arm here.

max_insert templates are primary proper pairs from the aligner BAM. A 500 bp
arm is not a filter of a 700 bp BAM. Newly admitted counts and the
high-confidence short-template fraction come from
``compare_cuttag_insert_stability.py``: MAPQ is the minimum of the mates,
length bins are that module's bins, and a template missing from the 700 arm
stays in the denominator without counting as a coordinate change.

Length-binned coverage correlation uses 100 kb windows on the occupied union
only. Shared empty windows are not treated as quality evidence.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
import subprocess
from pathlib import Path
from typing import NamedTuple


WINDOW_BP = 100_000
FRIP_RELATIVE_DROP = 0.10
OFFPEAK_RELATIVE_RISE = 0.05
SHORT_TEMPLATE_CHANGE_LIMIT = 0.01
COVERAGE_LENGTH_BINS: tuple[tuple[str, int, int | None], ...] = (
    ("1-100", 1, 100),
    ("101-200", 101, 200),
    ("201-300", 201, 300),
    ("301-500", 301, 500),
    ("501-700", 501, 700),
    ("701+", 701, None),
)
_PROPER = 0x2
_UNMAPPED = 0x4
_SECONDARY = 0x100
_SUPPLEMENTARY = 0x800
_READ1 = 0x40
_READ2 = 0x80


class Mate(NamedTuple):
    chrom: str
    pos: int
    mapq: int


class Template(NamedTuple):
    qname: str
    r1: Mate
    r2: Mate
    length: int


def _duplicate_module():
    path = Path(__file__).with_name("compare_cuttag_duplicate_arms.py")
    spec = importlib.util.spec_from_file_location("compare_cuttag_duplicate_arms", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _insert_module():
    path = Path(__file__).with_name("compare_cuttag_insert_stability.py")
    spec = importlib.util.spec_from_file_location(
        "compare_cuttag_insert_stability", path
    )
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def coverage_length_bin(length: int) -> str | None:
    if length <= 0:
        return None
    for name, low, high in COVERAGE_LENGTH_BINS:
        if high is None and length >= low:
            return name
        if high is not None and low <= length <= high:
            return name
    return None


def _flag_skipped(flag: int) -> bool:
    return bool(flag & (_UNMAPPED | _SECONDARY | _SUPPLEMENTARY)) or not (
        flag & _PROPER
    )


def templates_from_sam_text(text: str) -> tuple[dict[str, Template], list[str]]:
    """Return primary proper templates and qnames that never formed a pair."""

    mates: dict[str, dict[str, Mate]] = {}
    lengths: dict[str, dict[str, int]] = {}
    templates: dict[str, Template] = {}
    for line_number, line in enumerate(text.splitlines(), 1):
        if not line or line.startswith("@"):
            continue
        cells = line.split("\t")
        if len(cells) < 11:
            raise ValueError(f"SAM line {line_number} has fewer than 11 fields")
        qname = cells[0]
        flag = int(cells[1])
        if _flag_skipped(flag):
            continue
        if qname in templates:
            raise ValueError(f"duplicate primary alignment for {qname}")
        end = "r1" if flag & _READ1 else "r2" if flag & _READ2 else ""
        if not end:
            raise ValueError(f"SAM line {line_number} is proper but not R1 or R2")
        mate_slot = mates.setdefault(qname, {})
        length_slot = lengths.setdefault(qname, {})
        if end in mate_slot:
            raise ValueError(f"duplicate primary {end} for {qname}")
        mate_slot[end] = Mate(cells[2], int(cells[3]), int(cells[4]))
        length_slot[end] = abs(int(cells[8]))
        if "r1" in mate_slot and "r2" in mate_slot:
            r1_length = length_slot["r1"]
            r2_length = length_slot["r2"]
            if r1_length and r2_length and r1_length != r2_length:
                raise ValueError(f"{qname} has disagreeing template lengths")
            templates[qname] = Template(
                qname,
                mate_slot["r1"],
                mate_slot["r2"],
                r1_length or r2_length,
            )
            del mates[qname]
            del lengths[qname]
    return templates, sorted(mates)


def load_templates(path: Path) -> tuple[dict[str, Template], list[str]]:
    if path.suffix == ".sam":
        return templates_from_sam_text(path.read_text(encoding="utf-8"))
    completed = subprocess.run(
        ["samtools", "view", "-f", "2", "-F", "2304", str(path)],
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        raise RuntimeError(completed.stderr.strip() or f"samtools failed for {path}")
    return templates_from_sam_text(completed.stdout)


def _stability_record(template: Template) -> dict[str, object]:
    return {
        "qname": template.qname,
        "r1_chrom": template.r1.chrom,
        "r1_pos": template.r1.pos,
        "r2_chrom": template.r2.chrom,
        "r2_pos": template.r2.pos,
        "length": template.length,
        "mapq": min(template.r1.mapq, template.r2.mapq),
    }


def max_insert_template_report(
    arm_500: dict[str, Template],
    arm_700: dict[str, Template],
    incomplete_500: list[str],
    incomplete_700: list[str],
) -> dict[str, object]:
    stability = _insert_module()
    contrast = stability.comparison_record(
        [_stability_record(template) for template in arm_500.values()],
        [_stability_record(template) for template in arm_700.values()],
    )
    fraction = contrast["short_template_stability"]["fraction"]
    return {
        "stage": "primary proper pairs from the aligner BAM or SAM",
        "not_a_bam_filter": True,
        "metric_source": "compare_cuttag_insert_stability.comparison_record",
        "templates_500": len(arm_500),
        "templates_700": len(arm_700),
        "incomplete_500": len(incomplete_500),
        "incomplete_700": len(incomplete_700),
        "qnames_only_in_500": len(set(arm_500) - set(arm_700)),
        "insert_contrast": contrast,
        "short_template_trigger": {
            "limit_exclusive": SHORT_TEMPLATE_CHANGE_LIMIT,
            "definition": (
                "changed high-confidence short templates divided by all "
                "such templates in the 500 arm; absence is not a coordinate change"
            ),
            "triggered": bool(
                fraction is not None and fraction > SHORT_TEMPLATE_CHANGE_LIMIT
            ),
            "status": "not_comparable" if fraction is None else "measured",
        },
    }


def relative_change(candidate: float, other: float) -> float | None:
    if other == 0:
        return None
    return (candidate - other) / other


def signal_regression_triggered(candidate_frip: float, other_frip: float) -> dict:
    """True only when both pre-registered inequalities are strict."""

    duplicate = _duplicate_module()
    candidate_off = duplicate.offpeak_fraction(candidate_frip)
    other_off = duplicate.offpeak_fraction(other_frip)
    frip_change = relative_change(candidate_frip, other_frip)
    off_change = relative_change(candidate_off, other_off)
    comparable = frip_change is not None and off_change is not None
    triggered = bool(
        comparable
        and frip_change < -FRIP_RELATIVE_DROP
        and off_change > OFFPEAK_RELATIVE_RISE
    )
    return {
        "candidate_frip": candidate_frip,
        "other_frip": other_frip,
        "relative_frip_change": frip_change,
        "relative_offpeak_change": off_change,
        "triggered": triggered,
        "status": "measured" if comparable else "not_comparable",
    }


def peak_jaccard(
    left: list[tuple[str, int, int]], right: list[tuple[str, int, int]]
) -> dict[str, float | int | None]:
    duplicate = _duplicate_module()
    raw = duplicate.bp_jaccard(left, right)
    return {
        "left_merged_bp": raw["retain_merged_bp"],
        "right_merged_bp": raw["exclude_merged_bp"],
        "intersection_bp": raw["intersection_bp"],
        "union_bp": raw["union_bp"],
        "jaccard": raw["jaccard"],
        "left_peaks_touched": duplicate.overlap_fraction(left, right),
        "right_peaks_touched": duplicate.overlap_fraction(right, left),
    }


def pearson(left: list[float], right: list[float]) -> float | None:
    if len(left) != len(right) or len(left) < 3:
        return None
    mean_left = math.fsum(left) / len(left)
    mean_right = math.fsum(right) / len(right)
    diff_left = [value - mean_left for value in left]
    diff_right = [value - mean_right for value in right]
    numerator = math.fsum(
        one * other for one, other in zip(diff_left, diff_right, strict=True)
    )
    left_scale = math.sqrt(math.fsum(value * value for value in diff_left))
    right_scale = math.sqrt(math.fsum(value * value for value in diff_right))
    if left_scale == 0 or right_scale == 0:
        return None
    return numerator / (left_scale * right_scale)


def _window_key(template: Template, window_bp: int) -> tuple[str, int] | None:
    if template.r1.chrom != template.r2.chrom or template.length <= 0:
        return None
    start = min(template.r1.pos, template.r2.pos) - 1
    if start < 0:
        return None
    return template.r1.chrom, start // window_bp


def length_binned_coverage(
    left: dict[str, Template],
    right: dict[str, Template],
    window_bp: int = WINDOW_BP,
) -> dict[str, object]:
    """Pearson correlation on windows where either replicate has a template."""

    def counts(templates: dict[str, Template]) -> dict[tuple[str, str, int], int]:
        counted: dict[tuple[str, str, int], int] = {}
        for template in templates.values():
            bin_name = coverage_length_bin(template.length)
            window = _window_key(template, window_bp)
            if bin_name is None or window is None:
                continue
            key = (bin_name, window[0], window[1])
            counted[key] = counted.get(key, 0) + 1
        return counted

    left_counts = counts(left)
    right_counts = counts(right)
    bins: dict[str, dict[str, float | int | None | bool]] = {}
    for name, _low, _high in COVERAGE_LENGTH_BINS:
        left_keys = {key[1:] for key in left_counts if key[0] == name}
        right_keys = {key[1:] for key in right_counts if key[0] == name}
        keys = left_keys | right_keys
        left_values = [float(left_counts.get((name, *key), 0)) for key in sorted(keys)]
        right_values = [
            float(right_counts.get((name, *key), 0)) for key in sorted(keys)
        ]
        bins[name] = {
            "occupied_windows": len(keys),
            "empty_windows_included": False,
            "pearson": pearson(left_values, right_values),
            "quality_evidence": False,
        }
    return {"window_bp": window_bp, "bins": bins}


def assert_same_policy(left: dict[str, str], right: dict[str, str]) -> None:
    duplicate = _duplicate_module()
    for field in ("duplicate_analysis", *duplicate.FIXED_QC_FIELDS):
        if field not in left or field not in right:
            raise ValueError(f"fixed field {field} is missing")
        if left[field] != right[field]:
            raise ValueError(
                f"fixed field {field} differs: {left[field]} vs {right[field]}"
            )


def assert_single_variable_max_insert(
    arm_500: dict[str, str], arm_700: dict[str, str]
) -> None:
    """Require retain-or-same duplicate policy and max_insert 500 versus 700."""

    duplicate = _duplicate_module()
    if arm_500.get("duplicate_analysis") != arm_700.get("duplicate_analysis"):
        raise ValueError("duplicate_analysis differs between max_insert arms")
    if arm_500.get("resolved_max_insert") != "500":
        raise ValueError("500 arm resolved_max_insert is not 500")
    if arm_700.get("resolved_max_insert") != "700":
        raise ValueError("700 arm resolved_max_insert is not 700")
    for field in duplicate.FIXED_QC_FIELDS:
        if field == "resolved_max_insert":
            continue
        if field not in arm_500 or field not in arm_700:
            raise ValueError(f"fixed field {field} is missing")
        if arm_500[field] != arm_700[field]:
            raise ValueError(
                f"fixed field {field} differs: {arm_500[field]} vs {arm_700[field]}"
            )


def evaluate_freeze(
    *,
    hard_failure: bool,
    signal_triggered_on_every_independent_replicate: bool | None,
    independent_replicate_count: int,
    confirmation_evaluated: bool,
    short_template_triggered: bool | None,
) -> dict[str, object]:
    """Apply the pre-registered freeze gate. Missing evidence stays closed."""

    reasons: list[str] = []
    if hard_failure:
        reasons.append("hard_failure")
    if independent_replicate_count < 2:
        reasons.append("fewer_than_two_independent_replicates")
    elif signal_triggered_on_every_independent_replicate is None:
        reasons.append("signal_trigger_not_evaluated")
    elif signal_triggered_on_every_independent_replicate:
        reasons.append("signal_background_regression")
    if not confirmation_evaluated:
        reasons.append("confirmation_not_evaluated")
    if short_template_triggered:
        reasons.append("short_template_coordinate_change")
    return {"freeze_default": not reasons, "blocking_reasons": reasons}


def _write_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def _selected_alignment(sam: Path | None, bam: Path | None, label: str) -> Path:
    selected = [path for path in (sam, bam) if path is not None]
    if len(selected) != 1:
        raise SystemExit(f"provide exactly one SAM or BAM for the {label} arm")
    return selected[0]


def _signal_main(args: argparse.Namespace) -> int:
    duplicate = _duplicate_module()
    arm_500 = duplicate.load_qc_summary(args.qc_500, args.sample)
    arm_700 = duplicate.load_qc_summary(args.qc_700, args.sample)
    assert_single_variable_max_insert(arm_500, arm_700)
    peaks_500 = duplicate.parse_narrow_peaks(args.peaks_500)
    peaks_700 = duplicate.parse_narrow_peaks(args.peaks_700)
    payload: dict[str, object] = {
        "sample": args.sample,
        "candidate_arm": "max_insert_700",
        "other_arm": "max_insert_500",
        "duplicate_analysis": arm_700["duplicate_analysis"],
        "peak_jaccard": peak_jaccard(peaks_500, peaks_700),
        "peak_count": {
            "max_insert_500": len(peaks_500),
            "max_insert_700": len(peaks_700),
        },
    }
    if args.bam_500 or args.bam_700:
        if not (args.bam_500 and args.bam_700 and args.chrom_sizes):
            raise SystemExit("common-universe FRiP needs both BAMs and chrom.sizes")
        universe = args.output.with_suffix(".common.bed")
        merged = duplicate.merge_intervals(peaks_500 + peaks_700)
        duplicate.write_bed(universe, merged)
        frip_500 = duplicate.count_common_frip(
            args.frip_script,
            args.sample,
            args.bam_500,
            universe,
            args.output.with_suffix(".500.frip.tsv"),
        )
        frip_700 = duplicate.count_common_frip(
            args.frip_script,
            args.sample,
            args.bam_700,
            universe,
            args.output.with_suffix(".700.frip.tsv"),
        )
        value_500 = float(frip_500["frip"])
        value_700 = float(frip_700["frip"])
        payload["common_universe_frip"] = {
            "epsilon": duplicate.FRIP_EPSILON,
            "max_insert_500": frip_500,
            "max_insert_700": frip_700,
            "effects_500_versus_700": duplicate.frip_effects(value_700, value_500),
            "offpeak_background_ratio_500_over_700": (
                None
                if duplicate.offpeak_fraction(value_700) == 0
                else duplicate.offpeak_fraction(value_500)
                / duplicate.offpeak_fraction(value_700)
            ),
        }
        payload["candidate_regression"] = signal_regression_triggered(
            value_700, value_500
        )
        genome_bp = duplicate.load_chrom_sizes(args.chrom_sizes)
        peak_bp = duplicate.interval_bp(merged)
        payload["signal_background"] = {
            "genome_bp": genome_bp,
            "max_insert_500": duplicate.density_enrichment(
                int(frip_500["reads_in_peaks"]),
                int(frip_500["total_reads"]),
                peak_bp,
                genome_bp,
            ),
            "max_insert_700": duplicate.density_enrichment(
                int(frip_700["reads_in_peaks"]),
                int(frip_700["total_reads"]),
                peak_bp,
                genome_bp,
            ),
        }
    _write_json(args.output, payload)
    return 0


def _templates_main(args: argparse.Namespace) -> int:
    path_500 = _selected_alignment(args.sam_500, args.bam_500, "500")
    path_700 = _selected_alignment(args.sam_700, args.bam_700, "700")
    arm_500, incomplete_500 = load_templates(path_500)
    arm_700, incomplete_700 = load_templates(path_700)
    _write_json(
        args.output,
        max_insert_template_report(arm_500, arm_700, incomplete_500, incomplete_700),
    )
    return 0


def _replicates_main(args: argparse.Namespace) -> int:
    duplicate = _duplicate_module()
    if args.left_qc or args.right_qc:
        if not (
            args.left_qc and args.right_qc and args.left_sample and args.right_sample
        ):
            raise SystemExit(
                "replicate QC comparison needs both QC files and both sample names"
            )
        assert_same_policy(
            duplicate.load_qc_summary(args.left_qc, args.left_sample),
            duplicate.load_qc_summary(args.right_qc, args.right_sample),
        )
    left_peaks = duplicate.parse_narrow_peaks(args.left_peaks)
    right_peaks = duplicate.parse_narrow_peaks(args.right_peaks)
    payload: dict[str, object] = {
        "left_label": args.left_label,
        "right_label": args.right_label,
        "biological_independence": "not_established_by_this_script",
        "peak_jaccard": peak_jaccard(left_peaks, right_peaks),
        "peak_count": {"left": len(left_peaks), "right": len(right_peaks)},
    }
    if args.left_bam or args.right_bam:
        if not (args.left_bam and args.right_bam):
            raise SystemExit("provide both replicate BAMs for coverage correlation")
        left_templates, left_incomplete = load_templates(args.left_bam)
        right_templates, right_incomplete = load_templates(args.right_bam)
        payload["coverage"] = length_binned_coverage(left_templates, right_templates)
        payload["incomplete_templates"] = {
            "left": len(left_incomplete),
            "right": len(right_incomplete),
        }
    _write_json(args.output, payload)
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    templates = commands.add_parser("max-insert")
    templates.add_argument("--sam-500", type=Path)
    templates.add_argument("--bam-500", type=Path)
    templates.add_argument("--sam-700", type=Path)
    templates.add_argument("--bam-700", type=Path)
    templates.add_argument("--output", type=Path, required=True)
    templates.set_defaults(func=_templates_main)
    signal = commands.add_parser("signal")
    signal.add_argument("--sample", required=True)
    signal.add_argument("--qc-500", type=Path, required=True)
    signal.add_argument("--qc-700", type=Path, required=True)
    signal.add_argument("--peaks-500", type=Path, required=True)
    signal.add_argument("--peaks-700", type=Path, required=True)
    signal.add_argument("--bam-500", type=Path)
    signal.add_argument("--bam-700", type=Path)
    signal.add_argument("--chrom-sizes", type=Path)
    signal.add_argument(
        "--frip-script",
        type=Path,
        default=Path(__file__).with_name("calc_frip.py"),
    )
    signal.add_argument("--output", type=Path, required=True)
    signal.set_defaults(func=_signal_main)
    replicates = commands.add_parser("replicates")
    replicates.add_argument("--left-sample")
    replicates.add_argument("--right-sample")
    replicates.add_argument("--left-label", required=True)
    replicates.add_argument("--right-label", required=True)
    replicates.add_argument("--left-peaks", type=Path, required=True)
    replicates.add_argument("--right-peaks", type=Path, required=True)
    replicates.add_argument("--left-qc", type=Path)
    replicates.add_argument("--right-qc", type=Path)
    replicates.add_argument("--left-bam", type=Path)
    replicates.add_argument("--right-bam", type=Path)
    replicates.add_argument("--output", type=Path, required=True)
    replicates.set_defaults(func=_replicates_main)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
