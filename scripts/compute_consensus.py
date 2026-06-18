#!/usr/bin/env python3
"""Compute replicate-validated consensus peaks from per-biorep peak calls.

Stage 54 consensus engine. Pure Python stdlib, no BEDTools dependency.

Algorithm:
  1. Parse per-biorep peak files into source peaks.
  2. Build undirected graph: edges = pairwise reciprocal overlap >= threshold.
  3. Find connected components.
  4. Filter: support_count = distinct bioreps in component,
     keep if support_count >= min_replicates.
  5. Output consensus peaks (merged intervals) + summary TSV.

Usage:
  python3 scripts/compute_consensus.py \\
    --peaks rep1.narrowPeak rep2.narrowPeak \\
    --bioreps 1 2 \\
    --format narrowPeak \\
    --output consensus.narrowPeak \\
    --summary consensus.summary.tsv \\
    --experiment exp1 --assay chipseq --caller macs3 --peak-mode narrow
"""

import argparse
import json
import os
import sys
from collections import defaultdict

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

NARROWPEAK_COLS = 10
BROADPEAK_COLS = 9
BED_MIN_COLS = 3
BED_OUTPUT_COLS = 6

VALID_FORMATS = {"narrowPeak", "broadPeak", "bed"}

HEADER_PREFIXES = ("#", "track", "browser")

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Compute replicate-validated consensus peaks (Stage 54)"
    )
    parser.add_argument(
        "--peaks", nargs="+", required=True,
        help="Per-biorep peak files",
    )
    parser.add_argument(
        "--bioreps", nargs="+", required=True,
        help="Biological replicate labels (one per --peaks file)",
    )
    parser.add_argument(
        "--format", required=True,
        choices=sorted(VALID_FORMATS),
        help="Peak file format",
    )
    parser.add_argument(
        "--min-replicates", type=int, default=2,
        help="Minimum supporting biological replicates (default: 2)",
    )
    parser.add_argument(
        "--reciprocal-overlap", type=float, default=0.5,
        help="Reciprocal overlap threshold in (0, 1] (default: 0.5)",
    )
    parser.add_argument(
        "--output", required=True,
        help="Output consensus peak file path",
    )
    parser.add_argument(
        "--summary", required=True,
        help="Output summary TSV path",
    )
    parser.add_argument(
        "--experiment", required=True,
        help="Experiment identifier",
    )
    parser.add_argument(
        "--assay", required=True,
        help="Assay type (e.g. chipseq, cuttag, atac)",
    )
    parser.add_argument(
        "--caller", required=True,
        help="Peak caller (e.g. macs3, seacr)",
    )
    parser.add_argument(
        "--peak-mode", required=True,
        help="Peak mode (e.g. narrow, broad, stringent)",
    )
    parser.add_argument(
        "--final-method", default="",
        help="Primary final reproducibility method for this mode",
    )
    parser.add_argument(
        "--final-output", default="",
        help="Path to primary final validated peak file",
    )
    return parser.parse_args(argv)


# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------


def validate_inputs(args):
    """Validate CLI arguments. Exits with message on failure."""
    errors = []

    n_peaks = len(args.peaks)
    n_bioreps = len(args.bioreps)

    if n_peaks < 2:
        errors.append(
            "at least two peak files / biological replicates are required"
        )
    if n_peaks != n_bioreps:
        errors.append(
            f"{n_peaks} peak file(s) but {n_bioreps} biorep label(s)"
        )

    seen_bioreps = set()
    for label in args.bioreps:
        if label in seen_bioreps:
            errors.append(f"duplicate biorep label: {label!r}")
        seen_bioreps.add(label)

    if args.min_replicates < 2:
        errors.append("min-replicates must be >= 2")
    elif args.min_replicates > n_bioreps:
        errors.append(
            f"min-replicates ({args.min_replicates}) exceeds number of "
            f"biological replicates ({n_bioreps})"
        )

    if not (0 < args.reciprocal_overlap <= 1):
        errors.append(
            f"reciprocal-overlap must be in (0, 1], "
            f"got {args.reciprocal_overlap}"
        )

    for path in args.peaks:
        if not os.path.isfile(path):
            errors.append(f"cannot read: {path}")

    if errors:
        for msg in errors:
            print(f"ERROR: {msg}", file=sys.stderr)
        sys.exit(1)


# ---------------------------------------------------------------------------
# Peak parsing
# ---------------------------------------------------------------------------


def parse_peak_file(path, biorep_label, fmt):
    """Parse a peak file and return list of source peak tuples.

    Each tuple: (chrom, start, end, biorep_label, fields_dict)

    fields_dict keys vary by format:
      narrowPeak: name, score, strand, signalValue, pValue, qValue, peak
      broadPeak:  name, score, strand, signalValue, pValue, qValue
      bed:        name, score, strand (synthesized defaults if absent)

    Raises SystemExit on malformed lines.
    """
    peaks = []
    with open(path) as fh:
        for line_no, raw_line in enumerate(fh, start=1):
            line = raw_line.rstrip("\n")

            # Skip headers and blank lines
            if not line or line.startswith(HEADER_PREFIXES):
                continue

            fields = line.split("\t")

            try:
                peak = _parse_peak_line(fields, fmt, path, line_no)
            except ValueError as e:
                print(
                    f"ERROR: {path}: line {line_no}: {e}",
                    file=sys.stderr,
                )
                sys.exit(1)

            chrom, start, end, parsed_fields = peak
            peaks.append((chrom, start, end, biorep_label, parsed_fields))

    return peaks


def _parse_peak_line(fields, fmt, path, line_no):
    """Parse a single peak line. Returns (chrom, start, end, fields_dict).

    Raises ValueError for malformed input.
    """
    if fmt == "narrowPeak":
        if len(fields) < NARROWPEAK_COLS:
            raise ValueError(
                f"expected {NARROWPEAK_COLS} narrowPeak columns, "
                f"got {len(fields)}"
            )
        chrom = fields[0]
        start = _parse_coord(fields[1], "start")
        end = _parse_coord(fields[2], "end")
        _check_interval(start, end)
        score = _parse_int_safe(fields[4], "score")
        strand_val = fields[5] if len(fields) > 5 else "."
        signal_value = _parse_float_safe(fields[6], "signalValue")
        pvalue = fields[7] if len(fields) > 7 else "-1"
        qvalue = fields[8] if len(fields) > 8 else "-1"
        peak_val = _parse_int_or_none(fields[9])
        return chrom, start, end, {
            "name": fields[3] if len(fields) > 3 else ".",
            "score": score,
            "strand": strand_val,
            "signalValue": signal_value,
            "pValue": pvalue,
            "qValue": qvalue,
            "peak": peak_val,
        }

    elif fmt == "broadPeak":
        if len(fields) < BROADPEAK_COLS:
            raise ValueError(
                f"expected {BROADPEAK_COLS} broadPeak columns, "
                f"got {len(fields)}"
            )
        chrom = fields[0]
        start = _parse_coord(fields[1], "start")
        end = _parse_coord(fields[2], "end")
        _check_interval(start, end)
        score = _parse_int_safe(fields[4], "score")
        strand_val = fields[5] if len(fields) > 5 else "."
        signal_value = _parse_float_safe(fields[6], "signalValue")
        pvalue = fields[7] if len(fields) > 7 else "-1"
        qvalue = fields[8] if len(fields) > 8 else "-1"
        return chrom, start, end, {
            "name": fields[3] if len(fields) > 3 else ".",
            "score": score,
            "strand": strand_val,
            "signalValue": signal_value,
            "pValue": pvalue,
            "qValue": qvalue,
        }

    elif fmt == "bed":
        if len(fields) < BED_MIN_COLS:
            raise ValueError(
                f"expected at least {BED_MIN_COLS} BED columns, "
                f"got {len(fields)}"
            )
        chrom = fields[0]
        start = _parse_coord(fields[1], "start")
        end = _parse_coord(fields[2], "end")
        _check_interval(start, end)
        name = fields[3] if len(fields) > 3 else "."
        score = _parse_int_safe(fields[4], "score") if len(fields) > 4 else 0
        strand_val = fields[5] if len(fields) > 5 else "."
        signal_value = _parse_float_safe(
            fields[6], "signalValue"
        ) if len(fields) > 6 else 0.0
        return chrom, start, end, {
            "name": name,
            "score": score,
            "strand": strand_val,
            "signalValue": signal_value,
        }

    else:
        raise ValueError(f"unknown format: {fmt}")


def _parse_coord(value, name):
    """Parse a coordinate field as non-negative integer."""
    try:
        val = int(value)
    except (ValueError, TypeError):
        raise ValueError(f"{name} must be an integer, got {value!r}")
    if val < 0:
        raise ValueError(f"{name} must be non-negative, got {val}")
    return val


def _check_interval(start, end):
    """Raise ValueError if end <= start."""
    if end <= start:
        raise ValueError(
            f"end ({end}) must be greater than start ({start})"
        )


def _parse_int_safe(value, name):
    """Parse int or return 0."""
    try:
        return int(value)
    except (ValueError, TypeError):
        return 0


def _parse_int_or_none(value):
    """Parse int or return None for fields whose invalid value has meaning."""
    try:
        return int(value)
    except (ValueError, TypeError):
        return None


def _parse_float_safe(value, name):
    """Parse float or return 0.0."""
    try:
        return float(value)
    except (ValueError, TypeError):
        return 0.0


# ---------------------------------------------------------------------------
# Overlap graph
# ---------------------------------------------------------------------------


def build_overlap_graph(peaks, reciprocal_overlap):
    """Build adjacency list for peaks overlapping by reciprocal_overlap.

    peaks: list of (chrom, start, end, biorep_label, fields_dict)
    reciprocal_overlap: float in (0, 1]

    Returns adjacency dict: node_idx -> set of adjacent node indices.
    """
    n = len(peaks)
    graph = defaultdict(set)

    # Sort by (chrom, start) for efficient scanning
    indexed = sorted(enumerate(peaks), key=lambda x: (x[1][0], x[1][1]))

    for i in range(n):
        idx_a, (chrom_a, start_a, end_a, _, _) = indexed[i]
        len_a = end_a - start_a

        # Scan forward until other.start > end_a
        for j in range(i + 1, n):
            idx_b, (chrom_b, start_b, end_b, _, _) = indexed[j]
            if chrom_b != chrom_a:
                break
            if start_b > end_a:
                break

            len_b = end_b - start_b
            overlap = min(end_a, end_b) - max(start_a, start_b)
            if overlap <= 0:
                continue

            # Reciprocal overlap check
            if (overlap / len_a >= reciprocal_overlap and
                    overlap / len_b >= reciprocal_overlap):
                graph[idx_a].add(idx_b)
                graph[idx_b].add(idx_a)

    return graph


# ---------------------------------------------------------------------------
# Connected components
# ---------------------------------------------------------------------------


def find_components(graph, n_nodes):
    """Find connected components via DFS.

    Returns list of lists of node indices.
    Nodes not in graph are isolated (singleton components).
    """
    visited = [False] * n_nodes
    components = []

    for node in range(n_nodes):
        if visited[node]:
            continue
        # DFS
        stack = [node]
        visited[node] = True
        component = []
        while stack:
            current = stack.pop()
            component.append(current)
            for neighbor in graph.get(current, ()):
                if not visited[neighbor]:
                    visited[neighbor] = True
                    stack.append(neighbor)
        components.append(component)

    return components


# ---------------------------------------------------------------------------
# Filtering
# ---------------------------------------------------------------------------


def filter_by_support(components, peaks, min_replicates):
    """Filter components by distinct biorep support count.

    Returns list of (component_nodes, support_count, distinct_bioreps).
    """
    retained = []
    for comp in components:
        bioreps = set()
        for node_idx in comp:
            _, _, _, biorep_label, _ = peaks[node_idx]
            bioreps.add(biorep_label)
        support_count = len(bioreps)
        if support_count >= min_replicates:
            retained.append((comp, support_count, sorted(bioreps)))
    return retained


# ---------------------------------------------------------------------------
# Consensus peak construction
# ---------------------------------------------------------------------------


def compute_consensus_peaks(retained_components, peaks, n_bioreps, fmt):
    """Build consensus peak rows from retained components.

    Returns list of tab-separated strings.
    """
    # Sort components by consensus coordinates, then deterministic source order.
    def component_key(item):
        comp_nodes, _, bioreps = item
        chroms = set()
        min_start = float("inf")
        max_end = -1
        first_source_order = min(comp_nodes)
        for node_idx in comp_nodes:
            chrom, start, end, _, _ = peaks[node_idx]
            chroms.add(chrom)
            if start < min_start:
                min_start = start
            if end > max_end:
                max_end = end
        return (
            sorted(chroms)[0],
            min_start,
            max_end,
            tuple(bioreps),
            first_source_order,
        )

    sorted_components = sorted(retained_components, key=component_key)

    rows = []
    for peak_idx, (comp_nodes, support_count, _supporting_bioreps) in enumerate(
        sorted_components, start=1
    ):
        name = f"consensus_peak_{peak_idx}"
        row = _build_consensus_row(
            comp_nodes, peaks, support_count, n_bioreps, name, fmt
        )
        rows.append(row)

    return rows


def _build_consensus_row(comp_nodes, peaks, support_count, n_bioreps, name, fmt):
    """Build a single consensus peak row string."""
    chrom = None
    min_start = float("inf")
    max_end = -1

    for node_idx in comp_nodes:
        c, start, end, _, fields = peaks[node_idx]
        if chrom is None:
            chrom = c
        if start < min_start:
            min_start = start
        if end > max_end:
            max_end = end

    score = int(1000 * support_count / n_bioreps + 0.5)
    strand = "."

    # Max signalValue
    max_signal = 0.0
    best_signal_peak = None
    for node_idx in comp_nodes:
        _, start, end, biorep, fields = peaks[node_idx]
        sv = fields.get("signalValue", 0.0)
        if isinstance(sv, str):
            try:
                sv = float(sv)
            except (ValueError, TypeError):
                sv = 0.0
        if sv > max_signal:
            max_signal = sv
            best_signal_peak = (start, end, biorep, fields, node_idx)
        elif sv == max_signal and best_signal_peak is not None:
            # Tie-break: smallest (chrom, start, end), then biorep, then index
            old_start, old_end, old_biorep, _, old_idx = best_signal_peak
            if (start < old_start or
                    (start == old_start and end < old_end) or
                    (start == old_start and end == old_end and biorep < old_biorep) or
                    (start == old_start and end == old_end and biorep == old_biorep
                     and node_idx < old_idx)):
                best_signal_peak = (start, end, biorep, fields, node_idx)

    signal_value_str = str(max_signal) if max_signal != int(max_signal) else str(int(max_signal))

    if chrom is None:
        chrom = "."

    if fmt == "narrowPeak":
        summit = _derive_summit(best_signal_peak, min_start, max_end)
        return "\t".join([
            str(chrom), str(int(min_start)), str(int(max_end)), str(name), str(score),
            str(strand), str(signal_value_str), "-1", "-1", str(summit),
        ])

    elif fmt == "broadPeak":
        return "\t".join([
            str(chrom), str(int(min_start)), str(int(max_end)), str(name), str(score),
            str(strand), str(signal_value_str), "-1", "-1",
        ])

    elif fmt == "bed":
        return "\t".join([
            str(chrom), str(int(min_start)), str(int(max_end)), str(name), str(score), str(strand),
        ])

    return ""


def _derive_summit(best_signal_peak, consensus_start, consensus_end):
    """Derive consensus summit from the max-signal source peak.

    Falls back to midpoint if source summit is invalid.
    """
    if best_signal_peak is None:
        return (consensus_end - consensus_start) // 2

    source_start, source_end, _biorep, fields, _idx = best_signal_peak
    peak_val = fields.get("peak")

    if peak_val is not None and isinstance(peak_val, (int, float)) and peak_val >= 0:
        peak_val = int(peak_val)
        summit_abs = source_start + peak_val
        if summit_abs < source_end:
            consensus_summit = summit_abs - consensus_start
            # Clamp to valid range
            span = consensus_end - consensus_start
            if consensus_summit < 0:
                consensus_summit = 0
            elif consensus_summit >= span:
                consensus_summit = span - 1
            return consensus_summit

    # Fallback: midpoint
    return (consensus_end - consensus_start) // 2


# ---------------------------------------------------------------------------
# Summary TSV
# ---------------------------------------------------------------------------


def build_support_distribution(retained_components):
    """Build support_distribution dict from retained components."""
    dist = defaultdict(int)
    for _, support_count, _ in retained_components:
        dist[str(support_count)] += 1
    return dict(sorted(dist.items(), key=lambda x: int(x[0])))


def write_summary_tsv(
    path, args, n_bioreps, retained_components, support_distribution, output_file
):
    """Write consensus summary TSV."""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)

    consensus_peak_count = len(retained_components)

    header = [
        "experiment", "assay", "peak_mode", "caller", "n_bioreps",
        "min_replicates", "reciprocal_overlap", "consensus_peak_count",
        "support_distribution", "biorep_labels", "source_peak_files",
        "final_method", "final_output",
    ]

    biorep_labels_str = ",".join(sorted(args.bioreps))
    source_peak_files_json = json.dumps([os.path.abspath(p) for p in args.peaks])
    support_dist_json = json.dumps(support_distribution)

    row = [
        args.experiment, args.assay, args.peak_mode, args.caller,
        str(n_bioreps), str(args.min_replicates),
        str(args.reciprocal_overlap), str(consensus_peak_count),
        support_dist_json, biorep_labels_str, source_peak_files_json,
        args.final_method, args.final_output,
    ]

    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        fh.write("\t".join(row) + "\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main(argv=None):
    args = parse_args(argv)

    # 1. Validate inputs
    validate_inputs(args)
    n_bioreps = len(args.bioreps)

    # 2. Parse all peak files
    all_peaks = []
    for path, biorep in zip(args.peaks, args.bioreps):
        file_peaks = parse_peak_file(path, biorep, args.format)
        all_peaks.extend(file_peaks)

    # 3. Build overlap graph
    graph = build_overlap_graph(all_peaks, args.reciprocal_overlap)

    # 4. Find connected components
    components = find_components(graph, len(all_peaks))

    # 5. Filter by support count
    retained = filter_by_support(components, all_peaks, args.min_replicates)

    # 6. Compute consensus peaks
    consensus_rows = compute_consensus_peaks(
        retained, all_peaks, n_bioreps, args.format
    )

    # 7. Write output peak file
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, "w") as fh:
        for row in consensus_rows:
            fh.write(row + "\n")
        # If no rows, file is empty (0 bytes) — deliberate

    # 8. Write summary TSV
    support_dist = build_support_distribution(retained)
    write_summary_tsv(
        args.summary, args, n_bioreps, retained, support_dist, args.output
    )

    return 0


if __name__ == "__main__":
    sys.exit(main())
