#!/usr/bin/env python3
"""Calculate NRF/PBC library complexity metrics from a BAM file.

Stage 3c-1: BAM-derived NRF (Non-Redundant Fraction) and PBC (PCR
Bottleneck Coefficient) metrics.  Uses samtools view + stdlib counting.
Does NOT require preseq, phantompeakqualtools, or R packages.

PE fragment key (per ENCODE-style convention):
  - Only first-in-pair records from proper pairs.
  - Fragment span = (min(POS, PNEXT), min(POS, PNEXT) + abs(TLEN) - 1).
  - Falls back to (RNAME, POS, strand) when TLEN is zero / unusable.

SE fragment key:
  - (RNAME, POS, strand) from primary mapped reads.

Output columns:
    sample  total_fragments  distinct_fragments  one_read_fragments
    two_read_fragments  nrf  pbc1  pbc2

Usage:
    python3 scripts/calc_nrf_pbc.py \\
        --sample SAMPLE \\
        --bam sample.final.bam \\
        --output sample.nrf_pbc.tsv
"""

import argparse
import collections
import os
import subprocess
import sys

_NA = "NA"

# SAM flag bits
_FPAIRED = 0x1
_FPROPER_PAIR = 0x2
_FUNMAP = 0x4
_FMUNMAP = 0x8
_FREVERSE = 0x10
_FFIRST = 0x40
_FSECONDARY = 0x100
_FSUPP = 0x800

# PE: paired + proper_pair + first_in_pair
_PE_INCLUDE = _FPAIRED | _FPROPER_PAIR | _FFIRST
# Exclude: unmapped + mate_unmapped + secondary + supplementary
_PE_EXCLUDE = _FUNMAP | _FMUNMAP | _FSECONDARY | _FSUPP

# SE: exclude unmapped + secondary + supplementary
_SE_EXCLUDE = _FUNMAP | _FSECONDARY | _FSUPP

# Columns in samtools view output (zero-indexed)
_COL_FLAG = 1
_COL_RNAME = 2
_COL_POS = 3
_COL_PNEXT = 7
_COL_TLEN = 8


def _read_fragments(bam_path: str, layout: str):
    """Yield fragment keys from a BAM.

    For PE, each fragment is represented by its first-in-pair record with
    a span-aware key.  For SE, each fragment is a single alignment.
    """
    if layout == "PE":
        include = _PE_INCLUDE
        exclude = _PE_EXCLUDE
    else:
        include = None
        exclude = _SE_EXCLUDE

    cmd = ["samtools", "view"]
    if include is not None:
        cmd.extend(["-f", str(include)])
    if exclude != 0:
        cmd.extend(["-F", str(exclude)])
    cmd.append(bam_path)

    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    if proc.stdout is None:
        return

    for line in proc.stdout:
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 11:
            continue
        try:
            flag = int(cols[_COL_FLAG])
        except ValueError:
            continue

        rname = cols[_COL_RNAME]
        pos = cols[_COL_POS]

        if layout == "PE":
            yield _pe_key(cols, flag, rname, pos)
        else:
            strand = "-" if (flag & _FREVERSE) else "+"
            yield (rname, pos, strand)

    proc.stdout.close()
    rc = proc.wait()
    if rc != 0 and proc.stderr:
        err = proc.stderr.read().strip()
        if err:
            print(f"ERROR: samtools view failed: {err}", file=sys.stderr)
            sys.exit(1)


def _pe_key(cols: list[str], flag: int, rname: str, pos: str) -> tuple:
    """Return a span-aware PE fragment key.

    When TLEN is usable (non-zero), the key is (RNAME, start, end)
    where start = min(POS, PNEXT) and end = start + abs(TLEN) - 1.

    Falls back to (RNAME, POS, strand) when TLEN is zero or PNEXT is
    missing / non-numeric.
    """
    tlen_str = cols[_COL_TLEN]
    pnext_str = cols[_COL_PNEXT]

    try:
        tlen = int(tlen_str)
    except (ValueError, TypeError):
        tlen = 0
    try:
        pnext = int(pnext_str)
    except (ValueError, TypeError):
        pnext = 0

    if tlen != 0 and pnext > 0:
        pos_int = int(pos)
        start = min(pos_int, pnext)
        end = start + abs(tlen) - 1
        return (rname, str(start), str(end))

    # Fallback for unusable TLEN
    strand = "-" if (flag & _FREVERSE) else "+"
    return (rname, pos, strand)


def _compute_metrics(fragments) -> dict[str, str]:
    """Count fragment keys and compute NRF/PBC metrics.

    Returns a dict with all output columns.
    """
    counter = collections.Counter(fragments)
    total = 0
    for count in counter.values():
        total += count
    distinct = len(counter)
    one_read = sum(1 for c in counter.values() if c == 1)
    two_read = sum(1 for c in counter.values() if c == 2)

    def _ratio(num, denom):
        if denom == 0:
            return _NA
        return f"{num / denom:.6f}"

    return {
        "total_fragments": str(total),
        "distinct_fragments": str(distinct),
        "one_read_fragments": str(one_read),
        "two_read_fragments": str(two_read),
        "nrf": _ratio(distinct, total),
        "pbc1": _ratio(one_read, distinct),
        "pbc2": _ratio(one_read, two_read),
    }


def _detect_layout(bam_path: str) -> str:
    """Heuristic: stream at most 100 mapped reads to guess PE vs SE."""
    proc = subprocess.Popen(
        ["samtools", "view", "-F", str(_FUNMAP | _FSECONDARY | _FSUPP),
         bam_path],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    layout = "SE"
    try:
        count = 0
        for line in proc.stdout or []:
            if count >= 100:
                break
            count += 1
            cols = line.split("\t")
            if len(cols) >= 2:
                try:
                    flag = int(cols[1])
                except ValueError:
                    continue
                if flag & _FPAIRED:
                    layout = "PE"
                    break
    finally:
        proc.stdout.close() if proc.stdout else None
        proc.terminate()
        proc.wait()
    return layout


def main():
    parser = argparse.ArgumentParser(
        description="Calculate NRF/PBC library complexity from a BAM"
    )
    parser.add_argument("--sample", required=True, help="Sample ID")
    parser.add_argument("--bam", required=True, help="Path to BAM file")
    parser.add_argument("--output", required=True, help="Path to output TSV")
    args = parser.parse_args()

    if not os.path.isfile(args.bam):
        print(f"ERROR: BAM file not found: {args.bam}", file=sys.stderr)
        sys.exit(1)

    layout = _detect_layout(args.bam)

    metrics = _compute_metrics(_read_fragments(args.bam, layout))

    columns = [
        "sample",
        "total_fragments",
        "distinct_fragments",
        "one_read_fragments",
        "two_read_fragments",
        "nrf",
        "pbc1",
        "pbc2",
    ]

    with open(args.output, "w") as fh:
        fh.write("\t".join(columns) + "\n")
        row = [args.sample]
        for col in columns[1:]:
            row.append(metrics.get(col, _NA))
        fh.write("\t".join(row) + "\n")


if __name__ == "__main__":
    main()
