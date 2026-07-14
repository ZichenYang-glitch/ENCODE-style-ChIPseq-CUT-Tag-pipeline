#!/usr/bin/env python3
"""CUT&Tag fragment-size QC from final.bam.

For PE samples, computes fragment-length statistics from properly paired
read1 alignments using abs(TLEN). For SE samples, writes NA values with
`layout_not_supported` status. Python stdlib only; samtools required on PATH
(provided by the workflow Conda environment).

Usage:
    python3 scripts/calc_cuttag_fragment_size.py \
        --sample <sample_id> \
        --bam <final.bam> \
        --layout PE|SE \
        --output <output.tsv>
"""

import argparse
import os
import statistics
import subprocess
import sys

_SAMTOOLS = os.environ.get("SAMTOOLS", "samtools")

PE_INCLUDE_FLAGS = 0x1 | 0x2  # paired + properly paired
PE_EXCLUDE_FLAGS = 0x4 | 0x100 | 0x800  # unmapped, secondary, supplementary
READ1_FLAG = 0x40  # first in pair


def _compute_pe_stats(fragments):
    n = len(fragments)
    if n == 0:
        return {
            "fragment_count": "0",
            "fragment_mean": "NA",
            "fragment_median": "NA",
            "fragment_mode": "NA",
            "fragment_min": "NA",
            "fragment_max": "NA",
            "fraction_lt_150": "NA",
            "fraction_150_300": "NA",
            "fraction_300_500": "NA",
            "fraction_ge_500": "NA",
            "fraction_lt_120": "NA",
            "status": "no_fragments",
        }

    mean_val = statistics.mean(fragments)
    median_val = statistics.median(fragments)

    freq = {}
    for v in fragments:
        freq[v] = freq.get(v, 0) + 1
    max_count = max(freq.values())
    mode_val = min(v for v, c in freq.items() if c == max_count)

    min_val = min(fragments)
    max_val = max(fragments)

    lt_150 = sum(1 for v in fragments if v < 150)
    r150_300 = sum(1 for v in fragments if 150 <= v < 300)
    r300_500 = sum(1 for v in fragments if 300 <= v < 500)
    ge_500 = sum(1 for v in fragments if v >= 500)
    lt_120 = sum(1 for v in fragments if v < 120)

    return {
        "fragment_count": str(n),
        "fragment_mean": "{:.1f}".format(mean_val),
        "fragment_median": "{:.1f}".format(median_val),
        "fragment_mode": str(mode_val),
        "fragment_min": str(min_val),
        "fragment_max": str(max_val),
        "fraction_lt_150": "{:.4f}".format(lt_150 / n),
        "fraction_150_300": "{:.4f}".format(r150_300 / n),
        "fraction_300_500": "{:.4f}".format(r300_500 / n),
        "fraction_ge_500": "{:.4f}".format(ge_500 / n),
        "fraction_lt_120": "{:.4f}".format(lt_120 / n),
        "status": "ok",
    }


def _se_row(sample_id):
    return {
        "sample": sample_id,
        "layout": "SE",
        "fragment_count": "NA",
        "fragment_mean": "NA",
        "fragment_median": "NA",
        "fragment_mode": "NA",
        "fragment_min": "NA",
        "fragment_max": "NA",
        "fraction_lt_150": "NA",
        "fraction_150_300": "NA",
        "fraction_300_500": "NA",
        "fraction_ge_500": "NA",
        "fraction_lt_120": "NA",
        "status": "layout_not_supported",
    }


def main():
    parser = argparse.ArgumentParser(description="CUT&Tag fragment-size QC")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--layout", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    if args.layout == "SE":
        print(
            "CUT&Tag fragment-size QC: SE layout not supported; "
            "insert size cannot be reliably computed from SE alignments. "
            "Output will contain NA values.",
            file=sys.stderr,
        )
        row = _se_row(args.sample)
    else:
        try:
            subprocess.run([_SAMTOOLS, "--version"], capture_output=True, check=True)
        except (subprocess.CalledProcessError, FileNotFoundError):
            print("ERROR: samtools not found on PATH", file=sys.stderr)
            sys.exit(1)

        view_proc = subprocess.Popen(
            [_SAMTOOLS, "view", args.bam],
            stdout=subprocess.PIPE,
        )

        fragments = []
        for line in view_proc.stdout:
            line_str = line.decode("utf-8")
            fields = line_str.split("\t")
            if len(fields) < 10:
                continue
            try:
                flag = int(fields[1])
            except ValueError:
                continue
            if (flag & PE_INCLUDE_FLAGS) != PE_INCLUDE_FLAGS:
                continue
            if (flag & PE_EXCLUDE_FLAGS) != 0:
                continue
            if (flag & READ1_FLAG) == 0:
                continue
            try:
                tlen = int(fields[8])
            except ValueError:
                continue
            fragments.append(abs(tlen))

        view_proc.wait()
        if view_proc.returncode != 0:
            print("ERROR: samtools view failed", file=sys.stderr)
            sys.exit(1)

        stats = _compute_pe_stats(fragments)
        row = {
            "sample": args.sample,
            "layout": "PE",
            **stats,
        }

    header = [
        "sample",
        "layout",
        "fragment_count",
        "fragment_mean",
        "fragment_median",
        "fragment_mode",
        "fragment_min",
        "fragment_max",
        "fraction_lt_150",
        "fraction_150_300",
        "fraction_300_500",
        "fraction_ge_500",
        "fraction_lt_120",
        "status",
    ]

    with open(args.output, "w") as fh:
        fh.write("\t".join(header) + "\n")
        fh.write("\t".join(row[h] for h in header) + "\n")


if __name__ == "__main__":
    main()
