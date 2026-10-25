#!/usr/bin/env python3
"""Calculate organelle read metrics from pre-duplicate, pre-filter MAPQ BAM.

Consumes:
- Pre-duplicate, pre-organelle, pre-blacklist MAPQ BAM ({sample}.{MAPQ_TAG}.bam)
- Resolves candidate organelle contigs against BAM contigs from samtools idxstats.

Invariants:
- If organelle_mode == 'exclude' and no candidates match: exits non-zero before analysis outputs.
- If organelle_mode == 'report' and no candidates match: records status 'unresolved_contigs',
  fraction 'NA', organelle_reads 'NA' (does not fabricate 0.0).
- If matching contigs exist but total mapped reads == 0: records status 'no_mapped_reads',
  fraction 'NA'.
- If matching contigs exist, mapped reads > 0, and organelle reads == 0: records status 'ok',
  fraction '0.000000'.
- If organelle_mode == 'keep': does not enforce exclusion or fail on missing contigs.
"""

import argparse
import os
import subprocess
import sys


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", required=True, help="Input MAPQ BAM path")
    parser.add_argument("--sample", required=True, help="Sample identifier")
    parser.add_argument(
        "--mode",
        required=True,
        choices=("auto", "exclude", "report", "keep"),
        help="Effective organelle policy mode",
    )
    parser.add_argument(
        "--chroms",
        nargs="*",
        default=[],
        help="Candidate organelle contig names",
    )
    parser.add_argument(
        "--is-explicit",
        action="store_true",
        default=False,
        help="True if organelle-chroms was explicitly configured",
    )
    parser.add_argument("--output", required=True, help="Output TSV path")
    return parser.parse_args(argv)


def run_idxstats(bam: str) -> dict[str, int]:
    """Run samtools idxstats on a BAM and return dict of contig to mapped reads."""
    process = subprocess.Popen(
        ["samtools", "idxstats", bam],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    stdout, stderr = process.communicate()
    if process.returncode != 0:
        raise RuntimeError(f"samtools idxstats failed on {bam}: {stderr.strip()}")

    mapped_counts: dict[str, int] = {}
    for line in stdout.splitlines():
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) < 4:
            continue
        contig = fields[0]
        if contig == "*":
            continue
        try:
            mapped = int(fields[2])
        except ValueError:
            continue
        mapped_counts[contig] = mapped
    return mapped_counts


def main(argv=None):
    args = parse_args(argv)

    mapped_counts = run_idxstats(args.bam)
    total_mapped = sum(mapped_counts.values())

    candidate_chroms: list[str] = []
    for item in args.chroms:
        candidate_chroms.extend(item.split())

    if args.is_explicit:
        missing = [c for c in candidate_chroms if c not in mapped_counts]
        if missing:
            sys.stderr.write(
                f"Error: Sample {args.sample}: explicit organelle contig(s) {missing} "
                f"not found in BAM reference contigs\n"
            )
            return 1
        matched_contigs = [c for c in candidate_chroms if c in mapped_counts]
    else:
        matched_contigs = [c for c in candidate_chroms if c in mapped_counts]

    # Resolve metrics and status
    if args.mode == "keep":
        if not matched_contigs:
            status = "unresolved_contigs"
            organelle_str = "NA"
            fraction_str = "NA"
        elif total_mapped == 0:
            status = "no_mapped_reads"
            organelle_str = "0"
            fraction_str = "NA"
        else:
            organelle_reads = sum(mapped_counts[c] for c in matched_contigs)
            frac = organelle_reads / total_mapped
            status = "ok"
            fraction_str = f"{frac:.6f}"
            organelle_str = str(organelle_reads)
    elif args.mode == "exclude":
        if not matched_contigs:
            sys.stderr.write(
                f"Error: Sample {args.sample}: organelle_mode is 'exclude' but none of "
                f"candidate contigs {candidate_chroms} exist in BAM contigs. "
                f"Configure genome_resources.<genome>.organelle_chroms or set organelle_mode: keep.\n"
            )
            return 1
        organelle_reads = sum(mapped_counts[c] for c in matched_contigs)
        if total_mapped == 0:
            status = "no_mapped_reads"
            organelle_str = "0"
            fraction_str = "NA"
        else:
            status = "ok"
            organelle_str = str(organelle_reads)
            frac = organelle_reads / total_mapped
            fraction_str = f"{frac:.6f}"
    elif args.mode == "report":
        if not matched_contigs:
            status = "unresolved_contigs"
            organelle_str = "NA"
            fraction_str = "NA"
        elif total_mapped == 0:
            status = "no_mapped_reads"
            organelle_str = "0"
            fraction_str = "NA"
        else:
            status = "ok"
            organelle_reads = sum(mapped_counts[c] for c in matched_contigs)
            organelle_str = str(organelle_reads)
            frac = organelle_reads / total_mapped
            fraction_str = f"{frac:.6f}"
    else:
        sys.stderr.write(f"Error: Unknown organelle mode: {args.mode!r}\n")
        return 1

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, "w", encoding="utf-8") as handle:
        handle.write(
            "sample\torganelle_reads\ttotal_mapped_reads\tfraction\tmode\tstatus\n"
        )
        handle.write(
            f"{args.sample}\t{organelle_str}\t{total_mapped}\t{fraction_str}\t{args.mode}\t{status}\n"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
