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
  fraction 'NA', organelle_reads 'NA' (missing is not a biological zero).
- If matching contigs exist, mapped reads > 0, and organelle reads == 0: records status 'ok',
  fraction '0.000000'.
- If organelle_mode == 'keep': does not enforce exclusion or fail on missing contigs.
"""

from __future__ import annotations

import argparse
import subprocess
import sys

from qc_common import MITO_QC_V2_COLUMNS, NA, write_exact_tsv


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


def _resolve(mapped_counts, candidate_chroms, is_explicit, mode, sample):
    if is_explicit:
        missing = [c for c in candidate_chroms if c not in mapped_counts]
        if missing:
            sys.stderr.write(
                f"Error: Sample {sample}: explicit organelle contig(s) {missing} "
                f"not found in BAM reference contigs\n"
            )
            return None
        matched_contigs = list(candidate_chroms)
    else:
        matched_contigs = [c for c in candidate_chroms if c in mapped_counts]

    if mode == "exclude" and not matched_contigs:
        sys.stderr.write(
            f"Error: Sample {sample}: organelle_mode is 'exclude' but none of "
            f"candidate contigs {candidate_chroms} exist in BAM contigs. "
            f"Configure genome_resources.<genome>.organelle_chroms or set organelle_mode: keep.\n"
        )
        return None
    return matched_contigs


def main(argv=None):
    args = parse_args(argv)

    mapped_counts = run_idxstats(args.bam)
    total_mapped = sum(mapped_counts.values())

    candidate_chroms: list[str] = []
    for item in args.chroms:
        candidate_chroms.extend(item.split())

    matched_contigs = _resolve(
        mapped_counts,
        candidate_chroms,
        args.is_explicit,
        args.mode,
        args.sample,
    )
    if matched_contigs is None:
        return 1

    if not matched_contigs:
        status = "unresolved_contigs"
        organelle_str = NA
        fraction_str = NA
    elif total_mapped == 0:
        status = "no_mapped_reads"
        organelle_str = NA
        fraction_str = NA
    else:
        organelle_reads = sum(mapped_counts[c] for c in matched_contigs)
        status = "ok"
        organelle_str = str(organelle_reads)
        fraction_str = f"{organelle_reads / total_mapped:.6f}"

    row = {
        "schema_version": "2",
        "sample": args.sample,
        "organelle_reads": organelle_str,
        "total_mapped_reads": str(total_mapped),
        "fraction": fraction_str,
        "mode": args.mode,
        "status": status,
        "input_stage": "pre_duplicate_mapq_bam",
        "matched_contigs": ",".join(matched_contigs) if matched_contigs else NA,
        "organelle_chroms_explicit": "true" if args.is_explicit else "false",
    }
    write_exact_tsv(args.output, MITO_QC_V2_COLUMNS, row)
    return 0


if __name__ == "__main__":
    sys.exit(main())
