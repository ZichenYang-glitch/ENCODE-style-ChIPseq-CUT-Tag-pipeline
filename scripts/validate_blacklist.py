#!/usr/bin/env python3
"""Validate a blacklist BED against a BAM reference before sidecar filtering.

The analysis branch validates the BED inside filter_blacklist_templates.py.
The copies-mode sidecar rules use bedtools directly, so they call this helper
first to fail closed on missing files, malformed intervals, unknown contigs,
or intervals beyond the BAM @SQ length.
"""

import argparse
import sys

from filter_blacklist_templates import read_bam_header, validate_and_load_blacklist


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", required=True, help="BAM used to validate contigs")
    parser.add_argument("--blacklist", required=True, help="Blacklist BED path")
    parser.add_argument("--sample", required=True, help="Sample identifier")
    args = parser.parse_args(argv)

    sq = read_bam_header(args.bam)
    validate_and_load_blacklist(args.blacklist, sq, args.sample)
    return 0


if __name__ == "__main__":
    sys.exit(main())
