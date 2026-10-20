#!/usr/bin/env python3
"""Remove complete templates overlapping blacklist intervals from a BAM.

The updated analysis view contains only proper pairs (PE) or unpaired reads
(SE). For PE, the fragment span is ``[min(POS, PNEXT), min + abs(TLEN) - 1]``;
for SE, the alignment span is used. If any span overlaps a blacklist interval,
the entire template is removed. Non-proper/orphan PE records fail closed.

The implementation is a two-pass streaming filter:
1. scan ``samtools view`` text output and collect QNAMEs to drop;
2. stream headers and retained records into ``samtools view -b``.
"""

import argparse
import subprocess
import sys

from blacklist_filter import Blacklist


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", required=True)
    parser.add_argument("--blacklist", required=True)
    parser.add_argument("--layout", required=True, choices=("PE", "SE"))
    parser.add_argument("--output", required=True)
    parser.add_argument("--dropped-count", required=True)
    parser.add_argument("--sample", required=True)
    return parser.parse_args(argv)


def cigar_ref_length(cigar):
    """Return the reference-consuming length of a CIGAR string."""
    length = 0
    digits = ""
    for char in cigar:
        if char.isdigit():
            digits += char
            continue
        if not digits:
            continue
        count = int(digits)
        digits = ""
        if char in ("M", "D", "N", "=", "X"):
            length += count
    return length


def template_span(fields, layout):
    """Return (chrom, start0, end0, drop_reason) for one SAM record."""
    flag = int(fields[1])
    chrom = fields[2]
    pos = int(fields[3])
    cigar = fields[5]
    pnext = int(fields[7])
    tlen = int(fields[8])

    if flag & 0x4 or pos <= 0:
        return chrom, 0, 0, "unmapped"
    if layout == "SE":
        if flag & 0x1:
            return chrom, 0, 0, "paired record in SE layout"
        length = cigar_ref_length(cigar)
        if length <= 0:
            return chrom, 0, 0, "empty CIGAR span"
        return chrom, pos - 1, pos - 1 + length, None

    if not flag & 0x1:
        return chrom, 0, 0, "unpaired record in PE layout"
    if not flag & 0x2 or flag & 0x8 or tlen == 0 or pnext <= 0:
        return chrom, 0, 0, "non-proper or orphan pair"
    start = min(pos, pnext)
    length = abs(tlen)
    if length <= 0:
        return chrom, 0, 0, "non-positive TLEN"
    return chrom, start - 1, start - 1 + length, None


def scan_drop_names(bam, layout, blacklist):
    """Return (drop_names, dropped_alignments, total_alignments)."""
    process = subprocess.Popen(
        ["samtools", "view", bam],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    assert process.stdout is not None
    drop = set()
    total = 0
    for line in process.stdout:
        if not line or line.startswith("@"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 11:
            continue
        total += 1
        qname = fields[0]
        chrom, start, end, reason = template_span(fields, layout)
        if reason is not None:
            drop.add(qname)
            continue
        if blacklist.overlaps(chrom, start, end):
            drop.add(qname)
    stderr = process.stderr.read() if process.stderr is not None else ""
    returncode = process.wait()
    if returncode != 0:
        raise RuntimeError(f"samtools view failed: {stderr.strip()}")
    return drop, total


def write_filtered_bam(bam, output, drop):
    """Stream the BAM through a QNAME drop filter."""
    reader = subprocess.Popen(
        ["samtools", "view", "-h", bam],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    writer = subprocess.Popen(
        ["samtools", "view", "-b", "-o", output, "-"],
        stdin=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    assert reader.stdout is not None
    assert writer.stdin is not None
    dropped_alignments = 0
    kept_alignments = 0
    try:
        for line in reader.stdout:
            if line.startswith("@"):
                writer.stdin.write(line)
                continue
            qname = line.split("\t", 1)[0]
            if qname in drop:
                dropped_alignments += 1
            else:
                kept_alignments += 1
                writer.stdin.write(line)
    finally:
        writer.stdin.close()
    reader_stderr = reader.stderr.read() if reader.stderr is not None else ""
    writer_stderr = writer.stderr.read() if writer.stderr is not None else ""
    reader_code = reader.wait()
    writer_code = writer.wait()
    if reader_code != 0:
        raise RuntimeError(f"samtools view -h failed: {reader_stderr.strip()}")
    if writer_code != 0:
        raise RuntimeError(f"samtools view -b failed: {writer_stderr.strip()}")
    return dropped_alignments, kept_alignments


def main(argv=None):
    args = parse_args(argv)
    blacklist = Blacklist(args.blacklist)
    drop, total = scan_drop_names(args.bam, args.layout, blacklist)
    dropped_alignments, kept_alignments = write_filtered_bam(
        args.bam, args.output, drop
    )
    with open(args.dropped_count, "w", encoding="utf-8") as handle:
        handle.write(
            "sample\tlayout\tdropped_templates\tdropped_alignments\t"
            "kept_alignments\ttotal_alignments\n"
        )
        handle.write(
            f"{args.sample}\t{args.layout}\t{len(drop)}\t"
            f"{dropped_alignments}\t{kept_alignments}\t{total}\n"
        )
    print(
        f"{args.sample}: dropped {len(drop)} templates "
        f"({dropped_alignments} alignments); kept {kept_alignments}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
