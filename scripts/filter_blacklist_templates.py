#!/usr/bin/env python3
"""Canonical analysis BAM template closure, organelle, and blacklist filter.

Applies:
1. Complete template closure:
   - For PE: proper-pair (0x2), same contig, paired R1/R2, mutually consistent
     mate coordinates (pos1==pnext2, pos2==pnext1), mutually consistent non-zero
     TLEN (tlen1==-tlen2), legal negative TLEN kept as abs(TLEN). Incomplete
     templates or orphan mates fail closed and are excluded as a pair.
   - For SE: single mapped record with positive reference span.
2. Duplicate policy:
   - When --exclude-dup is set, any template with 0x400 duplicate flags is excluded.
3. Organelle policy:
   - When --organelle-mode is 'exclude', verifies organelle contig candidates
     against BAM @SQ and excludes whole templates with any mate on an organelle contig.
   - If explicit candidate contigs are provided, all must exist in @SQ.
   - If candidate contigs have no matches in @SQ and mode is 'exclude', fails closed.
4. Blacklist policy:
   - When --blacklist is provided, verifies BED coordinates (0 <= start < end <= contig_length)
     against BAM @SQ and excludes whole templates overlapping blacklist intervals.
   - When --blacklist is omitted, blacklist overlap is skipped while template closure,
     duplicate, and organelle filtering remain active.
"""

import argparse
import json
import os
import subprocess
import sys

from blacklist_filter import Blacklist


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", required=True, help="Input BAM path")
    parser.add_argument("--output", required=True, help="Output filtered BAM path")
    parser.add_argument(
        "--dropped-count", required=True, help="Output blacklist dropped count TSV"
    )
    parser.add_argument("--sample", required=True, help="Sample identifier")
    parser.add_argument("--layout", required=True, choices=("PE", "SE"))
    parser.add_argument(
        "--exclude-dup",
        action="store_true",
        default=False,
        help="Exclude duplicate (0x400) templates",
    )
    parser.add_argument(
        "--organelle-mode",
        choices=("auto", "exclude", "report", "keep"),
        default="keep",
        help="Organelle filtering mode",
    )
    parser.add_argument(
        "--organelle-chroms",
        nargs="*",
        default=[],
        help="Candidate organelle contig names",
    )
    parser.add_argument(
        "--organelle-chroms-explicit",
        action="store_true",
        default=False,
        help="True if organelle-chroms was explicitly configured",
    )
    parser.add_argument(
        "--blacklist",
        default=None,
        help="Optional path to blacklist BED file",
    )
    parser.add_argument(
        "--chrom-sizes",
        default=None,
        help="Optional path to chrom_sizes file for contig length verification",
    )
    parser.add_argument(
        "--fai",
        default=None,
        help="Optional path to FASTA index (.fai) for contig length verification",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Compression/decompression threads for samtools",
    )
    parser.add_argument(
        "--stage-metrics",
        default=None,
        help="Optional path to output stage metrics JSON",
    )
    return parser.parse_args(argv)


def cigar_ref_length(cigar: str) -> int:
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


def read_bam_header(bam: str) -> dict[str, int]:
    """Return dictionary of contig name to length from BAM @SQ headers."""
    process = subprocess.Popen(
        ["samtools", "view", "-H", bam],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    stdout, stderr = process.communicate()
    if process.returncode != 0:
        raise RuntimeError(f"samtools view -H failed on {bam}: {stderr.strip()}")
    sq: dict[str, int] = {}
    for line in stdout.splitlines():
        if not line.startswith("@SQ\t"):
            continue
        fields = dict(
            field.split(":", 1) for field in line.split("\t")[1:] if ":" in field
        )
        name = fields.get("SN")
        length = fields.get("LN")
        if name and length and length.isdigit():
            sq[name] = int(length)
    return sq


def validate_reference_contig_lengths(
    sq: dict[str, int],
    chrom_sizes_path: str | None = None,
    fai_path: str | None = None,
    sample: str = "",
) -> set[str]:
    """Verify configured reference contigs against BAM @SQ.

    Every provided reference source must share at least one contig with the
    BAM; otherwise the source belongs to a different reference role/assembly.
    Shared contig lengths are checked. BAM-only alt/decoy contigs are returned
    as unchecked rather than rejected, because the reference resources may
    intentionally describe only primary contigs.
    """
    checked_contigs: set[str] = set()
    for source_path, source_name in [
        (fai_path, "FAI"),
        (chrom_sizes_path, "chrom_sizes"),
    ]:
        if not (source_path and os.path.isfile(source_path)):
            continue
        ref_lengths: dict[str, int] = {}
        with open(source_path, encoding="utf-8") as handle:
            for line in handle:
                parts = line.strip().split()
                if len(parts) >= 2 and parts[1].isdigit():
                    ref_lengths[parts[0]] = int(parts[1])
        if not ref_lengths:
            raise RuntimeError(
                f"Sample {sample}: reference {source_name} ({source_path}) "
                "contains no valid contig lengths"
            )
        common = set(sq).intersection(ref_lengths)
        if not common:
            raise RuntimeError(
                f"Sample {sample}: BAM @SQ and reference {source_name} "
                f"({source_path}) have no contigs in common"
            )
        for chrom in sorted(common):
            if sq[chrom] != ref_lengths[chrom]:
                raise RuntimeError(
                    f"Sample {sample}: BAM @SQ contig length mismatch for {chrom!r}: "
                    f"BAM has {sq[chrom]} bp, reference {source_name} ({source_path}) has {ref_lengths[chrom]} bp"
                )
        checked_contigs.update(common)
    return set(sq) - checked_contigs


def validate_organelle_contigs(
    sq: dict[str, int],
    candidate_chroms: list[str],
    is_explicit: bool,
    mode: str,
    sample: str,
) -> set[str]:
    """Resolve matched organelle contigs and enforce fail-closed verification."""
    if is_explicit:
        for c in candidate_chroms:
            if c not in sq:
                raise RuntimeError(
                    f"Sample {sample}: explicit organelle contig {c!r} not found "
                    f"in BAM @SQ reference contigs"
                )
        matched = set(candidate_chroms)
    else:
        matched = set(candidate_chroms).intersection(sq.keys())

    if mode == "keep":
        return set()

    if mode == "exclude" and not matched:
        raise RuntimeError(
            f"Sample {sample}: organelle_mode is 'exclude' but none of candidate "
            f"contigs {candidate_chroms} exist in BAM @SQ contigs. Explicitly configure "
            f"genome_resources.<genome>.organelle_chroms or set organelle_mode: keep."
        )
    return matched


def validate_and_load_blacklist(
    blacklist_path: str, sq: dict[str, int], sample: str
) -> Blacklist:
    """Validate blacklist BED against BAM contigs and return Blacklist instance."""
    if not os.path.isfile(blacklist_path):
        raise RuntimeError(
            f"Sample {sample}: blacklist file not found: {blacklist_path}"
        )

    valid_intervals = 0
    with open(blacklist_path, encoding="utf-8") as handle:
        for line_no, raw_line in enumerate(handle, 1):
            line = raw_line.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                raise RuntimeError(
                    f"Sample {sample}: blacklist line {line_no} has fewer than 3 fields: {raw_line.rstrip()!r}"
                )
            chrom = fields[0]
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError as exc:
                raise RuntimeError(
                    f"Sample {sample}: blacklist line {line_no} non-integer coordinates: {raw_line.rstrip()!r}"
                ) from exc
            if start < 0 or end <= start:
                raise RuntimeError(
                    f"Sample {sample}: blacklist line {line_no} invalid interval start={start} >= end={end}"
                )
            if chrom not in sq:
                raise RuntimeError(
                    f"Sample {sample}: blacklist line {line_no} contig {chrom!r} not found in BAM @SQ contigs"
                )
            if end > sq[chrom]:
                raise RuntimeError(
                    f"Sample {sample}: blacklist line {line_no} interval {chrom}:{start}-{end} exceeds contig length {sq[chrom]}"
                )
            valid_intervals += 1

    bl = Blacklist(blacklist_path)
    if valid_intervals == 0 or bl.total_intervals == 0:
        raise RuntimeError(
            f"Sample {sample}: blacklist BED {blacklist_path} is empty or contains no valid intervals"
        )
    return bl


def scan_drop_qnames(
    bam: str,
    layout: str,
    exclude_dup: bool,
    matched_organelle: set[str],
    blacklist: Blacklist | None,
    threads: int = 1,
) -> tuple[set[str], dict[str, str], int]:
    """Scan BAM and identify QNAMEs to drop, tracking reasons and total alignments."""
    view_cmd = ["samtools", "view"]
    if threads > 1:
        view_cmd.extend(["-@", str(threads)])
    view_cmd.append(bam)

    process = subprocess.Popen(
        view_cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    assert process.stdout is not None

    drop_qnames: set[str] = set()
    drop_reasons: dict[str, str] = {}
    total_alignments = 0

    if layout == "SE":
        for line in process.stdout:
            if not line or line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                continue
            total_alignments += 1
            qname = fields[0]
            flag = int(fields[1])
            chrom = fields[2]
            pos = int(fields[3])
            cigar = fields[5]

            # Invalid template: paired, unmapped, non-positive pos, empty CIGAR, secondary/supplementary
            ref_len = cigar_ref_length(cigar)
            if (flag & 0x901) or (flag & 0x4) or pos <= 0 or ref_len <= 0:
                drop_qnames.add(qname)
                drop_reasons[qname] = "invalid"
                continue

            # Duplicate exclusion
            if exclude_dup and (flag & 0x400):
                drop_qnames.add(qname)
                drop_reasons[qname] = "duplicate"
                continue

            # Organelle exclusion
            if chrom in matched_organelle:
                drop_qnames.add(qname)
                drop_reasons[qname] = "organelle"
                continue

            # Blacklist overlap
            if blacklist is not None:
                start0 = pos - 1
                end0 = start0 + ref_len
                if blacklist.overlaps(chrom, start0, end0):
                    drop_qnames.add(qname)
                    drop_reasons[qname] = "blacklist"
                    continue

    else:  # PE
        open_pairs: dict[str, tuple] = {}
        seen_counts: dict[str, int] = {}

        for line in process.stdout:
            if not line or line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                continue
            total_alignments += 1
            qname = fields[0]
            flag = int(fields[1])
            chrom = fields[2]
            pos = int(fields[3])
            cigar = fields[5]
            rnext = fields[6]
            pnext = int(fields[7])
            tlen = int(fields[8])

            count = seen_counts.get(qname, 0) + 1
            seen_counts[qname] = count

            if count == 1:
                open_pairs[qname] = (flag, chrom, pos, cigar, rnext, pnext, tlen)
            elif count == 2:
                mate1 = open_pairs.pop(qname)
                flag1, chrom1, pos1, cigar1, rnext1, pnext1, tlen1 = mate1
                flag2, chrom2, pos2, cigar2, rnext2, pnext2, tlen2 = (
                    flag,
                    chrom,
                    pos,
                    cigar,
                    rnext,
                    pnext,
                    tlen,
                )

                # Proper template closure verification:
                # 1. Not secondary (0x100) or supplementary (0x800)
                # 2. Both paired (0x1)
                # 3. Both mapped (not 0x4) and mate mapped (not 0x8)
                # 4. Both proper-pair (0x2)
                # 5. Exactly one R1 and one R2 (0x40 vs 0x80)
                # 6. Same contig
                # 7. Coordinate reciprocity: pos1 == pnext2 and pos2 == pnext1
                # 8. TLEN reciprocity: tlen1 == -tlen2 and tlen1 != 0
                is_valid = True
                if (flag1 & 0x900) or (flag2 & 0x900):
                    is_valid = False
                elif not (flag1 & 0x1) or not (flag2 & 0x1):
                    is_valid = False
                elif (flag1 & 0x4) or (flag2 & 0x4) or (flag1 & 0x8) or (flag2 & 0x8):
                    is_valid = False
                elif not (flag1 & 0x2) or not (flag2 & 0x2):
                    is_valid = False
                elif not (
                    (flag1 & 0x40 and flag2 & 0x80) or (flag1 & 0x80 and flag2 & 0x40)
                ):
                    is_valid = False
                elif chrom1 != chrom2:
                    is_valid = False
                elif rnext1 not in ("=", chrom1) or rnext2 not in ("=", chrom2):
                    is_valid = False
                elif pos1 <= 0 or pos2 <= 0:
                    is_valid = False
                elif cigar_ref_length(cigar1) <= 0 or cigar_ref_length(cigar2) <= 0:
                    is_valid = False
                elif pos1 != pnext2 or pos2 != pnext1:
                    is_valid = False
                elif tlen1 == 0 or tlen2 == 0 or tlen1 != -tlen2:
                    is_valid = False

                if not is_valid:
                    drop_qnames.add(qname)
                    drop_reasons[qname] = "invalid"
                    continue

                # Duplicate exclusion
                if exclude_dup and ((flag1 & 0x400) or (flag2 & 0x400)):
                    drop_qnames.add(qname)
                    drop_reasons[qname] = "duplicate"
                    continue

                # Organelle exclusion
                if chrom1 in matched_organelle or chrom2 in matched_organelle:
                    drop_qnames.add(qname)
                    drop_reasons[qname] = "organelle"
                    continue

                # Blacklist overlap
                if blacklist is not None:
                    # Template fragment span: [min(pos1, pos2) - 1, min(pos1, pos2) - 1 + abs(tlen)]
                    start0 = min(pos1, pos2) - 1
                    end0 = start0 + abs(tlen1)
                    len1 = cigar_ref_length(cigar1)
                    len2 = cigar_ref_length(cigar2)
                    m1_start = pos1 - 1
                    m1_end = m1_start + max(len1, 1)
                    m2_start = pos2 - 1
                    m2_end = m2_start + max(len2, 1)

                    if (
                        blacklist.overlaps(chrom1, start0, end0)
                        or blacklist.overlaps(chrom1, m1_start, m1_end)
                        or blacklist.overlaps(chrom2, m2_start, m2_end)
                    ):
                        drop_qnames.add(qname)
                        drop_reasons[qname] = "blacklist"
                        continue

            else:  # count > 2: >2 records for same QNAME
                drop_qnames.add(qname)
                drop_reasons[qname] = "invalid"

        # Any unpartnered mate left in open_pairs at EOF is an incomplete template / orphan
        for q in open_pairs:
            drop_qnames.add(q)
            drop_reasons[q] = "invalid"

    stderr = process.stderr.read() if process.stderr is not None else ""
    ret = process.wait()
    if ret != 0:
        raise RuntimeError(f"samtools view failed on {bam}: {stderr.strip()}")

    return drop_qnames, drop_reasons, total_alignments


def write_filtered_bam(
    bam: str,
    output: str,
    drop_qnames: set[str],
    drop_reasons: dict[str, str] | None = None,
    threads: int = 1,
) -> tuple[int, int, dict[str, int]]:
    """Stream BAM records, writing only retained records to output BAM."""
    os.makedirs(os.path.dirname(os.path.abspath(output)), exist_ok=True)
    reader_cmd = ["samtools", "view", "-h"]
    writer_cmd = ["samtools", "view", "-b", "-o", output, "-"]
    if threads > 1:
        reader_cmd.extend(["-@", str(threads)])
        writer_cmd.extend(["-@", str(threads)])
    reader_cmd.append(bam)

    reader = subprocess.Popen(
        reader_cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    writer = subprocess.Popen(
        writer_cmd,
        stdin=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    assert reader.stdout is not None
    assert writer.stdin is not None

    dropped_alignments = 0
    kept_alignments = 0
    dropped_by_reason: dict[str, int] = {}

    try:
        for line in reader.stdout:
            if line.startswith("@"):
                writer.stdin.write(line)
                continue
            qname = line.split("\t", 1)[0]
            if qname in drop_qnames:
                dropped_alignments += 1
                if drop_reasons is not None:
                    reason = drop_reasons.get(qname, "unknown")
                    dropped_by_reason[reason] = dropped_by_reason.get(reason, 0) + 1
            else:
                kept_alignments += 1
                writer.stdin.write(line)
    finally:
        writer.stdin.close()

    reader_err = reader.stderr.read() if reader.stderr is not None else ""
    writer_err = writer.stderr.read() if writer.stderr is not None else ""
    r_code = reader.wait()
    w_code = writer.wait()
    if r_code != 0:
        raise RuntimeError(f"samtools view -h failed: {reader_err.strip()}")
    if w_code != 0:
        raise RuntimeError(f"samtools view -b failed: {writer_err.strip()}")

    return dropped_alignments, kept_alignments, dropped_by_reason


def main(argv=None):
    args = parse_args(argv)
    sq = read_bam_header(args.bam)
    validate_reference_contig_lengths(
        sq=sq,
        chrom_sizes_path=args.chrom_sizes,
        fai_path=args.fai,
        sample=args.sample,
    )

    # Flatten candidate organelle contigs
    candidate_chroms: list[str] = []
    for item in args.organelle_chroms:
        candidate_chroms.extend(item.split())

    matched_organelle = validate_organelle_contigs(
        sq=sq,
        candidate_chroms=candidate_chroms,
        is_explicit=args.organelle_chroms_explicit,
        mode=args.organelle_mode,
        sample=args.sample,
    )

    # Only exclude organelle contigs when organelle_mode is 'exclude'
    organelle_to_exclude = (
        matched_organelle if args.organelle_mode == "exclude" else set()
    )

    blacklist = None
    if args.blacklist and args.blacklist != "None":
        blacklist = validate_and_load_blacklist(args.blacklist, sq, args.sample)

    drop_qnames, drop_reasons, total_alignments = scan_drop_qnames(
        bam=args.bam,
        layout=args.layout,
        exclude_dup=args.exclude_dup,
        matched_organelle=organelle_to_exclude,
        blacklist=blacklist,
        threads=args.threads,
    )

    dropped_alignments, kept_alignments, dropped_by_reason = write_filtered_bam(
        bam=args.bam,
        output=args.output,
        drop_qnames=drop_qnames,
        drop_reasons=drop_reasons,
        threads=args.threads,
    )

    # Build stage counts for blacklist_dropped.tsv (measuring specifically the blacklist stage)
    dropped_bl_templates = sum(1 for r in drop_reasons.values() if r == "blacklist")
    dropped_bl_alignments = dropped_by_reason.get("blacklist", 0)
    entering_alignments = kept_alignments + dropped_bl_alignments
    if blacklist is not None:
        dropped_tpl_str = str(dropped_bl_templates)
        dropped_aln_str = str(dropped_bl_alignments)
        kept_str = str(kept_alignments)
        total_str = str(entering_alignments)
    else:
        dropped_tpl_str = "0"
        dropped_aln_str = "0"
        kept_str = str(kept_alignments)
        total_str = str(kept_alignments)

    os.makedirs(os.path.dirname(os.path.abspath(args.dropped_count)), exist_ok=True)
    with open(args.dropped_count, "w", encoding="utf-8") as handle:
        handle.write(
            "sample\tlayout\tdropped_templates\tdropped_alignments\t"
            "kept_alignments\ttotal_alignments\n"
        )
        handle.write(
            f"{args.sample}\t{args.layout}\t{dropped_tpl_str}\t"
            f"{dropped_aln_str}\t{kept_str}\t{total_str}\n"
        )

    stage_metrics_path = (
        args.stage_metrics
        if args.stage_metrics
        else os.path.splitext(args.dropped_count)[0] + ".stage_metrics.json"
    )
    os.makedirs(os.path.dirname(os.path.abspath(stage_metrics_path)), exist_ok=True)

    dropped_templates_by_reason = {
        "invalid": 0,
        "duplicate": 0,
        "organelle": 0,
        "blacklist": 0,
    }
    for reason in drop_reasons.values():
        dropped_templates_by_reason[reason] = (
            dropped_templates_by_reason.get(reason, 0) + 1
        )

    dropped_alignments_by_reason = {
        "invalid": 0,
        "duplicate": 0,
        "organelle": 0,
        "blacklist": 0,
    }
    for reason, count in dropped_by_reason.items():
        dropped_alignments_by_reason[reason] = count

    stage_metrics = {
        "sample": args.sample,
        "layout": args.layout,
        "organelle_mode": args.organelle_mode,
        "exclude_dup": args.exclude_dup,
        "has_blacklist": blacklist is not None,
        "total_alignments": total_alignments,
        "kept_alignments": kept_alignments,
        "dropped_alignments": dropped_alignments,
        "dropped_alignments_by_reason": dropped_alignments_by_reason,
        "dropped_by_reason": dropped_alignments_by_reason,
        "dropped_templates": len(drop_qnames),
        "dropped_templates_by_reason": dropped_templates_by_reason,
    }
    with open(stage_metrics_path, "w", encoding="utf-8") as handle:
        json.dump(stage_metrics, handle, indent=2)

    print(
        f"{args.sample}: filtered BAM wrote {kept_alignments} alignments "
        f"(dropped {dropped_alignments} across all filters; total {total_alignments})"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
