#!/usr/bin/env python3
"""Deterministically split a BAM into two complementary pseudoreplicates.

Each read/template is assigned to pr1 or pr2 by a stable sha256 hash of
(seed + canonical read name). Paired-end mates always go to the same
pseudoreplicate. Streams records directly from samtools view into two
samtools sort subprocess stdin pipes — no temporary uncompressed files.

Usage:
    python3 scripts/split_pseudoreps.py \
        --input in.bam --out1 in.pr1.bam --out2 in.pr2.bam \
        --seed 42 --threads 4
"""

import argparse
import hashlib
import os
import subprocess
import sys

_SAMTOOLS = os.environ.get("SAMTOOLS", "samtools")


def _canonical_name(query_name):
    """Strip trailing /1 or /2 from a read query name for stable hashing."""
    if query_name.endswith("/1") or query_name.endswith("/2"):
        return query_name[:-2]
    return query_name


def _assign_to_pseudorep(seed, name):
    """Return 1 or 2 for a deterministic pseudoreplicate assignment."""
    h = hashlib.sha256(f"{seed}:{name}".encode()).hexdigest()
    val = int(h, 16)
    return 1 if (val % 2 == 0) else 2


def _count_reads(bam_path):
    """Return read count from a BAM file via samtools view -c."""
    result = subprocess.run(
        [_SAMTOOLS, "view", "-c", bam_path],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        print(f"ERROR: samtools view -c failed for {bam_path}: {result.stderr}",
              file=sys.stderr)
        sys.exit(1)
    try:
        return int(result.stdout.strip())
    except ValueError:
        print(f"ERROR: samtools view -c returned non-integer: {result.stdout}",
              file=sys.stderr)
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description="Deterministic BAM pseudoreplicate splitter"
    )
    parser.add_argument("--input", required=True, help="Input BAM file")
    parser.add_argument("--out1", required=True, help="Output pr1 BAM")
    parser.add_argument("--out2", required=True, help="Output pr2 BAM")
    parser.add_argument("--seed", type=int, required=True, help="Seed for hash")
    parser.add_argument("--threads", type=int, default=4,
                        help="Threads for samtools sort")
    args = parser.parse_args()

    # Check samtools is callable
    try:
        subprocess.run([_SAMTOOLS, "--version"], capture_output=True, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        print(f"ERROR: {_SAMTOOLS} not found on PATH", file=sys.stderr)
        sys.exit(1)

    # Open input SAM stream
    view_proc = subprocess.Popen(
        [_SAMTOOLS, "view", "-h", args.input],
        stdout=subprocess.PIPE,
    )

    # Start two samtools sort subprocesses with stdin pipes
    sort1_proc = subprocess.Popen(
        [_SAMTOOLS, "sort", "-@", str(args.threads), "-o", args.out1, "-"],
        stdin=subprocess.PIPE,
    )
    sort2_proc = subprocess.Popen(
        [_SAMTOOLS, "sort", "-@", str(args.threads), "-o", args.out2, "-"],
        stdin=subprocess.PIPE,
    )

    assert sort1_proc.stdin is not None
    assert sort2_proc.stdin is not None

    try:
        for line in view_proc.stdout:
            line_str = line.decode("utf-8")
            if line_str.startswith("@"):
                sort1_proc.stdin.write(line)
                sort2_proc.stdin.write(line)
            else:
                fields = line_str.split("\t")
                if len(fields) < 1:
                    continue
                qname = fields[0]
                can = _canonical_name(qname)
                assign = _assign_to_pseudorep(args.seed, can)
                if assign == 1:
                    sort1_proc.stdin.write(line)
                else:
                    sort2_proc.stdin.write(line)

        sort1_proc.stdin.close()
        sort2_proc.stdin.close()

        view_proc.wait()
        if view_proc.returncode != 0:
            print("ERROR: samtools view failed", file=sys.stderr)
            sys.exit(1)

        if sort1_proc.wait() != 0 or sort2_proc.wait() != 0:
            print("ERROR: samtools sort failed", file=sys.stderr)
            sys.exit(1)

        subprocess.run(
            [_SAMTOOLS, "index", "-@", str(args.threads), args.out1], check=True)
        subprocess.run(
            [_SAMTOOLS, "index", "-@", str(args.threads), args.out2], check=True)

        n_input = _count_reads(args.input)
        n_pr1 = _count_reads(args.out1)
        n_pr2 = _count_reads(args.out2)

        if n_pr1 + n_pr2 != n_input:
            print(
                f"ERROR: complementarity failed: "
                f"pr1={n_pr1} + pr2={n_pr2} != input={n_input}",
                file=sys.stderr,
            )
            sys.exit(1)

        print(f"Split complete: pr1={n_pr1}, pr2={n_pr2}, "
              f"total={n_pr1 + n_pr2}, input={n_input}")

    finally:
        view_proc.wait()
        for p in (sort1_proc, sort2_proc):
            if p.poll() is None:
                p.stdin.close()
                p.terminate()
                p.wait()


if __name__ == "__main__":
    main()
