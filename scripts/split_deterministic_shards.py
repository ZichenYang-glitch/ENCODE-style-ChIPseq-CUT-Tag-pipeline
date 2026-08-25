#!/usr/bin/env python3
"""Assign the deterministic test files to mutually exclusive, exhaustive shards.

The partition is a pure function of the checked-out tree: each collected
``test_*.py`` file's relative POSIX path is hashed with SHA-256 and assigned
to ``int(digest) % shard_count``.  No pytest plugins, timing history, or
network services are involved, so every checkout of the same commit yields
the exact same shards on any machine.
"""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
import sys


def shard_index(relative_posix_path: str, shard_count: int) -> int:
    """Return the zero-based shard index for one relative POSIX path."""
    if not relative_posix_path or "\\" in relative_posix_path:
        raise ValueError("shard paths must be nonempty relative POSIX paths")
    digest = hashlib.sha256(relative_posix_path.encode("utf-8")).digest()
    return int.from_bytes(digest, "big") % shard_count


def list_test_files(root: Path) -> list[str]:
    """Return every collected test file below root as sorted POSIX paths."""
    if not root.is_dir() or root.is_symlink():
        raise ValueError(f"shard root is not a directory: {root}")
    files = []
    for path in root.rglob("test_*.py"):
        if not path.is_file() or path.is_symlink():
            continue
        files.append(path.relative_to(root).as_posix())
    unique = sorted(set(files))
    if not unique:
        raise ValueError(f"shard root contains no test_*.py files: {root}")
    return unique


def assign_shards(root: Path, shard_count: int) -> list[list[str]]:
    """Partition the tree's test files into ``shard_count`` ordered shards."""
    if shard_count < 1:
        raise ValueError("shard count must be a positive integer")
    shards: list[list[str]] = [[] for _ in range(shard_count)]
    for relative in list_test_files(root):
        shards[shard_index(relative, shard_count)].append(relative)
    return shards


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("test"))
    parser.add_argument("--shard", type=int, required=True, help="1-based shard")
    parser.add_argument("--of", type=int, required=True, help="total shard count")
    args = parser.parse_args(argv)
    if args.of < 1 or not 1 <= args.shard <= args.of:
        parser.error("--shard must satisfy 1 <= shard <= --of")
    try:
        shards = assign_shards(args.root, args.of)
    except ValueError as error:
        print(f"shard assignment failed: {error}", file=sys.stderr)
        return 2
    for relative in shards[args.shard - 1]:
        print(f"{args.root.as_posix()}/{relative}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
