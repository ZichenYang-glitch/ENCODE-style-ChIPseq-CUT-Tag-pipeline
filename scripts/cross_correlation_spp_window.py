#!/usr/bin/env python3
"""Print the locked run_spp.R exclusion window when the tagAlign is ready."""

from __future__ import annotations

import argparse
import sys

from qc_common import read_cross_correlation_metadata


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", required=True)
    args = parser.parse_args()
    try:
        metadata = read_cross_correlation_metadata(args.metadata)
    except (OSError, ValueError) as exc:
        print(f"ERROR: cross-correlation metadata is invalid: {exc}", file=sys.stderr)
        return 1
    if metadata.get("status") != "ready":
        print(
            f"ERROR: cross-correlation tagAlign is not ready: "
            f"{metadata.get('status')!r}",
            file=sys.stderr,
        )
        return 1
    lower = metadata.get("window_lower")
    upper = metadata.get("window_upper")
    if type(lower) is not int or type(upper) is not int or lower > upper:
        print("ERROR: cross-correlation window is invalid", file=sys.stderr)
        return 1
    print(f"{lower}:{upper}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
