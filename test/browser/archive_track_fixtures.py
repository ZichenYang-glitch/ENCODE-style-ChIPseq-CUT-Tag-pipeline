"""Generate the two synthetic registered-result sources for the archive browser test."""

from __future__ import annotations

import argparse
from pathlib import Path
import random

import pyBigWig


SIZE_BYTES = 17_995_618


def write_signal(path: Path, value: float) -> None:
    path.parent.mkdir(parents=True, exist_ok=False)
    filler = random.Random(20261125)
    with pyBigWig.open(str(path), "w") as writer:
        writer.addHeader([("chr1", 248956422), ("chr7", 159345973)], maxZooms=0)
        # Keep the file large with real, indexed data away from the test window.
        for chunk in range(100):
            writer.addEntries(
                "chr1",
                chunk * 500_000,
                values=[filler.random() for _ in range(50_000)],
                step=10,
                span=3,
            )
        writer.addEntries(["chr7"], [55_104_400], ends=[55_105_900], values=[value])

    actual_size = path.stat().st_size
    if actual_size != SIZE_BYTES:
        raise RuntimeError(
            f"{path.parent.name} signal.bw has {actual_size} bytes, expected {SIZE_BYTES}. "
            "Use the pyBigWig build in workflow/envs/deeptools.lock."
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dest", required=True, type=Path)
    args = parser.parse_args()
    args.dest.mkdir(parents=True, exist_ok=False)
    for version, value in [(1, 4.5), (2, 9.0)]:
        write_signal(args.dest / f"v{version}" / "signal.bw", value)


if __name__ == "__main__":
    main()
