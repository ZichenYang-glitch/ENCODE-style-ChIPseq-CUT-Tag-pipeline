#!/usr/bin/env python3
"""broad_idr_summary.py — Thin wrapper around idr_reproducibility_summary.py.

Provides the Stage 65 broad-peak IDR summary CLI by delegating to the unified
summary script with broad peak-mode default. The assay (chipseq or cuttag) must
still be supplied by the caller.
"""

import sys

from idr_reproducibility_summary import main


if __name__ == "__main__":
    argv = sys.argv[1:]
    if "--peak-mode" not in argv:
        argv = argv + ["--peak-mode", "broad"]
    main(argv)
