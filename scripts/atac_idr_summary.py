#!/usr/bin/env python3
"""atac_idr_summary.py — Thin wrapper around idr_reproducibility_summary.py.

Provides the Stage 55 ATAC IDR summary CLI by delegating to the unified
summary script with ATAC narrow defaults.
"""

import sys

from idr_reproducibility_summary import main


if __name__ == "__main__":
    argv = sys.argv[1:]
    if "--assay" not in argv:
        argv = argv + ["--assay", "atac"]
    if "--peak-mode" not in argv:
        argv = argv + ["--peak-mode", "narrow"]
    main(argv)
