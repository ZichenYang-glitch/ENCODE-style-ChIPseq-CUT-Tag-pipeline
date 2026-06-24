#!/usr/bin/env python3
"""cuttag_idr_summary.py — Thin wrapper around idr_reproducibility_summary.py.

Provides the Stage 64 CUT&Tag narrow IDR summary CLI by delegating to the
unified summary script with CUT&Tag narrow defaults.
"""

import sys

from idr_reproducibility_summary import main


if __name__ == "__main__":
    argv = sys.argv[1:]
    if "--assay" not in argv:
        argv = argv + ["--assay", "cuttag"]
    if "--peak-mode" not in argv:
        argv = argv + ["--peak-mode", "narrow"]
    main(argv)
