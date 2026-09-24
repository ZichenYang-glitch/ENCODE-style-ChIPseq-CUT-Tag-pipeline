#!/usr/bin/env python3
"""Private Hi-TrAC qualification, using the canonical source-provenance bootstrap."""

from pathlib import Path
import runpy
import sys

if not sys.flags.isolated or not sys.flags.no_site:
    raise SystemExit("invoke with Python -I -S")
sys.dont_write_bytecode = True
root = Path(__file__).resolve().parents[1]
provenance = runpy.run_path(str(root / "scripts/source_provenance.py"))
provenance["bootstrap_checkout"](root)
from encode_pipeline.adapters.hitrac_preprocess.entrypoint import main  # noqa: E402

if __name__ == "__main__":
    raise SystemExit(main())
