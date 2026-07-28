#!/usr/bin/env python3
"""Checkout-verified wrapper for the encode_pipeline validator.

The implementation has moved to src/encode_pipeline/. This file preserves
`python3 scripts/validate_samples.py --config config/config.yaml`.
"""

from pathlib import Path
import runpy
import sys

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT / "src"))
_require_checkout = runpy.run_path(
    str(Path(__file__).with_name("source_provenance.py"))
)["require_checkout"]
_require_checkout(REPOSITORY_ROOT)

from encode_pipeline.cli.validate import main  # noqa: E402

if __name__ == "__main__":
    main()
