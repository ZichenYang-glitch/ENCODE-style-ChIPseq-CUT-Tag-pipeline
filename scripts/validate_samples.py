#!/usr/bin/env python3
"""Backward-compatible wrapper for the encode_pipeline validator.

The implementation has moved to src/encode_pipeline/. This file preserves
`python3 scripts/validate_samples.py --config config/config.yaml`.
"""

import os
import sys

# Local-source fallback for environments where encode-pipeline is not installed.
try:
    from encode_pipeline.cli.validate import main
except ImportError:
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
    from encode_pipeline.cli.validate import main

if __name__ == "__main__":
    main()
