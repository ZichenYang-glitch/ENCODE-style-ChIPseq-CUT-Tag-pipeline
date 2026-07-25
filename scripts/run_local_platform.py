#!/usr/bin/env python3
"""Compatibility wrapper for the canonical HelixWeave local platform CLI."""

from __future__ import annotations

import importlib
import sys
from pathlib import Path


SOURCE_ROOT = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(SOURCE_ROOT))
_local_platform = importlib.import_module("encode_pipeline.cli.local_platform")

module_path = Path(_local_platform.__file__ or "").resolve()
if not module_path.is_relative_to(SOURCE_ROOT):
    raise RuntimeError("local platform CLI did not resolve from this source checkout")
main = _local_platform.main


if __name__ == "__main__":
    raise SystemExit(main())
