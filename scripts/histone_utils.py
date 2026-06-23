"""Backward-compatible wrapper for histone target classification.

The implementation has moved to src/encode_pipeline/qc/histone.py.
This file preserves `from histone_utils import classify_histone_target`.
"""

import os
import sys

try:
    from encode_pipeline.qc.histone import classify_histone_target
except ImportError:
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
    from encode_pipeline.qc.histone import classify_histone_target

__all__ = ["classify_histone_target"]
