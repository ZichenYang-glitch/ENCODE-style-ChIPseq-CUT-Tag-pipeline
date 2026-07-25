"""Compatibility import for the package-owned deterministic demo fixture."""

from __future__ import annotations

import importlib
import sys
from pathlib import Path


SOURCE_ROOT = Path(__file__).resolve().parents[1] / "src"
sys.path.insert(0, str(SOURCE_ROOT))
_fixture = importlib.import_module("encode_pipeline.cli.results_visibility_fixture")

module_path = Path(_fixture.__file__ or "").resolve()
if not module_path.is_relative_to(SOURCE_ROOT):
    raise RuntimeError("demo fixture did not resolve from this source checkout")
FIXTURE_SENTINEL_CONTENT = _fixture.FIXTURE_SENTINEL_CONTENT
FIXTURE_SENTINEL_NAME = _fixture.FIXTURE_SENTINEL_NAME
ResultsVisibilityInputs = _fixture.ResultsVisibilityInputs
prepare_results_visibility_fixture = _fixture.prepare_results_visibility_fixture


__all__ = [
    "FIXTURE_SENTINEL_CONTENT",
    "FIXTURE_SENTINEL_NAME",
    "ResultsVisibilityInputs",
    "prepare_results_visibility_fixture",
]
