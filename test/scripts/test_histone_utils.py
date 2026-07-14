"""Behavior tests for histone target classification."""

import pytest

from scripts.histone_utils import classify_histone_target


@pytest.mark.parametrize(
    ("target", "peak_mode", "histone_class", "status"),
    [
        ("H3K27me3", "broad", "broad_like", "ok"),
        ("H3K27me3", "narrow", "broad_like", "mismatch"),
        ("H3K27ac", "broad", "context_dependent", "ok"),
        ("H3K27ac", "narrow", "context_dependent", "ok"),
        ("H3K4me3", "narrow", "narrow_like", "ok"),
        ("CTCF", "narrow", "unknown", "unknown"),
        ("h3k27me3", "broad", "broad_like", "ok"),
        ("H2A.Z", "narrow", "narrow_like", "ok"),
    ],
)
def test_classify_histone_target(target, peak_mode, histone_class, status):
    result = classify_histone_target(target, peak_mode)

    assert result["inferred_histone_class"] == histone_class
    assert result["peak_mode_status"] == status
