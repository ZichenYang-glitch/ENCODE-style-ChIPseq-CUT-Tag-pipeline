#!/usr/bin/env python3
"""Unit tests for workflow/rules/idr_paths.smk helpers.

These helpers only perform string interpolation against OUTDIR.  We
monkey-patch OUTDIR here to avoid importing the full Snakefile namespace.
"""

import sys
from pathlib import Path

import pytest


# The helpers under test rely on a module-level OUTDIR global.
IDR_PATHS = Path(__file__).parent.parent / "workflow" / "rules" / "idr_paths.smk"


def _load_idr_paths_module(outdir):
    """Load idr_paths.smk with OUTDIR set to a test value."""
    code = IDR_PATHS.read_text()
    namespace = {"OUTDIR": outdir}
    exec(compile(code, str(IDR_PATHS), "exec"), namespace)
    return namespace


@pytest.fixture
def h():
    return _load_idr_paths_module("results")


@pytest.mark.parametrize(
    "exp,br,expected",
    [
        (
            "EXP1",
            1,
            "results/experiments/EXP1/02_align/biorep1.final.bam",
        ),
        (
            "EXP2",
            2,
            "results/experiments/EXP2/02_align/biorep2.final.bam",
        ),
    ],
)
def test_idr_biorep_bam(h, exp, br, expected):
    assert h["idr_biorep_bam"](exp, br) == expected


@pytest.mark.parametrize(
    "exp,br,expected",
    [
        (
            "EXP1",
            1,
            "results/experiments/EXP1/02_align/biorep1.final.bam.bai",
        ),
    ],
)
def test_idr_biorep_bai(h, exp, br, expected):
    assert h["idr_biorep_bai"](exp, br) == expected


def test_idr_pooled_treatment_bam(h):
    assert (
        h["idr_pooled_treatment_bam"]("EXP1")
        == "results/experiments/EXP1/02_align/EXP1.pooled.final.bam"
    )


def test_idr_pooled_control_bam(h):
    assert (
        h["idr_pooled_control_bam"]("EXP1")
        == "results/experiments/EXP1/02_align/EXP1.pooled.control.final.bam"
    )


@pytest.mark.parametrize(
    "source,pr,expected",
    [
        (
            "biorep1",
            1,
            "results/experiments/EXP1/05_pseudorep/EXP1_biorep1.pr1.bam",
        ),
        (
            "pooled",
            2,
            "results/experiments/EXP1/05_pseudorep/EXP1_pooled.pr2.bam",
        ),
    ],
)
def test_idr_pseudorep_bam(h, source, pr, expected):
    assert h["idr_pseudorep_bam"]("EXP1", source, pr) == expected


def test_idr_pseudorep_bai(h):
    assert (
        h["idr_pseudorep_bai"]("EXP1", "biorep2", 1)
        == "results/experiments/EXP1/05_pseudorep/EXP1_biorep2.pr1.bam.bai"
    )


def test_idr_biorep_bam_old_literal_match():
    """Ensure the helper returns the exact string previously inlined."""
    h = _load_idr_paths_module("results")
    assert h["idr_biorep_bam"]("EXP", 1) == "results/experiments/EXP/02_align/biorep1.final.bam"
