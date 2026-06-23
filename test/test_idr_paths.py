#!/usr/bin/env python3
"""Unit tests for workflow/rules/idr_paths.smk helpers.

These helpers perform string interpolation against OUTDIR or dispatch on
sample metadata. We monkey-patch the required Snakefile globals here to
avoid importing the full Snakefile namespace.
"""

import sys
from pathlib import Path

import pytest


# The helpers under test rely on a module-level OUTDIR global.
IDR_PATHS = Path(__file__).parent.parent / "workflow" / "rules" / "idr_paths.smk"


def _load_idr_paths_module(
    outdir,
    pooled_control_experiments=None,
    sample_map=None,
    treatment_samples_by_experiment=None,
):
    """Load idr_paths.smk with Snakefile globals monkey-patched."""
    code = IDR_PATHS.read_text()
    namespace = {
        "OUTDIR": outdir,
        "POOLED_CONTROL_EXPERIMENTS": set(pooled_control_experiments or []),
        "SAMPLE_MAP": sample_map or {},
        "TREATMENT_SAMPLES_BY_EXPERIMENT": treatment_samples_by_experiment or {},
        "_normalize_genome": lambda genome: {"hg38": "hs", "mm10": "mm"}.get(genome, genome),
        "_tool_param": lambda tool, key, default: {
            ("idr_macs3", "pvalue"): 0.1,
            ("idr_macs3", "extra_args"): "",
            ("macs3", "broad_cutoff"): 0.1,
        }.get((tool, key), default),
    }
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
        # Assay-prefixed sources used by ATAC, CUT&Tag, and broad variants.
        (
            "atac_biorep1",
            1,
            "results/experiments/EXP1/05_pseudorep/EXP1_atac_biorep1.pr1.bam",
        ),
        (
            "cuttag_pooled",
            2,
            "results/experiments/EXP1/05_pseudorep/EXP1_cuttag_pooled.pr2.bam",
        ),
        (
            "broad_chipseq_biorep2",
            1,
            "results/experiments/EXP1/05_pseudorep/EXP1_broad_chipseq_biorep2.pr1.bam",
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


def test_idr_inputs_match_legacy_inlined_paths(h):
    """Contract test: helper-built inputs equal the old inlined literals.

    This guards against accidental signature drift in the helpers used by
    idr.smk, idr_atac.smk, idr_cuttag.smk, and idr_broad.smk.
    """
    exp = "EXP"
    br = 1
    assert h["idr_biorep_bam"](exp, br) == f"results/experiments/{exp}/02_align/biorep{br}.final.bam"
    assert h["idr_biorep_bai"](exp, br) == f"results/experiments/{exp}/02_align/biorep{br}.final.bam.bai"
    assert (
        h["idr_pooled_treatment_bam"](exp)
        == f"results/experiments/{exp}/02_align/{exp}.pooled.final.bam"
    )
    assert (
        h["idr_pooled_control_bam"](exp)
        == f"results/experiments/{exp}/02_align/{exp}.pooled.control.final.bam"
    )
    assert (
        h["idr_pseudorep_bam"](exp, "biorep1", 1)
        == f"results/experiments/{exp}/05_pseudorep/{exp}_biorep1.pr1.bam"
    )
    assert (
        h["idr_pseudorep_bam"](exp, "atac_pooled", 2)
        == f"results/experiments/{exp}/05_pseudorep/{exp}_atac_pooled.pr2.bam"
    )


def test_idr_biorep_peaks_inputs(h):
    assert h["idr_biorep_peaks_inputs"]("EXP1", 1) == [
        "results/experiments/EXP1/02_align/biorep1.final.bam",
        "results/experiments/EXP1/02_align/biorep1.final.bam.bai",
    ]


def test_idr_biorep_peaks_inputs_with_control():
    h = _load_idr_paths_module("results", pooled_control_experiments=["EXP1"])
    assert h["idr_biorep_peaks_inputs"]("EXP1", 2) == [
        "results/experiments/EXP1/02_align/biorep2.final.bam",
        "results/experiments/EXP1/02_align/biorep2.final.bam.bai",
        "results/experiments/EXP1/02_align/EXP1.pooled.control.final.bam",
    ]


@pytest.mark.parametrize(
    "source,expected",
    [
        ("pooled", "results/experiments/EXP1/02_align/EXP1.pooled.final.bam"),
        ("biorep2", "results/experiments/EXP1/02_align/biorep2.final.bam"),
    ],
)
def test_idr_split_input(h, source, expected):
    assert h["idr_split_input"]("EXP1", source) == expected


def test_idr_pseudorep_peaks_inputs(h):
    assert h["idr_pseudorep_peaks_inputs"]("EXP1", "biorep2", "1") == [
        "results/experiments/EXP1/05_pseudorep/EXP1_biorep2.pr1.bam",
        "results/experiments/EXP1/05_pseudorep/EXP1_biorep2.pr1.bam.bai",
    ]


def test_idr_pseudorep_peaks_inputs_with_prefix():
    h = _load_idr_paths_module("results", pooled_control_experiments=["EXP1"])
    assert h["idr_pseudorep_peaks_inputs"]("EXP1", "pooled", "2", source_prefix="atac_") == [
        "results/experiments/EXP1/05_pseudorep/EXP1_atac_pooled.pr2.bam",
        "results/experiments/EXP1/05_pseudorep/EXP1_atac_pooled.pr2.bam.bai",
        "results/experiments/EXP1/02_align/EXP1.pooled.control.final.bam",
    ]


@pytest.fixture
def macs3_h():
    """Fixture for idr_macs3_args with a single mocked PE treatment sample."""
    return _load_idr_paths_module(
        "results",
        sample_map={
            "S1": {
                "layout": "PE",
                "genome": "hg38",
            },
        },
        treatment_samples_by_experiment={"EXP1": ["S1"]},
    )


@pytest.mark.parametrize(
    "assay,peak_mode,expected_subset",
    [
        ("chipseq", "narrow", ["-f BAMPE", "-g hs", "-p 0.1"]),
        ("atac", "narrow", ["--nomodel", "--shift -100", "--extsize 200"]),
        ("cuttag", "narrow", ["--nomodel", "--shift -100", "--extsize 200"]),
        ("chipseq", "broad", ["--broad", "--broad-cutoff 0.1"]),
        ("cuttag", "broad", ["--broad", "--broad-cutoff 0.1"]),
    ],
)
def test_idr_macs3_args_dispatch(macs3_h, assay, peak_mode, expected_subset):
    args = macs3_h["idr_macs3_args"]("EXP1", assay, peak_mode)
    for fragment in expected_subset:
        assert fragment in args, f"{fragment!r} not in {args!r}"


def test_idr_macs3_args_se_layout():
    h = _load_idr_paths_module(
        "results",
        sample_map={"S1": {"layout": "SE", "genome": "mm10"}},
        treatment_samples_by_experiment={"EXP2": ["S1"]},
    )
    assert "-f BAM" in h["idr_macs3_args"]("EXP2", "chipseq", "narrow")
    assert "-g mm" in h["idr_macs3_args"]("EXP2", "chipseq", "narrow")


def test_idr_macs3_args_empty_experiment():
    h = _load_idr_paths_module("results")
    assert h["idr_macs3_args"]("NOEXP", "chipseq", "narrow") == ""
