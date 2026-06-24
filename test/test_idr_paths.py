#!/usr/bin/env python3
"""Unit tests for workflow/rules/idr_paths.smk helpers.

These helpers perform string interpolation against OUTDIR or dispatch on
sample metadata. We monkey-patch the required Snakefile globals here to
avoid importing the full Snakefile namespace.
"""

import pytest


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
def test_idr_biorep_bam(idr_paths_namespace, exp, br, expected):
    assert idr_paths_namespace["idr_biorep_bam"](exp, br) == expected


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
def test_idr_biorep_bai(idr_paths_namespace, exp, br, expected):
    assert idr_paths_namespace["idr_biorep_bai"](exp, br) == expected


def test_idr_pooled_treatment_bam(idr_paths_namespace):
    assert (
        idr_paths_namespace["idr_pooled_treatment_bam"]("EXP1")
        == "results/experiments/EXP1/02_align/EXP1.pooled.final.bam"
    )


def test_idr_pooled_control_bam(idr_paths_namespace):
    assert (
        idr_paths_namespace["idr_pooled_control_bam"]("EXP1")
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
def test_idr_pseudorep_bam(idr_paths_namespace, source, pr, expected):
    assert idr_paths_namespace["idr_pseudorep_bam"]("EXP1", source, pr) == expected


def test_idr_pseudorep_bai(idr_paths_namespace):
    assert (
        idr_paths_namespace["idr_pseudorep_bai"]("EXP1", "biorep2", 1)
        == "results/experiments/EXP1/05_pseudorep/EXP1_biorep2.pr1.bam.bai"
    )


def test_idr_biorep_bam_old_literal_match(idr_paths_namespace):
    """Ensure the helper returns the exact string previously inlined."""
    assert idr_paths_namespace["idr_biorep_bam"]("EXP", 1) == "results/experiments/EXP/02_align/biorep1.final.bam"


def test_idr_inputs_match_legacy_inlined_paths(idr_paths_namespace):
    """Contract test: helper-built inputs equal the old inlined literals.

    This guards against accidental signature drift in the helpers used by
    idr.smk, idr_atac.smk, idr_cuttag.smk, and idr_broad.smk.
    """
    exp = "EXP"
    br = 1
    assert idr_paths_namespace["idr_biorep_bam"](exp, br) == f"results/experiments/{exp}/02_align/biorep{br}.final.bam"
    assert idr_paths_namespace["idr_biorep_bai"](exp, br) == f"results/experiments/{exp}/02_align/biorep{br}.final.bam.bai"
    assert (
        idr_paths_namespace["idr_pooled_treatment_bam"](exp)
        == f"results/experiments/{exp}/02_align/{exp}.pooled.final.bam"
    )
    assert (
        idr_paths_namespace["idr_pooled_control_bam"](exp)
        == f"results/experiments/{exp}/02_align/{exp}.pooled.control.final.bam"
    )
    assert (
        idr_paths_namespace["idr_pseudorep_bam"](exp, "biorep1", 1)
        == f"results/experiments/{exp}/05_pseudorep/{exp}_biorep1.pr1.bam"
    )
    assert (
        idr_paths_namespace["idr_pseudorep_bam"](exp, "atac_pooled", 2)
        == f"results/experiments/{exp}/05_pseudorep/{exp}_atac_pooled.pr2.bam"
    )


def test_idr_biorep_peaks_inputs(idr_paths_namespace):
    assert idr_paths_namespace["idr_biorep_peaks_inputs"]("EXP1", 1) == [
        "results/experiments/EXP1/02_align/biorep1.final.bam",
        "results/experiments/EXP1/02_align/biorep1.final.bam.bai",
    ]


def test_idr_biorep_peaks_inputs_with_control(idr_paths_namespace_with_control):
    assert idr_paths_namespace_with_control["idr_biorep_peaks_inputs"]("EXP1", 2) == [
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
def test_idr_split_input(idr_paths_namespace, source, expected):
    assert idr_paths_namespace["idr_split_input"]("EXP1", source) == expected


def test_idr_pseudorep_peaks_inputs(idr_paths_namespace):
    assert idr_paths_namespace["idr_pseudorep_peaks_inputs"]("EXP1", "biorep2", "1") == [
        "results/experiments/EXP1/05_pseudorep/EXP1_biorep2.pr1.bam",
        "results/experiments/EXP1/05_pseudorep/EXP1_biorep2.pr1.bam.bai",
    ]


def test_idr_pseudorep_peaks_inputs_with_prefix(idr_paths_namespace_with_control):
    assert idr_paths_namespace_with_control["idr_pseudorep_peaks_inputs"]("EXP1", "pooled", "2", source_prefix="atac_") == [
        "results/experiments/EXP1/05_pseudorep/EXP1_atac_pooled.pr2.bam",
        "results/experiments/EXP1/05_pseudorep/EXP1_atac_pooled.pr2.bam.bai",
        "results/experiments/EXP1/02_align/EXP1.pooled.control.final.bam",
    ]


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
def test_idr_macs3_args_dispatch(idr_paths_namespace_chipseq_pe, assay, peak_mode, expected_subset):
    args = idr_paths_namespace_chipseq_pe["idr_macs3_args"]("EXP1", assay, peak_mode)
    for fragment in expected_subset:
        assert fragment in args, f"{fragment!r} not in {args!r}"


def test_idr_macs3_args_se_layout(idr_paths_namespace_chipseq_se):
    assert "-f BAM" in idr_paths_namespace_chipseq_se["idr_macs3_args"]("EXP2", "chipseq", "narrow")
    assert "-g mm" in idr_paths_namespace_chipseq_se["idr_macs3_args"]("EXP2", "chipseq", "narrow")


def test_idr_macs3_args_empty_experiment(idr_paths_namespace):
    assert idr_paths_namespace["idr_macs3_args"]("NOEXP", "chipseq", "narrow") == ""


def test_idr_biorep_labels(idr_paths_namespace):
    assert idr_paths_namespace["idr_biorep_labels"]("EXP1") == ("1", "2")


# -----------------------------------------------------------------------------
# Legacy contract tests: prove the shared helpers produce the exact strings that
# the four IDR rule files previously inlined. These are intentionally verbose so
# any accidental signature drift is caught immediately.
# -----------------------------------------------------------------------------


@pytest.mark.parametrize(
    "exp,br",
    [("EXP", 1), ("EXP", 2)],
)
def test_legacy_idr_biorep_peaks_inputs_match_inlined(idr_paths_namespace, exp, br):
    """Legacy idr.smk / idr_atac.smk / idr_cuttag.smk / idr_broad.smk inlined:

        inputs = [
            idr_biorep_bam(exp, br),
            idr_biorep_bai(exp, br),
        ]
    """
    expected = [
        f"results/experiments/{exp}/02_align/biorep{br}.final.bam",
        f"results/experiments/{exp}/02_align/biorep{br}.final.bam.bai",
    ]
    assert idr_paths_namespace["idr_biorep_peaks_inputs"](exp, br) == expected


def test_legacy_idr_biorep_peaks_inputs_with_control_match_inlined(
    idr_paths_namespace_with_control,
):
    """Legacy inlined logic appended pooled control when present."""
    expected = [
        "results/experiments/EXP1/02_align/biorep2.final.bam",
        "results/experiments/EXP1/02_align/biorep2.final.bam.bai",
        "results/experiments/EXP1/02_align/EXP1.pooled.control.final.bam",
    ]
    assert idr_paths_namespace_with_control["idr_biorep_peaks_inputs"]("EXP1", 2) == expected


@pytest.mark.parametrize(
    "source,expected",
    [
        ("pooled", "results/experiments/EXP1/02_align/EXP1.pooled.final.bam"),
        ("biorep2", "results/experiments/EXP1/02_align/biorep2.final.bam"),
    ],
)
def test_legacy_idr_split_input_match_inlined(idr_paths_namespace, source, expected):
    """Legacy _split_input / _atac_split_input / _cuttag_split_input / _broad_split_input."""
    assert idr_paths_namespace["idr_split_input"]("EXP1", source) == expected


@pytest.mark.parametrize(
    "source,pr,source_prefix,expected",
    [
        (
            "biorep2",
            "1",
            "",
            [
                "results/experiments/EXP1/05_pseudorep/EXP1_biorep2.pr1.bam",
                "results/experiments/EXP1/05_pseudorep/EXP1_biorep2.pr1.bam.bai",
            ],
        ),
        (
            "pooled",
            "2",
            "atac_",
            [
                "results/experiments/EXP1/05_pseudorep/EXP1_atac_pooled.pr2.bam",
                "results/experiments/EXP1/05_pseudorep/EXP1_atac_pooled.pr2.bam.bai",
            ],
        ),
        (
            "biorep1",
            "2",
            "cuttag_",
            [
                "results/experiments/EXP1/05_pseudorep/EXP1_cuttag_biorep1.pr2.bam",
                "results/experiments/EXP1/05_pseudorep/EXP1_cuttag_biorep1.pr2.bam.bai",
            ],
        ),
        (
            "pooled",
            "1",
            "broad_chipseq_",
            [
                "results/experiments/EXP1/05_pseudorep/EXP1_broad_chipseq_pooled.pr1.bam",
                "results/experiments/EXP1/05_pseudorep/EXP1_broad_chipseq_pooled.pr1.bam.bai",
            ],
        ),
    ],
)
def test_legacy_idr_pseudorep_peaks_inputs_match_inlined(
    idr_paths_namespace, source, pr, source_prefix, expected
):
    """Legacy pseudorep input helpers for chipseq, atac, cuttag, and broad."""
    assert (
        idr_paths_namespace["idr_pseudorep_peaks_inputs"](
            "EXP1", source, pr, source_prefix=source_prefix
        )
        == expected
    )


def test_legacy_idr_pseudorep_peaks_inputs_with_control_match_inlined(
    idr_paths_namespace_with_control,
):
    """Legacy pseudorep input helpers appended pooled control when present."""
    expected = [
        "results/experiments/EXP1/05_pseudorep/EXP1_broad_cuttag_pooled.pr2.bam",
        "results/experiments/EXP1/05_pseudorep/EXP1_broad_cuttag_pooled.pr2.bam.bai",
        "results/experiments/EXP1/02_align/EXP1.pooled.control.final.bam",
    ]
    assert (
        idr_paths_namespace_with_control["idr_pseudorep_peaks_inputs"](
            "EXP1", "pooled", "2", source_prefix="broad_cuttag_"
        )
        == expected
    )


@pytest.mark.parametrize(
    "assay,peak_mode,expected_parts",
    [
        (
            "chipseq",
            "narrow",
            ["-f BAMPE", "-g hs", "-p 0.1"],
        ),
        (
            "atac",
            "narrow",
            ["-f BAMPE", "-g hs", "-p 0.1", "--nomodel", "--shift -100", "--extsize 200"],
        ),
        (
            "cuttag",
            "narrow",
            ["-f BAMPE", "-g hs", "-p 0.1", "--nomodel", "--shift -100", "--extsize 200"],
        ),
        (
            "chipseq",
            "broad",
            ["-f BAMPE", "-g hs", "-p 0.1", "--broad", "--broad-cutoff 0.1"],
        ),
        (
            "cuttag",
            "broad",
            ["-f BAMPE", "-g hs", "-p 0.1", "--broad", "--broad-cutoff 0.1"],
        ),
    ],
)
def test_legacy_idr_macs3_args_match_inlined(
    idr_paths_namespace_chipseq_pe, assay, peak_mode, expected_parts
):
    """Legacy _idr_macs3_args / _atac_idr_macs3_args / _cuttag_idr_macs3_args / _broad_idr_macs3_args."""
    args = idr_paths_namespace_chipseq_pe["idr_macs3_args"]("EXP1", assay, peak_mode)
    for part in expected_parts:
        assert part in args, f"{part!r} not in {args!r}"


def test_legacy_idr_biorep_labels_match_inlined(idr_paths_namespace):
    """Legacy idr.smk / idr_atac.smk / idr_cuttag.smk / idr_broad.smk summary rules
    previously inlined:

        bio_rep_a = str(sorted(_bioreps_for(wc.experiment, "treatment"))[0])
        bio_rep_b = str(sorted(_bioreps_for(wc.experiment, "treatment"))[1])

    or via private helpers:

        bioreps = sorted(_bioreps_for(experiment, "treatment"))
        return str(bioreps[0]), str(bioreps[1])
    """
    exp = "EXP1"
    expected_a = str(sorted(["1", "2"])[0])
    expected_b = str(sorted(["1", "2"])[1])
    labels = idr_paths_namespace["idr_biorep_labels"](exp)
    assert labels == (expected_a, expected_b)
    assert labels == ("1", "2")
