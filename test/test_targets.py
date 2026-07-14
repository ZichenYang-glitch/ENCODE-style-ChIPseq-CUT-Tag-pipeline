"""Contract tests for workflow/rules/targets.smk IDR target builders.

These tests lock the output of _idr_targets, _atac_idr_targets,
_cuttag_idr_targets, and _broad_idr_targets. They run the target builders
in a mock Snakefile namespace so the full workflow is not imported.
"""

import itertools
from pathlib import Path

import pytest


@pytest.fixture(scope="session")
def targets_file(repo_root):
    """Return the path to workflow/rules/targets.smk."""
    return Path(repo_root) / "workflow" / "rules" / "targets.smk"


def _expand(pattern, zip=False, **kwargs):
    """Lightweight reimplementation of Snakemake expand for target tests."""
    builtins_zip = (
        __builtins__["zip"] if isinstance(__builtins__, dict) else __builtins__.zip
    )

    def _as_list(value):
        if isinstance(value, str):
            return [value]
        try:
            iter(value)
        except TypeError:
            return [value]
        return list(value)

    if not kwargs:
        return [pattern]
    keys = list(kwargs.keys())
    values = [_as_list(v) for v in kwargs.values()]
    if zip:
        target_len = max(len(v) for v in values)
        values = [v * target_len if len(v) == 1 else v for v in values]
        combos = builtins_zip(*values)
    else:
        combos = itertools.product(*values)
    return [pattern.format(**dict(builtins_zip(keys, combo))) for combo in combos]


def _load_targets_namespace(targets_file, idr_globals=None):
    """Load targets.smk with IDR-related globals mocked.

    Only the IDR target builders and the globals they directly reference are
    provided. Other target builders (base, blacklist, qc, ...) are stubbed so
    _manifest_dependency_targets can be exercised if desired, but they are not
    the focus of these contract tests.
    """
    code = targets_file.read_text()
    namespace = {
        "OUTDIR": "results",
        "expand": _expand,
        "VALIDATED_CONFIG": {},
        "TREATMENT_SAMPLE_IDS": [],
        "ACTIVE_SAMPLE_IDS": [],
        "PEAK_SAMPLE_IDS": [],
        "MULTIQC": False,
        "QC_CONFIG": {},
        "BLACKLIST_SAMPLES": [],
        "SAMPLE_MAP": {},
        "SIGNAL_BW_SAMPLE_IDS": [],
        "SIGNAL_BW_EXPERIMENTS": [],
        "STAGE4B": False,
        "PEAK_MULTI_BIOREP_EXPERIMENTS": [],
        "SEACR_ENABLED": False,
        "SEACR_SAMPLE_IDS": [],
        "SEACR_MODE": "stringent",
        "CUTTAG_ACTIVE_SAMPLE_IDS": [],
        "MNASE_SAMPLE_IDS": [],
        "MNASE_MULTI_BIOREP_EXPERIMENTS": [],
        "TSS_GENOMES": [],
        "TSS_SAMPLE_IDS": [],
        "MULTI_BIOREP_EXPERIMENTS": [],
        "_EXP_LIST": [],
        "_BR_LIST": [],
        "POOLED_CONTROL_EXPERIMENTS": [],
        "CONSENSUS_ENABLED": False,
        "CONSENSUS_EXPERIMENTS": {},
        "SEACR_CONSENSUS_EXPERIMENTS": [],
        "BROAD_CHIPSEQ_IDR_EXPERIMENTS": [],
        "BROAD_CUTTAG_IDR_EXPERIMENTS": [],
        # IDR config gating and experiment lists are injected by callers.
        "STAGE5": False,
        "IDR_EXPERIMENTS": [],
        "IDR_BIOREP_EXP_LIST": [],
        "IDR_BIOREP_LIST": [],
        "IDR_PR_PEAK_EXP": [],
        "IDR_PR_PEAK_SRC": [],
        "IDR_PR_PEAK_PR": [],
        "IDR_SELF_EXP": [],
        "IDR_SELF_BR": [],
        "ATAC_IDR_ENABLED": False,
        "ATAC_IDR_EXPERIMENTS": [],
        "ATAC_IDR_BIOREP_EXP_LIST": [],
        "ATAC_IDR_BIOREP_LIST": [],
        "ATAC_IDR_PR_PEAK_EXP": [],
        "ATAC_IDR_PR_PEAK_SRC": [],
        "ATAC_IDR_PR_PEAK_PR": [],
        "ATAC_IDR_SELF_EXP": [],
        "ATAC_IDR_SELF_BR": [],
        "CUTTAG_IDR_ENABLED": False,
        "CUTTAG_IDR_EXPERIMENTS": [],
        "CUTTAG_IDR_BIOREP_EXP_LIST": [],
        "CUTTAG_IDR_BIOREP_LIST": [],
        "CUTTAG_IDR_PR_PEAK_EXP": [],
        "CUTTAG_IDR_PR_PEAK_SRC": [],
        "CUTTAG_IDR_PR_PEAK_PR": [],
        "CUTTAG_IDR_SELF_EXP": [],
        "CUTTAG_IDR_SELF_BR": [],
        "BROAD_CHIPSEQ_IDR_ENABLED": False,
        "BROAD_CUTTAG_IDR_ENABLED": False,
        "BROAD_IDR_ENABLED": False,
        "BROAD_IDR_EXPERIMENTS": [],
        "BROAD_IDR_EXPERIMENT_ASSAY": [],
        "BROAD_IDR_BIOREP_EXP_LIST": [],
        "BROAD_IDR_BIOREP_LIST": [],
        "BROAD_IDR_BIOREP_ASSAY": [],
        "BROAD_IDR_PR_PEAK_EXP": [],
        "BROAD_IDR_PR_PEAK_SRC": [],
        "BROAD_IDR_PR_PEAK_ASSAY": [],
        "BROAD_IDR_PR_PEAK_PR": [],
        "BROAD_IDR_SELF_EXP": [],
        "BROAD_IDR_SELF_BR": [],
        "BROAD_IDR_SELF_ASSAY": [],
    }
    # Stub non-IDR target builders so _manifest_dependency_targets stays loadable.
    for name in [
        "_base_targets",
        "_blacklist_targets",
        "_single_sample_qc_targets",
        "_signal_targets",
        "_cuttag_targets",
        "_mnase_targets",
        "_advanced_qc_targets",
        "_tss_targets",
        "_replicate_targets",
        "_consensus_targets",
    ]:
        namespace[name] = lambda: []
    if idr_globals:
        namespace.update(idr_globals)
    exec(compile(code, str(targets_file), "exec"), namespace)
    return namespace


@pytest.fixture
def targets_namespace_idr_enabled(targets_file):
    """Return a namespace where legacy IDR target builders are active."""
    return _load_targets_namespace(
        targets_file,
        idr_globals={
            "STAGE5": True,
            "IDR_EXPERIMENTS": ["EXP1", "EXP2"],
            "IDR_BIOREP_EXP_LIST": ["EXP1", "EXP1", "EXP2", "EXP2"],
            "IDR_BIOREP_LIST": ["1", "2", "1", "2"],
            "IDR_PR_PEAK_EXP": [
                "EXP1",
                "EXP1",
                "EXP1",
                "EXP1",
                "EXP1",
                "EXP1",
                "EXP2",
                "EXP2",
                "EXP2",
                "EXP2",
                "EXP2",
                "EXP2",
            ],
            "IDR_PR_PEAK_SRC": [
                "biorep1",
                "biorep1",
                "biorep2",
                "biorep2",
                "pooled",
                "pooled",
                "biorep1",
                "biorep1",
                "biorep2",
                "biorep2",
                "pooled",
                "pooled",
            ],
            "IDR_PR_PEAK_PR": [
                "1",
                "2",
                "1",
                "2",
                "1",
                "2",
                "1",
                "2",
                "1",
                "2",
                "1",
                "2",
            ],
            "IDR_SELF_EXP": ["EXP1", "EXP1", "EXP2", "EXP2"],
            "IDR_SELF_BR": ["1", "2", "1", "2"],
        },
    )


@pytest.fixture
def targets_namespace_atac_idr_enabled(targets_file):
    """Return a namespace where ATAC IDR target builders are active."""
    return _load_targets_namespace(
        targets_file,
        idr_globals={
            "ATAC_IDR_ENABLED": True,
            "ATAC_IDR_EXPERIMENTS": ["ATAC1"],
            "ATAC_IDR_BIOREP_EXP_LIST": ["ATAC1", "ATAC1"],
            "ATAC_IDR_BIOREP_LIST": ["1", "2"],
            "ATAC_IDR_PR_PEAK_EXP": ["ATAC1"] * 6,
            "ATAC_IDR_PR_PEAK_SRC": [
                "biorep1",
                "biorep1",
                "biorep2",
                "biorep2",
                "pooled",
                "pooled",
            ],
            "ATAC_IDR_PR_PEAK_PR": ["1", "2", "1", "2", "1", "2"],
            "ATAC_IDR_SELF_EXP": ["ATAC1", "ATAC1"],
            "ATAC_IDR_SELF_BR": ["1", "2"],
        },
    )


@pytest.fixture
def targets_namespace_cuttag_idr_enabled(targets_file):
    """Return a namespace where CUT&Tag IDR target builders are active."""
    return _load_targets_namespace(
        targets_file,
        idr_globals={
            "CUTTAG_IDR_ENABLED": True,
            "CUTTAG_IDR_EXPERIMENTS": ["CUT1"],
            "CUTTAG_IDR_BIOREP_EXP_LIST": ["CUT1", "CUT1"],
            "CUTTAG_IDR_BIOREP_LIST": ["1", "2"],
            "CUTTAG_IDR_PR_PEAK_EXP": ["CUT1"] * 6,
            "CUTTAG_IDR_PR_PEAK_SRC": [
                "biorep1",
                "biorep1",
                "biorep2",
                "biorep2",
                "pooled",
                "pooled",
            ],
            "CUTTAG_IDR_PR_PEAK_PR": ["1", "2", "1", "2", "1", "2"],
            "CUTTAG_IDR_SELF_EXP": ["CUT1", "CUT1"],
            "CUTTAG_IDR_SELF_BR": ["1", "2"],
        },
    )


@pytest.fixture
def targets_namespace_broad_idr_enabled(targets_file):
    """Return a namespace where broad IDR target builders are active."""
    return _load_targets_namespace(
        targets_file,
        idr_globals={
            "BROAD_CHIPSEQ_IDR_ENABLED": True,
            "BROAD_CUTTAG_IDR_ENABLED": True,
            "BROAD_IDR_ENABLED": True,
            "BROAD_IDR_EXPERIMENTS": ["BROAD1", "BROAD2"],
            "BROAD_IDR_EXPERIMENT_ASSAY": ["chipseq", "cuttag"],
            "BROAD_IDR_BIOREP_EXP_LIST": ["BROAD1", "BROAD1", "BROAD2", "BROAD2"],
            "BROAD_IDR_BIOREP_LIST": ["1", "2", "1", "2"],
            "BROAD_IDR_BIOREP_ASSAY": ["chipseq", "chipseq", "cuttag", "cuttag"],
            "BROAD_IDR_PR_PEAK_EXP": ["BROAD1"] * 6 + ["BROAD2"] * 6,
            "BROAD_IDR_PR_PEAK_SRC": [
                "biorep1",
                "biorep1",
                "biorep2",
                "biorep2",
                "pooled",
                "pooled",
                "biorep1",
                "biorep1",
                "biorep2",
                "biorep2",
                "pooled",
                "pooled",
            ],
            "BROAD_IDR_PR_PEAK_ASSAY": ["chipseq"] * 6 + ["cuttag"] * 6,
            "BROAD_IDR_PR_PEAK_PR": ["1", "2"] * 6,
            "BROAD_IDR_SELF_EXP": ["BROAD1", "BROAD1", "BROAD2", "BROAD2"],
            "BROAD_IDR_SELF_BR": ["1", "2", "1", "2"],
            "BROAD_IDR_SELF_ASSAY": ["chipseq", "chipseq", "cuttag", "cuttag"],
        },
    )


def _legacy_idr_expected():
    """Return expected legacy IDR target list for the idr_enabled fixture."""
    outdir = "results"
    expected = []
    for exp, br in zip(["EXP1", "EXP1", "EXP2", "EXP2"], ["1", "2", "1", "2"]):
        expected.append(
            f"{outdir}/experiments/{exp}/04_peaks/idr/"
            f"{exp}_biorep{br}_idr_peaks.narrowPeak"
        )
    for exp in ["EXP1", "EXP2"]:
        expected.append(f"{outdir}/experiments/{exp}/06_idr/true_replicates/idr.txt")
    for exp in ["EXP1", "EXP2"]:
        expected.append(
            f"{outdir}/experiments/{exp}/06_idr/true_replicates/"
            f"idr.thresholded.narrowPeak"
        )
    pr_data = zip(
        ["EXP1"] * 6 + ["EXP2"] * 6,
        ["biorep1", "biorep1", "biorep2", "biorep2", "pooled", "pooled"] * 2,
        ["1", "2"] * 6,
    )
    for exp, src, pr in pr_data:
        expected.append(
            f"{outdir}/experiments/{exp}/04_peaks/idr/"
            f"{exp}_{src}_pr{pr}_idr_peaks.narrowPeak"
        )
    for exp, br in zip(["EXP1", "EXP1", "EXP2", "EXP2"], ["1", "2", "1", "2"]):
        expected.append(
            f"{outdir}/experiments/{exp}/06_idr/self_pseudoreplicates/"
            f"biorep{br}.idr.txt"
        )
    for exp, br in zip(["EXP1", "EXP1", "EXP2", "EXP2"], ["1", "2", "1", "2"]):
        expected.append(
            f"{outdir}/experiments/{exp}/06_idr/self_pseudoreplicates/"
            f"biorep{br}.idr.thresholded.narrowPeak"
        )
    for exp in ["EXP1", "EXP2"]:
        expected.append(
            f"{outdir}/experiments/{exp}/06_idr/pooled_pseudoreplicates/idr.txt"
        )
    for exp in ["EXP1", "EXP2"]:
        expected.append(
            f"{outdir}/experiments/{exp}/06_idr/pooled_pseudoreplicates/"
            f"idr.thresholded.narrowPeak"
        )
    for exp in ["EXP1", "EXP2"]:
        expected.append(
            f"{outdir}/experiments/{exp}/06_idr/final/reproducibility_summary.tsv"
        )
    for exp in ["EXP1", "EXP2"]:
        expected.append(
            f"{outdir}/experiments/{exp}/06_idr/final/conservative.narrowPeak"
        )
    for exp in ["EXP1", "EXP2"]:
        expected.append(f"{outdir}/experiments/{exp}/06_idr/final/optimal.narrowPeak")
    return expected


def test_idr_targets_returns_expected_list(targets_namespace_idr_enabled):
    targets = targets_namespace_idr_enabled["_idr_targets"]()
    expected = _legacy_idr_expected()
    assert isinstance(targets, list)
    assert len(targets) == len(expected)
    assert targets == expected
    assert len(targets) == len(set(targets))


def test_idr_targets_key_markers(targets_namespace_idr_enabled):
    targets = set(targets_namespace_idr_enabled["_idr_targets"]())
    assert (
        "results/experiments/EXP1/04_peaks/idr/EXP1_biorep1_idr_peaks.narrowPeak"
        in targets
    )
    assert "results/experiments/EXP1/06_idr/true_replicates/idr.txt" in targets
    assert "results/experiments/EXP1/06_idr/final/optimal.narrowPeak" in targets


def test_atac_idr_targets_returns_expected_list(targets_namespace_atac_idr_enabled):
    targets = targets_namespace_atac_idr_enabled["_atac_idr_targets"]()
    outdir = "results"
    expected = []
    expected.append(
        f"{outdir}/experiments/ATAC1/06_reproducibility/idr/idr_peaks/"
        f"ATAC1_atac_biorep1_idr.narrowPeak"
    )
    expected.append(
        f"{outdir}/experiments/ATAC1/06_reproducibility/idr/idr_peaks/"
        f"ATAC1_atac_biorep2_idr.narrowPeak"
    )
    expected.append(
        f"{outdir}/experiments/ATAC1/06_reproducibility/idr/"
        f"true_replicates/ATAC1_atac_idr.txt"
    )
    expected.append(
        f"{outdir}/experiments/ATAC1/06_reproducibility/idr/"
        f"true_replicates/ATAC1_atac_idr.thresholded.narrowPeak"
    )
    pr_pairs = [
        ("biorep1", "1"),
        ("biorep1", "2"),
        ("biorep2", "1"),
        ("biorep2", "2"),
        ("pooled", "1"),
        ("pooled", "2"),
    ]
    for src, pr in pr_pairs:
        expected.append(
            f"{outdir}/experiments/ATAC1/06_reproducibility/idr/idr_peaks/"
            f"ATAC1_atac_{src}_pr{pr}_idr.narrowPeak"
        )
    for br in ["1", "2"]:
        expected.append(
            f"{outdir}/experiments/ATAC1/06_reproducibility/idr/"
            f"self_pseudoreplicates/ATAC1_atac_biorep{br}_idr.txt"
        )
    for br in ["1", "2"]:
        expected.append(
            f"{outdir}/experiments/ATAC1/06_reproducibility/idr/"
            f"self_pseudoreplicates/ATAC1_atac_biorep{br}_idr.thresholded.narrowPeak"
        )
    expected.append(
        f"{outdir}/experiments/ATAC1/06_reproducibility/idr/"
        f"pooled_pseudoreplicates/ATAC1_atac_idr.txt"
    )
    expected.append(
        f"{outdir}/experiments/ATAC1/06_reproducibility/idr/"
        f"pooled_pseudoreplicates/ATAC1_atac_idr.thresholded.narrowPeak"
    )
    expected.append(
        f"{outdir}/experiments/ATAC1/06_reproducibility/final/"
        f"ATAC1.atac.macs3.narrow.replicate_validated.idr.narrowPeak"
    )
    expected.append(
        f"{outdir}/experiments/ATAC1/06_reproducibility/final/"
        f"reproducibility_summary.tsv"
    )
    assert isinstance(targets, list)
    assert targets == expected
    assert len(targets) == len(set(targets))


def test_cuttag_idr_targets_returns_expected_list(targets_namespace_cuttag_idr_enabled):
    targets = targets_namespace_cuttag_idr_enabled["_cuttag_idr_targets"]()
    outdir = "results"
    expected = []
    expected.append(
        f"{outdir}/experiments/CUT1/06_reproducibility/idr/idr_peaks/"
        f"CUT1_cuttag_biorep1_idr.narrowPeak"
    )
    expected.append(
        f"{outdir}/experiments/CUT1/06_reproducibility/idr/idr_peaks/"
        f"CUT1_cuttag_biorep2_idr.narrowPeak"
    )
    expected.append(
        f"{outdir}/experiments/CUT1/06_reproducibility/idr/"
        f"true_replicates/CUT1_cuttag_idr.txt"
    )
    expected.append(
        f"{outdir}/experiments/CUT1/06_reproducibility/idr/"
        f"true_replicates/CUT1_cuttag_idr.thresholded.narrowPeak"
    )
    pr_pairs = [
        ("biorep1", "1"),
        ("biorep1", "2"),
        ("biorep2", "1"),
        ("biorep2", "2"),
        ("pooled", "1"),
        ("pooled", "2"),
    ]
    for src, pr in pr_pairs:
        expected.append(
            f"{outdir}/experiments/CUT1/06_reproducibility/idr/idr_peaks/"
            f"CUT1_cuttag_{src}_pr{pr}_idr.narrowPeak"
        )
    for br in ["1", "2"]:
        expected.append(
            f"{outdir}/experiments/CUT1/06_reproducibility/idr/"
            f"self_pseudoreplicates/CUT1_cuttag_biorep{br}_idr.txt"
        )
    for br in ["1", "2"]:
        expected.append(
            f"{outdir}/experiments/CUT1/06_reproducibility/idr/"
            f"self_pseudoreplicates/CUT1_cuttag_biorep{br}_idr.thresholded.narrowPeak"
        )
    expected.append(
        f"{outdir}/experiments/CUT1/06_reproducibility/idr/"
        f"pooled_pseudoreplicates/CUT1_cuttag_idr.txt"
    )
    expected.append(
        f"{outdir}/experiments/CUT1/06_reproducibility/idr/"
        f"pooled_pseudoreplicates/CUT1_cuttag_idr.thresholded.narrowPeak"
    )
    expected.append(
        f"{outdir}/experiments/CUT1/06_reproducibility/final/"
        f"CUT1.cuttag.macs3.narrow.replicate_validated.idr.narrowPeak"
    )
    expected.append(
        f"{outdir}/experiments/CUT1/06_reproducibility/final/"
        f"reproducibility_summary.tsv"
    )
    assert isinstance(targets, list)
    assert targets == expected
    assert len(targets) == len(set(targets))


def test_broad_idr_targets_returns_expected_list(targets_namespace_broad_idr_enabled):
    targets = targets_namespace_broad_idr_enabled["_broad_idr_targets"]()
    outdir = "results"
    expected = []
    for exp, assay, br in zip(
        ["BROAD1", "BROAD1", "BROAD2", "BROAD2"],
        ["chipseq", "chipseq", "cuttag", "cuttag"],
        ["1", "2", "1", "2"],
    ):
        expected.append(
            f"{outdir}/experiments/{exp}/06_reproducibility/idr/idr_peaks/"
            f"{exp}_broad_{assay}_biorep{br}_idr.broadPeak"
        )
    for exp, assay in zip(["BROAD1", "BROAD2"], ["chipseq", "cuttag"]):
        expected.append(
            f"{outdir}/experiments/{exp}/06_reproducibility/idr/"
            f"true_replicates/{exp}_broad_{assay}_idr.txt"
        )
    for exp, assay in zip(["BROAD1", "BROAD2"], ["chipseq", "cuttag"]):
        expected.append(
            f"{outdir}/experiments/{exp}/06_reproducibility/idr/"
            f"true_replicates/{exp}_broad_{assay}_idr.thresholded.broadPeak"
        )
    pr_data = zip(
        ["BROAD1"] * 6 + ["BROAD2"] * 6,
        ["chipseq"] * 6 + ["cuttag"] * 6,
        ["biorep1", "biorep1", "biorep2", "biorep2", "pooled", "pooled"] * 2,
        ["1", "2"] * 6,
    )
    for exp, assay, src, pr in pr_data:
        expected.append(
            f"{outdir}/experiments/{exp}/06_reproducibility/idr/idr_peaks/"
            f"{exp}_broad_{assay}_{src}_pr{pr}_idr.broadPeak"
        )
    for exp, assay, br in zip(
        ["BROAD1", "BROAD1", "BROAD2", "BROAD2"],
        ["chipseq", "chipseq", "cuttag", "cuttag"],
        ["1", "2", "1", "2"],
    ):
        expected.append(
            f"{outdir}/experiments/{exp}/06_reproducibility/idr/"
            f"self_pseudoreplicates/{exp}_broad_{assay}_biorep{br}_idr.txt"
        )
    for exp, assay, br in zip(
        ["BROAD1", "BROAD1", "BROAD2", "BROAD2"],
        ["chipseq", "chipseq", "cuttag", "cuttag"],
        ["1", "2", "1", "2"],
    ):
        expected.append(
            f"{outdir}/experiments/{exp}/06_reproducibility/idr/"
            f"self_pseudoreplicates/{exp}_broad_{assay}_biorep{br}_idr.thresholded.broadPeak"
        )
    for exp, assay in zip(["BROAD1", "BROAD2"], ["chipseq", "cuttag"]):
        expected.append(
            f"{outdir}/experiments/{exp}/06_reproducibility/idr/"
            f"pooled_pseudoreplicates/{exp}_broad_{assay}_idr.txt"
        )
    for exp, assay in zip(["BROAD1", "BROAD2"], ["chipseq", "cuttag"]):
        expected.append(
            f"{outdir}/experiments/{exp}/06_reproducibility/idr/"
            f"pooled_pseudoreplicates/{exp}_broad_{assay}_idr.thresholded.broadPeak"
        )
    for exp, assay in zip(["BROAD1", "BROAD2"], ["chipseq", "cuttag"]):
        expected.append(
            f"{outdir}/experiments/{exp}/06_reproducibility/final/"
            f"{exp}.{assay}.macs3.broad.replicate_validated.idr.broadPeak"
        )
    for exp in ["BROAD1", "BROAD2"]:
        expected.append(
            f"{outdir}/experiments/{exp}/06_reproducibility/final/"
            f"reproducibility_summary.tsv"
        )
    assert isinstance(targets, list)
    assert targets == expected
    assert len(targets) == len(set(targets))


def test_idr_targets_empty_when_disabled(targets_file):
    namespace = _load_targets_namespace(targets_file)
    assert namespace["_idr_targets"]() == []
    assert namespace["_atac_idr_targets"]() == []
    assert namespace["_cuttag_idr_targets"]() == []
    assert namespace["_broad_idr_targets"]() == []
