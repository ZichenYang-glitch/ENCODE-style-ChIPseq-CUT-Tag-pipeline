"""Contracts between MNase path helpers and the artifact catalog."""

from pathlib import Path

import pytest

from encode_pipeline.artifacts import artifacts_by_id, load_catalog


REPO_ROOT = Path(__file__).resolve().parents[2]
INVENTORY_PATH = REPO_ROOT / "docs" / "architecture" / "artifact-inventory.yaml"
PATH_HELPERS = REPO_ROOT / "workflow" / "rules" / "paths.smk"


@pytest.fixture(scope="module")
def mnase_path_helpers():
    namespace = {"OUTDIR": "results"}
    source = PATH_HELPERS.read_text(encoding="utf-8")
    exec(compile(source, PATH_HELPERS, "exec"), namespace)
    return namespace


@pytest.mark.parametrize(
    "artifact_id,helper_name,args",
    [
        ("mnase_mono_bam", "mnase_fragment_bam", ("SAMPLE", "mono")),
        ("mnase_mono_bai", "mnase_fragment_bai", ("SAMPLE", "mono")),
        ("mnase_sub_bam", "mnase_fragment_bam", ("SAMPLE", "sub")),
        ("mnase_sub_bai", "mnase_fragment_bai", ("SAMPLE", "sub")),
        ("mnase_di_bam", "mnase_fragment_bam", ("SAMPLE", "di")),
        ("mnase_di_bai", "mnase_fragment_bai", ("SAMPLE", "di")),
        ("mnase_dyad_bigwig", "mnase_signal_bw", ("SAMPLE", "dyad")),
        ("mnase_mono_bigwig", "mnase_signal_bw", ("SAMPLE", "mono")),
        ("mnase_qc_summary", "mnase_qc_summary_tsv", ("SAMPLE",)),
        (
            "pooled_mnase_mono_bam",
            "mnase_pooled_fragment_bam",
            ("EXPERIMENT", "mono"),
        ),
        (
            "pooled_mnase_mono_bai",
            "mnase_pooled_fragment_bai",
            ("EXPERIMENT", "mono"),
        ),
        (
            "pooled_mnase_dyad_bigwig",
            "mnase_pooled_signal_bw",
            ("EXPERIMENT", "dyad"),
        ),
        (
            "pooled_mnase_mono_bigwig",
            "mnase_pooled_signal_bw",
            ("EXPERIMENT", "mono"),
        ),
    ],
)
def test_mnase_path_helper_matches_catalog_template(
    mnase_path_helpers, artifact_id, helper_name, args
):
    catalog = artifacts_by_id(load_catalog(str(INVENTORY_PATH)))
    artifact = catalog[artifact_id]
    helper = mnase_path_helpers[helper_name]
    expected = artifact.path_template.replace("<sample>", "SAMPLE").replace(
        "<experiment>", "EXPERIMENT"
    )

    assert callable(helper)
    assert helper(*args) == expected
