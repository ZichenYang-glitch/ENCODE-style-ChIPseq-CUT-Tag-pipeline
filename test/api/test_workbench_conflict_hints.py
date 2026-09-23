"""Real HTTP validation projection, reference resolution and correction hints."""

from copy import deepcopy
from pathlib import Path
import runpy

import pytest

from api_test_client import ApiTestClient
from conftest import seed_test_authentication


@pytest.fixture
def validation_app(tmp_path, monkeypatch):
    fixture = runpy.run_path(
        str(
            Path(__file__).resolve().parents[1]
            / "browser/workbench_validation_runtime.py"
        )
    )
    config = fixture["prepare_reference"](tmp_path)
    monkeypatch.setenv("ENCODE_PIPELINE_REFERENCE_PROFILE_CONFIG", str(config))
    app = fixture["create_validation_app"](tmp_path)
    seed_test_authentication(app)
    try:
        yield app
    finally:
        app.state.run_queue.close()
        app.state.persistence.close()


@pytest.mark.parametrize(
    "config,path,hint_part,correction",
    [
        (
            {
                "standard": {
                    "umi": {"enabled": False},
                    "outputs": {"umi_intermediates": True},
                }
            },
            "config.standard.outputs",
            "Turn off standard.outputs.umi_intermediates",
            "umi_intermediates",
        ),
        (
            {
                "standard": {
                    "trimming": {"enabled": False, "tool": "trimgalore"},
                    "outputs": {"trimmed_reads": True},
                }
            },
            "config.standard.outputs",
            "Turn off standard.outputs.trimmed_reads",
            "trimmed_reads",
        ),
        (
            {
                "standard": {"qc": {"enabled": False}},
                "advanced": {"rseqc_modules": "bam_stat"},
            },
            "config.advanced.rseqc_modules",
            "enable both standard.qc.enabled and standard.qc.rseqc",
            "rseqc_modules",
        ),
        (
            {"standard": {"qc": {"enabled": False}}, "advanced": {"deseq2_vst": False}},
            "config.advanced.deseq2_vst",
            "Setting deseq2_vst to false does not remove it",
            "deseq2_vst",
        ),
        (
            {
                "standard": {"qc": {"enabled": False}},
                "advanced": {"featurecounts_group_type": "gene_type"},
            },
            "config.advanced",
            "Remove advanced.featurecounts_group_type",
            "featurecounts_group_type",
        ),
        (
            {
                "standard": {"trimming": {"enabled": False, "tool": "trimgalore"}},
                "advanced": {"min_trimmed_reads": 1},
            },
            "config.advanced.min_trimmed_reads",
            "enable standard.trimming.enabled",
            "min_trimmed_reads",
        ),
    ],
)
def test_real_http_conflict_hints_and_explicit_correction(
    validation_app,
    tmp_path,
    config,
    path,
    hint_part,
    correction,
):
    before = deepcopy(config)
    payload = {
        "config": config,
        "samples": [
            {
                "sample": "S1",
                "library": "lib1",
                "lane": "L001",
                "layout": "SE",
                "fastq_1": str(tmp_path / "reads.fastq.gz"),
                "strandedness": "reverse",
                "platform": "ILLUMINA",
            }
        ],
        "options": {},
        "reference_profile_revision_id": validation_app.state.test_reference_profile.revision_id,
    }
    with ApiTestClient(validation_app) as client:
        response = client.post("/api/v1/workflows/bulk-rnaseq/validate", json=payload)
        assert response.status_code == 200
        result = response.json()
        assert result["ok"] is False
        expected_code = (
            "BULK_RNASEQ_OUTPUT_CONFLICT"
            if "outputs" in path
            else "BULK_RNASEQ_ADVANCED_CONTEXT_CONFLICT"
        )
        assert [issue["code"] for issue in result["issues"]] == [expected_code]
        assert result["issues"][0]["path"] == path
        assert hint_part in result["issues"][0]["hint"]
        assert str(tmp_path) not in response.text
        assert config == before
        corrected = deepcopy(payload)
        if "outputs" in path:
            corrected["config"]["standard"]["outputs"][correction] = False
        else:
            del corrected["config"]["advanced"][correction]
        accepted = client.post("/api/v1/workflows/bulk-rnaseq/validate", json=corrected)
        assert accepted.status_code == 200
        assert accepted.json()["ok"] is True, accepted.json()
        # A validation-only fixture must not pretend to issue an executable snapshot.
        assert accepted.json()["snapshot"] is None
