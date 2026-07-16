"""Identity and scope tests for the pinned bulk RNA-seq results contract."""

from __future__ import annotations

import hashlib
from pathlib import Path

from encode_pipeline.adapters.bulk_rnaseq.results_contract import (
    RESULTS_CONTRACT_FILE,
    RESULTS_CONTRACT_SHA256,
    RESULTS_CONTRACT_SIZE,
    effective_rseqc_modules,
    load_bulk_rnaseq_results_contract,
)


PROJECT_ROOT = Path(__file__).resolve().parents[2]


def test_results_contract_has_exact_identity_and_fresh_loads():
    path = (
        PROJECT_ROOT
        / "src/encode_pipeline/contracts/nfcore_rnaseq"
        / RESULTS_CONTRACT_FILE
    )
    content = path.read_bytes()

    assert len(content) == RESULTS_CONTRACT_SIZE
    assert hashlib.sha256(content).hexdigest() == RESULTS_CONTRACT_SHA256
    first = load_bulk_rnaseq_results_contract()
    second = load_bulk_rnaseq_results_contract()
    assert first == second
    first["artifact_families"].clear()
    assert load_bulk_rnaseq_results_contract() == second


def test_results_contract_records_fixed_route_and_explicit_exclusions():
    contract = load_bulk_rnaseq_results_contract()

    assert contract["upstream"] == {
        "project": "nf-core/rnaseq",
        "release": "3.26.0",
        "commit": "e7ca46272c8f9d5ceee3f71759f4ba551d3217a4",
        "route": "star_salmon",
    }
    exclusions = {entry["path"] for entry in contract["public_exclusions"]}
    assert "multiqc/star_salmon/multiqc_report.html" in exclusions
    assert "multiqc/star_salmon/multiqc_report_data/multiqc_data.json" in exclusions
    assert "star_salmon/<sample>/cmd_info.json" in exclusions
    assert "trimgalore/*_trimming_report.txt" in exclusions
    assert "star_salmon/featurecounts/<sample>.featureCounts.tsv" in exclusions
    assert "workspace/logs/nextflow.log" in exclusions
    assert contract["explicit_profile_policy"]["save_reference"].startswith(
        "fail_closed"
    )
    assert contract["explicit_profile_policy"]["stringtie"].startswith("fail_closed")
    assert (
        contract["explicit_profile_policy"]["deseq2_qc"]
        == "excluded_qc_visualization_not_differential_expression"
    )
    assert (
        "exact_final_retained_bases_require_a_future_sanitized_derivative"
        in contract["explicit_profile_policy"]["trimgalore_metrics"]
    )
    assert "159" in contract["explicit_profile_policy"]["multiqc_sample_identity"]
    assert contract["multiqc_sample_identity_sources"] == {
        "multiqc_config_defaults_sha256": (
            "037598900d99fb4a5a32aeea2afffa4702fb2b4c309a05819fd5b1c655ca55de"
        ),
        "nfcore_multiqc_config_path": (
            "workflows/rnaseq/assets/multiqc/multiqc_config.yml"
        ),
        "nfcore_multiqc_config_size": 8940,
        "nfcore_multiqc_config_sha256": (
            "d859d137f2aa7800f7e12a84dec6658b1f558197e08f3041d05cd69fa6aa0b7e"
        ),
    }
    assert contract["bounded_limits"] == {
        "artifact_candidates": 128128,
        "audited_namespace_entries": 128128,
        "sample_status_table_bytes": 194128,
        "qc_source_files": 16016,
        "qc_source_bytes_each": 16777216,
        "qc_source_bytes_total": 268435456,
        "qc_metrics": 128000,
        "json_depth": 16,
        "json_nodes": 50000,
        "json_string_characters": 8192,
    }


def test_csi_effective_rseqc_policy_removes_fixed_incompatible_modules():
    default_modules = effective_rseqc_modules({"bam_csi_index": True})
    explicit_modules = effective_rseqc_modules(
        {
            "bam_csi_index": True,
            "rseqc_modules": "bam_stat,inner_distance,read_distribution,tin",
        }
    )

    assert {"inner_distance", "read_distribution", "tin"}.isdisjoint(default_modules)
    assert explicit_modules == ("bam_stat",)
