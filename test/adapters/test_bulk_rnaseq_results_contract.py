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
from encode_pipeline.platform.adapters import MAX_SAMPLE_ROWS


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
    assert len(exclusions) == len(contract["public_exclusions"])
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
        "future sanitized derivative"
        in contract["explicit_profile_policy"]["trimgalore_metrics"]
    )
    assert "159" in contract["explicit_profile_policy"]["multiqc_sample_identity"]
    assert (
        contract["explicit_profile_policy"]["rseqc_tin_sample_identity"]
        == "effective TIN rejects case-sensitive lowercase bam in sample IDs before workspace planning; authored IDs are never rewritten"
    )
    inner_distance = contract["artifact_semantics"]["rseqc_inner_distance"]
    assert inner_distance["distance"]["output_type"] == (
        "bulk_rnaseq.rseqc.inner_distance.distance"
    )
    assert inner_distance["distance"]["semantic"] == "per_read_pair_detail"
    assert inner_distance["distance"]["layout"] == "paired_end_only"
    assert inner_distance["mean"]["path"] == (
        "star_salmon/rseqc/inner_distance/txt/<sample>.inner_distance_mean.txt"
    )
    assert inner_distance["mean"]["semantic"] == (
        "mean_median_standard_deviation_summary"
    )
    assert inner_distance["single_end"] == {
        "path": "star_salmon/rseqc/inner_distance/txt/<sample>.inner_distance.txt",
        "semantic": "fixed_upstream_unsupported_placeholder",
        "public_policy": "excluded_but_namespace_owner_audited",
    }
    salmon = contract["metric_semantics"]["salmon_meta_info"]
    assert salmon["num_processed"]["metric_key"] == "salmon.processed_fragments"
    assert salmon["num_processed"]["display_name"] == "Salmon processed fragments"
    assert salmon["num_processed"]["unit"] == "count"
    assert salmon["num_mapped"]["metric_key"] == "salmon.mapped_fragments"
    assert salmon["num_mapped"]["display_name"] == "Salmon mapped fragments"
    assert salmon["mapping_fraction"] == {
        "metric_key": "salmon.mapping_fraction",
        "display_name": "Salmon mapped fragment fraction",
        "semantic": "mapped_fragments_divided_by_processed_fragments",
        "unit": "fraction",
    }
    assert salmon["paired_end_fragment"] == "one_read_pair"
    assert salmon["single_end_fragment"] == "one_single_read"
    assert contract["upstream_behavior_evidence"]["rseqc_tin"]["version"] == "5.0.4"
    assert contract["upstream_behavior_evidence"]["rseqc_tin"]["stdev_max"] == 50
    tin_evidence = contract["upstream_behavior_evidence"]["rseqc_tin"]
    assert tin_evidence["sdist_source"] == {
        "path": "scripts/tin.py",
        "size": 11563,
        "sha256": ("f2225494b2ea71461f969b6934ce9bdb8a813894c909e5989151c91db0f5e40b"),
    }
    assert tin_evidence["wheel_script"] == {
        "path": "RSeQC-5.0.4.data/scripts/tin.py",
        "size": 11550,
        "sha256": ("fe8514d74a5a9b467105823ed2792b91818c37de68ddd4eaae68edb712dc2fa6"),
    }
    star = contract["metric_semantics"]["star_log_final"]
    assert star["paired_end_template"] == "one_read_pair"
    assert star["metrics"]["pure_unmapped_fraction"] == {
        "metric_key": "star.pure_unmapped_template_fraction",
        "display_name": "STAR pure-unmapped template fraction",
        "unit": "fraction",
    }
    star_evidence = contract["upstream_behavior_evidence"]["star_log_final"]
    assert star_evidence["stats_source"]["sha256"] == (
        "8120977869f2360d1486f3f5f3e52e6032f53ce255e569ac86026047fa12f92f"
    )
    assert star_evidence["read_alignment_source"]["sha256"] == (
        "76aa9321c3cf63cc45acd8c16d0b841faec1b9b0aac3e78bbdea865867da9538"
    )
    featurecounts = contract["metric_semantics"]["featurecounts_summary"]
    assert featurecounts["count_unit"] == "read"
    assert featurecounts["metrics"] == {
        "assigned_reads": {
            "metric_key": "featurecounts.assigned_reads",
            "display_name": "featureCounts assigned reads",
            "unit": "count",
        },
        "classified_reads": {
            "metric_key": "featurecounts.classified_reads",
            "display_name": "featureCounts classified reads",
            "unit": "count",
        },
        "assigned_read_fraction": {
            "metric_key": "featurecounts.assigned_read_fraction",
            "display_name": "featureCounts assigned read fraction",
            "unit": "fraction",
        },
    }
    featurecounts_evidence = contract["upstream_behavior_evidence"][
        "featurecounts_summary"
    ]
    assert featurecounts_evidence["nfcore_config_sha256"] == (
        "cd2f70a786aa0c5e4a26e4cbf7342f044fe483324a231dc9011445db5cb7b880"
    )
    assert featurecounts_evidence["subread_source_sha256"] == (
        "470ffb8e1acbd6255a4ca09ec0f4ce587e3a6071f441ccd6c884ead238dff938"
    )
    salmon_evidence = contract["upstream_behavior_evidence"]["salmon_meta_info"]
    assert salmon_evidence["counter_source_sha256"] == (
        "de8ffc555de5856d8b0afb079841bc616a5babd7d220bcaef679c825df016129"
    )
    assert contract["upstream_behavior_evidence"]["sample_status_filtering"][
        "cutadapt_percent_tolerance"
    ] == {
        "value": "0.000000001",
        "unit": "percentage_points",
        "comparison": "absolute_reported_minus_exact_base_ratio",
    }
    assert contract["multiqc_sample_identity_sources"] == {
        "multiqc_version": "1.33",
        "multiqc_package_sha256": (
            "06e2b82fd9bfa458a79d8da869cb7d88260f624c03297120a70446a49ce852ed"
        ),
        "multiqc_config_defaults_sha256": (
            "037598900d99fb4a5a32aeea2afffa4702fb2b4c309a05819fd5b1c655ca55de"
        ),
        "filename_cleaning_order": (
            "prepended_extra_fn_clean_exts_then_default_fn_clean_exts_then_"
            "fn_clean_trim_then_global_exact_name_replacement"
        ),
        "fn_clean_trim_count": 26,
        "nfcore_multiqc_config_path": (
            "workflows/rnaseq/assets/multiqc/multiqc_config.yml"
        ),
        "nfcore_multiqc_config_size": 8940,
        "nfcore_multiqc_config_sha256": (
            "d859d137f2aa7800f7e12a84dec6658b1f558197e08f3041d05cd69fa6aa0b7e"
        ),
        "nfcore_multiqc_subworkflow_path": (
            "subworkflows/local/multiqc_rnaseq/main.nf"
        ),
        "nfcore_multiqc_subworkflow_size": 19937,
        "nfcore_multiqc_subworkflow_sha256": (
            "e9a3ba3e2823c726d01421f6b3f30a5479d3d289e59573c1a6f541c09d42df34"
        ),
        "multiqc_base_module_path": "multiqc/base_module.py",
        "multiqc_base_module_size": 59366,
        "multiqc_base_module_sha256": (
            "abc1fa1b848e8513e3e2ab4affcea4e6b4eec034ec2f5b30b5a4c1aa5f3073e4"
        ),
        "multiqc_table_object_path": "multiqc/plots/table_object.py",
        "multiqc_table_object_size": 60417,
        "multiqc_table_object_sha256": (
            "8d152780a1d761cf5409e74bce62b2add16e16648fe9f42258697f60813752a3"
        ),
        "multiqc_report_path": "multiqc/report.py",
        "multiqc_report_size": 48562,
        "multiqc_report_sha256": (
            "7d30a943b4611b2eebc14422f3c11aeb3a4fba83783bf665ca91ad7d6da3c046"
        ),
        "general_stats_grouping": (
            "exact per-sample union of canonical S and, only on fixed paired "
            "routes with authored raw FastQC or paired Trim Galore evidence, "
            "grouped S Read 1 and S Read 2 rows; threshold failures do not "
            "remove pre-threshold identities"
        ),
    }
    assert contract["bounded_limits"] == {
        "artifact_candidates": 128128,
        "audited_namespace_entries": 128128,
        "sample_status_table_bytes": 194128,
        "published_multiqc_table_data_rows": MAX_SAMPLE_ROWS * 3,
        "published_multiqc_table_bytes_each": 16777216,
        "published_multiqc_table_bytes_total": 67108864,
        "sample_status_evidence_bytes_each": 16777216,
        "sample_status_evidence_bytes_total": 268435456,
        "sample_status_evidence_line_characters": 8192,
        "sample_status_numeric_characters": 64,
        "qc_source_files": 16016,
        "qc_source_bytes_each": 16777216,
        "qc_source_bytes_total": 268435456,
        "qc_metrics": 128000,
        "json_depth": 16,
        "json_nodes": 50000,
        "json_string_characters": 8192,
    }
    assert set(contract["audited_published_multiqc_tables"]) == {
        "multiqc_general_stats.txt",
        "multiqc_star.txt",
        "multiqc_salmon.txt",
        "multiqc_cutadapt.txt",
        "multiqc_fastp.txt",
        "multiqc_fastqc_fastqc_raw.txt",
        "multiqc_fastqc_fastqc_trimmed.txt",
        "multiqc_fastqc_fastqc_filtered.txt",
        "multiqc_picard_dups.txt",
        "multiqc_featurecounts_biotype_plot.txt",
        "multiqc_rseqc_bam_stat.txt",
        "multiqc_rseqc_infer_experiment.txt",
        "multiqc_rseqc_read_distribution.txt",
        "multiqc_rseqc_tin.txt",
    }


def test_results_contract_artifact_family_membership_matches_public_catalog():
    contract = load_bulk_rnaseq_results_contract()
    membership = contract["artifact_family_output_types"]
    declared = set(contract["artifact_families"])

    assert set(membership) == declared
    assert all(values for values in membership.values())
    flattened = [
        output_type for values in membership.values() for output_type in values
    ]
    assert len(flattened) == len(set(flattened))
    for required_prefix in (
        "bulk_rnaseq.fastq.merged.",
        "bulk_rnaseq.star.unaligned.",
        "bulk_rnaseq.samtools.",
        "bulk_rnaseq.umi.",
    ):
        assert any(value.startswith(required_prefix) for value in flattened)


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
