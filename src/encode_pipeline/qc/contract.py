"""Versioned QC output contracts for the ENCODE-style workflow adapter.

The workflow writers and the adapter reader share these exact headers so a
new writer cannot silently drift from the parser.  Version 2 is additive:
the version column comes first, the complete version 1 column order is
preserved, and new fields follow.  Version 1 documents remain readable and
are never rewritten during extraction.
"""

from __future__ import annotations

SAMPLE_QC_SCHEMA_VERSION = "2"
MNASE_QC_SCHEMA_VERSION = "2"
FRAGMENT_QC_SCHEMA_VERSION = "2"
MITO_QC_SCHEMA_VERSION = "2"
ANALYSIS_POLICY_SCHEMA_VERSION = "1"

SAMPLE_QC_V1_COLUMNS = (
    "sample",
    "assay",
    "target",
    "genome",
    "layout",
    "peak_mode",
    "use_control",
    "control_type",
    "final_bam",
    "peaks",
    "blacklist",
    "blacklist_filtered_bam",
    "blacklist_filtered_peaks",
    "total_reads",
    "reads_in_peaks",
    "frip",
    "peak_count",
    "blacklist_filtered_peak_count",
    "metrics_source",
    "unpaired_reads_examined",
    "read_pairs_examined",
    "secondary_or_supplementary_reads",
    "unmapped_reads",
    "unpaired_read_duplicates",
    "read_pair_duplicates",
    "read_pair_optical_duplicates",
    "percent_duplication",
    "estimated_library_size",
    "total_reads_examined",
    "duplicate_reads_estimate",
    "total_fragments",
    "distinct_fragments",
    "one_read_fragments",
    "two_read_fragments",
    "nrf",
    "pbc1",
    "pbc2",
)

SAMPLE_QC_V2_NEW_COLUMNS = (
    "analysis_bam",
    "analysis_input_stage",
    "analysis_valid_fragments",
    "analysis_invalid_templates",
    "analysis_tlen_zero_templates",
    "analysis_orphan_templates",
    "analysis_cross_contig_templates",
    "analysis_inconsistent_templates",
    "analysis_upstream_dropped_templates",
    "analysis_upstream_dropped_alignments",
    "raw_peak_count",
    "delivered_peak_count",
    "duplicate_analysis",
    "resolved_max_insert",
    "blacklist_mode",
    "blacklist_applied",
    "organelle_mode",
    "organelle_status",
    "organelle_reads",
    "organelle_total_mapped_reads",
    "organelle_fraction",
    "matched_organelle_contigs",
    "cross_correlation_nsc",
    "cross_correlation_rsc",
    "cross_correlation_status",
    "scientific_policy_revision",
)

SAMPLE_QC_V2_COLUMNS = (
    ("schema_version",) + SAMPLE_QC_V1_COLUMNS + SAMPLE_QC_V2_NEW_COLUMNS
)

MNASE_QC_V1_COLUMNS = (
    "sample",
    "assay",
    "peak_mode",
    "sub_min",
    "sub_max",
    "mono_min",
    "mono_max",
    "di_min",
    "di_max",
    "dyad_min",
    "dyad_max",
    "sub_bam",
    "mono_bam",
    "di_bam",
    "dyad_bigwig",
    "mono_bigwig",
    "insert_size_metrics",
    "caller_danpos3_enabled",
    "caller_inps_enabled",
    "caller_sem_enabled",
    "sub_reads",
    "mono_reads",
    "di_reads",
)

MNASE_QC_V2_NEW_COLUMNS = (
    "analysis_bam",
    "analysis_input_stage",
    "analysis_status",
    "analysis_total_templates",
    "analysis_valid_templates",
    "analysis_invalid_templates",
    "analysis_tlen_zero_templates",
    "analysis_orphan_templates",
    "analysis_cross_contig_templates",
    "analysis_inconsistent_templates",
    "analysis_classified_templates",
    "analysis_unclassified_templates",
    "analysis_fraction_classified",
    "window_overlap",
    "sub_fragments",
    "mono_fragments",
    "di_fragments",
    "upstream_dropped_templates",
    "upstream_dropped_alignments",
)

MNASE_QC_V2_COLUMNS = (
    ("schema_version",) + MNASE_QC_V1_COLUMNS + MNASE_QC_V2_NEW_COLUMNS
)

FRAGMENT_QC_V1_COLUMNS = (
    "sample",
    "layout",
    "fragment_count",
    "fragment_mean",
    "fragment_median",
    "fragment_mode",
    "fragment_min",
    "fragment_max",
    "fraction_lt_150",
    "fraction_150_300",
    "fraction_300_500",
    "fraction_ge_500",
    "fraction_lt_120",
    "status",
)

FRAGMENT_QC_V2_NEW_COLUMNS = (
    "input_stage",
    "duplicate_analysis",
    "blacklist_applied",
    "templates_total",
    "templates_valid",
    "templates_invalid",
    "invalid_tlen_zero",
    "invalid_orphan",
    "invalid_cross_contig",
    "invalid_inconsistent",
    "alignments_total",
    "upstream_dropped_templates",
    "upstream_dropped_alignments",
    "upstream_alignments_total",
    "upstream_alignments_kept",
)

FRAGMENT_QC_V2_COLUMNS = (
    ("schema_version",) + FRAGMENT_QC_V1_COLUMNS + FRAGMENT_QC_V2_NEW_COLUMNS
)

MITO_QC_V1_COLUMNS = (
    "sample",
    "organelle_reads",
    "total_mapped_reads",
    "fraction",
    "mode",
    "status",
)

MITO_QC_V2_COLUMNS = (
    ("schema_version",)
    + MITO_QC_V1_COLUMNS
    + (
        "input_stage",
        "matched_contigs",
        "organelle_chroms_explicit",
    )
)

POOLED_QC_COLUMNS = (
    "experiment",
    "assay",
    "target",
    "inferred_histone_class",
    "expected_peak_mode",
    "configured_peak_mode",
    "peak_mode_status",
    "biological_replicates",
    "biological_replicate_labels",
    "pooled_bam",
    "pooled_peaks",
    "pooled_peak_count",
    "pooled_FE_bdg",
    "pooled_ppois_bdg",
    "signal_tracks_status",
)

ANALYSIS_POLICY_KEYS = (
    "schema_version",
    "scientific_policy_revision",
    "sample",
    "assay",
    "profile",
    "resolved_duplicate_analysis",
    "resolved_max_insert",
    "organelle_mode",
    "organelle_status",
    "matched_organelle_contigs",
    "blacklist_mode",
    "blacklist_applied",
    "input_stages",
    "stage_counts",
    "reference_summary",
    "qc_schema_versions",
)

ANALYSIS_POLICY_INPUT_STAGE_KEYS = (
    "fragment_size",
    "mnase_classification",
    "duplicate_metrics",
    "nrf_pbc",
    "organelle",
    "blacklist",
    "peaks",
)

POLICY_STAGE_COUNT_KEYS = (
    "analysis_input_alignments",
    "analysis_kept_alignments",
    "analysis_valid_templates",
    "analysis_invalid_templates",
    "upstream_dropped_templates",
    "upstream_dropped_alignments",
    "blacklist_dropped_templates",
    "blacklist_dropped_alignments",
    "organelle_reads",
    "organelle_total_mapped_reads",
)

POLICY_ORGANELLE_STATUSES = frozenset(
    {
        "ok",
        "unresolved_contigs",
        "no_mapped_reads",
        "not_applicable",
        "missing",
        "parse_failed",
    }
)
POLICY_CROSS_CORRELATION_STATUSES = frozenset(
    {"reported", "not_enabled", "missing", "parse_failed"}
)
POLICY_DUPLICATE_POLICIES = frozenset({"exclude", "retain"})
POLICY_BLACKLIST_MODES = frozenset({"off", "copies", "analysis"})
POLICY_ORGANELLE_MODES = frozenset({"auto", "exclude", "report", "keep"})
POLICY_MAX_INSERT_NOT_APPLICABLE = "not_applicable"


def is_positive_integer_text(value: str) -> bool:
    """Return True for a canonical positive decimal integer string."""
    return isinstance(value, str) and value.isdigit() and int(value) > 0
