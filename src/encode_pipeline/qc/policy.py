"""Strict reader for the versioned per-sample analysis policy record."""

from __future__ import annotations

import json
import re

from encode_pipeline.qc.contract import (
    ANALYSIS_POLICY_INPUT_STAGE_KEYS,
    ANALYSIS_POLICY_KEYS,
    ANALYSIS_POLICY_V1_KEYS,
    ANALYSIS_POLICY_V2_INPUT_STAGE_KEYS,
    POLICY_ANALYSIS_TEMPLATE_STATUSES,
    POLICY_BLACKLIST_MODES,
    POLICY_DUPLICATE_POLICIES,
    POLICY_MAX_INSERT_NOT_APPLICABLE,
    POLICY_ORGANELLE_MODES,
    POLICY_ORGANELLE_STATUSES,
    POLICY_STAGE_COUNT_KEYS,
)

ANALYSIS_POLICY_VERSION = 2
_ALLOWED_INPUT_STAGES = frozenset(
    {
        "analysis_bam",
        "pre_duplicate_markdup_bam",
        "pre_duplicate_mapq_bam",
        "final_bam_to_analysis_bam",
        "not_applicable",
    }
)
_SAFE_TOKEN = re.compile(r"^[^\s/\\\x00]+$")
_ABSOLUTE_PATH = re.compile(r"^(?:/|[A-Za-z]:[\\/])")


class AnalysisPolicyError(ValueError):
    """The policy document does not satisfy the fixed versioned contract."""


def _reject_duplicate_keys(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise AnalysisPolicyError("analysis policy has duplicate keys")
        result[key] = value
    return result


def _require_exact_keys(value, expected, label):
    if not isinstance(value, dict) or set(value) != set(expected):
        raise AnalysisPolicyError(f"{label} keys are invalid")


def _require_str(value, label, *, token=False):
    if not isinstance(value, str) or not value:
        raise AnalysisPolicyError(f"{label} is invalid")
    if _ABSOLUTE_PATH.match(value) or "\x00" in value:
        raise AnalysisPolicyError(f"{label} leaks a path")
    if token and _SAFE_TOKEN.fullmatch(value) is None:
        raise AnalysisPolicyError(f"{label} is invalid")
    return value


def _optional_count(value, label):
    if value is None:
        return None
    if type(value) is not int or value < 0:
        raise AnalysisPolicyError(f"{label} is invalid")
    return value


def _validate_analysis_template_status(value, counts, status):
    """Validate the version 2 reason behind analysis-template counts."""
    if status not in POLICY_ANALYSIS_TEMPLATE_STATUSES:
        raise AnalysisPolicyError("analysis template status is invalid")
    valid = counts["analysis_valid_templates"]
    invalid = counts["analysis_invalid_templates"]
    if status in ("not_applicable", "not_collected", "parse_failed"):
        if valid is not None or invalid is not None:
            raise AnalysisPolicyError(
                "analysis template counts must be NA for the stated status"
            )
    else:
        if valid is None or invalid is None:
            raise AnalysisPolicyError("analysis template counts are missing")
        if status == "ok" and valid == 0:
            raise AnalysisPolicyError("analysis template ok status has no valid count")
        if status == "no_valid_templates" and valid != 0:
            raise AnalysisPolicyError(
                "analysis template zero status has a nonzero valid count"
            )
    stage = value["input_stages"]["analysis_template_metrics"]
    if status in ("ok", "no_valid_templates") and stage != "analysis_bam":
        raise AnalysisPolicyError("analysis template input stage is inconsistent")
    if status in ("not_applicable", "not_collected") and stage != "not_applicable":
        raise AnalysisPolicyError("analysis template input stage is inconsistent")
    if status == "parse_failed" and stage not in ("analysis_bam", "not_applicable"):
        raise AnalysisPolicyError("analysis template input stage is inconsistent")


def parse_analysis_policy_bytes(content: bytes) -> dict:
    """Parse and validate one exact analysis-policy JSON document."""
    try:
        text = content.decode("utf-8", errors="strict")
    except UnicodeDecodeError as exc:
        raise AnalysisPolicyError("analysis policy is not UTF-8") from exc
    try:
        value = json.loads(text, object_pairs_hook=_reject_duplicate_keys)
    except (json.JSONDecodeError, TypeError) as exc:
        raise AnalysisPolicyError("analysis policy is not valid JSON") from exc
    if not isinstance(value, dict):
        raise AnalysisPolicyError("analysis policy must be an object")
    schema_version = value.get("schema_version")
    if schema_version == 1:
        expected_keys = ANALYSIS_POLICY_V1_KEYS
        input_stage_keys = ANALYSIS_POLICY_INPUT_STAGE_KEYS
        template_status = None
    elif schema_version == ANALYSIS_POLICY_VERSION:
        expected_keys = ANALYSIS_POLICY_KEYS
        input_stage_keys = ANALYSIS_POLICY_V2_INPUT_STAGE_KEYS
        template_status = value.get("analysis_template_status")
    else:
        raise AnalysisPolicyError("analysis policy version is unsupported")
    _require_exact_keys(value, expected_keys, "analysis policy")
    revision = value["scientific_policy_revision"]
    if type(revision) is not int or revision <= 0:
        raise AnalysisPolicyError("scientific policy revision is invalid")
    _require_str(value["sample"], "sample", token=True)
    _require_str(value["assay"], "assay", token=True)
    if value["profile"] not in ("legacy", "updated"):
        raise AnalysisPolicyError("analysis profile is invalid")
    if value["resolved_duplicate_analysis"] not in POLICY_DUPLICATE_POLICIES:
        raise AnalysisPolicyError("duplicate policy is invalid")
    max_insert = value["resolved_max_insert"]
    if max_insert != POLICY_MAX_INSERT_NOT_APPLICABLE and (
        type(max_insert) is not int or max_insert <= 0
    ):
        raise AnalysisPolicyError("resolved max insert is invalid")
    if value["organelle_mode"] not in POLICY_ORGANELLE_MODES:
        raise AnalysisPolicyError("organelle mode is invalid")
    if value["organelle_status"] not in POLICY_ORGANELLE_STATUSES:
        raise AnalysisPolicyError("organelle status is invalid")
    matched = value["matched_organelle_contigs"]
    if not isinstance(matched, list):
        raise AnalysisPolicyError("matched organelle contigs are invalid")
    for contig in matched:
        _require_str(contig, "matched organelle contig", token=True)
    if len(set(matched)) != len(matched):
        raise AnalysisPolicyError("matched organelle contigs are duplicated")
    if value["blacklist_mode"] not in POLICY_BLACKLIST_MODES:
        raise AnalysisPolicyError("blacklist mode is invalid")
    if type(value["blacklist_applied"]) is not bool:
        raise AnalysisPolicyError("blacklist applied flag is invalid")

    _require_exact_keys(value["input_stages"], input_stage_keys, "input stages")
    for key, stage in value["input_stages"].items():
        if stage not in _ALLOWED_INPUT_STAGES:
            raise AnalysisPolicyError(f"input stage {key} is invalid")

    _require_exact_keys(value["stage_counts"], POLICY_STAGE_COUNT_KEYS, "stage counts")
    counts = {
        key: _optional_count(item, f"stage count {key}")
        for key, item in value["stage_counts"].items()
    }
    if value["organelle_status"] != "ok":
        if counts["organelle_reads"] is not None:
            raise AnalysisPolicyError("organelle reads must be NA when not ok")
    if counts["organelle_total_mapped_reads"] is None:
        if counts["organelle_reads"] is not None:
            raise AnalysisPolicyError("organelle counts are inconsistent")

    if template_status is not None:
        _validate_analysis_template_status(value, counts, template_status)

    _require_exact_keys(
        value["reference_summary"],
        ("genome", "organelle_chroms_explicit"),
        "reference summary",
    )
    _require_str(value["reference_summary"]["genome"], "reference genome", token=True)
    if type(value["reference_summary"]["organelle_chroms_explicit"]) is not bool:
        raise AnalysisPolicyError("organelle explicit flag is invalid")

    _require_exact_keys(
        value["qc_schema_versions"],
        ("sample", "mnase", "fragment", "mito"),
        "QC schema versions",
    )
    if any(version != "2" for version in value["qc_schema_versions"].values()):
        raise AnalysisPolicyError("QC schema version is invalid")
    return value
