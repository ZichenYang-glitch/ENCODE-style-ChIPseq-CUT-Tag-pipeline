"""Build ENCODE sample rows from an explicit authoring plan.

The browser supplies revision ids and business fields. This module assigns
the server-resolved read paths. It does not infer role, pairing, or controls
from file names.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
import re
from typing import Any

from encode_pipeline.config.defaults import SAMPLE_ID_RE
from encode_pipeline.platform.authoring import WORKFLOW_SAMPLE_PLAN_INVALID
from encode_pipeline.platform.results import Issue, Result


_SAMPLE_REVISION_ID = re.compile(r"^smpr_[0-9a-f]{32}$")
_SAMPLE_TOKEN = re.compile(r"^[A-Za-z0-9_.-]{1,64}$")
_FIELDS = {
    "assay",
    "control_sample_revision_id",
    "layout",
    "peak_mode",
    "reads",
    "role",
    "sample_revision_id",
    "target",
}
_ASSAYS = {"chipseq", "cuttag", "atac", "mnase"}
_LAYOUTS = {"PE", "SE"}
_ROLES = {"treatment", "control"}
_PEAK_MODES = {"narrow", "broad", "nucleosome"}


def materialize_encode_sample_rows(
    plans: Sequence[Mapping[str, Any]],
    *,
    sample_names: Mapping[str, str],
    read_paths: Mapping[str, str],
) -> Result[list[dict[str, str]]]:
    """Return workflow rows for one already-normalized plan."""

    if not plans:
        return _invalid()
    rows: list[dict[str, str]] = []
    names: dict[str, str] = {}
    for plan in plans:
        revision_id = str(plan["sample_revision_id"])
        sample_name = sample_names.get(revision_id)
        if sample_name is None or SAMPLE_ID_RE.fullmatch(sample_name) is None:
            return _invalid()
        if sample_name in names.values():
            return _invalid()
        names[revision_id] = sample_name
        reads = {
            int(item["read"]): str(item["input_file_revision_id"])
            for item in plan["reads"]
        }
        fastq_1 = read_paths.get(reads.get(1, ""))
        fastq_2 = ""
        if 2 in reads:
            fastq_2 = read_paths.get(reads[2], "")
        if not fastq_1 or (plan["layout"] == "PE" and not fastq_2):
            return _invalid()
        if not _private_absolute(fastq_1) or (
            fastq_2 and not _private_absolute(fastq_2)
        ):
            return _invalid()
        rows.append(
            {
                "assay": str(plan["assay"]),
                "control_sample": "",
                "fastq_1": fastq_1,
                "fastq_2": fastq_2,
                "layout": str(plan["layout"]),
                "peak_mode": str(plan["peak_mode"]),
                "role": str(plan["role"]),
                "sample": sample_name,
                "target": str(plan["target"]),
            }
        )
    by_revision = {
        str(plan["sample_revision_id"]): row
        for plan, row in zip(plans, rows, strict=True)
    }
    for plan, row in zip(plans, rows, strict=True):
        control = plan["control_sample_revision_id"]
        if control is None:
            continue
        control_row = by_revision.get(str(control))
        if control_row is None or control_row["role"] != "control":
            return _invalid()
        row["control_sample"] = control_row["sample"]
    return Result.success(rows)


def normalize_encode_workflow_sample_plans(
    plans: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    """Apply ENCODE sample shape and assay rules to an identity-checked plan."""

    if not isinstance(plans, Sequence) or isinstance(plans, (str, bytes)):
        raise ValueError("workflow_samples must be a bounded list")
    normalized: list[dict[str, Any]] = []
    for item in plans:
        required = _FIELDS - {"control_sample_revision_id"}
        if not isinstance(item, Mapping) or not required <= set(item) <= _FIELDS:
            raise ValueError("workflow sample plan has an unexpected shape")
        sample_revision_id = item["sample_revision_id"]
        assay = item["assay"]
        layout = item["layout"]
        role = item["role"]
        peak_mode = item["peak_mode"]
        target = item["target"]
        control = item.get("control_sample_revision_id")
        if (
            not isinstance(sample_revision_id, str)
            or _SAMPLE_REVISION_ID.fullmatch(sample_revision_id) is None
            or assay not in _ASSAYS
            or layout not in _LAYOUTS
            or role not in _ROLES
            or peak_mode not in _PEAK_MODES
            or not isinstance(target, str)
            or _SAMPLE_TOKEN.fullmatch(target) is None
        ):
            raise ValueError("workflow sample plan is missing an explicit field")
        if (assay == "atac" and peak_mode != "narrow") or (
            assay == "mnase" and (layout != "PE" or peak_mode != "nucleosome")
        ):
            raise ValueError("workflow sample assay and layout do not agree")
        if assay != "mnase" and peak_mode == "nucleosome":
            raise ValueError("workflow sample assay and peak mode do not agree")
        if role == "control" and control is not None:
            raise ValueError("a control sample cannot name another control")
        if control is not None and (
            not isinstance(control, str)
            or _SAMPLE_REVISION_ID.fullmatch(control) is None
            or control == sample_revision_id
        ):
            raise ValueError("control sample revision is not explicit")
        reads = item["reads"]
        if not isinstance(reads, list):
            raise ValueError("workflow sample reads must be explicit")
        expected = (1, 2) if layout == "PE" else (1,)
        numbers = tuple(sorted(int(read["read"]) for read in reads))
        if numbers != expected:
            raise ValueError("workflow sample reads do not match the layout")
        normalized.append(
            {
                "assay": assay,
                "control_sample_revision_id": control,
                "layout": layout,
                "peak_mode": peak_mode,
                "reads": [dict(read) for read in reads],
                "role": role,
                "sample_revision_id": sample_revision_id,
                "target": target,
            }
        )
    controls = {
        item["sample_revision_id"] for item in normalized if item["role"] == "control"
    }
    for item in normalized:
        control = item["control_sample_revision_id"]
        if control is not None and control not in controls:
            raise ValueError("control sample revision is not a control row")
    return normalized


def _private_absolute(value: str) -> bool:
    return value.startswith("/") and "\n" not in value and "\t" not in value


def _invalid() -> Result[list[dict[str, str]]]:
    return Result.failure(
        [
            Issue(
                code=WORKFLOW_SAMPLE_PLAN_INVALID,
                message=(
                    "The workflow sample plan does not match the selected revisions."
                ),
                path="workflow_samples",
                source="authoring",
            )
        ]
    )
