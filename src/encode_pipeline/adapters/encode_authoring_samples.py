"""Build ENCODE sample rows from an explicit authoring plan.

The browser supplies revision ids and business fields. This module assigns
the server-resolved read paths. It does not infer role, pairing, or controls
from file names.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

from encode_pipeline.config.defaults import SAMPLE_ID_RE
from encode_pipeline.platform.authoring import WORKFLOW_SAMPLE_PLAN_INVALID
from encode_pipeline.platform.results import Issue, Result


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
