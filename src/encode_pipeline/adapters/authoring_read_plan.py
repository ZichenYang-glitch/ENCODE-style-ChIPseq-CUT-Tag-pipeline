"""The read-list format used by existing managed-input drafts."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
import re
from typing import Any

_INPUT_FILE_REVISION_ID = re.compile(r"^inpfr_[0-9a-f]{32}$")


def normalize_read_plan(
    plans: Sequence[Mapping[str, Any]],
) -> list[dict[str, Any]]:
    """Normalize the existing read coordinates without applying any assay rules."""

    normalized = []
    seen_revisions: set[str] = set()
    for plan in plans:
        if not isinstance(plan, Mapping):
            raise ValueError("workflow sample plan has an unexpected shape")
        reads = plan.get("reads")
        if not isinstance(reads, list) or not 1 <= len(reads) <= 2:
            raise ValueError("workflow sample reads must be explicit")
        by_number = {}
        for read in reads:
            if not isinstance(read, Mapping) or set(read) != {
                "read",
                "input_file_revision_id",
            }:
                raise ValueError("workflow sample read has an unexpected shape")
            number = read["read"]
            revision = read["input_file_revision_id"]
            if (
                number not in (1, 2)
                or not isinstance(revision, str)
                or _INPUT_FILE_REVISION_ID.fullmatch(revision) is None
                or number in by_number
                or revision in seen_revisions
            ):
                raise ValueError("workflow sample read is missing or repeated")
            by_number[number] = dict(read)
            seen_revisions.add(revision)
        normalized.append({**plan, "reads": [by_number[n] for n in sorted(by_number)]})
    return normalized
