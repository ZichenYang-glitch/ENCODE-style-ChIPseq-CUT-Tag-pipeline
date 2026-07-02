from __future__ import annotations

import re
from dataclasses import dataclass


_SAFE_MESSAGE = (
    "The assistant response was filtered because it contained execution-like "
    "instructions. Review the validation issues directly or ask for a "
    "non-executing explanation."
)

_COMMAND_PREFIX = r"(?:^\s*|[\.\n:]\s+|\$\s*)"


@dataclass(frozen=True)
class FilterResult:
    text: str
    filtered: bool = False
    matched_patterns: tuple[str, ...] = ()


class OutputFilter:
    _PATTERNS: tuple[tuple[str, re.Pattern[str]], ...] = (
        ("snakemake_command", re.compile(_COMMAND_PREFIX + r"snakemake\s+", re.IGNORECASE)),
        ("rm_command", re.compile(_COMMAND_PREFIX + r"rm\s+(?:-[a-zA-Z]*f|--force)", re.IGNORECASE)),
        ("python_command", re.compile(_COMMAND_PREFIX + r"python3?(?:\.\d+)?\s+\S+\.py", re.IGNORECASE)),
        ("bash_command", re.compile(_COMMAND_PREFIX + r"bash\s+", re.IGNORECASE)),
        ("sbatch_command", re.compile(_COMMAND_PREFIX + r"sbatch\s+", re.IGNORECASE)),
        ("qsub_command", re.compile(_COMMAND_PREFIX + r"qsub\s+", re.IGNORECASE)),
        ("execution_directive", re.compile(r"\b(?:run|execute)\s+(?:this\s+)?(?:shell\s+)?command\b", re.IGNORECASE)),
        ("delete_directive", re.compile(r"\bdelete\s+this\s+file\b", re.IGNORECASE)),
    )

    def filter_text(self, text: str) -> FilterResult:
        matched = [
            name for name, pattern in self._PATTERNS
            if pattern.search(text)
        ]
        if matched:
            return FilterResult(
                text=_SAFE_MESSAGE,
                filtered=True,
                matched_patterns=tuple(matched),
            )
        return FilterResult(text=text)
