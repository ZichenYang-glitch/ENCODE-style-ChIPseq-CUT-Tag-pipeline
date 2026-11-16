"""Navigation scaffolds built into WashU, separate from reference identity."""

from __future__ import annotations


# Chromosome lists compiled into wuepgg 62.0.2. An assembly label selects one
# of those lists for navigation only. It does not prove the run's reference
# files are the same build WashU published.
_SCAFFOLDS = {
    "hg19": "hg19",
    "grch37": "hg19",
    "hg38": "hg38",
    "grch38": "hg38",
    "mm9": "mm9",
    "mm10": "mm10",
    "grcm38": "mm10",
    "mm39": "mm39",
    "grcm39": "mm39",
}


def viewer_scaffold(assembly_label: str | None) -> str | None:
    """Return a built-in WashU genome name for navigation, when one exists."""

    if not isinstance(assembly_label, str):
        return None
    return _SCAFFOLDS.get(assembly_label.strip().lower())
