"""Histone target classification for ChIP-seq QC.

Pure stdlib module. No file I/O, no Snakemake dependency. Used by the
Stage 6b pooled experiment QC summary rule to classify targets and
report peak_mode compatibility.

Usage:
    from histone_utils import classify_histone_target
    info = classify_histone_target("H3K27me3", "broad")
"""

# Target groups (uppercase, de-punctuated for matching)
BROAD_LIKE = {
    "H3K27ME3", "H3K36ME3", "H3K9ME3",
    "H3K79ME2", "H3K79ME3",
    "H4K20ME1", "H4K20ME3",
}

NARROW_LIKE = {
    "H3K4ME3",
    "H3K9AC", "H3K14AC", "H3K18AC", "H3K23AC", "H3K56AC",
    "H4K5AC", "H4K8AC", "H4K12AC", "H4K16AC",
    # H2A.Z variants
    "H2AZ", "H2AFZ",
}

CONTEXT_DEPENDENT = {
    "H3K27AC", "H3K4ME1", "H3K4ME2",
}


def _normalize(target):
    """Normalize a target string for matching: uppercase, strip, de-punctuate."""
    return target.strip().upper().replace("-", "").replace("_", "").replace(".", "")


def classify_histone_target(target, configured_peak_mode=None):
    """Classify a histone target for QC purposes.

    Args:
        target: the target name from the sample sheet (e.g. "H3K27me3")
        configured_peak_mode: the peak_mode from the sample sheet
            ("broad" or "narrow"), or None

    Returns a dict with keys:
        target_normalized       — uppercase, de-punctuated target
        inferred_histone_class  — "broad_like" | "narrow_like" |
                                   "context_dependent" | "unknown"
        expected_peak_mode      — "broad" | "narrow" | "broad_or_narrow" |
                                   "unknown"
        compatible_peak_modes   — ["broad"], ["narrow"], ["broad","narrow"], []
        peak_mode_status        — "ok" | "mismatch" | "unknown"

    Does NOT raise exceptions. Mismatches are reported via status only.
    """
    norm = _normalize(target)

    if norm in BROAD_LIKE:
        cls = "broad_like"
        expected = "broad"
        compatible = ["broad"]
    elif norm in NARROW_LIKE:
        cls = "narrow_like"
        expected = "narrow"
        compatible = ["narrow"]
    elif norm in CONTEXT_DEPENDENT:
        cls = "context_dependent"
        expected = "broad_or_narrow"
        compatible = ["broad", "narrow"]
    else:
        cls = "unknown"
        expected = "unknown"
        compatible = []

    # peak_mode_status
    if configured_peak_mode is None:
        status = "unknown"
    elif cls == "unknown" or expected == "unknown":
        status = "unknown"
    elif configured_peak_mode in compatible:
        status = "ok"
    else:
        status = "mismatch"

    return {
        "target_normalized": norm,
        "inferred_histone_class": cls,
        "expected_peak_mode": expected,
        "compatible_peak_modes": compatible,
        "peak_mode_status": status,
    }
