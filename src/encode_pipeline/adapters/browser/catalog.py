"""Map public artifact identity to one WashU track description."""

from __future__ import annotations

from dataclasses import dataclass


_HITRAC_WORKFLOW = "hitrac-preprocess"
_ENCODE_WORKFLOW = "encode-style-chipseq-cuttag-atac-mnase"
_PET_OUTPUTS = frozenset({"hitrac_bedpe_all", "hitrac_bedpe_no_bg"})
_SIGNAL_OUTPUTS = frozenset(
    {
        "cpm_bigwig",
        "fragment_cpm_bigwig",
        "mnase_fragment_cpm_bigwig",
        "mnase_dyad_bigwig",
        "mnase_mono_bigwig",
        "pooled_mnase_dyad_bigwig",
        "pooled_mnase_mono_bigwig",
        "macs3_fe_bw",
        "macs3_ppois_bw",
        "pooled_fe_bw",
        "pooled_ppois_bw",
    }
)
_INTERVAL_OUTPUTS = frozenset(
    {
        "macs3_peak",
        "pooled_macs3_peak",
        "pooled_filtered_peak",
        "chipseq_macs3_narrow_consensus_peak",
        "chipseq_macs3_broad_consensus_peak",
        "cuttag_macs3_narrow_consensus_peak",
        "cuttag_macs3_broad_consensus_peak",
        "atac_macs3_narrow_consensus_peak",
        "cuttag_seacr_consensus_peak",
        "cuttag_macs3_narrow_idr_final_peak",
        "atac_macs3_narrow_idr_final_peak",
        "chipseq_macs3_broad_idr_final_peak",
        "cuttag_macs3_broad_idr_final_peak",
        "chipseq_macs3_broad_consensus_final_peak",
        "cuttag_macs3_narrow_consensus_final_peak",
        "cuttag_macs3_broad_consensus_final_peak",
        "cuttag_seacr_consensus_final_peak",
    }
)
_ANNOTATION_OUTPUTS = frozenset({"gene_annotation_bed", "refbed"})


@dataclass(frozen=True)
class BrowserArtifactDescription:
    """Scientific reading of one artifact for the WashU view."""

    wash_u_type: str
    role: str
    kind: str
    source_form: str
    label: str
    detail: str
    score_meaning: str
    available: bool
    unavailable_reason: str | None = None


def describe_browser_artifact(
    *,
    workflow_id: str,
    output_type: str,
    name: str,
) -> BrowserArtifactDescription | None:
    """Return a track description, or None when the artifact is not a track."""

    if not isinstance(output_type, str) or not isinstance(name, str):
        return None
    lowered = name.lower()
    if output_type.endswith("_summary") or output_type == "hitrac_summary":
        return None
    if lowered.endswith((".cool", ".mcool")) or output_type in {"cool", "mcool"}:
        return _unavailable(
            "hic",
            "interaction",
            "contact_matrix",
            "Cooler matrices are not loaded. This view does not run HiGlass.",
            "unsupported_contact_matrix",
        )
    if lowered.endswith(".hic") or output_type in {"hic", "contact_matrix_hic"}:
        return _original(
            "hic",
            "interaction",
            "contact_matrix",
            "Contact matrix",
            "Values are matrix contacts from the .hic file.",
            "matrix_contact",
        )
    if output_type in _PET_OUTPUTS or (
        workflow_id == _HITRAC_WORKFLOW and lowered.endswith(".bedpe.gz")
    ):
        return _derivative(
            "longrange",
            "interaction",
            "pet_records",
            "PET records",
            "Each source row is one paired-end tag, not a significant loop "
            "and not a contact count. The display score is a flag on each "
            "anchor and must not be added.",
            "display_flag",
        )
    if lowered.endswith((".bedpe", ".bedpe.gz")):
        return _derivative(
            "longrange",
            "interaction",
            "paired_intervals",
            "Paired intervals",
            "These are paired genomic intervals. Source score and MAPQ "
            "columns are not treated as interaction strength.",
            "display_flag",
        )
    if output_type in _SIGNAL_OUTPUTS or lowered.endswith(
        (".bw", ".bigwig", ".bigbed")
    ):
        if lowered.endswith(".bigbed"):
            return _original(
                "bigbed",
                "interval",
                "interval",
                "Intervals",
                "Intervals from the original bigBed file.",
                "source_interval",
            )
        return _original(
            "bigwig",
            "signal",
            "signal",
            "Signal",
            "Numeric values are signal from the original bigWig, including zeros.",
            "signal_value",
        )
    if output_type in _ANNOTATION_OUTPUTS or lowered.endswith(
        (".refbed.gz", ".genes.bed.gz")
    ):
        return _derivative(
            "bed",
            "annotation",
            "annotation",
            "Annotation",
            "Intervals from an annotation artifact for this run.",
            "annotation_label",
        )
    if (
        workflow_id == _ENCODE_WORKFLOW and output_type in _INTERVAL_OUTPUTS
    ) or _is_interval_name(lowered):
        return _derivative(
            "bed",
            "interval",
            "peak_intervals",
            "Peaks",
            "Intervals from the peak or BED artifact. A copied score is that "
            "file's own score, not a contact count.",
            "source_interval_score",
        )
    return None


def _is_interval_name(name: str) -> bool:
    if ".bedpe" in name:
        return False
    return name.endswith(
        (
            ".bed",
            ".bed.gz",
            ".narrowpeak",
            ".narrowpeak.gz",
            ".broadpeak",
            ".broadpeak.gz",
        )
    )


def _original(
    wash_u_type: str,
    role: str,
    kind: str,
    label: str,
    detail: str,
    score_meaning: str,
) -> BrowserArtifactDescription:
    return BrowserArtifactDescription(
        wash_u_type=wash_u_type,
        role=role,
        kind=kind,
        source_form="original",
        label=label,
        detail=detail,
        score_meaning=score_meaning,
        available=True,
    )


def _derivative(
    wash_u_type: str,
    role: str,
    kind: str,
    label: str,
    detail: str,
    score_meaning: str,
) -> BrowserArtifactDescription:
    return BrowserArtifactDescription(
        wash_u_type=wash_u_type,
        role=role,
        kind=kind,
        source_form="display_derivative",
        label=label,
        detail=detail,
        score_meaning=score_meaning,
        available=True,
    )


def _unavailable(
    wash_u_type: str,
    role: str,
    kind: str,
    detail: str,
    reason: str,
) -> BrowserArtifactDescription:
    return BrowserArtifactDescription(
        wash_u_type=wash_u_type,
        role=role,
        kind=kind,
        source_form="original",
        label="Unsupported",
        detail=detail,
        score_meaning="not_applicable",
        available=False,
        unavailable_reason=reason,
    )
