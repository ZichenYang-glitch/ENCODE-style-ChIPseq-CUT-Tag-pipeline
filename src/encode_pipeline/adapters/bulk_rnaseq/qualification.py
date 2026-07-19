"""Server-owned execution modes for bulk RNA-seq qualification runs."""

from __future__ import annotations

from collections.abc import Mapping
from enum import Enum
from typing import Any


class BulkRnaSeqExecutionMode(Enum):
    """Closed deployment policy; workflow inputs cannot select these modes."""

    STANDARD = "standard-v1"
    RAPID_QUANT = "rapid-quant-v1"


# nf-core/rnaseq 3.26.0's immutable ``conf/rapid_quant.config`` owns these
# values. Omitting every key from params.json is necessary because params-file
# values would otherwise override the selected profile.
RAPID_QUANT_PROFILE_OWNED_PARAMETERS = frozenset(
    {
        "pseudo_aligner",
        "skip_alignment",
        "skip_bigwig",
        "skip_biotype_qc",
        "skip_deseq2_qc",
        "skip_dupradar",
        "skip_markduplicates",
        "skip_multiqc",
        "skip_preseq",
        "skip_qualimap",
        "skip_quantification_merge",
        "skip_rseqc",
        "skip_stringtie",
    }
)

# The pinned profile relies on the base-config default ``false`` for this
# upstream branch guard. The authoring snapshot deliberately sets it to true
# for standard runs, so qualification mode must omit it alongside the values
# written directly by rapid_quant.config or Salmon would never be scheduled.
RAPID_QUANT_MODE_OWNED_PARAMETERS = RAPID_QUANT_PROFILE_OWNED_PARAMETERS | {
    "skip_pseudo_alignment"
}

_NEXTFLOW_PROFILE_BY_MODE = {
    BulkRnaSeqExecutionMode.STANDARD: None,
    BulkRnaSeqExecutionMode.RAPID_QUANT: "rapid_quant",
}


def execution_mode_identity(mode: BulkRnaSeqExecutionMode) -> str:
    """Return the stable identity coordinate for one server-owned mode."""
    if not isinstance(mode, BulkRnaSeqExecutionMode):
        raise ValueError("bulk RNA-seq execution mode is invalid")
    return mode.value


def nextflow_profile_for_mode(mode: BulkRnaSeqExecutionMode) -> str | None:
    """Return the exact scientific profile, if this mode requires one.

    Docker is enabled only by the platform hard config. Selecting nf-core's
    ``docker`` profile would replace the platform-owned run options after the
    hard config is evaluated.
    """
    try:
        return _NEXTFLOW_PROFILE_BY_MODE[mode]
    except (KeyError, TypeError) as error:
        raise ValueError("bulk RNA-seq execution mode is invalid") from error


def runtime_params_for_mode(
    params: Mapping[str, Any],
    *,
    mode: BulkRnaSeqExecutionMode,
) -> dict[str, Any]:
    """Remove only immutable-profile-owned params for qualification mode."""
    execution_mode_identity(mode)
    excluded = (
        RAPID_QUANT_MODE_OWNED_PARAMETERS
        if mode is BulkRnaSeqExecutionMode.RAPID_QUANT
        else frozenset()
    )
    return {name: value for name, value in params.items() if name not in excluded}
