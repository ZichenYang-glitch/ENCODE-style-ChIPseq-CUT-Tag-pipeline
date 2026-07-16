"""Immutable identity loader for the pinned bulk RNA-seq results contract."""

from __future__ import annotations

from copy import deepcopy
from hashlib import sha256
from importlib import resources
import json
from typing import Any, Mapping


RESULTS_CONTRACT_FILE = "results-contract-3.26.0.json"
RESULTS_CONTRACT_SIZE = 10_586
RESULTS_CONTRACT_SHA256 = (
    "2e2a64389b66a6fcfb59d83281562b8973bbf3846542a4f20ddaa3bda3e0fe9f"
)
_CONTRACT_PACKAGE = "encode_pipeline.contracts.nfcore_rnaseq"
_DEFAULT_RSEQC_MODULES = (
    "bam_stat",
    "inner_distance",
    "infer_experiment",
    "junction_annotation",
    "junction_saturation",
    "read_distribution",
    "read_duplication",
)
_CSI_INCOMPATIBLE_RSEQC_MODULES = frozenset(
    {"inner_distance", "read_distribution", "tin"}
)


def effective_downstream_layout(
    original_layout: str,
    params: Mapping[str, object],
) -> str:
    """Return the fixed post-UMI-extraction layout used by downstream tools."""
    if original_layout not in {"SE", "PE"}:
        raise ValueError("invalid bulk RNA-seq sample layout")
    if (
        original_layout == "PE"
        and params.get("with_umi") is True
        and params.get("skip_umi_extract") is False
        and type(params.get("umi_discard_read")) is int
        and params.get("umi_discard_read") in {1, 2}
    ):
        return "SE"
    return original_layout


def effective_rseqc_modules(params: Mapping[str, object]) -> tuple[str, ...]:
    """Mirror the pinned nf-core RSeQC module removal for CSI-indexed BAMs."""
    configured = params.get("rseqc_modules")
    modules = (
        tuple(str(configured).split(","))
        if configured is not None
        else _DEFAULT_RSEQC_MODULES
    )
    if params.get("bam_csi_index") is True:
        modules = tuple(
            module
            for module in modules
            if module not in _CSI_INCOMPATIBLE_RSEQC_MODULES
        )
    return modules


def load_bulk_rnaseq_results_contract() -> dict[str, Any]:
    """Return a fresh verified copy of the exact results contract document."""
    content = (
        resources.files(_CONTRACT_PACKAGE).joinpath(RESULTS_CONTRACT_FILE).read_bytes()
    )
    if (
        len(content) != RESULTS_CONTRACT_SIZE
        or sha256(content).hexdigest() != RESULTS_CONTRACT_SHA256
    ):
        raise ValueError("bulk RNA-seq results contract identity mismatch")
    try:
        value = json.loads(content)
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ValueError("bulk RNA-seq results contract is invalid") from exc
    if (
        not isinstance(value, dict)
        or value.get("schema_version") != "1.0.0"
        or value.get("workflow_id") != "bulk-rnaseq"
        or not isinstance(value.get("upstream"), dict)
        or value["upstream"].get("release") != "3.26.0"
        or value["upstream"].get("commit") != "e7ca46272c8f9d5ceee3f71759f4ba551d3217a4"
        or value["upstream"].get("route") != "star_salmon"
    ):
        raise ValueError("bulk RNA-seq results contract is invalid")
    return deepcopy(value)
