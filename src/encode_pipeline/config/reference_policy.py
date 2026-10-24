"""Pure reference-policy resolution shared by workflow and Python entry points.

This module intentionally has no filesystem or validator imports. It is safe
for the samples layer to import and is consumed by config validation, the
manifest/adapter path, and the Snakemake workflow.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence

from encode_pipeline.config import defaults
from encode_pipeline.errors import ValidationError


def _strict_bool(value, *, name: str) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        text = value.strip().lower()
        if text in ("true", "yes", "1"):
            return True
        if text in ("false", "no", "0"):
            return False
    raise ValidationError(f"{name} must be true or false, got {value!r}")


def _choice(value, *, name: str, allowed: Sequence[str]) -> str:
    if value not in allowed:
        raise ValidationError(f"{name} must be one of {sorted(allowed)}, got {value!r}")
    return str(value)


def resolve_blacklist_mode(config: Mapping, qc_raw: Mapping | None = None) -> str:
    """Resolve canonical/legacy blacklist mode without filesystem checks."""
    if qc_raw is None:
        qc_raw = config.get("qc")

    raw_mode = config.get("blacklist_mode")
    if raw_mode is False or (
        isinstance(raw_mode, str) and raw_mode.lower() in ("false", "no", "off")
    ):
        raw_mode = "off"
    elif raw_mode is True or (
        isinstance(raw_mode, str) and raw_mode.lower() in ("true", "yes")
    ):
        raw_mode = "copies"
    alias_present = isinstance(qc_raw, Mapping) and "blacklist_filter" in qc_raw
    alias_value = None
    if alias_present:
        alias_value = _strict_bool(
            qc_raw["blacklist_filter"], name="qc.blacklist_filter"
        )
        mapped = "copies" if alias_value else "off"
    else:
        mapped = None

    if raw_mode is None:
        return mapped if mapped is not None else "analysis"

    mode = _choice(
        raw_mode,
        name="blacklist_mode",
        allowed=defaults.BLACKLIST_MODES,
    )
    if mapped is not None:
        if mode == "analysis":
            if alias_value is True:
                raise ValidationError(
                    "config: blacklist_mode=analysis conflicts with "
                    "qc.blacklist_filter=true (copies)."
                )
        elif mode != mapped:
            raise ValidationError(
                "config: blacklist_mode conflicts with qc.blacklist_filter. "
                "Use one canonical key or make them consistent."
            )
    return mode


def assay_organelle_mode(assay: str) -> str:
    return "exclude" if assay in ("atac", "cuttag") else "report"


def control_context_assays(
    sample_id: str,
    sample_map: Mapping[str, Mapping],
    treatment_ids: Sequence[str],
) -> list[str]:
    return sorted(
        {
            str(sample_map[sid].get("assay", ""))
            for sid in treatment_ids
            if sample_map[sid].get("control_sample") == sample_id
        }
    )


def resolve_organelle_mode(
    sample: Mapping,
    *,
    config_mode: str,
    sample_map: Mapping[str, Mapping],
    treatment_ids: Sequence[str],
    error_cls=ValueError,
) -> str:
    """Resolve auto/explicit organelle mode for a treatment or control.

    Controls inherit the effective policy of the treatments that reference
    them. Assays that resolve to the same policy are interchangeable; different
    resolved policies fail closed.
    """
    if config_mode != "auto":
        return str(config_mode)
    if sample.get("role") != "control":
        return assay_organelle_mode(str(sample.get("assay", "")))

    sample_id = str(sample.get("id", ""))
    contexts = control_context_assays(sample_id, sample_map, treatment_ids)
    if not contexts:
        return assay_organelle_mode(str(sample.get("assay", "")))
    modes = sorted({assay_organelle_mode(assay) for assay in contexts})
    if len(modes) > 1:
        raise error_cls(
            f"Control {sample_id!r}: referenced treatments resolve to "
            f"different organelle policies: {modes}"
        )
    return modes[0]
