"""MNase-seq config validation helpers."""

from encode_pipeline.config import defaults

__all__ = ["validate_mnase_config"]


def _parse_range_int(value, label, original, error_cls=ValueError):
    """Parse one range endpoint as an integer, rejecting bools and floats."""
    if isinstance(value, bool):
        raise error_cls(
            f"{label} elements must be integers, got {original!r}"
        )
    if isinstance(value, int):
        return value
    if isinstance(value, str):
        text = value.strip()
        if text.isdigit():
            return int(text)
    raise error_cls(
        f"{label} elements must be integers, got {original!r}"
    )


def _validate_range_pair(r, label, error_cls=ValueError):
    """Validate a fragment/dyad range pair [min, max].

    *r* may be a list or tuple. Must have exactly 2 positive ints,
    with min < max. Raises *error_cls* on failure.
    """
    if not isinstance(r, (list, tuple)):
        raise error_cls(
            f"{label} must be a list of 2 positive ints "
            f"(min < max), got {type(r).__name__}"
        )
    if len(r) != 2:
        raise error_cls(
            f"{label} must have exactly 2 elements "
            f"(min, max), got {len(r)}"
        )
    lo, hi = (
        _parse_range_int(r[0], label, r, error_cls=error_cls),
        _parse_range_int(r[1], label, r, error_cls=error_cls),
    )
    if lo <= 0 or hi <= 0:
        raise error_cls(
            f"{label} values must be positive, got [{lo}, {hi}]"
        )
    if lo >= hi:
        raise error_cls(
            f"{label}: min must be < max, got [{lo}, {hi}]"
        )


def validate_mnase_config(mnase: dict, error_cls=ValueError) -> dict:
    """Validate the mnase config block. Returns normalized dict.

    Absent block -> all defaults. Stage 39 keys + Stage 40 fragments,
    dyad_range, and callers.

    Fragment range precedence: fragments.mono > mono_range > [140, 200].
    """
    if not isinstance(mnase, dict):
        raise error_cls(
            f"mnase must be a mapping, got {type(mnase).__name__}"
        )

    known = defaults.MNASE_TOP_KEYS
    for key in mnase:
        if key not in known:
            raise error_cls(
                f"mnase: unknown key {key!r}. Known: {sorted(known)}"
            )

    # --- mono_range (Stage 39, deprecated in favor of fragments.mono) ---

    mono_range = mnase.get(
        "mono_range", defaults.MNASE_MONO_RANGE_DEFAULT
    )
    _validate_range_pair(mono_range, "mnase.mono_range", error_cls=error_cls)

    # --- fragments (Stage 40) ---

    fragments_raw = mnase.get("fragments", {})
    if not isinstance(fragments_raw, dict):
        raise error_cls(
            f"mnase.fragments must be a mapping, "
            f"got {type(fragments_raw).__name__}"
        )

    fragments_known = defaults.MNASE_FRAGMENTS_KEYS
    for key in fragments_raw:
        if key not in fragments_known:
            raise error_cls(
                f"mnase.fragments: unknown key {key!r}. "
                f"Known: {sorted(fragments_known)}"
            )

    def _resolve_fragment_range(class_name, legacy_mono, hard_default):
        """Resolve a single fragment range with the backward-compat chain.

        class_name   -> key in fragments_raw
        legacy_mono  -> value from mono_range (only used for "mono")
        hard_default -> [lo, hi] fallback
        """
        if class_name in fragments_raw:
            r = fragments_raw[class_name]
            _validate_range_pair(
                r, f"mnase.fragments.{class_name}", error_cls=error_cls
            )
            return [int(r[0]), int(r[1])]
        if class_name == "mono":
            return [int(legacy_mono[0]), int(legacy_mono[1])]
        return [int(hard_default[0]), int(hard_default[1])]

    fragments = {
        "sub": _resolve_fragment_range(
            "sub", mono_range, defaults.MNASE_FRAGMENT_DEFAULTS["sub"]
        ),
        "mono": _resolve_fragment_range(
            "mono", mono_range, defaults.MNASE_FRAGMENT_DEFAULTS["mono"]
        ),
        "di": _resolve_fragment_range(
            "di", mono_range, defaults.MNASE_FRAGMENT_DEFAULTS["di"]
        ),
    }

    # --- dyad_range (Stage 40) ---

    dyad_range = mnase.get("dyad_range", defaults.MNASE_DYAD_RANGE_DEFAULT)
    _validate_range_pair(dyad_range, "mnase.dyad_range", error_cls=error_cls)

    # --- callers (Stage 40, execution deferred) ---

    callers_raw = mnase.get("callers", {})
    if not isinstance(callers_raw, dict):
        raise error_cls(
            f"mnase.callers must be a mapping, "
            f"got {type(callers_raw).__name__}"
        )

    callers_known = defaults.MNASE_CALLERS
    callers = {}
    for ck in callers_raw:
        if ck not in callers_known:
            raise error_cls(
                f"mnase.callers: unknown key {ck!r}. "
                f"Known: {sorted(callers_known)}"
            )
        val = callers_raw[ck]
        if not isinstance(val, bool):
            raise error_cls(
                f"mnase.callers.{ck} must be boolean true or false, "
                f"got {val!r}"
            )
        enabled = val
        if enabled:
            raise error_cls(
                f"mnase.callers.{ck}: caller execution is not implemented "
                f"in v0.2; set {ck}: false"
            )
        callers[ck] = False

    # fill defaults for absent caller keys
    for ck in callers_known:
        callers.setdefault(ck, False)

    return {
        "mono_range": [int(mono_range[0]), int(mono_range[1])],
        "fragments": fragments,
        "dyad_range": [int(dyad_range[0]), int(dyad_range[1])],
        "callers": callers,
    }
