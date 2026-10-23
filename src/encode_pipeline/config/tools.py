"""Tool parameter validation helpers."""

import math
import shlex

from encode_pipeline.config import defaults

__all__ = ["validate_tool_params"]


_KEEP_DUP_CANONICAL = "--keep-dup"


def _is_keep_dup_option(name):
    """True when ``name`` selects MACS3's ``--keep-dup`` option.

    The pinned MACS3 uses argparse, which accepts any unambiguous prefix of
    a long option; ``--keep-dup`` is the only callpeak option starting with
    ``--k``, so every prefix from ``--k`` upwards (e.g. ``--keep-d``) is the
    same flag. Longer spellings such as ``--keep-duplicates`` are not
    prefixes of the canonical option and are left alone.
    """
    return len(name) >= 3 and _KEEP_DUP_CANONICAL.startswith(name)


def _migrate_keep_dup_extra_args(tool, extra, error_cls):
    """Normalize legacy ``--keep-dup`` tokens in MACS3 extra_args.

    Duplicate policy is executed at the BAM layer and the workflow renders
    ``--keep-dup all`` itself. A synonymous ``--keep-dup all`` (including
    argparse-accepted abbreviations and the ``=`` form) is removed here so
    the final command carries the flag exactly once; any other value
    conflicts with that contract and is rejected with migration guidance
    instead of being decided by argument order.

    When no managed flag is present the original string is returned
    byte-for-byte so quoting of unrelated arguments is preserved. When a
    migration does happen, the remaining tokens are re-quoted with
    shlex.join so shell argument boundaries (spaces, apostrophes) survive,
    and a second validation pass is a stable no-op.
    """
    try:
        tokens = shlex.split(extra)
    except ValueError as exc:
        raise error_cls(
            f"tool_parameters.{tool}.extra_args: invalid shell quoting "
            f"in {extra!r}: {exc}"
        ) from exc
    kept = []
    migrated = False
    index = 0
    while index < len(tokens):
        token = tokens[index]
        name, sep, attached = token.partition("=")
        if not _is_keep_dup_option(name):
            kept.append(token)
            index += 1
            continue
        migrated = True
        if sep:
            value = attached
            index += 1
        else:
            value = tokens[index + 1] if index + 1 < len(tokens) else ""
            index += 2
        if value != "all":
            raise error_cls(
                f"tool_parameters.{tool}.extra_args: {name} "
                f"{value or '<missing>'!r} conflicts with the BAM-layer "
                "duplicate policy. Remove --keep-dup from extra_args (the "
                "workflow always renders --keep-dup all) or use 'all'."
            )
    if not migrated:
        return extra
    return shlex.join(kept)


def validate_tool_params(tool_params, error_cls=ValueError) -> dict:
    """Validate and normalize the Stage 4c tool_parameters config block."""
    if tool_params is None:
        return {}

    if isinstance(tool_params, bool):
        raise error_cls("tool_parameters must be a mapping, got boolean")

    if not isinstance(tool_params, dict):
        raise error_cls(
            f"tool_parameters must be a mapping, got {type(tool_params).__name__}"
        )

    known_tools = defaults.TOOL_PARAMETERS_TOOLS
    known_keys = defaults.TOOL_PARAMETERS_KEYS

    def _normalize_bool(tool, key, raw):
        if isinstance(raw, bool):
            return raw
        val = str(raw).lower()
        if val in ("true", "false"):
            return val == "true"
        raise error_cls(
            f"tool_parameters.{tool}.{key} must be true or false, got {raw!r}"
        )

    def _normalize_int(tool, key, raw, minimum=0):
        if isinstance(raw, bool):
            raise error_cls(
                f"tool_parameters.{tool}.{key} must be an integer, got {raw!r}"
            )
        if isinstance(raw, int):
            parsed = raw
        elif isinstance(raw, str):
            if raw.strip() == "":
                return ""
            if not raw.strip().isdigit():
                raise error_cls(
                    f"tool_parameters.{tool}.{key} must be an integer, got {raw!r}"
                )
            parsed = int(raw.strip())
        else:
            raise error_cls(
                f"tool_parameters.{tool}.{key} must be an integer, got {raw!r}"
            )
        if parsed < minimum:
            raise error_cls(
                f"tool_parameters.{tool}.{key} must be non-negative, got {parsed}"
            )
        return parsed

    def _normalize_probability(tool, key, raw):
        label = f"tool_parameters.{tool}.{key}"
        if isinstance(raw, bool):
            raise error_cls(f"{label} must be a probability in (0, 1], got {raw!r}")
        if isinstance(raw, (int, float)):
            val = float(raw)
        elif isinstance(raw, str):
            try:
                val = float(raw.strip())
            except (ValueError, TypeError):
                raise error_cls(f"{label} must be a probability in (0, 1], got {raw!r}")
        else:
            raise error_cls(f"{label} must be a probability in (0, 1], got {raw!r}")
        if not math.isfinite(val) or val <= 0 or val > 1:
            raise error_cls(f"{label} must be a probability in (0, 1], got {val}")
        return val

    def _normalize_filter_flags(key, raw):
        if raw is None or raw == "":
            return ""
        if isinstance(raw, bool):
            raise error_cls(
                f"tool_parameters.samtools_filter.{key} must be "
                f"a positive integer or hex string like 0x904, "
                f"got {raw!r}"
            )
        if isinstance(raw, int):
            if raw <= 0:
                raise error_cls(
                    f"tool_parameters.samtools_filter.{key} must be positive, got {raw}"
                )
            return raw
        if isinstance(raw, str):
            text = raw.strip()
            if text.startswith("0x") or text.startswith("0X"):
                try:
                    parsed = int(text, 16)
                except ValueError:
                    raise error_cls(
                        f"tool_parameters.samtools_filter.{key}: "
                        f"invalid hex value {raw!r}"
                    )
                if parsed <= 0:
                    raise error_cls(
                        f"tool_parameters.samtools_filter.{key} must be "
                        f"positive, got {parsed} (from {raw!r})"
                    )
                return parsed
            if text.isdigit():
                parsed = int(text)
                if parsed <= 0:
                    raise error_cls(
                        f"tool_parameters.samtools_filter.{key} must be "
                        f"positive, got {parsed}"
                    )
                return parsed
            raise error_cls(
                f"tool_parameters.samtools_filter.{key} must be "
                f"a positive integer or hex string like 0x904, "
                f"got {raw!r}"
            )
        raise error_cls(
            f"tool_parameters.samtools_filter.{key} must be "
            f"a positive integer or hex string like 0x904, "
            f"got {raw!r}"
        )

    for tool in tool_params:
        if tool not in known_tools:
            raise error_cls(
                f"tool_parameters: unknown tool block {tool!r}. "
                f"Known: {sorted(known_tools)}"
            )

    normalized = {}

    for tool in known_tools:
        block = tool_params.get(tool, {})
        if not isinstance(block, dict):
            raise error_cls(
                f"tool_parameters.{tool} must be a mapping, got {type(block).__name__}"
            )

        for key in block:
            if key not in known_keys[tool]:
                raise error_cls(
                    f"tool_parameters.{tool}: unknown key {key!r}. "
                    f"Known: {sorted(known_keys[tool])}"
                )

        norm: dict = {}

        if tool == "fastqc":
            extra = block.get("extra_args", "")
            if not isinstance(extra, str):
                raise error_cls("tool_parameters.fastqc.extra_args must be a string")
            norm["extra_args"] = extra

        elif tool == "trim_galore":
            norm["quality"] = _normalize_int(tool, "quality", block.get("quality", ""))
            norm["length"] = _normalize_int(tool, "length", block.get("length", ""))
            norm["stringency"] = _normalize_int(
                tool, "stringency", block.get("stringency", "")
            )
            extra = block.get("extra_args", "")
            if not isinstance(extra, str):
                raise error_cls(
                    "tool_parameters.trim_galore.extra_args must be a string"
                )
            norm["extra_args"] = extra

        elif tool == "bowtie2":
            mode = block.get("mode", "")
            if mode != "" and not isinstance(mode, str):
                raise error_cls("tool_parameters.bowtie2.mode must be a string")
            if mode not in defaults.BOWTIE2_MODES:
                raise error_cls(
                    f"tool_parameters.bowtie2.mode must be one of: "
                    "very-fast, fast, sensitive, very-sensitive, "
                    f"or empty. Got {mode!r}"
                )
            norm["mode"] = mode
            norm["dovetail"] = _normalize_bool(
                tool, "dovetail", block.get("dovetail", False)
            )
            norm["no_mixed"] = _normalize_bool(
                tool, "no_mixed", block.get("no_mixed", False)
            )
            norm["no_discordant"] = _normalize_bool(
                tool, "no_discordant", block.get("no_discordant", False)
            )
            extra = block.get("extra_args", "")
            if not isinstance(extra, str):
                raise error_cls("tool_parameters.bowtie2.extra_args must be a string")
            norm["extra_args"] = extra

        elif tool == "samtools_filter":
            norm["filter_flags"] = _normalize_filter_flags(
                "filter_flags", block.get("filter_flags", "")
            )
            extra = block.get("extra_args", "")
            if not isinstance(extra, str):
                raise error_cls(
                    "tool_parameters.samtools_filter.extra_args must be a string"
                )
            norm["extra_args"] = extra

        elif tool == "picard_markduplicates":
            norm["optical_duplicate_pixel_distance"] = _normalize_int(
                tool,
                "optical_duplicate_pixel_distance",
                block.get("optical_duplicate_pixel_distance", ""),
            )
            extra = block.get("extra_args", "")
            if not isinstance(extra, str):
                raise error_cls(
                    "tool_parameters.picard_markduplicates.extra_args must be a string"
                )
            norm["extra_args"] = extra

        elif tool == "bamcoverage":
            norm_using = block.get("normalize_using", "CPM")
            if not isinstance(norm_using, str):
                raise error_cls(
                    "tool_parameters.bamcoverage.normalize_using must be a string"
                )
            allowed_norm = defaults.BAMCOVERAGE_NORMALIZE_USING
            if norm_using not in allowed_norm:
                raise error_cls(
                    f"tool_parameters.bamcoverage.normalize_using must be "
                    f"one of: {sorted(allowed_norm)}. "
                    f"RPGC is not supported in this slice. Got {norm_using!r}"
                )
            norm["normalize_using"] = norm_using
            norm["smooth_length"] = _normalize_int(
                tool, "smooth_length", block.get("smooth_length", ""), minimum=1
            )
            extra = block.get("extra_args", "")
            if not isinstance(extra, str):
                raise error_cls(
                    "tool_parameters.bamcoverage.extra_args must be a string"
                )
            norm["extra_args"] = extra

        elif tool == "macs3":
            norm["qvalue"] = _normalize_probability(
                tool, "qvalue", block.get("qvalue", 0.01)
            )
            norm["broad_cutoff"] = _normalize_probability(
                tool, "broad_cutoff", block.get("broad_cutoff", 0.1)
            )
            extra = block.get("extra_args", "")
            if not isinstance(extra, str):
                raise error_cls("tool_parameters.macs3.extra_args must be a string")
            norm["extra_args"] = _migrate_keep_dup_extra_args(tool, extra, error_cls)

        elif tool == "multiqc":
            title = block.get("title", "")
            if not isinstance(title, str):
                raise error_cls("tool_parameters.multiqc.title must be a string")
            norm["title"] = title
            extra = block.get("extra_args", "")
            if not isinstance(extra, str):
                raise error_cls("tool_parameters.multiqc.extra_args must be a string")
            norm["extra_args"] = extra

        elif tool == "idr_macs3":
            norm["pvalue"] = _normalize_probability(
                tool, "pvalue", block.get("pvalue", 0.1)
            )
            extra = block.get("extra_args", "")
            if not isinstance(extra, str):
                raise error_cls("tool_parameters.idr_macs3.extra_args must be a string")
            norm["extra_args"] = _migrate_keep_dup_extra_args(tool, extra, error_cls)

        normalized[tool] = norm

    return normalized
