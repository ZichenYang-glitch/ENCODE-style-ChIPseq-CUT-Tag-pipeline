"""Tool parameter validation helpers."""

import math
import re

from encode_pipeline.config import defaults

__all__ = [
    "extract_bowtie2_min_insert",
    "validate_bowtie2_insert_bounds",
    "validate_tool_params",
]


_KEEP_DUP_CANONICAL = "--keep-dup"
_BOWTIE2_MAXINS_CANONICAL = "--maxins"
_BOWTIE2_MININS_CANONICAL = "--minins"
_BUNDLED_BOWTIE2_X_RE = re.compile(r"^-[bfqzhcradtC]+X")
_BOWTIE2_SHORT_I_RE = re.compile(r"^-[bfqzhcradtC]*I(.*)$")

_EXTRA_ARGS_WHITESPACE = " \t\n\r"


def _is_keep_dup_option(name):
    """True when ``name`` selects MACS3's ``--keep-dup`` option.

    The pinned MACS3 uses argparse, which accepts any unambiguous prefix of
    a long option; ``--keep-dup`` is the only callpeak option starting with
    ``--k``, so every prefix from ``--k`` upwards (e.g. ``--keep-d``) is the
    same flag. Longer spellings such as ``--keep-duplicates`` are not
    prefixes of the canonical option and are left alone.
    """
    return len(name) >= 3 and _KEEP_DUP_CANONICAL.startswith(name)


def _is_bowtie2_maxins_option(name):
    """True when ``name`` selects Bowtie2's ``--maxins`` long option.

    In the pinned Bowtie2 (GNU getopt_long in ``bowtie2-align-s``), ``--ma``
    is an exact option (match bonus), while every prefix from ``--max``
    through ``--maxins`` unambiguously resolves to ``--maxins``.
    """
    return len(name) >= 5 and _BOWTIE2_MAXINS_CANONICAL.startswith(name)


def _is_bowtie2_minins_option(name):
    """True when ``name`` selects Bowtie2's ``--minins`` long option.

    In the pinned Bowtie2, ``--min`` is ambiguous between ``--minins`` and
    ``--min-score``, while every prefix from ``--mini`` through ``--minins``
    unambiguously resolves to ``--minins``.
    """
    return len(name) >= 6 and _BOWTIE2_MININS_CANONICAL.startswith(name)


def _scan_extra_args_tokens(tool, extra, error_cls):
    """Tokenize extra_args into ``(value, start, end)`` spans.

    Follows POSIX quoting for token boundaries: single quotes quote
    literally, unquoted backslashes escape the next character except that
    ``\\`` + newline is removed as a line continuation, and double quotes
    group while still folding ``\\`` + newline continuations and ``$``,
    `` ` ``, ``"``, ``\\`` escapes. ``value`` is the dequoted token used for
    flag comparisons; ``start``/``end`` delimit the raw source text of each
    token (excluding any trailing unquoted line continuation) so untouched
    tokens keep their original quoting and escaped spaces byte-for-byte.
    Unterminated quotes raise a field-named validation error instead of
    leaking parser exceptions.
    """

    def _invalid_quoting():
        return error_cls(
            f"tool_parameters.{tool}.extra_args: invalid shell quoting "
            f"in {extra!r}: No closing quotation"
        )

    tokens = []
    index = 0
    length = len(extra)
    while index < length:
        if extra[index] in _EXTRA_ARGS_WHITESPACE:
            index += 1
            continue
        if extra[index] == "\\" and index + 1 < length and extra[index + 1] == "\n":
            index += 2
            continue
        start = index
        token_end = index
        chars = []
        while index < length and extra[index] not in _EXTRA_ARGS_WHITESPACE:
            char = extra[index]
            if char == "\\":
                if index + 1 < length and extra[index + 1] == "\n":
                    index += 2
                elif index + 1 < length:
                    chars.append(extra[index + 1])
                    index += 2
                    token_end = index
                else:
                    chars.append(char)
                    index += 1
                    token_end = index
            elif char == "'":
                closing = extra.find("'", index + 1)
                if closing == -1:
                    raise _invalid_quoting()
                chars.append(extra[index + 1 : closing])
                index = closing + 1
                token_end = index
            elif char == '"':
                cursor = index + 1
                while cursor < length and extra[cursor] != '"':
                    if extra[cursor] == "\\" and cursor + 1 < length:
                        next_char = extra[cursor + 1]
                        if next_char == "\n":
                            cursor += 2
                            continue
                        if next_char in '$`"\\':
                            chars.append(next_char)
                            cursor += 2
                            continue
                    chars.append(extra[cursor])
                    cursor += 1
                if cursor >= length:
                    raise _invalid_quoting()
                index = cursor + 1
                token_end = index
            else:
                chars.append(char)
                index += 1
                token_end = index
        tokens.append(("".join(chars), start, token_end))
    return tokens


def _migrate_keep_dup_extra_args(tool, extra, error_cls):
    """Normalize legacy ``--keep-dup`` tokens in MACS3 extra_args.

    Duplicate policy is executed at the BAM layer and the workflow renders
    ``--keep-dup all`` itself. A synonymous ``--keep-dup all`` (including
    argparse-accepted abbreviations and the ``=`` form) is removed here so
    the final command carries the flag exactly once; any other value
    conflicts with that contract and is rejected with migration guidance
    instead of being decided by argument order.

    Only the raw source spans of the managed tokens are removed; every
    other argument keeps its original quoting and shell expansions
    byte-for-byte (e.g. ``--tempdir "$TMPDIR"`` still expands), and a
    second validation pass is a stable no-op.
    """
    tokens = _scan_extra_args_tokens(tool, extra, error_cls)
    removed = set()
    migrated = False
    index = 0
    while index < len(tokens):
        value, _start, _end = tokens[index]
        name, sep, attached = value.partition("=")
        if not _is_keep_dup_option(name):
            index += 1
            continue
        migrated = True
        removed.add(index)
        if sep:
            flag_value = attached
            index += 1
        elif index + 1 < len(tokens):
            flag_value = tokens[index + 1][0]
            removed.add(index + 1)
            index += 2
        else:
            flag_value = ""
            index += 1
        if flag_value != "all":
            raise error_cls(
                f"tool_parameters.{tool}.extra_args: {name} "
                f"{flag_value or '<missing>'!r} conflicts with the BAM-layer "
                "duplicate policy. Remove --keep-dup from extra_args (the "
                "workflow always renders --keep-dup all) or use 'all'."
            )
    if not migrated:
        return extra
    return " ".join(
        extra[start:end]
        for token_index, (_value, start, end) in enumerate(tokens)
        if token_index not in removed
    )


def _normalize_bowtie2_max_insert(raw, error_cls):
    """Validate and normalize ``tool_parameters.bowtie2.max_insert``.

    Accepts ``"auto"`` or a positive integer (``int`` or ASCII digit string).
    Rejects booleans, zero, negative integers, decimals/floats, empty strings,
    and other values.
    """
    label = "tool_parameters.bowtie2.max_insert"
    if isinstance(raw, bool):
        raise error_cls(f"{label} must be 'auto' or a positive integer, got {raw!r}")
    if isinstance(raw, int):
        if raw < 1:
            raise error_cls(
                f"{label} must be 'auto' or a positive integer, got {raw!r}"
            )
        return raw
    if isinstance(raw, str):
        if raw == "auto":
            return "auto"
        text = raw.strip()
        if text.isascii() and text.isdigit():
            parsed = int(text)
            if parsed >= 1:
                return parsed
        raise error_cls(f"{label} must be 'auto' or a positive integer, got {raw!r}")
    raise error_cls(f"{label} must be 'auto' or a positive integer, got {raw!r}")


def _migrate_bowtie2_extra_args(max_insert, extra, error_cls):
    """Normalize ``-X``/``--maxins`` and validate ``-I``/``--minins`` in Bowtie2 extra_args.

    - Legacy ``-X`` / ``--maxins`` (including attached ``-X<N>``, ``=`` form,
      and unambiguous ``getopt_long`` abbreviations ``--max``, ``--maxi``,
      ``--maxin``, ``--maxins``) is migrated into ``max_insert`` when
      ``max_insert`` is ``"auto"``, or stripped when it matches an explicit
      ``max_insert``.
    - Conflicting ``-X``/``--maxins`` values (either against ``max_insert`` or
      multiple distinct values in ``extra_args``), bundled ``-<flags>X``, or
      invalid values (``<= 0``, ``-X=N``, non-integers, missing argument) are
      rejected with migration guidance instead of relying on argument order.
    - ``-I`` / ``--minins`` (including attached ``-I<N>``, bundled
      ``-<flags>I``, ``=`` form, and unambiguous ``--mini``, ``--minin``,
      ``--minins``) is validated as a non-negative integer and checked against
      ``effective_max_insert`` while remaining in ``extra_args``.
    """
    tokens = _scan_extra_args_tokens("bowtie2", extra, error_cls)
    removed = set()
    migrated_max_insert = None
    seen_min_insert = None
    index = 0
    while index < len(tokens):
        value, _start, _end = tokens[index]
        name, sep, attached = value.partition("=")

        if _BUNDLED_BOWTIE2_X_RE.match(value):
            raise error_cls(
                f"tool_parameters.bowtie2.extra_args: bundled {value!r} "
                "cannot be combined with workflow-managed max_insert "
                "(-X/--maxins). Configure tool_parameters.bowtie2.max_insert "
                "and separate short flags."
            )

        is_short_x = value == "-X" or (
            value.startswith("-X") and not value.startswith("--")
        )
        is_long_maxins = _is_bowtie2_maxins_option(name)
        if is_short_x or is_long_maxins:
            removed.add(index)
            if is_short_x:
                opt_label = "-X"
                if value == "-X":
                    if index + 1 < len(tokens):
                        flag_value = tokens[index + 1][0]
                        removed.add(index + 1)
                        index += 2
                    else:
                        flag_value = ""
                        index += 1
                else:
                    flag_value = value[2:]
                    index += 1
            else:
                opt_label = name
                if sep:
                    flag_value = attached
                    index += 1
                elif index + 1 < len(tokens):
                    flag_value = tokens[index + 1][0]
                    removed.add(index + 1)
                    index += 2
                else:
                    flag_value = ""
                    index += 1

            if not (
                flag_value.isascii() and flag_value.isdigit() and int(flag_value) >= 1
            ):
                raise error_cls(
                    f"tool_parameters.bowtie2.extra_args: {opt_label} requires "
                    f"a positive integer, got {flag_value or '<missing>'!r}. "
                    "Set tool_parameters.bowtie2.max_insert to a positive "
                    "integer or 'auto' and remove -X/--maxins from extra_args."
                )
            parsed_max = int(flag_value)
            if migrated_max_insert is not None and migrated_max_insert != parsed_max:
                raise error_cls(
                    f"tool_parameters.bowtie2.extra_args: conflicting "
                    f"{opt_label} values ({migrated_max_insert} vs "
                    f"{parsed_max}). Set tool_parameters.bowtie2.max_insert "
                    "in one place and remove -X/--maxins from extra_args."
                )
            migrated_max_insert = parsed_max
            continue

        if name == "--min":
            raise error_cls(
                "tool_parameters.bowtie2.extra_args: '--min' is ambiguous in "
                "Bowtie2 (matches --minins and --min-score); use --minins or -I."
            )

        short_i_match = (
            None if value.startswith("--") else _BOWTIE2_SHORT_I_RE.match(value)
        )
        is_long_minins = _is_bowtie2_minins_option(name)
        if short_i_match is not None or is_long_minins:
            if short_i_match is not None:
                opt_label = "-I"
                attached_i = short_i_match.group(1)
                if attached_i == "":
                    if index + 1 < len(tokens):
                        flag_value = tokens[index + 1][0]
                        index += 2
                    else:
                        flag_value = ""
                        index += 1
                else:
                    flag_value = attached_i
                    index += 1
            else:
                opt_label = name
                if sep:
                    flag_value = attached
                    index += 1
                elif index + 1 < len(tokens):
                    flag_value = tokens[index + 1][0]
                    index += 2
                else:
                    flag_value = ""
                    index += 1

            if not (flag_value.isascii() and flag_value.isdigit()):
                raise error_cls(
                    f"tool_parameters.bowtie2.extra_args: {opt_label} requires "
                    f"a non-negative integer, got {flag_value or '<missing>'!r}."
                )
            parsed_min = int(flag_value)
            if seen_min_insert is not None and seen_min_insert != parsed_min:
                raise error_cls(
                    f"tool_parameters.bowtie2.extra_args: conflicting "
                    f"{opt_label} values ({seen_min_insert} vs {parsed_min})."
                )
            seen_min_insert = parsed_min
            continue

        index += 1

    effective_max_insert = max_insert
    if migrated_max_insert is not None:
        if max_insert == "auto":
            effective_max_insert = migrated_max_insert
        elif max_insert != migrated_max_insert:
            raise error_cls(
                f"tool_parameters.bowtie2.extra_args: -X/--maxins "
                f"{migrated_max_insert} conflicts with "
                f"tool_parameters.bowtie2.max_insert={max_insert}. "
                "Configure tool_parameters.bowtie2.max_insert in one place "
                "and remove -X/--maxins from extra_args."
            )

    if seen_min_insert is not None:
        if isinstance(effective_max_insert, int):
            if seen_min_insert > effective_max_insert:
                raise error_cls(
                    f"tool_parameters.bowtie2.extra_args: minimum insert "
                    f"length ({seen_min_insert}) cannot exceed max_insert "
                    f"({effective_max_insert})."
                )
        else:
            max_auto = max(defaults.BOWTIE2_MAX_INSERT_BY_ASSAY.values())
            if seen_min_insert > max_auto:
                raise error_cls(
                    f"tool_parameters.bowtie2.extra_args: minimum insert "
                    f"length ({seen_min_insert}) exceeds maximum auto "
                    f"max_insert ({max_auto})."
                )

    if not removed:
        normalized_extra = extra
    else:
        normalized_extra = " ".join(
            extra[start:end]
            for token_index, (_value, start, end) in enumerate(tokens)
            if token_index not in removed
        )
    return effective_max_insert, normalized_extra


def extract_bowtie2_min_insert(extra, error_cls=ValueError):
    """Return the validated ``-I``/``--minins`` integer from ``extra``, or None."""
    if not extra:
        return None
    tokens = _scan_extra_args_tokens("bowtie2", extra, error_cls)
    seen_min_insert = None
    index = 0
    while index < len(tokens):
        value, _start, _end = tokens[index]
        name, sep, attached = value.partition("=")
        short_i_match = (
            None if value.startswith("--") else _BOWTIE2_SHORT_I_RE.match(value)
        )
        is_long_minins = _is_bowtie2_minins_option(name)
        if short_i_match is not None:
            attached_i = short_i_match.group(1)
            if attached_i == "":
                flag_value = tokens[index + 1][0] if index + 1 < len(tokens) else ""
                index += 2
            else:
                flag_value = attached_i
                index += 1
            seen_min_insert = int(flag_value)
            continue
        if is_long_minins:
            if sep:
                flag_value = attached
                index += 1
            else:
                flag_value = tokens[index + 1][0] if index + 1 < len(tokens) else ""
                index += 2
            seen_min_insert = int(flag_value)
            continue
        index += 1
    return seen_min_insert


def validate_bowtie2_insert_bounds(
    validated_config: dict,
    samples: list[dict],
    error_cls=ValueError,
) -> None:
    """Ensure Bowtie2 ``-I``/``--minins`` does not exceed any PE sample's resolved max_insert."""
    bowtie2_cfg = validated_config.get("tool_parameters", {}).get("bowtie2", {})
    extra = bowtie2_cfg.get("extra_args", "")
    min_insert = extract_bowtie2_min_insert(extra, error_cls=error_cls)
    if min_insert is None:
        return
    max_insert = bowtie2_cfg.get("max_insert", "auto")
    for sample in samples:
        if sample.get("layout") != "PE":
            continue
        assay = sample.get("assay", "")
        resolved_max = (
            max_insert
            if isinstance(max_insert, int)
            else defaults.BOWTIE2_MAX_INSERT_BY_ASSAY[assay]
        )
        if min_insert > resolved_max:
            sample_id = sample.get("id") or sample.get("sample") or "<unknown>"
            raise error_cls(
                f"tool_parameters.bowtie2.extra_args: minimum insert length "
                f"({min_insert}) cannot exceed resolved max_insert "
                f"({resolved_max}) for PE sample {sample_id!r} (assay={assay!r})."
            )


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
            max_insert = _normalize_bowtie2_max_insert(
                block.get("max_insert", "auto"), error_cls
            )
            extra = block.get("extra_args", "")
            if not isinstance(extra, str):
                raise error_cls("tool_parameters.bowtie2.extra_args must be a string")
            norm["max_insert"], norm["extra_args"] = _migrate_bowtie2_extra_args(
                max_insert, extra, error_cls
            )

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
