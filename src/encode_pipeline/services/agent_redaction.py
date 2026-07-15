from __future__ import annotations

import re
from collections.abc import Mapping
from dataclasses import dataclass


SENSITIVE_KEYS: frozenset[str] = frozenset(
    {
        "password",
        "token",
        "secret",
        "api_key",
        "access_key",
        "private_key",
        "credential",
    }
)

_PATH_TERMINATOR = r"(?:\s+(?=\w+(?:[.,;:!?](?=\s|[()\[\]{}\"\'\`<>]|$)|\s|[()\[\]{}\"\'\`<>]|$))|\s+(?=[()\[\]{}\"\'\`<>])|(?=[()\[\]{}\"\'\`<>])|[.,;:!?]\s|[.,;:!?]$|$)"

# Combined absolute-path regex. Unix paths start with "/" and are not preceded
# by a word char, another "/", or a colon (so URL schemes like https:// are
# preserved). A leading "//" is treated as a protocol-relative URL because the second "//" character is excluded.
# The "/api" exclusion is case-insensitive and uses a word boundary
# so "/api/...", "/api?...", and bare "/api" are preserved while "/apiary" is
# still redacted. Windows paths are drive-letter (C:\...) or UNC (\\server...).
# The match continues across spaces when the following token looks like a path
# component (contains "/", ".", or "-") and stops before trailing sentence
# punctuation or common bracket/quote delimiters.
_PATH_RE = re.compile(
    r"(?<![\w/:.~])/(?!/)(?!api\b).*?(?=" + _PATH_TERMINATOR + r")"
    r"|(?<!\w)(?:[A-Za-z]:[\\/](?!/)|\\\\).*?(?=" + _PATH_TERMINATOR + r")",
    re.IGNORECASE,
)

# Maximum container nesting depth handled by recursive redaction. Kept well
# below Python's default recursion limit so the depth guard itself does not crash.
_MAX_DEPTH = 200


def _is_sensitive_key(key: object) -> bool:
    return isinstance(key, str) and key.lower() in SENSITIVE_KEYS


@dataclass(frozen=True)
class RedactionResult:
    value: object
    sensitive_keys_redacted: int = 0
    absolute_paths_redacted: int = 0


class _RedactionContext:
    def __init__(self) -> None:
        self.seen_ids: set[int] = set()


class RedactionPolicy:
    def redact_value(self, value: object) -> RedactionResult:
        ctx = _RedactionContext()
        return self._redact(value, ctx, depth=0)

    def redact_text(self, text: str) -> RedactionResult:
        matches = _PATH_RE.findall(text)
        redacted = _PATH_RE.sub("<FILE_PATH>", text)
        return RedactionResult(
            value=redacted,
            absolute_paths_redacted=len(matches),
        )

    def redact_issue(self, issue: Mapping[str, object]) -> RedactionResult:
        preserved = {"code", "severity", "source"}
        redacted: dict[str, object] = {}
        total_sensitive = 0
        total_paths = 0
        for key, val in issue.items():
            if key in preserved:
                redacted[key] = val
                continue
            if key == "path" and isinstance(val, str):
                path_result = self.redact_text(val)
                redacted[key] = path_result.value
                total_paths += path_result.absolute_paths_redacted
                continue
            result = self.redact_value(val)
            redacted[key] = result.value
            total_sensitive += result.sensitive_keys_redacted
            total_paths += result.absolute_paths_redacted
        return RedactionResult(
            value=redacted,
            sensitive_keys_redacted=total_sensitive,
            absolute_paths_redacted=total_paths,
        )

    def redact_issues(self, issues: list[Mapping[str, object]]) -> RedactionResult:
        redacted: list[object] = []
        total_sensitive = 0
        total_paths = 0
        for issue in issues:
            result = self.redact_issue(issue)
            redacted.append(result.value)
            total_sensitive += result.sensitive_keys_redacted
            total_paths += result.absolute_paths_redacted
        return RedactionResult(
            value=redacted,
            sensitive_keys_redacted=total_sensitive,
            absolute_paths_redacted=total_paths,
        )

    def _redact(
        self, value: object, ctx: _RedactionContext, depth: int
    ) -> RedactionResult:
        if depth > _MAX_DEPTH:
            return RedactionResult(value="<MAX_DEPTH_REACHED>")
        if isinstance(value, str):
            return self.redact_text(value)
        if isinstance(value, Mapping):
            obj_id = id(value)
            if obj_id in ctx.seen_ids:
                return RedactionResult(value="<CYCLIC_REFERENCE>")
            ctx.seen_ids.add(obj_id)
            try:
                return self._redact_mapping(value, ctx, depth)
            finally:
                ctx.seen_ids.discard(obj_id)
        if isinstance(value, list):
            obj_id = id(value)
            if obj_id in ctx.seen_ids:
                return RedactionResult(value="<CYCLIC_REFERENCE>")
            ctx.seen_ids.add(obj_id)
            try:
                return self._redact_list(value, ctx, depth)
            finally:
                ctx.seen_ids.discard(obj_id)
        return RedactionResult(value=value)

    def _redact_mapping(
        self, mapping: Mapping[object, object], ctx: _RedactionContext, depth: int
    ) -> RedactionResult:
        redacted: dict[object, object] = {}
        total_sensitive = 0
        total_paths = 0
        for key, val in mapping.items():
            if _is_sensitive_key(key):
                redacted[key] = "<REDACTED>"
                total_sensitive += 1
            else:
                result = self._redact(val, ctx, depth + 1)
                redacted[key] = result.value
                total_sensitive += result.sensitive_keys_redacted
                total_paths += result.absolute_paths_redacted
        return RedactionResult(
            value=redacted,
            sensitive_keys_redacted=total_sensitive,
            absolute_paths_redacted=total_paths,
        )

    def _redact_list(
        self, items: list[object], ctx: _RedactionContext, depth: int
    ) -> RedactionResult:
        redacted: list[object] = []
        total_sensitive = 0
        total_paths = 0
        for item in items:
            result = self._redact(item, ctx, depth + 1)
            redacted.append(result.value)
            total_sensitive += result.sensitive_keys_redacted
            total_paths += result.absolute_paths_redacted
        return RedactionResult(
            value=redacted,
            sensitive_keys_redacted=total_sensitive,
            absolute_paths_redacted=total_paths,
        )
