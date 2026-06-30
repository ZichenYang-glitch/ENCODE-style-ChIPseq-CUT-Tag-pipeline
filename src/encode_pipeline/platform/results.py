"""Workflow-platform Result and Issue primitives."""

from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Generic, TypeVar


T = TypeVar("T")


class IssueSeverity(str, Enum):
    """Severity levels for workflow-platform issues."""

    ERROR = "error"
    WARNING = "warning"
    INFO = "info"


@dataclass(frozen=True)
class Issue:
    """A structured user-facing issue produced by platform boundaries."""

    code: str
    message: str
    severity: IssueSeverity | str = IssueSeverity.ERROR
    path: str | None = None
    source: str | None = None
    technical_message: str | None = None
    hint: str | None = None
    context: Mapping[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        code = _normalize_required_string(self.code, "code")
        message = _normalize_required_string(self.message, "message")
        severity = _normalize_severity(self.severity)
        context = _copy_context(self.context)

        object.__setattr__(self, "code", code)
        object.__setattr__(self, "message", message)
        object.__setattr__(self, "severity", severity)
        object.__setattr__(self, "context", context)

    @classmethod
    def from_exception(
        cls,
        exc: Exception,
        *,
        code: str = "VALIDATION_ERROR",
        source: str | None = None,
        path: str | None = None,
        hint: str | None = None,
        context: Mapping[str, object] | None = None,
    ) -> "Issue":
        """Build an error issue from any exception without importing adapters."""
        message = str(exc)
        return cls(
            code=code,
            message=message,
            source=source,
            path=path,
            technical_message=message,
            hint=hint,
            context=context or {},
        )

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-ready dict with a fresh context copy."""
        return {
            "code": self.code,
            "message": self.message,
            "severity": self.severity.value,
            "path": self.path,
            "source": self.source,
            "technical_message": self.technical_message,
            "hint": self.hint,
            "context": dict(self.context),
        }


@dataclass(frozen=True)
class Result(Generic[T]):
    """A value plus structured issues from platform/service boundaries."""

    value: T | None
    issues: tuple[Issue, ...] = ()

    def __post_init__(self) -> None:
        object.__setattr__(self, "issues", _normalize_issues(self.issues))

    @classmethod
    def success(
        cls,
        value: T,
        issues: Iterable[Issue] = (),
    ) -> "Result[T]":
        """Return a successful result, optionally with warning/info issues."""
        issue_tuple = _normalize_issues(issues)
        if any(issue.severity is IssueSeverity.ERROR for issue in issue_tuple):
            raise ValueError("Result.success cannot include error issues")
        return cls(value=value, issues=issue_tuple)

    @classmethod
    def failure(cls, issues: Iterable[Issue]) -> "Result[T]":
        """Return a failed result requiring at least one error issue."""
        issue_tuple = _normalize_issues(issues)
        if not any(issue.severity is IssueSeverity.ERROR for issue in issue_tuple):
            raise ValueError("Result.failure requires at least one error issue")
        return cls(value=None, issues=issue_tuple)

    @property
    def is_success(self) -> bool:
        """Return True when no error-severity issues are present."""
        return not self.is_failure

    @property
    def is_failure(self) -> bool:
        """Return True when any error-severity issue is present."""
        return any(issue.severity is IssueSeverity.ERROR for issue in self.issues)

    @property
    def errors(self) -> tuple[Issue, ...]:
        """Return error-severity issues."""
        return tuple(
            issue for issue in self.issues
            if issue.severity is IssueSeverity.ERROR
        )

    @property
    def warnings(self) -> tuple[Issue, ...]:
        """Return warning-severity issues."""
        return tuple(
            issue for issue in self.issues
            if issue.severity is IssueSeverity.WARNING
        )

    def to_dict(
        self,
        *,
        value_serializer: Callable[[T | None], Any] | None = None,
    ) -> dict[str, Any]:
        """Return a JSON-ready dict for API/frontend/agent consumers."""
        value: Any = self.value
        if value_serializer is not None:
            value = value_serializer(self.value)
        return {
            "ok": self.is_success,
            "value": value,
            "issues": [issue.to_dict() for issue in self.issues],
        }


def _normalize_required_string(value: str, name: str) -> str:
    if not isinstance(value, str):
        raise ValueError(f"Issue {name} must be a string")
    normalized = value.strip()
    if not normalized:
        raise ValueError(f"Issue {name} must be non-empty")
    return normalized


def _normalize_severity(value: IssueSeverity | str) -> IssueSeverity:
    if isinstance(value, IssueSeverity):
        return value
    if isinstance(value, str):
        try:
            return IssueSeverity(value)
        except ValueError as exc:
            raise ValueError(f"Invalid issue severity: {value!r}") from exc
    raise ValueError(f"Invalid issue severity: {value!r}")


def _copy_context(context: Mapping[str, object]) -> dict[str, object]:
    if not isinstance(context, Mapping):
        raise ValueError("Issue context must be a mapping")
    return dict(context)


def _normalize_issues(issues: Iterable[Issue]) -> tuple[Issue, ...]:
    issue_tuple = tuple(issues)
    for issue in issue_tuple:
        if not isinstance(issue, Issue):
            raise ValueError("Result issues must contain only Issue entries")
    return issue_tuple
