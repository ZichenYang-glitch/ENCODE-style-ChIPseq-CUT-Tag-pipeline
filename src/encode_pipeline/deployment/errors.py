"""Stable, redacted failures for supported deployment operations."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class DeploymentIssue:
    """One public operator issue without paths, secrets, or exception text."""

    code: str
    message: str
    component: str | None = None
    recoverable: bool = False

    def to_dict(self) -> dict[str, object]:
        value: dict[str, object] = {
            "code": self.code,
            "message": self.message,
            "recoverable": self.recoverable,
        }
        if self.component is not None:
            value["component"] = self.component
        return value


class DeploymentError(RuntimeError):
    """Internal exception carrying only a stable public issue."""

    def __init__(self, issue: DeploymentIssue) -> None:
        self.issue = issue
        super().__init__(issue.code)


def fail(
    code: str,
    message: str,
    *,
    component: str | None = None,
    recoverable: bool = False,
) -> DeploymentError:
    return DeploymentError(
        DeploymentIssue(
            code=code,
            message=message,
            component=component,
            recoverable=recoverable,
        )
    )


__all__ = ["DeploymentError", "DeploymentIssue", "fail"]
