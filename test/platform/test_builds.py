"""Tests for durable workflow build identity primitives."""

from datetime import datetime, timezone

import pytest

from encode_pipeline.platform.builds import WorkflowBuildIdentity


def _identity(**changes) -> WorkflowBuildIdentity:
    values = {
        "workflow_id": "workflow",
        "adapter_version": "1.0.0",
        "scheme": "sha256-tree-v1",
        "logical_entrypoint": "workflow/Snakefile",
        "digest": "a" * 64,
        "captured_at": datetime.now(timezone.utc),
    }
    values.update(changes)
    return WorkflowBuildIdentity(**values)


def test_workflow_build_identity_matches_without_timestamp_identity():
    first = _identity(captured_at=datetime(2026, 1, 1, tzinfo=timezone.utc))
    second = _identity(captured_at=datetime(2026, 1, 2, tzinfo=timezone.utc))

    assert first != second
    assert first.matches(second)


@pytest.mark.parametrize(
    ("field_name", "value"),
    [
        ("workflow_id", " "),
        ("adapter_version", ""),
        ("scheme", ""),
        ("logical_entrypoint", "/workflow/Snakefile"),
        ("logical_entrypoint", "workflow/../Snakefile"),
        ("digest", "A" * 64),
        ("digest", "a" * 63),
    ],
)
def test_workflow_build_identity_rejects_invalid_fields(field_name, value):
    with pytest.raises(ValueError):
        _identity(**{field_name: value})
