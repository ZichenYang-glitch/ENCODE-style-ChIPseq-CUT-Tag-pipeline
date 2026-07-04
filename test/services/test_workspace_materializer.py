"""Tests for the WorkspaceMaterializer filesystem materialization boundary."""

from pathlib import Path


def test_workspace_materializer_refuses_invalid_plan_type(tmp_path):
    from encode_pipeline.services.materialization import WorkspaceMaterializer

    materializer = WorkspaceMaterializer()
    result = materializer.materialize("not-a-plan", tmp_path.resolve())

    assert result.is_failure is True
    assert len(result.issues) == 1
    issue = result.issues[0]
    assert issue.code == "WORKSPACE_MATERIALIZATION_INVALID_PLAN"
    assert issue.severity.value == "error"
    assert issue.source == "workspace_materializer"
    assert issue.path == "plan"
