"""Tests for the WorkspaceMaterializer filesystem materialization boundary."""

from pathlib import Path


def test_workspace_materializer_refuses_relative_base_dir(tmp_path):
    from encode_pipeline.platform.adapters import WorkspacePlan
    from encode_pipeline.services.materialization import WorkspaceMaterializer

    materializer = WorkspaceMaterializer()
    result = materializer.materialize(WorkspacePlan(), Path("relative/path"))

    assert result.is_failure is True
    issue = result.issues[0]
    assert issue.code == "WORKSPACE_MATERIALIZATION_BASE_DIR_RELATIVE"
    assert issue.path == "base_dir"


def test_workspace_materializer_refuses_symlink_base_dir(tmp_path):
    from encode_pipeline.platform.adapters import WorkspacePlan
    from encode_pipeline.services.materialization import WorkspaceMaterializer

    real_dir = tmp_path / "real"
    real_dir.mkdir()
    symlink_dir = tmp_path / "link"
    symlink_dir.symlink_to(real_dir)

    materializer = WorkspaceMaterializer()
    result = materializer.materialize(WorkspacePlan(), symlink_dir)

    assert result.is_failure is True
    issue = result.issues[0]
    assert issue.code == "WORKSPACE_MATERIALIZATION_SYMLINK"
    assert issue.path == "base_dir"


def test_workspace_materializer_refuses_symlink_parent_of_base_dir(tmp_path):
    from encode_pipeline.platform.adapters import WorkspacePlan
    from encode_pipeline.services.materialization import WorkspaceMaterializer

    real_dir = tmp_path / "real"
    real_dir.mkdir()
    symlink_dir = tmp_path / "link"
    symlink_dir.symlink_to(real_dir)
    base_dir = symlink_dir / "workspace"

    materializer = WorkspaceMaterializer()
    result = materializer.materialize(WorkspacePlan(), base_dir)

    assert result.is_failure is True
    issue = result.issues[0]
    assert issue.code == "WORKSPACE_MATERIALIZATION_SYMLINK"
    assert issue.path == "base_dir"


def test_workspace_materializer_refuses_broken_symlink_base_dir(tmp_path):
    from encode_pipeline.platform.adapters import WorkspacePlan
    from encode_pipeline.services.materialization import WorkspaceMaterializer

    broken_link = tmp_path / "broken_link"
    broken_link.symlink_to(tmp_path / "does_not_exist")

    materializer = WorkspaceMaterializer()
    result = materializer.materialize(WorkspacePlan(), broken_link)

    assert result.is_failure is True
    issue = result.issues[0]
    assert issue.code == "WORKSPACE_MATERIALIZATION_SYMLINK"
    assert issue.path == "base_dir"


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
