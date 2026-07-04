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


def test_workspace_materializer_creates_directories(tmp_path):
    from encode_pipeline.platform.adapters import WorkspacePlan
    from encode_pipeline.services.materialization import WorkspaceMaterializer

    base_dir = tmp_path.resolve()
    plan = WorkspacePlan(directories=("logs", "results/peaks"))

    materializer = WorkspaceMaterializer()
    result = materializer.materialize(plan, base_dir)

    assert result.is_success is True
    assert (base_dir / "logs").is_dir()
    assert (base_dir / "results" / "peaks").is_dir()


def test_workspace_materializer_existing_directories_are_ok(tmp_path):
    from encode_pipeline.platform.adapters import WorkspacePlan
    from encode_pipeline.services.materialization import WorkspaceMaterializer

    base_dir = tmp_path.resolve()
    (base_dir / "logs").mkdir()
    plan = WorkspacePlan(directories=("logs",))

    result = WorkspaceMaterializer().materialize(plan, base_dir)

    assert result.is_success is True


def test_workspace_materializer_refuses_wrong_type_for_directory(tmp_path):
    from encode_pipeline.platform.adapters import WorkspacePlan
    from encode_pipeline.services.materialization import WorkspaceMaterializer

    base_dir = tmp_path.resolve()
    (base_dir / "logs").write_text("not a directory", encoding="utf-8")
    plan = WorkspacePlan(directories=("logs",))

    result = WorkspaceMaterializer().materialize(plan, base_dir)

    assert result.is_failure is True
    issue = result.issues[0]
    assert issue.code == "WORKSPACE_MATERIALIZATION_WRONG_TYPE"
    assert issue.path == "workspace_plan.directories[0]"


def test_workspace_materializer_refuses_symlink_target_directory(tmp_path):
    from encode_pipeline.platform.adapters import WorkspacePlan
    from encode_pipeline.services.materialization import WorkspaceMaterializer

    base_dir = tmp_path.resolve()
    real = base_dir / "real"
    real.mkdir()
    link = base_dir / "link"
    link.symlink_to(real)
    plan = WorkspacePlan(directories=("link",))

    result = WorkspaceMaterializer().materialize(plan, base_dir)

    assert result.is_failure is True
    issue = result.issues[0]
    assert issue.code == "WORKSPACE_MATERIALIZATION_SYMLINK"
    assert issue.path == "workspace_plan.directories[0]"


def test_workspace_materializer_refuses_symlink_parent_directory_under_base_dir(tmp_path):
    from encode_pipeline.platform.adapters import WorkspacePlan
    from encode_pipeline.services.materialization import WorkspaceMaterializer

    base_dir = tmp_path.resolve()
    real = base_dir / "real"
    real.mkdir()
    link = base_dir / "link"
    link.symlink_to(real)
    plan = WorkspacePlan(directories=("link/nested",))

    result = WorkspaceMaterializer().materialize(plan, base_dir)

    assert result.is_failure is True
    issue = result.issues[0]
    assert issue.code == "WORKSPACE_MATERIALIZATION_SYMLINK"
    assert issue.path == "workspace_plan.directories[0]"
    assert issue.message == "Planned directory path has a symlinked parent under base_dir."


def test_workspace_materializer_path_policy_errors_use_safe_messages(tmp_path):
    from encode_pipeline.platform.adapters import WorkspacePlan
    from encode_pipeline.services.materialization import WorkspaceMaterializer

    base_dir = tmp_path.resolve()

    materializer = WorkspaceMaterializer()

    absolute_result = materializer.materialize(
        WorkspacePlan(directories=("/absolute/path",)), base_dir
    )
    assert absolute_result.is_failure is True
    absolute_issue = absolute_result.issues[0]
    assert absolute_issue.code == "WORKSPACE_MATERIALIZATION_PATH_POLICY_PATH_ABSOLUTE"
    assert absolute_issue.message == "Planned path must not be absolute."
    assert "/absolute/path" not in absolute_issue.message

    traversal_result = materializer.materialize(
        WorkspacePlan(directories=("../escape",)), base_dir
    )
    assert traversal_result.is_failure is True
    traversal_issue = traversal_result.issues[0]
    assert traversal_issue.code == "WORKSPACE_MATERIALIZATION_PATH_POLICY_PATH_TRAVERSAL"
    assert traversal_issue.message == "Planned path must not contain '.' or '..' components."
    assert "../escape" not in traversal_issue.message

    invalid_result = materializer.materialize(
        WorkspacePlan(directories=("has space",)), base_dir
    )
    assert invalid_result.is_failure is True
    invalid_issue = invalid_result.issues[0]
    assert invalid_issue.code == "WORKSPACE_MATERIALIZATION_PATH_POLICY_PATH_INVALID"
    assert invalid_issue.message == "Planned path is invalid."
    assert "has space" not in invalid_issue.message


def test_workspace_materializer_creates_files(tmp_path):
    from encode_pipeline.platform.adapters import WorkspacePlan
    from encode_pipeline.services.materialization import WorkspaceMaterializer

    base_dir = tmp_path.resolve()
    plan = WorkspacePlan(
        directories=("logs",),
        files=(("logs/run.log", b"started\n"),),
    )

    result = WorkspaceMaterializer().materialize(plan, base_dir)

    assert result.is_success is True
    assert (base_dir / "logs" / "run.log").read_bytes() == b"started\n"


def test_workspace_materializer_creates_parent_directories_for_files(tmp_path):
    from encode_pipeline.platform.adapters import WorkspacePlan
    from encode_pipeline.services.materialization import WorkspaceMaterializer

    base_dir = tmp_path.resolve()
    plan = WorkspacePlan(files=(("results/peaks/out.narrowPeak", b""),))

    result = WorkspaceMaterializer().materialize(plan, base_dir)

    assert result.is_success is True
    assert (base_dir / "results" / "peaks" / "out.narrowPeak").is_file()


def test_workspace_materializer_refuses_existing_file(tmp_path):
    from encode_pipeline.platform.adapters import WorkspacePlan
    from encode_pipeline.services.materialization import WorkspaceMaterializer

    base_dir = tmp_path.resolve()
    (base_dir / "existing.txt").write_text("old", encoding="utf-8")
    plan = WorkspacePlan(files=(("existing.txt", b"new"),))

    result = WorkspaceMaterializer().materialize(plan, base_dir)

    assert result.is_failure is True
    issue = result.issues[0]
    assert issue.code == "WORKSPACE_MATERIALIZATION_ALREADY_EXISTS"
    assert issue.path == "workspace_plan.files[0]"
    assert (base_dir / "existing.txt").read_text(encoding="utf-8") == "old"


def test_workspace_materializer_refuses_wrong_type_for_file(tmp_path):
    from encode_pipeline.platform.adapters import WorkspacePlan
    from encode_pipeline.services.materialization import WorkspaceMaterializer

    base_dir = tmp_path.resolve()
    (base_dir / "is_dir").mkdir()
    plan = WorkspacePlan(files=(("is_dir", b"x"),))

    result = WorkspaceMaterializer().materialize(plan, base_dir)

    assert result.is_failure is True
    issue = result.issues[0]
    assert issue.code == "WORKSPACE_MATERIALIZATION_WRONG_TYPE"
    assert issue.path == "workspace_plan.files[0]"


def test_workspace_materializer_refuses_symlink_target_file(tmp_path):
    from encode_pipeline.platform.adapters import WorkspacePlan
    from encode_pipeline.services.materialization import WorkspaceMaterializer

    base_dir = tmp_path.resolve()
    real = base_dir / "real.txt"
    real.write_text("x", encoding="utf-8")
    link = base_dir / "link.txt"
    link.symlink_to(real)
    plan = WorkspacePlan(files=(("link.txt", b"y"),))

    result = WorkspaceMaterializer().materialize(plan, base_dir)

    assert result.is_failure is True
    issue = result.issues[0]
    assert issue.code == "WORKSPACE_MATERIALIZATION_SYMLINK"
    assert issue.path == "workspace_plan.files[0]"


def test_workspace_materializer_refuses_symlink_parent_under_base_dir(tmp_path):
    from encode_pipeline.platform.adapters import WorkspacePlan
    from encode_pipeline.services.materialization import WorkspaceMaterializer

    base_dir = tmp_path.resolve()
    real = base_dir / "real"
    real.mkdir()
    link = base_dir / "link"
    link.symlink_to(real)
    plan = WorkspacePlan(files=(("link/file.txt", b"x"),))

    result = WorkspaceMaterializer().materialize(plan, base_dir)

    assert result.is_failure is True
    issue = result.issues[0]
    assert issue.code == "WORKSPACE_MATERIALIZATION_SYMLINK"
    assert issue.path == "workspace_plan.files[0]"


def test_workspace_materializer_refuses_plan_path_traversal(tmp_path):
    from encode_pipeline.platform.adapters import WorkspacePlan
    from encode_pipeline.services.materialization import WorkspaceMaterializer

    base_dir = tmp_path.resolve()
    plan = WorkspacePlan(directories=("../escape",))

    result = WorkspaceMaterializer().materialize(plan, base_dir)

    assert result.is_failure is True
    issue = result.issues[0]
    assert issue.code == "WORKSPACE_MATERIALIZATION_PATH_POLICY_PATH_TRAVERSAL"
    assert issue.path == "workspace_plan.directories[0]"


def test_workspace_materializer_refuses_plan_path_absolute(tmp_path):
    from encode_pipeline.platform.adapters import WorkspacePlan
    from encode_pipeline.services.materialization import WorkspaceMaterializer

    base_dir = tmp_path.resolve()
    plan = WorkspacePlan(files=(("/etc/passwd", b"x"),))

    result = WorkspaceMaterializer().materialize(plan, base_dir)

    assert result.is_failure is True
    issue = result.issues[0]
    assert issue.code == "WORKSPACE_MATERIALIZATION_PATH_POLICY_PATH_ABSOLUTE"
    assert issue.path == "workspace_plan.files[0]"


def test_workspace_materializer_preflight_blocks_all_writes_on_any_failure(tmp_path):
    from encode_pipeline.platform.adapters import WorkspacePlan
    from encode_pipeline.services.materialization import WorkspaceMaterializer

    base_dir = tmp_path.resolve()
    plan = WorkspacePlan(
        directories=("logs",),
        files=(("../escape", b"x"),),
    )

    result = WorkspaceMaterializer().materialize(plan, base_dir)

    assert result.is_failure is True
    assert not (base_dir / "logs").exists()


def test_workspace_materializer_import_boundary() -> None:
    import ast
    from pathlib import Path

    source_path = (
        Path(__file__).resolve().parents[2]
        / "src/encode_pipeline/services/materialization.py"
    )
    source = source_path.read_text(encoding="utf-8")
    tree = ast.parse(source)

    forbidden_modules = {
        "encode_pipeline.api",
        "encode_pipeline.frontend",
        "snakemake",
        "subprocess",
    }
    imported_modules: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported_modules.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module is not None:
            imported_modules.add(node.module)

    assert not any(
        module == forbidden or module.startswith(f"{forbidden}.")
        for module in imported_modules
        for forbidden in forbidden_modules
    )
