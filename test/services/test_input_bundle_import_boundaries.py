"""Static and process boundaries for Omics Intake Bundle consumption."""

from __future__ import annotations

import ast
import os
from pathlib import Path
import subprocess
import sys
import textwrap
import tomllib


REPO_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = REPO_ROOT / "src"
BOUNDARY_MODULES = (
    SRC_ROOT / "encode_pipeline/platform/input_bundles.py",
    SRC_ROOT / "encode_pipeline/services/input_bundle_imports.py",
)
ENCODE_ADAPTER_MODULE = SRC_ROOT / "encode_pipeline/adapters/encode.py"


def _imported_names(path: Path) -> set[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"))
    names: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            names.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module is not None:
            names.add(node.module)
    return names


def _dotted_name(node: ast.expr) -> str:
    if isinstance(node, ast.Name):
        return node.id
    if isinstance(node, ast.Attribute):
        parent = _dotted_name(node.value)
        return f"{parent}.{node.attr}" if parent else node.attr
    return ""


def test_boundary_modules_do_not_import_producer_or_execution_internals() -> None:
    forbidden = {
        "omics_intake",
        "subprocess",
        "snakemake",
        "encode_pipeline.adapters",
        "encode_pipeline.config",
        "encode_pipeline.samples",
        "encode_pipeline.persistence",
        "encode_pipeline.workers",
    }

    for path in BOUNDARY_MODULES:
        imports = _imported_names(path)
        violations = sorted(
            name
            for name in imports
            if any(name == item or name.startswith(item + ".") for item in forbidden)
        )
        assert not violations, f"{path.name} imports forbidden modules: {violations}"


def test_generic_service_contains_no_encode_private_authoring_fields() -> None:
    source = (BOUNDARY_MODULES[1]).read_text(encoding="utf-8")

    for token in ("fastq", "bowtie2_index", "control_bam", "peak_mode"):
        assert token not in source


def test_project_declares_schema_validator_but_no_producer_dependency() -> None:
    project = tomllib.loads((REPO_ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    dependencies = project["project"]["dependencies"]
    normalized = tuple(item.lower().replace("_", "-") for item in dependencies)

    assert any(item.startswith("jsonschema") for item in normalized)
    assert not any("omics-intake" in item for item in normalized)


def test_importing_boundary_does_not_load_producer_or_workflow_modules() -> None:
    code = """
        import sys
        import encode_pipeline.services.input_bundle_imports

        forbidden = (
            "omics_intake",
            "encode_pipeline.adapters.encode",
            "encode_pipeline.config.validator",
            "encode_pipeline.samples",
            "snakemake",
        )
        for name in forbidden:
            loaded = name in sys.modules or any(
                module.startswith(name + ".") for module in sys.modules
            )
            print(f"{name}={loaded}")
    """
    environment = dict(os.environ)
    environment["PYTHONPATH"] = str(SRC_ROOT)
    environment["PYTHONDONTWRITEBYTECODE"] = "1"

    completed = subprocess.run(
        [sys.executable, "-c", textwrap.dedent(code)],
        capture_output=True,
        text=True,
        env=environment,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr
    assert set(completed.stdout.splitlines()) == {
        "omics_intake=False",
        "encode_pipeline.adapters.encode=False",
        "encode_pipeline.config.validator=False",
        "encode_pipeline.samples=False",
        "snakemake=False",
    }


def test_encode_bundle_mapping_reachable_graph_has_no_filesystem_probes() -> None:
    tree = ast.parse(ENCODE_ADAPTER_MODULE.read_text(encoding="utf-8"))
    functions = {
        node.name: node
        for node in tree.body
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
    }
    adapter_class = next(
        node
        for node in tree.body
        if isinstance(node, ast.ClassDef) and node.name == "EncodeStyleWorkflowAdapter"
    )
    methods = {
        node.name: node
        for node in adapter_class.body
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
    }
    pending = [methods["import_input_bundle"], functions["_map_input_bundle"]]
    visited: set[str] = set()
    violations: set[str] = set()
    forbidden_attributes = {
        "exists",
        "glob",
        "is_file",
        "isfile",
        "iterdir",
        "lstat",
        "open",
        "read_bytes",
        "read_text",
        "resolve",
        "rglob",
        "stat",
    }
    forbidden_bare_calls = {"open", "stat", "lstat"}

    while pending:
        function = pending.pop()
        if function.name in visited:
            continue
        visited.add(function.name)
        for node in ast.walk(function):
            if isinstance(node, (ast.Import, ast.ImportFrom)):
                names = (
                    [alias.name for alias in node.names]
                    if isinstance(node, ast.Import)
                    else [node.module or ""]
                )
                if any(name == "os" or name == "stat" for name in names):
                    violations.add(f"{function.name}:import")
            if not isinstance(node, ast.Call):
                continue
            dotted = _dotted_name(node.func)
            if isinstance(node.func, ast.Name):
                if node.func.id in forbidden_bare_calls:
                    violations.add(f"{function.name}:{dotted}")
                helper = functions.get(node.func.id)
                if helper is not None:
                    pending.append(helper)
            elif isinstance(node.func, ast.Attribute):
                if (
                    node.func.attr in forbidden_attributes
                    or dotted == "os.path"
                    or dotted.startswith("os.path.")
                    or dotted in {"os.open", "os.stat", "os.lstat"}
                ):
                    violations.add(f"{function.name}:{dotted}")
                if (
                    isinstance(node.func.value, ast.Name)
                    and node.func.value.id == "self"
                ):
                    helper = methods.get(node.func.attr)
                    if helper is not None:
                        pending.append(helper)

    assert "_map_bowtie2_index" in visited
    assert not violations, sorted(violations)


def test_full_bundle_inspection_does_not_load_producer_or_execution_stacks() -> None:
    code = """
        import sys
        from pathlib import Path
        from tempfile import TemporaryDirectory

        from encode_pipeline.adapters.encode import EncodeStyleWorkflowAdapter
        from encode_pipeline.platform.registry import WorkflowRegistry
        from encode_pipeline.services.input_bundle_imports import InputBundleImportService
        from input_bundle_support import create_runnable_bundle

        with TemporaryDirectory() as directory:
            bundle_path = create_runnable_bundle(Path(directory) / "bundle")
            adapter = EncodeStyleWorkflowAdapter()
            result = InputBundleImportService(
                WorkflowRegistry((adapter,))
            ).inspect(bundle_path, adapter.metadata.workflow_id)
            if result.is_failure:
                raise SystemExit(result.errors[0].code)

        forbidden = (
            "omics_intake",
            "snakemake",
            "encode_pipeline.workers",
            "encode_pipeline.persistence",
        )
        for prefix in forbidden:
            loaded = any(
                name == prefix or name.startswith(prefix + ".")
                for name in sys.modules
            )
            print(f"{prefix}={loaded}")
    """
    environment = dict(os.environ)
    environment["PYTHONPATH"] = os.pathsep.join(
        (str(SRC_ROOT), str(REPO_ROOT / "test"))
    )
    environment["PYTHONDONTWRITEBYTECODE"] = "1"

    completed = subprocess.run(
        [sys.executable, "-c", textwrap.dedent(code)],
        capture_output=True,
        text=True,
        env=environment,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr
    assert set(completed.stdout.splitlines()) == {
        "omics_intake=False",
        "snakemake=False",
        "encode_pipeline.workers=False",
        "encode_pipeline.persistence=False",
    }
