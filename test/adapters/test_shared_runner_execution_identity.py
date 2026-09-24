"""The real default runner's execution authority must be identity controlled.

Candidate modules are loaded only to inspect the allowlist; no candidate program,
scientific process, Docker endpoint or Redis instance is executed.
"""

from __future__ import annotations

import ast
import hashlib
import importlib.util
from pathlib import Path
import shutil
import sys

import pytest

from encode_pipeline.adapters.bulk_rnaseq.adapter import BulkRnaSeqWorkflowAdapter
from encode_pipeline.adapters.bulk_rnaseq.execution_identity import (
    EXECUTION_IMPLEMENTATION_PATHS,
    build_execution_implementation_manifest,
    canonical_execution_manifest_bytes,
    verify_execution_implementation,
)
from encode_pipeline.adapters.hitrac_preprocess import deployment, execution
from encode_pipeline.adapters.hitrac_preprocess.adapter import HiTracPreprocessAdapter
from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.services.defaults import create_default_process_runner
from encode_pipeline.workers.settings import WorkerSettings


ROOT = Path(__file__).resolve().parents[2]
CONTROL_MODULES = (
    "__init__",
    "deployment",
    "adapter",
    "execution",
    "admission",
    "qualification",
    "calls",
)
DEPLOYMENT = "src/encode_pipeline/adapters/hitrac_preprocess/deployment.py"


def _controlled_copy(destination):
    for relative in EXECUTION_IMPLEMENTATION_PATHS:
        target = destination / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / relative, target)
    return destination


def _runner(registry, tmp_path):
    settings = WorkerSettings(
        database_url=f"sqlite:///{tmp_path / 'unused.db'}",
        redis_url="redis://unused.invalid/0",
        queue_name="identity-contract",
        workspace_root=tmp_path / "unused-workspaces",
    )
    return create_default_process_runner(registry=registry, settings=settings)


def test_changed_runner_authority_cannot_keep_the_same_bulk_identity(
    tmp_path, monkeypatch
):
    registry = WorkflowRegistry(adapters=[BulkRnaSeqWorkflowAdapter()])
    original_allowlist = _runner(registry, tmp_path)._allowed_executables
    copied = _controlled_copy(tmp_path / "implementation")
    before = build_execution_implementation_manifest(copied)

    original = Path(deployment.__file__).read_text()
    function = next(
        item
        for item in ast.parse(original).body
        if isinstance(item, ast.FunctionDef)
        and item.name == "local_execution_executable"
    )
    unexecuted = str(tmp_path / "unreviewed-executable-never-created-or-run")
    lines = original.splitlines(keepends=True)
    candidate = copied / DEPLOYMENT
    candidate.parent.mkdir(parents=True, exist_ok=True)
    candidate.write_text(
        "".join(lines[: function.lineno - 1])
        + f"def local_execution_executable(adapter):\n    return {unexecuted!r}\n"
        + "".join(lines[function.end_lineno :])
    )
    spec = importlib.util.spec_from_file_location(deployment.__name__, candidate)
    altered = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(altered)
    monkeypatch.setitem(sys.modules, deployment.__name__, altered)
    changed_allowlist = _runner(registry, tmp_path)._allowed_executables
    assert original_allowlist == ("snakemake",)
    assert changed_allowlist == ("snakemake", unexecuted)
    assert not Path(unexecuted).exists()

    after = build_execution_implementation_manifest(copied)
    # Same real factory and Bulk-only registry, different authority: old identity
    # must not continue to describe this executable policy.
    assert before["aggregate_sha256"] != after["aggregate_sha256"]
    verified = verify_execution_implementation(
        manifest_bytes=canonical_execution_manifest_bytes(before),
        package_root=copied / "src/encode_pipeline",
    )
    assert verified.is_failure
    assert verified.issues[0].code == "BULK_RNASEQ_EXECUTION_IMPLEMENTATION_INVALID"


@pytest.mark.parametrize("module", CONTROL_MODULES)
def test_every_hitrac_runner_control_dependency_is_identity_bound(tmp_path, module):
    relative = f"src/encode_pipeline/adapters/hitrac_preprocess/{module}.py"
    assert relative in EXECUTION_IMPLEMENTATION_PATHS
    copied = _controlled_copy(tmp_path / "implementation")
    before = build_execution_implementation_manifest(copied)
    path = copied / relative
    path.write_bytes(path.read_bytes() + b"\n# controlled byte-drift probe\n")
    after = build_execution_implementation_manifest(copied)
    assert after["aggregate_sha256"] != before["aggregate_sha256"]
    verified = verify_execution_implementation(
        manifest_bytes=canonical_execution_manifest_bytes(before),
        package_root=copied / "src/encode_pipeline",
    )
    assert verified.is_failure


def test_runtime_lock_drift_cannot_change_shared_runner_admission(
    tmp_path, monkeypatch
):
    # Only the expensive scientific runtime verification is replaced. This test
    # isolates the shared execution authority's fixed runtime-policy input.
    monkeypatch.setattr(execution.RuntimeAdmission, "verify", lambda self: None)
    runtime = execution.RuntimeAdmission(
        binding=tmp_path / "unused-binding.json",
        binding_sha256="a" * 64,
        python=Path(sys.executable).absolute(),
        python_sha256=hashlib.sha256(Path(sys.executable).read_bytes()).hexdigest(),
        implementation_sha256="b" * 64,
    )
    adapter = HiTracPreprocessAdapter(runtime=runtime)
    registry = WorkflowRegistry(adapters=[adapter])
    copied = tmp_path / "implementation"
    lock = copied / "config/hitrac_preprocess/tools.lock.json"
    lock.parent.mkdir(parents=True)
    lock.write_bytes((ROOT / "config/hitrac_preprocess/tools.lock.json").read_bytes())
    monkeypatch.setattr(deployment, "__file__", str(copied / DEPLOYMENT))
    assert _runner(registry, tmp_path)._allowed_executables == (
        "snakemake",
        str(runtime.python),
    )
    lock.write_bytes(lock.read_bytes() + b"\n")
    assert _runner(registry, tmp_path)._allowed_executables == ("snakemake",)
    assert adapter.execution_availability().execution == "not_configured"
