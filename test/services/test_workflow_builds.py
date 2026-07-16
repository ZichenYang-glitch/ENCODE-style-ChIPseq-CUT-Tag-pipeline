"""Tests for content-addressed workflow build identities."""

from __future__ import annotations

from datetime import datetime, timezone
import os
from pathlib import Path

from encode_pipeline.platform.adapters import (
    CommandSpec,
    WorkflowCapabilities,
    WorkflowMetadata,
)
from encode_pipeline.platform.builds import WorkflowBuildIdentity
from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.platform.results import Issue, Result
from encode_pipeline.services.defaults import create_default_workflow_registry
from encode_pipeline.services.workflow_builds import WorkflowBuildIdentityProvider


WORKFLOW_ID = "encode-style-chipseq-cuttag-atac-mnase"


class _IdentityAdapter:
    metadata = WorkflowMetadata(
        workflow_id="identity-adapter",
        name="Identity adapter",
        version="2.0.0",
        engines=("opaque-engine",),
    )
    capabilities = WorkflowCapabilities(supports=("validation",))

    def __init__(self, callback):
        self._callback = callback

    def schema(self): ...
    def validate(self, inputs): ...
    def preview_dag(self, inputs): ...
    def plan_workspace(self, inputs, workspace): ...

    def build_command(self, plan, workspace):
        return Result.success(CommandSpec(argv=("unused",)))

    def extract_artifacts(self, inputs, workspace): ...

    def capture_build_identity(self):
        return self._callback()


def _adapter_identity(**overrides) -> WorkflowBuildIdentity:
    values = {
        "workflow_id": "identity-adapter",
        "adapter_version": "2.0.0",
        "scheme": "adapter-lock-v1",
        "logical_entrypoint": "main.nf",
        "digest": "a" * 64,
        "captured_at": datetime.now(timezone.utc),
    }
    values.update(overrides)
    return WorkflowBuildIdentity(**values)


def _project(root: Path, *, marker: str = "same") -> Path:
    files = {
        "pyproject.toml": f"[project]\nname = {marker!r}\n",
        "docs/architecture/artifact-inventory.yaml": "artifacts: []\n",
        "src/encode_pipeline/runtime.py": f"VALUE = {marker!r}\n",
        "src/encode_pipeline/adapters/encode_qc.py": f"CATALOG = {marker!r}\n",
        "workflow/Snakefile": f"# {marker}\nrule all:\n    input: []\n",
        "workflow/rules/main.smk": f"# {marker}\n",
        "workflow/schemas/config.schema.json": '{"type":"object"}\n',
        "workflow/envs/tool.lock": f"# {marker}\n",
        "profiles/default/config.yaml": "cores: 1\n",
        "scripts/tool.py": f"VALUE = {marker!r}\n",
        "scripts/tool.sh": f"#!/bin/sh\n# {marker}\n",
    }
    for relative_path, contents in files.items():
        path = root / relative_path
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(contents, encoding="utf-8")
    return root


def _capture(root: Path):
    return WorkflowBuildIdentityProvider(
        create_default_workflow_registry(),
        project_root=root,
    ).capture(WORKFLOW_ID)


def test_build_identity_delegates_to_adapter_without_reading_project_root(tmp_path):
    identity = _adapter_identity()
    adapter = _IdentityAdapter(lambda: Result.success(identity))
    provider = WorkflowBuildIdentityProvider(
        WorkflowRegistry([adapter]),
        project_root=(tmp_path / "absent").resolve(),
    )

    result = provider.capture(adapter.metadata.workflow_id)

    assert result.is_success
    assert result.value is identity


def test_build_identity_rejects_mismatched_adapter_identity(tmp_path):
    adapter = _IdentityAdapter(
        lambda: Result.success(_adapter_identity(workflow_id="different"))
    )
    result = WorkflowBuildIdentityProvider(
        WorkflowRegistry([adapter]),
        project_root=tmp_path.resolve(),
    ).capture(adapter.metadata.workflow_id)

    assert result.is_failure
    assert result.issues[0].code == "WORKFLOW_BUILD_SOURCE_UNAVAILABLE"


def test_build_identity_sanitizes_adapter_failure(tmp_path):
    secret = str(tmp_path / "private" / "asset-lock")
    adapter = _IdentityAdapter(
        lambda: Result.failure(
            [
                Issue(
                    code="PRIVATE_IDENTITY_FAILURE",
                    message=secret,
                    technical_message=secret,
                )
            ]
        )
    )
    result = WorkflowBuildIdentityProvider(
        WorkflowRegistry([adapter]),
        project_root=tmp_path.resolve(),
    ).capture(adapter.metadata.workflow_id)

    assert result.is_failure
    assert [issue.code for issue in result.issues] == [
        "WORKFLOW_BUILD_SOURCE_UNAVAILABLE"
    ]
    assert secret not in str(result.issues[0].to_dict())


def test_build_identity_sanitizes_adapter_exception(tmp_path):
    secret = str(tmp_path / "private" / "asset-lock")

    def fail():
        raise RuntimeError(secret)

    adapter = _IdentityAdapter(fail)
    result = WorkflowBuildIdentityProvider(
        WorkflowRegistry([adapter]),
        project_root=tmp_path.resolve(),
    ).capture(adapter.metadata.workflow_id)

    assert result.is_failure
    assert [issue.code for issue in result.issues] == [
        "WORKFLOW_BUILD_SOURCE_UNAVAILABLE"
    ]
    assert secret not in str(result.issues[0].to_dict())


def test_encode_build_digest_remains_byte_for_byte_compatible(tmp_path):
    result = _capture(_project(tmp_path / "project"))

    assert result.is_success
    assert result.value.digest == (
        "e51ab94092d85f50baf11ec67056b034b98c05e101fcf0c9c30bd0d3bfdcbd07"
    )


def test_build_digest_is_stable_across_absolute_roots_and_mtimes(tmp_path):
    first_root = _project(tmp_path / "first")
    second_root = _project(tmp_path / "second")
    os.utime(second_root / "workflow" / "Snakefile", (1, 1))

    first = _capture(first_root)
    second = _capture(second_root)

    assert first.is_success
    assert second.is_success
    assert first.value.matches(second.value)


def test_build_digest_changes_when_controlled_source_changes(tmp_path):
    root = _project(tmp_path / "project")
    before = _capture(root)
    (root / "scripts" / "tool.py").write_text("VALUE = 'changed'\n", encoding="utf-8")
    after = _capture(root)

    assert before.is_success
    assert after.is_success
    assert not before.value.matches(after.value)


def test_build_digest_changes_when_artifact_inventory_changes(tmp_path):
    root = _project(tmp_path / "project")
    before = _capture(root)
    (root / "docs/architecture/artifact-inventory.yaml").write_text(
        "artifacts:\n- id: changed\n",
        encoding="utf-8",
    )
    after = _capture(root)

    assert before.is_success
    assert after.is_success
    assert not before.value.matches(after.value)


def test_build_digest_changes_when_qc_parser_catalog_changes(tmp_path):
    root = _project(tmp_path / "project")
    before = _capture(root)
    (root / "src/encode_pipeline/adapters/encode_qc.py").write_text(
        "CATALOG = 'changed'\n",
        encoding="utf-8",
    )
    after = _capture(root)

    assert before.is_success
    assert after.is_success
    assert not before.value.matches(after.value)


def test_build_digest_covers_workflow_files_without_known_suffixes(tmp_path):
    root = _project(tmp_path / "project")
    source = root / "workflow" / "scripts" / "analysis.R"
    source.parent.mkdir()
    source.write_text("VALUE <- 1\n", encoding="utf-8")
    before = _capture(root)
    source.write_text("VALUE <- 2\n", encoding="utf-8")
    after = _capture(root)

    assert before.is_success
    assert after.is_success
    assert not before.value.matches(after.value)


def test_build_digest_ignores_cache_and_pyc_files(tmp_path):
    root = _project(tmp_path / "project")
    before = _capture(root)
    cache = root / "src" / "encode_pipeline" / "__pycache__"
    cache.mkdir()
    (cache / "runtime.cpython-313.pyc").write_bytes(b"generated")
    after = _capture(root)

    assert before.is_success
    assert after.is_success
    assert before.value.matches(after.value)


def test_build_digest_includes_all_default_profile_files(tmp_path):
    root = _project(tmp_path / "project")
    profile_marker = root / "profiles" / "default" / "runtime-flags"
    profile_marker.write_text("first\n", encoding="utf-8")
    before = _capture(root)
    profile_marker.write_text("second\n", encoding="utf-8")
    after = _capture(root)

    assert before.is_success
    assert after.is_success
    assert not before.value.matches(after.value)


def test_build_digest_uses_only_top_level_execution_scripts(tmp_path):
    root = _project(tmp_path / "project")
    nested_script = root / "scripts" / "internal" / "helper.py"
    nested_script.parent.mkdir()
    nested_script.write_text("VALUE = 'first'\n", encoding="utf-8")
    before = _capture(root)
    nested_script.write_text("VALUE = 'second'\n", encoding="utf-8")
    after = _capture(root)

    assert before.is_success
    assert after.is_success
    assert before.value.matches(after.value)


def test_build_capture_rejects_missing_required_source(tmp_path):
    root = _project(tmp_path / "project")
    (root / "workflow" / "Snakefile").unlink()

    result = _capture(root)

    assert result.is_failure
    assert result.issues[0].code == "WORKFLOW_BUILD_SOURCE_UNAVAILABLE"
    assert str(root) not in str(result.issues[0].to_dict())


def test_build_capture_rejects_symlinked_source(tmp_path):
    root = _project(tmp_path / "project")
    target = root / "scripts" / "tool.py"
    link = root / "scripts" / "linked.py"
    link.symlink_to(target)

    result = _capture(root)

    assert result.is_failure
    assert result.issues[0].code == "WORKFLOW_BUILD_SOURCE_UNAVAILABLE"
