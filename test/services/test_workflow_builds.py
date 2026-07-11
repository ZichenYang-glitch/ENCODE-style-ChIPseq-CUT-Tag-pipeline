"""Tests for content-addressed workflow build identities."""

from __future__ import annotations

import os
from pathlib import Path

from encode_pipeline.services.defaults import create_default_workflow_registry
from encode_pipeline.services.workflow_builds import WorkflowBuildIdentityProvider


WORKFLOW_ID = "encode-style-chipseq-cuttag-atac-mnase"


def _project(root: Path, *, marker: str = "same") -> Path:
    files = {
        "pyproject.toml": f"[project]\nname = {marker!r}\n",
        "src/encode_pipeline/runtime.py": f"VALUE = {marker!r}\n",
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
