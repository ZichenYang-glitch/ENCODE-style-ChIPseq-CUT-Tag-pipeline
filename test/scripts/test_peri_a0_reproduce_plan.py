"""Relocation preflight for the frozen Peri A0 manifest."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import shutil
import subprocess

import pytest

_REPO = Path(__file__).resolve().parents[2]
_SCRIPT = _REPO / "scripts/peri_a0_reproduce/plan.py"
_SPEC = importlib.util.spec_from_file_location("peri_a0_reproduce_plan", _SCRIPT)
plan = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(plan)

_REAL_MANIFEST = _REPO / "scripts/peri_a0_closeout/evidence/build-testfix.json"
_REAL_PATCH = _REPO / "scripts/peri_a0_closeout/evidence/candidate-testfix.patch"
_OLD_ROOTS = (
    "/tmp/helixweave-peri-a0-20261109",
    "/tmp/hw-peri-a0-close-20261109",
)
_ACCEPTED_PATCH = "06272e1a66423a9e562509a38d044cd8506a1c10e61ef1fd3fe239107680587b"
_ACCEPTED_LOCK = "0692683359fb75073eca3d14b364071bc02058be7b11e0335ccaf2b0f314f5dd"
_ACCEPTED_COMMIT = "c7d154c9985cbe08ba25e964cb7468d15014e6d0"


def _executable(path: Path, version_line: str) -> None:
    path.write_text(f"#!/bin/sh\necho '{version_line}'\n", encoding="utf-8")
    path.chmod(0o755)


def _git(source: Path, *args: str) -> None:
    subprocess.run(
        [
            "git",
            "-C",
            str(source),
            "-c",
            "user.email=peri-a0-plan@example.com",
            "-c",
            "user.name=Peri A0 Plan Test",
            "-c",
            "commit.gpgsign=false",
            *args,
        ],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )


def _source_tree(path: Path, lock_text: str = "fixture-lock\n") -> None:
    path.mkdir(parents=True)
    (path / "Cargo.lock").write_text(lock_text, encoding="utf-8")
    (path / "rust-toolchain.toml").write_text(
        '[toolchain]\nchannel = "1.99.0"\n',
        encoding="utf-8",
    )
    script = path / "scripts"
    script.mkdir()
    (script / "cargo-rmcp-patched.sh").write_text("#!/bin/sh\n", encoding="utf-8")
    _git(path, "init", "-q")
    _git(
        path,
        "add",
        "Cargo.lock",
        "rust-toolchain.toml",
        "scripts/cargo-rmcp-patched.sh",
    )
    _git(path, "commit", "-q", "-m", "fixture")


def _toolchain(path: Path, rustc_line: str = "rustc 1.99.0 (fixture)") -> None:
    bindir = path / "bin"
    bindir.mkdir(parents=True, exist_ok=True)
    _executable(bindir / "rustc", rustc_line)
    _executable(bindir / "cargo", "cargo 1.99.0 (fixture)")


def _cache(path: Path) -> None:
    (path / "registry").mkdir(parents=True)


def _manifest(path: Path, environment: dict, command: list[str] | None = None) -> None:
    payload = {
        "source_commit": plan.PERI_COMMIT,
        "patch_sha256": plan.PATCH_SHA256,
        "cargo_lock_before": plan.CARGO_LOCK_SHA256,
        "environment": environment,
        "command": command
        or ["/tmp/peri-a0-historical-moved/scripts/cargo-rmcp-patched.sh", "build"],
        "binary": "/tmp/peri-a0-historical-moved/build/target/release/peri",
        "cache_copy_from": "/tmp/peri-a0-historical-moved",
        "toolchain_reused_read_only": "/tmp/peri-a0-historical-moved/toolchain",
    }
    path.write_text(json.dumps(payload), encoding="utf-8")


def _old_environment() -> dict[str, str]:
    root = "/tmp/peri-a0-historical-moved"
    return {
        "PATH": f"{root}/toolchain/bin:/usr/bin:/bin",
        "HOME": f"{root}/build-home",
        "CARGO_HOME": f"{root}/cargo-home",
        "CARGO_TARGET_DIR": f"{root}/build/target",
        "CARGO_BUILD_JOBS": "99",
        "CARGO_NET_OFFLINE": "false",
        "LANG": "C.UTF-8",
        "GIT_CONFIG_GLOBAL": "/dev/null",
        "GIT_CONFIG_NOSYSTEM": "1",
        "CARGO_NET_GIT_FETCH_WITH_CLI": "true",
        "GIT_CONFIG_COUNT": "1",
        "GIT_CONFIG_KEY_0": "url.ssh://git@github.com/.insteadOf",
        "GIT_CONFIG_VALUE_0": "https://github.com/",
        "GIT_SSH_COMMAND": f"ssh -i {root}/id -o BatchMode=yes",
        "LD_PRELOAD": f"{root}/libpreload.so",
        "RUSTC": f"{root}/toolchain/bin/rustc",
    }


def _layout(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    source = tmp_path / "peri"
    _source_tree(source)
    toolchain = tmp_path / "toolchain"
    _toolchain(toolchain)
    cache = tmp_path / "cache"
    _cache(cache)
    manifest = tmp_path / "historical.json"
    _manifest(manifest, _old_environment())
    task = tmp_path / "task"
    monkeypatch.setattr(plan, "source_commit", lambda _source: plan.PERI_COMMIT)
    monkeypatch.setattr(plan, "lock_sha256", lambda _source: plan.CARGO_LOCK_SHA256)
    return source, toolchain, cache, manifest, task


def _argv(source, toolchain, cache, manifest, task, patch=_REAL_PATCH):
    return [
        "--source",
        str(source),
        "--patch",
        str(patch),
        "--toolchain",
        str(toolchain),
        "--cache",
        str(cache),
        "--task-root",
        str(task),
        "--historical-manifest",
        str(manifest),
    ]


def _assert_relocated(result: dict, toolchain: Path, task: Path, old: str) -> None:
    environment = result["environment"]
    command = result["command"]
    rendered = json.dumps({"command": command, "environment": environment})
    assert old not in rendered
    assert environment["PATH"] == f"{toolchain.resolve() / 'bin'}:/usr/bin:/bin"
    assert environment["HOME"] == str(task.resolve() / "build-home")
    assert environment["CARGO_HOME"] == str(task.resolve() / "cargo-home")
    assert environment["CARGO_TARGET_DIR"] == str(task.resolve() / "build" / "target")
    assert environment["CARGO_HOME"] != str(result["cache_source"])
    assert environment["CARGO_BUILD_JOBS"] == "4"
    assert environment["CARGO_NET_OFFLINE"] == "true"
    assert command == list(plan.BUILD_COMMAND)
    assert result["starts_build"] is False
    assert result["online_locked_supported"] is False
    assert not (task / "cargo-home").exists()
    assert not (task / "build").exists()
    assert not (task / "build-home").exists()


def test_accepted_identity_constants_match_the_archive():
    assert plan.PERI_COMMIT == _ACCEPTED_COMMIT
    assert plan.PATCH_SHA256 == _ACCEPTED_PATCH
    assert plan.CARGO_LOCK_SHA256 == _ACCEPTED_LOCK
    assert hashlib.sha256(_REAL_PATCH.read_bytes()).hexdigest() == _ACCEPTED_PATCH


def test_relocated_inputs_do_not_inherit_historical_paths(tmp_path, monkeypatch):
    source, toolchain, cache, manifest, task = _layout(tmp_path, monkeypatch)
    code = plan.main(_argv(source, toolchain, cache, manifest, task))
    written = json.loads((task / "reproduce-plan.json").read_text(encoding="utf-8"))

    assert code == 0
    _assert_relocated(written, toolchain, task, "/tmp/peri-a0-historical-moved")
    assert "GIT_SSH_COMMAND" in written["dropped_historical_keys"]
    assert "LD_PRELOAD" in written["dropped_historical_keys"]
    assert "RUSTC" in written["dropped_historical_keys"]
    assert "PATH" in written["replaced_historical_keys"]
    assert written["environment"]["GIT_CONFIG_VALUE_0"] == "https://github.com/"
    assert written["patch_application"] == "apply-before-build"
    assert "/tmp/peri-a0-historical-moved" not in json.dumps(written)


def test_real_closeout_manifest_paths_are_replaced(tmp_path, monkeypatch):
    source, toolchain, cache, _manifest_path, task = _layout(tmp_path, monkeypatch)
    historical = json.loads(_REAL_MANIFEST.read_text(encoding="utf-8"))
    assert (
        "/tmp/helixweave-peri-a0-20261109/toolchain/bin"
        in historical["environment"]["PATH"]
    )
    result = plan.run(
        _argv(source, toolchain, cache, _REAL_MANIFEST, task),
    )

    rendered = json.dumps(
        {"command": result["command"], "environment": result["environment"]}
    )
    for root in _OLD_ROOTS:
        assert root not in rendered
    assert result["environment"]["PATH"] != historical["environment"]["PATH"]
    _assert_relocated(result, toolchain, task, _OLD_ROOTS[0])
    assert _OLD_ROOTS[1] not in json.dumps(result["environment"])


def test_online_locked_is_not_a_plan_argument():
    with pytest.raises(SystemExit) as caught:
        plan.main(["--online-locked"])
    assert caught.value.code == 2


def test_original_experiment_directory_is_not_created(tmp_path, monkeypatch):
    source, toolchain, cache, manifest, _task = _layout(tmp_path, monkeypatch)
    forbidden = Path("/tmp/hw-peri-a0-close-20261109/not-created-by-a0-plan")
    assert not forbidden.exists()

    with pytest.raises(plan.PlanRefusal, match="original experiment"):
        plan.run(_argv(source, toolchain, cache, manifest, forbidden))

    assert not forbidden.exists()


def test_shared_environment_and_shared_target_are_refused(tmp_path, monkeypatch):
    source, toolchain, cache, manifest, _task = _layout(tmp_path, monkeypatch)
    shared = tmp_path / "shared-env"
    (shared / "conda-meta").mkdir(parents=True)
    with pytest.raises(plan.PlanRefusal, match="shared environment"):
        plan.run(_argv(source, toolchain, cache, manifest, shared / "task"))

    historical_target = tmp_path / "old-target"
    historical_target.mkdir()
    environment = _old_environment()
    environment["CARGO_TARGET_DIR"] = str(historical_target)
    target_manifest = tmp_path / "target.json"
    _manifest(target_manifest, environment)
    with pytest.raises(plan.PlanRefusal, match="historical local path"):
        plan.run(_argv(source, toolchain, cache, target_manifest, historical_target))
    assert list(historical_target.iterdir()) == []


def test_cache_inside_task_output_is_refused(tmp_path, monkeypatch):
    source, toolchain, _existing_cache, manifest, task = _layout(tmp_path, monkeypatch)
    task.mkdir()
    nested_cache = task / "cache"
    _cache(nested_cache)

    with pytest.raises(plan.PlanRefusal, match="cache source"):
        plan.run(_argv(source, toolchain, nested_cache, manifest, task))

    assert not (task / "reproduce-plan.json").exists()


def test_wrong_commit_patch_lock_and_toolchain_are_refused(tmp_path, monkeypatch):
    source, toolchain, cache, manifest, task = _layout(tmp_path, monkeypatch)
    monkeypatch.undo()

    with pytest.raises(plan.PlanRefusal, match="source commit"):
        plan.run(_argv(source, toolchain, cache, manifest, task))
    assert not task.exists()

    monkeypatch.setattr(plan, "source_commit", lambda _source: plan.PERI_COMMIT)
    bad_patch = tmp_path / "bad.patch"
    bad_patch.write_text("not-the-patch\n", encoding="utf-8")
    with pytest.raises(plan.PlanRefusal, match="patch file"):
        plan.run(_argv(source, toolchain, cache, manifest, task, patch=bad_patch))

    monkeypatch.setattr(plan, "lock_sha256", lambda _source: "0" * 64)
    with pytest.raises(plan.PlanRefusal, match="Cargo.lock"):
        plan.run(_argv(source, toolchain, cache, manifest, task))

    monkeypatch.setattr(plan, "lock_sha256", lambda _source: plan.CARGO_LOCK_SHA256)
    _toolchain(toolchain, rustc_line="rustc 1.98.0 (wrong)")
    with pytest.raises(plan.PlanRefusal, match="rustc is not 1.99.0"):
        plan.run(_argv(source, toolchain, cache, manifest, task))


def test_cache_toolchain_channel_and_dirty_source_are_refused(tmp_path, monkeypatch):
    source, toolchain, cache, manifest, task = _layout(tmp_path, monkeypatch)
    (cache / "registry").rmdir()
    with pytest.raises(plan.PlanRefusal, match="registry nor git"):
        plan.run(_argv(source, toolchain, cache, manifest, task))
    _cache(cache)

    (source / "rust-toolchain.toml").write_text("[toolchain]\n", encoding="utf-8")
    with pytest.raises(plan.PlanRefusal, match="no channel"):
        plan.run(_argv(source, toolchain, cache, manifest, task))
    (source / "rust-toolchain.toml").write_text(
        '[toolchain]\nchannel = "1.98.0"\n',
        encoding="utf-8",
    )
    with pytest.raises(plan.PlanRefusal, match="channel"):
        plan.run(_argv(source, toolchain, cache, manifest, task))
    (source / "rust-toolchain.toml").write_text(
        '[toolchain]\nchannel = "1.99.0"\n',
        encoding="utf-8",
    )

    (source / "Cargo.lock").write_text("changed\n", encoding="utf-8")
    with pytest.raises(plan.PlanRefusal, match="source diff"):
        plan.run(_argv(source, toolchain, cache, manifest, task))
    assert not (task / "reproduce-plan.json").exists()


def test_nonempty_task_and_failed_toolchain_command_are_refused(tmp_path, monkeypatch):
    source, toolchain, cache, manifest, task = _layout(tmp_path, monkeypatch)
    task.mkdir()
    (task / "existing").write_text("keep\n", encoding="utf-8")
    with pytest.raises(plan.PlanRefusal, match="not empty"):
        plan.run(_argv(source, toolchain, cache, manifest, task))
    assert (task / "existing").read_text(encoding="utf-8") == "keep\n"

    rustc = toolchain / "bin" / "rustc"
    rustc.write_text("#!/bin/sh\nexit 1\n", encoding="utf-8")
    rustc.chmod(0o755)
    with pytest.raises(plan.PlanRefusal, match="rustc -V failed"):
        plan.run(_argv(source, toolchain, cache, manifest, tmp_path / "other-task"))
    rustc.chmod(0o644)
    with pytest.raises(plan.PlanRefusal, match="not executable"):
        plan.run(_argv(source, toolchain, cache, manifest, tmp_path / "third-task"))


def test_already_applied_diff_is_recorded(tmp_path, monkeypatch):
    source, toolchain, cache, _manifest_path, task = _layout(tmp_path, monkeypatch)
    task.mkdir()
    nested = task / "peri"
    os.rename(source, nested)
    (nested / "Cargo.lock").write_text("fixture-lock\nchanged\n", encoding="utf-8")
    diff = subprocess.check_output(["git", "-C", str(nested), "diff", "HEAD"])
    digest = hashlib.sha256(diff).hexdigest()
    monkeypatch.setattr(plan, "PATCH_SHA256", digest)
    patch = tmp_path / "applied.patch"
    patch.write_bytes(diff)
    manifest = tmp_path / "applied.json"
    _manifest(manifest, _old_environment())

    result = plan.run(_argv(nested, toolchain, cache, manifest, task, patch=patch))

    assert result["patch_application"] == "already-applied"
    assert result["cwd"] == str(nested.resolve())
    _assert_relocated(result, toolchain, task, "/tmp/peri-a0-historical-moved")


def test_missing_git_repository_is_refused(tmp_path, monkeypatch):
    _source, toolchain, cache, manifest, task = _layout(tmp_path, monkeypatch)
    monkeypatch.undo()
    plain = tmp_path / "plain"
    _source_tree(plain)
    shutil.rmtree(plain / ".git")

    with pytest.raises(plan.PlanRefusal, match="not a git repository"):
        plan.run(_argv(plain, toolchain, cache, manifest, task))


def test_more_refusals_leave_outputs_unwritten(tmp_path, monkeypatch, capsys):
    source, toolchain, cache, manifest, task = _layout(tmp_path, monkeypatch)
    code = plan.main(
        _argv(
            source,
            toolchain,
            cache,
            manifest,
            Path("/tmp/helixweave-peri-design-20261108/not-created-by-a0-plan"),
        )
    )
    assert code == 2
    assert "original experiment" in capsys.readouterr().err

    experiment = tmp_path / "experiment"
    experiment.mkdir()
    placed = experiment / "peri"
    os.rename(source, placed)
    monkeypatch.setattr(plan, "ORIGINAL_EXPERIMENT_ROOTS", (experiment,))
    with pytest.raises(plan.PlanRefusal, match="source is inside an original"):
        plan.run(_argv(placed, toolchain, cache, manifest, task))
    os.rename(placed, source)
    monkeypatch.undo()
    monkeypatch.setattr(plan, "source_commit", lambda _source: plan.PERI_COMMIT)
    monkeypatch.setattr(plan, "lock_sha256", lambda _source: plan.CARGO_LOCK_SHA256)

    historical = _old_environment()
    historical["CARGO_TARGET_DIR"] = str(source)
    overlap = tmp_path / "overlap.json"
    _manifest(overlap, historical)
    with pytest.raises(plan.PlanRefusal, match="source overlaps"):
        plan.run(_argv(source, toolchain, cache, overlap, task))

    task.mkdir()
    nested_toolchain = task / "toolchain"
    os.rename(toolchain, nested_toolchain)
    with pytest.raises(plan.PlanRefusal, match="toolchain must stay outside"):
        plan.run(_argv(source, nested_toolchain, cache, manifest, task))
    os.rename(nested_toolchain, toolchain)
    task.rmdir()

    venv = tmp_path / "venv"
    venv.mkdir()
    (venv / "pyvenv.cfg").write_text("home = fixture\n", encoding="utf-8")
    with pytest.raises(plan.PlanRefusal, match="shared environment"):
        plan.run(_argv(source, toolchain, cache, manifest, venv / "task"))

    (source / "scripts" / "cargo-rmcp-patched.sh").unlink()
    with pytest.raises(plan.PlanRefusal, match="cargo-rmcp-patched.sh"):
        plan.run(_argv(source, toolchain, cache, manifest, task))
    (source / "scripts" / "cargo-rmcp-patched.sh").write_text(
        "#!/bin/sh\n", encoding="utf-8"
    )

    (source / "Cargo.lock").unlink()
    with pytest.raises(plan.PlanRefusal, match="missing Cargo.lock"):
        plan.run(_argv(source, toolchain, cache, manifest, task))
    (source / "Cargo.lock").write_text("fixture-lock\n", encoding="utf-8")

    with pytest.raises(plan.PlanRefusal, match="parent does not exist"):
        plan.run(
            _argv(source, toolchain, cache, manifest, tmp_path / "missing" / "task")
        )

    nested_root = tmp_path / "nested-task"
    (nested_root / "child").mkdir(parents=True)
    nested_source = nested_root / "child" / "peri"
    os.rename(source, nested_source)
    with pytest.raises(plan.PlanRefusal, match="direct child"):
        plan.run(_argv(nested_source, toolchain, cache, manifest, nested_root))
    os.rename(nested_source, source)

    broken = tmp_path / "broken.json"
    broken.write_text("{", encoding="utf-8")
    with pytest.raises(plan.PlanRefusal, match="not JSON"):
        plan.run(_argv(source, toolchain, cache, broken, task))
    broken.write_text("[]", encoding="utf-8")
    with pytest.raises(plan.PlanRefusal, match="not an object"):
        plan.run(_argv(source, toolchain, cache, broken, task))
    broken.write_text(
        json.dumps(
            {
                "source_commit": "0" * 40,
                "patch_sha256": plan.PATCH_SHA256,
                "cargo_lock_before": plan.CARGO_LOCK_SHA256,
                "environment": {},
            }
        ),
        encoding="utf-8",
    )
    with pytest.raises(
        plan.PlanRefusal, match="source commit is not the accepted commit"
    ):
        plan.run(_argv(source, toolchain, cache, broken, task))

    result = plan.run(_argv(source, toolchain, cache, manifest, task))
    with pytest.raises(plan.PlanRefusal, match="not empty"):
        plan.run(_argv(source, toolchain, cache, manifest, task))
    with pytest.raises(plan.PlanRefusal, match="overwrite"):
        plan._write_plan(task / "reproduce-plan.json", result)
    assert not (task / "cargo-home").exists()


def test_whitespace_attribute_is_limited_to_frozen_evidence():
    attributes = (_REPO / ".gitattributes").read_text(encoding="utf-8")
    assert "* -whitespace" not in attributes
    assert "scripts/peri_a0_reproduce" not in attributes
    listed = []
    for line in attributes.splitlines():
        if not line or line.startswith("#"):
            continue
        path, rule = line.split(" ", 1)
        assert rule == "whitespace=-blank-at-eol,-blank-at-eof"
        assert path.startswith(
            ("scripts/peri_a0_closeout/evidence/", "scripts/peri_a0_revision/evidence/")
        )
        assert (_REPO / path).is_file()
        listed.append(path)
    assert "scripts/peri_a0_closeout/evidence/candidate-testfix.patch" in listed

    output = subprocess.check_output(
        [
            "git",
            "check-attr",
            "whitespace",
            "--",
            "scripts/peri_a0_reproduce/plan.py",
            "docs/development/peri-a0-archive.md",
            "scripts/peri_a0_closeout/evidence/candidate-testfix.patch",
        ],
        cwd=_REPO,
        text=True,
    )
    assert "scripts/peri_a0_reproduce/plan.py: whitespace: unspecified" in output
    assert "docs/development/peri-a0-archive.md: whitespace: unspecified" in output
    assert (
        "scripts/peri_a0_closeout/evidence/candidate-testfix.patch: "
        "whitespace: -blank-at-eol,-blank-at-eof"
    ) in output


def test_ruff_exclude_lists_frozen_python_files_only():
    text = (_REPO / "pyproject.toml").read_text(encoding="utf-8")
    section = text.split("[tool.ruff]", 1)[1].split("[tool.ruff.lint]", 1)[0]
    array = section.split("extend-exclude = [", 1)[1].split("]", 1)[0]
    paths = re.findall(r'"([^"]+)"', array)
    assert paths
    assert "scripts/peri_a0" not in paths
    assert "scripts/peri_a0_revision" not in paths
    assert "scripts/peri_a0_closeout" not in paths
    assert not any(item.startswith("scripts/peri_a0_reproduce") for item in paths)
    for relative in paths:
        assert relative.endswith(".py")
        assert relative.startswith(
            (
                "scripts/peri_a0/",
                "scripts/peri_a0_revision/",
                "scripts/peri_a0_closeout/",
            )
        )
        assert (_REPO / relative).is_file()
