"""Plan a relocated Peri A0 build without starting Cargo.

The frozen ``build_candidate.py`` copies a historical manifest environment and
rewrites only ``HOME``, ``CARGO_HOME``, and ``CARGO_TARGET_DIR``. This planner
does not call it. Path settings come from the current toolchain, cache source,
and task directory. ``--online-locked`` belongs to ``verify_patch.py`` and is
not accepted here.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

PERI_COMMIT = "c7d154c9985cbe08ba25e964cb7468d15014e6d0"
PATCH_SHA256 = "06272e1a66423a9e562509a38d044cd8506a1c10e61ef1fd3fe239107680587b"
CARGO_LOCK_SHA256 = "0692683359fb75073eca3d14b364071bc02058be7b11e0335ccaf2b0f314f5dd"
RUST_VERSION = "1.99.0"
BUILD_COMMAND = (
    "./scripts/cargo-rmcp-patched.sh",
    "build",
    "--locked",
    "-p",
    "peri-tui",
    "--release",
    "--bin",
    "peri",
)
ORIGINAL_EXPERIMENT_ROOTS = (
    Path("/tmp/hw-peri-a0-close-20261109"),
    Path("/tmp/helixweave-peri-a0-20261109"),
    Path("/tmp/helixweave-peri-design-20261108"),
)
_SYSTEM_ROOTS = (
    "/bin",
    "/dev",
    "/etc",
    "/lib",
    "/lib64",
    "/proc",
    "/usr/bin",
    "/usr/lib",
    "/usr/lib64",
)
_CONSTRUCTED_KEYS = {
    "CARGO_BUILD_JOBS",
    "CARGO_HOME",
    "CARGO_NET_OFFLINE",
    "CARGO_TARGET_DIR",
    "GIT_CONFIG_GLOBAL",
    "GIT_CONFIG_NOSYSTEM",
    "HOME",
    "LANG",
    "PATH",
}
_COPY_IF_SAFE = {
    "CARGO_NET_GIT_FETCH_WITH_CLI",
    "GIT_CONFIG_COUNT",
    "GIT_CONFIG_KEY_0",
    "GIT_CONFIG_VALUE_0",
    "GIT_SSH_COMMAND",
}
_HISTORICAL_PATH_FIELDS = (
    "binary",
    "cache_copy_from",
    "toolchain_reused_read_only",
)


class PlanRefusal(Exception):
    """The preflight stopped before writing a plan or starting a build."""


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def is_system_path(token: str) -> bool:
    return any(token == root or token.startswith(root + "/") for root in _SYSTEM_ROOTS)


def absolute_tokens(value: str) -> list[str]:
    tokens = []
    for chunk in value.split():
        for part in chunk.split(":"):
            if part.startswith("/") and not part.startswith("//"):
                tokens.append(part)
    return tokens


def historical_local_paths(manifest: dict) -> list[str]:
    values = []
    environment = manifest.get("environment")
    if isinstance(environment, dict):
        values.extend(str(item) for item in environment.values())
    command = manifest.get("command")
    if isinstance(command, list):
        values.extend(str(item) for item in command)
    for key in _HISTORICAL_PATH_FIELDS:
        if key in manifest:
            values.append(str(manifest[key]))
    found = []
    for value in values:
        for token in absolute_tokens(value):
            if not is_system_path(token):
                found.append(token)
    return list(dict.fromkeys(found))


def _resolve(path: Path) -> Path:
    return path.expanduser().resolve()


def _is_relative_to(path: Path, root: Path) -> bool:
    try:
        path.resolve().relative_to(root.resolve())
    except ValueError:
        return False
    return True


def _conflicts(path: Path, root: Path) -> bool:
    return _is_relative_to(path, root) or _is_relative_to(root, path)


def _require_file(path: Path, label: str) -> None:
    if not path.is_file():
        raise PlanRefusal(f"{label} is not a file")


def _require_dir(path: Path, label: str) -> None:
    if not path.is_dir():
        raise PlanRefusal(f"{label} is not a directory")


def _shared_environment(path: Path) -> bool:
    current = path.resolve()
    for parent in (current, *current.parents):
        if (parent / "conda-meta").is_dir() or (parent / "pyvenv.cfg").is_file():
            return True
    return False


def _git_output(source: Path, *args: str, text: bool) -> str | bytes:
    try:
        return subprocess.check_output(
            ["git", "-C", str(source), *args],
            stderr=subprocess.PIPE,
            text=text,
        )
    except (OSError, subprocess.CalledProcessError) as exc:
        detail = exc.stderr if isinstance(exc, subprocess.CalledProcessError) else ""
        if isinstance(detail, bytes):
            detail = detail.decode("utf-8", "replace")
        message = f"git {' '.join(args)} failed while checking the Peri source"
        if str(detail).strip():
            message = message + ": " + str(detail).strip()
        raise PlanRefusal(message) from exc


def source_commit(source: Path) -> str:
    output = _git_output(source, "rev-parse", "HEAD", text=True)
    if not isinstance(output, str):
        raise PlanRefusal("git rev-parse returned an unexpected result")
    return output.strip()


def source_diff(source: Path) -> bytes:
    output = _git_output(source, "diff", "HEAD", text=False)
    if not isinstance(output, bytes):
        raise PlanRefusal("git diff returned an unexpected result")
    return output


def lock_sha256(source: Path) -> str:
    return file_sha256(source / "Cargo.lock")


def toolchain_channel(source: Path) -> str:
    channel = None
    text = (source / "rust-toolchain.toml").read_text(encoding="utf-8")
    for raw in text.splitlines():
        line = raw.split("#", 1)[0].strip()
        if line.startswith("channel") and "=" in line:
            channel = line.split("=", 1)[1].strip().strip('"').strip("'")
    if channel is None:
        raise PlanRefusal("rust-toolchain.toml has no channel")
    return channel


def _command_version(binary: Path) -> None:
    if not os.access(binary, os.X_OK):
        raise PlanRefusal(f"toolchain binary is not executable: {binary.name}")
    try:
        output = subprocess.check_output(
            [str(binary), "-V"],
            env={"PATH": f"{binary.parent}:/usr/bin:/bin", "LANG": "C.UTF-8"},
            text=True,
            stderr=subprocess.PIPE,
            timeout=30,
        )
    except subprocess.TimeoutExpired as exc:
        raise PlanRefusal(f"{binary.name} -V timed out") from exc
    except (OSError, subprocess.CalledProcessError) as exc:
        raise PlanRefusal(f"{binary.name} -V failed") from exc
    first = output.splitlines()[0] if output else ""
    expected = f"{binary.name} {RUST_VERSION}"
    if not first.startswith(expected):
        raise PlanRefusal(f"{binary.name} is not {RUST_VERSION}")


def _load_manifest(path: Path) -> dict:
    try:
        manifest = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise PlanRefusal("historical manifest is not JSON") from exc
    if not isinstance(manifest, dict):
        raise PlanRefusal("historical manifest is not an object")
    if manifest.get("source_commit") != PERI_COMMIT:
        raise PlanRefusal(
            "historical manifest source commit is not the accepted commit"
        )
    if manifest.get("patch_sha256") != PATCH_SHA256:
        raise PlanRefusal("historical manifest patch digest is not the accepted patch")
    if manifest.get("cargo_lock_before") != CARGO_LOCK_SHA256:
        raise PlanRefusal("historical manifest lock digest is not the accepted lock")
    if not isinstance(manifest.get("environment"), dict):
        raise PlanRefusal("historical manifest has no environment object")
    return manifest


def _refuse_layout(
    *,
    source: Path,
    cache: Path,
    toolchain: Path,
    task_root: Path,
    historical_paths: list[str],
) -> None:
    for root in ORIGINAL_EXPERIMENT_ROOTS:
        if _conflicts(task_root, root):
            raise PlanRefusal("task output is inside an original experiment directory")
        if _conflicts(source, root):
            raise PlanRefusal("source is inside an original experiment directory")
    if _shared_environment(task_root):
        raise PlanRefusal("task output is inside a shared environment")
    if _conflicts(task_root, cache):
        raise PlanRefusal("cache source must stay outside the task output")
    if _conflicts(task_root, toolchain):
        raise PlanRefusal("toolchain must stay outside the task output")
    for historical in historical_paths:
        if _conflicts(task_root, Path(historical)):
            raise PlanRefusal("task output overlaps a historical local path")
        if _conflicts(source, Path(historical)):
            raise PlanRefusal("source overlaps a historical local path")


def _path_allowed(token: str, toolchain: Path, task_root: Path) -> bool:
    if is_system_path(token):
        return True
    path = Path(token)
    return _is_relative_to(path, toolchain) or _is_relative_to(path, task_root)


def build_environment(
    manifest: dict,
    *,
    toolchain: Path,
    task_root: Path,
) -> tuple[dict[str, str], list[str], list[str]]:
    """Build the Cargo environment from the current roots, not historical paths."""
    environment = {
        "PATH": f"{toolchain / 'bin'}:/usr/bin:/bin",
        "HOME": str(task_root / "build-home"),
        "CARGO_HOME": str(task_root / "cargo-home"),
        "CARGO_TARGET_DIR": str(task_root / "build" / "target"),
        "CARGO_BUILD_JOBS": "4",
        "CARGO_NET_OFFLINE": "true",
        "LANG": "C.UTF-8",
        "GIT_CONFIG_GLOBAL": "/dev/null",
        "GIT_CONFIG_NOSYSTEM": "1",
    }
    replaced = []
    dropped = []
    historical = manifest["environment"]
    for key, value in historical.items():
        text = str(value)
        if key in _CONSTRUCTED_KEYS:
            replaced.append(key)
            continue
        unsafe = [
            token
            for token in absolute_tokens(text)
            if not _path_allowed(token, toolchain, task_root)
        ]
        if key not in _COPY_IF_SAFE or unsafe:
            dropped.append(str(key))
            continue
        environment[key] = text
    for key, value in environment.items():
        for token in absolute_tokens(value):
            if not _path_allowed(token, toolchain, task_root):
                raise PlanRefusal(
                    f"constructed environment keeps an unrelated path: {token}"
                )
    return environment, sorted(replaced), sorted(dropped)


def _patch_application(source: Path) -> str:
    diff = source_diff(source)
    if not diff:
        return "apply-before-build"
    if hashlib.sha256(diff).hexdigest() == PATCH_SHA256:
        return "already-applied"
    raise PlanRefusal("source diff is not the accepted patch")


def _verify_source(source: Path) -> str:
    if source_commit(source) != PERI_COMMIT:
        raise PlanRefusal("source commit is not the accepted Peri commit")
    if not (source / "scripts" / "cargo-rmcp-patched.sh").is_file():
        raise PlanRefusal("source is missing scripts/cargo-rmcp-patched.sh")
    if not (source / "Cargo.lock").is_file():
        raise PlanRefusal("source is missing Cargo.lock")
    if lock_sha256(source) != CARGO_LOCK_SHA256:
        raise PlanRefusal("Cargo.lock is not the accepted A0 lock")
    if toolchain_channel(source) != RUST_VERSION:
        raise PlanRefusal("rust-toolchain.toml channel is not 1.99.0")
    return _patch_application(source)


def _verify_cache(cache: Path) -> None:
    if (cache / "registry").is_dir() or (cache / "git").is_dir():
        return
    raise PlanRefusal("cache source has neither registry nor git")


def _verify_toolchain(toolchain: Path) -> None:
    bindir = toolchain / "bin"
    _require_dir(bindir, "toolchain bin")
    _command_version(bindir / "rustc")
    _command_version(bindir / "cargo")


def _prepare_output(task_root: Path, source: Path) -> None:
    allowed: set[str] = set()
    if task_root.exists():
        if not task_root.is_dir():
            raise PlanRefusal("task output is not a directory")
        if _is_relative_to(source, task_root):
            if source.parent != task_root:
                raise PlanRefusal(
                    "source inside the task output must be a direct child"
                )
            allowed.add(source.name)
        if any(child.name not in allowed for child in task_root.iterdir()):
            raise PlanRefusal("task output is not empty")
        return
    if not task_root.parent.is_dir():
        raise PlanRefusal("task output parent does not exist")
    task_root.mkdir()


def _guard_command(command: list[str], toolchain: Path, task_root: Path) -> None:
    for item in command:
        for token in absolute_tokens(item):
            if not _path_allowed(token, toolchain, task_root):
                raise PlanRefusal(
                    f"constructed command keeps an unrelated path: {token}"
                )


def assemble_plan(
    manifest: dict,
    *,
    source: Path,
    patch: Path,
    cache: Path,
    toolchain: Path,
    task_root: Path,
    historical_manifest: Path,
    patch_application: str,
) -> dict:
    environment, replaced, dropped = build_environment(
        manifest, toolchain=toolchain, task_root=task_root
    )
    command = list(BUILD_COMMAND)
    _guard_command(command, toolchain, task_root)
    return {
        "cache_source": str(cache),
        "cargo_lock_sha256": CARGO_LOCK_SHA256,
        "command": command,
        "cwd": str(source),
        "dropped_historical_keys": dropped,
        "environment": environment,
        "historical_manifest": str(historical_manifest),
        "inherited_non_path_keys": sorted(set(environment) - _CONSTRUCTED_KEYS),
        "online_locked_entrypoint": "scripts/peri_a0_closeout/verify_patch.py",
        "online_locked_supported": False,
        "patch": str(patch),
        "patch_application": patch_application,
        "patch_sha256": PATCH_SHA256,
        "replaced_historical_keys": replaced,
        "source_commit": PERI_COMMIT,
        "starts_build": False,
        "task_root": str(task_root),
        "toolchain": str(toolchain),
        "toolchain_version": RUST_VERSION,
        "writes": [
            str(task_root / "reproduce-plan.json"),
            str(task_root / "build-home"),
            str(task_root / "cargo-home"),
            str(task_root / "build" / "target"),
        ],
    }


def _write_plan(path: Path, plan: dict) -> None:
    payload = json.dumps(plan, indent=2, sort_keys=True) + "\n"
    try:
        descriptor = os.open(path, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o644)
    except FileExistsError as exc:
        raise PlanRefusal(
            "refusing to overwrite an existing reproduction plan"
        ) from exc
    with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
        handle.write(payload)


def run(argv: list[str] | None = None) -> dict:
    parser = argparse.ArgumentParser(
        description="Plan a relocated Peri A0 build without starting Cargo."
    )
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--patch", type=Path, required=True)
    parser.add_argument("--toolchain", type=Path, required=True)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--task-root", type=Path, required=True)
    parser.add_argument("--historical-manifest", type=Path, required=True)
    args = parser.parse_args(argv)
    source = _resolve(args.source)
    patch = _resolve(args.patch)
    cache = _resolve(args.cache)
    toolchain = _resolve(args.toolchain)
    task_root = _resolve(args.task_root)
    historical_manifest = _resolve(args.historical_manifest)
    _require_dir(source, "source")
    _require_file(patch, "patch")
    _require_dir(cache, "cache source")
    _require_dir(toolchain, "toolchain")
    _require_file(historical_manifest, "historical manifest")
    manifest = _load_manifest(historical_manifest)
    _refuse_layout(
        source=source,
        cache=cache,
        toolchain=toolchain,
        task_root=task_root,
        historical_paths=historical_local_paths(manifest),
    )
    if file_sha256(patch) != PATCH_SHA256:
        raise PlanRefusal("patch file is not the accepted A0 patch")
    _verify_cache(cache)
    _verify_toolchain(toolchain)
    patch_application = _verify_source(source)
    plan = assemble_plan(
        manifest,
        source=source,
        patch=patch,
        cache=cache,
        toolchain=toolchain,
        task_root=task_root,
        historical_manifest=historical_manifest,
        patch_application=patch_application,
    )
    _prepare_output(task_root, source)
    _write_plan(task_root / "reproduce-plan.json", plan)
    return plan


def main(argv: list[str] | None = None) -> int:
    try:
        plan = run(argv)
    except PlanRefusal as error:
        print(error, file=sys.stderr)
        return 2
    print(json.dumps(plan, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
