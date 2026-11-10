"""Build an explicitly patched, pinned candidate using private copies of old caches."""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time

ROOT = Path("/tmp/helixweave-peri-a0-revision-20261109")
SOURCE = ROOT / "peri-candidate"
OLD = Path("/tmp/helixweave-peri-a0-20261109")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--task-root", type=Path, default=ROOT)
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--reuse-root", type=Path, default=OLD)
    parser.add_argument("--label", default="final")
    args = parser.parse_args()
    if args.label not in {"final"}:
        raise ValueError("unknown build label")
    root, source, old = (
        args.task_root.resolve(),
        args.source.resolve(),
        args.reuse_root.resolve(),
    )
    manifest = root / "evidence/build-final.json"
    assert not manifest.exists(), "do not overwrite build evidence"
    original = json.loads((old / "evidence/build.json").read_text())
    assert (
        subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=source, text=True
        ).strip()
        == original["source_commit"]
    )
    patch = subprocess.check_output(["git", "diff", "HEAD"], cwd=source)
    assert patch
    (root / "evidence/candidate-final.patch").write_bytes(patch)
    env = dict(original["environment"])
    for key in ["HOME", "CARGO_HOME", "CARGO_TARGET_DIR"]:
        env[key] = env[key].replace(str(old), str(root))
    env["CARGO_NET_OFFLINE"] = "true"
    env["CARGO_BUILD_JOBS"] = "4"
    command = original["command"]
    record = {
        **original,
        "binary": str(root / "build/target/release/peri"),
        "command": command,
        "environment": env,
        "local_patch": patch.decode(),
        "patch_sha256": hashlib.sha256(patch).hexdigest(),
        "success": False,
        "started_at_unix": time.time(),
        "cargo_lock_before": digest(source / "Cargo.lock"),
        "cache_copy_from": str(old),
        "toolchain_reused_read_only": str(old / "toolchain"),
    }
    manifest = root / "evidence/build-final.json"
    assert not manifest.exists(), "do not overwrite build evidence"
    try:
        with (root / "evidence/native-build-final.log").open("xb") as log:
            result = subprocess.run(
                command,
                cwd=source,
                env=env,
                stdout=log,
                stderr=subprocess.STDOUT,
                timeout=900,
            )
        record["exit_code"] = result.returncode
        if result.returncode == 0:
            record["binary_sha256"] = digest(Path(record["binary"]))
            record["success"] = True
    finally:
        record["finished_at_unix"] = time.time()
        record["cargo_lock_after"] = digest(source / "Cargo.lock")
        assert subprocess.check_output(["git", "diff", "HEAD"], cwd=source) == patch
        manifest.write_text(json.dumps(record, indent=2) + "\n")
    print(
        json.dumps(
            {
                k: record.get(k)
                for k in ["success", "exit_code", "binary_sha256", "patch_sha256"]
            }
        ),
        flush=True,
    )
    return 0 if record["success"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
