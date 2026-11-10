"""Reproduce a task-local build of the pinned source; never changes global config."""

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import tarfile
import time
import tomllib
import urllib.request

from configuration import PERI_COMMIT

RUST_VERSION = "1.99.0"
TARGET = "x86_64-unknown-linux-gnu"


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def download(url, path, expected=None):
    if not path.exists():
        temporary = path.with_suffix(path.suffix + ".partial")
        with (
            urllib.request.urlopen(url, timeout=120) as response,
            temporary.open("wb") as output,
        ):
            while data := response.read(1024 * 1024):
                output.write(data)
        temporary.replace(path)
    if expected is not None and digest(path) != expected:
        raise RuntimeError("archive checksum mismatch: " + path.name)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--task-root", type=Path, required=True)
    parser.add_argument(
        "--git-ssh",
        action="store_true",
        help="Use process-local Git CLI/SSH URL mapping for pinned Git dependencies",
    )
    args = parser.parse_args()
    source, root = args.source.resolve(), args.task_root.resolve()
    identity = subprocess.check_output(
        ["git", "-C", str(source), "rev-parse", "HEAD"], text=True
    ).strip()
    if identity != PERI_COMMIT or subprocess.check_output(
        ["git", "-C", str(source), "diff", "HEAD"]
    ):
        raise RuntimeError("source identity mismatch or unrecorded patch")
    for directory in ["downloads", "build", "evidence", "cargo-home", "build-home"]:
        (root / directory).mkdir(parents=True, exist_ok=True)
    pin = tomllib.loads((source / "rust-toolchain.toml").read_text())["toolchain"][
        "channel"
    ]
    if pin != RUST_VERSION:
        raise RuntimeError("toolchain pin changed")
    manifest = root / "downloads" / f"channel-rust-{pin}.toml"
    download(f"https://static.rust-lang.org/dist/channel-rust-{pin}.toml", manifest)
    packages = tomllib.loads(manifest.read_text())["pkg"]
    records = []
    for package in ["rustc", "cargo", "rust-std"]:
        item = packages[package]["target"][TARGET]
        archive = root / "downloads" / item["xz_url"].rsplit("/", 1)[1]
        download(item["xz_url"], archive, item["xz_hash"])
        records.append(
            {"name": package, "url": item["xz_url"], "sha256": digest(archive)}
        )
        marker = root / "toolchain" / (".a0-" + package)
        if not marker.exists():
            with tarfile.open(archive) as tar:
                tar.extractall(root / "build", filter="data")
            installer = (
                root / "build" / archive.name.removesuffix(".tar.xz") / "install.sh"
            )
            subprocess.run(
                [
                    "sh",
                    str(installer),
                    "--prefix=" + str(root / "toolchain"),
                    "--disable-ldconfig",
                ],
                check=True,
            )
            marker.write_text(item["xz_hash"] + "\n")
    submodule_pin = subprocess.check_output(
        ["git", "-C", str(source), "rev-parse", "HEAD:peri-cool"], text=True
    ).strip()
    submodule_actual = subprocess.check_output(
        ["git", "-C", str(source / "peri-cool"), "rev-parse", "HEAD"], text=True
    ).strip()
    if submodule_actual != submodule_pin:
        raise RuntimeError("peri-cool submodule identity mismatch")
    env = {
        "PATH": str(root / "toolchain" / "bin") + ":/usr/bin:/bin",
        "HOME": str(root / "build-home"),
        "CARGO_HOME": str(root / "cargo-home"),
        "CARGO_TARGET_DIR": str(root / "build" / "target"),
        "CARGO_BUILD_JOBS": "4",
        "LANG": "C.UTF-8",
        "GIT_CONFIG_GLOBAL": "/dev/null",
        "GIT_CONFIG_NOSYSTEM": "1",
    }
    if args.git_ssh:
        env.update(
            CARGO_NET_GIT_FETCH_WITH_CLI="true",
            GIT_CONFIG_COUNT="1",
            GIT_CONFIG_KEY_0="url.ssh://git@github.com/.insteadOf",
            GIT_CONFIG_VALUE_0="https://github.com/",
            GIT_SSH_COMMAND="ssh -o BatchMode=yes",
        )
    command = [
        "./scripts/cargo-rmcp-patched.sh",
        "build",
        "--locked",
        "-p",
        "peri-tui",
        "--release",
        "--bin",
        "peri",
    ]
    lock_digest = digest(source / "Cargo.lock")
    record = {
        "source_commit": identity,
        "local_patch": "",
        "rust_version": pin,
        "toolchain_packages": records,
        "command": command,
        "environment": env,
        "cargo_lock_before": lock_digest,
        "started_at_unix": time.time(),
    }
    record["rustc_vv"] = subprocess.check_output(["rustc", "-Vv"], env=env, text=True)
    record["cargo_version"] = subprocess.check_output(
        ["cargo", "-V"], env=env, text=True
    )
    record["submodule_gitlink"] = subprocess.check_output(
        ["git", "-C", str(source), "ls-tree", "HEAD", "peri-cool"], text=True
    ).strip()
    try:
        with (root / "evidence" / "native-build.log").open("wb") as log:
            subprocess.run(
                command,
                cwd=source,
                env=env,
                stdout=log,
                stderr=subprocess.STDOUT,
                check=True,
            )
        binary = root / "build" / "target" / "release" / "peri"
        record.update(binary=str(binary), binary_sha256=digest(binary), success=True)
    except BaseException:
        record["success"] = False
        raise
    finally:
        record["cargo_lock_after"] = digest(source / "Cargo.lock")
        record["finished_at_unix"] = time.time()
        record["local_patch"] = subprocess.check_output(
            ["git", "-C", str(source), "diff", "HEAD"], text=True
        )
        (root / "evidence" / "build.json").write_text(
            json.dumps(record, indent=2) + "\n"
        )
    if record["cargo_lock_after"] != lock_digest or record["local_patch"]:
        raise RuntimeError("build changed source or lock; do not run this binary")


if __name__ == "__main__":
    main()
