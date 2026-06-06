#!/usr/bin/env python3
"""Guard against workstation-specific paths and stale env references."""

import re
import sys
from pathlib import Path

import yaml


REPO_ROOT = Path(__file__).resolve().parents[1]

RUNTIME_ROOTS = [
    "scripts",
    "workflow",
    "config",
    "test",
]

ACTIVE_DOCS = [
    "README.md",
    "KNOWN_ISSUES.md",
    "docs/configuration.md",
    "docs/environments.md",
    "docs/quickstart.md",
    "docs/sample-sheet.md",
]

RUNTIME_FORBIDDEN = [
    ("/home/", "home-directory absolute path"),
    ("/Users/", "macOS home-directory absolute path"),
    ("/data/", "site-specific data path"),
    ("/mnt/c/", "WSL user path"),
    ("C:\\", "Windows absolute path"),
]

DOC_FORBIDDEN = [
    ("/home/irenadler", "author workstation path"),
    ("miniconda3/envs/chipseq", "author Conda env path"),
]


def _iter_runtime_files():
    suffixes = {".py", ".sh", ".smk", ".yaml", ".yml", ".tsv"}
    skip = {
        Path("test/test_no_hardcoded_paths.py"),
        Path("test/test_stage43_artifact_inventory.py"),
    }
    for root in RUNTIME_ROOTS:
        base = REPO_ROOT / root
        if not base.exists():
            continue
        for path in base.rglob("*"):
            rel = path.relative_to(REPO_ROOT)
            if rel in skip or "__pycache__" in rel.parts:
                continue
            if path.is_file() and path.suffix in suffixes:
                yield path


def _check_forbidden(files, patterns):
    failures = []
    for path in files:
        rel = path.relative_to(REPO_ROOT)
        text = path.read_text(errors="replace")
        for lineno, line in enumerate(text.splitlines(), start=1):
            for pattern, reason in patterns:
                if pattern in line:
                    failures.append(f"{rel}:{lineno}: {reason}: {pattern}")
    return failures


def _check_rule_env_references():
    failures = []
    env_ref = re.compile(r'"(\.\./envs/[^"]+\.yml)"')
    for rule_file in (REPO_ROOT / "workflow" / "rules").glob("*.smk"):
        text = rule_file.read_text()
        for match in env_ref.finditer(text):
            ref = match.group(1)
            env_path = (rule_file.parent / ref).resolve()
            if not env_path.exists():
                failures.append(
                    f"{rule_file.relative_to(REPO_ROOT)} references missing env {ref}"
                )
            if env_path.name == "chipseq.yml":
                failures.append(
                    f"{rule_file.relative_to(REPO_ROOT)} references monolithic env {ref}"
                )
    return failures


def _check_env_files():
    failures = []
    for env_file in sorted((REPO_ROOT / "workflow" / "envs").glob("*.yml")):
        data = yaml.safe_load(env_file.read_text())
        channels = data.get("channels", [])
        deps = data.get("dependencies", [])
        rel = env_file.relative_to(REPO_ROOT)
        if "nodefaults" not in channels:
            failures.append(f"{rel}: missing nodefaults channel guard")
        if not deps:
            failures.append(f"{rel}: dependencies list is empty")
    return failures


def main():
    failures = []
    failures.extend(_check_forbidden(_iter_runtime_files(), RUNTIME_FORBIDDEN))
    failures.extend(
        _check_forbidden((REPO_ROOT / p for p in ACTIVE_DOCS), DOC_FORBIDDEN)
    )
    failures.extend(_check_rule_env_references())
    failures.extend(_check_env_files())

    if failures:
        print("No-hardcoding guard failed:\n")
        for failure in failures:
            print(f"FAIL: {failure}")
        return 1

    print("PASS: no hardcoded workstation paths or stale env references")
    return 0


if __name__ == "__main__":
    sys.exit(main())
