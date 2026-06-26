#!/usr/bin/env python3
"""Compare current `snakemake --lint` output to a known-good baseline.

Run with --update-baseline to regenerate
`docs/operations/snakemake-lint-warnings.txt` after deliberate changes.
"""

import argparse
import difflib
import os
import subprocess
import sys
from pathlib import Path

# Allow importing sibling test helpers.
REPO_ROOT = Path(__file__).resolve().parent.parent
BASELINE = REPO_ROOT / "docs" / "operations" / "snakemake-lint-warnings.txt"


def _snakemake_command():
    """Resolve the snakemake executable using the same helper tests use."""
    sys.path.insert(0, str(REPO_ROOT / "test"))
    from _tool_resolver import resolve_tool

    return resolve_tool("snakemake", "SNAKEMAKE")


def _normalize(text, repo_root):
    """Remove machine-specific absolute paths and trailing whitespace."""
    root = str(repo_root)
    # Match both '/root/' prefix and '/root' at end of string/line.
    text = text.replace(root + os.sep, "")
    text = text.replace(root, ".")
    # Strip trailing whitespace from each line and remove trailing blank lines
    # so the baseline stays stable regardless of how snakemake terminates its
    # output.
    lines = [line.rstrip() for line in text.splitlines()]
    while lines and lines[-1] == "":
        lines.pop()
    return "\n".join(lines) + "\n"


def _run_lint(repo_root):
    """Run `snakemake --lint` and return normalized stdout+stderr."""
    cmd = [_snakemake_command(), "--lint", "--workflow-profile", "none"]
    env = os.environ.copy()
    if not env.get("XDG_CACHE_HOME"):
        env["XDG_CACHE_HOME"] = "/tmp/encode-pipeline-snakemake-cache"
    result = subprocess.run(
        cmd,
        cwd=repo_root,
        capture_output=True,
        text=True,
        env=env,
    )
    combined = result.stdout + result.stderr
    return _normalize(combined, repo_root)


def main():
    parser = argparse.ArgumentParser(description="Check snakemake --lint baseline")
    parser.add_argument(
        "--update-baseline",
        action="store_true",
        help="Write the current lint output as the new baseline",
    )
    args = parser.parse_args()

    output = _run_lint(REPO_ROOT)

    if args.update_baseline:
        BASELINE.parent.mkdir(parents=True, exist_ok=True)
        BASELINE.write_text(output)
        print(f"Updated baseline: {BASELINE}")
        return 0

    if not BASELINE.exists():
        print(f"Baseline missing: {BASELINE}")
        print("Run with --update-baseline to create it.")
        return 1

    baseline_text = BASELINE.read_text()
    if output == baseline_text:
        print("snakemake --lint output matches baseline.")
        return 0

    diff = difflib.unified_diff(
        baseline_text.splitlines(keepends=True),
        output.splitlines(keepends=True),
        fromfile=str(BASELINE),
        tofile="current snakemake --lint output",
    )
    sys.stdout.writelines(diff)
    print(
        "\nERROR: snakemake --lint produced new or changed warnings.\n"
        "If these warnings are expected, run:\n\n"
        "  python3 test/check_snakemake_lint.py --update-baseline\n"
    )
    return 1


if __name__ == "__main__":
    sys.exit(main())
