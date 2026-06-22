"""Stage 59 env pinning enforcement test.

Scans workflow/envs/*.yml and requires every top-level core CLI tool
dependency to have a version constraint (=, ==, >=, <=, <, ~=).
Does not require transitive or non-core dependencies to be pinned.
"""

import os
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ENVS_DIR = os.path.join(REPO_ROOT, "workflow", "envs")

# Core top-level CLI tools that must be version-constrained
CORE_TOOLS = {
    "fastqc",
    "trim-galore",
    "cutadapt",
    "bowtie2",
    "samtools",
    "bedtools",
    "macs3",
    "deeptools",
    "idr",
    "multiqc",
    "picard",
    "phantompeakqualtools",
    "preseq",
    "seacr",
    "ucsc-bedgraphtobigwig",
}

# Tools that appear in env files but are NOT required to be pinned
# (e.g., pure Python packages, conda infrastructure)
NON_CORE = {
    "python",
    "snakemake",
    "snakemake-minimal",
    "openjdk",
    "pyyaml",
}

VERSION_CONSTRAINT_CHARS = set("=<>~")


def parse_env_file(path):
    """Parse a conda env YAML and return list of (line_number, raw_line, package_name).

    Handles the simple subset of YAML used in these env files:
    dependencies as a list of strings under 'dependencies:'.
    """
    results = []
    with open(path) as f:
        lines = f.readlines()

    in_deps = False
    indent = 0
    for i, _raw_line in enumerate(lines, start=1):
        line = _raw_line.rstrip("\n")
        stripped = line.strip()

        if stripped.startswith("dependencies:"):
            in_deps = True
            indent = len(line) - len(line.lstrip())
            continue

        if not in_deps:
            continue

        # Detect end of dependencies block
        if stripped and not stripped.startswith("- ") and not stripped.startswith("#"):
            current_indent = len(line) - len(line.lstrip())
            if current_indent <= indent and not stripped.startswith("-"):
                in_deps = False
                continue

        # Parse dependency line
        if stripped.startswith("- "):
            dep_text = stripped[2:].strip()
            # Handle inline comments
            if "#" in dep_text and not dep_text.startswith("#"):
                dep_text = dep_text.split("#")[0].strip()
            if dep_text:
                # Extract package name (first word or everything before version spec)
                pkg_name = _extract_package_name(dep_text)
                results.append((i, _raw_line.rstrip("\n"), pkg_name, dep_text))

    return results


def _extract_package_name(dep_text):
    """Extract the package name from a dependency line.

    Examples:
      "samtools >=1.16,<2" -> "samtools"
      "samtools" -> "samtools"
      "python >=3.10" -> "python"
    """
    text = dep_text.strip()
    for sep in (" ", "\t"):
        if sep in text:
            return text.split(sep)[0].strip()
    return text


def has_version_constraint(dep_text):
    """Check if a dependency text has a version constraint after the package name."""
    text = dep_text.strip()
    pkg = _extract_package_name(text)
    rest = text[len(pkg):].strip()
    if not rest:
        return False
    return any(c in VERSION_CONSTRAINT_CHARS for c in rest)


def main():
    errors = []
    ok_count = 0
    total_core = 0

    for fname in sorted(os.listdir(ENVS_DIR)):
        if not fname.endswith(".yml"):
            continue
        fpath = os.path.join(ENVS_DIR, fname)
        relpath = os.path.relpath(fpath, REPO_ROOT)

        for line_no, raw_line, pkg_name, dep_text in parse_env_file(fpath):
            if pkg_name not in CORE_TOOLS:
                continue
            total_core += 1
            if has_version_constraint(dep_text):
                ok_count += 1
            else:
                errors.append(
                    f"  {relpath}: line {line_no}: {pkg_name!r} has no "
                    f"version constraint"
                )

    if errors:
        print(f"FAIL: {len(errors)} unpinned core tool(s):\n")
        for err in errors:
            print(err)
        print(f"\n  {ok_count}/{total_core} core tool entries are pinned.")
        print(f"  Add a version constraint (e.g., '>=major.minor,<next-major')")
        print(f"  to each unpinned entry.")
        return 1

    print(f"PASS: All {total_core} core tool entries have version constraints")
    return 0


if __name__ == "__main__":
    sys.exit(main())
