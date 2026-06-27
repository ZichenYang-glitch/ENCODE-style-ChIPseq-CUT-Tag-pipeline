"""CLI entry point: encode-dag.

Provides snapshot/diff subcommands for smoke-test DAG dry-runs.
"""

import argparse
import csv
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from encode_pipeline.cli._logging import emit, make_context


REQUIRED_REPO_FILES = ("workflow" / Path("Snakefile"), "test" / Path("profiles"))


def _resolve_repo_root(repo_root=None):
    """Resolve the pipeline repository root.

    Resolution order:
    1. Explicit ``repo_root`` argument.
    2. ``ENCODE_PIPELINE_REPO_ROOT`` environment variable.
    3. Walk upward from the current working directory looking for
       ``workflow/Snakefile`` and ``test/profiles``.

    Raises:
        RuntimeError: if the repository root cannot be found.
    """
    if repo_root:
        return Path(repo_root).resolve()

    env_root = os.environ.get("ENCODE_PIPELINE_REPO_ROOT")
    if env_root:
        return Path(env_root).resolve()

    current = Path.cwd().resolve()
    for parent in [current, *current.parents]:
        if all((parent / rel).exists() for rel in REQUIRED_REPO_FILES):
            return parent

    raise RuntimeError(
        "Could not locate encode-pipeline repository root. "
        "Use --repo-root or set ENCODE_PIPELINE_REPO_ROOT."
    )


def _discover_placeholders(samples_tsv_path):
    """Parse samples.tsv and return placeholder file paths needed for dry-run."""
    paths = set()
    with open(samples_tsv_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            for key in ("fastq_1", "fastq_2", "control_bam"):
                val = (row.get(key) or "").strip()
                if val:
                    paths.add(val)
    return paths


def _rewrite_profile_config(profile_config_path, workdir):
    """Rewrite the samples path in a profile config to point into workdir."""
    with open(profile_config_path) as fh:
        content = fh.read()
    abs_samples = os.path.join(workdir, "samples.tsv")
    return re.sub(
        r"^samples:.*$",
        f'samples: "{abs_samples}"',
        content,
        flags=re.MULTILINE,
    )


def _prepare_profile_workdir(profile, profiles_dir):
    """Create a temporary workdir with placeholders and rewritten config."""
    profile_dir = profiles_dir / profile
    workdir = tempfile.mkdtemp(prefix=f"dag_{profile}_", dir="/tmp")
    samples_tsv_src = profile_dir / "samples.tsv"
    config_yaml_src = profile_dir / "config.yaml"

    for rel_path in _discover_placeholders(samples_tsv_src):
        placeholder_path = Path(workdir) / os.path.basename(rel_path)
        with open(placeholder_path, "w"):
            pass

    dest_samples = Path(workdir) / "samples.tsv"
    shutil.copy2(samples_tsv_src, dest_samples)

    rewritten_config = _rewrite_profile_config(config_yaml_src, workdir)
    dest_config = Path(workdir) / "config.yaml"
    with open(dest_config, "w") as fh:
        fh.write(rewritten_config)

    return workdir, str(dest_config)


def _extract_rule_names(stdout, stderr):
    """Parse scheduled rule names from snakemake -n output."""
    output = stdout + stderr
    rules = set()
    for line in output.splitlines():
        match = re.match(r"^\s*(?:local)?rule\s+(\w+)\s*:", line)
        if match:
            rules.add(match.group(1))
    return sorted(rules)


def _find_snakemake():
    """Resolve snakemake executable from env, PATH, or common conda paths."""
    env_value = os.environ.get("SNAKEMAKE", "")
    if env_value:
        return env_value
    path_value = shutil.which("snakemake")
    if path_value:
        return path_value
    candidates = [
        Path.home() / "miniconda3" / "envs" / "chipseq" / "bin" / "snakemake",
        Path.home() / "miniconda3" / "bin" / "snakemake",
        Path.home() / "micromamba" / "envs" / "chipseq" / "bin" / "snakemake",
    ]
    for candidate in candidates:
        if candidate.exists() and os.access(candidate, os.X_OK):
            return str(candidate)
    return None


def _run_snakemake_dryrun(config_path, snakefile):
    """Run snakemake -n and return sorted scheduled rule names."""
    snakemake = _find_snakemake()
    if not snakemake:
        raise RuntimeError("snakemake executable not found on PATH")

    workdir = os.path.dirname(config_path)
    env = os.environ.copy()
    if not env.get("XDG_CACHE_HOME"):
        env["XDG_CACHE_HOME"] = "/tmp/encode-pipeline-snakemake-cache"
    result = subprocess.run(
        [
            snakemake,
            "-s",
            str(snakefile),
            "--configfile",
            config_path,
            "-n",
        ],
        cwd=workdir,
        capture_output=True,
        text=True,
        env={**env, "PYTHONDONTWRITEBYTECODE": "1"},
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"snakemake -n failed:\n{result.stdout}\n{result.stderr}"
        )
    return _extract_rule_names(result.stdout, result.stderr)


def _snapshot(profile, repo_root):
    snakefile = repo_root / "workflow" / "Snakefile"
    profiles_dir = repo_root / "test" / "profiles"
    workdir, config_path = _prepare_profile_workdir(profile, profiles_dir)
    try:
        return _run_snakemake_dryrun(config_path, snakefile)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


def _read_snapshot(snapshot_path):
    with open(snapshot_path) as fh:
        return [
            line.strip()
            for line in fh
            if line.strip() and not line.strip().startswith("#")
        ]


def _cmd_snapshot(args):
    rules = _snapshot(args.profile, args.repo_root)
    lines = ["# DAG snapshot", "# Generated by: encode-dag snapshot", ""] + rules
    output = "\n".join(lines) + "\n"
    if args.output:
        with open(args.output, "w") as fh:
            fh.write(output)
        print(f"Snapshot written: {args.output}")
    else:
        print(output, end="")


def _cmd_diff(args):
    rules = _snapshot(args.profile, args.repo_root)
    snapshot_path = args.repo_root / "test" / "fixtures" / "dag_snapshots" / f"{args.profile}.txt"
    if not snapshot_path.exists():
        raise FileNotFoundError(f"Snapshot missing: {snapshot_path}")
    expected = _read_snapshot(snapshot_path)

    if rules == expected:
        print(f"No differences for {args.profile}")
        return 0

    print(f"DAG diff for {args.profile}:")
    for rule in sorted(set(expected) - set(rules)):
        print(f"- {rule}")
    for rule in sorted(set(rules) - set(expected)):
        print(f"+ {rule}")
    return 1


def main(argv=None):
    parser = argparse.ArgumentParser(description="DAG snapshot and diff utility")
    sub = parser.add_subparsers(dest="command", required=True)

    snap = sub.add_parser("snapshot", help="Generate a DAG rule snapshot")
    snap.add_argument("--profile", required=True, help="Smoke-test profile name")
    snap.add_argument("--output", default=None, help="Output file (default: stdout)")
    snap.add_argument("--repo-root", default=None, help="Path to pipeline repository root")

    diff = sub.add_parser("diff", help="Diff current DAG against stored snapshot")
    diff.add_argument("--profile", required=True, help="Smoke-test profile name")
    diff.add_argument("--repo-root", default=None, help="Path to pipeline repository root")

    args = parser.parse_args(argv)

    try:
        args.repo_root = _resolve_repo_root(args.repo_root)
    except RuntimeError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2

    run_id, version, config_hash = make_context()
    emit(run_id, version, config_hash)

    if args.command == "snapshot":
        _cmd_snapshot(args)
        return 0
    if args.command == "diff":
        return _cmd_diff(args)
    return 2


if __name__ == "__main__":
    sys.exit(main())
