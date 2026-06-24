"""Shared pytest fixtures for the test suite."""

import csv
import glob
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import pytest


# Ignore legacy test_stage*.py files that have their own main()/__main__ guard;
# they are executed via test_stage_shim.py instead. Files without an entry point
# (e.g. test_stage8_smoke_profiles.py) remain normal pytest modules.
_TEST_DIR = os.path.dirname(os.path.abspath(__file__))
_LEGACY_DIR = os.path.join(_TEST_DIR, "legacy")
collect_ignore = []
for _search_dir in (_TEST_DIR, _LEGACY_DIR):
    if not os.path.isdir(_search_dir):
        continue
    for _f in os.listdir(_search_dir):
        if _f.startswith("test_stage") and _f.endswith(".py") and _f != "test_stage_shim.py":
            _path = os.path.join(_search_dir, _f)
            try:
                with open(_path, encoding="utf-8") as _fh:
                    _src = _fh.read()
            except OSError:
                continue
            if "def main(" in _src or 'if __name__ == "__main__":' in _src:
                collect_ignore.append(_f)


def pytest_addoption(parser):
    """Add custom command-line options."""
    parser.addoption(
        "--update-snapshots",
        action="store_true",
        default=False,
        help="Update DAG snapshot fixtures from current dry-run output",
    )


@pytest.fixture(scope="session")
def repo_root():
    """Return the repository root directory."""
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


@pytest.fixture(scope="session")
def snakefile(repo_root):
    """Return the path to the workflow Snakefile."""
    return os.path.join(repo_root, "workflow", "Snakefile")


@pytest.fixture(scope="session")
def profiles_dir(repo_root):
    """Return the directory containing test profiles."""
    return os.path.join(repo_root, "test", "profiles")


@pytest.fixture(scope="session")
def snapshots_dir(repo_root):
    """Return the directory for DAG snapshot fixtures."""
    return os.path.join(repo_root, "test", "fixtures", "dag_snapshots")


SMOKE_PROFILES = [
    "chipseq_se_noctrl",
    "chipseq_pe_noctrl",
    "chipseq_pe_ctrlsample",
    "cuttag_pe_noctrl",
    "cuttag_pe_seacr",
    "chipseq_idr_dryrun",
    "chipseq_pe_external_ctrlbam",
    "mnase_pe_noctrl",
]


@pytest.fixture(scope="session")
def smoke_profiles():
    """Return the list of smoke-test profile names."""
    return SMOKE_PROFILES


@pytest.fixture(scope="session")
def idr_paths_file(repo_root):
    """Return the path to workflow/rules/idr_paths.smk."""
    return Path(repo_root) / "workflow" / "rules" / "idr_paths.smk"


def _discover_placeholders(samples_tsv_path):
    """Parse samples.tsv and return placeholder file paths needed for dry-run."""
    paths = set()
    with open(samples_tsv_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            fq1 = (row.get("fastq_1") or "").strip()
            fq2 = (row.get("fastq_2") or "").strip()
            cb = (row.get("control_bam") or "").strip()
            if fq1:
                paths.add(fq1)
            if fq2:
                paths.add(fq2)
            if cb:
                paths.add(cb)
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


def prepare_profile_workdir(profile_dir):
    """Create a temporary workdir with placeholders and rewritten config.

    Returns (workdir, dest_config). Caller is responsible for cleaning up
    workdir when done.
    """
    workdir = tempfile.mkdtemp(prefix="profile_", dir="/tmp")
    samples_tsv_src = os.path.join(profile_dir, "samples.tsv")
    config_yaml_src = os.path.join(profile_dir, "config.yaml")

    for rel_path in _discover_placeholders(samples_tsv_src):
        placeholder_path = os.path.join(workdir, os.path.basename(rel_path))
        with open(placeholder_path, "w"):
            pass

    dest_samples = os.path.join(workdir, "samples.tsv")
    shutil.copy2(samples_tsv_src, dest_samples)

    rewritten_config = _rewrite_profile_config(config_yaml_src, workdir)
    dest_config = os.path.join(workdir, "config.yaml")
    with open(dest_config, "w") as fh:
        fh.write(rewritten_config)

    return workdir, dest_config


@pytest.fixture(scope="session")
def test_data_dir(repo_root):
    """Return the test data directory."""
    return Path(repo_root) / "test" / "data"


@pytest.fixture(scope="session")
def valid_config_path(repo_root):
    """Return the path to a valid example config YAML."""
    return Path(repo_root) / "config" / "config.yaml"


@pytest.fixture(scope="session")
def valid_samples_path(repo_root):
    """Return the path to a valid example samples TSV."""
    return Path(repo_root) / "config" / "samples.tsv"


@pytest.fixture(scope="session")
def validator_script(repo_root):
    """Return the path to scripts/validate_samples.py."""
    return os.path.join(repo_root, "scripts", "validate_samples.py")


@pytest.fixture(scope="session")
def snakemake_executable():
    """Resolve the snakemake executable used by dry-run tests."""
    sys.path.insert(0, _TEST_DIR)
    from _tool_resolver import resolve_tool

    return resolve_tool("snakemake", "SNAKEMAKE")


@pytest.fixture
def tmp_config(tmp_path):
    """Return a helper that writes a temporary config + samples TSV pair.

    Usage::

        workdir, config_path, samples_path = tmp_config(
            config={"samples": "...", "use_control": False},
            samples="sample\tfastq_1\nS1\tR1.fq\n",
        )

    The returned *workdir* is the tmp_path directory; caller need not clean up.
    """

    def _make(config, samples="", placeholders=None):
        workdir = tmp_path
        samples_path = workdir / "samples.tsv"
        samples_path.write_text(samples, encoding="utf-8")
        config_path = workdir / "config.yaml"

        resolved = dict(config)
        if "samples" not in resolved:
            resolved["samples"] = str(samples_path)

        with open(config_path, "w", encoding="utf-8") as fh:
            _write_yaml(fh, resolved)

        for name in placeholders or []:
            p = workdir / name
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_text("", encoding="utf-8")

        return workdir, str(config_path), str(samples_path)

    return _make


def _write_yaml(fh, data, indent=0):
    """Write a nested YAML mapping to a file handle."""
    prefix = "  " * indent
    for key, value in data.items():
        if isinstance(value, dict):
            fh.write(f"{prefix}{key}:\n")
            _write_yaml(fh, value, indent + 1)
        elif isinstance(value, bool):
            fh.write(f"{prefix}{key}: {str(value).lower()}\n")
        elif isinstance(value, str):
            if value:
                fh.write(f'{prefix}{key}: "{value}"\n')
            else:
                fh.write(f'{prefix}{key}: ""\n')
        elif isinstance(value, list):
            fh.write(f"{prefix}{key}:\n")
            for item in value:
                if isinstance(item, dict):
                    fh.write(f"{prefix}  -\n")
                    _write_yaml(fh, item, indent + 2)
                else:
                    fh.write(f"{prefix}  - {item}\n")
        else:
            fh.write(f"{prefix}{key}: {value}\n")


@pytest.fixture
def run_validator(validator_script):
    """Return a helper that runs validate_samples.py on a config path.

    Returns a ``Result`` object with ``rc``, ``stdout``, ``stderr`` attributes.
    """

    def _run(config_path, strict_inputs=False):
        cmd = [sys.executable, validator_script, "--config", str(config_path)]
        if strict_inputs:
            cmd.append("--strict-inputs")
        p = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
        )
        return p

    return _run


@pytest.fixture
def run_snakemake(snakemake_executable, snakefile):
    """Return a helper that runs snakemake -n for a config file.

    Returns a ``Result`` object with ``rc``, ``stdout``, ``stderr`` attributes.
    """

    def _run(config_path, extra_args=None, quiet=True):
        if extra_args:
            quiet = False
        cmd = [
            snakemake_executable,
            "-s",
            snakefile,
            "--configfile",
            str(config_path),
            "--dry-run",
        ]
        if quiet:
            cmd.append("--quiet")
        if extra_args:
            cmd.extend(extra_args)
        return subprocess.run(cmd, capture_output=True, text=True)

    return _run


def _load_idr_paths_namespace(
    idr_paths_file,
    outdir="results",
    pooled_control_experiments=None,
    sample_map=None,
    treatment_samples_by_experiment=None,
):
    """Load idr_paths.smk helpers with Snakefile globals monkey-patched.

    This lets unit tests exercise pure string helpers without importing the
    full Snakefile namespace.
    """
    code = idr_paths_file.read_text()
    namespace = {
        "OUTDIR": outdir,
        "POOLED_CONTROL_EXPERIMENTS": set(pooled_control_experiments or []),
        "SAMPLE_MAP": sample_map or {},
        "TREATMENT_SAMPLES_BY_EXPERIMENT": treatment_samples_by_experiment or {},
        "_normalize_genome": lambda genome: {"hg38": "hs", "mm10": "mm"}.get(
            genome, genome
        ),
        "_tool_param": lambda tool, key, default: {
            ("idr_macs3", "pvalue"): 0.1,
            ("idr_macs3", "extra_args"): "",
            ("macs3", "broad_cutoff"): 0.1,
        }.get((tool, key), default),
        "_bioreps_for": lambda experiment, sample_type: ["1", "2"],
    }
    exec(compile(code, str(idr_paths_file), "exec"), namespace)
    return namespace


@pytest.fixture
def idr_paths_namespace(idr_paths_file):
    """Return a namespace with the default idr_paths.smk helpers loaded."""
    return _load_idr_paths_namespace(idr_paths_file)


@pytest.fixture
def idr_paths_namespace_with_control(idr_paths_file):
    """Return a namespace where EXP1 has a pooled control experiment."""
    return _load_idr_paths_namespace(
        idr_paths_file,
        pooled_control_experiments=["EXP1"],
    )


@pytest.fixture
def idr_paths_namespace_chipseq_pe(idr_paths_file):
    """Return a namespace with a single ChIP-seq PE treatment sample."""
    return _load_idr_paths_namespace(
        idr_paths_file,
        sample_map={"S1": {"layout": "PE", "genome": "hg38"}},
        treatment_samples_by_experiment={"EXP1": ["S1"]},
    )


@pytest.fixture
def idr_paths_namespace_chipseq_se(idr_paths_file):
    """Return a namespace with a single ChIP-seq SE treatment sample."""
    return _load_idr_paths_namespace(
        idr_paths_file,
        sample_map={"S1": {"layout": "SE", "genome": "mm10"}},
        treatment_samples_by_experiment={"EXP2": ["S1"]},
    )
