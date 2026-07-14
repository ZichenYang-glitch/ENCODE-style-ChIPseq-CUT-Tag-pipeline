"""Shared pytest fixtures for the test suite."""

import asyncio
import csv
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from threading import Event, Thread

import pytest


_TEST_DIR = os.path.dirname(os.path.abspath(__file__))


def pytest_addoption(parser):
    """Add custom command-line options."""
    parser.addoption(
        "--update-snapshots",
        action="store_true",
        default=False,
        help="Update DAG snapshot fixtures from current dry-run output",
    )


@pytest.fixture(autouse=True)
def isolate_api_database(request, monkeypatch):
    """Give API tests a disposable file-backed database without shadowing conftest."""
    api_test_dir = Path(__file__).resolve().parent / "api"
    if api_test_dir not in Path(request.fspath).resolve().parents:
        yield
        return

    tmp_path = request.getfixturevalue("tmp_path")
    monkeypatch.setenv(
        "ENCODE_PIPELINE_DATABASE_URL",
        f"sqlite:///{tmp_path / 'platform.db'}",
    )
    yield


async def _run_in_joined_test_thread(function, *args, **kwargs):
    """Run sync API routes without leaking this environment's executor threads."""
    completed = Event()
    results = []
    exceptions = []

    def invoke():
        try:
            results.append(function(*args, **kwargs))
        except BaseException as exc:
            exceptions.append(exc)
        finally:
            completed.set()

    thread = Thread(target=invoke)
    thread.start()
    try:
        while not completed.is_set():
            await asyncio.sleep(0.001)
    finally:
        thread.join(timeout=3)
        if thread.is_alive():  # pragma: no cover - test deadlock guard
            raise RuntimeError("test threadpool call did not terminate")
    if exceptions:
        raise exceptions[0]
    return results[0]


@pytest.fixture(autouse=True)
def joined_api_test_threadpool(request, monkeypatch):
    """Exercise synchronous FastAPI routes with a deterministic joined thread."""
    api_test_dir = Path(__file__).resolve().parent / "api"
    if api_test_dir not in Path(request.fspath).resolve().parents:
        yield
        return

    import fastapi.routing

    monkeypatch.setattr(
        fastapi.routing,
        "run_in_threadpool",
        _run_in_joined_test_thread,
    )
    yield


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
def run_validator():
    """Return a helper that runs the config validator CLI entry point.

    Returns a ``Result`` object with ``returncode``, ``stdout``, ``stderr``
    attributes, matching the subprocess interface used by legacy tests.
    """
    from encode_pipeline.config import validator

    class _Result:
        def __init__(self, returncode, stdout, stderr):
            self.returncode = returncode
            self.stdout = stdout
            self.stderr = stderr

    def _run(config_path, strict_inputs=False):
        import io

        old_stdout = sys.stdout
        old_stderr = sys.stderr
        old_argv = sys.argv
        stdout_capture = io.StringIO()
        stderr_capture = io.StringIO()
        sys.stdout = stdout_capture
        sys.stderr = stderr_capture
        sys.argv = ["validate_samples.py", "--config", str(config_path)]
        if strict_inputs:
            sys.argv.append("--strict-inputs")
        try:
            validator.main()
            return _Result(0, stdout_capture.getvalue(), stderr_capture.getvalue())
        except SystemExit as exc:
            code = exc.code if exc.code is not None else 0
            return _Result(code, stdout_capture.getvalue(), stderr_capture.getvalue())
        finally:
            sys.stdout = old_stdout
            sys.stderr = old_stderr
            sys.argv = old_argv

    return _run


@pytest.fixture
def run_snakemake(snakemake_executable, snakefile):
    """Return a helper that runs snakemake -n for a config file.

    Returns a ``Result`` object with ``rc``, ``stdout``, ``stderr`` attributes.
    """

    def _run(config_path, extra_args=None, quiet=True):
        config_path = Path(config_path).resolve()
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
        env = os.environ.copy()
        env["PYTHONDONTWRITEBYTECODE"] = "1"
        env["XDG_CACHE_HOME"] = str(config_path.parent / ".cache")
        return subprocess.run(
            cmd,
            cwd=config_path.parent,
            capture_output=True,
            text=True,
            env=env,
        )

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
