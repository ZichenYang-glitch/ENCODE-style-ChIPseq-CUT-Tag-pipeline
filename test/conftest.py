"""Shared pytest fixtures for the test suite."""

import asyncio
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import runpy
import shutil
import subprocess
import sys
import tempfile
from threading import Event, Thread

_REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
_require_checkout = runpy.run_path(
    str(_REPOSITORY_ROOT / "scripts" / "source_provenance.py")
)["require_checkout"]
_require_checkout(_REPOSITORY_ROOT)

import pytest  # noqa: E402


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


@pytest.fixture
def reference_ready_app(tmp_path, monkeypatch):
    """Return an API app with one enabled, operator-prepared tiny ENCODE profile."""
    from encode_pipeline.api.main import create_app

    reference_root = tmp_path / "operator-reference"
    resources = {}
    for name, suffix, content in (
        ("reference_fasta", "fa", b">chr1\nACGT\n"),
        ("gtf", "gtf", b'chr1\ttest\texon\t1\t4\t.\t+\t.\tgene_id "g1";\n'),
        ("chrom_sizes", "sizes", b"chr1\t4\n"),
        ("blacklist", "bed", b"chr1\t1\t2\n"),
    ):
        path = reference_root / f"{name}.{suffix}"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(content)
        resources[name] = {
            "path": str(path),
            "sha256": hashlib.sha256(content).hexdigest(),
        }
    prefix = reference_root / "bowtie2" / "GRCh38"
    index_files = {}
    for suffix in (
        ".1.bt2",
        ".2.bt2",
        ".3.bt2",
        ".4.bt2",
        ".rev.1.bt2",
        ".rev.2.bt2",
    ):
        path = Path(f"{prefix}{suffix}")
        path.parent.mkdir(parents=True, exist_ok=True)
        content = suffix.encode("ascii")
        path.write_bytes(content)
        index_files[suffix] = hashlib.sha256(content).hexdigest()
    binding = {
        "schema_version": "encode-reference-binding-v1",
        "assembly": "GRCh38",
        "effective_genome_size": "hs",
        "genome_resources": resources,
        "bowtie2_index": {"prefix": str(prefix), "files": index_files},
    }
    private_config = tmp_path / "operator-reference-profiles.json"
    private_config.write_text(
        json.dumps(
            {
                "schema_version": "helixweave-reference-profiles-v1",
                "profiles": {
                    "grch38-test": {
                        "bindings": {"encode-style-chipseq-cuttag-atac-mnase": binding}
                    }
                },
            }
        ),
        encoding="utf-8",
    )
    monkeypatch.setenv(
        "ENCODE_PIPELINE_REFERENCE_PROFILE_CONFIG",
        str(private_config),
    )
    app = create_app()
    summary = app.state.reference_profile_service.register(
        safe_key="grch38-test",
        display_name="GRCh38 tiny",
        organism="Homo sapiens",
        assembly="GRCh38",
        config_key="grch38-test",
    )
    enabled = app.state.reference_profile_service.enable(
        summary.profile_id,
        revision_id=summary.revision_id,
    )
    app.state.test_reference_profile = enabled
    seed_test_authentication(app)
    return app


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


@pytest.fixture
def bulk_rnaseq_qualifications():
    """Return explicit synthetic qualification for composition-only tests."""
    from encode_pipeline.adapters.bulk_rnaseq.execution_identity import (
        ExecutionImplementationQualification,
    )

    return {
        "implementation_qualification": ExecutionImplementationQualification(
            manifest_sha256="1" * 64,
            aggregate_sha256="2" * 64,
            file_count=1,
            persistence_contract_version="1.0.0",
            persistence_contract_sha256="3" * 64,
        ),
    }


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


TEST_AUTH_PASSWORD_HASH = (
    "$argon2id$v=19$m=65536,t=3,p=4$"
    "pGLQVK/BOQLC0oPSA8RQTg$TdiQWUP9gwBAI8iAXUT7oEtjOPqKjJhqyW0JS8ye/Ag"
)
TEST_AUTH_SESSION_TOKEN = "a" * 43
TEST_AUTH_CSRF_TOKEN = "b" * 43
TEST_AUTH_USER_ID = "usr_" + "0" * 32


def seed_test_authentication(app, monkeypatch=None):
    """Seed one administrator and session so API tests keep their own focus."""
    from datetime import datetime, timedelta, timezone

    from encode_pipeline.platform.authentication import (
        UserAccount,
        UserRole,
        UserStatus,
    )
    from encode_pipeline.services.authentication import (
        SessionSecrets,
        new_session_record,
    )

    repository = app.state.authentication_repository
    now = datetime.now(timezone.utc)
    account = UserAccount(
        user_id=TEST_AUTH_USER_ID,
        username="test-admin",
        role=UserRole.ADMINISTRATOR,
        status=UserStatus.ENABLED,
        password_hash=TEST_AUTH_PASSWORD_HASH,
        created_at=now,
        updated_at=now,
        password_changed_at=now,
    )
    from encode_pipeline.services.authentication_repositories import (
        AuthenticationAccountConflictError,
    )

    try:
        repository.create_account(account)
    except AuthenticationAccountConflictError:
        pass
    secrets = SessionSecrets._from_generated(
        TEST_AUTH_SESSION_TOKEN,
        TEST_AUTH_CSRF_TOKEN,
    )
    record = new_session_record(
        user_id=account.user_id,
        secrets=secrets,
        created_at=now,
        lifetime=timedelta(hours=8),
    )
    try:
        repository.get_session(record.session_digest)
    except KeyError:
        repository.create_session(record)
    app.state.test_auth_tokens = (TEST_AUTH_SESSION_TOKEN, TEST_AUTH_CSRF_TOKEN)
    return account
