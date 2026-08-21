#!/usr/bin/env python3
"""Prepare a controlled tiny workflow and exec the real local platform stack."""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from contextlib import redirect_stderr, redirect_stdout
from copy import deepcopy
import csv
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import shutil
import stat
import sys
from uuid import uuid4

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from encode_pipeline.cli.results_visibility_fixture import (  # noqa: E402
    ResultsVisibilityInputs,
    prepare_results_visibility_fixture,
)
from encode_pipeline.cli import app as helixweave_cli  # noqa: E402
from encode_pipeline.services.private_reference_profiles import (  # noqa: E402
    PRIVATE_REFERENCE_PROFILE_SCHEMA_VERSION,
)
from encode_pipeline.workers.settings import (  # noqa: E402
    REFERENCE_PROFILE_CONFIG_ENV,
)
from bulk_product_runtime import (  # noqa: E402
    prepare_bulk_product_browser_runtime,
)


OWNERSHIP_SENTINEL = ".encode-platform-playwright-owned"
PRODUCT_BULK_GATE_ENV = "HELIXWEAVE_REQUIRE_BULK_RNASEQ_PRODUCT_GATE"
ENCODE_WORKFLOW_ID = "encode-style-chipseq-cuttag-atac-mnase"
BULK_WORKFLOW_ID = "bulk-rnaseq"
_PRIVATE_REFERENCE_FIELDS = frozenset({"genome", "bowtie2_index"})
_BOWTIE2_SUFFIXES = (
    ".1.bt2",
    ".2.bt2",
    ".3.bt2",
    ".4.bt2",
    ".rev.1.bt2",
    ".rev.2.bt2",
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare_bulk_authoring_fixture(runtime_root: Path) -> dict[str, object]:
    """Create tiny synthetic authoring inputs that are never used for execution."""
    fixture_root = runtime_root / "bulk-authoring"
    reads_root = fixture_root / "reads"
    rrna_root = fixture_root / "rrna"
    sortmerna_index_root = rrna_root / "sortmerna-index"
    for directory in (reads_root, sortmerna_index_root):
        directory.mkdir(parents=True, exist_ok=True)

    read_names = (
        "pairedA_1.lib1.L001.fastq.gz",
        "pairedA_2.lib1.L001.fastq.gz",
        "pairedA_1.lib1.L002.fastq.gz",
        "pairedA_2.lib1.L002.fastq.gz",
        "singleB.lib1.L001.fastq.gz",
    )
    synthetic_fastq = b"@synthetic-read\nACGTACGT\n+\nIIIIIIII\n"
    for read_name in read_names:
        (reads_root / read_name).write_bytes(gzip.compress(synthetic_fastq, mtime=0))

    rrna_fasta_path = rrna_root / "synthetic-rrna.fa"
    rrna_fasta_path.write_text(">rrnaSynthetic\nACGTACGT\n", encoding="utf-8")
    rrna_manifest_path = rrna_root / "database-manifest.json"
    rrna_manifest_path.write_text(
        json.dumps(
            {
                "database": rrna_fasta_path.name,
                "sha256": _sha256(rrna_fasta_path),
            },
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    sortmerna_manifest_path = sortmerna_index_root / "identity.json"
    sortmerna_manifest_path.write_text(
        json.dumps(
            {
                "database_manifest_sha256": _sha256(rrna_manifest_path),
                "fixture": "authoring-only",
            },
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )

    samples_path = fixture_root / "samples.tsv"
    header = "sample\tlibrary\tlane\tlayout\tfastq_1\tfastq_2\tstrandedness\tplatform"
    rows = (
        (
            "pairedA",
            "lib1",
            "L001",
            "PE",
            reads_root / read_names[0],
            reads_root / read_names[1],
            "reverse",
            "ILLUMINA",
        ),
        (
            "pairedA",
            "lib1",
            "L002",
            "PE",
            reads_root / read_names[2],
            reads_root / read_names[3],
            "reverse",
            "ILLUMINA",
        ),
        (
            "singleB",
            "lib1",
            "L001",
            "SE",
            reads_root / read_names[4],
            "",
            "auto",
            "ILLUMINA",
        ),
    )
    samples_path.write_text(
        "\n".join([header, *("\t".join(map(str, row)) for row in rows)]) + "\n",
        encoding="utf-8",
    )

    config = {
        "standard": {
            "analysis": {"alignment": "star", "quantification": "salmon"},
            "trimming": {"enabled": True, "tool": "fastp"},
            "ribosomal_rna_removal": {
                "enabled": True,
                "tool": "sortmerna",
                "save_filtered_reads": True,
                "database_manifest": {
                    "path": str(rrna_manifest_path),
                    "identity_sha256": _sha256(rrna_manifest_path),
                },
                "sortmerna_index": {
                    "path": str(sortmerna_index_root),
                    "identity_sha256": _sha256(sortmerna_manifest_path),
                },
            },
            "umi": {
                "enabled": True,
                "mode": "read_name",
                "deduplication_tool": "umitools",
                "read_name_separator": ":",
                "grouping_method": "directional",
                "emit_dedup_stats": False,
                "primary_alignments_only": False,
            },
            "outputs": {"trimmed_reads": True, "umi_intermediates": True},
        },
        "advanced": {"min_trimmed_reads": 1},
    }
    return {
        "bulkWorkflowId": "bulk-rnaseq",
        "bulkSamplesPath": str(samples_path),
        "bulkConfig": config,
        "bulkOptions": {},
        "bulkExpectedExecution": "not_configured",
        "bulkRequiredArtifactOutputTypes": [],
        "bulkRequiredQcMetricKeys": [],
    }


def prepare_bulk_browser_fixture(
    runtime_root: Path,
    environ: dict[str, str] | None = None,
) -> tuple[dict[str, object], dict[str, str], Mapping[str, object] | None]:
    """Select ordinary authoring or the explicitly protected product fixture."""
    source = dict(os.environ if environ is None else environ)
    requested = source.get(PRODUCT_BULK_GATE_ENV)
    if requested is None:
        return prepare_bulk_authoring_fixture(runtime_root), {}, None
    if requested != "1":
        raise ValueError(f"{PRODUCT_BULK_GATE_ENV} must be exactly 1 when configured")
    prepared = prepare_bulk_product_browser_runtime(runtime_root, source)
    return (
        dict(prepared.manifest_fields),
        dict(prepared.deployment_environment),
        prepared.reference_profile_binding,
    )


def prepare_reference_profile_runtime(
    runtime_root: Path,
    inputs: ResultsVisibilityInputs,
    *,
    bulk_binding: Mapping[str, object] | None,
    cli_command: Callable[[Sequence[str] | None], int] = helixweave_cli.main,
) -> tuple[ResultsVisibilityInputs, dict[str, str]]:
    """Register exact test-private references through ``helixweave admin``."""
    public_inputs = _strip_encode_reference_fields(inputs)
    operator_root = runtime_root / "operator-references"
    operator_root.mkdir(mode=0o700)
    configured_profiles: dict[str, object] = {}
    registrations = [
        (
            "browser-grch38",
            "GRCh38 browser tiny",
            "Homo sapiens",
            "GRCh38",
            "browser-grch38-private",
            ENCODE_WORKFLOW_ID,
            _write_encode_reference_binding(operator_root / "grch38", "GRCh38", "hs"),
        ),
        (
            "browser-mm10",
            "mm10 browser tiny",
            "Mus musculus",
            "mm10",
            "browser-mm10-private",
            ENCODE_WORKFLOW_ID,
            _write_encode_reference_binding(operator_root / "mm10", "mm10", "mm"),
        ),
    ]
    if bulk_binding is not None:
        registrations.append(
            (
                "browser-bulk-rnaseq",
                "Bulk RNA-seq browser tiny",
                "Synthetic organism",
                "tiny",
                "browser-bulk-private",
                BULK_WORKFLOW_ID,
                dict(bulk_binding),
            )
        )
    for _, _, _, _, config_key, workflow_id, binding in registrations:
        configured_profiles[config_key] = {"bindings": {workflow_id: binding}}

    config_path = operator_root / "reference-profiles.json"
    config_descriptor = os.open(
        config_path,
        os.O_WRONLY
        | os.O_CREAT
        | os.O_EXCL
        | getattr(os, "O_CLOEXEC", 0)
        | getattr(os, "O_NOFOLLOW", 0),
        0o600,
    )
    with os.fdopen(config_descriptor, "w", encoding="utf-8") as handle:
        json.dump(
            {
                "schema_version": PRIVATE_REFERENCE_PROFILE_SCHEMA_VERSION,
                "profiles": configured_profiles,
            },
            handle,
            ensure_ascii=True,
            sort_keys=True,
            separators=(",", ":"),
        )
        handle.write("\n")

    database_url = f"sqlite:///{runtime_root / 'platform.db'}"
    for safe_key, display_name, organism, assembly, config_key, _, _ in registrations:
        registered = _run_reference_profile_admin(
            [
                "--database-url",
                database_url,
                "--reference-profile-config",
                str(config_path),
                "reference-profile",
                "register",
                "--safe-key",
                safe_key,
                "--display-name",
                display_name,
                "--organism",
                organism,
                "--assembly",
                assembly,
                "--config-key",
                config_key,
            ],
            cli_command=cli_command,
        )
        profile_id = registered.get("profile_id")
        revision_id = registered.get("revision_id")
        if not isinstance(profile_id, str) or not isinstance(revision_id, str):
            raise RuntimeError(
                "reference profile registration returned invalid identity"
            )
        _run_reference_profile_admin(
            [
                "--database-url",
                database_url,
                "--reference-profile-config",
                str(config_path),
                "reference-profile",
                "verify",
                revision_id,
            ],
            cli_command=cli_command,
        )
        _run_reference_profile_admin(
            [
                "--database-url",
                database_url,
                "--reference-profile-config",
                str(config_path),
                "reference-profile",
                "enable",
                profile_id,
                "--revision-id",
                revision_id,
            ],
            cli_command=cli_command,
        )
    return public_inputs, {REFERENCE_PROFILE_CONFIG_ENV: str(config_path)}


def _run_reference_profile_admin(
    arguments: Sequence[str],
    *,
    cli_command: Callable[[Sequence[str] | None], int],
) -> dict[str, object]:
    stdout = io.StringIO()
    stderr = io.StringIO()
    with redirect_stdout(stdout), redirect_stderr(stderr):
        exit_code = cli_command(["admin", *arguments])
    if exit_code != 0:
        raise RuntimeError("reference profile administrator command failed")
    try:
        result = json.loads(stdout.getvalue())
    except (TypeError, ValueError):
        raise RuntimeError(
            "reference profile administrator output is invalid"
        ) from None
    if not isinstance(result, dict) or stderr.getvalue():
        raise RuntimeError("reference profile administrator output is invalid")
    return result


def _strip_encode_reference_fields(
    inputs: ResultsVisibilityInputs,
) -> ResultsVisibilityInputs:
    with inputs.samples_path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = tuple(
            field
            for field in (reader.fieldnames or ())
            if field not in _PRIVATE_REFERENCE_FIELDS
        )
        rows = [
            {key: value for key, value in row.items() if key in fieldnames}
            for row in reader
        ]
    if not rows or not fieldnames:
        raise ValueError("controlled ENCODE sample fixture is invalid")
    temporary = inputs.samples_path.with_suffix(".public.tmp")
    with temporary.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(inputs.samples_path)

    def public_config(value: Mapping[str, object]) -> dict[str, object]:
        result = deepcopy(dict(value))
        result.pop("genome_resources", None)
        return result

    return ResultsVisibilityInputs(
        samples_path=inputs.samples_path,
        results_config=public_config(inputs.results_config),
        cancel_config=public_config(inputs.cancel_config),
        empty_config=public_config(inputs.empty_config),
        malformed_config=public_config(inputs.malformed_config),
        expected_qc_summary=inputs.expected_qc_summary,
    )


def _write_encode_reference_binding(
    root: Path,
    assembly: str,
    effective_genome_size: str,
) -> dict[str, object]:
    root.mkdir(parents=True)
    resources: dict[str, object] = {}
    for name, suffix, content in (
        ("reference_fasta", "fa", f">{assembly}\nACGTACGT\n".encode()),
        (
            "gtf",
            "gtf",
            (
                f"{assembly}\tHelixWeave\texon\t1\t8\t.\t+\t.\t"
                'gene_id "g1"; transcript_id "t1";\n'
            ).encode(),
        ),
        ("chrom_sizes", "sizes", f"{assembly}\t8\n".encode()),
        ("blacklist", "bed", f"{assembly}\t1\t2\n".encode()),
    ):
        path = root / f"{name}.{suffix}"
        path.write_bytes(content)
        resources[name] = {"path": str(path), "sha256": _sha256(path)}
    index_prefix = root / "bowtie2" / assembly
    index_prefix.parent.mkdir()
    index_files: dict[str, str] = {}
    for suffix in _BOWTIE2_SUFFIXES:
        path = Path(f"{index_prefix}{suffix}")
        path.write_bytes(f"{assembly}{suffix}\n".encode())
        index_files[suffix] = _sha256(path)
    return {
        "schema_version": "encode-reference-binding-v1",
        "assembly": assembly,
        "effective_genome_size": effective_genome_size,
        "genome_resources": resources,
        "bowtie2_index": {"prefix": str(index_prefix), "files": index_files},
    }


def write_manifest(
    runtime_root: Path,
    queue_name: str,
    inputs: ResultsVisibilityInputs,
    bulk_authoring: dict[str, object],
) -> None:
    """Write the bounded browser-visible fixture contract."""
    manifest = {
        "workflowId": "encode-style-chipseq-cuttag-atac-mnase",
        "samplesPath": str(inputs.samples_path),
        "resultsConfig": inputs.results_config,
        "cancelConfig": inputs.cancel_config,
        "emptyConfig": inputs.empty_config,
        "malformedConfig": inputs.malformed_config,
        "expectedQcSummary": inputs.expected_qc_summary,
        "runtimeRoot": str(runtime_root),
        "workspaceRoot": str(runtime_root / "workspaces"),
        "markerRoot": str(runtime_root / "tmp"),
        "queueName": queue_name,
    }
    manifest.update(bulk_authoring)
    (runtime_root / "runtime.json").write_text(
        json.dumps(manifest, sort_keys=True), encoding="utf-8"
    )


def _validated_owned_runtime_root(
    runtime_value: str, runtime_owner: str
) -> tuple[Path, int]:
    """Resolve one invocation-owned runtime without following its root entry."""
    requested_root = Path(runtime_value)
    if not requested_root.is_absolute():
        raise ValueError("E2E runtime root must be an absolute path")
    lexical_root = Path(os.path.abspath(runtime_value))
    runtime_root = lexical_root.resolve(strict=True)
    if runtime_root != lexical_root:
        raise ValueError("E2E runtime root must not traverse symbolic links")
    if runtime_root in {
        Path("/"),
        Path("/tmp").resolve(),
        Path.home().resolve(),
        REPOSITORY_ROOT.resolve(),
    }:
        raise ValueError("E2E runtime root must be a dedicated temporary directory")
    runtime_info = runtime_root.lstat()
    if not stat.S_ISDIR(runtime_info.st_mode) or runtime_info.st_uid != os.geteuid():
        raise ValueError("E2E runtime root is not an owned directory")
    sentinel = runtime_root / OWNERSHIP_SENTINEL
    sentinel_info = sentinel.lstat()
    if (
        not stat.S_ISREG(sentinel_info.st_mode)
        or sentinel_info.st_nlink != 1
        or sentinel_info.st_uid != os.geteuid()
        or sentinel.read_text(encoding="utf-8") != runtime_owner
    ):
        raise ValueError("E2E runtime root is not owned by this Playwright invocation")
    return runtime_root, runtime_info.st_dev


def _make_owned_directories_removable(path: Path, expected_device: int) -> None:
    """Restore owner access on directories without following links or mounts."""
    info = path.lstat()
    if (
        not stat.S_ISDIR(info.st_mode)
        or info.st_uid != os.geteuid()
        or info.st_dev != expected_device
        or path.is_mount()
    ):
        raise ValueError("E2E runtime contains an unsafe directory")
    path.chmod(info.st_mode | stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR)
    with os.scandir(path) as entries:
        children = tuple(entries)
    for child in children:
        child_info = child.stat(follow_symlinks=False)
        if child_info.st_uid != os.geteuid():
            raise ValueError("E2E runtime contains an unowned entry")
        if stat.S_ISLNK(child_info.st_mode):
            continue
        child_path = Path(child.path)
        if child_path.is_mount():
            raise ValueError("E2E runtime contains an unexpected mount")
        if stat.S_ISDIR(child_info.st_mode):
            _make_owned_directories_removable(child_path, expected_device)


def cleanup_owned_runtime_root(runtime_value: str, runtime_owner: str) -> None:
    """Remove only a sentinel-owned runtime, including copied read-only trees."""
    runtime_root, expected_device = _validated_owned_runtime_root(
        runtime_value, runtime_owner
    )
    _make_owned_directories_removable(runtime_root, expected_device)
    shutil.rmtree(runtime_root)
    if os.path.lexists(runtime_root):
        raise OSError("E2E runtime root still exists after cleanup")


def prepare_owned_runtime_root(runtime_value: str, runtime_owner: str) -> Path:
    """Reset only the invocation-owned temporary runtime directory."""
    runtime_root = Path(runtime_value).resolve()
    cleanup_owned_runtime_root(runtime_value, runtime_owner)
    runtime_root.mkdir(mode=0o700)
    (runtime_root / OWNERSHIP_SENTINEL).write_text(runtime_owner, encoding="utf-8")
    return runtime_root


E2E_ADMIN_USERNAME = "e2e-admin"
E2E_ADMIN_PASSWORD = "e2e playwright admin password"


def prepare_auth_admin(runtime_root: Path) -> None:
    """Create the browser runtime's unique administrator before services start."""
    from encode_pipeline.persistence.database import (
        create_database_engine,
        create_session_factory,
    )
    from encode_pipeline.persistence.authentication import (
        SqlAlchemyAuthenticationRepository,
    )
    from encode_pipeline.persistence.migrations import upgrade_database
    from encode_pipeline.services.authentication_service import (
        AccountAdministrationService,
    )

    database_url = f"sqlite:///{runtime_root / 'platform.db'}"
    upgrade_database(database_url)
    engine = create_database_engine(database_url)
    try:
        repository = SqlAlchemyAuthenticationRepository(create_session_factory(engine))
        AccountAdministrationService(
            repository=repository
        ).bootstrap_initial_administrator(
            E2E_ADMIN_USERNAME,
            E2E_ADMIN_PASSWORD,
        )
    finally:
        engine.dispose()


def main() -> None:
    runtime_value = os.environ.get("ENCODE_PIPELINE_E2E_ROOT")
    runtime_owner = os.environ.get("ENCODE_PIPELINE_E2E_OWNER")
    if not runtime_value or not runtime_owner:
        raise ValueError("Playwright runtime root and owner token are required")
    runtime_root = prepare_owned_runtime_root(runtime_value, runtime_owner)
    project_root = runtime_root / "project"
    inputs = prepare_results_visibility_fixture(project_root)
    (
        bulk_authoring,
        deployment_environment,
        bulk_reference_binding,
    ) = prepare_bulk_browser_fixture(runtime_root)
    inputs, reference_environment = prepare_reference_profile_runtime(
        runtime_root,
        inputs,
        bulk_binding=bulk_reference_binding,
    )
    deployment_environment.update(reference_environment)
    queue_name = f"encode-pipeline-browser-{uuid4().hex}"
    prepare_auth_admin(runtime_root)
    write_manifest(runtime_root, queue_name, inputs, bulk_authoring)
    redis_url = os.environ.get(
        "ENCODE_PIPELINE_E2E_REDIS_URL", "redis://127.0.0.1:6380/0"
    )
    argv = [
        sys.executable,
        str(REPOSITORY_ROOT / "scripts" / "run_local_platform.py"),
        "--project-root",
        str(project_root),
        "--frontend-root",
        str(REPOSITORY_ROOT / "frontend"),
        "--runtime-root",
        str(runtime_root),
        "--redis-url",
        redis_url,
        "--queue-name",
        queue_name,
        "--api-port",
        "8010",
        "--frontend-port",
        "4173",
        "--readiness-timeout",
        "300" if deployment_environment else "120",
    ]
    child_environment = os.environ.copy()
    child_environment.update(deployment_environment)
    os.execvpe(sys.executable, argv, child_environment)


def _cleanup_owned_runtime_from_environment() -> int:
    runtime_value = os.environ.get("ENCODE_PIPELINE_E2E_ROOT")
    runtime_owner = os.environ.get("ENCODE_PIPELINE_E2E_OWNER")
    if not runtime_value or not runtime_owner:
        print("owned Playwright runtime cleanup was not configured", file=sys.stderr)
        return 2
    try:
        cleanup_owned_runtime_root(runtime_value, runtime_owner)
    except (OSError, ValueError):
        print("owned Playwright runtime cleanup failed", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    if sys.argv[1:] == ["--cleanup-owned-runtime"]:
        raise SystemExit(_cleanup_owned_runtime_from_environment())
    main()
