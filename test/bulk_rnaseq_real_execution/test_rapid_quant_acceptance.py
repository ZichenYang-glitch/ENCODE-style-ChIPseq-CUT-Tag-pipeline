"""Explicit tiny rapid_quant qualification through the production runtime."""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
import gzip
from hashlib import sha256
import json
import os
from pathlib import Path
import re
import signal
import stat
import time

import pytest

from encode_pipeline.adapters.bulk_rnaseq import (
    BulkRnaSeqExecutionBinding,
    BulkRnaSeqRapidQuantQualificationAdapter,
    BulkRnaSeqTranscriptomeBinding,
    RuntimeAssetBinding,
)
from encode_pipeline.adapters.bulk_rnaseq.qualification import (
    RAPID_QUANT_MODE_OWNED_PARAMETERS,
    RAPID_QUANT_PROFILE_OWNED_PARAMETERS,
)
from encode_pipeline.platform.adapters import CommandSpec, WorkflowInputs
from encode_pipeline.platform.managed_containers import MANAGED_CONTAINER_SCOPE_LABEL
from encode_pipeline.platform.results import Issue
from encode_pipeline.services.managed_containers import ManagedContainerCleaner
from encode_pipeline.services.materialization import WorkspaceMaterializer
from encode_pipeline.services.process_runner import ProcessResult, ProcessRunner

from .support import (
    AcceptanceFixture,
    GateSettings,
    assert_no_managed_containers,
    load_acceptance_fixture,
    require_gate_settings,
)


_RAPID_PROFILE = "rapid_quant"
_RAPID_MODE = "rapid-quant-v1"
_EXECUTION_TIMEOUT_SECONDS = 3_600
_MAX_GATE_FILE_BYTES = 64 * 1024 * 1024
_SHA256_IMAGE = re.compile(r"sha256:[0-9a-f]{64}")


@dataclass(frozen=True)
class RapidQuantRun:
    """Private, path-free evidence from one bounded qualification attempt."""

    process_result: ProcessResult
    process_issues: tuple[Issue, ...]
    quant_surface: tuple[tuple[str, str], ...]
    trace_processes: tuple[str, ...]
    container_script_count: int


def test_rapid_quant_input_is_single_sample_minimal_and_does_not_mutate_fixture():
    fixture = AcceptanceFixture(
        workflow_inputs=WorkflowInputs(
            config={
                "standard": {
                    "reference": {"reference_id": "tiny"},
                    "trimming": {"enabled": True, "tool": "trimgalore"},
                    "ribosomal_rna_removal": {
                        "enabled": True,
                        "tool": "sortmerna",
                    },
                }
            },
            samples=[
                {
                    "sample": "PE1",
                    "library": "pe",
                    "lane": "L001",
                    "layout": "PE",
                    "fastq_1": "/fixture/pe-r1.fastq.gz",
                    "fastq_2": "/fixture/pe-r2.fastq.gz",
                    "strandedness": "auto",
                    "platform": "ILLUMINA",
                },
                {
                    "sample": "SE1",
                    "library": "se",
                    "lane": "L001",
                    "layout": "SE",
                    "fastq_1": "/fixture/se.fastq.gz",
                    "strandedness": "auto",
                    "platform": "ILLUMINA",
                },
            ],
            options={},
        ),
        transcriptome=BulkRnaSeqTranscriptomeBinding(
            reference_id="tiny",
            fasta_sha256="a" * 64,
            gtf_sha256="b" * 64,
            transcript_fasta=Path("/fixture/transcripts.fa"),
            transcript_fasta_sha256="c" * 64,
        ),
        acceptance_manifest_sha256="d" * 64,
        source_manifest_sha256="e" * 64,
        source_identity_sha256="f" * 64,
        index_provenance_manifest_sha256="0" * 64,
        index_provenance_identity_sha256="1" * 64,
        required_artifact_output_types=(),
        required_qc_metric_keys=(),
        required_sample_ids=("PE1", "SE1"),
        required_artifact_sample_output_types=(),
        required_qc_sample_metric_keys=(),
        required_qc_sample_metric_values=(),
    )
    before = fixture.workflow_inputs.to_dict()

    inputs = _rapid_quant_inputs(fixture)

    assert inputs.samples == [before["samples"][1]]
    assert inputs.config["standard"]["trimming"] == {
        "enabled": False,
        "tool": "trimgalore",
    }
    assert inputs.config["standard"]["ribosomal_rna_removal"] == {"enabled": False}
    assert fixture.workflow_inputs.to_dict() == before


def test_rapid_quant_command_contract_is_fail_closed(tmp_path: Path):
    workspace = (tmp_path / "workspace").resolve()
    command = _unit_command(workspace)

    _assert_rapid_command_contract(command, workspace)

    argv = list(command.argv)
    argv[argv.index(_RAPID_PROFILE)] = "test"
    weakened = CommandSpec(
        argv=tuple(argv),
        cwd=command.cwd,
        env=command.env,
        preflight_argv=command.preflight_argv,
        managed_container_scope=command.managed_container_scope,
        managed_container_endpoint_identity=(
            command.managed_container_endpoint_identity
        ),
    )
    with pytest.raises(AssertionError, match="rapid_quant"):
        _assert_rapid_command_contract(weakened, workspace)


def test_rapid_quant_surface_uses_relative_deterministic_file_identity(
    tmp_path: Path,
):
    first = (tmp_path / "first").resolve()
    second = (tmp_path / "second").resolve()
    for workspace in (first, second):
        quant = workspace / "results/SE1/salmon/SE1/quant.sf"
        quant.parent.mkdir(parents=True)
        quant.write_bytes(
            b"Name\tLength\tEffectiveLength\tTPM\tNumReads\nTX\t10\t5\t1\t1\n"
        )

    assert _quant_surface(first) == _quant_surface(second)
    assert _quant_surface(first) == (
        (
            "SE1/salmon/SE1/quant.sf",
            "35cefed904a582162febaccf7c3717a0d6e99d83a0ee27228202f971d4ebced9",
        ),
    )


@pytest.mark.bulk_rnaseq_real_execution
def test_controlled_tiny_rapid_quant_nextflow_container_gate(
    tmp_path: Path,
) -> None:
    """Prove execution mechanics only; the synthetic fixture is not biology."""
    settings = require_gate_settings()
    fixture = load_acceptance_fixture(settings.fixture_manifest)
    inputs = _rapid_quant_inputs(fixture)
    binding = BulkRnaSeqExecutionBinding(
        assets=RuntimeAssetBinding(
            root=settings.runtime_root,
            docker_executable=settings.docker_executable,
            docker_socket=settings.docker_socket,
        ),
        transcriptome=fixture.transcriptome,
    )
    adapter = BulkRnaSeqRapidQuantQualificationAdapter(execution=binding)

    doctor = adapter.doctor()
    assert doctor.ready is True
    assert doctor.issues == ()
    build_identity = adapter.capture_build_identity()
    assert build_identity.is_success

    successful = tuple(
        _execute_rapid_quant(
            adapter=adapter,
            inputs=inputs,
            settings=settings,
            workspace=(tmp_path / f"success-{index}").resolve(),
        )
        for index in (1, 2)
    )
    for run in successful:
        assert run.process_result.exit_code == 0
        assert not run.process_issues
        assert tuple(path for path, _digest in run.quant_surface) == (
            "SE1/salmon/SE1/quant.sf",
        )
        assert any(name.endswith(":FASTQC") for name in run.trace_processes)
        assert any(name.endswith(":SALMON_QUANT") for name in run.trace_processes)
        assert not any(name.endswith(":STAR_ALIGN") for name in run.trace_processes)
        assert not any(
            name.endswith((":MAKE_TRANSCRIPTS_FASTA", ":RSEM_PREPAREREFERENCE"))
            for name in run.trace_processes
        )
        assert run.container_script_count >= 2
    assert successful[0].quant_surface == successful[1].quant_surface
    assert successful[0].trace_processes == successful[1].trace_processes

    broken_inputs = _broken_fastq_inputs(
        inputs,
        destination=(tmp_path / "broken.fastq.gz").resolve(),
    )
    failed = _execute_rapid_quant(
        adapter=adapter,
        inputs=broken_inputs,
        settings=settings,
        workspace=(tmp_path / "failure").resolve(),
        expect_success=False,
    )
    assert failed.process_result.exit_code != 0
    assert [issue.code for issue in failed.process_issues] == [
        "PROCESS_RUNNER_NONZERO_EXIT"
    ]
    assert failed.container_script_count >= 1


def _rapid_quant_inputs(fixture: AcceptanceFixture) -> WorkflowInputs:
    """Select the shortest controlled SE route while retaining real indexes."""
    document = fixture.workflow_inputs.to_dict()
    samples = document["samples"]
    if not isinstance(samples, list):
        raise AssertionError("rapid_quant fixture samples must be inline")
    sample_ids = sorted(
        {
            row["sample"]
            for row in samples
            if row.get("layout") == "SE" and isinstance(row.get("sample"), str)
        }
    )
    if not sample_ids:
        raise AssertionError("rapid_quant fixture requires one controlled SE sample")
    selected_id = sample_ids[0]
    selected = [deepcopy(row) for row in samples if row.get("sample") == selected_id]
    if len(selected) != 1:
        raise AssertionError("rapid_quant route requires one single-lane SE sample")
    standard = deepcopy(document["config"]["standard"])
    trimming = standard.get("trimming")
    if not isinstance(trimming, dict) or not isinstance(trimming.get("tool"), str):
        raise AssertionError("rapid_quant fixture requires a typed trimming policy")
    standard["trimming"] = {"enabled": False, "tool": trimming["tool"]}
    standard["ribosomal_rna_removal"] = {"enabled": False}
    return WorkflowInputs(
        config={"standard": standard},
        samples=selected,
        options=deepcopy(document["options"]),
    )


def _broken_fastq_inputs(
    inputs: WorkflowInputs,
    *,
    destination: Path,
) -> WorkflowInputs:
    """Replace only the controlled read with deterministic corrupt gzip."""
    if destination.exists() or not destination.is_absolute():
        raise AssertionError("broken FASTQ destination must be a new absolute path")
    compressed = bytearray(gzip.compress(b"@broken\nACGT\n+\nIIII\n", mtime=0))
    compressed[-8] ^= 0xFF
    destination.write_bytes(compressed)
    document = inputs.to_dict()
    samples = document["samples"]
    assert isinstance(samples, list) and len(samples) == 1
    samples[0]["fastq_1"] = str(destination)
    return WorkflowInputs(
        config=document["config"],
        samples=samples,
        options=document["options"],
    )


def _execute_rapid_quant(
    *,
    adapter: BulkRnaSeqRapidQuantQualificationAdapter,
    inputs: WorkflowInputs,
    settings: GateSettings,
    workspace: Path,
    expect_success: bool = True,
) -> RapidQuantRun:
    validation = adapter.validate(inputs)
    assert validation.is_success
    planned = adapter.plan_workspace(inputs, workspace)
    assert planned.is_success
    materialized = WorkspaceMaterializer().materialize(planned.value, workspace)
    assert materialized.is_success
    command_result = adapter.build_command(planned.value, workspace)
    assert command_result.is_success
    command = command_result.value
    _assert_rapid_command_contract(command, workspace)

    cleaner = ManagedContainerCleaner(
        executable=settings.docker_executable,
        unix_socket=settings.docker_socket,
    )
    runner = ProcessRunner(
        allowed_executables=(command.argv[0],),
        timeout_seconds=_EXECUTION_TIMEOUT_SECONDS,
        managed_container_cleaner=cleaner,
    )
    preflight = runner.run(_preflight_spec(command))
    assert preflight.is_success
    assert preflight.value.exit_code == 0
    assert not preflight.issues

    executed = runner.run(command)
    assert executed.is_success
    if expect_success:
        assert executed.value.exit_code == 0
    else:
        assert executed.value.exit_code != 0
    assert command.managed_container_scope is not None
    try:
        assert_no_managed_containers(cleaner, command.managed_container_scope)
    finally:
        if cleaner.cleanup(command.managed_container_scope).is_failure:
            raise AssertionError("rapid_quant hygiene cleanup failed")
    assert_no_managed_containers(cleaner, command.managed_container_scope)
    _assert_workspace_processes_reaped(workspace)
    _assert_nonempty_nextflow_log(workspace)

    scripts = _container_scripts(workspace, command.managed_container_scope)
    quant_surface = _quant_surface(workspace)
    trace_processes = _trace_processes(workspace) if expect_success else ()
    return RapidQuantRun(
        process_result=executed.value,
        process_issues=executed.issues,
        quant_surface=quant_surface,
        trace_processes=trace_processes,
        container_script_count=len(scripts),
    )


def _preflight_spec(command: CommandSpec) -> CommandSpec:
    assert command.preflight_argv is not None
    return CommandSpec(
        argv=command.preflight_argv,
        cwd=command.cwd,
        env=command.env,
        redaction_values=command.redaction_values,
    )


def _assert_rapid_command_contract(command: CommandSpec, workspace: Path) -> None:
    assert command.cwd == str(workspace / "engine/launch")
    assert command.argv.count("-C") == 1
    assert command.argv.count("-offline") == 1
    assert "-c" not in command.argv
    assert "-resume" not in command.argv
    assert command.argv[command.argv.index("-profile") + 1] == _RAPID_PROFILE
    assert command.preflight_argv is not None
    assert command.preflight_argv.count("-C") == 1
    assert (
        command.preflight_argv[command.preflight_argv.index("-profile") + 1]
        == _RAPID_PROFILE
    )
    assert command.env["NXF_OFFLINE"] == "true"
    assert command.env["NXF_DISABLE_CHECK_LATEST"] == "true"
    assert command.managed_container_scope is not None
    assert command.managed_container_endpoint_identity is not None

    config = _read_bounded(workspace / "config/platform.nextflow.config")
    assert "docker.registry = ''" in config
    assert "docker.remove = true" in config
    assert "wave.enabled = false" in config
    assert "--pull=never" in config
    assert "--network=none" in config
    assert (
        f"--label={MANAGED_CONTAINER_SCOPE_LABEL}={command.managed_container_scope}"
    ) in config

    params = json.loads(_read_bounded(workspace / "config/params.json"))
    assert set(params).isdisjoint(RAPID_QUANT_PROFILE_OWNED_PARAMETERS)
    assert set(params).isdisjoint(RAPID_QUANT_MODE_OWNED_PARAMETERS)
    identity = json.loads(_read_bounded(workspace / "config/execution-identity.json"))
    assert identity["execution_mode"] == _RAPID_MODE
    assert identity["resume_enabled"] is False


def _container_scripts(workspace: Path, scope: str) -> tuple[Path, ...]:
    scripts = tuple(sorted((workspace / "engine/work").rglob(".command.run")))
    if not scripts:
        raise AssertionError("rapid_quant did not materialize a container command")
    expected_label = f"--label={MANAGED_CONTAINER_SCOPE_LABEL}={scope}"
    for script in scripts:
        content = _read_bounded(script)
        assert "docker run" in content
        assert "--pull=never" in content
        assert "--network=none" in content
        assert expected_label in content
        assert _SHA256_IMAGE.search(content) is not None
        assert "quay.io/" not in content
        assert "docker.io/" not in content
        assert "community.wave.seqera.io/" not in content
    return scripts


def _quant_surface(workspace: Path) -> tuple[tuple[str, str], ...]:
    results_root = workspace / "results"
    quant_files = tuple(sorted(results_root.glob("*/salmon/*/quant.sf")))
    surface = []
    for path in quant_files:
        relative = path.relative_to(results_root).as_posix()
        surface.append((relative, _regular_file_sha256(path)))
    return tuple(surface)


def _trace_processes(workspace: Path) -> tuple[str, ...]:
    lines = _read_bounded(workspace / "reports/trace.txt").splitlines()
    if len(lines) < 2:
        raise AssertionError("rapid_quant trace is empty")
    header = lines[0].split("\t")
    try:
        name_index = header.index("name")
        status_index = header.index("status")
    except ValueError:
        raise AssertionError("rapid_quant trace fields are incomplete") from None
    rows = [line.split("\t") for line in lines[1:] if line]
    if any(len(row) != len(header) for row in rows):
        raise AssertionError("rapid_quant trace rows are malformed")
    if any(row[status_index] != "COMPLETED" for row in rows):
        raise AssertionError("rapid_quant successful trace is not complete")
    return tuple(sorted(row[name_index].split(" (")[0] for row in rows))


def _assert_nonempty_nextflow_log(workspace: Path) -> None:
    log = workspace / "logs/nextflow.log"
    if not _read_bounded(log).strip():
        raise AssertionError("rapid_quant Nextflow log is empty")


def _assert_workspace_processes_reaped(workspace: Path) -> None:
    marker = os.fsencode(str(workspace))
    deadline = time.monotonic() + 5
    while time.monotonic() < deadline:
        if not _workspace_processes(marker):
            return
        time.sleep(0.05)
    _terminate_workspace_process_groups(_workspace_processes(marker))
    raise AssertionError("rapid_quant left a workspace process alive")


def _workspace_processes(marker: bytes) -> tuple[int, ...]:
    observed = []
    try:
        candidates = tuple(Path("/proc").iterdir())
    except OSError:
        raise AssertionError("rapid_quant process cleanup cannot be audited") from None
    for candidate in candidates:
        if not candidate.name.isdigit():
            continue
        try:
            command = (candidate / "cmdline").read_bytes()
        except OSError:
            continue
        if marker in command:
            observed.append(int(candidate.name))
    return tuple(sorted(observed))


def _terminate_workspace_process_groups(process_ids: tuple[int, ...]) -> None:
    own_group = os.getpgrp()
    groups: set[int] = set()
    for process_id in process_ids:
        try:
            group = os.getpgid(process_id)
        except (OSError, ProcessLookupError):
            continue
        if group > 0 and group != own_group:
            groups.add(group)
    for requested_signal in (signal.SIGTERM, signal.SIGKILL):
        for group in groups:
            try:
                os.killpg(group, requested_signal)
            except (OSError, ProcessLookupError):
                continue
        if requested_signal is signal.SIGTERM and groups:
            time.sleep(0.25)


def _read_bounded(path: Path) -> str:
    try:
        status = path.lstat()
        if (
            not stat.S_ISREG(status.st_mode)
            or status.st_size <= 0
            or status.st_size > _MAX_GATE_FILE_BYTES
        ):
            raise OSError
        return path.read_text(encoding="utf-8")
    except (OSError, UnicodeError):
        raise AssertionError("rapid_quant gate file is unavailable") from None


def _regular_file_sha256(path: Path) -> str:
    try:
        status = path.lstat()
        if (
            not stat.S_ISREG(status.st_mode)
            or status.st_size <= 0
            or status.st_size > _MAX_GATE_FILE_BYTES
        ):
            raise OSError
        digest = sha256()
        with path.open("rb") as stream:
            while chunk := stream.read(1024 * 1024):
                digest.update(chunk)
        return digest.hexdigest()
    except OSError:
        raise AssertionError("rapid_quant output file is unavailable") from None


def _unit_command(workspace: Path) -> CommandSpec:
    config = workspace / "config/platform.nextflow.config"
    config.parent.mkdir(parents=True)
    scope = "a" * 64
    config.write_text(
        "\n".join(
            (
                "docker.registry = ''",
                "docker.remove = true",
                "wave.enabled = false",
                (
                    "docker.runOptions = '--pull=never --network=none "
                    f"--label={MANAGED_CONTAINER_SCOPE_LABEL}={scope}'"
                ),
            )
        ),
        encoding="utf-8",
    )
    (workspace / "config/params.json").write_text("{}\n", encoding="utf-8")
    (workspace / "config/execution-identity.json").write_text(
        json.dumps({"execution_mode": _RAPID_MODE, "resume_enabled": False}) + "\n",
        encoding="utf-8",
    )
    argv = (
        "/runtime/nextflow",
        "-C",
        str(config),
        "run",
        "/runtime/source",
        "-profile",
        _RAPID_PROFILE,
        "-offline",
    )
    return CommandSpec(
        argv=argv,
        cwd=str(workspace / "engine/launch"),
        env={"NXF_OFFLINE": "true", "NXF_DISABLE_CHECK_LATEST": "true"},
        preflight_argv=(
            "/runtime/nextflow",
            "-C",
            str(config),
            "config",
            "/runtime/source",
            "-profile",
            _RAPID_PROFILE,
        ),
        managed_container_scope=scope,
        managed_container_endpoint_identity="b" * 64,
    )
