"""Focused composition contracts for the private real-execution runner."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from encode_pipeline.platform.adapters import WorkflowInputs
from encode_pipeline.services.managed_containers import ManagedContainerCleaner
from encode_pipeline.workers.timeouts import WorkerHardTimeout

from .support import (
    AcceptanceFixture,
    GateSettings,
    build_acceptance_process_runner,
)
from .test_platform_acceptance import _assert_full_trace_contract


def _docker_coordinates(tmp_path: Path) -> tuple[Path, Path]:
    executable = tmp_path / "bin" / "docker"
    executable.parent.mkdir()
    executable.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
    executable.chmod(0o755)
    socket_path = tmp_path / "docker.sock"
    socket_path.write_text("test-only socket coordinate\n", encoding="utf-8")
    return executable, socket_path


def _write_full_trace(path: Path, processes: tuple[str, ...]) -> None:
    path.write_text(
        "name\tstatus\n"
        + "".join(f"NFCORE_RNASEQ:{name}\tCOMPLETED\n" for name in processes),
        encoding="utf-8",
    )


def test_full_trace_accepts_tagged_required_processes(tmp_path: Path) -> None:
    trace = tmp_path / "trace.txt"
    _write_full_trace(
        trace,
        (
            "STAR_ALIGN (PE1)",
            "SALMON_QUANT (PE1)",
            "FASTQC (PE1)",
            "SORTMERNA (PE1)",
        ),
    )

    _assert_full_trace_contract(trace)


@pytest.mark.parametrize(
    "process",
    ("STAR_GENOMEGENERATE", "SALMON_INDEX", "SORTMERNA_INDEX"),
)
def test_full_trace_rejects_runtime_index_builds(
    tmp_path: Path,
    process: str,
) -> None:
    trace = tmp_path / "trace.txt"
    _write_full_trace(
        trace,
        (
            "STAR_ALIGN (SE1)",
            "SALMON_QUANT (SE1)",
            "FASTQC (SE1)",
            "SORTMERNA (SE1)",
            f"{process} (tiny)",
        ),
    )

    with pytest.raises(AssertionError):
        _assert_full_trace_contract(trace)


def test_acceptance_process_runner_binds_unshare_docker_and_timeout(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    docker_executable, docker_socket = _docker_coordinates(tmp_path)
    monkeypatch.setattr(
        ManagedContainerCleaner,
        "_endpoint_identities",
        lambda _self: ((1, 2, 3), (4, 5, 6)),
    )
    settings = GateSettings(
        runtime_root=(tmp_path / "runtime").resolve(),
        fixture_manifest=(tmp_path / "fixture.json").resolve(),
        redis_url="redis://127.0.0.1:6379/15",
        docker_executable=docker_executable,
        docker_socket=docker_socket,
    )
    binding = SimpleNamespace(
        assets=SimpleNamespace(
            network_isolation_executable=Path("/usr/bin/unshare"),
            docker_executable=docker_executable,
            docker_socket=docker_socket,
        )
    )

    runner = build_acceptance_process_runner(
        settings=settings,
        binding=binding,
        timeout_seconds=37,
        passthrough_exceptions=(WorkerHardTimeout,),
    )

    assert runner._allowed_executables == ("/usr/bin/unshare",)
    assert runner._timeout_seconds == 37.0
    assert runner._passthrough_exceptions == (WorkerHardTimeout,)
    cleaner = runner._managed_container_cleaner
    assert cleaner is not None
    assert cleaner.executable == docker_executable
    assert cleaner.unix_socket == docker_socket


def test_acceptance_process_runner_rejects_docker_binding_drift(
    tmp_path: Path,
) -> None:
    docker_executable, docker_socket = _docker_coordinates(tmp_path)
    settings = GateSettings(
        runtime_root=(tmp_path / "runtime").resolve(),
        fixture_manifest=(tmp_path / "fixture.json").resolve(),
        redis_url="redis://127.0.0.1:6379/15",
        docker_executable=docker_executable,
        docker_socket=docker_socket,
    )
    binding = SimpleNamespace(
        assets=SimpleNamespace(
            network_isolation_executable=Path("/usr/bin/unshare"),
            docker_executable=tmp_path / "other" / "docker",
            docker_socket=docker_socket,
        )
    )

    try:
        build_acceptance_process_runner(
            settings=settings,
            binding=binding,
            timeout_seconds=37,
        )
    except ValueError as error:
        assert (
            str(error) == "execution binding Docker endpoint differs from gate settings"
        )
    else:
        raise AssertionError("Docker binding drift was accepted")


def test_platform_submission_injects_the_acceptance_process_runner(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from . import platform_harness

    docker_executable, docker_socket = _docker_coordinates(tmp_path)
    settings = GateSettings(
        runtime_root=(tmp_path / "runtime").resolve(),
        fixture_manifest=(tmp_path / "fixture.json").resolve(),
        redis_url="redis://127.0.0.1:6379/15",
        docker_executable=docker_executable,
        docker_socket=docker_socket,
    )
    binding = SimpleNamespace(assets=SimpleNamespace())
    composition = SimpleNamespace(
        binding=binding,
        registry=object(),
        build_identity_provider=object(),
    )
    monkeypatch.setattr(
        platform_harness,
        "build_results_composition",
        lambda *_args, **_kwargs: composition,
    )
    harness = platform_harness.PlatformAcceptanceHarness(
        gate_settings=settings,
        repository_root=tmp_path.resolve(),
        temporary_root=(tmp_path / "acceptance").resolve(),
        job_timeout_seconds=41,
    )
    harness._run_queue = object()

    runner = object()
    runner_arguments: dict[str, object] = {}

    def build_runner(**kwargs):
        runner_arguments.update(kwargs)
        return runner

    captured_runtime_arguments: dict[str, object] = {}
    runtime = SimpleNamespace(
        persistence=SimpleNamespace(repository=object()),
        preflight_service=SimpleNamespace(
            preflight=lambda _run_id: SimpleNamespace(
                is_failure=False,
                issues=(),
            )
        ),
        run_service=SimpleNamespace(
            get_execution_assignment=lambda _run_id: SimpleNamespace(job_id="job-1")
        ),
        build_identity_provider=composition.build_identity_provider,
    )

    class RuntimeContext:
        def __enter__(self):
            return runtime

        def __exit__(self, *_args):
            return None

    def open_runtime(*args, **kwargs):
        captured_runtime_arguments["args"] = args
        captured_runtime_arguments["kwargs"] = kwargs
        return RuntimeContext()

    class FakeValidatedInputService:
        def __init__(self, **_kwargs):
            pass

        def validate(self, _workflow_id, _workflow_inputs):
            return SimpleNamespace(
                is_failure=False,
                value=SimpleNamespace(
                    snapshot_id="snapshot-1",
                    payload_digest="a" * 64,
                ),
                issues=(),
            )

    class FakeValidatedRunCreationService:
        def __init__(self, **_kwargs):
            pass

        def create_run(self, _workflow_id, _snapshot_id):
            return SimpleNamespace(record=SimpleNamespace(run_id="run-1"))

    class FakeRunSubmissionService:
        def __init__(self, **_kwargs):
            pass

        def start_run(self, _run_id):
            return SimpleNamespace(status=SimpleNamespace(value="queued"))

    monkeypatch.setattr(
        platform_harness, "build_acceptance_process_runner", build_runner
    )
    monkeypatch.setattr(platform_harness, "open_worker_runtime", open_runtime)
    monkeypatch.setattr(
        platform_harness, "ValidationService", lambda **_kwargs: object()
    )
    monkeypatch.setattr(
        platform_harness,
        "ValidatedInputService",
        FakeValidatedInputService,
    )
    monkeypatch.setattr(
        platform_harness,
        "ValidatedRunCreationService",
        FakeValidatedRunCreationService,
    )
    monkeypatch.setattr(
        platform_harness,
        "RunSubmissionService",
        FakeRunSubmissionService,
    )
    fixture = AcceptanceFixture(
        workflow_inputs=WorkflowInputs(config={}, samples=None, options={}),
        transcriptome=object(),
        acceptance_manifest_sha256="a" * 64,
        source_manifest_sha256="b" * 64,
        source_identity_sha256="c" * 64,
        index_provenance_manifest_sha256="d" * 64,
        index_provenance_identity_sha256="e" * 64,
        required_artifact_output_types=(),
        required_qc_metric_keys=(),
        required_sample_ids=(),
        required_artifact_sample_output_types=(),
        required_qc_sample_metric_keys=(),
        required_qc_sample_metric_values=(),
    )
    monkeypatch.setattr(
        platform_harness, "load_acceptance_fixture", lambda _path: fixture
    )

    submitted = harness._submit(fixture)

    assert submitted.run_id == "run-1"
    assert runner_arguments == {
        "settings": settings,
        "binding": binding,
        "timeout_seconds": 41,
        "passthrough_exceptions": (WorkerHardTimeout,),
    }
    assert captured_runtime_arguments["args"] == (harness.worker_settings,)
    assert captured_runtime_arguments["kwargs"] == {
        "registry": composition.registry,
        "build_identity_provider": composition.build_identity_provider,
        "process_runner": runner,
    }


def test_worker_entry_injects_a_fresh_hard_timeout_runner(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from encode_pipeline.workers import jobs as worker_jobs

    from . import worker_entry

    docker_executable, docker_socket = _docker_coordinates(tmp_path)
    monkeypatch.setattr(
        ManagedContainerCleaner,
        "_endpoint_identities",
        lambda _self: ((1, 2, 3), (4, 5, 6)),
    )
    gate_settings = GateSettings(
        runtime_root=(tmp_path / "runtime").resolve(),
        fixture_manifest=(tmp_path / "fixture.json").resolve(),
        redis_url="redis://127.0.0.1:6379/15",
        docker_executable=docker_executable,
        docker_socket=docker_socket,
    )
    binding = SimpleNamespace(
        assets=SimpleNamespace(
            network_isolation_executable=Path("/usr/bin/unshare"),
            docker_executable=docker_executable,
            docker_socket=docker_socket,
        )
    )
    composition = SimpleNamespace(
        binding=binding,
        registry=object(),
        build_identity_provider=object(),
    )
    worker_settings = SimpleNamespace(job_timeout_seconds=53)
    captured: dict[str, object] = {}
    restored_runtime_factory = worker_jobs.open_worker_runtime

    def open_runtime(*args, **kwargs):
        captured["args"] = args
        captured["kwargs"] = kwargs
        return object()

    def run_worker(arguments):
        captured["worker_arguments"] = arguments
        captured["runtime"] = worker_jobs.open_worker_runtime()
        return 0

    monkeypatch.setattr(worker_entry, "require_gate_settings", lambda: gate_settings)
    monkeypatch.setattr(
        worker_entry,
        "build_results_composition",
        lambda *_args, **_kwargs: composition,
    )
    monkeypatch.setattr(worker_entry, "load_worker_settings", lambda: worker_settings)
    monkeypatch.setattr(worker_entry, "open_worker_runtime", open_runtime)
    monkeypatch.setattr(worker_entry, "worker_main", run_worker)

    assert worker_entry.main() == 0

    assert captured["worker_arguments"] == ("--burst",)
    assert captured["args"] == (worker_settings,)
    runtime_kwargs = captured["kwargs"]
    assert runtime_kwargs["registry"] is composition.registry
    assert (
        runtime_kwargs["build_identity_provider"] is composition.build_identity_provider
    )
    runner = runtime_kwargs["process_runner"]
    assert runner._allowed_executables == ("/usr/bin/unshare",)
    assert runner._timeout_seconds == 53.0
    assert runner._passthrough_exceptions == (WorkerHardTimeout,)
    assert runner._managed_container_cleaner.executable == docker_executable
    assert runner._managed_container_cleaner.unix_socket == docker_socket
    assert worker_jobs.open_worker_runtime is restored_runtime_factory
