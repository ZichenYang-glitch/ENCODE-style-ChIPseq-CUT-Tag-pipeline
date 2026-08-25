"""Focused composition contracts for the private real-execution runner."""

from __future__ import annotations

import json
from pathlib import Path
import stat
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


def _acceptance_fixture(tmp_path: Path) -> AcceptanceFixture:
    reference = {
        "reference_id": "tiny",
        "fasta": str(tmp_path / "reference.fa"),
        "fasta_sha256": "1" * 64,
        "gtf": str(tmp_path / "genes.gtf"),
        "gtf_sha256": "2" * 64,
        "annotation_style": "ensembl",
        "star_index": {
            "path": str(tmp_path / "star"),
            "identity_sha256": "3" * 64,
        },
        "salmon_index": {
            "path": str(tmp_path / "salmon"),
            "identity_sha256": "4" * 64,
        },
    }
    return AcceptanceFixture(
        workflow_inputs=WorkflowInputs(
            config={"standard": {"reference": reference}},
            samples=[{"sample": "tiny"}],
            options={},
        ),
        transcriptome=SimpleNamespace(
            reference_id="tiny",
            fasta_sha256="1" * 64,
            gtf_sha256="2" * 64,
            transcript_fasta=tmp_path / "transcripts.fa",
            transcript_fasta_sha256="5" * 64,
        ),
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


def _private_config_harness(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    from . import platform_harness

    docker_executable, docker_socket = _docker_coordinates(tmp_path)
    monkeypatch.setattr(
        platform_harness,
        "build_results_composition",
        lambda *_args, **_kwargs: SimpleNamespace(
            binding=SimpleNamespace(assets=SimpleNamespace()),
            registry=SimpleNamespace(get=lambda _workflow_id: object()),
            build_identity_provider=object(),
        ),
    )
    harness = platform_harness.PlatformAcceptanceHarness(
        gate_settings=GateSettings(
            runtime_root=(tmp_path / "runtime").resolve(),
            fixture_manifest=(tmp_path / "fixture.json").resolve(),
            redis_url="redis://127.0.0.1:6379/15",
            docker_executable=docker_executable,
            docker_socket=docker_socket,
        ),
        repository_root=tmp_path.resolve(),
        temporary_root=(tmp_path / "acceptance").resolve(),
        job_timeout_seconds=41,
    )
    # Mirror the harness execution-stage entry: prepare the fresh database
    # once before any existing-only persistence access.
    harness.temporary_root.mkdir(parents=True, exist_ok=True)
    platform_harness.prepare_acceptance_database(harness.database_url)
    return harness


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


def test_execution_activity_refreshes_claimed_job_before_observing_worker(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from encode_pipeline.platform.runs import RunStatus
    from encode_pipeline.services import runs as runs_module
    from . import platform_harness

    harness = _private_config_harness(tmp_path, monkeypatch)
    submitted = platform_harness.SubmittedAcceptanceRun(
        run_id="run-observe-claimed-worker",
        job_id="run-execution-observe-claimed-worker",
        validated_snapshot_id="vsnap_observe_claimed_worker",
        fixture_acceptance_manifest_sha256="a" * 64,
    )
    full_refreshes: list[None] = []
    status_refreshes: list[bool] = []
    worker_reads: list[int] = []

    class FakeJob:
        def __init__(self) -> None:
            self._status = platform_harness.JobStatus.QUEUED
            self._worker_name: str | None = None

        @property
        def worker_name(self) -> str | None:
            worker_reads.append(len(full_refreshes))
            return self._worker_name

        def refresh(self) -> None:
            full_refreshes.append(None)
            self._status = platform_harness.JobStatus.STARTED
            self._worker_name = "worker-observe-claimed-worker"

        def get_status(self, *, refresh: bool = False):
            status_refreshes.append(refresh)
            if refresh:
                self._status = platform_harness.JobStatus.STARTED
            return self._status

    class FakeProcess:
        pid = 987_654_321

        @staticmethod
        def poll() -> None:
            return None

    record = SimpleNamespace(status=RunStatus.RUNNING)
    assignment = SimpleNamespace(
        dispatched_at=object(),
        claimed_at=object(),
    )

    class FakeRunService:
        def __init__(self, *_args, **_kwargs) -> None:
            return None

        @staticmethod
        def get_run(run_id):
            assert run_id == submitted.run_id
            return record

        @staticmethod
        def get_execution_assignment(run_id):
            assert run_id == submitted.run_id
            return assignment

    job = FakeJob()
    process = FakeProcess()
    harness._submitted.append(submitted)
    harness._worker_processes.append(process)
    harness._run_queue = SimpleNamespace(
        _queue=SimpleNamespace(
            fetch_job=lambda job_id: job if job_id == submitted.job_id else None
        )
    )
    monkeypatch.setattr(runs_module, "RunService", FakeRunService)
    monkeypatch.setattr(
        platform_harness,
        "open_run_persistence",
        lambda _database_url: SimpleNamespace(
            repository=object(),
            close=lambda: None,
        ),
    )
    monkeypatch.setattr(
        platform_harness.ManagedContainerCleaner,
        "_endpoint_identities",
        lambda _self: ((1, 2, 3), (4, 5, 6)),
    )
    monkeypatch.setattr(
        platform_harness.ManagedContainerCleaner,
        "verify_endpoint",
        lambda _self: SimpleNamespace(is_failure=False),
    )
    monkeypatch.setattr(
        platform_harness,
        "_worker_session_nextflow_processes",
        lambda *_args, **_kwargs: (123,),
    )
    monkeypatch.setattr(
        platform_harness,
        "_worker_session_process_groups",
        lambda *_args, **_kwargs: (456,),
    )

    def reject_another_poll(_seconds: float) -> None:
        raise AssertionError("activity observer did not discover the claimed worker")

    monkeypatch.setattr(platform_harness.time, "sleep", reject_another_poll)

    evidence = harness.wait_for_execution_activity(
        submitted,
        process,
        require_managed_container=False,
        timeout_seconds=1,
    )

    assert evidence == platform_harness.ExecutionActivityEvidence(
        worker_session_id=process.pid,
        process_group_count=1,
        nextflow_observed=True,
        managed_container_observed=False,
    )
    assert full_refreshes == [None]
    assert status_refreshes == [False]
    assert worker_reads == [1]


@pytest.mark.parametrize(
    "cleanup_observation",
    ["clean", "residual", "constructor-error"],
)
@pytest.mark.parametrize(
    (
        "terminal_status_name",
        "rq_observation",
        "expected_rq_status",
        "error_code",
        "error_reason_code",
    ),
    [
        (
            "FAILED",
            "failed",
            "failed",
            "RUN_EXECUTION_FAILED",
            "WORKER_JOB_TIMEOUT",
        ),
        ("SUCCEEDED", "finished", "finished", None, None),
        (
            "FAILED",
            "missing",
            "unavailable",
            "RUN_EXECUTION_FAILED",
            "WORKER_JOB_TIMEOUT",
        ),
        (
            "FAILED",
            "status-error",
            "unavailable",
            "RUN_EXECUTION_FAILED",
            "WORKER_JOB_TIMEOUT",
        ),
    ],
)
def test_terminal_before_activity_publishes_path_free_reason_evidence(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    cleanup_observation: str,
    terminal_status_name: str,
    rq_observation: str,
    expected_rq_status: str,
    error_code: str | None,
    error_reason_code: str | None,
) -> None:
    from encode_pipeline.platform.runs import RunStatus
    from encode_pipeline.services import runs as runs_module
    from . import platform_harness

    harness = _private_config_harness(tmp_path, monkeypatch)
    submitted = platform_harness.SubmittedAcceptanceRun(
        run_id="run-terminal-before-activity",
        job_id="run-execution-terminal-before-activity",
        validated_snapshot_id="vsnap_terminal_before_activity",
        fixture_acceptance_manifest_sha256="a" * 64,
    )
    private_detail = str(tmp_path / "private-terminal-detail")
    cleanup_confirmed = cleanup_observation == "clean"

    terminal_status = getattr(RunStatus, terminal_status_name)
    rq_status = {
        "failed": platform_harness.JobStatus.FAILED,
        "finished": platform_harness.JobStatus.FINISHED,
    }.get(rq_observation)
    status_reads: list[None] = []
    evidence_path = (
        harness.temporary_root
        / "evidence"
        / f"terminal-lifecycle-{submitted.run_id}.json"
    )

    class FakeJob:
        worker_name = "worker-terminal-before-activity"

        @staticmethod
        def get_status(*, refresh: bool = False):
            assert refresh is True
            status_reads.append(None)
            if rq_observation == "status-error":
                raise RuntimeError(private_detail)
            assert rq_status is not None
            return rq_status

    class FakeProcess:
        pid = 987_654_321
        returncode: int | None = None

        def poll(self):
            return self.returncode

        def communicate(self, *, timeout):
            assert timeout > 0
            assert evidence_path.is_file()
            self.returncode = 0
            return ("", "")

    record = SimpleNamespace(
        status=terminal_status,
        cancellation_reason=private_detail,
        error=(
            None
            if error_code is None
            else SimpleNamespace(
                code=error_code,
                context={
                    "reason_code": error_reason_code,
                    "private_detail": private_detail,
                },
            )
        ),
    )
    assignment = SimpleNamespace(
        job_id=submitted.job_id,
        dispatched_at=object(),
        claimed_at=object(),
        cancellation_requested_at=None,
        cancellation_acknowledged_at=None,
    )
    events = (
        SimpleNamespace(event_type="status_changed", status=RunStatus.QUEUED),
        SimpleNamespace(event_type="status_changed", status=RunStatus.RUNNING),
        SimpleNamespace(event_type="status_changed", status=terminal_status),
    )
    result_state = SimpleNamespace(
        artifact_revision=0,
        artifact_attempt_id=None,
        artifact_attempt_status=None,
        qc_revision=0,
        qc_attempt_id=None,
        qc_attempt_status=None,
    )

    class FakeRunService:
        def __init__(self, *_args, **_kwargs) -> None:
            return None

        @staticmethod
        def get_run(run_id):
            assert run_id == submitted.run_id
            return record

        @staticmethod
        def get_execution_assignment(run_id):
            assert run_id == submitted.run_id
            return assignment

        @staticmethod
        def list_events(run_id, *, limit):
            assert run_id == submitted.run_id
            assert limit == 1000
            return events

        @staticmethod
        def get_result_state(run_id):
            assert run_id == submitted.run_id
            return result_state

        @staticmethod
        def list_artifacts(run_id):
            assert run_id == submitted.run_id
            return ()

        @staticmethod
        def list_qc_metrics(run_id):
            assert run_id == submitted.run_id
            return ()

    persistence = SimpleNamespace(repository=object(), close=lambda: None)
    job = FakeJob()
    process = FakeProcess()
    harness._submitted.append(submitted)
    harness._worker_processes.append(process)

    def fetch_job(job_id):
        assert job_id == submitted.job_id
        if rq_observation == "missing":
            return None
        return job

    harness._run_queue = SimpleNamespace(_queue=SimpleNamespace(fetch_job=fetch_job))
    monkeypatch.setattr(runs_module, "RunService", FakeRunService)
    monkeypatch.setattr(
        platform_harness,
        "open_run_persistence",
        lambda _database_url: persistence,
    )
    monkeypatch.setattr(
        platform_harness.ManagedContainerCleaner,
        "_endpoint_identities",
        lambda _self: ((1, 2, 3), (4, 5, 6)),
    )
    monkeypatch.setattr(
        platform_harness.ManagedContainerCleaner,
        "verify_endpoint",
        lambda _self: SimpleNamespace(is_failure=False),
    )

    if cleanup_observation == "constructor-error":

        def fail_cleaner_construction(**_kwargs):
            raise RuntimeError(private_detail)

        monkeypatch.setattr(
            platform_harness,
            "ManagedContainerCleaner",
            fail_cleaner_construction,
        )

    def assert_cleanup(*_args, **_kwargs) -> None:
        if cleanup_observation == "residual":
            raise AssertionError(private_detail)

    monkeypatch.setattr(
        platform_harness,
        "assert_no_managed_containers",
        assert_cleanup,
    )

    expected_error_code = error_code or "UNAVAILABLE"
    expected_error_reason_code = error_reason_code or "UNAVAILABLE"
    with pytest.raises(AssertionError) as raised:
        harness.wait_for_execution_activity(
            submitted,
            process,
            require_managed_container=False,
            timeout_seconds=1,
        )
    assert str(raised.value) == (
        "TERMINAL_BEFORE_REQUIRED_ACTIVITY "
        f"error_code={expected_error_code} "
        f"error_reason_code={expected_error_reason_code}"
    )

    rendered = evidence_path.read_text(encoding="utf-8")
    payload = json.loads(rendered)
    assert payload["assertion_reason_code"] == "TERMINAL_BEFORE_REQUIRED_ACTIVITY"
    assert payload["lifecycle_status"] == terminal_status.value
    assert payload["lifecycle_history"] == [
        "queued",
        "running",
        terminal_status.value,
    ]
    assert payload["error_code"] == error_code
    assert payload["error_reason_code"] == error_reason_code
    assert payload["assignment_dispatched"] is True
    assert payload["assignment_claimed"] is True
    assert payload["rq_status"] == expected_rq_status
    assert payload["rq_failed"] is (
        expected_rq_status == platform_harness.JobStatus.FAILED.value
    )
    assert payload["rq_stopped"] is False
    assert payload["rq_finished"] is (
        expected_rq_status == platform_harness.JobStatus.FINISHED.value
    )
    assert len(status_reads) == (0 if rq_observation == "missing" else 1)
    assert payload["cleanup_confirmed"] is cleanup_confirmed
    assert "/" not in rendered
    assert str(tmp_path) not in rendered
    assert private_detail not in rendered
    assert process not in harness._worker_processes


def test_rq_terminal_metadata_stabilization_is_bounded_without_sleeping() -> None:
    from . import platform_harness

    clock = SimpleNamespace(value=0.0)
    calls: list[float] = []

    class FakeJob:
        @staticmethod
        def get_status(*, refresh: bool):
            assert refresh is True
            calls.append(clock.value)
            return platform_harness.JobStatus.STARTED

    status = platform_harness._wait_for_rq_terminal_status(
        FakeJob(),
        timeout_seconds=0.2,
        monotonic=lambda: clock.value,
        sleep=lambda seconds: setattr(clock, "value", clock.value + seconds),
    )

    assert status is platform_harness.JobStatus.STARTED
    assert clock.value == pytest.approx(0.2)
    assert len(calls) == 5


def test_rq_terminal_metadata_stabilization_observes_failed_state() -> None:
    from . import platform_harness

    statuses = iter(
        (
            platform_harness.JobStatus.STARTED,
            platform_harness.JobStatus.FAILED,
        )
    )
    clock = SimpleNamespace(value=0.0)

    class FakeJob:
        @staticmethod
        def get_status(*, refresh: bool):
            assert refresh is True
            return next(statuses)

    status = platform_harness._wait_for_rq_terminal_status(
        FakeJob(),
        timeout_seconds=1,
        monotonic=lambda: clock.value,
        sleep=lambda seconds: setattr(clock, "value", clock.value + seconds),
    )

    assert status is platform_harness.JobStatus.FAILED
    assert clock.value == pytest.approx(platform_harness._RQ_TERMINAL_POLL_SECONDS)


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
        registry=SimpleNamespace(get=lambda _workflow_id: object()),
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
    harness.temporary_root.mkdir(parents=True, exist_ok=True)
    platform_harness.prepare_acceptance_database(harness.database_url)
    harness._run_queue = object()

    runner = object()
    runner_arguments: dict[str, object] = {}

    def build_runner(**kwargs):
        runner_arguments.update(kwargs)
        return runner

    captured_runtime_arguments: dict[str, object] = {}
    reference_binding = SimpleNamespace(
        revision_id="refpr_gate_tiny",
        revision_public_identity_sha256="f" * 64,
    )
    binding_service = object()
    reference_resolver = object()
    run_number = 0

    def create_run(_workflow_id, _snapshot_id):
        nonlocal run_number
        run_number += 1
        return SimpleNamespace(record=SimpleNamespace(run_id=f"run-{run_number}"))

    run_service = SimpleNamespace(
        get_execution_assignment=lambda run_id: SimpleNamespace(
            job_id=f"job-{run_id.removeprefix('run-')}"
        ),
        get_validated_reference_binding=lambda _snapshot_id: reference_binding,
        get_run_reference_binding=lambda _run_id: reference_binding,
    )
    runtime = SimpleNamespace(
        persistence=SimpleNamespace(
            repository=object(),
            reference_profile_repository=object(),
        ),
        preflight_service=SimpleNamespace(
            preflight=lambda _run_id: SimpleNamespace(
                is_failure=False,
                issues=(),
            )
        ),
        run_service=run_service,
        reference_profile_binding_service=binding_service,
        reference_profile_resolver=reference_resolver,
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

    service_arguments: dict[str, list[dict[str, object]]] = {
        "catalog": [],
        "validation": [],
        "creation": [],
        "submission": [],
    }
    registration_calls: list[dict[str, object]] = []
    enable_calls: list[tuple[str, str | None]] = []

    class FakeReferenceProfileService:
        def __init__(self, **kwargs):
            service_arguments["catalog"].append(kwargs)

        def register(self, **kwargs):
            registration_calls.append(kwargs)
            return SimpleNamespace(
                profile_id="refp_gate_tiny",
                revision_id=reference_binding.revision_id,
                public_identity_sha256=(
                    reference_binding.revision_public_identity_sha256
                ),
            )

        def enable(self, profile_id, *, revision_id=None):
            enable_calls.append((profile_id, revision_id))
            return SimpleNamespace(
                profile_id=profile_id,
                revision_id=revision_id,
                public_identity_sha256=(
                    reference_binding.revision_public_identity_sha256
                ),
                enabled=True,
            )

        def get_revision_summary(self, revision_id):
            assert revision_id == reference_binding.revision_id
            return SimpleNamespace(
                profile_id="refp_gate_tiny",
                revision_id=revision_id,
                public_identity_sha256=(
                    reference_binding.revision_public_identity_sha256
                ),
                enabled=True,
            )

    validated_revision_ids: list[str] = []
    validated_inputs: list[WorkflowInputs] = []

    class FakeValidatedInputService:
        def __init__(self, **kwargs):
            service_arguments["validation"].append(kwargs)

        def validate(
            self,
            _workflow_id,
            workflow_inputs,
            *,
            reference_profile_revision_id,
        ):
            validated_revision_ids.append(reference_profile_revision_id)
            validated_inputs.append(workflow_inputs)
            return SimpleNamespace(
                is_failure=False,
                value=SimpleNamespace(
                    snapshot_id=f"snapshot-{len(validated_revision_ids)}",
                    payload_digest="a" * 64,
                ),
                issues=(),
            )

    class FakeValidatedRunCreationService:
        def __init__(self, **kwargs):
            service_arguments["creation"].append(kwargs)

        def create_run(self, workflow_id, snapshot_id):
            return create_run(workflow_id, snapshot_id)

    class FakeRunSubmissionService:
        def __init__(self, **kwargs):
            service_arguments["submission"].append(kwargs)

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
        "ReferenceProfileService",
        FakeReferenceProfileService,
        raising=False,
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
    fixture = _acceptance_fixture(tmp_path)
    reference = fixture.workflow_inputs.config["standard"]["reference"]
    monkeypatch.setattr(
        platform_harness, "load_acceptance_fixture", lambda _path: fixture
    )

    first = harness._submit(fixture)
    second = harness._submit(fixture)

    assert (first.run_id, second.run_id) == ("run-1", "run-2")
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
    config_path = harness.worker_settings.reference_profile_config
    assert config_path is not None
    assert stat.S_IMODE(config_path.parent.stat().st_mode) == 0o700
    assert stat.S_IMODE(config_path.stat().st_mode) == 0o600
    private_config = json.loads(config_path.read_text(encoding="utf-8"))
    private_binding = private_config["profiles"]["bulk-rnaseq-gate-tiny-private"][
        "bindings"
    ]["bulk-rnaseq"]
    assert private_binding["reference"] == reference
    assert private_binding["transcriptome"] == {
        "reference_id": "tiny",
        "fasta_sha256": "1" * 64,
        "gtf_sha256": "2" * 64,
        "transcript_fasta": str(tmp_path / "transcripts.fa"),
        "transcript_fasta_sha256": "5" * 64,
    }
    assert registration_calls == [
        {
            "safe_key": "bulk-rnaseq-gate-tiny",
            "display_name": "Bulk RNA-seq protected tiny",
            "organism": "Synthetic organism",
            "assembly": "tiny",
            "config_key": "bulk-rnaseq-gate-tiny-private",
        }
    ]
    assert enable_calls == [("refp_gate_tiny", "refpr_gate_tiny")]
    assert validated_revision_ids == ["refpr_gate_tiny", "refpr_gate_tiny"]
    assert all(
        "reference" not in inputs.config["standard"] for inputs in validated_inputs
    )
    assert all(
        arguments["reference_profile_binding_service"] is binding_service
        and arguments["reference_profile_catalog"] is not None
        for arguments in service_arguments["validation"]
    )
    assert all(
        arguments["reference_profile_binding_service"] is binding_service
        for arguments in service_arguments["creation"]
    )
    assert all(
        arguments["reference_profile_resolver"] is reference_resolver
        for arguments in service_arguments["submission"]
    )
    assert harness._worker_environment()[
        platform_harness.REFERENCE_PROFILE_CONFIG_ENV
    ] == str(config_path)
    config_directory = config_path.parent
    assert harness._cleanup_reference_profile_config() is True
    assert not config_path.exists()
    assert not config_directory.exists()


@pytest.mark.parametrize("failure_stage", ("document", "open"))
def test_private_reference_config_failure_is_redacted_and_cleans_owned_directory(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    failure_stage: str,
) -> None:
    from . import platform_harness

    harness = _private_config_harness(tmp_path, monkeypatch)
    fixture = _acceptance_fixture(tmp_path)
    private_path = str(harness.reference_profile_config_path)
    if failure_stage == "document":
        monkeypatch.setattr(
            platform_harness,
            "_private_reference_profile_document",
            lambda _fixture: (_ for _ in ()).throw(OSError(private_path)),
        )
    else:
        original_open = platform_harness.os.open

        def fail_config_open(path, flags, mode=0o777):
            if Path(path) == harness.reference_profile_config_path:
                raise OSError(private_path)
            return original_open(path, flags, mode)

        monkeypatch.setattr(platform_harness.os, "open", fail_config_open)

    with pytest.raises(AssertionError) as captured:
        harness._prepare_reference_profile_config(fixture)

    assert (
        str(captured.value) == "private reference profile config could not be prepared"
    )
    assert private_path not in str(captured.value)
    assert not harness.reference_profile_config_path.exists()
    assert not harness.reference_profile_config_path.parent.exists()


def test_private_reference_config_collision_fails_closed_without_false_cleanup(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from . import platform_harness

    harness = _private_config_harness(tmp_path, monkeypatch)
    fixture = _acceptance_fixture(tmp_path)
    original_document = platform_harness._private_reference_profile_document

    def collide(current_fixture):
        harness.reference_profile_config_path.write_text(
            "untrusted claimant\n",
            encoding="utf-8",
        )
        return original_document(current_fixture)

    monkeypatch.setattr(
        platform_harness,
        "_private_reference_profile_document",
        collide,
    )

    with pytest.raises(AssertionError) as captured:
        harness._prepare_reference_profile_config(fixture)

    assert str(captured.value) == (
        "private reference profile config cleanup could not be confirmed"
    )
    assert str(harness.reference_profile_config_path) not in str(captured.value)
    assert harness._cleanup_reference_profile_config() is False
    harness.reference_profile_config_path.unlink()
    harness.reference_profile_config_path.parent.rmdir()


def test_private_reference_config_cleanup_rejects_identity_replacement(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    harness = _private_config_harness(tmp_path, monkeypatch)
    harness._prepare_reference_profile_config(_acceptance_fixture(tmp_path))
    config_path = harness.reference_profile_config_path
    replacement_path = config_path.parent / "replacement.json"
    replacement_path.write_text("replacement\n", encoding="utf-8")
    replacement_path.chmod(0o600)
    config_path.unlink()
    replacement_path.rename(config_path)

    assert harness._cleanup_reference_profile_config() is False
    assert config_path.read_text(encoding="utf-8") == "replacement\n"
    config_path.unlink()
    config_path.parent.rmdir()


def test_private_reference_config_cleanup_rejects_directory_replacement(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    harness = _private_config_harness(tmp_path, monkeypatch)
    harness._prepare_reference_profile_config(_acceptance_fixture(tmp_path))
    config_directory = harness.reference_profile_config_path.parent
    original_directory = config_directory.with_name("original-reference-profile")
    config_directory.rename(original_directory)
    config_directory.mkdir(mode=0o700)

    assert harness._cleanup_reference_profile_config() is False
    assert config_directory.is_dir()
    assert (original_directory / "reference-profiles.json").is_file()
    config_directory.rmdir()
    (original_directory / "reference-profiles.json").unlink()
    original_directory.rmdir()


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
