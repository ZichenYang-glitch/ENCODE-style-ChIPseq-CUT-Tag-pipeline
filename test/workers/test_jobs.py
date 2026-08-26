"""Integration tests for the PR123 worker ownership handshake."""

from __future__ import annotations

from dataclasses import replace
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import shlex
import traceback
from types import SimpleNamespace

import fakeredis
import pytest
from rq import SimpleWorker
from rq.serializers import JSONSerializer
from rq.timeouts import JobTimeoutException

import encode_pipeline.workers.jobs as worker_jobs
from encode_pipeline.adapters.bulk_rnaseq import (
    BulkRnaSeqExecutionBinding,
    BulkRnaSeqResultsWorkflowAdapter,
    BulkRnaSeqTranscriptomeBinding,
    RuntimeAssetBinding,
)
from encode_pipeline.adapters.bulk_rnaseq.execution_identity import (
    verify_execution_implementation,
)
from encode_pipeline.adapters.bulk_rnaseq.qualification import (
    load_default_execution_qualification,
)
from encode_pipeline.adapters.bulk_rnaseq.runtime_assets import (
    RuntimeAssetAdmission,
    VerifiedRuntimeAssets,
)
from encode_pipeline.persistence.runtime import open_run_persistence
from encode_pipeline.platform.adapters import WorkflowInputs
from encode_pipeline.platform.execution import RunExecutionAssignment, RunExecutionClaim
from encode_pipeline.platform.managed_containers import (
    managed_container_endpoint_identity,
)
from encode_pipeline.platform.notifications import DisabledTerminalRunNotifier
from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.platform.results import Issue, Result
from encode_pipeline.platform.result_generations import (
    validate_artifact_generation,
    validate_qc_generation,
    validate_result_attempt_id,
)
from encode_pipeline.platform.runs import RunStatus
from encode_pipeline.services.defaults import (
    create_default_run_service,
    create_default_workflow_registry,
)
from encode_pipeline.services.artifact_extraction import ArtifactExtractionService
from encode_pipeline.services.local_execution import LocalExecutionService
from encode_pipeline.services.private_reference_profiles import (
    load_private_reference_profile_config,
)
from encode_pipeline.services.process_runner import ProcessRunner
from encode_pipeline.services.qc_summary_indexing import QcSummaryIndexingService
from encode_pipeline.services.reference_profile_runtime import (
    ReferenceProfileBindingService,
)
from encode_pipeline.services.reference_profiles import ReferenceProfileService
from encode_pipeline.services.run_repositories import RunEventDraft
from encode_pipeline.services.runs import RunService
from encode_pipeline.services.validated_inputs import (
    ValidatedInputService,
    ValidatedRunCreationService,
)
from encode_pipeline.services.validation import ValidationService
from encode_pipeline.services.workflow_builds import WorkflowBuildIdentityProvider
from encode_pipeline.workers.jobs import (
    handle_execution_stopped,
    handle_work_horse_killed,
)
from encode_pipeline.workers.rq_queue import RqRunQueue
from encode_pipeline.workers.runtime import open_worker_runtime
from encode_pipeline.workers.settings import (
    MANAGED_DOCKER_EXECUTABLE_ENV,
    MANAGED_DOCKER_SOCKET_ENV,
    QUEUE_NAME_ENV,
    REFERENCE_PROFILE_CONFIG_ENV,
    REDIS_URL_ENV,
    WORKSPACE_ROOT_ENV,
)
from encode_pipeline.workers.timeouts import WorkerHardTimeout

from .conftest import create_planned_run, worker_settings


ARTIFACT_GENERATION = f"artifactgen-{'a' * 64}"
QC_GENERATION = f"qcgen-{'b' * 64}"


def _acquired_claim(run_id: str = "run-1") -> RunExecutionClaim:
    claimed_at = datetime.now(timezone.utc)
    return RunExecutionClaim(
        assignment=RunExecutionAssignment(
            run_id=run_id,
            job_id=f"job-{run_id}",
            backend="rq",
            queue_name="test-queue",
            created_at=claimed_at,
            dispatched_at=claimed_at,
            claimed_at=claimed_at,
        ),
        acquired=True,
    )


def _post_success_runtime(**services):
    terminal_notifier = services.pop(
        "terminal_notifier",
        DisabledTerminalRunNotifier(),
    )
    return SimpleNamespace(
        terminal_notifier=terminal_notifier,
        **services,
    )


class _ResultStateTracker:
    def __init__(self):
        self.artifact_attempt_id = None
        self.artifact_attempt_status = None
        self.artifact_outcome = None
        self.artifact_generation = None
        self.qc_attempt_id = None
        self.qc_attempt_status = None
        self.qc_attempt_artifact_generation = None
        self.qc_artifact_generation = None
        self.qc_generation = None
        self.qc_outcome = None

    def complete_artifact(self, attempt_id, result, *, status="succeeded"):
        self.artifact_attempt_id = attempt_id
        self.artifact_attempt_status = status
        self.artifact_outcome = status
        if status == "succeeded":
            self.artifact_generation = ARTIFACT_GENERATION
        return result

    def complete_qc(
        self, attempt_id, artifact_generation, result=None, *, status="succeeded"
    ):
        self.qc_attempt_id = attempt_id
        self.qc_attempt_status = status
        self.qc_attempt_artifact_generation = artifact_generation
        self.qc_artifact_generation = artifact_generation
        self.qc_generation = QC_GENERATION
        self.qc_outcome = status
        return result

    def get_result_state(self, _run_id):
        return self


def _verified_bulk_runtime_assets(root: Path) -> VerifiedRuntimeAssets:
    return VerifiedRuntimeAssets(
        root=root,
        source_tree=root / "source/rnaseq",
        nextflow_executable=root / "nextflow/nextflow-25.04.3-dist",
        jdk_archive=root / "jdk/corretto.tar.gz",
        jdk_tree=root / "jdk/corretto",
        java_executable=root / "jdk/corretto/bin/java",
        plugin_root=root / "plugins",
        plugin_archive=root / "plugins/nf-schema-2.5.1.zip",
        plugin_meta=root / "plugins/nf-schema-2.5.1-meta.json",
        plugin_tree=root / "plugins/nf-schema-2.5.1",
        container_lock=root / "containers/availability-lock.json",
        containers=(),
        source_tree_sha256="1" * 64,
        runtime_identity_sha256="2" * 64,
        nextflow_sha256="3" * 64,
        jdk_archive_sha256="8" * 64,
        jdk_tree_sha256="9" * 64,
        java_executable_sha256="a" * 64,
        plugin_archive_sha256="4" * 64,
        plugin_tree_sha256="5" * 64,
        container_inventory_sha256="6" * 64,
        container_lock_sha256="7" * 64,
    )


def _bulk_execution_adapter(
    runtime_root: Path,
    transcript_fasta: Path,
    *,
    reference_fasta_sha256: str,
    reference_gtf_sha256: str,
    implementation_qualification,
) -> BulkRnaSeqResultsWorkflowAdapter:
    transcript_contents = b">tx1\nACGT\n"
    transcript_fasta.write_bytes(transcript_contents)
    binding = BulkRnaSeqExecutionBinding(
        assets=RuntimeAssetBinding(root=runtime_root),
        transcriptome=BulkRnaSeqTranscriptomeBinding(
            reference_id="tiny",
            fasta_sha256=reference_fasta_sha256,
            gtf_sha256=reference_gtf_sha256,
            transcript_fasta=transcript_fasta,
            transcript_fasta_sha256=hashlib.sha256(transcript_contents).hexdigest(),
        ),
        implementation_qualification=implementation_qualification,
    )
    return BulkRnaSeqResultsWorkflowAdapter(execution=binding)


def _bulk_execution_inputs(root: Path) -> tuple[WorkflowInputs, str, str]:
    input_root = root / "bulk-inputs"
    input_root.mkdir(parents=True)
    reference_fasta = (input_root / "reference.fa").resolve()
    reference_gtf = (input_root / "reference.gtf").resolve()
    fastq = (input_root / "S1.fastq.gz").resolve()
    reference_fasta.write_bytes(b">chr1\nACGT\n")
    reference_gtf.write_bytes(b'chr1\ttest\texon\t1\t4\t.\t+\t.\tgene_id "g1";\n')
    fastq.write_bytes(b"controlled-fastq")
    fasta_sha256 = hashlib.sha256(reference_fasta.read_bytes()).hexdigest()
    gtf_sha256 = hashlib.sha256(reference_gtf.read_bytes()).hexdigest()
    return (
        WorkflowInputs(
            config={
                "standard": {
                    "reference": {
                        "reference_id": "tiny",
                        "fasta": str(reference_fasta),
                        "fasta_sha256": fasta_sha256,
                        "gtf": str(reference_gtf),
                        "gtf_sha256": gtf_sha256,
                    }
                }
            },
            samples=[
                {
                    "sample": "S1",
                    "library": "lib1",
                    "lane": "L001",
                    "layout": "SE",
                    "fastq_1": str(fastq),
                    "strandedness": "auto",
                    "platform": "ILLUMINA",
                }
            ],
            options={},
        ),
        fasta_sha256,
        gtf_sha256,
    )


def _configure_worker_environment(
    monkeypatch,
    configured,
    *,
    write_qc_summary: bool = False,
):
    monkeypatch.setenv("ENCODE_PIPELINE_DATABASE_URL", configured.database_url)
    monkeypatch.setenv(REDIS_URL_ENV, configured.redis_url)
    monkeypatch.setenv(QUEUE_NAME_ENV, configured.queue_name)
    monkeypatch.setenv(WORKSPACE_ROOT_ENV, str(configured.workspace_root))
    assert configured.reference_profile_config is not None
    monkeypatch.setenv(
        REFERENCE_PROFILE_CONFIG_ENV,
        str(configured.reference_profile_config),
    )
    executable_dir = configured.workspace_root.parent / "test-bin"
    executable_dir.mkdir(parents=True, exist_ok=True)
    qc_commands = ""
    if write_qc_summary:
        from encode_pipeline.adapters.encode_qc import _QC_HEADER

        row = {column: "NA" for column in _QC_HEADER}
        row.update(
            {
                "sample": "S1",
                "assay": "chipseq",
                "total_reads": "1000",
                "frip": "0.125",
                "peak_count": "50",
                "percent_duplication": "0.2",
                "estimated_library_size": "9007199254740993",
                "nrf": "0.8",
                "pbc1": "0.75",
                "pbc2": "3.0",
            }
        )
        contents = "\t".join(_QC_HEADER) + "\n"
        contents += "\t".join(row[column] for column in _QC_HEADER) + "\n"
        qc_commands = (
            "mkdir -p results/S1/01_qc\n"
            f"printf %s {shlex.quote(contents)} > "
            "results/S1/01_qc/S1.qc_summary.tsv\n"
        )
    snakemake = executable_dir / "snakemake"
    snakemake.write_text(
        "#!/bin/sh\n"
        'while [ "$#" -gt 0 ]; do\n'
        '  if [ "$1" = "--directory" ]; then shift; cd "$1"; break; fi\n'
        "  shift\n"
        "done\n"
        "mkdir -p results/multiqc\n"
        "printf 'output_type\\tpath\\n' > results/multiqc/result_manifest.tsv\n"
        f"{qc_commands}"
        "printf 'worker stdout\\n'\n"
        "printf 'worker stderr\\n' >&2\n",
        encoding="utf-8",
    )
    snakemake.chmod(0o755)
    monkeypatch.setenv("PATH", f"{executable_dir}{os.pathsep}{os.environ['PATH']}")


def _run_burst(connection, run_queue):
    worker = SimpleWorker(
        [run_queue._queue],
        connection=connection,
        serializer=JSONSerializer,
    )
    return worker.work(burst=True, logging_level="WARNING")


def _read_run(configured, run_id):
    persistence = open_run_persistence(configured.database_url)
    service = create_default_run_service(
        registry=create_default_workflow_registry(),
        repository=persistence.repository,
    )
    try:
        return (
            service.get_run(run_id),
            service.list_events(run_id, limit=100),
            service.get_execution_assignment(run_id),
        )
    finally:
        persistence.close()


def _prepare_requeue(configured, assignment):
    persistence = open_run_persistence(configured.database_url)
    service = RunService(
        create_default_workflow_registry(),
        repository=persistence.repository,
    )
    try:
        service.mark_execution_dispatched(
            assignment.run_id,
            job_id=assignment.job_id,
        )
        service.queue_dispatched_run(
            assignment.run_id,
            job_id=assignment.job_id,
            backend=assignment.backend,
            queue_name=assignment.queue_name,
        )
        current = service.get_execution_assignment(assignment.run_id)
        assert current is not None
        return persistence.repository.prepare_execution_requeue(
            assignment.run_id,
            expected_status=RunStatus.QUEUED,
            expected_assignment=current,
            requested_at=datetime.now(timezone.utc),
            event=RunEventDraft(
                event_type="run_requeue_requested_by_admin",
                message="Run requeue requested by an administrator.",
                status=RunStatus.QUEUED,
                stage="execution",
                context={"reason_code": "RUN_REQUEUED_BY_ADMIN_RECOVERY"},
            ),
        ).assignment
    finally:
        persistence.close()


def _read_results(configured, run_id):
    persistence = open_run_persistence(configured.database_url)
    service = create_default_run_service(
        registry=create_default_workflow_registry(),
        repository=persistence.repository,
    )
    try:
        return (
            service.get_result_state(run_id),
            service.list_artifacts(run_id),
            service.list_qc_metrics(run_id),
        )
    finally:
        persistence.close()


def _complete_workflow_successfully(execution_service, run_id, _claim):
    execution_service._run_service.transition_run(
        run_id,
        RunStatus.RUNNING,
        stage="execution",
        message="Worker started local workflow execution.",
        context={"reason_code": "LOCAL_EXECUTION_STARTED"},
    )
    completed = execution_service._run_service.transition_run(
        run_id,
        RunStatus.SUCCEEDED,
        stage="execution",
        message="Local workflow execution completed successfully.",
        context={"reason_code": "LOCAL_EXECUTION_SUCCEEDED"},
    )
    return Result.success(completed)


def _copy_controlled_project(destination: Path) -> None:
    source = Path(__file__).resolve().parents[2]
    destination.mkdir()
    shutil.copy2(source / "pyproject.toml", destination / "pyproject.toml")
    for relative in (
        Path("docs/architecture/artifact-inventory.yaml"),
        Path("src/encode_pipeline"),
        Path("workflow"),
        Path("profiles/default"),
        Path("scripts"),
    ):
        if (source / relative).is_dir():
            shutil.copytree(source / relative, destination / relative)
        else:
            (destination / relative).parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source / relative, destination / relative)


def _fail_if_process_starts(*_args, **_kwargs):
    raise AssertionError("ProcessRunner must not run before build verification")


def test_rq_worker_rebuilds_dependencies_and_persists_handshake_event(
    tmp_path, monkeypatch
):
    configured = worker_settings(tmp_path)
    assignment = create_planned_run(
        configured,
        "handshake-run",
        assign_queue=configured.queue_name,
    )
    _configure_worker_environment(monkeypatch, configured)
    connection = fakeredis.FakeRedis()
    run_queue = RqRunQueue(configured, connection=connection)
    assert assignment is not None
    run_queue.enqueue_execution(assignment)

    assert _run_burst(connection, run_queue) is True

    record, events, persisted_assignment = _read_run(configured, "handshake-run")
    handshake = [
        event for event in events if event.event_type == "worker_dependencies_rebuilt"
    ]
    assert record.status is RunStatus.SUCCEEDED
    assert persisted_assignment is not None
    assert persisted_assignment.dispatched_at is not None
    assert persisted_assignment.claimed_at is not None
    assert len(handshake) == 1
    assert handshake[0].context == {
        "backend": "rq",
        "job_id": assignment.job_id,
        "queue_name": configured.queue_name,
        "workflow_id": record.workflow_id,
    }
    persistence = open_run_persistence(configured.database_url)
    try:
        service = create_default_run_service(
            registry=create_default_workflow_registry(),
            repository=persistence.repository,
        )
        assert service.list_logs("handshake-run", "stdout")[0].lines == (
            "worker stdout",
        )
        assert service.list_logs("handshake-run", "stderr")[0].lines == (
            "worker stderr",
        )
        artifacts = service.list_artifacts("handshake-run")
        assert len(artifacts) == 1
        assert artifacts[0].metadata["output_type"] == "result_manifest"
        assert artifacts[0].metadata["relative_path"] == (
            "results/multiqc/result_manifest.tsv"
        )
        assert service.list_qc_metrics("handshake-run") == ()
        qc_events = [
            event for event in events if event.event_type == "qc_metrics_indexed"
        ]
        assert len(qc_events) == 1
        assert qc_events[0].context["metric_count"] == 0
        validate_result_attempt_id(qc_events[0].context["attempt_id"])
        validate_artifact_generation(qc_events[0].context["artifact_generation"])
        validate_qc_generation(qc_events[0].context["qc_generation"])
    finally:
        persistence.close()

    job = run_queue._queue.fetch_job(assignment.job_id)
    assert job is not None
    job.delete()
    run_queue.enqueue_execution(assignment)
    assert _run_burst(connection, run_queue) is True

    _, events_after_retry, assignment_after_retry = _read_run(
        configured,
        "handshake-run",
    )
    assert assignment_after_retry == persisted_assignment
    assert [event.event_type for event in events_after_retry].count(
        "worker_dependencies_rebuilt"
    ) == 1


def test_rq_worker_does_not_resume_preclaimed_legacy_null_cleanup_assignment(
    tmp_path,
    monkeypatch,
):
    configured = worker_settings(tmp_path)
    assignment = create_planned_run(
        configured,
        "preclaimed-legacy-run",
        assign_queue=configured.queue_name,
    )
    assert assignment is not None
    persistence = open_run_persistence(configured.database_url)
    try:
        service = create_default_run_service(
            registry=create_default_workflow_registry(),
            repository=persistence.repository,
        )
        service.mark_execution_dispatched(
            assignment.run_id,
            job_id=assignment.job_id,
        )
        service.queue_dispatched_run(
            assignment.run_id,
            job_id=assignment.job_id,
            backend=assignment.backend,
            queue_name=assignment.queue_name,
        )
        claim = service.claim_execution_assignment(
            assignment.run_id,
            job_id=assignment.job_id,
            backend=assignment.backend,
            queue_name=assignment.queue_name,
        )
        assert claim.acquired is True
        assert claim.assignment.managed_container_scope is None
        assert claim.assignment.managed_container_endpoint_identity is None
    finally:
        persistence.close()

    before_record, before_events, before_assignment = _read_run(
        configured,
        assignment.run_id,
    )

    def fail_if_execution_starts(*_args, **_kwargs):
        raise AssertionError("pre-claimed jobs must stop before execution rebuild")

    monkeypatch.setattr(LocalExecutionService, "execute", fail_if_execution_starts)
    monkeypatch.setattr(ProcessRunner, "run", _fail_if_process_starts)
    _configure_worker_environment(monkeypatch, configured)
    connection = fakeredis.FakeRedis()
    run_queue = RqRunQueue(configured, connection=connection)
    run_queue.enqueue_execution(claim.assignment)

    assert _run_burst(connection, run_queue) is True

    job = run_queue._queue.fetch_job(assignment.job_id)
    record, events, persisted_assignment = _read_run(
        configured,
        assignment.run_id,
    )
    assert job is not None
    assert job.is_finished
    assert record == before_record
    assert events == before_events
    assert persisted_assignment == before_assignment


def test_rq_worker_persists_nonempty_qc_metrics_for_new_sqlite_reader(
    tmp_path,
    monkeypatch,
):
    configured = worker_settings(tmp_path)
    samples = tmp_path / "qc-samples.tsv"
    samples.write_text(
        "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\tpeak_mode\tgenome\t"
        "bowtie2_index\texperiment\tcondition\treplicate\tbiological_replicate\t"
        "technical_replicate\trole\tcontrol_sample\tcontrol_bam\n"
        "S1\t/tmp/R1.fq.gz\t/tmp/R2.fq.gz\tPE\tchipseq\tCTCF\tnarrow\ths\t"
        "/tmp/index\tEXP1\ttreatment\t1\t1\t1\ttreatment\t\t\n",
        encoding="utf-8",
    )
    assignment = create_planned_run(
        configured,
        "qc-run",
        assign_queue=configured.queue_name,
        enable_qc_summary=True,
        samples_path=samples,
    )
    assert assignment is not None
    _configure_worker_environment(monkeypatch, configured, write_qc_summary=True)
    connection = fakeredis.FakeRedis()
    run_queue = RqRunQueue(configured, connection=connection)
    run_queue.enqueue_execution(assignment)

    assert _run_burst(connection, run_queue) is True

    persistence = open_run_persistence(configured.database_url)
    try:
        service = create_default_run_service(
            registry=create_default_workflow_registry(),
            repository=persistence.repository,
        )
        record = service.get_run("qc-run")
        metrics = service.list_qc_metrics("qc-run")
        events = service.list_events("qc-run", limit=100)
    finally:
        persistence.close()

    assert record.status is RunStatus.SUCCEEDED
    assert {metric.metric_key for metric in metrics} == {
        "library.estimated_size",
        "library.nrf",
        "library.pbc1",
        "library.pbc2",
        "library.percent_duplication",
        "peaks.count",
        "peaks.frip",
        "sequencing.total_reads",
    }
    assert (
        next(
            metric
            for metric in metrics
            if metric.metric_key == "library.estimated_size"
        ).value
        == 9007199254740993
    )
    qc_event = next(
        event for event in events if event.event_type == "qc_metrics_indexed"
    )
    assert qc_event.context["metric_count"] == 8
    validate_result_attempt_id(qc_event.context["attempt_id"])
    validate_artifact_generation(qc_event.context["artifact_generation"])
    validate_qc_generation(qc_event.context["qc_generation"])


def test_bulk_binding_survives_durable_worker_reconstruction_and_admission(
    tmp_path,
    monkeypatch,
):
    configured = worker_settings(tmp_path)
    runtime_root = (tmp_path / "bulk-runtime").resolve()
    transcript_fasta = (tmp_path / "transcripts.fa").resolve()
    inputs, reference_fasta_sha256, reference_gtf_sha256 = _bulk_execution_inputs(
        tmp_path
    )
    verified_assets = _verified_bulk_runtime_assets(runtime_root)
    implementation = verify_execution_implementation()
    assert implementation.is_success
    qualification = load_default_execution_qualification(implementation.value)
    assert qualification.is_success
    monkeypatch.setattr(
        RuntimeAssetAdmission,
        "acquire",
        lambda _self: Result.success(verified_assets),
    )
    api_adapter = _bulk_execution_adapter(
        runtime_root,
        transcript_fasta,
        reference_fasta_sha256=reference_fasta_sha256,
        reference_gtf_sha256=reference_gtf_sha256,
        implementation_qualification=qualification.value.implementation,
    )
    worker_adapter = _bulk_execution_adapter(
        runtime_root,
        transcript_fasta,
        reference_fasta_sha256=reference_fasta_sha256,
        reference_gtf_sha256=reference_gtf_sha256,
        implementation_qualification=qualification.value.implementation,
    )
    api_registry = WorkflowRegistry((api_adapter,))
    worker_registry = WorkflowRegistry((worker_adapter,))
    api_provider = WorkflowBuildIdentityProvider(api_registry)
    worker_provider = WorkflowBuildIdentityProvider(worker_registry)
    api_identity = api_provider.capture_executable("bulk-rnaseq")
    worker_identity = worker_provider.capture_executable("bulk-rnaseq")
    assert api_identity.is_success
    assert worker_identity.is_success
    assert api_identity.value.matches(worker_identity.value)
    reference = dict(inputs.config["standard"]["reference"])
    transcriptome = {
        "reference_id": reference["reference_id"],
        "fasta_sha256": reference["fasta_sha256"],
        "gtf_sha256": reference["gtf_sha256"],
        "transcript_fasta": str(transcript_fasta),
        "transcript_fasta_sha256": hashlib.sha256(
            transcript_fasta.read_bytes()
        ).hexdigest(),
    }
    private_config_path = (tmp_path / "bulk-reference-profile.json").resolve()
    private_config_path.write_text(
        json.dumps(
            {
                "schema_version": "helixweave-reference-profiles-v1",
                "profiles": {
                    "bulk-test-reference": {
                        "bindings": {
                            "bulk-rnaseq": {
                                "schema_version": "bulk-rnaseq-reference-binding-v1",
                                "reference": reference,
                                "transcriptome": transcriptome,
                            }
                        }
                    }
                },
            }
        ),
        encoding="utf-8",
    )
    configured = replace(
        configured,
        reference_profile_config=private_config_path,
    )
    persistence = open_run_persistence(configured.database_url)
    try:

        def private_config_provider():
            return load_private_reference_profile_config(private_config_path)

        profiles = ReferenceProfileService(
            repository=persistence.reference_profile_repository,
            private_config_provider=private_config_provider,
            adapter_provider=api_registry.get,
        )
        registered = profiles.register(
            safe_key="bulk-test-reference",
            display_name="Bulk test reference",
            organism="Test organism",
            assembly="tiny",
            config_key="bulk-test-reference",
        )
        enabled = profiles.enable(
            registered.profile_id,
            revision_id=registered.revision_id,
        )
        bindings = ReferenceProfileBindingService(
            repository=persistence.reference_profile_repository,
            private_config_provider=private_config_provider,
            adapter_provider=api_registry.get,
        )
        service = RunService(
            api_registry,
            id_factory=lambda: "bulk-worker-run",
            repository=persistence.repository,
        )
        standard = dict(inputs.config["standard"])
        standard.pop("reference")
        unbound_inputs = WorkflowInputs(
            config={**inputs.config, "standard": standard},
            samples=inputs.samples,
            options=inputs.options,
        )
        snapshot_result = ValidatedInputService(
            registry=api_registry,
            validation_service=ValidationService(api_registry),
            build_identity_provider=api_provider,
            repository=persistence.repository,
            reference_profile_binding_service=bindings,
            reference_profile_catalog=profiles,
        ).validate(
            "bulk-rnaseq",
            unbound_inputs,
            reference_profile_revision_id=enabled.revision_id,
        )
        assert snapshot_result.is_success and snapshot_result.value is not None
        record = (
            ValidatedRunCreationService(
                run_service=service,
                build_identity_provider=api_provider,
                reference_profile_binding_service=bindings,
            )
            .create_run("bulk-rnaseq", snapshot_result.value.snapshot_id)
            .record
        )
        service.transition_run(record.run_id, RunStatus.VALIDATING)
        service.complete_preflight(
            record.run_id,
            snapshot_result.value.workflow_build_identity,
        )
        assignment = service.ensure_execution_assignment(
            record.run_id,
            queue_name=configured.queue_name,
        )
    finally:
        persistence.close()

    worker_runner = ProcessRunner(allowed_executables=("/usr/bin/unshare",))
    with open_worker_runtime(
        configured,
        registry=worker_registry,
        build_identity_provider=worker_provider,
        process_runner=worker_runner,
    ) as runtime:
        reconstructed = runtime.build_identity_provider.capture_executable(
            "bulk-rnaseq"
        )
        assert reconstructed.is_success
        assert api_identity.value.matches(reconstructed.value)
        claim = worker_jobs._initialize_execution_with_runtime(
            runtime,
            SimpleNamespace(
                id=assignment.job_id,
                origin=configured.queue_name,
            ),
            assignment.run_id,
        )
        assert claim is not None
        assert claim.acquired is True
        assert runtime.run_service.get_run(assignment.run_id).status is RunStatus.QUEUED
        persisted_assignment = runtime.run_service.get_execution_assignment(
            assignment.run_id
        )
        assert persisted_assignment is not None
        assert persisted_assignment.claimed_at is not None
        assert claim.assignment == persisted_assignment
        events = runtime.run_service.list_events(assignment.run_id)
        assert events[-1].event_type == "worker_dependencies_rebuilt"
        assert events[-1].context["workflow_id"] == "bulk-rnaseq"
        durable_plan = runtime.execution_planner.plan_run(assignment.run_id)
        assert durable_plan.is_success
        workspace = (tmp_path / "worker-rebuilt-workspace").resolve()
        workspace_plan = runtime.workspace_planner.plan_workspace(
            durable_plan.value,
            base_dir=workspace,
        )
        assert workspace_plan.is_success
        command = runtime.command_builder.build_command(
            workspace_plan.value,
            workspace,
        )
        assert command.is_success
        assert command.value.command_spec.argv[0] == "/usr/bin/unshare"
        assert worker_runner._admit_executable(command.value.command_spec).is_success


def test_rq_worker_rejects_stale_job_identity_without_writing_handshake(
    tmp_path, monkeypatch
):
    configured = worker_settings(tmp_path)
    assignment = create_planned_run(
        configured,
        "stale-run",
        assign_queue=configured.queue_name,
    )
    _configure_worker_environment(monkeypatch, configured)
    connection = fakeredis.FakeRedis()
    run_queue = RqRunQueue(configured, connection=connection)
    assert assignment is not None
    stale_assignment = replace(assignment, job_id="stale-job-id")
    run_queue.enqueue_execution(stale_assignment)

    assert _run_burst(connection, run_queue) is True

    job = run_queue._queue.fetch_job("stale-job-id")
    _, events, _ = _read_run(configured, "stale-run")
    assert job is not None
    assert job.is_failed
    assert all(event.event_type != "worker_dependencies_rebuilt" for event in events)


def test_worker_entry_durably_confirms_pending_requeue_before_build_admission(
    tmp_path,
    monkeypatch,
):
    configured = worker_settings(tmp_path)
    assignment = create_planned_run(
        configured,
        "pending-requeue-entry",
        assign_queue=configured.queue_name,
    )
    assert assignment is not None
    prepared = _prepare_requeue(configured, assignment)
    monkeypatch.setattr(
        worker_jobs,
        "_require_matching_workflow_build",
        lambda *_args: False,
    )

    with open_worker_runtime(configured) as runtime:
        claim = worker_jobs._initialize_execution_with_runtime(
            runtime,
            SimpleNamespace(id=prepared.job_id, origin=prepared.queue_name),
            prepared.run_id,
        )

        persisted = runtime.run_service.get_execution_assignment(prepared.run_id)
        events = runtime.run_service.list_events(prepared.run_id, limit=100)

    assert claim is None
    assert persisted is not None
    assert persisted.claimed_at is None
    assert persisted.requeue_requested_at == prepared.requeue_requested_at
    assert persisted.requeue_confirmed_at is not None
    assert events[-1].event_type == "run_requeue_delivery_observed_by_worker"
    assert events[-1].context == {
        "reason_code": "RUN_REQUEUED_BY_ADMIN_RECOVERY",
        "confirmation_source": "worker_entry",
    }


def test_worker_requeue_confirmation_cas_failure_stops_before_build_or_claim(
    tmp_path,
    monkeypatch,
):
    configured = worker_settings(tmp_path)
    assignment = create_planned_run(
        configured,
        "pending-requeue-cas-failure",
        assign_queue=configured.queue_name,
    )
    assert assignment is not None
    prepared = _prepare_requeue(configured, assignment)

    with open_worker_runtime(configured) as runtime:
        monkeypatch.setattr(
            runtime.run_service,
            "confirm_execution_requeue_observed",
            lambda *_args, **_kwargs: (_ for _ in ()).throw(
                RuntimeError("private persistence failure")
            ),
        )
        monkeypatch.setattr(
            worker_jobs,
            "_require_matching_workflow_build",
            lambda *_args: (_ for _ in ()).throw(
                AssertionError("build admission must not run")
            ),
        )

        with pytest.raises(RuntimeError, match="private persistence failure"):
            worker_jobs._initialize_execution_with_runtime(
                runtime,
                SimpleNamespace(id=prepared.job_id, origin=prepared.queue_name),
                prepared.run_id,
            )

        persisted = runtime.run_service.get_execution_assignment(prepared.run_id)
    assert persisted is not None
    assert persisted.claimed_at is None
    assert persisted.requeue_confirmed_at is None


def test_rq_worker_rejects_missing_durable_assignment(tmp_path, monkeypatch):
    configured = worker_settings(tmp_path)
    create_planned_run(configured, "unassigned-run")
    _configure_worker_environment(monkeypatch, configured)
    connection = fakeredis.FakeRedis()
    run_queue = RqRunQueue(configured, connection=connection)
    orphan_assignment = RunExecutionAssignment(
        run_id="unassigned-run",
        job_id="orphan-job-id",
        backend="rq",
        queue_name=configured.queue_name,
        created_at=datetime.now(timezone.utc),
    )
    run_queue.enqueue_execution(orphan_assignment)

    assert _run_burst(connection, run_queue) is True

    job = run_queue._queue.fetch_job("orphan-job-id")
    _, events, _ = _read_run(configured, "unassigned-run")
    assert job is not None
    assert job.is_failed
    assert all(event.event_type != "worker_dependencies_rebuilt" for event in events)


def test_rq_worker_fails_legacy_planned_run_without_build_before_claim(
    tmp_path,
    monkeypatch,
):
    configured = worker_settings(tmp_path)
    assignment = create_planned_run(
        configured,
        "missing-build-run",
        assign_queue=configured.queue_name,
        bind_build_identity=False,
    )
    assert assignment is not None
    _configure_worker_environment(monkeypatch, configured)
    monkeypatch.setattr(ProcessRunner, "run", _fail_if_process_starts)
    connection = fakeredis.FakeRedis()
    run_queue = RqRunQueue(configured, connection=connection)
    run_queue.enqueue_execution(assignment)
    before_record, before_events, before_assignment = _read_run(
        configured,
        assignment.run_id,
    )

    assert _run_burst(connection, run_queue) is True

    job = run_queue._queue.fetch_job(assignment.job_id)
    record, events, persisted_assignment = _read_run(
        configured,
        assignment.run_id,
    )
    assert job is not None
    assert job.is_failed
    assert record == before_record
    assert events == before_events
    assert persisted_assignment == before_assignment
    assert all(event.event_type != "worker_dependencies_rebuilt" for event in events)


def test_rq_worker_rejects_legacy_reference_run_without_revision_evidence(
    tmp_path,
    monkeypatch,
):
    configured = worker_settings(tmp_path)
    persistence = open_run_persistence(configured.database_url)
    try:
        registry = create_default_workflow_registry()
        service = RunService(
            registry,
            id_factory=lambda: "legacy-reference-run",
            repository=persistence.repository,
        )
        record = service.create_run(
            "encode-style-chipseq-cuttag-atac-mnase",
            WorkflowInputs(config={}),
        )
        service.transition_run(record.run_id, RunStatus.VALIDATING)
        build = WorkflowBuildIdentityProvider(registry).capture_executable(
            record.workflow_id
        )
        assert build.is_success and build.value is not None
        service.complete_preflight(record.run_id, build.value)
        assignment = service.ensure_execution_assignment(
            record.run_id,
            queue_name=configured.queue_name,
        )
    finally:
        persistence.close()

    with open_worker_runtime(configured) as runtime:
        admission = worker_jobs._capture_current_workflow_build(runtime, record)
    assert admission.is_failure
    assert admission.issues[0].code == "RUN_WORKFLOW_BUILD_IDENTITY_UNAVAILABLE"

    _configure_worker_environment(monkeypatch, configured)
    monkeypatch.setattr(ProcessRunner, "run", _fail_if_process_starts)
    connection = fakeredis.FakeRedis()
    run_queue = RqRunQueue(configured, connection=connection)
    run_queue.enqueue_execution(assignment)

    assert _run_burst(connection, run_queue) is True

    job = run_queue._queue.fetch_job(assignment.job_id)
    failed, events, persisted_assignment = _read_run(
        configured,
        assignment.run_id,
    )
    assert job is not None and job.is_failed
    assert failed.status is RunStatus.PLANNED
    assert failed.error is None
    assert persisted_assignment is not None
    assert persisted_assignment.claimed_at is None
    assert all(event.event_type != "worker_dependencies_rebuilt" for event in events)


def test_rq_worker_rejects_project_a_build_on_project_b_before_process(
    tmp_path,
    monkeypatch,
):
    configured = worker_settings(tmp_path)
    project_a = tmp_path / "project-a"
    _copy_controlled_project(project_a)
    snakefile_a = project_a / "workflow" / "Snakefile"
    snakefile_a.write_text(
        snakefile_a.read_text(encoding="utf-8") + "\n# project A drift\n",
        encoding="utf-8",
    )
    registry = create_default_workflow_registry()
    identity_result = WorkflowBuildIdentityProvider(
        registry,
        project_root=project_a,
    ).capture("encode-style-chipseq-cuttag-atac-mnase")
    assert identity_result.is_success
    assignment = create_planned_run(
        configured,
        "mismatched-build-run",
        assign_queue=configured.queue_name,
        build_identity=identity_result.value,
    )
    assert assignment is not None
    _configure_worker_environment(monkeypatch, configured)
    monkeypatch.setattr(ProcessRunner, "run", _fail_if_process_starts)
    connection = fakeredis.FakeRedis()
    run_queue = RqRunQueue(configured, connection=connection)
    run_queue.enqueue_execution(assignment)
    before_record, before_events, before_assignment = _read_run(
        configured,
        assignment.run_id,
    )

    assert _run_burst(connection, run_queue) is True

    job = run_queue._queue.fetch_job(assignment.job_id)
    record, events, persisted_assignment = _read_run(
        configured,
        assignment.run_id,
    )
    assert job is not None
    assert job.is_failed
    assert record == before_record
    assert events == before_events
    assert persisted_assignment == before_assignment
    assert all(event.event_type != "worker_dependencies_rebuilt" for event in events)


def test_rq_worker_rejects_adapter_version_drift_before_any_worker_side_effect(
    tmp_path,
    monkeypatch,
):
    configured = worker_settings(tmp_path)
    registry = create_default_workflow_registry()
    identity_result = WorkflowBuildIdentityProvider(registry).capture(
        "encode-style-chipseq-cuttag-atac-mnase"
    )
    assert identity_result.is_success
    current_identity = identity_result.value
    assignment = create_planned_run(
        configured,
        "adapter-version-drift-run",
        assign_queue=configured.queue_name,
        build_identity=replace(
            current_identity,
            adapter_version=f"{current_identity.adapter_version}-stale",
        ),
    )
    assert assignment is not None
    _configure_worker_environment(monkeypatch, configured)

    before_record, before_events, before_assignment = _read_run(
        configured,
        assignment.run_id,
    )
    assert before_record.status is RunStatus.PLANNED
    assert before_assignment == assignment
    assert before_assignment.dispatched_at is None
    assert before_assignment.claimed_at is None

    def fail_if_execution_rebuild_starts(*_args, **_kwargs):
        raise AssertionError(
            "LocalExecutionService must not rebuild or inspect the workspace "
            "before adapter admission"
        )

    monkeypatch.setattr(
        LocalExecutionService,
        "execute",
        fail_if_execution_rebuild_starts,
    )
    monkeypatch.setattr(ProcessRunner, "run", _fail_if_process_starts)
    connection = fakeredis.FakeRedis()
    run_queue = RqRunQueue(configured, connection=connection)
    run_queue.enqueue_execution(assignment)

    assert _run_burst(connection, run_queue) is True

    job = run_queue._queue.fetch_job(assignment.job_id)
    record, events, persisted_assignment = _read_run(
        configured,
        assignment.run_id,
    )
    assert job is not None
    assert job.is_failed
    assert record == before_record
    assert events == before_events
    assert persisted_assignment == before_assignment
    assert all(event.event_type != "worker_dependencies_rebuilt" for event in events)


def test_stale_worker_cannot_fail_queued_run_before_assignment_check(
    tmp_path,
    monkeypatch,
):
    configured = worker_settings(tmp_path)
    registry = create_default_workflow_registry()
    identity_result = WorkflowBuildIdentityProvider(registry).capture_executable(
        "encode-style-chipseq-cuttag-atac-mnase"
    )
    assert identity_result.is_success
    assignment = create_planned_run(
        configured,
        "stale-build-drift-run",
        assign_queue=configured.queue_name,
        build_identity=replace(
            identity_result.value,
            adapter_version=f"{identity_result.value.adapter_version}-stale",
        ),
    )
    assert assignment is not None
    persistence = open_run_persistence(configured.database_url)
    try:
        service = create_default_run_service(
            registry=create_default_workflow_registry(),
            repository=persistence.repository,
        )
        service.mark_execution_dispatched(
            assignment.run_id,
            job_id=assignment.job_id,
        )
        service.queue_dispatched_run(
            assignment.run_id,
            job_id=assignment.job_id,
            backend=assignment.backend,
            queue_name=assignment.queue_name,
        )
    finally:
        persistence.close()
    _configure_worker_environment(monkeypatch, configured)
    monkeypatch.setattr(ProcessRunner, "run", _fail_if_process_starts)
    connection = fakeredis.FakeRedis()
    run_queue = RqRunQueue(configured, connection=connection)
    stale_assignment = replace(assignment, job_id="stale-build-drift-job")
    run_queue.enqueue_execution(stale_assignment)
    before_record, before_events, before_assignment = _read_run(
        configured,
        assignment.run_id,
    )
    assert before_record.status is RunStatus.QUEUED

    assert _run_burst(connection, run_queue) is True

    job = run_queue._queue.fetch_job(stale_assignment.job_id)
    record, events, persisted_assignment = _read_run(
        configured,
        assignment.run_id,
    )
    assert job is not None
    assert job.is_failed
    assert record == before_record
    assert events == before_events
    assert persisted_assignment == before_assignment


def test_rq_worker_fails_closed_when_local_build_cannot_be_fingerprinted(
    tmp_path,
    monkeypatch,
):
    configured = worker_settings(tmp_path)
    assignment = create_planned_run(
        configured,
        "unavailable-build-run",
        assign_queue=configured.queue_name,
    )
    assert assignment is not None
    _configure_worker_environment(monkeypatch, configured)

    def unavailable(_self, _adapter):
        return Result.failure(
            [
                Issue(
                    code="WORKFLOW_BUILD_SOURCE_UNAVAILABLE",
                    message="Controlled source is unavailable.",
                    severity="error",
                    source="test",
                )
            ]
        )

    monkeypatch.setattr(
        WorkflowBuildIdentityProvider,
        "capture_resolved_executable",
        unavailable,
    )
    monkeypatch.setattr(ProcessRunner, "run", _fail_if_process_starts)
    connection = fakeredis.FakeRedis()
    run_queue = RqRunQueue(configured, connection=connection)
    run_queue.enqueue_execution(assignment)
    before_record, before_events, before_assignment = _read_run(
        configured,
        assignment.run_id,
    )

    assert _run_burst(connection, run_queue) is True

    job = run_queue._queue.fetch_job(assignment.job_id)
    record, events, persisted_assignment = _read_run(
        configured,
        assignment.run_id,
    )
    assert job is not None
    assert job.is_failed
    assert record == before_record
    assert events == before_events
    assert persisted_assignment == before_assignment
    assert all(event.event_type != "worker_dependencies_rebuilt" for event in events)


def test_rq_worker_treats_queued_cancellation_as_clean_noop_without_process(
    tmp_path, monkeypatch
):
    configured = worker_settings(tmp_path)
    assignment = create_planned_run(
        configured,
        "cancelled-run",
        assign_queue=configured.queue_name,
    )
    assert assignment is not None
    _configure_worker_environment(monkeypatch, configured)
    connection = fakeredis.FakeRedis()
    run_queue = RqRunQueue(configured, connection=connection)
    run_queue.enqueue_execution(assignment)

    persistence = open_run_persistence(configured.database_url)
    try:
        run_service = create_default_run_service(
            registry=create_default_workflow_registry(),
            repository=persistence.repository,
        )
        run_service.mark_execution_dispatched(
            assignment.run_id,
            job_id=assignment.job_id,
        )
        queued = run_service.queue_dispatched_run(
            assignment.run_id,
            job_id=assignment.job_id,
            backend=assignment.backend,
            queue_name=assignment.queue_name,
        )
        assert queued.status is RunStatus.QUEUED
        run_service.cancel_run("cancelled-run", reason="Cancelled before worker claim.")
    finally:
        persistence.close()

    def fail_if_process_starts(*_args, **_kwargs):
        raise AssertionError("ProcessRunner must not run for a cancelled queued job")

    monkeypatch.setattr(ProcessRunner, "run", fail_if_process_starts)

    assert _run_burst(connection, run_queue) is True

    job = run_queue._queue.fetch_job(assignment.job_id)
    record, events, _ = _read_run(configured, "cancelled-run")
    assert job is not None
    assert job.is_finished
    assert record.status is RunStatus.CANCELLED
    assert all(event.event_type != "worker_dependencies_rebuilt" for event in events)
    persistence = open_run_persistence(configured.database_url)
    try:
        run_service = create_default_run_service(
            registry=create_default_workflow_registry(),
            repository=persistence.repository,
        )
        assert run_service.list_logs("cancelled-run", "stdout") == ()
        assert run_service.list_logs("cancelled-run", "stderr") == ()
    finally:
        persistence.close()


def test_rq_worker_persists_nonzero_failure_when_notification_hits_deadline(
    tmp_path,
    monkeypatch,
):
    configured = worker_settings(tmp_path)
    assignment = create_planned_run(
        configured,
        "nonzero-run",
        assign_queue=configured.queue_name,
    )
    assert assignment is not None
    _configure_worker_environment(monkeypatch, configured)
    notification_calls = []

    class DeadlineNotifier:
        def notify_terminal_run(self, run_id, status, *, include_qc=False):
            notification_calls.append((run_id, status, include_qc))
            raise WorkerHardTimeout("deadline during failed-run terminal email")

    monkeypatch.setattr(
        "encode_pipeline.workers.runtime.compose_terminal_run_notifier",
        lambda **_kwargs: DeadlineNotifier(),
    )
    snakemake = configured.workspace_root.parent / "test-bin" / "snakemake"
    snakemake.write_text(
        "#!/bin/sh\nprintf 'before failure\\n'\nprintf 'safe error\\n' >&2\nexit 9\n",
        encoding="utf-8",
    )
    snakemake.chmod(0o755)
    connection = fakeredis.FakeRedis()
    run_queue = RqRunQueue(configured, connection=connection)
    run_queue.enqueue_execution(assignment)

    assert _run_burst(connection, run_queue) is True

    job = run_queue._queue.fetch_job(assignment.job_id)
    record, _events, persisted_assignment = _read_run(configured, "nonzero-run")
    assert job is not None
    assert job.is_failed
    failure = job.latest_result()
    assert failure is not None
    assert "WorkerExecutionError" in (failure.exc_string or "")
    assert "WorkerHardTimeout" not in (failure.exc_string or "")
    assert record.status is RunStatus.FAILED
    assert record.error is not None
    assert record.error.code == "RUN_EXECUTION_FAILED"
    assert record.error.context == {"reason_code": "LOCAL_RUN_EXECUTION_FAILED"}
    assert persisted_assignment is not None
    assert persisted_assignment.claimed_at is not None
    assert notification_calls == [
        (assignment.run_id, RunStatus.FAILED, False),
    ]


def test_rq_worker_sanitizes_unexpected_execution_exception(tmp_path, monkeypatch):
    configured = worker_settings(tmp_path)
    assignment = create_planned_run(
        configured,
        "unexpected-run",
        assign_queue=configured.queue_name,
    )
    assert assignment is not None
    _configure_worker_environment(monkeypatch, configured)

    def explode(_self, _run_id, _claim):
        raise RuntimeError("redis://private-password@internal:6379/0")

    monkeypatch.setattr(LocalExecutionService, "execute", explode)
    connection = fakeredis.FakeRedis()
    run_queue = RqRunQueue(configured, connection=connection)
    run_queue.enqueue_execution(assignment)

    assert _run_burst(connection, run_queue) is True

    job = run_queue._queue.fetch_job(assignment.job_id)
    record, _events, _assignment = _read_run(configured, "unexpected-run")
    assert job is not None
    assert job.is_failed
    failure = job.latest_result()
    assert failure is not None
    assert "private-password" not in (failure.exc_string or "")
    assert record.status is RunStatus.FAILED
    assert record.error is not None
    assert record.error.code == "RUN_WORKER_FAILED"
    assert "private-password" not in str(record.error.to_dict())


def test_rq_worker_persists_and_reraises_job_timeout(tmp_path, monkeypatch):
    configured = worker_settings(tmp_path)
    assignment = create_planned_run(
        configured,
        "timeout-run",
        assign_queue=configured.queue_name,
    )
    assert assignment is not None
    _configure_worker_environment(monkeypatch, configured)

    def timeout(_self, _run_id, _claim):
        raise WorkerHardTimeout("RQ deadline reached")

    monkeypatch.setattr(LocalExecutionService, "execute", timeout)
    connection = fakeredis.FakeRedis()
    run_queue = RqRunQueue(configured, connection=connection)
    run_queue.enqueue_execution(assignment)

    assert _run_burst(connection, run_queue) is True

    job = run_queue._queue.fetch_job(assignment.job_id)
    record, _events, _assignment = _read_run(configured, "timeout-run")
    assert job is not None
    assert job.is_failed
    failure = job.latest_result()
    assert failure is not None
    assert "WorkerHardTimeout" in (failure.exc_string or "")
    assert record.status is RunStatus.FAILED
    assert record.error is not None
    assert record.error.code == "RUN_WORKER_FAILED"
    assert record.error.context == {"reason_code": "WORKER_JOB_TIMEOUT"}


def test_rq_worker_reraises_timeout_after_atomic_artifact_success(
    tmp_path,
    monkeypatch,
):
    configured = worker_settings(tmp_path)
    assignment = create_planned_run(
        configured,
        "artifact-finalization-timeout",
        assign_queue=configured.queue_name,
    )
    assert assignment is not None
    _configure_worker_environment(monkeypatch, configured)
    monkeypatch.setattr(
        LocalExecutionService,
        "execute",
        _complete_workflow_successfully,
    )

    def commit_then_time_out(self, run_id, *, attempt_id):
        self._run_service.begin_artifact_result_attempt(
            run_id,
            attempt_id=attempt_id,
        )
        self._run_service.replace_artifacts(
            run_id,
            (),
            attempt_id=attempt_id,
        )
        raise WorkerHardTimeout("deadline after artifact finalization commit")

    monkeypatch.setattr(
        ArtifactExtractionService,
        "extract",
        commit_then_time_out,
    )
    connection = fakeredis.FakeRedis()
    run_queue = RqRunQueue(configured, connection=connection)
    run_queue.enqueue_execution(assignment)

    assert _run_burst(connection, run_queue) is True

    job = run_queue._queue.fetch_job(assignment.job_id)
    record, events, _assignment = _read_run(
        configured,
        "artifact-finalization-timeout",
    )
    result_state, artifacts, qc_metrics = _read_results(
        configured,
        "artifact-finalization-timeout",
    )
    assert job is not None
    assert job.is_failed
    failure = job.latest_result()
    assert failure is not None
    assert "WorkerHardTimeout" in (failure.exc_string or "")
    assert record.status is RunStatus.SUCCEEDED
    assert record.error is None
    assert result_state.artifact_attempt_status == "succeeded"
    assert result_state.artifact_outcome == "succeeded"
    assert result_state.qc_outcome is None
    assert artifacts == ()
    assert qc_metrics == ()
    outcome_types = [
        event.event_type
        for event in events
        if event.event_type
        in {
            "artifacts_indexed",
            "artifact_extraction_failed",
            "qc_metrics_indexed",
            "qc_metrics_indexing_failed",
        }
    ]
    assert outcome_types == ["artifacts_indexed"]
    assert "/private/" not in str([event.to_dict() for event in events])


def test_rq_worker_reraises_timeout_after_atomic_qc_success(
    tmp_path,
    monkeypatch,
):
    configured = worker_settings(tmp_path)
    assignment = create_planned_run(
        configured,
        "qc-finalization-timeout",
        assign_queue=configured.queue_name,
    )
    assert assignment is not None
    _configure_worker_environment(monkeypatch, configured)
    monkeypatch.setattr(
        LocalExecutionService,
        "execute",
        _complete_workflow_successfully,
    )

    def commit_artifacts(self, run_id, *, attempt_id):
        self._run_service.begin_artifact_result_attempt(
            run_id,
            attempt_id=attempt_id,
        )
        artifacts = self._run_service.replace_artifacts(
            run_id,
            (),
            attempt_id=attempt_id,
        )
        return Result.success(artifacts)

    def commit_then_time_out(
        self,
        run_id,
        artifacts,
        *,
        attempt_id,
        expected_artifact_generation,
    ):
        self._run_service.begin_qc_result_attempt(
            run_id,
            expected_artifact_generation=expected_artifact_generation,
            expected_artifacts=artifacts,
            attempt_id=attempt_id,
        )
        self._run_service.replace_qc_metrics(
            run_id,
            (),
            expected_artifacts=artifacts,
            attempt_id=attempt_id,
            expected_artifact_generation=expected_artifact_generation,
        )
        raise WorkerHardTimeout("deadline after QC finalization commit")

    monkeypatch.setattr(
        ArtifactExtractionService,
        "extract",
        commit_artifacts,
    )
    monkeypatch.setattr(
        QcSummaryIndexingService,
        "index",
        commit_then_time_out,
    )
    connection = fakeredis.FakeRedis()
    run_queue = RqRunQueue(configured, connection=connection)
    run_queue.enqueue_execution(assignment)

    assert _run_burst(connection, run_queue) is True

    job = run_queue._queue.fetch_job(assignment.job_id)
    record, events, _assignment = _read_run(
        configured,
        "qc-finalization-timeout",
    )
    result_state, artifacts, qc_metrics = _read_results(
        configured,
        "qc-finalization-timeout",
    )
    assert job is not None
    assert job.is_failed
    failure = job.latest_result()
    assert failure is not None
    assert "WorkerHardTimeout" in (failure.exc_string or "")
    assert record.status is RunStatus.SUCCEEDED
    assert record.error is None
    assert result_state.artifact_attempt_status == "succeeded"
    assert result_state.artifact_outcome == "succeeded"
    assert result_state.qc_attempt_status == "succeeded"
    assert result_state.qc_outcome == "succeeded"
    assert artifacts == ()
    assert qc_metrics == ()
    outcome_types = [
        event.event_type
        for event in events
        if event.event_type
        in {
            "artifacts_indexed",
            "artifact_extraction_failed",
            "qc_metrics_indexed",
            "qc_metrics_indexing_failed",
        }
    ]
    assert outcome_types == ["artifacts_indexed", "qc_metrics_indexed"]
    assert "/private/" not in str([event.to_dict() for event in events])


def test_rq_worker_generalizes_failure_when_durable_mapping_itself_fails(
    tmp_path,
    monkeypatch,
):
    configured = worker_settings(tmp_path)
    assignment = create_planned_run(
        configured,
        "mapping-failure-run",
        assign_queue=configured.queue_name,
    )
    assert assignment is not None
    _configure_worker_environment(monkeypatch, configured)

    def execution_failure(_self, _run_id, _claim):
        raise RuntimeError("redis://execution-secret@internal:6379/0")

    def mapping_failure(*_args, **_kwargs):
        raise RuntimeError("sqlite:///mapping-secret/private.db")

    monkeypatch.setattr(LocalExecutionService, "execute", execution_failure)
    monkeypatch.setattr(
        worker_jobs,
        "_record_unexpected_failure",
        mapping_failure,
    )
    connection = fakeredis.FakeRedis()
    run_queue = RqRunQueue(configured, connection=connection)
    run_queue.enqueue_execution(assignment)

    assert _run_burst(connection, run_queue) is True

    job = run_queue._queue.fetch_job(assignment.job_id)
    record, _events, _assignment = _read_run(configured, "mapping-failure-run")
    assert job is not None
    assert job.is_failed
    failure = job.latest_result()
    assert failure is not None
    assert "execution-secret" not in (failure.exc_string or "")
    assert "mapping-secret" not in (failure.exc_string or "")
    assert "WorkerExecutionError" in (failure.exc_string or "")
    assert record.status is RunStatus.QUEUED


def test_unexpected_failure_mapping_swallows_repository_errors():
    class BrokenRunService:
        def get_run(self, _run_id):
            raise RuntimeError("sqlite:///private-backend.db")

    worker_jobs._record_unexpected_failure(BrokenRunService(), "run-1")


def test_post_success_artifact_exception_does_not_fail_execution_job():
    recorded: list[tuple[str, str]] = []

    class Extraction:
        def extract(self, _run_id, *, attempt_id):
            validate_result_attempt_id(attempt_id)
            raise RuntimeError("/private/workspace")

        def record_unexpected_failure(self, run_id, *, attempt_id):
            recorded.append((run_id, attempt_id))

    runtime = _post_success_runtime(
        local_execution_service=SimpleNamespace(
            execute=lambda _run_id, _claim: Result.success(object())
        ),
        artifact_extraction_service=Extraction(),
    )

    worker_jobs._execute_claimed_run(runtime, "run-1", _acquired_claim())

    assert len(recorded) == 1
    assert recorded[0][0] == "run-1"
    validate_result_attempt_id(recorded[0][1])


def test_post_success_artifact_hard_timeout_is_recorded_and_reraised():
    recorded: list[tuple[str, str]] = []

    class Extraction:
        def extract(self, _run_id, *, attempt_id):
            validate_result_attempt_id(attempt_id)
            raise WorkerHardTimeout("artifact indexing exceeded the RQ deadline")

        def record_unexpected_failure(self, run_id, *, attempt_id):
            recorded.append((run_id, attempt_id))

    runtime = _post_success_runtime(
        local_execution_service=SimpleNamespace(
            execute=lambda _run_id, _claim: Result.success(object())
        ),
        artifact_extraction_service=Extraction(),
    )

    with pytest.raises(
        WorkerHardTimeout,
        match="artifact indexing exceeded the RQ deadline",
    ):
        worker_jobs._execute_claimed_run(runtime, "run-1", _acquired_claim())

    assert len(recorded) == 1
    assert recorded[0][0] == "run-1"
    validate_result_attempt_id(recorded[0][1])


@pytest.mark.parametrize(
    ("callback_error", "expected_message"),
    (
        (
            WorkerHardTimeout("deadline during artifact timeout evidence"),
            "deadline during artifact timeout evidence",
        ),
        (
            RuntimeError("/private/artifact-timeout-evidence"),
            "artifact indexing exceeded the RQ deadline",
        ),
    ),
)
def test_post_success_artifact_timeout_survives_failure_callback_fault(
    callback_error,
    expected_message,
):
    class Extraction:
        def extract(self, _run_id, *, attempt_id):
            validate_result_attempt_id(attempt_id)
            raise WorkerHardTimeout("artifact indexing exceeded the RQ deadline")

        def record_unexpected_failure(self, _run_id, *, attempt_id):
            validate_result_attempt_id(attempt_id)
            raise callback_error

    runtime = _post_success_runtime(
        local_execution_service=SimpleNamespace(
            execute=lambda _run_id, _claim: Result.success(object())
        ),
        artifact_extraction_service=Extraction(),
    )

    with pytest.raises(WorkerHardTimeout, match=expected_message) as raised:
        worker_jobs._execute_claimed_run(runtime, "run-1", _acquired_claim())

    formatted = "".join(traceback.format_exception(raised.value))
    assert "/private/artifact-timeout-evidence" not in formatted


def test_post_success_artifact_failure_callback_reraises_hard_timeout():
    class Extraction:
        def extract(self, _run_id, *, attempt_id):
            validate_result_attempt_id(attempt_id)
            raise RuntimeError("/private/workspace")

        def record_unexpected_failure(self, _run_id, *, attempt_id):
            validate_result_attempt_id(attempt_id)
            raise WorkerHardTimeout("deadline during artifact failure evidence")

    runtime = _post_success_runtime(
        local_execution_service=SimpleNamespace(
            execute=lambda _run_id, _claim: Result.success(object())
        ),
        artifact_extraction_service=Extraction(),
    )

    with pytest.raises(
        WorkerHardTimeout,
        match="deadline during artifact failure evidence",
    ) as raised:
        worker_jobs._execute_claimed_run(runtime, "run-1", _acquired_claim())

    formatted = "".join(traceback.format_exception(raised.value))
    assert "/private/workspace" not in formatted


def test_post_success_artifact_failure_callback_exception_remains_nonfatal():
    notification_calls = []

    class Extraction:
        def extract(self, _run_id, *, attempt_id):
            validate_result_attempt_id(attempt_id)
            raise RuntimeError("/private/workspace")

        def record_unexpected_failure(self, _run_id, *, attempt_id):
            validate_result_attempt_id(attempt_id)
            raise RuntimeError("/private/artifact-failure-evidence")

    runtime = _post_success_runtime(
        local_execution_service=SimpleNamespace(
            execute=lambda _run_id, _claim: Result.success(object())
        ),
        artifact_extraction_service=Extraction(),
        terminal_notifier=SimpleNamespace(
            notify_terminal_run=lambda *_args, **_kwargs: notification_calls.append(
                "notify"
            )
        ),
    )

    assert worker_jobs._execute_claimed_run(runtime, "run-1", _acquired_claim()) is None
    assert notification_calls == []


def test_post_success_result_state_hard_timeout_stops_before_qc():
    qc_calls = []

    runtime = _post_success_runtime(
        local_execution_service=SimpleNamespace(
            execute=lambda _run_id, _claim: Result.success(object())
        ),
        artifact_extraction_service=SimpleNamespace(
            extract=lambda _run_id, **_kwargs: Result.success(())
        ),
        run_service=SimpleNamespace(
            get_result_state=lambda _run_id: (_ for _ in ()).throw(
                WorkerHardTimeout("deadline while reading artifact generation")
            )
        ),
        qc_summary_indexing_service=SimpleNamespace(
            index=lambda *_args, **_kwargs: qc_calls.append("indexed")
        ),
    )

    with pytest.raises(
        WorkerHardTimeout,
        match="deadline while reading artifact generation",
    ):
        worker_jobs._execute_claimed_run(runtime, "run-1", _acquired_claim())

    assert qc_calls == []


def test_post_success_result_state_exception_stops_before_qc():
    qc_calls = []

    runtime = _post_success_runtime(
        local_execution_service=SimpleNamespace(
            execute=lambda _run_id, _claim: Result.success(object())
        ),
        artifact_extraction_service=SimpleNamespace(
            extract=lambda _run_id, **_kwargs: Result.success(())
        ),
        run_service=SimpleNamespace(
            get_result_state=lambda _run_id: (_ for _ in ()).throw(
                RuntimeError("/private/result-state")
            )
        ),
        qc_summary_indexing_service=SimpleNamespace(
            index=lambda *_args, **_kwargs: qc_calls.append("indexed")
        ),
    )

    assert worker_jobs._execute_claimed_run(runtime, "run-1", _acquired_claim()) is None
    assert qc_calls == []


def test_post_success_result_wrapper_hard_timeout_stops_before_qc():
    state = _ResultStateTracker()

    class TimedOutArtifactResult:
        @property
        def is_failure(self):
            raise WorkerHardTimeout("deadline while reading artifact result")

    runtime = _post_success_runtime(
        local_execution_service=SimpleNamespace(
            execute=lambda _run_id, _claim: Result.success(object())
        ),
        artifact_extraction_service=SimpleNamespace(
            extract=lambda _run_id, attempt_id: state.complete_artifact(
                attempt_id,
                TimedOutArtifactResult(),
            )
        ),
        run_service=state,
    )

    with pytest.raises(
        WorkerHardTimeout,
        match="deadline while reading artifact result",
    ):
        worker_jobs._execute_claimed_run(runtime, "run-1", _acquired_claim())


def test_post_success_qc_receives_only_successful_complete_artifact_set():
    artifacts = (object(), object())
    calls = []
    state = _ResultStateTracker()

    runtime = _post_success_runtime(
        local_execution_service=SimpleNamespace(
            execute=lambda _run_id, _claim: Result.success(object())
        ),
        artifact_extraction_service=SimpleNamespace(
            extract=lambda _run_id, attempt_id: state.complete_artifact(
                attempt_id,
                Result.success(artifacts),
            )
        ),
        run_service=state,
        qc_summary_indexing_service=SimpleNamespace(
            index=lambda run_id, values, **kwargs: (
                calls.append((run_id, values, kwargs)),
                state.complete_qc(
                    kwargs["attempt_id"],
                    kwargs["expected_artifact_generation"],
                ),
            )[1]
        ),
    )

    worker_jobs._execute_claimed_run(runtime, "run-1", _acquired_claim())

    assert len(calls) == 1
    assert calls[0][:2] == ("run-1", artifacts)
    assert calls[0][2]["expected_artifact_generation"] == ARTIFACT_GENERATION
    validate_result_attempt_id(calls[0][2]["attempt_id"])


def test_success_notification_runs_only_after_qc_finalizer_returns():
    calls = []
    state = _ResultStateTracker()

    class Notifier:
        def notify_terminal_run(self, run_id, status, *, include_qc=False):
            calls.append(("notify", run_id, status, include_qc))

    runtime = _post_success_runtime(
        local_execution_service=SimpleNamespace(
            execute=lambda _run_id, _claim: Result.success(object())
        ),
        artifact_extraction_service=SimpleNamespace(
            extract=lambda _run_id, attempt_id: (
                calls.append(("artifact",)),
                state.complete_artifact(attempt_id, Result.success(())),
            )[1]
        ),
        run_service=state,
        qc_summary_indexing_service=SimpleNamespace(
            index=lambda *_args, **kwargs: (
                calls.append(("qc",)),
                state.complete_qc(
                    kwargs["attempt_id"],
                    kwargs["expected_artifact_generation"],
                ),
            )[1]
        ),
        terminal_notifier=Notifier(),
    )

    worker_jobs._execute_claimed_run(runtime, "run-1", _acquired_claim())

    assert calls == [
        ("artifact",),
        ("qc",),
        ("notify", "run-1", RunStatus.SUCCEEDED, True),
    ]


def test_artifact_failure_without_durable_outcome_skips_success_notification():
    state = _ResultStateTracker()
    notifications = []

    def unresolved_failure(_run_id, attempt_id):
        state.artifact_attempt_id = attempt_id
        state.artifact_attempt_status = "pending"
        return Result.failure(
            [
                Issue(
                    code="ARTIFACT_INDEXING_FAILED",
                    message="Artifact indexing failed.",
                    source="test",
                )
            ]
        )

    runtime = _post_success_runtime(
        local_execution_service=SimpleNamespace(
            execute=lambda _run_id, _claim: Result.success(object())
        ),
        artifact_extraction_service=SimpleNamespace(extract=unresolved_failure),
        run_service=state,
        terminal_notifier=SimpleNamespace(
            notify_terminal_run=lambda *_args, **_kwargs: notifications.append("notify")
        ),
    )

    worker_jobs._execute_claimed_run(runtime, "run-1", _acquired_claim())

    assert notifications == []


def test_qc_failure_without_durable_outcome_skips_success_notification():
    state = _ResultStateTracker()
    notifications = []

    def unresolved_failure(_run_id, _artifacts, **kwargs):
        state.qc_attempt_id = kwargs["attempt_id"]
        state.qc_attempt_status = "pending"
        state.qc_attempt_artifact_generation = kwargs["expected_artifact_generation"]
        return Result.failure(
            [
                Issue(
                    code="QC_INDEXING_FAILED",
                    message="QC indexing failed.",
                    source="test",
                )
            ]
        )

    runtime = _post_success_runtime(
        local_execution_service=SimpleNamespace(
            execute=lambda _run_id, _claim: Result.success(object())
        ),
        artifact_extraction_service=SimpleNamespace(
            extract=lambda _run_id, attempt_id: state.complete_artifact(
                attempt_id,
                Result.success(()),
            )
        ),
        run_service=state,
        qc_summary_indexing_service=SimpleNamespace(index=unresolved_failure),
        terminal_notifier=SimpleNamespace(
            notify_terminal_run=lambda *_args, **_kwargs: notifications.append("notify")
        ),
    )

    worker_jobs._execute_claimed_run(runtime, "run-1", _acquired_claim())

    assert notifications == []


def test_losing_success_cas_skips_result_indexing_and_notification():
    calls = []
    runtime = _post_success_runtime(
        local_execution_service=SimpleNamespace(
            execute=lambda _run_id, _claim: Result.success(
                SimpleNamespace(terminal_transition_won=False)
            )
        ),
        artifact_extraction_service=SimpleNamespace(
            extract=lambda *_args, **_kwargs: calls.append("artifact")
        ),
        terminal_notifier=SimpleNamespace(
            notify_terminal_run=lambda *_args, **_kwargs: calls.append("notify")
        ),
    )

    worker_jobs._execute_claimed_run(runtime, "run-1", _acquired_claim())

    assert calls == []


def test_success_notification_hard_timeout_does_not_change_worker_result():
    calls = []
    state = _ResultStateTracker()

    class Notifier:
        def notify_terminal_run(self, _run_id, _status, *, include_qc=False):
            assert include_qc is True
            calls.append("notify")
            raise WorkerHardTimeout("deadline during terminal email")

    runtime = _post_success_runtime(
        local_execution_service=SimpleNamespace(
            execute=lambda _run_id, _claim: Result.success(object())
        ),
        artifact_extraction_service=SimpleNamespace(
            extract=lambda _run_id, attempt_id: (
                calls.append("artifact"),
                state.complete_artifact(
                    attempt_id,
                    Result.failure(
                        [
                            Issue(
                                code="ARTIFACT_INDEXING_FAILED",
                                message="Artifact indexing failed.",
                                source="test",
                            )
                        ]
                    ),
                    status="failed",
                ),
            )[1]
        ),
        run_service=state,
        terminal_notifier=Notifier(),
    )

    assert worker_jobs._execute_claimed_run(runtime, "run-1", _acquired_claim()) is None
    assert calls == ["artifact", "notify"]


def test_post_success_qc_exception_is_recorded_without_failing_execution():
    recorded = []
    state = _ResultStateTracker()

    class QcIndexing:
        def index(self, _run_id, _artifacts, **_kwargs):
            raise RuntimeError("/private/qc/workspace")

        def record_unexpected_failure(self, run_id, **kwargs):
            recorded.append((run_id, kwargs))
            state.complete_qc(
                kwargs["attempt_id"],
                kwargs["expected_artifact_generation"],
                status="failed",
            )

    runtime = _post_success_runtime(
        local_execution_service=SimpleNamespace(
            execute=lambda _run_id, _claim: Result.success(object())
        ),
        artifact_extraction_service=SimpleNamespace(
            extract=lambda _run_id, attempt_id: state.complete_artifact(
                attempt_id,
                Result.success(()),
            )
        ),
        run_service=state,
        qc_summary_indexing_service=QcIndexing(),
    )

    worker_jobs._execute_claimed_run(runtime, "run-1", _acquired_claim())

    assert len(recorded) == 1
    assert recorded[0][0] == "run-1"
    assert recorded[0][1]["expected_artifact_generation"] == ARTIFACT_GENERATION
    validate_result_attempt_id(recorded[0][1]["attempt_id"])


def test_post_success_qc_hard_timeout_is_recorded_and_reraised():
    recorded = []
    state = _ResultStateTracker()

    class QcIndexing:
        def index(self, _run_id, _artifacts, **_kwargs):
            raise WorkerHardTimeout("QC indexing exceeded the RQ deadline")

        def record_unexpected_failure(self, run_id, **kwargs):
            recorded.append((run_id, kwargs))

    runtime = _post_success_runtime(
        local_execution_service=SimpleNamespace(
            execute=lambda _run_id, _claim: Result.success(object())
        ),
        artifact_extraction_service=SimpleNamespace(
            extract=lambda _run_id, attempt_id: state.complete_artifact(
                attempt_id,
                Result.success(()),
            )
        ),
        run_service=state,
        qc_summary_indexing_service=QcIndexing(),
    )

    with pytest.raises(
        WorkerHardTimeout,
        match="QC indexing exceeded the RQ deadline",
    ):
        worker_jobs._execute_claimed_run(runtime, "run-1", _acquired_claim())

    assert len(recorded) == 1
    assert recorded[0][0] == "run-1"
    assert recorded[0][1]["expected_artifact_generation"] == ARTIFACT_GENERATION
    validate_result_attempt_id(recorded[0][1]["attempt_id"])


@pytest.mark.parametrize(
    ("callback_error", "expected_message"),
    (
        (
            WorkerHardTimeout("deadline during QC timeout evidence"),
            "deadline during QC timeout evidence",
        ),
        (
            RuntimeError("/private/qc-timeout-evidence"),
            "QC indexing exceeded the RQ deadline",
        ),
    ),
)
def test_post_success_qc_timeout_survives_failure_callback_fault(
    callback_error,
    expected_message,
):
    state = _ResultStateTracker()

    class QcIndexing:
        def index(self, _run_id, _artifacts, **_kwargs):
            raise WorkerHardTimeout("QC indexing exceeded the RQ deadline")

        def record_unexpected_failure(self, _run_id, **kwargs):
            validate_result_attempt_id(kwargs["attempt_id"])
            raise callback_error

    runtime = _post_success_runtime(
        local_execution_service=SimpleNamespace(
            execute=lambda _run_id, _claim: Result.success(object())
        ),
        artifact_extraction_service=SimpleNamespace(
            extract=lambda _run_id, attempt_id: state.complete_artifact(
                attempt_id,
                Result.success(()),
            )
        ),
        run_service=state,
        qc_summary_indexing_service=QcIndexing(),
    )

    with pytest.raises(WorkerHardTimeout, match=expected_message) as raised:
        worker_jobs._execute_claimed_run(runtime, "run-1", _acquired_claim())

    formatted = "".join(traceback.format_exception(raised.value))
    assert "/private/qc-timeout-evidence" not in formatted


def test_post_success_qc_failure_callback_reraises_hard_timeout():
    state = _ResultStateTracker()

    class QcIndexing:
        def index(self, _run_id, _artifacts, **_kwargs):
            raise RuntimeError("/private/qc/workspace")

        def record_unexpected_failure(self, _run_id, **kwargs):
            validate_result_attempt_id(kwargs["attempt_id"])
            raise WorkerHardTimeout("deadline during QC failure evidence")

    runtime = _post_success_runtime(
        local_execution_service=SimpleNamespace(
            execute=lambda _run_id, _claim: Result.success(object())
        ),
        artifact_extraction_service=SimpleNamespace(
            extract=lambda _run_id, attempt_id: state.complete_artifact(
                attempt_id,
                Result.success(()),
            )
        ),
        run_service=state,
        qc_summary_indexing_service=QcIndexing(),
    )

    with pytest.raises(
        WorkerHardTimeout,
        match="deadline during QC failure evidence",
    ) as raised:
        worker_jobs._execute_claimed_run(runtime, "run-1", _acquired_claim())

    formatted = "".join(traceback.format_exception(raised.value))
    assert "/private/qc/workspace" not in formatted


def test_post_success_qc_failure_callback_exception_remains_nonfatal():
    notification_calls = []
    state = _ResultStateTracker()

    class QcIndexing:
        def index(self, _run_id, _artifacts, **_kwargs):
            raise RuntimeError("/private/qc/workspace")

        def record_unexpected_failure(self, _run_id, **kwargs):
            validate_result_attempt_id(kwargs["attempt_id"])
            raise RuntimeError("/private/qc-failure-evidence")

    runtime = _post_success_runtime(
        local_execution_service=SimpleNamespace(
            execute=lambda _run_id, _claim: Result.success(object())
        ),
        artifact_extraction_service=SimpleNamespace(
            extract=lambda _run_id, attempt_id: state.complete_artifact(
                attempt_id,
                Result.success(()),
            )
        ),
        run_service=state,
        qc_summary_indexing_service=QcIndexing(),
        terminal_notifier=SimpleNamespace(
            notify_terminal_run=lambda *_args, **_kwargs: notification_calls.append(
                "notify"
            )
        ),
    )

    assert worker_jobs._execute_claimed_run(runtime, "run-1", _acquired_claim()) is None
    assert notification_calls == []


def test_failure_mapping_does_not_swallow_rq_timeout():
    class TimedOutRunService:
        def get_run(self, _run_id):
            raise WorkerHardTimeout("RQ deadline reached during failure mapping")

    with pytest.raises(WorkerHardTimeout, match="during failure mapping"):
        worker_jobs._record_unexpected_failure_safely(
            TimedOutRunService(),
            "run-1",
        )


def test_worker_job_preserves_migration_admission_reason_without_db_fallback(
    monkeypatch,
):
    monkeypatch.setattr(
        worker_jobs,
        "get_current_job",
        lambda: SimpleNamespace(id="job-1", origin="queue-1"),
    )
    monkeypatch.setattr(
        worker_jobs,
        "open_worker_runtime",
        lambda: (_ for _ in ()).throw(
            worker_jobs.MigrationAdmissionError("MIGRATION_REVISION_DIGEST_MISMATCH")
        ),
    )
    monkeypatch.setattr(
        worker_jobs,
        "_record_initialization_failure_fallback",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            AssertionError("migration rejection must not reopen SQLite")
        ),
    )

    with pytest.raises(worker_jobs.WorkerExecutionError) as raised:
        worker_jobs.run_execution_job("migration-rejected-run")

    assert raised.value.reason_code == "MIGRATION_REVISION_DIGEST_MISMATCH"
    assert str(raised.value) == (
        "migration execution admission failed [MIGRATION_REVISION_DIGEST_MISMATCH]"
    )
    assert raised.value.__cause__ is None


def test_worker_composition_failure_is_public_safe_and_durably_failed(
    tmp_path,
    monkeypatch,
):
    configured = worker_settings(tmp_path)
    assignment = create_planned_run(
        configured,
        "composition-failure-run",
        assign_queue=configured.queue_name,
    )
    assert assignment is not None
    prepared = _prepare_requeue(configured, assignment)
    _configure_worker_environment(monkeypatch, configured)
    monkeypatch.setattr(
        worker_jobs,
        "get_current_job",
        lambda: SimpleNamespace(
            id=prepared.job_id,
            origin=configured.queue_name,
        ),
    )

    def fail_composition():
        raise RuntimeError("sqlite:///private/path/platform.db")

    monkeypatch.setattr(worker_jobs, "open_worker_runtime", fail_composition)

    with pytest.raises(worker_jobs.WorkerExecutionError) as raised:
        worker_jobs.run_execution_job("composition-failure-run")

    assert "private/path" not in str(raised.value)
    assert raised.value.__cause__ is None
    record, events, persisted_assignment = _read_run(
        configured,
        "composition-failure-run",
    )
    assert record.status is RunStatus.FAILED
    assert record.error is not None
    assert record.error.context == {"reason_code": "WORKER_INITIALIZATION_FAILED"}
    assert persisted_assignment is not None
    assert persisted_assignment.dispatched_at is not None
    assert persisted_assignment.requeue_requested_at == prepared.requeue_requested_at
    assert persisted_assignment.requeue_confirmed_at is not None
    assert "run_requeue_delivery_observed_by_worker" in {
        event.event_type for event in events
    }


def test_invalid_enabled_notification_config_still_durably_fails_run(
    tmp_path,
    monkeypatch,
):
    configured = worker_settings(tmp_path)
    assignment = create_planned_run(
        configured,
        "invalid-notification-config-run",
        assign_queue=configured.queue_name,
    )
    assert assignment is not None
    prepared = _prepare_requeue(configured, assignment)
    _configure_worker_environment(monkeypatch, configured)
    monkeypatch.setenv("HELIXWEAVE_TERMINAL_EMAIL_ENABLED", "true")
    monkeypatch.setattr(
        worker_jobs,
        "get_current_job",
        lambda: SimpleNamespace(
            id=prepared.job_id,
            origin=configured.queue_name,
        ),
    )

    with pytest.raises(worker_jobs.WorkerExecutionError) as raised:
        worker_jobs.run_execution_job(assignment.run_id)

    assert "could not be initialized" in str(raised.value)
    record, _events, persisted_assignment = _read_run(
        configured,
        assignment.run_id,
    )
    assert record.status is RunStatus.FAILED
    assert record.error is not None
    assert record.error.context == {"reason_code": "WORKER_INITIALIZATION_FAILED"}
    assert persisted_assignment is not None
    assert persisted_assignment.requeue_confirmed_at is not None


@pytest.mark.parametrize("drift_field", ("job_id", "queue_name"))
def test_worker_composition_fallback_refuses_identity_drift_without_mutation(
    tmp_path,
    monkeypatch,
    drift_field,
):
    configured = worker_settings(tmp_path)
    assignment = create_planned_run(
        configured,
        f"identity-drift-{drift_field}",
        assign_queue=configured.queue_name,
    )
    assert assignment is not None
    prepared = _prepare_requeue(configured, assignment)
    _configure_worker_environment(monkeypatch, configured)
    current_job = SimpleNamespace(
        id="wrong-job" if drift_field == "job_id" else prepared.job_id,
        origin="wrong-queue" if drift_field == "queue_name" else configured.queue_name,
    )
    monkeypatch.setattr(worker_jobs, "get_current_job", lambda: current_job)
    monkeypatch.setattr(
        worker_jobs,
        "open_worker_runtime",
        lambda **_kwargs: (_ for _ in ()).throw(RuntimeError("composition failed")),
    )

    with pytest.raises(worker_jobs.WorkerExecutionError):
        worker_jobs.run_execution_job(assignment.run_id)

    record, events, persisted_assignment = _read_run(configured, assignment.run_id)
    assert record.status is RunStatus.QUEUED
    assert persisted_assignment == prepared
    assert persisted_assignment.requeue_confirmed_at is None
    assert all(
        event.event_type != "run_requeue_delivery_observed_by_worker"
        for event in events
    )
    assert all(event.status is not RunStatus.FAILED for event in events)


def test_worker_hard_timeout_during_composition_uses_durable_fallback(
    tmp_path,
    monkeypatch,
):
    configured = worker_settings(tmp_path)
    assignment = create_planned_run(
        configured,
        "composition-timeout-run",
        assign_queue=configured.queue_name,
    )
    assert assignment is not None
    _configure_worker_environment(monkeypatch, configured)
    monkeypatch.setattr(
        worker_jobs,
        "get_current_job",
        lambda: SimpleNamespace(
            id=assignment.job_id,
            origin=configured.queue_name,
        ),
    )
    monkeypatch.setattr(
        worker_jobs,
        "open_worker_runtime",
        lambda **_kwargs: (_ for _ in ()).throw(
            WorkerHardTimeout("RQ deadline reached")
        ),
    )

    with pytest.raises(WorkerHardTimeout, match="RQ deadline reached"):
        worker_jobs.run_execution_job(assignment.run_id)

    record, _events, persisted_assignment = _read_run(configured, assignment.run_id)
    assert record.status is RunStatus.FAILED
    assert record.error is not None
    assert record.error.context == {"reason_code": "WORKER_JOB_TIMEOUT"}
    assert persisted_assignment is not None
    assert persisted_assignment.dispatched_at is not None


def test_worker_claim_timeout_survives_timeout_in_durable_fallback(monkeypatch):
    run_id = "claim-timeout"
    runtime = SimpleNamespace(run_service=object())

    class RuntimeContext:
        def __enter__(self):
            return runtime

        def __exit__(self, *_args):
            return None

    monkeypatch.setattr(
        worker_jobs,
        "get_current_job",
        lambda: SimpleNamespace(id="job-1", origin="queue"),
    )
    monkeypatch.setattr(worker_jobs, "open_worker_runtime", RuntimeContext)
    monkeypatch.setattr(
        worker_jobs,
        "_initialize_execution_with_runtime",
        lambda *_args: (_ for _ in ()).throw(
            WorkerHardTimeout("original claim deadline")
        ),
    )
    monkeypatch.setattr(
        worker_jobs,
        "_cleanup_runtime_managed_containers",
        lambda *_args: True,
    )
    monkeypatch.setattr(
        worker_jobs,
        "_record_unexpected_failure_safely",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        worker_jobs,
        "_record_initialization_failure_fallback",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(
            WorkerHardTimeout("secondary fallback deadline")
        ),
    )

    with pytest.raises(WorkerHardTimeout, match="original claim deadline") as raised:
        worker_jobs.run_execution_job(run_id)

    formatted = "".join(traceback.format_exception(raised.value))
    assert "secondary fallback deadline" not in formatted


@pytest.mark.parametrize("cleanup_outcome", ("failure", "exception"))
def test_worker_hard_timeout_survives_cleanup_failure(
    monkeypatch,
    cleanup_outcome,
):
    run_id = f"cleanup-{cleanup_outcome}"
    runtime = SimpleNamespace(run_service=object())
    cleanup_events = []

    class RuntimeContext:
        def __enter__(self):
            return runtime

        def __exit__(self, *_args):
            return None

    monkeypatch.setattr(
        worker_jobs,
        "get_current_job",
        lambda: SimpleNamespace(id="job-1", origin="queue"),
    )
    monkeypatch.setattr(worker_jobs, "open_worker_runtime", RuntimeContext)
    monkeypatch.setattr(
        worker_jobs,
        "_initialize_execution_with_runtime",
        lambda *_args: True,
    )
    monkeypatch.setattr(
        worker_jobs,
        "_execute_claimed_run",
        lambda *_args: (_ for _ in ()).throw(WorkerHardTimeout("original RQ deadline")),
    )

    def cleanup(*_args):
        if cleanup_outcome == "exception":
            raise RuntimeError("/private/cleanup")
        return False

    monkeypatch.setattr(
        worker_jobs,
        "_cleanup_runtime_managed_containers",
        cleanup,
    )
    monkeypatch.setattr(
        worker_jobs,
        "_record_cleanup_failure_safely",
        lambda service, observed_run_id: cleanup_events.append(
            (service, observed_run_id)
        ),
    )
    monkeypatch.setattr(
        worker_jobs,
        "_record_unexpected_failure_safely",
        lambda *_args, **_kwargs: None,
    )
    monkeypatch.setattr(
        worker_jobs,
        "_record_initialization_failure_fallback",
        lambda *_args, **_kwargs: None,
    )

    with pytest.raises(WorkerHardTimeout, match="original RQ deadline") as raised:
        worker_jobs.run_execution_job(run_id)

    formatted = "".join(traceback.format_exception(raised.value))
    assert "/private/cleanup" not in formatted
    if cleanup_outcome == "exception":
        assert cleanup_events == [(runtime.run_service, run_id)]


@pytest.mark.parametrize(
    "secondary_fault",
    ("cleanup_timeout", "cleanup_evidence_timeout", "failure_evidence_timeout"),
)
def test_worker_hard_timeout_survives_secondary_timeout(
    monkeypatch,
    secondary_fault,
):
    run_id = f"secondary-{secondary_fault}"
    runtime = SimpleNamespace(run_service=object())

    class RuntimeContext:
        def __enter__(self):
            return runtime

        def __exit__(self, *_args):
            return None

    monkeypatch.setattr(
        worker_jobs,
        "get_current_job",
        lambda: SimpleNamespace(id="job-1", origin="queue"),
    )
    monkeypatch.setattr(worker_jobs, "open_worker_runtime", RuntimeContext)
    monkeypatch.setattr(
        worker_jobs,
        "_initialize_execution_with_runtime",
        lambda *_args: True,
    )
    monkeypatch.setattr(
        worker_jobs,
        "_execute_claimed_run",
        lambda *_args: (_ for _ in ()).throw(WorkerHardTimeout("original RQ deadline")),
    )

    def cleanup(*_args):
        if secondary_fault == "cleanup_timeout":
            raise WorkerHardTimeout("secondary cleanup deadline")
        if secondary_fault == "cleanup_evidence_timeout":
            raise RuntimeError("/private/cleanup")
        return True

    def record_cleanup_failure(*_args):
        if secondary_fault == "cleanup_evidence_timeout":
            raise WorkerHardTimeout("secondary cleanup evidence deadline")
        raise AssertionError("cleanup evidence callback was not expected")

    def record_failure(*_args, **_kwargs):
        if secondary_fault == "failure_evidence_timeout":
            raise WorkerHardTimeout("secondary failure evidence deadline")

    monkeypatch.setattr(
        worker_jobs,
        "_cleanup_runtime_managed_containers",
        cleanup,
    )
    monkeypatch.setattr(
        worker_jobs,
        "_record_cleanup_failure_safely",
        record_cleanup_failure,
    )
    monkeypatch.setattr(
        worker_jobs,
        "_record_unexpected_failure_safely",
        record_failure,
    )
    monkeypatch.setattr(
        worker_jobs,
        "_record_initialization_failure_fallback",
        lambda *_args, **_kwargs: None,
    )

    with pytest.raises(WorkerHardTimeout, match="original RQ deadline") as raised:
        worker_jobs.run_execution_job(run_id)

    formatted = "".join(traceback.format_exception(raised.value))
    assert "/private/cleanup" not in formatted
    assert "secondary cleanup deadline" not in formatted
    assert "secondary cleanup evidence deadline" not in formatted
    assert "secondary failure evidence deadline" not in formatted


def test_missing_worker_adapter_is_initialization_failure_not_identity_error(
    tmp_path,
    monkeypatch,
):
    configured = worker_settings(tmp_path)
    assignment = create_planned_run(
        configured,
        "missing-adapter-run",
        assign_queue=configured.queue_name,
    )
    assert assignment is not None
    _configure_worker_environment(monkeypatch, configured)
    monkeypatch.setattr(
        worker_jobs,
        "get_current_job",
        lambda: SimpleNamespace(
            id=assignment.job_id,
            origin=configured.queue_name,
        ),
    )
    monkeypatch.setattr(
        "encode_pipeline.services.defaults.create_default_workflow_registry",
        lambda **_kwargs: WorkflowRegistry(),
    )

    with pytest.raises(worker_jobs.WorkerExecutionError):
        worker_jobs.run_execution_job(assignment.run_id)

    record, _events, _assignment = _read_run(configured, assignment.run_id)
    assert record.status is RunStatus.FAILED
    assert record.error is not None
    assert record.error.context == {"reason_code": "WORKER_INITIALIZATION_FAILED"}


def test_work_horse_death_is_mapped_to_durable_failure(tmp_path, monkeypatch):
    configured = worker_settings(tmp_path)
    assignment = create_planned_run(
        configured,
        "killed-run",
        assign_queue=configured.queue_name,
    )
    assert assignment is not None
    _configure_worker_environment(monkeypatch, configured)

    class KilledJob:
        id = assignment.job_id
        func_name = "encode_pipeline.workers.jobs.run_execution_job"
        args = (assignment.run_id,)
        kwargs = {}
        origin = configured.queue_name

    handle_work_horse_killed(KilledJob(), 123, 9, object())

    record, _events, persisted_assignment = _read_run(configured, "killed-run")
    assert record.status is RunStatus.FAILED
    assert record.error is not None
    assert record.error.code == "RUN_WORKER_FAILED"
    assert record.error.context == {"reason_code": "WORKER_PROCESS_TERMINATED"}
    assert persisted_assignment is not None
    assert persisted_assignment.dispatched_at is not None


def test_work_horse_death_confirms_pending_requeue_before_terminal_failure(
    tmp_path,
    monkeypatch,
):
    configured = worker_settings(tmp_path)
    assignment = create_planned_run(
        configured,
        "killed-pending-requeue",
        assign_queue=configured.queue_name,
    )
    assert assignment is not None
    prepared = _prepare_requeue(configured, assignment)
    _configure_worker_environment(monkeypatch, configured)

    handle_work_horse_killed(_execution_job(prepared), 123, 9, object())

    record, events, persisted = _read_run(configured, prepared.run_id)
    assert record.status is RunStatus.FAILED
    assert persisted is not None
    assert persisted.claimed_at is None
    assert persisted.requeue_confirmed_at is not None
    assert "run_requeue_delivery_observed_by_worker" in {
        event.event_type for event in events
    }


def test_work_horse_cleanup_failure_does_not_fabricate_terminal_failure(
    tmp_path,
    monkeypatch,
):
    configured = worker_settings(tmp_path)
    assignment = create_planned_run(
        configured,
        "killed-cleanup-failure",
        assign_queue=configured.queue_name,
    )
    assert assignment is not None
    _configure_worker_environment(monkeypatch, configured)
    calls = []

    def cleanup(_runtime, run_id):
        calls.append(run_id)
        return False

    monkeypatch.setattr(worker_jobs, "_cleanup_runtime_managed_containers", cleanup)

    handle_work_horse_killed(
        _execution_job(assignment),
        123,
        9,
        object(),
    )

    record, _events, persisted = _read_run(configured, assignment.run_id)
    assert calls == [assignment.run_id]
    assert record.status is RunStatus.PLANNED
    assert record.error is None
    assert persisted == assignment


def test_stopped_callback_acknowledges_user_intent_after_horse_exit(
    tmp_path,
    monkeypatch,
):
    configured = worker_settings(tmp_path)
    assignment = create_planned_run(
        configured,
        "cancelled-running-run",
        assign_queue=configured.queue_name,
    )
    assert assignment is not None
    _configure_worker_environment(monkeypatch, configured)
    notification_calls = []

    class Notifier:
        def notify_terminal_run(self, run_id, status, *, include_qc=False):
            notification_calls.append((run_id, status, include_qc))

    monkeypatch.setattr(
        worker_jobs,
        "compose_terminal_run_notifier",
        lambda **_kwargs: Notifier(),
    )
    claimed = _prepare_running_assignment(configured, assignment, request_cancel=True)
    job = _execution_job(claimed)

    handle_execution_stopped(job, fakeredis.FakeRedis())

    record, events, persisted = _read_run(configured, assignment.run_id)
    assert record.status is RunStatus.CANCELLED
    assert record.ended_at is not None
    assert record.cancellation_reason == "User requested cancellation."
    assert persisted is not None
    assert persisted.cancellation_acknowledged_at is not None
    assert [event.event_type for event in events].count(
        "cancellation_acknowledged"
    ) == 1

    handle_execution_stopped(job, fakeredis.FakeRedis())
    _record, repeated_events, repeated_assignment = _read_run(
        configured,
        assignment.run_id,
    )
    assert repeated_events == events
    assert repeated_assignment == persisted
    assert notification_calls == [
        (assignment.run_id, RunStatus.CANCELLED, False),
    ]


def test_stopped_cleanup_failure_leaves_cancellation_unacknowledged(
    tmp_path,
    monkeypatch,
):
    configured = worker_settings(tmp_path)
    assignment = create_planned_run(
        configured,
        "cancel-cleanup-failure",
        assign_queue=configured.queue_name,
    )
    assert assignment is not None
    _configure_worker_environment(monkeypatch, configured)
    claimed = _prepare_running_assignment(configured, assignment, request_cancel=True)
    calls = []

    def cleanup(_settings, run_id, _run_service):
        calls.append(run_id)
        return False

    monkeypatch.setattr(worker_jobs, "_cleanup_settings_managed_containers", cleanup)

    handle_execution_stopped(_execution_job(claimed), fakeredis.FakeRedis())

    record, _events, persisted = _read_run(configured, assignment.run_id)
    assert calls == [assignment.run_id]
    assert record.status is RunStatus.RUNNING
    assert record.ended_at is None
    assert persisted is not None
    assert persisted.cancellation_requested_at is not None
    assert persisted.cancellation_acknowledged_at is None


@pytest.mark.parametrize("endpoint_matches", (True, False))
def test_stopped_callback_uses_only_the_durable_cleanup_binding(
    tmp_path,
    monkeypatch,
    endpoint_matches,
):
    configured = worker_settings(tmp_path)
    assignment = create_planned_run(
        configured,
        f"cancel-bound-cleanup-{endpoint_matches}",
        assign_queue=configured.queue_name,
    )
    assert assignment is not None
    _configure_worker_environment(monkeypatch, configured)
    executable = (tmp_path / "admin-bin" / "docker").resolve()
    socket = (tmp_path / "admin-docker.sock").resolve()
    current_endpoint = managed_container_endpoint_identity(executable, socket)
    durable_endpoint = current_endpoint if endpoint_matches else "f" * 64
    durable_scope = "e" * 64
    claimed = _prepare_running_assignment(
        configured,
        assignment,
        request_cancel=True,
        managed_container_scope=durable_scope,
        managed_container_endpoint_identity=durable_endpoint,
    )
    monkeypatch.setenv(MANAGED_DOCKER_EXECUTABLE_ENV, str(executable))
    monkeypatch.setenv(MANAGED_DOCKER_SOCKET_ENV, str(socket))
    cleanup_scopes = []

    class Cleaner:
        def __init__(self, *, executable, unix_socket):
            self.endpoint_identity = managed_container_endpoint_identity(
                executable,
                unix_socket,
            )

        def cleanup(self, scope):
            cleanup_scopes.append(scope)
            return Result.success(None)

    monkeypatch.setattr(worker_jobs, "ManagedContainerCleaner", Cleaner)
    before_record, before_events, before_assignment = _read_run(
        configured,
        assignment.run_id,
    )

    handle_execution_stopped(_execution_job(claimed), fakeredis.FakeRedis())

    record, events, persisted = _read_run(configured, assignment.run_id)
    if endpoint_matches:
        assert cleanup_scopes == [durable_scope]
        assert record.status is RunStatus.CANCELLED
        assert persisted is not None
        assert persisted.cancellation_acknowledged_at is not None
    else:
        assert cleanup_scopes == []
        assert record == before_record
        assert persisted == before_assignment
        assert [event.event_type for event in events[:-1]] == [
            event.event_type for event in before_events
        ]
        assert events[-1].event_type == "execution_cleanup_failed"
        assert all(event.event_type != "cancellation_acknowledged" for event in events)


def test_stopped_callback_without_user_intent_fails_truthfully(
    tmp_path,
    monkeypatch,
):
    configured = worker_settings(tmp_path)
    assignment = create_planned_run(
        configured,
        "unexpected-stopped-run",
        assign_queue=configured.queue_name,
    )
    assert assignment is not None
    _configure_worker_environment(monkeypatch, configured)
    claimed = _prepare_running_assignment(configured, assignment)

    handle_execution_stopped(_execution_job(claimed), fakeredis.FakeRedis())

    record, events, persisted = _read_run(configured, assignment.run_id)
    assert record.status is RunStatus.FAILED
    assert record.cancellation_reason is None
    assert record.error is not None
    assert record.error.code == "RUN_WORKER_STOPPED_UNEXPECTEDLY"
    assert record.error.context == {"reason_code": "WORKER_STOP_WITHOUT_CANCELLATION"}
    assert persisted is not None
    assert persisted.cancellation_acknowledged_at is None
    assert events[-1].event_type == "execution_stopped_unexpectedly"


@pytest.mark.parametrize("claimed", (False, True))
def test_stopped_callback_before_running_fails_truthfully(
    tmp_path,
    monkeypatch,
    claimed,
):
    configured = worker_settings(tmp_path)
    assignment = create_planned_run(
        configured,
        f"unexpected-queued-stop-{claimed}",
        assign_queue=configured.queue_name,
    )
    assert assignment is not None
    _configure_worker_environment(monkeypatch, configured)
    persistence = open_run_persistence(configured.database_url)
    try:
        service = create_default_run_service(
            create_default_workflow_registry(),
            repository=persistence.repository,
        )
        service.mark_execution_dispatched(
            assignment.run_id,
            job_id=assignment.job_id,
        )
        service.queue_dispatched_run(
            assignment.run_id,
            job_id=assignment.job_id,
            backend=assignment.backend,
            queue_name=assignment.queue_name,
        )
        if claimed:
            service.claim_execution_assignment(
                assignment.run_id,
                job_id=assignment.job_id,
                backend=assignment.backend,
                queue_name=assignment.queue_name,
            )
        queued_assignment = service.get_execution_assignment(assignment.run_id)
        assert queued_assignment is not None
    finally:
        persistence.close()

    handle_execution_stopped(
        _execution_job(queued_assignment),
        fakeredis.FakeRedis(),
    )

    record, events, persisted = _read_run(configured, assignment.run_id)
    assert record.status is RunStatus.FAILED
    assert record.error is not None
    assert record.error.code == "RUN_WORKER_STOPPED_UNEXPECTEDLY"
    assert persisted is not None
    assert persisted.cancellation_acknowledged_at is None
    assert events[-1].event_type == "execution_stopped_unexpectedly"
    assert events[-1].context["previous_status"] == "queued"


def test_stopped_callback_confirms_pending_requeue_before_terminal_failure(
    tmp_path,
    monkeypatch,
):
    configured = worker_settings(tmp_path)
    assignment = create_planned_run(
        configured,
        "stopped-pending-requeue",
        assign_queue=configured.queue_name,
    )
    assert assignment is not None
    prepared = _prepare_requeue(configured, assignment)
    _configure_worker_environment(monkeypatch, configured)

    handle_execution_stopped(_execution_job(prepared), fakeredis.FakeRedis())

    record, events, persisted = _read_run(configured, prepared.run_id)
    assert record.status is RunStatus.FAILED
    assert persisted is not None
    assert persisted.claimed_at is None
    assert persisted.requeue_confirmed_at is not None
    assert "run_requeue_delivery_observed_by_worker" in {
        event.event_type for event in events
    }


def test_stopped_callback_before_dispatch_fails_planned_truthfully(
    tmp_path,
    monkeypatch,
):
    configured = worker_settings(tmp_path)
    assignment = create_planned_run(
        configured,
        "unexpected-planned-stop",
        assign_queue=configured.queue_name,
    )
    assert assignment is not None
    _configure_worker_environment(monkeypatch, configured)

    handle_execution_stopped(
        _execution_job(assignment),
        fakeredis.FakeRedis(),
    )

    record, events, persisted = _read_run(configured, assignment.run_id)
    assert record.status is RunStatus.FAILED
    assert record.error is not None
    assert record.error.code == "RUN_WORKER_STOPPED_UNEXPECTEDLY"
    assert persisted is not None
    assert persisted.dispatched_at is None
    assert persisted.cancellation_acknowledged_at is None
    assert events[-1].event_type == "execution_stopped_unexpectedly"
    assert events[-1].context["previous_status"] == "planned"


def test_stopped_callback_preserves_natural_terminal_race(tmp_path, monkeypatch):
    configured = worker_settings(tmp_path)
    assignment = create_planned_run(
        configured,
        "naturally-finished-run",
        assign_queue=configured.queue_name,
    )
    assert assignment is not None
    _configure_worker_environment(monkeypatch, configured)
    claimed = _prepare_running_assignment(configured, assignment, request_cancel=True)
    persistence = open_run_persistence(configured.database_url)
    try:
        service = create_default_run_service(
            create_default_workflow_registry(),
            repository=persistence.repository,
        )
        succeeded = service.transition_run(assignment.run_id, RunStatus.SUCCEEDED)
        events_before = service.list_events(assignment.run_id, limit=100)
    finally:
        persistence.close()

    handle_execution_stopped(_execution_job(claimed), fakeredis.FakeRedis())

    record, events, persisted = _read_run(configured, assignment.run_id)
    assert record == succeeded
    assert events == events_before
    assert persisted is not None
    assert persisted.cancellation_acknowledged_at is None


@pytest.mark.parametrize("field", ["job_id", "queue_name", "run_id", "func_name"])
def test_stopped_callback_identity_drift_never_mutates_another_run(
    tmp_path,
    monkeypatch,
    field,
):
    configured = worker_settings(tmp_path)
    assignment = create_planned_run(
        configured,
        f"stopped-drift-{field}",
        assign_queue=configured.queue_name,
    )
    assert assignment is not None
    _configure_worker_environment(monkeypatch, configured)
    claimed = _prepare_running_assignment(configured, assignment, request_cancel=True)
    job = _execution_job(claimed)
    if field == "job_id":
        job.id = "wrong-job"
    elif field == "queue_name":
        job.origin = "wrong-queue"
    elif field == "run_id":
        job.args = ("another-run",)
    else:
        job.func_name = "other.worker"
    record_before, events_before, persisted_before = _read_run(
        configured,
        assignment.run_id,
    )

    handle_execution_stopped(job, fakeredis.FakeRedis())

    record, events, persisted = _read_run(configured, assignment.run_id)
    assert record == record_before
    assert events == events_before
    assert persisted == persisted_before


def test_stopped_callback_does_not_swallow_rq_callback_timeout(monkeypatch):
    monkeypatch.setattr(
        worker_jobs,
        "load_worker_settings",
        lambda: (_ for _ in ()).throw(JobTimeoutException("callback deadline")),
    )

    with pytest.raises(JobTimeoutException, match="callback deadline"):
        handle_execution_stopped(
            SimpleNamespace(
                id="job-1",
                origin="runs",
                func_name="encode_pipeline.workers.jobs.run_execution_job",
                args=("run-1",),
                kwargs={},
            ),
            fakeredis.FakeRedis(),
        )


def _prepare_running_assignment(
    configured,
    assignment,
    *,
    request_cancel=False,
    managed_container_scope=None,
    managed_container_endpoint_identity=None,
):
    persistence = open_run_persistence(configured.database_url)
    try:
        service = create_default_run_service(
            create_default_workflow_registry(),
            repository=persistence.repository,
        )
        service.mark_execution_dispatched(
            assignment.run_id,
            job_id=assignment.job_id,
        )
        service.queue_dispatched_run(
            assignment.run_id,
            job_id=assignment.job_id,
            backend=assignment.backend,
            queue_name=assignment.queue_name,
        )
        service.claim_execution_assignment(
            assignment.run_id,
            job_id=assignment.job_id,
            backend=assignment.backend,
            queue_name=assignment.queue_name,
        )
        claimed = service.get_execution_assignment(assignment.run_id)
        assert claimed is not None
        if (
            managed_container_scope is not None
            or managed_container_endpoint_identity is not None
        ):
            claimed = service.bind_execution_cleanup_identity(
                assignment.run_id,
                expected_assignment=claimed,
                managed_container_scope=managed_container_scope,
                managed_container_endpoint_identity=(
                    managed_container_endpoint_identity
                ),
            )
        service.transition_run(assignment.run_id, RunStatus.RUNNING)
        claimed = service.get_execution_assignment(assignment.run_id)
        assert claimed is not None
        if request_cancel:
            requested = service.request_execution_cancellation(
                assignment.run_id,
                job_id=claimed.job_id,
                backend=claimed.backend,
                queue_name=claimed.queue_name,
                reason="User requested cancellation.",
            )
            claimed = requested.assignment
        return claimed
    finally:
        persistence.close()


def _execution_job(assignment):
    return SimpleNamespace(
        id=assignment.job_id,
        origin=assignment.queue_name,
        func_name="encode_pipeline.workers.jobs.run_execution_job",
        args=(assignment.run_id,),
        kwargs={},
    )
