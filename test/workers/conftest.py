"""Shared worker-boundary test helpers."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import fakeredis
import pytest
import yaml

from encode_pipeline.persistence.runtime import open_run_persistence
from encode_pipeline.platform.adapters import WorkflowInputs
from encode_pipeline.platform.builds import WorkflowBuildIdentity
from encode_pipeline.platform.execution import RunExecutionAssignment
from encode_pipeline.platform.runs import RunStatus
from encode_pipeline.services.defaults import (
    create_default_workflow_build_identity_provider,
    create_default_workflow_registry,
)
from encode_pipeline.services.materialization import WorkspaceMaterializer
from encode_pipeline.services.planning import ExecutionPlanner, WorkspacePlanner
from encode_pipeline.services.private_reference_profiles import (
    load_private_reference_profile_config,
)
from encode_pipeline.services.reference_profile_runtime import (
    ReferenceProfileBindingService,
    ReferenceProfileRuntimeResolver,
)
from encode_pipeline.services.reference_profiles import ReferenceProfileService
from encode_pipeline.services.runs import RunService
from encode_pipeline.services.validated_inputs import (
    ValidatedInputService,
    ValidatedRunCreationService,
)
from encode_pipeline.services.validation import ValidationService
from encode_pipeline.workers.settings import WorkerSettings


WORKFLOW_ID = "encode-style-chipseq-cuttag-atac-mnase"
PROFILE_ROOT = Path(__file__).resolve().parents[1] / "profiles" / "platform_worker_tiny"
REFERENCE_CONFIG_KEY = "worker-test-reference"


def _write_reference_profile_config(tmp_path: Path) -> Path:
    reference_root = tmp_path / "operator-reference"
    reference_root.mkdir(exist_ok=True)
    resources: dict[str, dict[str, str]] = {}
    contents = {
        "blacklist": b"chr1\t0\t1\n",
        "chrom_sizes": b"chr1\t4\n",
        "gtf": b'chr1\ttest\texon\t1\t4\t.\t+\t.\tgene_id "g1";\n',
        "reference_fasta": b">chr1\nACGT\n",
    }
    for name, content in contents.items():
        path = (reference_root / name).resolve()
        path.write_bytes(content)
        resources[name] = {
            "path": str(path),
            "sha256": hashlib.sha256(content).hexdigest(),
        }
    prefix = (reference_root / "bowtie2" / "genome").resolve()
    prefix.parent.mkdir(exist_ok=True)
    index_files: dict[str, str] = {}
    for suffix in (
        ".1.bt2",
        ".2.bt2",
        ".3.bt2",
        ".4.bt2",
        ".rev.1.bt2",
        ".rev.2.bt2",
    ):
        content = suffix.encode("ascii")
        Path(f"{prefix}{suffix}").write_bytes(content)
        index_files[suffix] = hashlib.sha256(content).hexdigest()
    config_path = (tmp_path / "operator-reference-profiles.json").resolve()
    config_path.write_text(
        json.dumps(
            {
                "schema_version": "helixweave-reference-profiles-v1",
                "profiles": {
                    REFERENCE_CONFIG_KEY: {
                        "bindings": {
                            WORKFLOW_ID: {
                                "schema_version": "encode-reference-binding-v1",
                                "assembly": "GRCh38",
                                "effective_genome_size": "hs",
                                "genome_resources": resources,
                                "bowtie2_index": {
                                    "prefix": str(prefix),
                                    "files": index_files,
                                },
                            }
                        }
                    }
                },
            }
        ),
        encoding="utf-8",
    )
    return config_path


@pytest.fixture(autouse=True)
def _fakeredis_client_list_includes_real_redis_addresses(monkeypatch):
    """Supply the CLIENT LIST fields that real Redis guarantees to RQ."""
    client_list = fakeredis.FakeRedis.client_list

    def client_list_with_addresses(connection, *args, **kwargs):
        clients = client_list(connection, *args, **kwargs)
        for client in clients:
            client.setdefault("addr", "127.0.0.1:0")
            client.setdefault("laddr", "127.0.0.1:6379")
        return clients

    monkeypatch.setattr(
        fakeredis.FakeRedis,
        "client_list",
        client_list_with_addresses,
    )


def worker_settings(tmp_path: Path, queue_name: str = "worker-tests") -> WorkerSettings:
    """Return isolated file-backed settings for one worker test."""
    return WorkerSettings(
        database_url=f"sqlite:///{tmp_path / 'platform.db'}",
        redis_url="redis://unused.test/0",
        queue_name=queue_name,
        workspace_root=tmp_path / "workspaces",
        reference_profile_config=_write_reference_profile_config(tmp_path),
    )


def create_planned_run(
    settings: WorkerSettings,
    run_id: str,
    *,
    assign_queue: str | None = None,
    bind_build_identity: bool = True,
    build_identity: WorkflowBuildIdentity | None = None,
    enable_qc_summary: bool = False,
    samples_path: Path | None = None,
) -> RunExecutionAssignment | None:
    """Persist a PLANNED run and optionally its canonical execution job ID."""
    persistence = open_run_persistence(settings.database_url)
    try:
        registry = create_default_workflow_registry()
        assert settings.reference_profile_config is not None

        def private_config_provider():
            return load_private_reference_profile_config(
                settings.reference_profile_config
            )

        profile_service = ReferenceProfileService(
            repository=persistence.reference_profile_repository,
            private_config_provider=private_config_provider,
            adapter_provider=registry.get,
        )
        enabled_profiles = profile_service.list_enabled(WORKFLOW_ID)
        if enabled_profiles:
            reference_revision = enabled_profiles[0]
        else:
            registered = profile_service.register(
                safe_key="worker-test-reference",
                display_name="Worker test GRCh38",
                organism="Homo sapiens",
                assembly="GRCh38",
                config_key=REFERENCE_CONFIG_KEY,
            )
            reference_revision = profile_service.enable(
                registered.profile_id,
                revision_id=registered.revision_id,
            )
        binding_service = ReferenceProfileBindingService(
            repository=persistence.reference_profile_repository,
            private_config_provider=private_config_provider,
            adapter_provider=registry.get,
        )
        service = RunService(
            registry,
            id_factory=lambda: run_id,
            repository=persistence.repository,
        )
        config = yaml.safe_load(
            (PROFILE_ROOT / "config.yaml").read_text(encoding="utf-8")
        )
        configured_samples = (
            (PROFILE_ROOT / "samples.tsv").resolve()
            if samples_path is None
            else samples_path.resolve()
        )
        config.pop("samples", None)
        config.pop("genome_resources", None)
        if enable_qc_summary:
            config["qc"]["summary"] = True
        with configured_samples.open(encoding="utf-8", newline="") as handle:
            samples = list(csv.DictReader(handle, delimiter="\t"))
        for sample in samples:
            sample.pop("genome", None)
            sample.pop("bowtie2_index", None)
        inputs = WorkflowInputs(config=config, samples=samples)
        build_provider = create_default_workflow_build_identity_provider(
            registry=registry
        )
        validated = ValidatedInputService(
            registry=registry,
            validation_service=ValidationService(registry),
            build_identity_provider=build_provider,
            repository=persistence.repository,
            reference_profile_binding_service=binding_service,
            reference_profile_catalog=profile_service,
        ).validate(
            WORKFLOW_ID,
            inputs,
            reference_profile_revision_id=reference_revision.revision_id,
        )
        assert validated.is_success and validated.value is not None
        record = (
            ValidatedRunCreationService(
                run_service=service,
                build_identity_provider=build_provider,
                reference_profile_binding_service=binding_service,
            )
            .create_run(WORKFLOW_ID, validated.value.snapshot_id)
            .record
        )
        assert record.run_id == run_id
        service.transition_run(run_id, RunStatus.VALIDATING, stage="preflight")
        base_result = ExecutionPlanner(service).plan_run(run_id)
        assert base_result.is_success
        base_plan = base_result.value
        workspace_dir = settings.workspace_root / run_id
        reference_resolver = ReferenceProfileRuntimeResolver(
            persistence.repository,
            registry,
            binding_service,
        )
        workspace_result = WorkspacePlanner(
            registry,
            reference_profile_resolver=reference_resolver,
        ).plan_workspace(
            base_plan,
            workspace_dir,
        )
        assert workspace_result.is_success
        workspace_plan = workspace_result.value
        materialized = WorkspaceMaterializer().materialize(
            workspace_plan.workspace_plan,
            workspace_dir,
        )
        assert materialized.is_success
        if bind_build_identity:
            if build_identity is None:
                build_identity = validated.value.workflow_build_identity
            service.complete_preflight(
                run_id,
                build_identity,
                stage="preflight",
            )
        else:
            service.transition_run(run_id, RunStatus.PLANNED, stage="preflight")
        if assign_queue is None:
            return None
        assignment = service.ensure_execution_assignment(
            run_id,
            queue_name=assign_queue,
        )
        return assignment
    finally:
        persistence.close()
