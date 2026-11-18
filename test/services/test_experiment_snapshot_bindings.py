"""Memory freeze and replay for an experiment revision bound to a snapshot."""

from __future__ import annotations

from dataclasses import replace
from datetime import datetime, timedelta, timezone
from hashlib import sha256
from itertools import count
import json
from unittest.mock import patch

import pytest

from encode_pipeline.platform.adapters import WorkflowInputs
from encode_pipeline.platform.builds import WorkflowBuildIdentity
from encode_pipeline.api.catalog_models import CatalogExperimentRunListResponse
from encode_pipeline.platform.experiment_catalog import (
    CatalogError,
    ExperimentRevisionFreezeSource,
    association_pins_match,
    freeze_experiment_snapshot_binding,
)
from encode_pipeline.platform.runs import RunArtifactRef, RunRecord, RunStatus
from encode_pipeline.platform.data_registry import ProjectSampleSelection
from encode_pipeline.platform.snapshots import (
    PAYLOAD_DIGEST_SCHEME,
    VALIDATION_EVIDENCE_OUTCOME,
    ValidatedInputSnapshot,
    build_workflow_inputs_digest,
    canonical_workflow_inputs_json,
)
from encode_pipeline.services.data_registry import DataRegistryService
from encode_pipeline.services.data_registry_repositories import (
    InMemoryDataRegistryRepository,
)
from encode_pipeline.services.experiment_catalog import ExperimentCatalogService
from encode_pipeline.services.experiment_catalog_repositories import (
    InMemoryExperimentCatalogRepository,
)
from encode_pipeline.services.run_repositories import (
    InMemoryRunRepository,
    RunEventDraft,
    ValidatedSnapshotReplayConflictError,
)


NOW = datetime(2026, 11, 19, tzinfo=timezone.utc)


def _ids(prefix: str):
    serial = count(1)
    return lambda: f"{prefix}{next(serial):032x}"


def _write(operation: str, body: dict) -> dict:
    key = "".join(character if character.isalnum() else "-" for character in operation)
    return {
        "actor_user_id": "usr_" + "a" * 32,
        "idempotency_key": key[:120],
        "raw_body": json.dumps(body, sort_keys=True).encode(),
    }


def _world(*, attach_trace: bool = True):
    registry = InMemoryDataRegistryRepository(legacy_created_at=NOW)
    data = DataRegistryService(
        repository=registry,
        project_id_factory=_ids("prj_"),
        sample_id_factory=_ids("smp_"),
        sample_revision_id_factory=_ids("smpr_"),
        now_factory=lambda: NOW,
    )
    catalog_repository = InMemoryExperimentCatalogRepository(registry)
    catalog = ExperimentCatalogService(
        catalog_repository,
        now_factory=lambda: NOW,
        term_id_factory=_ids("term_"),
        term_revision_id_factory=_ids("trmv_"),
        experiment_id_factory=_ids("exp_"),
        experiment_revision_id_factory=_ids("expr_"),
        sample_revision_id_factory=_ids("smpr_"),
    )
    runs = InMemoryRunRepository(
        data_registry_repository=registry,
        experiment_freeze_lookup=catalog_repository,
    )
    if attach_trace:
        catalog_repository.attach_experiment_trace_evidence(runs)
    return data, catalog, runs, catalog_repository


def _term(catalog: ExperimentCatalogService, term_type: str, name: str) -> dict:
    body = {
        "term_type": term_type,
        "display_name": name,
        "namespace": None,
        "accession": None,
        "aliases": [],
    }
    return catalog.create_term(**_write(f"term-{name}", body), **body)


def _known(term: dict) -> dict:
    return {
        "status": "known",
        "term_id": term["term_id"],
        "term_revision_id": term["revision"]["term_revision_id"],
    }


def _snapshot(snapshot_id: str) -> ValidatedInputSnapshot:
    canonical_payload = canonical_workflow_inputs_json(
        WorkflowInputs(config={"threads": 1}, samples=[{"sample": "row"}])
    )
    return ValidatedInputSnapshot(
        snapshot_id=snapshot_id,
        workflow_id="workflow-a",
        adapter_version="1.0.0",
        schema_version="1.0.0",
        schema_dialect="https://json-schema.org/draft/2020-12/schema",
        workflow_build_identity=WorkflowBuildIdentity(
            workflow_id="workflow-a",
            adapter_version="1.0.0",
            scheme="sha256-tree-v1",
            logical_entrypoint="workflow/Snakefile",
            digest="a" * 64,
            captured_at=NOW,
        ),
        canonical_payload=canonical_payload,
        payload_digest_scheme=PAYLOAD_DIGEST_SCHEME,
        payload_digest=build_workflow_inputs_digest(canonical_payload),
        validation_outcome=VALIDATION_EVIDENCE_OUTCOME,
        validation_issue_codes=(),
        validated_at=NOW,
        expires_at=NOW + timedelta(minutes=30),
    )


def _record(run_id: str, snapshot: ValidatedInputSnapshot) -> RunRecord:
    return RunRecord(
        run_id=run_id,
        workflow_id="workflow-a",
        inputs=snapshot.to_workflow_inputs().to_dict(),
        status=RunStatus.CREATED,
        created_at=NOW + timedelta(minutes=1),
        updated_at=NOW + timedelta(minutes=1),
        started_at=None,
        ended_at=None,
        current_stage=None,
        cancellation_reason=None,
        error=None,
        tags={},
    )


def _event() -> RunEventDraft:
    return RunEventDraft(
        event_type="status_changed",
        message="Run created.",
        status=RunStatus.CREATED,
        context={"previous_status": None, "new_status": "created"},
    )


def _experiment(data, catalog, sample_count: int) -> dict:
    assay = _term(catalog, "assay", "CUT&Tag")
    target = _term(catalog, "target", "H3K27ac")
    batch = _term(catalog, "batch", "batch-1")
    role = _term(catalog, "replicate_role", "biological replicate")
    project = data.create_project("Pilot")
    imported = []
    for index in range(1, sample_count + 1):
        created_sample = data.import_sample(
            project.project_id,
            stable_key=f"sample-{index}",
            display_name=f"Sample {index}",
            attributes={"organism": "human"},
        )
        imported.append(created_sample.revision.sample_revision_id)
    body = {
        "stable_key": "cut-pair",
        "display_name": "cut-pair",
        "assay": _known(assay),
        "target": _known(target),
        "batch": _known(batch),
        "samples": [
            {"sample_revision_id": revision_id, "replicate_role": _known(role)}
            for revision_id in imported
        ],
    }
    created = catalog.create_experiment(
        **_write("experiment-cut-pair", body),
        project_id=project.project_id,
        **body,
    )
    return {
        "project_id": project.project_id,
        "sample_revision_ids": tuple(imported),
        "experiment": created,
    }


def test_memory_freeze_ignores_order_and_replay_compares_both_fields() -> None:
    data, catalog, runs, _repository = _world()
    seeded = _experiment(data, catalog, 2)
    left, right = seeded["sample_revision_ids"]
    snapshot = _snapshot("vsnap_" + "1" * 32)
    runs.create_validated_input_snapshot(
        snapshot,
        project_sample_selection=ProjectSampleSelection(
            project_id=seeded["project_id"],
            sample_revision_ids=(right, left),
        ),
        experiment_revision_id=seeded["experiment"]["experiment_revision_id"],
    )
    frozen = runs.get_validated_experiment_binding(snapshot.snapshot_id)
    assert frozen is not None
    assert (
        frozen.experiment_revision_id
        == (seeded["experiment"]["experiment_revision_id"])
    )
    assert len(frozen.payload_digest) == 64

    runs.consume_validated_input_snapshot(
        snapshot.snapshot_id,
        workflow_id="workflow-a",
        expected_build_identity=snapshot.workflow_build_identity,
        record=_record("run-bound", snapshot),
        consumed_at=NOW + timedelta(minutes=1),
        event=_event(),
    )
    copied = runs._run_experiment_bindings["run-bound"]
    assert copied == frozen

    object.__setattr__(copied, "payload_digest", "b" * 64)
    with pytest.raises(ValidatedSnapshotReplayConflictError):
        runs.consume_validated_input_snapshot(
            snapshot.snapshot_id,
            workflow_id="workflow-a",
            expected_build_identity=snapshot.workflow_build_identity,
            record=_record("run-bound", snapshot),
            consumed_at=NOW + timedelta(minutes=1),
            event=_event(),
        )
    object.__setattr__(copied, "payload_digest", frozen.payload_digest)
    object.__setattr__(copied, "experiment_revision_id", "expr_" + "c" * 32)
    with pytest.raises(ValidatedSnapshotReplayConflictError):
        runs.consume_validated_input_snapshot(
            snapshot.snapshot_id,
            workflow_id="workflow-a",
            expected_build_identity=snapshot.workflow_build_identity,
            record=_record("run-bound", snapshot),
            consumed_at=NOW + timedelta(minutes=1),
            event=_event(),
        )


def test_memory_mismatch_and_archive_leave_no_snapshot() -> None:
    data, catalog, runs, _repository = _world()
    seeded = _experiment(data, catalog, 2)
    only_one = seeded["sample_revision_ids"][:1]
    snapshot = _snapshot("vsnap_" + "2" * 32)
    with pytest.raises(CatalogError) as mismatch:
        runs.create_validated_input_snapshot(
            snapshot,
            project_sample_selection=ProjectSampleSelection(
                project_id=seeded["project_id"],
                sample_revision_ids=only_one,
            ),
            experiment_revision_id=seeded["experiment"]["experiment_revision_id"],
        )
    assert mismatch.value.code == "catalog_sample_set_mismatch"
    with pytest.raises(KeyError):
        runs.get_validated_input_snapshot(snapshot.snapshot_id)

    catalog.archive_experiment(
        **_write(
            "archive-cut-pair",
            {"expected_revision_number": 1},
        ),
        experiment_id=seeded["experiment"]["experiment_id"],
        expected_revision_number=1,
    )
    archived_snapshot = _snapshot("vsnap_" + "3" * 32)
    with pytest.raises(CatalogError) as archived:
        runs.create_validated_input_snapshot(
            archived_snapshot,
            project_sample_selection=ProjectSampleSelection(
                project_id=seeded["project_id"],
                sample_revision_ids=seeded["sample_revision_ids"],
            ),
            experiment_revision_id=seeded["experiment"]["experiment_revision_id"],
        )
    assert archived.value.code == "catalog_experiment_archived"
    with pytest.raises(KeyError):
        runs.get_validated_input_snapshot(archived_snapshot.snapshot_id)


def test_memory_omitted_experiment_stays_unbound_on_replay() -> None:
    data, catalog, runs, _repository = _world()
    seeded = _experiment(data, catalog, 1)
    snapshot = _snapshot("vsnap_" + "4" * 32)
    runs.create_validated_input_snapshot(
        snapshot,
        project_sample_selection=ProjectSampleSelection(
            project_id=seeded["project_id"],
            sample_revision_ids=seeded["sample_revision_ids"],
        ),
    )
    assert runs.get_validated_experiment_binding(snapshot.snapshot_id) is None
    runs.consume_validated_input_snapshot(
        snapshot.snapshot_id,
        workflow_id="workflow-a",
        expected_build_identity=snapshot.workflow_build_identity,
        record=_record("run-open", snapshot),
        consumed_at=NOW + timedelta(minutes=1),
        event=_event(),
    )
    assert "run-open" not in runs._run_experiment_bindings
    replay = runs.consume_validated_input_snapshot(
        snapshot.snapshot_id,
        workflow_id="workflow-a",
        expected_build_identity=snapshot.workflow_build_identity,
        record=_record("run-open", snapshot),
        consumed_at=NOW + timedelta(minutes=1),
        event=_event(),
    )
    assert replay.created is False
    trace = catalog.list_experiment_run_results(
        project_id=seeded["project_id"],
        experiment_revision_id=seeded["experiment"]["experiment_revision_id"],
        cursor=None,
        limit=50,
    )
    assert trace["trace_status"] == "unresolved"
    assert trace["reason"] == "experiment_snapshot_binding_absent"


def test_association_pins_match_requires_ordinal_group_and_role() -> None:
    first = (
        0,
        "smpr_" + "1" * 32,
        "a" * 64,
        None,
        "known",
        "term_" + "1" * 32,
        "trmv_" + "1" * 32,
    )
    second = (
        1,
        "smpr_" + "2" * 32,
        "b" * 64,
        "group-a",
        "known",
        "term_" + "2" * 32,
        "trmv_" + "2" * 32,
    )
    assert association_pins_match((first, second), (second, first))
    drifted = (
        (11, *first[1:]),
        (*first[:3], "drifted-group", *first[4:]),
        (*first[:4], "not_applicable", None, None),
    )
    for mirror in drifted:
        assert not association_pins_match((first, second), (mirror, second))


def test_freeze_rejects_a_mirror_that_drops_a_payload_pin() -> None:
    revision_id = "expr_" + "d" * 32
    sample_ids = ("smpr_" + "1" * 32, "smpr_" + "2" * 32)
    source = ExperimentRevisionFreezeSource(
        experiment_revision_id=revision_id,
        project_id="prj_" + "1" * 32,
        project_active=True,
        experiment_archived=False,
        payload_digest="a" * 64,
        sample_revision_ids=sample_ids,
        association_matches_payload=False,
    )
    with pytest.raises(CatalogError) as caught:
        freeze_experiment_snapshot_binding(
            source,
            project_id=source.project_id,
            sample_revision_ids=sample_ids,
        )
    assert caught.value.code == "catalog_sample_set_mismatch"


def test_memory_trace_reads_bindings_only_after_the_store_is_attached() -> None:
    data, catalog, runs, repository = _world(attach_trace=False)
    seeded = _experiment(data, catalog, 2)
    revision_id = seeded["experiment"]["experiment_revision_id"]

    def trace() -> dict:
        return catalog.list_experiment_run_results(
            project_id=seeded["project_id"],
            experiment_revision_id=revision_id,
            cursor=None,
            limit=50,
        )

    snapshot = _snapshot("vsnap_" + "5" * 32)
    runs.create_validated_input_snapshot(
        snapshot,
        project_sample_selection=ProjectSampleSelection(
            project_id=seeded["project_id"],
            sample_revision_ids=tuple(reversed(seeded["sample_revision_ids"])),
        ),
        experiment_revision_id=revision_id,
    )
    isolated = trace()
    assert isolated["trace_status"] == "unresolved"
    assert isolated["reason"] == "experiment_snapshot_binding_absent"

    repository.attach_experiment_trace_evidence(runs)
    frozen = trace()
    assert frozen["trace_status"] == "experiment_snapshot"
    assert "reason" not in frozen
    assert frozen["results"] == []

    runs.consume_validated_input_snapshot(
        snapshot.snapshot_id,
        workflow_id="workflow-a",
        expected_build_identity=snapshot.workflow_build_identity,
        record=_record("run-trace", snapshot),
        consumed_at=NOW + timedelta(minutes=1),
        event=_event(),
    )
    consumed = trace()
    assert consumed["trace_status"] == "experiment_snapshot"
    assert "reason" not in consumed
    assert consumed["results"] == []

    succeeded_at = NOW + timedelta(minutes=2)
    current = runs.get_run("run-trace")
    runs.update_run(
        replace(
            current,
            status=RunStatus.SUCCEEDED,
            updated_at=succeeded_at,
            ended_at=succeeded_at,
        ),
        expected_status=RunStatus.CREATED,
        event=RunEventDraft(
            event_type="status_changed",
            message="Scientific execution succeeded.",
            status=RunStatus.SUCCEEDED,
            context={"previous_status": "created", "new_status": "succeeded"},
        ),
    )
    attempt_id = "resultattempt-" + sha256(b"trace").hexdigest()
    runs.begin_artifact_result_attempt(
        "run-trace",
        attempt_id=attempt_id,
        expected_status=RunStatus.SUCCEEDED,
    )
    runs.replace_artifacts(
        "run-trace",
        (
            RunArtifactRef(
                artifact_id="peaks",
                run_id="run-trace",
                artifact_type="file",
                name="peaks.bed",
                uri="run://runs/run-trace/artifacts/peaks",
                mime_type="text/plain",
                produced_at=succeeded_at,
                revision="artifactrev-" + sha256(b"peaks").hexdigest(),
                metadata={"output_type": "peaks"},
            ),
        ),
        attempt_id=attempt_id,
        expected_status=RunStatus.SUCCEEDED,
        event=RunEventDraft(
            event_type="artifacts_indexed",
            message="Workflow artifacts indexed.",
            status=RunStatus.SUCCEEDED,
        ),
    )
    published = trace()
    assert published["trace_status"] == "experiment_snapshot"
    assert "reason" not in published
    assert len(published["results"]) == 1
    assert published["results"][0]["run_id"] == "run-trace"
    assert published["results"][0]["publication_evidence"] == "present"
    assert published["results"][0]["artifact_association_kind"] == "run_input_scope"

    copied = runs._run_experiment_bindings["run-trace"]
    object.__setattr__(copied, "payload_digest", "b" * 64)
    hidden = trace()
    assert hidden["trace_status"] == "experiment_snapshot"
    assert "reason" not in hidden
    assert hidden["results"] == []


def _publish(
    runs,
    run_id: str,
    stamp: datetime,
    artifact_ids: tuple[str, ...],
) -> None:
    current = runs.get_run(run_id)
    if current.status is RunStatus.CREATED:
        runs.update_run(
            replace(
                current,
                status=RunStatus.SUCCEEDED,
                updated_at=stamp,
                ended_at=stamp,
            ),
            expected_status=RunStatus.CREATED,
            event=RunEventDraft(
                event_type="status_changed",
                message="Scientific execution succeeded.",
                status=RunStatus.SUCCEEDED,
                context={"previous_status": "created", "new_status": "succeeded"},
            ),
        )
    attempt_id = (
        "resultattempt-"
        + sha256(f"{run_id}|{'|'.join(artifact_ids)}".encode()).hexdigest()
    )
    runs.begin_artifact_result_attempt(
        run_id,
        attempt_id=attempt_id,
        expected_status=RunStatus.SUCCEEDED,
    )

    class PublicationClock(datetime):
        @classmethod
        def now(cls, tz=None):
            return stamp.astimezone(tz)

    with patch("encode_pipeline.services.run_repositories.datetime", PublicationClock):
        runs.replace_artifacts(
            run_id,
            tuple(
                RunArtifactRef(
                    artifact_id=artifact_id,
                    run_id=run_id,
                    artifact_type="file",
                    name=f"{artifact_id}.bed",
                    uri=f"run://runs/{run_id}/artifacts/{artifact_id}",
                    mime_type="text/plain",
                    produced_at=stamp,
                    revision="artifactrev-" + sha256(artifact_id.encode()).hexdigest(),
                    metadata={"output_type": artifact_id},
                )
                for artifact_id in artifact_ids
            ),
            attempt_id=attempt_id,
            expected_status=RunStatus.SUCCEEDED,
            event=RunEventDraft(
                event_type="artifacts_indexed",
                message="Workflow artifacts indexed.",
                status=RunStatus.SUCCEEDED,
            ),
        )


def test_memory_trace_pages_by_publication_instant() -> None:
    data, catalog, runs, _ = _world()
    seeded = _experiment(data, catalog, 2)
    revision_id = seeded["experiment"]["experiment_revision_id"]
    selection = ProjectSampleSelection(
        project_id=seeded["project_id"],
        sample_revision_ids=tuple(seeded["sample_revision_ids"]),
    )
    for digit, run_id, micros in (("1", "run-a", 0), ("2", "run-b", 500_000)):
        snapshot = _snapshot("vsnap_" + digit * 32)
        runs.create_validated_input_snapshot(
            snapshot,
            project_sample_selection=selection,
            experiment_revision_id=revision_id,
        )
        runs.consume_validated_input_snapshot(
            snapshot.snapshot_id,
            workflow_id="workflow-a",
            expected_build_identity=snapshot.workflow_build_identity,
            record=_record(run_id, snapshot),
            consumed_at=NOW + timedelta(minutes=1),
            event=_event(),
        )
        _publish(
            runs,
            run_id,
            NOW + timedelta(minutes=2, microseconds=micros),
            ("peaks",),
        )
    cursor = None
    seen: list[tuple[str, str]] = []
    for _ in range(2):
        page = catalog.list_experiment_run_results(
            project_id=seeded["project_id"],
            experiment_revision_id=revision_id,
            cursor=cursor,
            limit=1,
        )
        CatalogExperimentRunListResponse.model_validate(page)
        assert page["trace_status"] == "experiment_snapshot"
        assert len(page["results"]) == 1
        row = page["results"][0]
        seen.append((row["run_id"], row["published_at"]))
        cursor = page["next_cursor"]
        assert (cursor is not None) == (len(seen) == 1)
    assert seen == [
        ("run-b", "2026-11-19T00:02:00.500000Z"),
        ("run-a", "2026-11-19T00:02:00Z"),
    ]
    assert cursor is None


def test_memory_trace_orders_the_same_instant_by_run_then_artifact() -> None:
    data, catalog, runs, _ = _world()
    seeded = _experiment(data, catalog, 2)
    revision_id = seeded["experiment"]["experiment_revision_id"]
    selection = ProjectSampleSelection(
        project_id=seeded["project_id"],
        sample_revision_ids=tuple(seeded["sample_revision_ids"]),
    )
    stamp = NOW + timedelta(minutes=2)
    for digit, run_id, artifact_ids in (
        ("3", "run-a", ("peaks",)),
        ("4", "run-b", ("peaks", "signal")),
    ):
        snapshot = _snapshot("vsnap_" + digit * 32)
        runs.create_validated_input_snapshot(
            snapshot,
            project_sample_selection=selection,
            experiment_revision_id=revision_id,
        )
        runs.consume_validated_input_snapshot(
            snapshot.snapshot_id,
            workflow_id="workflow-a",
            expected_build_identity=snapshot.workflow_build_identity,
            record=_record(run_id, snapshot),
            consumed_at=NOW + timedelta(minutes=1),
            event=_event(),
        )
        _publish(runs, run_id, stamp, artifact_ids)
    cursor = None
    seen: list[tuple[str, str]] = []
    for _ in range(3):
        page = catalog.list_experiment_run_results(
            project_id=seeded["project_id"],
            experiment_revision_id=revision_id,
            cursor=cursor,
            limit=1,
        )
        CatalogExperimentRunListResponse.model_validate(page)
        assert len(page["results"]) == 1
        row = page["results"][0]
        seen.append((row["run_id"], row["artifact_id"]))
        cursor = page["next_cursor"]
    assert seen == [
        ("run-b", "signal"),
        ("run-b", "peaks"),
        ("run-a", "peaks"),
    ]
    assert cursor is None
