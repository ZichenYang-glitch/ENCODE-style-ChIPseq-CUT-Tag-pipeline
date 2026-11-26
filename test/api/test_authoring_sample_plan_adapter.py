"""Public drafts keep ENCODE sample plans on the selected workflow adapter."""

from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import Mock

import pytest

from encode_pipeline.adapters.hitrac_preprocess.adapter import HiTracPreprocessAdapter
from encode_pipeline.api.main import create_app
from encode_pipeline.platform.authoring import WorkflowSamplePlanAdapter
from encode_pipeline.platform.input_registry import InputFileRevisionSelection
from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.platform.results import Result
from encode_pipeline.services.authoring import _sample_plan_adapter
from encode_pipeline.services.data_registry import DataRegistryService, SampleImportSpec
from api_test_client import ApiTestClient
from conftest import seed_test_authentication


WORKFLOW_ID = "encode-style-chipseq-cuttag-atac-mnase"
HITRAC_WORKFLOW_ID = "hitrac-preprocess"
BULK_WORKFLOW_ID = "bulk-rnaseq"


def test_declining_or_unrelated_adapters_do_not_own_sample_plans() -> None:
    assert _sample_plan_adapter(object()) is None
    assert _sample_plan_adapter(_DecliningSamplePlan()) is None


def test_public_encode_sample_plan_round_trip_reuses_the_snapshot(
    tmp_path, reference_ready_app, monkeypatch
):
    app = reference_ready_app
    root, project, sample_revision_id, reads = _managed_reads(app, tmp_path)
    plan = _encode_plan(
        sample_revision_id,
        [
            {"read": 2, "input_file_revision_id": reads["r2"]},
            {"read": 1, "input_file_revision_id": reads["r1"]},
        ],
    )
    validated_service = app.state.validated_input_service
    validate = Mock(wraps=validated_service.validate)
    monkeypatch.setattr(validated_service, "validate", validate)
    with ApiTestClient(app) as client:
        url = _session_url(client, project, WORKFLOW_ID)
        saved = client.request(
            "PUT",
            url + "/draft",
            json=_draft(
                0,
                sample_revision_id,
                app.state.test_reference_profile.revision_id,
                [plan],
            ),
        )
        assert saved.status_code == 200, saved.text
        opened = client.get(url)
        assert opened.status_code == 200, opened.text
        stored = opened.json()["draft"]["workflow_samples"][0]
        assert stored["assay"] == "chipseq"
        assert stored["layout"] == "PE"
        assert stored["role"] == "treatment"
        assert stored["target"] == "CTCF"
        assert stored["peak_mode"] == "narrow"
        assert stored["control_sample_revision_id"] is None
        assert [item["read"] for item in stored["reads"]] == [1, 2]
        assert [item["input_file_revision_id"] for item in stored["reads"]] == [
            reads["r1"],
            reads["r2"],
        ]
        assert str(root) not in opened.text
        resaved = client.request(
            "PUT",
            url + "/draft",
            json={
                "expected_revision": opened.json()["draft_revision"],
                **opened.json()["draft"],
            },
        )
        assert resaved.status_code == 200, resaved.text
        assert resaved.json()["changed"] is False
        assert resaved.json()["draft_revision"] == opened.json()["draft_revision"]
        validated = client.post(
            url + "/validate",
            json={"expected_revision": resaved.json()["draft_revision"]},
        )
        assert validated.status_code == 200, validated.text
        body = validated.json()
        assert body["ok"] is True, body
        snapshot_id = body["snapshot"]["snapshot_id"]
        binding = app.state.persistence.repository.get_validated_input_use_binding(
            snapshot_id
        )
        assert [
            (use.key, use.occurrence, use.members[0].input_file_revision_id)
            for use in binding.input_uses
        ] == [
            ("sample-a", 1, reads["r1"]),
            ("sample-a", 2, reads["r2"]),
        ]
        assert str(root) not in validated.text
        again = client.post(
            url + "/validate",
            json={"expected_revision": resaved.json()["draft_revision"]},
        )
        assert again.status_code == 200, again.text
        assert again.json()["ok"] is True
        assert again.json()["snapshot"]["snapshot_id"] == snapshot_id
        validate.assert_called_once()
        # Same revision IDs under a different use cannot reuse the old snapshot.
        with monkeypatch.context() as patch:
            patch.setattr(
                app.state.registry.get(WORKFLOW_ID),
                "describe_workflow_sample_inputs",
                lambda plans, *, sample_names: (
                    InputFileRevisionSelection("other", 1, (reads["r1"],)),
                    InputFileRevisionSelection("other", 2, (reads["r2"],)),
                ),
            )
            rebound = client.post(
                url + "/validate",
                json={"expected_revision": resaved.json()["draft_revision"]},
            )
        assert rebound.json()["ok"] is False
        assert rebound.json()["snapshot"] is None
        assert validate.call_count == 2
        mismatched = client.request(
            "PUT",
            url + "/draft",
            json=_draft(
                resaved.json()["draft_revision"],
                sample_revision_id,
                app.state.test_reference_profile.revision_id,
                [plan],
                input_selections=[
                    {
                        "input_use_key": "sample-a",
                        "occurrence": 2,
                        "input_file_revision_ids": [reads["r2"]],
                    },
                    {
                        "input_use_key": "sample-a",
                        "occurrence": 1,
                        "input_file_revision_ids": [reads["r1"]],
                    },
                ],
            ),
        )
        assert mismatched.status_code == 200, mismatched.text
        rejected = client.post(
            url + "/validate",
            json={"expected_revision": mismatched.json()["draft_revision"]},
        )
        assert rejected.json()["ok"] is False
        assert rejected.json()["snapshot"] is None
        assert rejected.json()["issues"][0]["code"] == "WORKFLOW_SAMPLE_PLAN_INVALID"
        assert validate.call_count == 2


def test_empty_materializer_result_is_a_validation_issue(
    tmp_path, reference_ready_app, monkeypatch
):
    from encode_pipeline.platform.results import Result

    app = reference_ready_app
    _, project, sample_id, reads = _managed_reads(app, tmp_path)
    plan = _encode_plan(
        sample_id,
        [
            {"read": 1, "input_file_revision_id": reads["r1"]},
            {"read": 2, "input_file_revision_id": reads["r2"]},
        ],
    )
    adapter = app.state.registry.get(WORKFLOW_ID)
    monkeypatch.setattr(
        type(adapter),
        "materialize_workflow_sample_rows",
        lambda *a, **k: Result.success(None),
    )
    with ApiTestClient(app) as client:
        session = client.post(
            "/api/v1/agent/authoring-sessions",
            json={"project_id": project.project_id, "workflow_id": WORKFLOW_ID},
        )
        sid = session.json()["sid"]
        draft = client.request(
            "PUT",
            f"/api/v1/agent/authoring-sessions/{sid}/draft",
            json=_draft(
                0, sample_id, app.state.test_reference_profile.revision_id, [plan]
            ),
        )
        assert draft.status_code == 200, draft.text
        response = client.post(
            f"/api/v1/agent/authoring-sessions/{sid}/validate",
            json={"expected_revision": draft.json()["draft_revision"]},
        )
        assert response.status_code == 200, response.text
        assert response.json()["snapshot"] is None
        assert response.json()["issues"][0]["code"] == "WORKFLOW_SAMPLE_PLAN_INVALID"


def test_public_draft_keeps_encode_plan_errors_off_other_workflows(
    reference_ready_app,
):
    app = reference_ready_app
    project, revision_id = _project_sample(app)
    encode = app.state.registry.get(WORKFLOW_ID)
    hitrac = app.state.registry.get(HITRAC_WORKFLOW_ID)
    bulk = app.state.registry.get(BULK_WORKFLOW_ID)
    assert isinstance(encode, WorkflowSamplePlanAdapter)
    assert encode.supports_workflow_sample_plan() is True
    assert not isinstance(hitrac, WorkflowSamplePlanAdapter)
    assert not isinstance(bulk, WorkflowSamplePlanAdapter)
    illegal = _encode_plan(
        revision_id,
        [{"read": 1, "input_file_revision_id": "inpfr_" + ("a" * 32)}],
        assay="rna",
        peak_mode="nucleosome",
        layout="SE",
    )
    with ApiTestClient(app) as client:
        url = _session_url(client, project, WORKFLOW_ID)
        rejected = client.request(
            "PUT",
            url + "/draft",
            json=_draft(
                0,
                revision_id,
                app.state.test_reference_profile.revision_id,
                [illegal],
            ),
        )
        assert rejected.status_code == 422, rejected.text
        assert rejected.json()["reason_code"] == "INVALID_REQUEST"
        opened = client.get(url)
        assert opened.json()["draft_revision"] == 0
        assert opened.json()["draft"]["workflow_samples"] == []
        missing = client.request(
            "PUT",
            url + "/draft",
            json=_draft(
                0,
                revision_id,
                app.state.test_reference_profile.revision_id,
                [],
            ),
        )
        assert missing.status_code == 200, missing.text
        checked = client.post(
            url + "/validate",
            json={"expected_revision": missing.json()["draft_revision"]},
        )
        assert checked.json()["ok"] is False
        assert checked.json()["snapshot"] is None
        assert checked.json()["issues"][0]["code"] == "WORKFLOW_SAMPLE_PLAN_REQUIRED"
        hitrac_url = _session_url(client, project, HITRAC_WORKFLOW_ID)
        saved = client.request(
            "PUT",
            hitrac_url + "/draft",
            json={
                **_draft(0, revision_id, None, [illegal]),
                "samples": [
                    {
                        "sample_id": "S1",
                        "fastq_1": "/data/S1.R1.fastq.gz",
                        "fastq_2": "/data/S1.R2.fastq.gz",
                    }
                ],
            },
        )
        assert saved.status_code == 200, saved.text
        opened_hitrac = client.get(hitrac_url)
        assert opened_hitrac.json()["draft"]["workflow_samples"][0]["assay"] == "rna"
        validated = client.post(
            hitrac_url + "/validate",
            json={"expected_revision": saved.json()["draft_revision"]},
        )
        assert validated.status_code == 409, validated.text
        assert validated.json()["issues"][0]["code"] == "REFERENCE_PROFILE_REQUIRED"


def test_registered_adapter_can_materialize_other_input_columns(tmp_path, monkeypatch):
    class FilePlanAdapter(HiTracPreprocessAdapter):
        def supports_workflow_sample_plan(self):
            return True

        def normalize_workflow_sample_plans(self, plans):
            return [dict(plan) for plan in plans]

        def describe_workflow_sample_inputs(self, plans, *, sample_names):
            return tuple(
                InputFileRevisionSelection(
                    "library:" + sample_names[plan["sample_revision_id"]],
                    0,
                    tuple(plan["files"]),
                )
                for plan in plans
            )

        def materialize_workflow_sample_rows(self, plans, *, sample_names, input_paths):
            return Result.success(
                [
                    {
                        "sample_id": sample_names[plan["sample_revision_id"]],
                        "fastq_1": input_paths[plan["files"][0]],
                        "fastq_2": input_paths[plan["files"][1]],
                    }
                    for plan in plans
                ]
            )

    app = create_app(registry=WorkflowRegistry([FilePlanAdapter()]))
    seed_test_authentication(app)
    root, project, sample, reads = _managed_reads(app, tmp_path)
    validated_service = app.state.validated_input_service
    validate = Mock(wraps=validated_service.validate)
    monkeypatch.setattr(validated_service, "validate", validate)
    plan = {"sample_revision_id": sample, "files": [reads["r1"], reads["r2"]]}
    with ApiTestClient(app) as client:
        url = _session_url(client, project, HITRAC_WORKFLOW_ID)
        saved = client.request(
            "PUT", url + "/draft", json=_draft(0, sample, None, [plan])
        )
        assert saved.status_code == 200, saved.text
        opened = client.get(url)
        assert opened.json()["draft"]["workflow_samples"] == [plan]
        validated = client.post(
            url + "/validate",
            json={"expected_revision": saved.json()["draft_revision"]},
        )
        assert validated.status_code == 409, validated.text
        assert validated.json()["issues"][0]["code"] == "REFERENCE_PROFILE_REQUIRED"
        assert str(root) not in opened.text + validated.text
        validate.assert_called_once()
        assert validate.call_args.args[1].samples == [
            {
                "sample_id": "sample-a",
                "fastq_1": str(root / "reads" / "s1.fastq.gz"),
                "fastq_2": str(root / "reads" / "s2.fastq.gz"),
            }
        ]
        assert validate.call_args.kwargs["input_file_revision_selections"] == (
            InputFileRevisionSelection(
                "library:sample-a", 0, (reads["r1"], reads["r2"])
            ),
        )
        other_project, other_sample = _project_sample(app)
        other_url = _session_url(client, other_project, HITRAC_WORKFLOW_ID)
        foreign_file_plan = {**plan, "sample_revision_id": other_sample}
        saved = client.request(
            "PUT",
            other_url + "/draft",
            json=_draft(0, other_sample, None, [foreign_file_plan]),
        )
        assert saved.status_code == 200, saved.text
        rejected = client.post(
            other_url + "/validate",
            json={"expected_revision": saved.json()["draft_revision"]},
        )
        assert rejected.json()["issues"][0]["code"] == "WORKFLOW_SAMPLE_PLAN_INVALID"
        assert rejected.json()["snapshot"] is None
        assert str(root) not in rejected.text
        validate.assert_called_once()


@pytest.mark.parametrize("workflow_id", [HITRAC_WORKFLOW_ID, BULK_WORKFLOW_ID])
def test_workflows_without_plan_support_keep_legacy_read_order(
    reference_ready_app, workflow_id
):
    app = reference_ready_app
    project, sample = _project_sample(app)
    read1 = {"read": 1, "input_file_revision_id": "inpfr_" + "a" * 32}
    read2 = {"read": 2, "input_file_revision_id": "inpfr_" + "b" * 32}
    plan = {"sample_revision_id": sample, "reads": [read2, read1]}
    with ApiTestClient(app) as client:
        url = _session_url(client, project, workflow_id)
        saved = client.request(
            "PUT", url + "/draft", json=_draft(0, sample, None, [plan])
        )
        assert saved.status_code == 200, saved.text
        assert saved.json()["draft"]["workflow_samples"] == [
            {**plan, "reads": [read1, read2]}
        ]
        revision = saved.json()["draft_revision"]
        invalid = client.request(
            "PUT",
            url + "/draft",
            json=_draft(
                revision,
                sample,
                None,
                [{**plan, "reads": [read1, read1]}],
            ),
        )
        assert invalid.status_code == 422, invalid.text
        assert invalid.json()["reason_code"] == "INVALID_REQUEST"
        assert client.get(url).json()["draft_revision"] == revision
        checked = client.post(url + "/validate", json={"expected_revision": revision})
        assert checked.status_code == 409, checked.text
        assert checked.json()["issues"][0]["code"] == "REFERENCE_PROFILE_REQUIRED"


def test_saved_draft_with_an_unregistered_workflow_returns_its_public_error(
    reference_ready_app,
):
    app = reference_ready_app
    project, sample = _project_sample(app)
    with ApiTestClient(app) as client:
        url = _session_url(client, project, HITRAC_WORKFLOW_ID)
        saved = client.request("PUT", url + "/draft", json=_draft(0, sample, None, []))
        assert saved.status_code == 200, saved.text
    restarted = create_app(registry=WorkflowRegistry(()))
    seed_test_authentication(restarted)
    with ApiTestClient(restarted) as client:
        checked = client.post(
            url + "/validate",
            json={"expected_revision": saved.json()["draft_revision"]},
        )
        assert checked.status_code == 404, checked.text
        assert checked.json()["reason_code"] == "WORKFLOW_UNAVAILABLE"


class _DecliningSamplePlan:
    def supports_workflow_sample_plan(self) -> bool:
        return False

    def normalize_workflow_sample_plans(self, plans):
        raise AssertionError(plans)

    def describe_workflow_sample_inputs(self, plans, *, sample_names):
        raise AssertionError((plans, sample_names))

    def materialize_workflow_sample_rows(self, plans, *, sample_names, input_paths):
        raise AssertionError((plans, sample_names, input_paths))


def _managed_reads(app, tmp_path: Path):
    from encode_pipeline.services.input_registry import (
        InputRegistryAdminService,
        InputRegistryService,
    )
    from encode_pipeline.services.private_storage_pools import (
        PRIVATE_STORAGE_POOL_SCHEMA_VERSION,
        load_private_storage_pool_config,
    )

    root = (tmp_path / "pool").resolve()
    (root / "reads").mkdir(parents=True)
    (root / "reads" / "s1.fastq.gz").write_bytes(b"@read\nACGT\n+\n!!!!\n")
    (root / "reads" / "s2.fastq.gz").write_bytes(b"@read\nTGCA\n+\n!!!!\n")
    config_path = tmp_path / "pools.json"
    config_path.write_text(
        json.dumps(
            {
                "schema_version": PRIVATE_STORAGE_POOL_SCHEMA_VERSION,
                "storage_pool_roots": {"authoring-inputs": str(root)},
            }
        ),
        encoding="utf-8",
    )
    app.state.authoring_service._storage_pool_config_path = config_path
    app.state.validated_input_service._storage_pool_config_path = config_path
    project, sample_revision_id = _project_sample(app)
    admin = InputRegistryAdminService(
        InputRegistryService(
            repository=app.state.persistence.input_registry_repository
        ),
        private_config=load_private_storage_pool_config(config_path),
    )
    pool = admin.register_storage_pool(
        display_name="Authoring inputs",
        config_key="authoring-inputs",
    )
    admin.bind_project_storage_pool(
        project_id=project.project_id,
        storage_pool_id=pool.storage_pool_id,
    )
    registered = {}
    for name, relative in (("r1", "reads/s1.fastq.gz"), ("r2", "reads/s2.fastq.gz")):
        created = admin.register_input_file(
            project_id=project.project_id,
            stable_key=f"sample-a-{name}",
            pool_relative_path=relative,
        )
        registered[name] = created["input_file_revision_id"]
    return root, project, sample_revision_id, registered


def _encode_plan(sample_revision_id: str, reads: list[dict[str, object]], **overrides):
    plan = {
        "sample_revision_id": sample_revision_id,
        "assay": "chipseq",
        "target": "CTCF",
        "peak_mode": "narrow",
        "layout": "PE",
        "role": "treatment",
        "control_sample_revision_id": None,
        "reads": reads,
    }
    plan.update(overrides)
    return plan


def _draft(
    expected_revision: int,
    sample_revision_id: str,
    reference_profile_revision_id: str | None,
    workflow_samples: list[dict[str, object]],
    *,
    input_selections: list[dict[str, object]] | None = None,
) -> dict[str, object]:
    return {
        "expected_revision": expected_revision,
        "config": {},
        "samples": None,
        "sample_revision_ids": [sample_revision_id],
        "reference_profile_revision_id": reference_profile_revision_id,
        "workflow_samples": workflow_samples,
        "input_selections": input_selections or [],
    }


def _project_sample(app):
    registry = DataRegistryService(
        repository=app.state.persistence.data_registry_repository
    )
    project = registry.create_project("Sample plan")
    imported = registry.import_samples(
        project.project_id,
        rows=(
            SampleImportSpec(
                stable_key="sample-a", display_name="Sample A", attributes={}
            ),
        ),
    )
    return project, imported[0].revision.sample_revision_id


def _session_url(client, project, workflow_id) -> str:
    response = client.post(
        "/api/v1/agent/authoring-sessions",
        json={"project_id": project.project_id, "workflow_id": workflow_id},
    )
    assert response.status_code == 200, response.text
    return "/api/v1/agent/authoring-sessions/" + response.json()["sid"]
