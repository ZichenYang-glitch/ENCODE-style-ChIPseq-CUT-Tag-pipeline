from __future__ import annotations

import copy

import pytest
from jsonschema import Draft202012Validator

from encode_pipeline.contracts.deployment import load_schema
from encode_pipeline.deployment.models import DeploymentState, PLATFORM

from .support import manifest_for


def test_bundle_manifest_matches_the_packaged_machine_schema() -> None:
    manifest, _payload = manifest_for(PLATFORM)
    validator = Draft202012Validator(load_schema("deployment-bundle-v1.schema.json"))

    assert list(validator.iter_errors(manifest.to_dict())) == []


def test_bundle_machine_schema_accepts_jvm_inner_class_paths() -> None:
    manifest, _payload = manifest_for(PLATFORM)
    document = copy.deepcopy(manifest.to_dict())
    document["files"][0]["path"] = "payload/runtime/plugins/Validation$_closure.class"
    validator = Draft202012Validator(load_schema("deployment-bundle-v1.schema.json"))

    assert list(validator.iter_errors(document)) == []


@pytest.mark.parametrize(
    ("record", "path"),
    (
        ("contract", "payload/not-a-contract.json"),
        ("contract", "payload/contracts/../escape.json"),
        ("file", "payload/../escape"),
    ),
)
def test_bundle_machine_schema_rejects_paths_rejected_by_the_parser(
    record: str,
    path: str,
) -> None:
    manifest, _payload = manifest_for(PLATFORM)
    document = copy.deepcopy(manifest.to_dict())
    if record == "contract":
        document["contracts"][0]["path"] = path
    else:
        document["files"][0]["path"] = path
    validator = Draft202012Validator(load_schema("deployment-bundle-v1.schema.json"))

    assert list(validator.iter_errors(document))


def test_deployment_state_matches_the_packaged_machine_schema() -> None:
    state = DeploymentState.initial()
    validator = Draft202012Validator(load_schema("deployment-state-v1.schema.json"))

    assert list(validator.iter_errors(state.to_dict())) == []
