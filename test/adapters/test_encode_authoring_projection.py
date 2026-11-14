"""ENCODE authoring projection exposes only approved resource fields."""

from encode_pipeline.adapters.encode_authoring import build_encode_authoring_schema
from encode_pipeline.adapters.encode_authoring_projection import project_encode_fields
from encode_pipeline.adapters.hitrac_preprocess.authoring import (
    build_hitrac_authoring_schema,
)
from encode_pipeline.platform.registry import ENCODE_WORKFLOW_ID


def test_encode_projection_keeps_resource_constraints_and_omits_schema_defaults():
    schema = build_encode_authoring_schema().to_dict()
    fields = project_encode_fields(ENCODE_WORKFLOW_ID, schema)

    assert fields is not None
    by_token = {field.token: field for field in fields}
    assert set(by_token) == {"config.threads", "options.cores"}
    assert by_token["config.threads"].constraint.minimum == 1
    assert by_token["config.threads"].constraint.maximum is None
    assert by_token["options.cores"].constraint.maximum == 1024
    rendered = str([field.public_dict() for field in fields])
    assert "default" not in rendered
    assert "description" not in rendered
    assert "8" not in rendered


def test_workflow_without_a_policy_has_no_projection():
    schema = build_hitrac_authoring_schema().to_dict()

    assert project_encode_fields("hitrac-preprocess", schema) is None
