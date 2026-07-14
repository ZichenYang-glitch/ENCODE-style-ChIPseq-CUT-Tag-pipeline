"""Contracts for the artifact catalog seam imported by Snakemake code."""

from dataclasses import asdict
from pathlib import Path

from encode_pipeline import artifacts as package_artifacts
from workflow.lib import artifact as workflow_artifacts


REPO_ROOT = Path(__file__).resolve().parents[2]
INVENTORY = REPO_ROOT / "docs" / "architecture" / "artifact-inventory.yaml"
WORKFLOW_EXPORTS = {
    "Artifact",
    "VALID_ASSAY_GATES",
    "VALID_LEVELS",
    "VALID_SCOPES",
    "artifacts_by_id",
    "artifacts_by_manifest_output_type",
    "filter_artifacts",
    "load_artifacts",
    "validate_artifact",
}


def test_workflow_artifact_seam_reexports_package_catalog_objects():
    assert set(workflow_artifacts.__all__) == WORKFLOW_EXPORTS
    for name in WORKFLOW_EXPORTS:
        assert getattr(workflow_artifacts, name) is getattr(package_artifacts, name)


def test_workflow_artifact_seam_loads_and_indexes_the_real_inventory():
    artifacts = workflow_artifacts.load_artifacts(str(INVENTORY))

    assert artifacts
    assert all(
        isinstance(artifact, workflow_artifacts.Artifact) for artifact in artifacts
    )
    assert all(
        workflow_artifacts.validate_artifact(asdict(artifact)) == []
        for artifact in artifacts
    )
    assert len(workflow_artifacts.artifacts_by_id(artifacts)) == len(artifacts)
