"""Tests for workflow-neutral input Bundle value objects."""

from __future__ import annotations

from dataclasses import FrozenInstanceError, fields, replace
import hashlib
from pathlib import Path

import pytest

from encode_pipeline.platform.adapters import WorkflowInputs
from encode_pipeline.platform.input_bundles import (
    ImportedWorkflowInputs,
    InputBundleArtifact,
    InputBundleFile,
    InputBundleFileSetAlternatives,
    InputBundleFileObservation,
    InputBundleIdentity,
    InputBundleMapping,
    WorkflowInputBundle,
    validate_input_bundle_relative_path,
)


SHA256 = "1" * 64


def _identity() -> InputBundleIdentity:
    return InputBundleIdentity(
        contract_id="urn:omics-intake:schema:intake-bundle:0.2",
        contract_version="0.2",
        schema_sha256=SHA256,
        bundle_sha256="2" * 64,
        producer_name="omics-intake",
        producer_version="0.2.0",
        canonical_project_identity="sha256:" + "3" * 64,
    )


def _artifact() -> InputBundleArtifact:
    content = b"sample\n"
    return InputBundleArtifact(
        role="sample_sheet",
        relative_path="samples.encode.tsv",
        media_type="text/tab-separated-values",
        size_bytes=len(content),
        sha256=hashlib.sha256(content).hexdigest(),
        content=content,
    )


def _file() -> InputBundleFile:
    return InputBundleFile(
        file_id="ena:SRR000001:1",
        relative_path="downloads/ena/SRR000001.fastq.gz",
        file_format="fastq",
        role="treatment",
        read_number=1,
        size_bytes=8,
        sha256="4" * 64,
    )


def test_identity_is_frozen_and_serializes_without_local_paths() -> None:
    identity = _identity()
    serialized = identity.to_dict()
    serialized["producer_name"] = "changed"

    with pytest.raises(FrozenInstanceError):
        identity.producer_name = "changed"

    assert identity.producer_name == "omics-intake"
    assert set(identity.to_dict()) == {
        "contract_id",
        "contract_version",
        "schema_sha256",
        "bundle_sha256",
        "producer_name",
        "producer_version",
        "canonical_project_identity",
    }


def test_bundle_copies_collections_and_only_resolves_valid_relative_paths(
    tmp_path: Path,
) -> None:
    artifacts = [_artifact()]
    files = [_file()]
    bundle = WorkflowInputBundle(
        identity=_identity(),
        workflow_name="encode-epigenomics",
        genome="mm10",
        render_contract="encode-render-v3",
        project_root=tmp_path.resolve(),
        artifacts=artifacts,
        files=files,
    )
    artifacts.clear()
    files.clear()

    assert bundle.artifact("sample_sheet") == _artifact()
    assert bundle.file_for_path("downloads/ena/SRR000001.fastq.gz") == _file()
    assert bundle.project_path("downloads/ena/SRR000001.fastq.gz") == (
        tmp_path.resolve() / "downloads/ena/SRR000001.fastq.gz"
    )
    with pytest.raises(ValueError):
        bundle.project_path("../outside.fastq.gz")


@pytest.mark.parametrize(
    "value",
    (
        "",
        "/absolute",
        "../outside",
        "nested/../outside",
        "nested//file",
        "nested\\file",
        ".hidden",
        "C:/private/file",
        "contains\x00nul",
        "contains\nnewline",
    ),
)
def test_relative_path_contract_rejects_non_portable_values(value: str) -> None:
    with pytest.raises(ValueError):
        validate_input_bundle_relative_path(value)


def test_mapping_and_import_result_copy_input_values() -> None:
    inputs = WorkflowInputs(
        config={"threads": 4},
        samples=[{"sample": "S1"}],
        options={"strict_inputs": True},
    )
    mapping = InputBundleMapping(
        inputs=inputs,
        required_project_files=("downloads/S1.fastq.gz",),
    )
    imported = ImportedWorkflowInputs(
        identity=_identity(),
        workflow_id="fake",
        inputs=WorkflowInputs(**mapping.inputs.to_dict()),
        source_files=(
            InputBundleFileObservation(
                relative_path="downloads/S1.fastq.gz",
                size_bytes=8,
                sha256="5" * 64,
                contract_bound=True,
            ),
        ),
    )
    serialized = imported.inputs.to_dict()
    serialized["config"]["threads"] = 8

    assert imported.inputs.config == {"threads": 4}
    assert imported.source_files[0].contract_bound is True


def test_file_set_alternatives_preserve_explicit_preference_without_io() -> None:
    large = ("indexes/genome.1.bt2l", "indexes/genome.2.bt2l")
    standard = ("indexes/genome.1.bt2", "indexes/genome.2.bt2")
    alternatives = InputBundleFileSetAlternatives(alternatives=(large, standard))
    mapping = InputBundleMapping(
        inputs=WorkflowInputs(config={}),
        required_project_files=("downloads/S1.fastq.gz",),
        required_project_file_sets=(alternatives,),
    )

    assert mapping.required_project_file_sets[0].alternatives == (large, standard)
    assert tuple(field.name for field in fields(InputBundleFileSetAlternatives)) == (
        "alternatives",
    )
    assert tuple(field.name for field in fields(InputBundleMapping)) == (
        "inputs",
        "required_project_files",
        "required_project_file_sets",
    )
    assert tuple(field.name for field in fields(InputBundleFileObservation)) == (
        "relative_path",
        "size_bytes",
        "sha256",
        "contract_bound",
    )


@pytest.mark.parametrize(
    "alternatives",
    (
        (),
        "indexes/a.bt2",
        ("indexes/a.bt2",),
        object(),
        ((),),
        (("indexes/z.bt2", "indexes/a.bt2"),),
        (("../outside.bt2",),),
        (("indexes/a.bt2",), ("indexes/a.bt2",)),
    ),
)
def test_file_set_alternatives_reject_invalid_declarations(alternatives) -> None:
    with pytest.raises(ValueError):
        InputBundleFileSetAlternatives(alternatives=alternatives)


@pytest.mark.parametrize(
    "factory",
    (
        lambda: InputBundleIdentity(
            contract_id="contract",
            contract_version="0.2",
            schema_sha256="A" * 64,
            bundle_sha256="2" * 64,
            producer_name="producer",
            producer_version="0.2.0",
            canonical_project_identity="sha256:" + "3" * 64,
        ),
        lambda: InputBundleArtifact(
            role="sample_sheet",
            relative_path="samples.tsv",
            media_type="text/tab-separated-values",
            size_bytes=1,
            sha256=hashlib.sha256(b"wrong").hexdigest(),
            content=b"x",
        ),
        lambda: InputBundleMapping(
            inputs=WorkflowInputs(config={}),
            required_project_files=(),
        ),
    ),
)
def test_bundle_values_reject_invalid_integrity_or_empty_mapping(factory) -> None:
    with pytest.raises(ValueError):
        factory()


@pytest.mark.parametrize(
    "factory",
    (
        lambda: replace(_identity(), producer_name=" producer"),
        lambda: replace(_identity(), canonical_project_identity="3" * 64),
        lambda: replace(_artifact(), role=""),
        lambda: replace(_artifact(), size_bytes=True),
        lambda: replace(_artifact(), sha256="A" * 64),
        lambda: replace(_artifact(), content="sample\n"),
        lambda: replace(_artifact(), size_bytes=1),
        lambda: replace(_file(), file_id=""),
        lambda: replace(_file(), read_number=3),
        lambda: replace(_file(), size_bytes=-1),
        lambda: replace(_file(), sha256="A" * 64),
        lambda: InputBundleFileObservation(
            relative_path="downloads/S1.fastq.gz",
            size_bytes=-1,
            sha256="5" * 64,
            contract_bound=True,
        ),
        lambda: InputBundleFileObservation(
            relative_path="downloads/S1.fastq.gz",
            size_bytes=1,
            sha256="A" * 64,
            contract_bound=True,
        ),
        lambda: InputBundleFileObservation(
            relative_path="downloads/S1.fastq.gz",
            size_bytes=1,
            sha256="5" * 64,
            contract_bound=1,
        ),
    ),
)
def test_bundle_value_fields_fail_closed(factory) -> None:
    with pytest.raises(ValueError):
        factory()


def test_bundle_collection_and_result_invariants_fail_closed(tmp_path: Path) -> None:
    common = {
        "identity": _identity(),
        "workflow_name": "encode-epigenomics",
        "genome": "mm10",
        "render_contract": "encode-render-v3",
        "project_root": tmp_path.resolve(),
        "artifacts": (_artifact(),),
        "files": (_file(),),
    }
    invalid_bundles = (
        {**common, "identity": object()},
        {**common, "workflow_name": ""},
        {**common, "project_root": Path("relative")},
        {**common, "artifacts": (object(),)},
        {**common, "files": (object(),)},
        {**common, "artifacts": (_artifact(), _artifact())},
        {
            **common,
            "files": (
                replace(_file(), relative_path="downloads/z.fastq.gz"),
                replace(_file(), file_id="ena:SRR000002:1"),
            ),
        },
    )
    for values in invalid_bundles:
        with pytest.raises(ValueError):
            WorkflowInputBundle(**values)

    with pytest.raises(ValueError):
        InputBundleMapping(
            inputs=object(),
            required_project_files=("downloads/S1.fastq.gz",),
        )
    with pytest.raises(ValueError):
        InputBundleMapping(
            inputs=WorkflowInputs(config={}),
            required_project_files=(
                "downloads/S1.fastq.gz",
                "downloads/S1.fastq.gz",
            ),
        )
    alternatives = InputBundleFileSetAlternatives(
        alternatives=(("downloads/S1.fastq.gz",),)
    )
    with pytest.raises(ValueError):
        InputBundleMapping(
            inputs=WorkflowInputs(config={}),
            required_project_files=("downloads/S1.fastq.gz",),
            required_project_file_sets=(alternatives,),
        )
    with pytest.raises(ValueError):
        InputBundleMapping(
            inputs=WorkflowInputs(config={}),
            required_project_files=(),
            required_project_file_sets=(object(),),
        )

    observation = InputBundleFileObservation(
        relative_path="downloads/S1.fastq.gz",
        size_bytes=1,
        sha256="5" * 64,
        contract_bound=True,
    )
    invalid_results = (
        {"identity": object(), "workflow_id": "fake", "inputs": WorkflowInputs({})},
        {"identity": _identity(), "workflow_id": "", "inputs": WorkflowInputs({})},
        {"identity": _identity(), "workflow_id": "fake", "inputs": object()},
    )
    for values in invalid_results:
        with pytest.raises(ValueError):
            ImportedWorkflowInputs(source_files=(observation,), **values)
    with pytest.raises(ValueError):
        ImportedWorkflowInputs(
            identity=_identity(),
            workflow_id="fake",
            inputs=WorkflowInputs({}),
            source_files=(object(),),
        )
    with pytest.raises(ValueError):
        ImportedWorkflowInputs(
            identity=_identity(),
            workflow_id="fake",
            inputs=WorkflowInputs({}),
            source_files=(observation, observation),
        )


def test_bundle_lookup_rejects_missing_or_ambiguous_roles_and_paths(
    tmp_path: Path,
) -> None:
    first = _file()
    second = replace(first, file_id="ena:SRR000002:1")
    bundle = WorkflowInputBundle(
        identity=_identity(),
        workflow_name="encode-epigenomics",
        genome="mm10",
        render_contract="encode-render-v3",
        project_root=tmp_path.resolve(),
        artifacts=(_artifact(),),
        files=(first, second),
    )

    with pytest.raises(KeyError):
        bundle.artifact("workflow_config")
    with pytest.raises(ValueError):
        bundle.file_for_path(first.relative_path)
