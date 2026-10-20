"""Reference-profile binding tests for the ENCODE-style adapter."""

from __future__ import annotations

from copy import deepcopy
import hashlib
import os
from pathlib import Path

import pytest

import encode_pipeline.adapters.bulk_rnaseq.resource_closure as resource_closure
from encode_pipeline.adapters.encode import EncodeStyleWorkflowAdapter
from encode_pipeline.platform.adapters import (
    ReferenceProfileBindingAdapter,
    WorkflowInputs,
)


def _sha256(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def _write(path: Path, value: bytes) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(value)
    return _sha256(value)


def _binding(root: Path) -> dict[str, object]:
    resources: dict[str, object] = {}
    for name, content in (
        ("reference_fasta", b">chr1\nACGT\n"),
        ("gtf", b'chr1\ttest\texon\t1\t4\t.\t+\t.\tgene_id "g1";\n'),
        ("chrom_sizes", b"chr1\t4\n"),
        ("blacklist", b"chr1\t1\t2\n"),
    ):
        suffix = {
            "reference_fasta": "fa",
            "gtf": "gtf",
            "chrom_sizes": "sizes",
            "blacklist": "bed",
        }[name]
        path = root / f"reference/{name}.{suffix}"
        resources[name] = {"path": str(path), "sha256": _write(path, content)}

    prefix = root / "reference/bowtie2/GRCh38"
    files: dict[str, str] = {}
    for suffix in (".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"):
        path = Path(f"{prefix}{suffix}")
        files[suffix] = _write(path, suffix.encode("ascii"))
    return {
        "schema_version": "encode-reference-binding-v1",
        "assembly": "GRCh38",
        "effective_genome_size": "hs",
        "genome_resources": resources,
        "bowtie2_index": {"prefix": str(prefix), "files": files},
    }


def _binding_v2(root: Path) -> dict[str, object]:
    payload = _binding(root)
    payload["schema_version"] = "encode-reference-binding-v2"
    resources = payload["genome_resources"]
    assert isinstance(resources, dict)
    fasta = Path(resources["reference_fasta"]["path"])
    fai = Path(str(fasta) + ".fai")
    resources["reference_fai"] = {
        "path": str(fai),
        "sha256": _write(fai, b"chr1\t4\t6\t4\t5\n"),
    }
    reference_dict = root / "reference/reference.dict"
    resources["reference_dict"] = {
        "path": str(reference_dict),
        "sha256": _write(reference_dict, b"@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:4\n"),
    }
    payload["organelle_chroms"] = ["chr1"]
    payload["canonical_chroms"] = ["chr1"]
    return payload


def _inputs(root: Path) -> WorkflowInputs:
    return WorkflowInputs(
        config={"threads": 1},
        samples=[
            {
                "sample": "S1",
                "fastq_1": str(root / "S1_R1.fastq.gz"),
                "fastq_2": str(root / "S1_R2.fastq.gz"),
                "layout": "PE",
                "assay": "chipseq",
                "target": "CTCF",
                "peak_mode": "narrow",
            }
        ],
        options={},
    )


def test_public_authoring_schema_hides_operator_reference_paths():
    schema = EncodeStyleWorkflowAdapter().schema().to_dict()

    assert "genome_resources" not in schema["config_schema"]["properties"]
    items = schema["sample_schema"]["items"]
    assert "bowtie2_index" not in items["properties"]
    assert "genome" not in items["properties"]
    assert "bowtie2_index" not in items["required"]
    assert "genome" not in items["required"]


def test_encode_binding_verifies_and_injects_exact_reference(tmp_path: Path):
    adapter = EncodeStyleWorkflowAdapter()
    payload = _binding(tmp_path)

    assert isinstance(adapter, ReferenceProfileBindingAdapter)
    verified = adapter.verify_reference_profile_binding(payload)
    bound = adapter.bind_reference_profile(_inputs(tmp_path), payload)

    assert verified.is_success
    assert bound.is_success
    assert bound.value.identity == verified.value
    assert verified.value.workflow_id == adapter.metadata.workflow_id
    assert verified.value.contract_version == "encode-reference-binding-v1"
    assert len(verified.value.identity_sha256) == 64
    assert str(tmp_path) not in repr(verified.value)
    assert str(tmp_path) not in repr(bound.value)
    resources = bound.value.inputs.config["genome_resources"]
    assert list(resources) == ["GRCh38"]
    assert resources["GRCh38"]["reference_fasta"].endswith("reference_fasta.fa")
    row = bound.value.inputs.samples[0]
    assert row["genome"] == "GRCh38"
    assert row["bowtie2_index"].endswith("bowtie2/GRCh38")
    assert bound.value.adapter.validate(bound.value.inputs).is_success
    direct_plan = adapter.plan_workspace(bound.value.inputs, tmp_path / "workspace")
    selected_plan = bound.value.adapter.plan_workspace(
        bound.value.inputs, tmp_path / "workspace"
    )
    assert direct_plan.is_success and selected_plan.is_success
    assert selected_plan.value == direct_plan.value


def test_encode_binding_v2_checks_contig_closure_and_injects_organelles(
    tmp_path: Path,
):
    adapter = EncodeStyleWorkflowAdapter()
    payload = _binding_v2(tmp_path)

    verified = adapter.verify_reference_profile_binding(payload)
    bound = adapter.bind_reference_profile(_inputs(tmp_path), payload)

    assert verified.is_success, verified.issues
    assert bound.is_success, bound.issues
    assert verified.value.contract_version == "encode-reference-binding-v2"
    resources = bound.value.inputs.config["genome_resources"]["GRCh38"]
    assert resources["organelle_chroms"] == ["chr1"]
    assert resources["reference_fasta"].endswith("reference_fasta.fa")


@pytest.mark.parametrize(
    ("resource", "content"),
    [
        ("chrom_sizes", b"chr1\t5\n"),
        ("gtf", b'chr2\ttest\texon\t1\t4\t.\t+\t.\tgene_id "g1";\n'),
        ("blacklist", b"chr2\t1\t2\n"),
        ("reference_fai", b"chr1\t5\t0\t5\t6\n"),
        ("reference_dict", b"@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:5\n"),
    ],
)
def test_encode_binding_v2_rejects_cross_resource_mismatch(
    tmp_path: Path, resource: str, content: bytes
):
    adapter = EncodeStyleWorkflowAdapter()
    payload = _binding_v2(tmp_path)
    resources = payload["genome_resources"]
    assert isinstance(resources, dict)
    path = Path(resources[resource]["path"])
    resources[resource]["sha256"] = _write(path, content)

    result = adapter.verify_reference_profile_binding(payload)

    assert result.is_failure


def test_encode_binding_v2_rejects_organelle_and_egs_mismatch(tmp_path: Path):
    adapter = EncodeStyleWorkflowAdapter()
    payload = _binding_v2(tmp_path)
    payload["organelle_chroms"] = ["chr2"]
    assert adapter.verify_reference_profile_binding(payload).is_failure

    payload = _binding_v2(tmp_path / "egs")
    payload["effective_genome_size"] = "mm"
    assert adapter.verify_reference_profile_binding(payload).is_failure


def test_encode_binding_v2_identity_is_path_free(tmp_path: Path):
    adapter = EncodeStyleWorkflowAdapter()
    first = adapter.verify_reference_profile_binding(_binding_v2(tmp_path / "a"))
    second = adapter.verify_reference_profile_binding(_binding_v2(tmp_path / "b"))

    assert first.is_success and second.is_success
    assert first.value == second.value


def test_encode_binding_identity_is_path_free_and_content_deterministic(
    tmp_path: Path,
):
    adapter = EncodeStyleWorkflowAdapter()
    first = adapter.verify_reference_profile_binding(_binding(tmp_path / "a"))
    second = adapter.verify_reference_profile_binding(_binding(tmp_path / "b"))

    assert first.is_success and second.is_success
    assert first.value == second.value


@pytest.mark.parametrize("caller_field", ["genome_resources", "sample_genome", "index"])
def test_encode_binding_rejects_caller_reference_conflicts(
    tmp_path: Path,
    caller_field: str,
):
    adapter = EncodeStyleWorkflowAdapter()
    payload = _binding(tmp_path)
    inputs = _inputs(tmp_path)
    config = deepcopy(dict(inputs.config))
    samples = deepcopy(inputs.samples)
    if caller_field == "genome_resources":
        config["genome_resources"] = {"mm10": {"effective_genome_size": "mm"}}
    elif caller_field == "sample_genome":
        samples[0]["genome"] = "mm10"
    else:
        samples[0]["bowtie2_index"] = "/operator/bypass"

    result = adapter.bind_reference_profile(
        WorkflowInputs(config=config, samples=samples), payload
    )

    assert result.is_failure
    assert result.errors[0].code == "ENCODE_REFERENCE_BINDING_CONFLICT"
    assert str(tmp_path) not in str(result.errors[0].to_dict())


def test_encode_binding_rejects_exact_private_values_from_caller(tmp_path: Path):
    adapter = EncodeStyleWorkflowAdapter()
    payload = _binding(tmp_path)
    inputs = _inputs(tmp_path)
    resources = payload["genome_resources"]
    config = {
        **inputs.config,
        "genome_resources": {
            payload["assembly"]: {
                "effective_genome_size": payload["effective_genome_size"],
                **{name: descriptor["path"] for name, descriptor in resources.items()},
            }
        },
    }

    result = adapter.bind_reference_profile(
        WorkflowInputs(config=config, samples=inputs.samples), payload
    )

    assert result.is_failure
    assert result.errors[0].code == "ENCODE_REFERENCE_BINDING_CONFLICT"
    assert str(tmp_path) not in str(result.errors[0].to_dict())


def test_encode_binding_rejects_mixed_assembly(tmp_path: Path):
    payload = _binding(tmp_path)
    inputs = _inputs(tmp_path)
    rows = deepcopy(inputs.samples)
    rows.append({**rows[0], "sample": "S2", "genome": "mm10"})

    result = EncodeStyleWorkflowAdapter().bind_reference_profile(
        WorkflowInputs(config=inputs.config, samples=rows), payload
    )

    assert result.is_failure
    assert result.errors[0].code == "ENCODE_REFERENCE_BINDING_CONFLICT"


def test_encode_binding_rejects_digest_mismatch_symlink_fifo_and_extra_index(
    tmp_path: Path,
):
    adapter = EncodeStyleWorkflowAdapter()

    mismatch = _binding(tmp_path / "mismatch")
    fasta = Path(mismatch["genome_resources"]["reference_fasta"]["path"])
    fasta.write_bytes(b"changed")
    assert adapter.verify_reference_profile_binding(mismatch).is_failure

    symlinked = _binding(tmp_path / "symlink")
    gtf = Path(symlinked["genome_resources"]["gtf"]["path"])
    original = gtf.with_suffix(".original")
    gtf.rename(original)
    gtf.symlink_to(original)
    assert adapter.verify_reference_profile_binding(symlinked).is_failure

    fifo_payload = _binding(tmp_path / "fifo")
    blacklist = Path(fifo_payload["genome_resources"]["blacklist"]["path"])
    blacklist.unlink()
    os.mkfifo(blacklist)
    try:
        assert adapter.verify_reference_profile_binding(fifo_payload).is_failure
    finally:
        blacklist.unlink()

    extra = _binding(tmp_path / "extra")
    prefix = Path(extra["bowtie2_index"]["prefix"])
    Path(f"{prefix}.extra.bt2").write_bytes(b"extra")
    result = adapter.verify_reference_profile_binding(extra)
    assert result.is_failure
    assert result.errors[0].code == "ENCODE_REFERENCE_BINDING_INVALID"


def test_encode_binding_rejects_unknown_payload_fields(tmp_path: Path):
    payload = _binding(tmp_path)
    payload["operator_note"] = "not contracted"

    result = EncodeStyleWorkflowAdapter().verify_reference_profile_binding(payload)

    assert result.is_failure
    assert result.errors[0].code == "ENCODE_REFERENCE_BINDING_INVALID"


def test_encode_binding_rejects_descriptor_replacement_race(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    payload = _binding(tmp_path)
    target = Path(payload["genome_resources"]["reference_fasta"]["path"])
    original = resource_closure._DescriptorRoot.reopen_identity
    replaced = False

    def replace_before_reopen(self, parts, *, is_directory):
        nonlocal replaced
        if parts == (target.name,) and not replaced:
            replacement = target.with_suffix(".replacement")
            replacement.write_bytes(b">changed\nTGCA\n")
            os.replace(replacement, target)
            replaced = True
        return original(self, parts, is_directory=is_directory)

    monkeypatch.setattr(
        resource_closure._DescriptorRoot,
        "reopen_identity",
        replace_before_reopen,
    )

    result = EncodeStyleWorkflowAdapter().verify_reference_profile_binding(payload)

    assert replaced is True
    assert result.is_failure
    assert result.errors[0].code == "ENCODE_REFERENCE_BINDING_INVALID"


def test_encode_binding_v2_rejects_replacement_between_hash_and_parse(
    tmp_path: Path, monkeypatch
):
    import encode_pipeline.adapters.encode_references as refs

    adapter = EncodeStyleWorkflowAdapter()
    payload = _binding_v2(tmp_path)
    resources = payload["genome_resources"]
    assert isinstance(resources, dict)
    fasta = Path(resources["reference_fasta"]["path"])
    original = refs._iter_verified_lines
    replaced = False

    def hooked(path, expected_sha256):
        nonlocal replaced
        if not replaced and str(path) == str(fasta):
            replaced = True
            fasta.write_bytes(b">chr1\nTTTT\n")
        return original(path, expected_sha256)

    monkeypatch.setattr(refs, "_iter_verified_lines", hooked)
    result = adapter.verify_reference_profile_binding(payload)

    assert result.is_failure


def test_encode_binding_v2_does_not_skip_track_like_data_contigs(
    tmp_path: Path,
):
    adapter = EncodeStyleWorkflowAdapter()
    payload = _binding_v2(tmp_path)
    resources = payload["genome_resources"]
    assert isinstance(resources, dict)
    gtf = Path(resources["gtf"]["path"])
    resources["gtf"]["sha256"] = _write(
        gtf,
        b'track_decoy\ttest\texon\t1\t4\t.\t+\t.\tgene_id "g1";\n',
    )
    assert adapter.verify_reference_profile_binding(payload).is_failure

    payload = _binding_v2(tmp_path / "bed")
    resources = payload["genome_resources"]
    assert isinstance(resources, dict)
    blacklist = Path(resources["blacklist"]["path"])
    resources["blacklist"]["sha256"] = _write(blacklist, b"browser_decoy\t1\t2\n")
    assert adapter.verify_reference_profile_binding(payload).is_failure


def test_encode_binding_v2_accepts_real_bed_headers(tmp_path: Path):
    adapter = EncodeStyleWorkflowAdapter()
    payload = _binding_v2(tmp_path)
    resources = payload["genome_resources"]
    assert isinstance(resources, dict)
    blacklist = Path(resources["blacklist"]["path"])
    resources["blacklist"]["sha256"] = _write(
        blacklist,
        b"#chrom\tstart\tend\n"
        b"track name=blacklist\n"
        b"browser position chr1:1-2\n"
        b"chr1\t1\t2\n",
    )
    gtf = Path(resources["gtf"]["path"])
    resources["gtf"]["sha256"] = _write(
        gtf,
        b'# annotation\nchr1\ttest\texon\t1\t4\t.\t+\t.\tgene_id "g1";\n',
    )

    assert adapter.verify_reference_profile_binding(payload).is_success
