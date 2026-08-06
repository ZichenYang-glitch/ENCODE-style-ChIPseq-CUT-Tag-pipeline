"""Reference-profile binding tests for the bulk RNA-seq adapter."""

from __future__ import annotations

from copy import deepcopy
import hashlib
from pathlib import Path

from encode_pipeline.adapters.bulk_rnaseq import (
    BulkRnaSeqExecutionBinding,
    BulkRnaSeqResultsWorkflowAdapter,
    BulkRnaSeqTranscriptomeBinding,
    RuntimeAssetBinding,
)
from encode_pipeline.platform.adapters import (
    ReferenceProfileBindingAdapter,
    WorkflowInputs,
)


def _write(path: Path, value: bytes) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(value)
    return hashlib.sha256(value).hexdigest()


def _payload(root: Path) -> dict[str, object]:
    fasta = root / "reference/genome.fa"
    gtf = root / "reference/genes.gtf"
    transcript = root / "reference/transcripts.fa"
    fasta_sha256 = _write(fasta, b">chr1\nACGT\n")
    gtf_sha256 = _write(gtf, b'chr1\ttest\texon\t1\t4\t.\t+\t.\tgene_id "g1";\n')
    transcript_sha256 = _write(transcript, b">tx1\nACGT\n")
    return {
        "schema_version": "bulk-rnaseq-reference-binding-v1",
        "reference": {
            "reference_id": "tiny-ref",
            "fasta": str(fasta),
            "fasta_sha256": fasta_sha256,
            "gtf": str(gtf),
            "gtf_sha256": gtf_sha256,
            "annotation_style": "ensembl",
        },
        "transcriptome": {
            "reference_id": "tiny-ref",
            "fasta_sha256": fasta_sha256,
            "gtf_sha256": gtf_sha256,
            "transcript_fasta": str(transcript),
            "transcript_fasta_sha256": transcript_sha256,
        },
    }


def _adapter(tmp_path: Path, bulk_rnaseq_qualifications):
    old_transcript = tmp_path / "old/transcripts.fa"
    old_sha = _write(old_transcript, b">old\nAAAA\n")
    execution = BulkRnaSeqExecutionBinding(
        assets=RuntimeAssetBinding(root=(tmp_path / "runtime").resolve()),
        transcriptome=BulkRnaSeqTranscriptomeBinding(
            reference_id="old-ref",
            fasta_sha256="a" * 64,
            gtf_sha256="b" * 64,
            transcript_fasta=old_transcript.resolve(),
            transcript_fasta_sha256=old_sha,
        ),
        **bulk_rnaseq_qualifications,
    )
    return BulkRnaSeqResultsWorkflowAdapter(execution=execution)


def _reference_unbound_adapter(tmp_path: Path, bulk_rnaseq_qualifications):
    return BulkRnaSeqResultsWorkflowAdapter(
        execution=BulkRnaSeqExecutionBinding(
            assets=RuntimeAssetBinding(root=(tmp_path / "runtime").resolve()),
            transcriptome=None,
            **bulk_rnaseq_qualifications,
        )
    )


def _inputs(root: Path) -> WorkflowInputs:
    fastq_1 = root / "reads/S1.fastq.gz"
    _write(fastq_1, b"reads")
    return WorkflowInputs(
        config={"standard": {}},
        samples=[
            {
                "sample": "S1",
                "library": "lib1",
                "lane": "L001",
                "layout": "SE",
                "fastq_1": str(fastq_1),
                "strandedness": "reverse",
                "platform": "ILLUMINA",
            }
        ],
    )


def test_public_bulk_schema_hides_operator_reference_paths():
    from encode_pipeline.adapters.bulk_rnaseq import BulkRnaSeqWorkflowAdapter

    document = BulkRnaSeqWorkflowAdapter().schema().to_dict()
    standard = document["config_schema"]["properties"]["standard"]
    assert "reference" not in standard["properties"]
    assert "reference" not in standard["required"]


def test_bulk_binding_injects_reference_and_preserves_adapter_subclass(
    tmp_path: Path,
    bulk_rnaseq_qualifications,
):
    adapter = _adapter(tmp_path, bulk_rnaseq_qualifications)
    payload = _payload(tmp_path)

    assert isinstance(adapter, ReferenceProfileBindingAdapter)
    verified = adapter.verify_reference_profile_binding(payload)
    bound = adapter.bind_reference_profile(_inputs(tmp_path), payload)

    assert verified.is_success
    assert bound.is_success
    assert bound.value.identity == verified.value
    assert str(tmp_path) not in repr(bound.value)
    assert isinstance(bound.value.adapter, BulkRnaSeqResultsWorkflowAdapter)
    assert bound.value.adapter is not adapter
    assert bound.value.inputs.config["standard"]["reference"] == payload["reference"]
    assert (
        bound.value.adapter.execution_binding.assets is adapter.execution_binding.assets
    )
    assert (
        bound.value.adapter.execution_binding.implementation_qualification
        == adapter.execution_binding.implementation_qualification
    )
    assert bound.value.adapter.validate(bound.value.inputs).is_success


def test_bulk_binding_supplies_transcriptome_to_reference_unbound_runtime(
    tmp_path: Path,
    bulk_rnaseq_qualifications,
):
    adapter = _reference_unbound_adapter(tmp_path, bulk_rnaseq_qualifications)
    payload = _payload(tmp_path)

    bound = adapter.bind_reference_profile(_inputs(tmp_path), payload)

    assert bound.is_success
    assert adapter.execution_binding.transcriptome is None
    assert bound.value.adapter.execution_binding.transcriptome is not None
    assert (
        bound.value.adapter.execution_binding.transcriptome.reference_id
        == payload["transcriptome"]["reference_id"]
    )


def test_bulk_binding_identity_is_path_free_and_content_deterministic(
    tmp_path: Path,
    bulk_rnaseq_qualifications,
):
    adapter = _adapter(tmp_path, bulk_rnaseq_qualifications)
    first = adapter.verify_reference_profile_binding(_payload(tmp_path / "a"))
    second = adapter.verify_reference_profile_binding(_payload(tmp_path / "b"))

    assert first.is_success and second.is_success
    assert first.value == second.value
    assert str(tmp_path) not in repr(first.value)


def test_bulk_binding_rejects_caller_reference_conflict(
    tmp_path: Path,
    bulk_rnaseq_qualifications,
):
    adapter = _adapter(tmp_path, bulk_rnaseq_qualifications)
    inputs = _inputs(tmp_path)
    config = deepcopy(dict(inputs.config))
    config["standard"]["reference"] = {
        "reference_id": "caller-bypass",
        "fasta": "/caller/bypass.fa",
    }

    result = adapter.bind_reference_profile(
        WorkflowInputs(config=config, samples=inputs.samples), _payload(tmp_path)
    )

    assert result.is_failure
    assert result.errors[0].code == "BULK_RNASEQ_REFERENCE_BINDING_CONFLICT"
    assert str(tmp_path) not in str(result.errors[0].to_dict())


def test_bulk_binding_rejects_exact_private_reference_from_caller(
    tmp_path: Path,
    bulk_rnaseq_qualifications,
):
    adapter = _adapter(tmp_path, bulk_rnaseq_qualifications)
    payload = _payload(tmp_path)
    inputs = _inputs(tmp_path)
    config = deepcopy(dict(inputs.config))
    config["standard"]["reference"] = deepcopy(payload["reference"])

    result = adapter.bind_reference_profile(
        WorkflowInputs(config=config, samples=inputs.samples), payload
    )

    assert result.is_failure
    assert result.errors[0].code == "BULK_RNASEQ_REFERENCE_BINDING_CONFLICT"
    assert str(tmp_path) not in str(result.errors[0].to_dict())


def test_bulk_binding_rejects_asset_and_transcriptome_identity_drift(
    tmp_path: Path,
    bulk_rnaseq_qualifications,
):
    adapter = _adapter(tmp_path, bulk_rnaseq_qualifications)
    payload = _payload(tmp_path)
    Path(payload["reference"]["fasta"]).write_bytes(b"changed")
    assert adapter.verify_reference_profile_binding(payload).is_failure

    payload = _payload(tmp_path / "second")
    payload["transcriptome"]["fasta_sha256"] = "f" * 64
    result = adapter.verify_reference_profile_binding(payload)
    assert result.is_failure
    assert result.errors[0].code == "BULK_RNASEQ_REFERENCE_BINDING_INVALID"


def test_bulk_binding_rejects_unknown_fields_and_unconfigured_execution(
    tmp_path: Path,
    bulk_rnaseq_qualifications,
):
    payload = _payload(tmp_path)
    payload["unknown"] = True
    configured = _adapter(tmp_path, bulk_rnaseq_qualifications)
    assert configured.verify_reference_profile_binding(payload).is_failure

    from encode_pipeline.adapters.bulk_rnaseq import BulkRnaSeqWorkflowAdapter

    result = BulkRnaSeqWorkflowAdapter().verify_reference_profile_binding(
        _payload(tmp_path / "unconfigured")
    )
    assert result.is_failure
    assert result.errors[0].code == "BULK_RNASEQ_REFERENCE_BINDING_UNAVAILABLE"
