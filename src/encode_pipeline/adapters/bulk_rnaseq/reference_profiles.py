"""Private Reference Profile binding for pinned bulk RNA-seq execution."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
import hashlib
from pathlib import Path
import re
from typing import TYPE_CHECKING, Any

from encode_pipeline.adapters.bulk_rnaseq.execution import (
    BulkRnaSeqExecutionBinding,
    BulkRnaSeqTranscriptomeBinding,
)
from encode_pipeline.adapters.bulk_rnaseq.reference_closure import (
    ReferenceClosure,
    verify_reference_closure,
)
from encode_pipeline.adapters.bulk_rnaseq.resource_closure import (
    VerifiedLocalFile,
    safe_regular_file_identity,
)
from encode_pipeline.adapters.bulk_rnaseq.validation import (
    validate_bulk_rnaseq_inputs,
)
from encode_pipeline.platform.adapters import WorkflowInputs
from encode_pipeline.platform.reference_profiles import (
    ADAPTER_REFERENCE_BINDING_IDENTITY_SCHEME,
    AdapterReferenceBindingIdentity,
)
from encode_pipeline.platform.results import Issue, Result

if TYPE_CHECKING:
    from encode_pipeline.adapters.bulk_rnaseq.adapter import BulkRnaSeqWorkflowAdapter


BULK_RNASEQ_REFERENCE_BINDING_CONTRACT = "bulk-rnaseq-reference-binding-v1"

_WORKFLOW_ID = "bulk-rnaseq"
_SHA256 = re.compile(r"^[0-9a-f]{64}$")
_REFERENCE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,127}$")
_TOP_LEVEL_FIELDS = {"reference", "schema_version", "transcriptome"}
_REFERENCE_REQUIRED = {"fasta", "fasta_sha256", "gtf", "gtf_sha256", "reference_id"}
_REFERENCE_OPTIONAL = {"annotation_style", "salmon_index", "star_index"}
_TRANSCRIPTOME_FIELDS = {
    "fasta_sha256",
    "gtf_sha256",
    "reference_id",
    "transcript_fasta",
    "transcript_fasta_sha256",
}


class _BindingInvalid(Exception):
    pass


@dataclass(frozen=True)
class _VerifiedBinding:
    reference_payload: dict[str, object]
    reference: ReferenceClosure
    transcriptome: BulkRnaSeqTranscriptomeBinding
    transcript_fasta: VerifiedLocalFile
    identity: AdapterReferenceBindingIdentity


def verify_bulk_rnaseq_reference_profile_binding(
    adapter: BulkRnaSeqWorkflowAdapter,
    payload: Mapping[str, object],
) -> Result[AdapterReferenceBindingIdentity]:
    """Verify an exact reference/transcriptome pair without returning paths."""
    if adapter.execution_binding is None:
        return _failure("BULK_RNASEQ_REFERENCE_BINDING_UNAVAILABLE")
    try:
        verified = _verify_binding(adapter.execution_binding, payload)
    except (OSError, TypeError, ValueError, UnicodeError, _BindingInvalid):
        return _failure("BULK_RNASEQ_REFERENCE_BINDING_INVALID")
    return Result.success(verified.identity)


def resolve_bulk_rnaseq_reference_profile(
    adapter: BulkRnaSeqWorkflowAdapter,
    inputs: WorkflowInputs,
    payload: Mapping[str, object],
) -> Result[
    tuple[WorkflowInputs, BulkRnaSeqExecutionBinding, AdapterReferenceBindingIdentity]
]:
    """Bind private reference assets to fresh inputs and a run-scoped runtime."""
    execution = adapter.execution_binding
    if execution is None:
        return _failure("BULK_RNASEQ_REFERENCE_BINDING_UNAVAILABLE")
    if not isinstance(inputs, WorkflowInputs):
        return _failure("BULK_RNASEQ_REFERENCE_BINDING_INVALID")
    try:
        config = dict(inputs.config)
        standard_value = config.get("standard")
        if not isinstance(standard_value, Mapping):
            raise _BindingInvalid
        standard = dict(standard_value)
        raw_reference = (
            payload.get("reference") if isinstance(payload, Mapping) else None
        )
        if "reference" in standard:
            return _failure("BULK_RNASEQ_REFERENCE_BINDING_CONFLICT")
        standard["reference"] = (
            dict(raw_reference) if isinstance(raw_reference, Mapping) else raw_reference
        )
        config["standard"] = standard
        resolved_inputs = WorkflowInputs(
            config=config,
            samples=inputs.samples,
            options=inputs.options,
        )
        validation = validate_bulk_rnaseq_inputs(resolved_inputs)
        if validation.is_failure:
            return Result.failure(validation.issues)
        if not isinstance(validation.value, Mapping) or not isinstance(
            validation.value.get("nfcore_params"), Mapping
        ):
            raise _BindingInvalid
        verified = _verify_binding(
            execution,
            payload,
            index_build_parameters=validation.value["nfcore_params"],
        )
        run_execution = BulkRnaSeqExecutionBinding(
            assets=execution.assets,
            transcriptome=verified.transcriptome,
            implementation_qualification=execution.implementation_qualification,
            container_uid=execution.container_uid,
            container_gid=execution.container_gid,
            resume_enabled=execution.resume_enabled,
            resource_policy=execution.resource_policy,
        )
        return Result.success((resolved_inputs, run_execution, verified.identity))
    except (OSError, TypeError, ValueError, UnicodeError, _BindingInvalid):
        return _failure("BULK_RNASEQ_REFERENCE_BINDING_INVALID")


def _verify_binding(
    execution: BulkRnaSeqExecutionBinding,
    payload: Mapping[str, object],
    *,
    index_build_parameters: Mapping[str, object] | None = None,
) -> _VerifiedBinding:
    if not isinstance(payload, Mapping) or set(payload) != _TOP_LEVEL_FIELDS:
        raise _BindingInvalid
    if payload.get("schema_version") != BULK_RNASEQ_REFERENCE_BINDING_CONTRACT:
        raise _BindingInvalid
    raw_reference = payload.get("reference")
    if not isinstance(raw_reference, Mapping):
        raise _BindingInvalid
    reference_fields = set(raw_reference)
    if not _REFERENCE_REQUIRED.issubset(
        reference_fields
    ) or reference_fields.difference(_REFERENCE_REQUIRED | _REFERENCE_OPTIONAL):
        raise _BindingInvalid
    reference = dict(raw_reference)
    reference_id = reference.get("reference_id")
    if (
        not isinstance(reference_id, str)
        or _REFERENCE_ID.fullmatch(reference_id) is None
    ):
        raise _BindingInvalid
    if reference.get("annotation_style", "ensembl") not in {"ensembl", "gencode"}:
        raise _BindingInvalid
    for path_name, digest_name in (
        ("fasta", "fasta_sha256"),
        ("gtf", "gtf_sha256"),
    ):
        path_value = reference.get(path_name)
        digest_value = reference.get(digest_name)
        if (
            not isinstance(path_value, str)
            or not isinstance(digest_value, str)
            or _SHA256.fullmatch(digest_value) is None
        ):
            raise _BindingInvalid
    for index_name in ("star_index", "salmon_index"):
        index = reference.get(index_name)
        if index is None:
            continue
        if not isinstance(index, Mapping) or set(index) != {
            "path",
            "identity_sha256",
        }:
            raise _BindingInvalid
        if (
            not isinstance(index.get("path"), str)
            or not isinstance(index.get("identity_sha256"), str)
            or _SHA256.fullmatch(index["identity_sha256"]) is None
        ):
            raise _BindingInvalid

    raw_transcriptome = payload.get("transcriptome")
    if (
        not isinstance(raw_transcriptome, Mapping)
        or set(raw_transcriptome) != _TRANSCRIPTOME_FIELDS
    ):
        raise _BindingInvalid
    transcript_path = raw_transcriptome.get("transcript_fasta")
    if not isinstance(transcript_path, str):
        raise _BindingInvalid
    transcript_reference_id = raw_transcriptome.get("reference_id")
    transcript_fasta_digest = raw_transcriptome.get("fasta_sha256")
    transcript_gtf_digest = raw_transcriptome.get("gtf_sha256")
    transcript_sequence_digest = raw_transcriptome.get("transcript_fasta_sha256")
    if (
        not isinstance(transcript_reference_id, str)
        or not isinstance(transcript_fasta_digest, str)
        or _SHA256.fullmatch(transcript_fasta_digest) is None
        or not isinstance(transcript_gtf_digest, str)
        or _SHA256.fullmatch(transcript_gtf_digest) is None
        or not isinstance(transcript_sequence_digest, str)
        or _SHA256.fullmatch(transcript_sequence_digest) is None
    ):
        raise _BindingInvalid
    transcriptome = BulkRnaSeqTranscriptomeBinding(
        reference_id=transcript_reference_id,
        fasta_sha256=transcript_fasta_digest,
        gtf_sha256=transcript_gtf_digest,
        transcript_fasta=Path(transcript_path),
        transcript_fasta_sha256=transcript_sequence_digest,
    )
    if (
        transcriptome.reference_id != reference_id
        or transcriptome.fasta_sha256 != reference.get("fasta_sha256")
        or transcriptome.gtf_sha256 != reference.get("gtf_sha256")
    ):
        raise _BindingInvalid
    transcript_result = safe_regular_file_identity(
        transcriptome.transcript_fasta,
        expected_sha256=transcriptome.transcript_fasta_sha256,
        policy=execution.resource_policy,
    )
    if transcript_result.is_failure:
        raise _BindingInvalid
    assert transcript_result.value is not None

    producer_images: Mapping[str, str] | None = None
    if "star_index" in reference or "salmon_index" in reference:
        assets = execution.runtime_admission.acquire()
        if assets.is_failure:
            raise _BindingInvalid
        assert assets.value is not None
        producer_images = {item.process: item.image for item in assets.value.containers}
    reference_result = verify_reference_closure(
        reference,
        producer_images=producer_images,
        index_build_parameters=index_build_parameters,
        transcript_fasta_sha256=transcriptome.transcript_fasta_sha256,
        policy=execution.resource_policy,
    )
    if reference_result.is_failure:
        raise _BindingInvalid
    assert reference_result.value is not None

    digest = hashlib.sha256()
    _frame(digest, _WORKFLOW_ID.encode())
    _frame(digest, BULK_RNASEQ_REFERENCE_BINDING_CONTRACT.encode())
    _frame(digest, bytes.fromhex(reference_result.value.identity_sha256))
    _frame(digest, str(reference.get("annotation_style", "ensembl")).encode())
    _frame(digest, str(transcript_result.value.size_bytes).encode())
    _frame(digest, bytes.fromhex(transcript_result.value.sha256))
    identity = AdapterReferenceBindingIdentity(
        workflow_id=_WORKFLOW_ID,
        contract_version=BULK_RNASEQ_REFERENCE_BINDING_CONTRACT,
        identity_scheme=ADAPTER_REFERENCE_BINDING_IDENTITY_SCHEME,
        identity_sha256=digest.hexdigest(),
    )
    return _VerifiedBinding(
        reference_payload=reference,
        reference=reference_result.value,
        transcriptome=transcriptome,
        transcript_fasta=transcript_result.value,
        identity=identity,
    )


def _frame(digest: Any, value: bytes) -> None:
    digest.update(len(value).to_bytes(8, "big"))
    digest.update(value)


def _failure(code: str) -> Result[Any]:
    return Result.failure(
        [
            Issue(
                code=code,
                message="Bulk RNA-seq reference binding is unavailable.",
                severity="error",
                path="reference_revision_id",
                source="adapter",
            )
        ]
    )
