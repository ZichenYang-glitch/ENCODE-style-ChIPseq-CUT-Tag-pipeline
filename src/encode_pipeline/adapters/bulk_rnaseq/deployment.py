"""Fail-closed operator composition for the default bulk RNA-seq product."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
import json
import os
from pathlib import Path
from typing import Any

from encode_pipeline.adapters.bulk_rnaseq.adapter import (
    BulkRnaSeqResultsWorkflowAdapter,
    BulkRnaSeqWorkflowAdapter,
)
from encode_pipeline.adapters.bulk_rnaseq.execution import (
    BulkRnaSeqExecutionBinding,
    BulkRnaSeqTranscriptomeBinding,
)
from encode_pipeline.adapters.bulk_rnaseq.runtime_assets import RuntimeAssetBinding
from encode_pipeline.platform.adapters import WorkflowAvailability


RUNTIME_ROOT_ENV = "HELIXWEAVE_BULK_RNASEQ_RUNTIME_ROOT"
TRANSCRIPTOME_BINDING_MANIFEST_ENV = (
    "HELIXWEAVE_BULK_RNASEQ_TRANSCRIPTOME_BINDING_MANIFEST"
)
MANAGED_DOCKER_EXECUTABLE_ENV = "ENCODE_PIPELINE_MANAGED_DOCKER_EXECUTABLE"
MANAGED_DOCKER_SOCKET_ENV = "ENCODE_PIPELINE_MANAGED_DOCKER_SOCKET"
TRANSCRIPTOME_BINDING_SCHEMA_VERSION = "1.0.0"

_COORDINATE_NAMES = (
    RUNTIME_ROOT_ENV,
    TRANSCRIPTOME_BINDING_MANIFEST_ENV,
    MANAGED_DOCKER_EXECUTABLE_ENV,
    MANAGED_DOCKER_SOCKET_ENV,
)
_TRANSCRIPTOME_FIELDS = {
    "schema_version",
    "reference_id",
    "fasta_sha256",
    "gtf_sha256",
    "transcript_fasta",
    "transcript_fasta_sha256",
}
_MAX_MANIFEST_BYTES = 64 * 1024


@dataclass(frozen=True)
class BulkRnaSeqLocalExecutionConfiguration:
    """Private local-process coordinates from an admitted adapter binding."""

    executable: Path
    docker_executable: Path
    docker_socket: Path


def load_default_bulk_rnaseq_adapter(
    environ: Mapping[str, str] | None = None,
) -> BulkRnaSeqWorkflowAdapter:
    """Return authoring-only or fully admitted execution composition.

    Missing configuration is a normal authoring-only deployment. Partial,
    malformed, or unadmitted configuration remains visible as unavailable and
    never exposes the rejected operator coordinates.
    """

    source = os.environ if environ is None else environ
    configured = tuple(name for name in _COORDINATE_NAMES if source.get(name))
    if not configured:
        return BulkRnaSeqWorkflowAdapter()
    if len(configured) != len(_COORDINATE_NAMES):
        return _unavailable_adapter()

    try:
        runtime_root = _canonical_absolute_path(source[RUNTIME_ROOT_ENV])
        manifest_path = _canonical_absolute_path(
            source[TRANSCRIPTOME_BINDING_MANIFEST_ENV]
        )
        docker_executable = _canonical_absolute_path(
            source[MANAGED_DOCKER_EXECUTABLE_ENV]
        )
        docker_socket = _canonical_absolute_path(source[MANAGED_DOCKER_SOCKET_ENV])
        transcriptome = _load_transcriptome_binding(manifest_path)
        binding = BulkRnaSeqExecutionBinding(
            assets=RuntimeAssetBinding(
                root=runtime_root,
                docker_executable=docker_executable,
                docker_socket=docker_socket,
            ),
            transcriptome=transcriptome,
        )
        adapter = BulkRnaSeqResultsWorkflowAdapter(execution=binding)
    except Exception:
        return _unavailable_adapter()
    return adapter


def local_execution_configuration(
    adapter: object,
) -> BulkRnaSeqLocalExecutionConfiguration | None:
    """Return private runner coordinates only for a composed results adapter."""
    if not isinstance(adapter, BulkRnaSeqResultsWorkflowAdapter):
        return None
    binding = adapter.execution_binding
    if binding is None:
        return None
    return BulkRnaSeqLocalExecutionConfiguration(
        executable=binding.assets.network_isolation_executable,
        docker_executable=binding.assets.docker_executable,
        docker_socket=binding.assets.docker_socket,
    )


def disable_local_execution(adapter: object) -> None:
    """Permanently remove an uncomposable private execution binding."""
    if isinstance(adapter, BulkRnaSeqWorkflowAdapter):
        adapter.disable_execution()


def _unavailable_adapter() -> BulkRnaSeqWorkflowAdapter:
    return BulkRnaSeqWorkflowAdapter(
        configured_availability=WorkflowAvailability(
            execution="unavailable",
            reason_code="WORKFLOW_EXECUTION_UNAVAILABLE",
        )
    )


def _canonical_absolute_path(value: object) -> Path:
    if not isinstance(value, str) or not value:
        raise ValueError("deployment coordinate is invalid")
    path = Path(value)
    rendered = str(path)
    if (
        not path.is_absolute()
        or rendered != value
        or rendered != str(Path(rendered))
        or any(part == ".." for part in path.parts)
        or any(character in rendered for character in ("\x00", "\n", "\r"))
    ):
        raise ValueError("deployment coordinate is invalid")
    return path


def _load_transcriptome_binding(path: Path) -> BulkRnaSeqTranscriptomeBinding:
    if path.is_symlink() or not path.is_file():
        raise ValueError("transcriptome binding manifest is unavailable")
    size = path.stat().st_size
    if size <= 0 or size > _MAX_MANIFEST_BYTES:
        raise ValueError("transcriptome binding manifest size is invalid")
    document = json.loads(
        path.read_bytes(),
        object_pairs_hook=_unique_object,
    )
    if not isinstance(document, dict) or set(document) != _TRANSCRIPTOME_FIELDS:
        raise ValueError("transcriptome binding manifest fields are invalid")
    if document["schema_version"] != TRANSCRIPTOME_BINDING_SCHEMA_VERSION:
        raise ValueError("transcriptome binding manifest version is unsupported")
    transcript_fasta = _canonical_absolute_path(document["transcript_fasta"])
    return BulkRnaSeqTranscriptomeBinding(
        reference_id=document["reference_id"],
        fasta_sha256=document["fasta_sha256"],
        gtf_sha256=document["gtf_sha256"],
        transcript_fasta=transcript_fasta,
        transcript_fasta_sha256=document["transcript_fasta_sha256"],
    )


def _unique_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    value: dict[str, Any] = {}
    for key, item in pairs:
        if key in value:
            raise ValueError("duplicate JSON object key")
        value[key] = item
    return value
