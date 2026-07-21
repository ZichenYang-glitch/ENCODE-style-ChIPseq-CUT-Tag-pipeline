"""Immutable reference and prebuilt-index identities for bulk RNA-seq."""

from __future__ import annotations

from contextlib import suppress
from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import stat
from typing import Any, Mapping

from encode_pipeline.adapters.bulk_rnaseq.resource_closure import (
    DEFAULT_RESOURCE_CLOSURE_POLICY,
    ResourceClosurePolicy,
    VerifiedLocalFile,
    _DescriptorRoot,
    _PathObservation,
    _ResourceFailure,
    _absolute_lexical_path,
    _directory_flags,
    _directory_identity,
    _file_identity,
    _safe_leaf_name,
    _safe_relative_parts,
    safe_regular_file_bytes,
    safe_regular_file_identity,
)
from encode_pipeline.platform.results import Issue, Result


REFERENCE_CLOSURE_SCHEME = "sha256-framed-reference-v1"
REFERENCE_INDEX_MANIFEST = "helixweave-reference-index.json"
REFERENCE_INDEX_SCHEMA_VERSION = "1.1.0"
MAXIMUM_INDEX_MANIFEST_BYTES = 2 * 1024 * 1024
MAXIMUM_INDEX_FILES = 16_384
MAXIMUM_INDEX_TOTAL_BYTES = 512 * 1024**3

_SHA256 = re.compile(r"^[0-9a-f]{64}$")
_IMMUTABLE_IMAGE = re.compile(r"^[^@\s]+@sha256:[0-9a-f]{64}$")
_INDEX_KINDS = frozenset({"star", "salmon"})
_INDEX_PRODUCERS = {
    "star": {
        "process": "STAR_GENOMEGENERATE",
        "tool": "star",
        "tool_version": "2.7.11b",
        "build_contract": "nfcore-rnaseq-3.26.0-star-genomegenerate-default-v1",
    },
    "salmon": {
        "process": "SALMON_INDEX",
        "tool": "salmon",
        "tool_version": "1.10.3",
        "build_contract": "nfcore-rnaseq-3.26.0-salmon-index-default-v1",
    },
}


@dataclass(frozen=True)
class ReferenceIndexFile:
    """One verified regular file in a prebuilt reference index."""

    relative_path: str
    path: Path
    size_bytes: int
    sha256: str


@dataclass(frozen=True)
class ReferenceIndexClosure:
    """A prebuilt index manifest bound to the exact FASTA/GTF identities."""

    kind: str
    path: Path
    manifest_sha256: str
    identity_sha256: str
    producer_process: str
    producer_tool: str
    producer_tool_version: str
    producer_container_image: str
    build_parameters: tuple[tuple[str, object], ...]
    files: tuple[ReferenceIndexFile, ...]


@dataclass(frozen=True)
class ReferenceClosure:
    """Verified primary reference files plus optional prebuilt indexes."""

    reference_id: str
    fasta: VerifiedLocalFile
    gtf: VerifiedLocalFile
    star_index: ReferenceIndexClosure | None
    salmon_index: ReferenceIndexClosure | None
    identity_sha256: str


@dataclass
class _WalkState:
    policy: ResourceClosurePolicy
    files: list[str]
    directories: list[str]
    observations: list[_PathObservation]
    total_bytes: int = 0


class _DuplicateJsonKey(ValueError):
    pass


def verify_reference_closure(
    reference: Mapping[str, object],
    *,
    producer_images: Mapping[str, str] | None = None,
    index_build_parameters: Mapping[str, object] | None = None,
    transcript_fasta_sha256: str | None = None,
    policy: ResourceClosurePolicy = DEFAULT_RESOURCE_CLOSURE_POLICY,
) -> Result[ReferenceClosure]:
    """Verify explicit FASTA/GTF bytes and every supplied index manifest."""
    try:
        if not isinstance(reference, Mapping) or not isinstance(
            policy, ResourceClosurePolicy
        ):
            raise _ResourceFailure
        reference_id = reference.get("reference_id")
        annotation_style = reference.get("annotation_style", "ensembl")
        fasta_sha256 = reference.get("fasta_sha256")
        gtf_sha256 = reference.get("gtf_sha256")
        if (
            not isinstance(reference_id, str)
            or not reference_id
            or annotation_style not in {"ensembl", "gencode"}
            or not _is_sha256(fasta_sha256)
            or not _is_sha256(gtf_sha256)
        ):
            raise _ResourceFailure
        fasta_result = safe_regular_file_identity(
            reference.get("fasta"),
            expected_sha256=fasta_sha256,
            policy=policy,
        )
        gtf_result = safe_regular_file_identity(
            reference.get("gtf"),
            expected_sha256=gtf_sha256,
            policy=policy,
        )
        if fasta_result.is_failure or gtf_result.is_failure:
            raise _ResourceFailure

        indexes: dict[str, ReferenceIndexClosure | None] = {
            "star": None,
            "salmon": None,
        }
        for kind in sorted(indexes):
            value = reference.get(f"{kind}_index")
            if value is None:
                continue
            if not isinstance(value, Mapping):
                raise _ResourceFailure
            result = verify_reference_index(
                value.get("path"),
                kind=kind,
                expected_manifest_sha256=value.get("identity_sha256"),
                fasta_sha256=fasta_sha256,
                gtf_sha256=gtf_sha256,
                transcript_fasta_sha256=transcript_fasta_sha256,
                annotation_style=annotation_style,
                index_build_parameters=index_build_parameters,
                expected_container_image=_expected_container_image(
                    kind,
                    producer_images,
                ),
                policy=policy,
            )
            if result.is_failure:
                raise _ResourceFailure
            indexes[kind] = result.value

        identity = _reference_identity(
            reference_id=reference_id,
            fasta=fasta_result.value,
            gtf=gtf_result.value,
            indexes=indexes,
        )
        return Result.success(
            ReferenceClosure(
                reference_id=reference_id,
                fasta=fasta_result.value,
                gtf=gtf_result.value,
                star_index=indexes["star"],
                salmon_index=indexes["salmon"],
                identity_sha256=identity,
            )
        )
    except (_ResourceFailure, OSError, TypeError, ValueError, UnicodeError):
        return _failure("BULK_RNASEQ_REFERENCE_INVALID")


def verify_reference_index(
    path: str | Path,
    *,
    kind: str,
    expected_manifest_sha256: str,
    fasta_sha256: str,
    gtf_sha256: str,
    transcript_fasta_sha256: str | None = None,
    annotation_style: str,
    expected_container_image: str,
    index_build_parameters: Mapping[str, object] | None = None,
    policy: ResourceClosurePolicy = DEFAULT_RESOURCE_CLOSURE_POLICY,
) -> Result[ReferenceIndexClosure]:
    """Verify an exact index file set described by a versioned local manifest."""
    try:
        if (
            kind not in _INDEX_KINDS
            or not _is_sha256(expected_manifest_sha256)
            or not _is_sha256(fasta_sha256)
            or not _is_sha256(gtf_sha256)
            or (kind == "salmon" and not _is_sha256(transcript_fasta_sha256))
            or annotation_style not in {"ensembl", "gencode"}
            or not isinstance(expected_container_image, str)
            or _IMMUTABLE_IMAGE.fullmatch(expected_container_image) is None
            or not isinstance(policy, ResourceClosurePolicy)
        ):
            raise _ResourceFailure
        root = _absolute_path(path)
        manifest_path = root / REFERENCE_INDEX_MANIFEST
        manifest_result = safe_regular_file_bytes(
            manifest_path,
            maximum_bytes=MAXIMUM_INDEX_MANIFEST_BYTES,
            expected_sha256=expected_manifest_sha256,
            policy=policy,
        )
        if manifest_result.is_failure:
            raise _ResourceFailure
        manifest = _parse_manifest(manifest_result.value.content)
        expected_producer = {
            **_INDEX_PRODUCERS[kind],
            "container_image": expected_container_image,
        }
        if (
            manifest["index_kind"] != kind
            or manifest["reference"]["fasta_sha256"] != fasta_sha256
            or manifest["reference"]["gtf_sha256"] != gtf_sha256
            or (
                kind == "salmon"
                and manifest["reference"]["transcript_fasta_sha256"]
                != transcript_fasta_sha256
            )
            or manifest["producer"] != expected_producer
            or not _valid_build_parameters(
                kind,
                manifest["build_parameters"],
                annotation_style=annotation_style,
                index_build_parameters=index_build_parameters,
            )
            or len(manifest["files"]) > policy.maximum_index_entries
        ):
            raise _ResourceFailure

        declared = manifest["files"]
        declared_names = tuple(item["path"] for item in declared)
        with _DescriptorRoot.open(root) as reader:
            state = _WalkState(
                policy=policy,
                files=[],
                directories=[],
                observations=[],
            )
            _walk_tree(reader, reader.descriptor, (), state)
            reader.verify_stable(state.observations)
        observed_names = tuple(
            name for name in sorted(state.files) if name != REFERENCE_INDEX_MANIFEST
        )
        if observed_names != declared_names:
            raise _ResourceFailure
        expected_directories: set[str] = set()
        for name in declared_names:
            parent = PurePosixPath(name).parent
            while parent.as_posix() != ".":
                expected_directories.add(parent.as_posix())
                parent = parent.parent
        if set(state.directories) != expected_directories:
            raise _ResourceFailure

        verified_files: list[ReferenceIndexFile] = []
        total_bytes = 0
        for item in declared:
            result = safe_regular_file_identity(
                root.joinpath(*PurePosixPath(item["path"]).parts),
                expected_sha256=item["sha256"],
                maximum_bytes=policy.maximum_index_file_bytes,
                policy=policy,
            )
            if result.is_failure or result.value.size_bytes != item["size_bytes"]:
                raise _ResourceFailure
            total_bytes += result.value.size_bytes
            if total_bytes > MAXIMUM_INDEX_TOTAL_BYTES:
                raise _ResourceFailure
            verified_files.append(
                ReferenceIndexFile(
                    relative_path=item["path"],
                    path=result.value.path,
                    size_bytes=result.value.size_bytes,
                    sha256=result.value.sha256,
                )
            )

        identity = _index_identity(
            kind=kind,
            manifest_sha256=manifest_result.value.sha256,
            files=verified_files,
        )
        return Result.success(
            ReferenceIndexClosure(
                kind=kind,
                path=root,
                manifest_sha256=manifest_result.value.sha256,
                identity_sha256=identity,
                producer_process=manifest["producer"]["process"],
                producer_tool=manifest["producer"]["tool"],
                producer_tool_version=manifest["producer"]["tool_version"],
                producer_container_image=manifest["producer"]["container_image"],
                build_parameters=_freeze_build_parameters(manifest["build_parameters"]),
                files=tuple(verified_files),
            )
        )
    except (_ResourceFailure, OSError, TypeError, ValueError, UnicodeError):
        return _failure("BULK_RNASEQ_REFERENCE_INDEX_INVALID")


def _walk_tree(
    reader: _DescriptorRoot,
    descriptor: int,
    prefix: tuple[str, ...],
    state: _WalkState,
) -> None:
    try:
        names = sorted(os.listdir(descriptor))
    except OSError as error:
        raise _ResourceFailure from error
    for raw_name in names:
        name = _safe_leaf_name(raw_name)
        parts = (*prefix, name)
        relative = PurePosixPath(*parts).as_posix()
        if (
            len(state.files) + len(state.directories)
            >= min(MAXIMUM_INDEX_FILES, state.policy.maximum_index_entries)
            or len(parts) > state.policy.maximum_index_depth
            or len(relative.encode("utf-8")) > state.policy.maximum_relative_path_bytes
        ):
            raise _ResourceFailure
        listed = os.stat(name, dir_fd=descriptor, follow_symlinks=False)
        if stat.S_ISDIR(listed.st_mode):
            child = -1
            try:
                child = os.open(name, _directory_flags(), dir_fd=descriptor)
                before = os.fstat(child)
                if _directory_identity(listed) != _directory_identity(before):
                    raise _ResourceFailure
                state.directories.append(relative)
                _walk_tree(reader, child, parts, state)
                after = os.fstat(child)
                if _directory_identity(before) != _directory_identity(after):
                    raise _ResourceFailure
                observation = reader.reopen_identity(parts, is_directory=True)
                if observation.identities[-1] != _directory_identity(after):
                    raise _ResourceFailure
                state.observations.append(observation)
            finally:
                if child >= 0:
                    with suppress(OSError):
                        os.close(child)
            continue
        if not stat.S_ISREG(listed.st_mode) or listed.st_nlink != 1:
            raise _ResourceFailure
        observed = reader.read_regular_file(
            parts,
            maximum_bytes=(
                MAXIMUM_INDEX_MANIFEST_BYTES
                if parts == (REFERENCE_INDEX_MANIFEST,)
                else state.policy.maximum_index_file_bytes
            ),
            capture=False,
        )
        if observed.path_observation.identities[-1] != _file_identity(listed):
            raise _ResourceFailure
        state.total_bytes += observed.size_bytes
        if state.total_bytes > MAXIMUM_INDEX_TOTAL_BYTES:
            raise _ResourceFailure
        state.files.append(relative)
        state.observations.append(observed.path_observation)


def _parse_manifest(content: bytes) -> dict[str, Any]:
    try:
        value = json.loads(content.decode("utf-8"), object_pairs_hook=_unique_object)
    except (UnicodeDecodeError, json.JSONDecodeError, _DuplicateJsonKey) as error:
        raise _ResourceFailure from error
    if not isinstance(value, dict) or set(value) != {
        "schema_version",
        "index_kind",
        "producer",
        "build_parameters",
        "reference",
        "files",
    }:
        raise _ResourceFailure
    if value["schema_version"] != REFERENCE_INDEX_SCHEMA_VERSION:
        raise _ResourceFailure
    if value["index_kind"] not in _INDEX_KINDS:
        raise _ResourceFailure
    producer = value["producer"]
    if (
        not isinstance(producer, dict)
        or set(producer)
        != {
            "process",
            "tool",
            "tool_version",
            "build_contract",
            "container_image",
        }
        or not all(isinstance(item, str) and item for item in producer.values())
        or _IMMUTABLE_IMAGE.fullmatch(producer["container_image"]) is None
    ):
        raise _ResourceFailure
    build_parameters = value["build_parameters"]
    if not isinstance(build_parameters, dict):
        raise _ResourceFailure
    reference = value["reference"]
    expected_reference_keys = {"fasta_sha256", "gtf_sha256"}
    if value["index_kind"] == "salmon":
        expected_reference_keys.add("transcript_fasta_sha256")
    if not isinstance(reference, dict) or set(reference) != expected_reference_keys:
        raise _ResourceFailure
    if not all(_is_sha256(item) for item in reference.values()):
        raise _ResourceFailure
    files = value["files"]
    if not isinstance(files, list) or not files or len(files) > MAXIMUM_INDEX_FILES:
        raise _ResourceFailure
    normalized: list[dict[str, object]] = []
    seen: set[str] = set()
    for item in files:
        if not isinstance(item, dict) or set(item) != {
            "path",
            "size_bytes",
            "sha256",
        }:
            raise _ResourceFailure
        path = item["path"]
        size = item["size_bytes"]
        digest = item["sha256"]
        if not isinstance(path, str):
            raise _ResourceFailure
        parts = _safe_relative_parts(path, maximum_bytes=4096)
        normalized_path = PurePosixPath(*parts).as_posix()
        if (
            normalized_path == REFERENCE_INDEX_MANIFEST
            or normalized_path in seen
            or isinstance(size, bool)
            or not isinstance(size, int)
            or size < 0
            or not _is_sha256(digest)
        ):
            raise _ResourceFailure
        seen.add(normalized_path)
        normalized.append(
            {"path": normalized_path, "size_bytes": size, "sha256": digest}
        )
    normalized.sort(key=lambda item: item["path"])
    if [item["path"] for item in files] != [item["path"] for item in normalized]:
        raise _ResourceFailure
    return {
        "schema_version": value["schema_version"],
        "index_kind": value["index_kind"],
        "producer": producer,
        "build_parameters": build_parameters,
        "reference": reference,
        "files": normalized,
    }


def _expected_container_image(
    kind: str,
    producer_images: Mapping[str, str] | None,
) -> str:
    if not isinstance(producer_images, Mapping):
        raise _ResourceFailure
    process = _INDEX_PRODUCERS[kind]["process"]
    image = producer_images.get(process)
    if not isinstance(image, str) or _IMMUTABLE_IMAGE.fullmatch(image) is None:
        raise _ResourceFailure
    return image


def _valid_build_parameters(
    kind: str,
    value: object,
    *,
    annotation_style: str,
    index_build_parameters: Mapping[str, object] | None,
) -> bool:
    if not isinstance(value, dict) or (
        index_build_parameters is not None
        and not isinstance(index_build_parameters, Mapping)
    ):
        return False
    expected = index_build_parameters or {}
    skip_gtf_filter = expected.get("skip_gtf_filter", False)
    skip_gtf_transcript_filter = expected.get("skip_gtf_transcript_filter", False)
    gffread_transcript_fasta = expected.get("gffread_transcript_fasta", False)
    if any(
        not isinstance(item, bool)
        for item in (
            skip_gtf_filter,
            skip_gtf_transcript_filter,
            gffread_transcript_fasta,
        )
    ):
        return False
    if kind == "star":
        genome_sa_index_nbases = value.get("genome_sa_index_nbases")
        return (
            set(value)
            == {
                "extra_args",
                "genome_sa_index_nbases",
                "genome_sa_index_nbases_strategy",
                "sjdb_gtf_feature",
                "sjdb_overhang",
                "skip_gtf_filter",
                "skip_gtf_transcript_filter",
            }
            and value["extra_args"] == []
            and isinstance(genome_sa_index_nbases, int)
            and not isinstance(genome_sa_index_nbases, bool)
            and 1 <= genome_sa_index_nbases <= 14
            and value["genome_sa_index_nbases_strategy"] == "nfcore-rnaseq-3.26.0-auto"
            and value["sjdb_gtf_feature"] == "exon"
            and value["sjdb_overhang"] == 100
            and value["skip_gtf_filter"] is skip_gtf_filter
            and value["skip_gtf_transcript_filter"] is skip_gtf_transcript_filter
        )
    return (
        set(value)
        == {
            "decoy_mode",
            "extra_args",
            "gencode",
            "gffread_transcript_fasta",
            "kmer_size",
            "skip_gtf_filter",
            "skip_gtf_transcript_filter",
        }
        and value["decoy_mode"] == "gentrome"
        and value["extra_args"] == []
        and value["gencode"] is (annotation_style == "gencode")
        and value["gffread_transcript_fasta"] is gffread_transcript_fasta
        and value["kmer_size"] == 31
        and value["skip_gtf_filter"] is skip_gtf_filter
        and value["skip_gtf_transcript_filter"] is skip_gtf_transcript_filter
    )


def _freeze_build_parameters(
    value: Mapping[str, object],
) -> tuple[tuple[str, object], ...]:
    return tuple(
        (name, tuple(item) if isinstance(item, list) else item)
        for name, item in sorted(value.items())
    )


def _reference_identity(
    *,
    reference_id: str,
    fasta: VerifiedLocalFile,
    gtf: VerifiedLocalFile,
    indexes: Mapping[str, ReferenceIndexClosure | None],
) -> str:
    digest = hashlib.sha256()
    _frame(digest, REFERENCE_CLOSURE_SCHEME.encode())
    _frame(digest, reference_id.encode())
    for name, file in (("fasta", fasta), ("gtf", gtf)):
        _frame(digest, name.encode())
        _frame(digest, str(file.size_bytes).encode())
        _frame(digest, bytes.fromhex(file.sha256))
    for kind in sorted(indexes):
        value = indexes[kind]
        if value is not None:
            _frame(digest, kind.encode())
            _frame(digest, bytes.fromhex(value.identity_sha256))
    return digest.hexdigest()


def _index_identity(
    *,
    kind: str,
    manifest_sha256: str,
    files: list[ReferenceIndexFile],
) -> str:
    digest = hashlib.sha256()
    _frame(digest, b"sha256-framed-reference-index-v1")
    _frame(digest, kind.encode())
    _frame(digest, bytes.fromhex(manifest_sha256))
    for file in files:
        _frame(digest, file.relative_path.encode())
        _frame(digest, str(file.size_bytes).encode())
        _frame(digest, bytes.fromhex(file.sha256))
    return digest.hexdigest()


def _unique_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    value: dict[str, Any] = {}
    for key, item in pairs:
        if key in value:
            raise _DuplicateJsonKey
        value[key] = item
    return value


def _is_sha256(value: object) -> bool:
    return isinstance(value, str) and _SHA256.fullmatch(value) is not None


def _absolute_path(value: object) -> Path:
    if not isinstance(value, (str, Path)) or isinstance(value, bool):
        raise _ResourceFailure
    return _absolute_lexical_path(Path(os.fspath(value)))


def _frame(digest: Any, value: bytes) -> None:
    digest.update(len(value).to_bytes(8, "big"))
    digest.update(value)


def _failure(code: str) -> Result[Any]:
    return Result.failure(
        [
            Issue(
                code=code,
                message="Bulk RNA-seq reference verification failed.",
                severity="error",
                path="config.standard.reference",
                source="adapter",
            )
        ]
    )
