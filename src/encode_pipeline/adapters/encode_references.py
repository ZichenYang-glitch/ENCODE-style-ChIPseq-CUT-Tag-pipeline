"""Private Reference Profile binding for the ENCODE-style adapter."""

from __future__ import annotations

from collections.abc import Mapping
import csv
import hashlib
import io
import os
from pathlib import Path
import re
from typing import Any

from encode_pipeline.adapters.bulk_rnaseq.resource_closure import (
    VerifiedLocalFile,
    safe_regular_file_bytes,
    safe_regular_file_identity,
)
from encode_pipeline.config.defaults import BT2_LARGE, BT2_STANDARD
from encode_pipeline.platform.adapters import (
    MAX_AUTHORING_REQUEST_BYTES,
    MAX_SAMPLE_CELL_LENGTH,
    MAX_SAMPLE_COLUMNS,
    MAX_SAMPLE_ROWS,
    WorkflowInputs,
)
from encode_pipeline.platform.reference_profiles import (
    ADAPTER_REFERENCE_BINDING_IDENTITY_SCHEME,
    AdapterReferenceBindingIdentity,
)
from encode_pipeline.platform.results import Issue, Result


ENCODE_REFERENCE_BINDING_CONTRACT = "encode-reference-binding-v1"

_WORKFLOW_ID = "encode-style-chipseq-cuttag-atac-mnase"
_ASSEMBLY = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,127}$")
_SHA256 = re.compile(r"^[0-9a-f]{64}$")
_RESOURCE_NAMES = ("blacklist", "chrom_sizes", "gtf", "reference_fasta")
_TOP_LEVEL_FIELDS = {
    "assembly",
    "bowtie2_index",
    "effective_genome_size",
    "genome_resources",
    "schema_version",
}
_MAX_BINDING_STRING_BYTES = 16 * 1024


class _BindingInvalid(Exception):
    pass


def verify_encode_reference_profile_binding(
    payload: Mapping[str, object],
) -> Result[AdapterReferenceBindingIdentity]:
    """Verify every private reference asset and return a path-free identity."""
    try:
        _resolved, identity = _verify_binding(payload)
    except (OSError, TypeError, ValueError, UnicodeError, _BindingInvalid):
        return _failure("ENCODE_REFERENCE_BINDING_INVALID")
    return Result.success(identity)


def resolve_encode_reference_profile(
    inputs: WorkflowInputs,
    payload: Mapping[str, object],
) -> Result[tuple[WorkflowInputs, AdapterReferenceBindingIdentity]]:
    """Inject a verified reference into fresh workflow inputs, never caller paths."""
    if not isinstance(inputs, WorkflowInputs):
        return _failure("ENCODE_REFERENCE_BINDING_INVALID")
    try:
        resolved, identity = _verify_binding(payload)
    except (OSError, TypeError, ValueError, UnicodeError, _BindingInvalid):
        return _failure("ENCODE_REFERENCE_BINDING_INVALID")
    try:
        config = dict(inputs.config)
        expected_resources = {
            resolved["assembly"]: {
                "effective_genome_size": resolved["effective_genome_size"],
                **{
                    name: str(resolved["resources"][name].path)
                    for name in _RESOURCE_NAMES
                },
            }
        }
        if "genome_resources" in config:
            raise _BindingInvalid
        config["genome_resources"] = expected_resources

        sample_rows = _sample_rows(inputs.samples)
        assembly = resolved["assembly"]
        index_prefix = str(resolved["index_prefix"])
        for row in sample_rows:
            if "genome" in row or "bowtie2_index" in row:
                raise _BindingInvalid
            row["genome"] = assembly
            row["bowtie2_index"] = index_prefix
        return Result.success(
            (
                WorkflowInputs(
                    config=config,
                    samples=sample_rows,
                    options=inputs.options,
                ),
                identity,
            )
        )
    except (OSError, TypeError, ValueError, UnicodeError, _BindingInvalid):
        return _failure("ENCODE_REFERENCE_BINDING_CONFLICT")


def _verify_binding(
    payload: Mapping[str, object],
) -> tuple[dict[str, Any], AdapterReferenceBindingIdentity]:
    if not isinstance(payload, Mapping) or set(payload) != _TOP_LEVEL_FIELDS:
        raise _BindingInvalid
    if payload.get("schema_version") != ENCODE_REFERENCE_BINDING_CONTRACT:
        raise _BindingInvalid
    assembly = payload.get("assembly")
    if not isinstance(assembly, str) or _ASSEMBLY.fullmatch(assembly) is None:
        raise _BindingInvalid
    effective_genome_size = payload.get("effective_genome_size")
    if isinstance(effective_genome_size, bool) or not (
        (isinstance(effective_genome_size, int) and effective_genome_size > 0)
        or effective_genome_size in {"hs", "mm"}
    ):
        raise _BindingInvalid

    resource_payload = payload.get("genome_resources")
    if not isinstance(resource_payload, Mapping) or set(resource_payload) != set(
        _RESOURCE_NAMES
    ):
        raise _BindingInvalid
    resources: dict[str, VerifiedLocalFile] = {}
    for name in _RESOURCE_NAMES:
        descriptor = resource_payload.get(name)
        if not isinstance(descriptor, Mapping) or set(descriptor) != {"path", "sha256"}:
            raise _BindingInvalid
        path_value = descriptor.get("path")
        expected = descriptor.get("sha256")
        if (
            not isinstance(path_value, str)
            or len(path_value.encode()) > _MAX_BINDING_STRING_BYTES
            or not isinstance(expected, str)
            or _SHA256.fullmatch(expected) is None
        ):
            raise _BindingInvalid
        result = safe_regular_file_identity(path_value, expected_sha256=expected)
        if result.is_failure:
            raise _BindingInvalid
        assert result.value is not None
        resources[name] = result.value

    index = payload.get("bowtie2_index")
    if not isinstance(index, Mapping) or set(index) != {"prefix", "files"}:
        raise _BindingInvalid
    prefix_value = index.get("prefix")
    if (
        not isinstance(prefix_value, str)
        or len(prefix_value.encode()) > _MAX_BINDING_STRING_BYTES
    ):
        raise _BindingInvalid
    prefix = Path(prefix_value)
    if not prefix.is_absolute() or str(prefix) != prefix_value or ".." in prefix.parts:
        raise _BindingInvalid
    file_digests = index.get("files")
    if not isinstance(file_digests, Mapping):
        raise _BindingInvalid
    standard = tuple(template.format(prefix="")[0:] for template in BT2_STANDARD)
    large = tuple(template.format(prefix="")[0:] for template in BT2_LARGE)
    supplied_suffixes = set(file_digests)
    if supplied_suffixes == set(standard):
        suffixes = standard
    elif supplied_suffixes == set(large):
        suffixes = large
    else:
        raise _BindingInvalid
    index_files: dict[str, VerifiedLocalFile] = {}
    for suffix in suffixes:
        expected = file_digests.get(suffix)
        if not isinstance(expected, str) or _SHA256.fullmatch(expected) is None:
            raise _BindingInvalid
        result = safe_regular_file_identity(
            Path(f"{prefix}{suffix}"), expected_sha256=expected
        )
        if result.is_failure:
            raise _BindingInvalid
        assert result.value is not None
        index_files[suffix] = result.value
    _verify_exact_index_set(prefix, suffixes)

    digest = hashlib.sha256()
    _frame(digest, _WORKFLOW_ID.encode())
    _frame(digest, ENCODE_REFERENCE_BINDING_CONTRACT.encode())
    _frame(digest, assembly.encode())
    _frame(digest, str(effective_genome_size).encode())
    for name in _RESOURCE_NAMES:
        resource = resources[name]
        _frame(digest, name.encode())
        _frame(digest, str(resource.size_bytes).encode())
        _frame(digest, bytes.fromhex(resource.sha256))
    for suffix in suffixes:
        file = index_files[suffix]
        _frame(digest, suffix.encode())
        _frame(digest, str(file.size_bytes).encode())
        _frame(digest, bytes.fromhex(file.sha256))
    identity = AdapterReferenceBindingIdentity(
        workflow_id=_WORKFLOW_ID,
        contract_version=ENCODE_REFERENCE_BINDING_CONTRACT,
        identity_scheme=ADAPTER_REFERENCE_BINDING_IDENTITY_SCHEME,
        identity_sha256=digest.hexdigest(),
    )
    return (
        {
            "assembly": assembly,
            "effective_genome_size": effective_genome_size,
            "resources": resources,
            "index_prefix": prefix,
        },
        identity,
    )


def _verify_exact_index_set(prefix: Path, suffixes: tuple[str, ...]) -> None:
    flags = os.O_RDONLY | getattr(os, "O_DIRECTORY", 0) | getattr(os, "O_CLOEXEC", 0)
    flags |= getattr(os, "O_NOFOLLOW", 0)
    descriptor = os.open(prefix.parent, flags)
    try:
        before = os.fstat(descriptor)
        names = tuple(os.listdir(descriptor))
        after = os.fstat(descriptor)
    finally:
        os.close(descriptor)
    if (before.st_dev, before.st_ino, before.st_ctime_ns) != (
        after.st_dev,
        after.st_ino,
        after.st_ctime_ns,
    ):
        raise _BindingInvalid
    expected = {f"{prefix.name}{suffix}" for suffix in suffixes}
    observed = {
        name
        for name in names
        if name.startswith(f"{prefix.name}.")
        and (name.endswith(".bt2") or name.endswith(".bt2l"))
    }
    if observed != expected:
        raise _BindingInvalid


def _sample_rows(samples: object) -> list[dict[str, str]]:
    if isinstance(samples, list):
        return [dict(row) for row in samples]
    if not isinstance(samples, (str, Path)):
        raise _BindingInvalid
    read = safe_regular_file_bytes(
        samples,
        maximum_bytes=MAX_AUTHORING_REQUEST_BYTES,
    )
    if read.is_failure:
        raise _BindingInvalid
    assert read.value is not None
    try:
        reader = csv.DictReader(
            io.StringIO(read.value.content.decode("utf-8")), delimiter="\t"
        )
        if reader.fieldnames is None or len(reader.fieldnames) > MAX_SAMPLE_COLUMNS:
            raise _BindingInvalid
        rows = [dict(row) for row in reader]
    except (csv.Error, UnicodeError):
        raise _BindingInvalid from None
    if (
        not rows
        or len(rows) > MAX_SAMPLE_ROWS
        or any(
            len(str(value)) > MAX_SAMPLE_CELL_LENGTH
            for row in rows
            for value in row.values()
        )
    ):
        raise _BindingInvalid
    return rows


def _frame(digest: Any, value: bytes) -> None:
    digest.update(len(value).to_bytes(8, "big"))
    digest.update(value)


def _failure(code: str) -> Result[Any]:
    return Result.failure(
        [
            Issue(
                code=code,
                message="ENCODE reference binding is unavailable.",
                severity="error",
                path="reference_revision_id",
                source="adapter",
            )
        ]
    )
