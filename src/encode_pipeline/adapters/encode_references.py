"""Private Reference Profile binding for the ENCODE-style adapter."""

from __future__ import annotations

from collections.abc import Iterable, Mapping
import csv
import hashlib
import io
import os
from pathlib import Path
import re
import stat
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
ENCODE_REFERENCE_BINDING_CONTRACT_V2 = "encode-reference-binding-v2"

_WORKFLOW_ID = "encode-style-chipseq-cuttag-atac-mnase"
_ASSEMBLY = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,127}$")
_SHA256 = re.compile(r"^[0-9a-f]{64}$")
_RESOURCE_NAMES = ("blacklist", "chrom_sizes", "gtf", "reference_fasta")
_V2_RESOURCE_NAMES = _RESOURCE_NAMES + ("reference_fai", "reference_dict")
_V2_REQUIRED_RESOURCES = _RESOURCE_NAMES + ("reference_fai",)
_TOP_LEVEL_FIELDS = {
    "assembly",
    "bowtie2_index",
    "effective_genome_size",
    "genome_resources",
    "schema_version",
}
_V2_OPTIONAL_FIELDS = {"organelle_chroms", "canonical_chroms"}
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
                **(
                    {"organelle_chroms": list(resolved["organelle_chroms"])}
                    if "organelle_chroms" in resolved
                    else {}
                ),
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
    if not isinstance(payload, Mapping):
        raise _BindingInvalid
    schema_version = payload.get("schema_version")
    if schema_version == ENCODE_REFERENCE_BINDING_CONTRACT:
        return _verify_binding_v1(payload)
    if schema_version == ENCODE_REFERENCE_BINDING_CONTRACT_V2:
        return _verify_binding_v2(payload)
    raise _BindingInvalid


def _verify_binding_v1(
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


_READ_CHUNK_BYTES = 1024 * 1024


def _iter_verified_raw_lines(
    path_value: str | Path,
    expected_sha256: str,
    *,
    maximum_bytes: int | None = None,
):
    """Yield (line, absolute_offset, newline_width) from one verified read."""
    path = Path(path_value)
    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0)
    descriptor = os.open(path, flags)
    try:
        before = os.fstat(descriptor)
        if not stat.S_ISREG(before.st_mode) or before.st_nlink != 1:
            raise _BindingInvalid
        if maximum_bytes is not None and before.st_size > maximum_bytes:
            raise _BindingInvalid
        digest = hashlib.sha256()
        buffer = bytearray()
        base_offset = 0
        remaining = before.st_size
        while remaining:
            chunk = os.read(descriptor, min(_READ_CHUNK_BYTES, remaining))
            if not chunk:
                raise _BindingInvalid
            digest.update(chunk)
            remaining -= len(chunk)
            buffer.extend(chunk)
            start = 0
            while True:
                newline = buffer.find(b"\n", start)
                if newline < 0:
                    break
                content = bytes(buffer[start:newline])
                newline_width = 1
                if content.endswith(b"\r"):
                    content = content[:-1]
                    newline_width = 2
                yield content, base_offset + start, newline_width
                start = newline + 1
            if start:
                del buffer[:start]
                base_offset += start
        if buffer:
            content = bytes(buffer)
            if content.endswith(b"\r"):
                content = content[:-1]
            yield content, base_offset, 0
        if os.read(descriptor, 1):
            raise _BindingInvalid
        after = os.fstat(descriptor)
        if (
            before.st_dev,
            before.st_ino,
            before.st_size,
            before.st_mtime_ns,
            before.st_ctime_ns,
        ) != (
            after.st_dev,
            after.st_ino,
            after.st_size,
            after.st_mtime_ns,
            after.st_ctime_ns,
        ):
            raise _BindingInvalid
        if digest.hexdigest() != expected_sha256:
            raise _BindingInvalid
    finally:
        os.close(descriptor)


def _iter_verified_lines(
    path_value: str | Path,
    expected_sha256: str,
    *,
    maximum_bytes: int | None = None,
):
    """Yield lines from one stable no-follow read and verify its digest."""
    for content, _offset, _newline_width in _iter_verified_raw_lines(
        path_value,
        expected_sha256,
        maximum_bytes=maximum_bytes,
    ):
        yield content


def _parse_fai_lines(
    lines: Iterable[bytes],
) -> list[tuple[str, int, int, int, int]]:
    rows = []
    for raw in lines:
        line = raw.decode("utf-8")
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) < 5:
            raise _BindingInvalid
        try:
            length = int(fields[1])
            offset = int(fields[2])
            line_bases = int(fields[3])
            line_width = int(fields[4])
        except ValueError as exc:
            raise _BindingInvalid from exc
        if (
            not fields[0]
            or length <= 0
            or offset < 0
            or line_bases <= 0
            or line_width < line_bases
        ):
            raise _BindingInvalid
        rows.append((fields[0], length, offset, line_bases, line_width))
    if not rows:
        raise _BindingInvalid
    return rows


def _parse_fasta_layout(
    lines: Iterable[tuple[bytes, int, int]],
) -> list[tuple[str, int, int, int, int]]:
    rows: list[tuple[str, int, int, int, int]] = []
    name = None
    length = 0
    first_offset = 0
    line_bases = 0
    line_width = 0
    pending: tuple[bytes, int] | None = None

    def flush() -> None:
        if name is None:
            return
        if pending is None:
            raise _BindingInvalid
        content, _newline_width = pending
        if not content or len(content) > line_bases:
            raise _BindingInvalid
        rows.append((name, length + len(content), first_offset, line_bases, line_width))

    for content, line_offset, newline_width in lines:
        if content.startswith(b">"):
            flush()
            header = content[1:].strip().split(None, 1)
            if not header:
                raise _BindingInvalid
            try:
                name = header[0].decode("utf-8")
            except UnicodeDecodeError as exc:
                raise _BindingInvalid from exc
            length = 0
            first_offset = 0
            line_bases = 0
            line_width = 0
            pending = None
            continue
        if name is None or not content:
            raise _BindingInvalid
        if pending is None:
            first_offset = line_offset
            line_bases = len(content)
            line_width = line_bases + (newline_width if newline_width else 1)
            pending = (content, newline_width)
            continue
        previous, previous_newline = pending
        if len(previous) != line_bases:
            raise _BindingInvalid
        if previous_newline != line_width - line_bases:
            raise _BindingInvalid
        length += len(previous)
        pending = (content, newline_width)
    flush()
    if not rows or any(row_length <= 0 for _, row_length, *_ in rows):
        raise _BindingInvalid
    return rows


def _parse_chrom_sizes_lines(lines: Iterable[bytes]) -> dict[str, int]:
    rows = []
    for raw in lines:
        line = raw.decode("utf-8")
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) < 2:
            raise _BindingInvalid
        try:
            length = int(fields[1])
        except ValueError as exc:
            raise _BindingInvalid from exc
        if not fields[0] or length <= 0:
            raise _BindingInvalid
        rows.append((fields[0], length))
    result = dict(rows)
    if not result or len(result) != len(rows):
        raise _BindingInvalid
    return result


def _parse_contig_column_lines(lines: Iterable[bytes], kind: str) -> set[str]:
    contigs = set()
    for raw in lines:
        line = raw.decode("utf-8")
        if not line.strip() or line.startswith("#"):
            continue
        if kind == "gtf":
            fields = line.split("\t")
            if len(fields) < 9:
                continue
        else:
            if re.match(r"^(track|browser)(\s|$)", line):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                continue
        contigs.add(fields[0])
    return contigs


def _parse_dict_contigs_lines(lines: Iterable[bytes]) -> dict[str, int]:
    contigs = {}
    for raw in lines:
        line = raw.decode("utf-8")
        if not line.startswith("@SQ"):
            continue
        fields = dict(
            field.split(":", 1) for field in line.split("\t")[1:] if ":" in field
        )
        name = fields.get("SN")
        length = fields.get("LN")
        if not name or not length:
            raise _BindingInvalid
        try:
            value = int(length)
        except ValueError as exc:
            raise _BindingInvalid from exc
        if value <= 0 or name in contigs:
            raise _BindingInvalid
        contigs[name] = value
    if not contigs:
        raise _BindingInvalid
    return contigs


def _verify_v2_contig_closure(
    resources: Mapping[str, VerifiedLocalFile],
    *,
    organelle_chroms: list[str],
    canonical_chroms: list[str],
) -> None:
    parsed = {
        "reference_fasta": _parse_fasta_layout(
            _iter_verified_raw_lines(
                resources["reference_fasta"].path,
                resources["reference_fasta"].sha256,
            )
        ),
        "reference_fai": _parse_fai_lines(
            _iter_verified_lines(
                resources["reference_fai"].path,
                resources["reference_fai"].sha256,
            )
        ),
        "chrom_sizes": _parse_chrom_sizes_lines(
            _iter_verified_lines(
                resources["chrom_sizes"].path,
                resources["chrom_sizes"].sha256,
            )
        ),
        "gtf": _parse_contig_column_lines(
            _iter_verified_lines(resources["gtf"].path, resources["gtf"].sha256),
            "gtf",
        ),
        "blacklist": _parse_contig_column_lines(
            _iter_verified_lines(
                resources["blacklist"].path, resources["blacklist"].sha256
            ),
            "bed",
        ),
    }
    fasta_rows = parsed["reference_fasta"]
    fai_rows = parsed["reference_fai"]
    if fasta_rows != fai_rows:
        raise _BindingInvalid
    lengths = {name: length for name, length, *_ in fai_rows}
    if len(lengths) != len(fai_rows):
        raise _BindingInvalid
    if parsed["chrom_sizes"] != lengths:
        raise _BindingInvalid
    fasta_names = set(lengths)
    for name in ("gtf", "blacklist"):
        if parsed[name] - fasta_names:
            raise _BindingInvalid
    if "reference_dict" in resources:
        dict_contigs = _parse_dict_contigs_lines(
            _iter_verified_lines(
                resources["reference_dict"].path,
                resources["reference_dict"].sha256,
            )
        )
        if dict_contigs != lengths:
            raise _BindingInvalid
    for chrom in organelle_chroms + canonical_chroms:
        if chrom not in fasta_names:
            raise _BindingInvalid


def _verify_v2_egs(assembly: str, effective_genome_size: object) -> None:
    if isinstance(effective_genome_size, bool) or not (
        (isinstance(effective_genome_size, int) and effective_genome_size > 0)
        or effective_genome_size in {"hs", "mm"}
    ):
        raise _BindingInvalid
    if effective_genome_size == "hs" and not re.match(r"^(GRCh|hg)", assembly):
        raise _BindingInvalid
    if effective_genome_size == "mm" and not re.match(r"^(GRCm|mm)", assembly):
        raise _BindingInvalid


def _verify_binding_v2(
    payload: Mapping[str, object],
) -> tuple[dict[str, Any], AdapterReferenceBindingIdentity]:
    required = _TOP_LEVEL_FIELDS
    allowed = required | _V2_OPTIONAL_FIELDS
    if not required.issubset(payload) or not set(payload).issubset(allowed):
        raise _BindingInvalid
    if payload.get("schema_version") != ENCODE_REFERENCE_BINDING_CONTRACT_V2:
        raise _BindingInvalid
    assembly = payload.get("assembly")
    if not isinstance(assembly, str) or _ASSEMBLY.fullmatch(assembly) is None:
        raise _BindingInvalid
    effective_genome_size = payload.get("effective_genome_size")
    _verify_v2_egs(assembly, effective_genome_size)

    resource_payload = payload.get("genome_resources")
    if not isinstance(resource_payload, Mapping):
        raise _BindingInvalid
    if not set(_V2_REQUIRED_RESOURCES).issubset(resource_payload):
        raise _BindingInvalid
    if not set(resource_payload).issubset(_V2_RESOURCE_NAMES):
        raise _BindingInvalid
    resources: dict[str, VerifiedLocalFile] = {}
    for name in _V2_REQUIRED_RESOURCES:
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
    if "reference_dict" in resource_payload:
        descriptor = resource_payload.get("reference_dict")
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
        resources["reference_dict"] = result.value

    organelle_chroms = _binding_chrom_list(
        payload.get("organelle_chroms", []), "organelle_chroms"
    )
    canonical_chroms = _binding_chrom_list(
        payload.get("canonical_chroms", []), "canonical_chroms"
    )
    _verify_v2_contig_closure(
        resources,
        organelle_chroms=organelle_chroms,
        canonical_chroms=canonical_chroms,
    )

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
    _frame(digest, ENCODE_REFERENCE_BINDING_CONTRACT_V2.encode())
    _frame(digest, assembly.encode())
    _frame(digest, str(effective_genome_size).encode())
    for name in _V2_RESOURCE_NAMES:
        resource = resources.get(name)
        if resource is None:
            continue
        _frame(digest, name.encode())
        _frame(digest, str(resource.size_bytes).encode())
        _frame(digest, bytes.fromhex(resource.sha256))
    for chrom in organelle_chroms:
        _frame(digest, b"organelle")
        _frame(digest, chrom.encode())
    for chrom in canonical_chroms:
        _frame(digest, b"canonical")
        _frame(digest, chrom.encode())
    for suffix in suffixes:
        file = index_files[suffix]
        _frame(digest, suffix.encode())
        _frame(digest, str(file.size_bytes).encode())
        _frame(digest, bytes.fromhex(file.sha256))
    identity = AdapterReferenceBindingIdentity(
        workflow_id=_WORKFLOW_ID,
        contract_version=ENCODE_REFERENCE_BINDING_CONTRACT_V2,
        identity_scheme=ADAPTER_REFERENCE_BINDING_IDENTITY_SCHEME,
        identity_sha256=digest.hexdigest(),
    )
    return (
        {
            "assembly": assembly,
            "effective_genome_size": effective_genome_size,
            "resources": resources,
            "index_prefix": prefix,
            "organelle_chroms": organelle_chroms,
            "canonical_chroms": canonical_chroms,
        },
        identity,
    )


def _binding_chrom_list(value: object, field: str) -> list[str]:
    if value is None:
        return []
    if not isinstance(value, (list, tuple)):
        raise _BindingInvalid
    result = []
    for item in value:
        if not isinstance(item, str) or not item or len(item.encode()) > 1024:
            raise _BindingInvalid
        result.append(item)
    if len(set(result)) != len(result):
        raise _BindingInvalid
    return list(result)


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
