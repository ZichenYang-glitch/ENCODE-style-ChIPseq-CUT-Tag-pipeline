"""Private Hi-TrAC qualification admission; no platform registry integration.

The runtime pin is repository-owned. Reference manifests are approved build
records: their digest is supplied separately, never inferred from the manifest.
Only staged, verified copies are passed to the upstream shell-based script.
"""

from __future__ import annotations

from dataclasses import dataclass
import gzip
import hashlib
from itertools import zip_longest
import json
import os
from pathlib import Path
import re
import stat
from typing import Iterator
import zlib


INDEX_SUFFIXES = (".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2")
SCRIPT_SHA256 = "c3c4ef4e6287fa4ade97a6f5980d20c81ea345b8e8e13c88ae3b7c4bd3a67aec"
_COMMIT = "de6cc732fa00b408551b9f4272933640c08447f1"
_SAFE_PATH = re.compile(r"/[A-Za-z0-9_./-]+\Z")
_SAMPLE = re.compile(r"[A-Za-z0-9][A-Za-z0-9_. -]{0,127}\Z")
_SEQUENCE = re.compile(r"[ACGTRYSWKMBDHVNacgtryswkmbdhvn]+\Z")
_SHA = re.compile(r"[0-9a-f]{64}\Z")


class AdmissionError(ValueError):
    """A controlled diagnostic; source paths and payloads are never in str()."""

    def __init__(self, reason_code: str, sample: str | None = None):
        super().__init__(reason_code)
        self.reason_code = reason_code
        self.sample = sample


@dataclass(frozen=True)
class Sample:
    id: str
    r1: Path
    r2: Path


@dataclass(frozen=True)
class RuntimeBinding:
    prefix: Path
    python: Path
    script: Path
    tools: dict[str, Path]
    lock_sha256: str
    binding_sha256: str


@dataclass(frozen=True)
class ReferenceBinding:
    fasta: Path
    prefix: Path
    files: dict[str, str]
    contigs: dict[str, int]
    binding_sha256: str


@dataclass(frozen=True)
class PreparedInputs:
    staged_fastq_dir: Path
    reference_prefix: Path
    samples: dict[str, dict]
    contigs: dict[str, int]
    input_identity: dict


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _json(path: Path) -> dict:
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise AdmissionError("binding_duplicate_key")
            result[key] = value
        return result

    try:
        value = json.loads(path.read_text(encoding="utf-8"), object_pairs_hook=pairs)
    except (OSError, ValueError, UnicodeError) as exc:
        if isinstance(exc, AdmissionError):
            raise
        raise AdmissionError("binding_invalid") from None
    if not isinstance(value, dict):
        raise AdmissionError("binding_invalid")
    return value


def _regular(path: Path) -> Path:
    """Reject links on all input-path components, including the final file."""
    if not path.is_absolute() or ".." in path.parts:
        raise AdmissionError("input_path_invalid")
    try:
        current = Path(path.anchor)
        mode = current.lstat().st_mode
        for part in path.parts[1:]:
            current /= part
            mode = current.lstat().st_mode
            if stat.S_ISLNK(mode):
                raise AdmissionError("input_symlink")
        if not stat.S_ISREG(mode):
            raise AdmissionError("input_not_regular")
    except OSError:
        raise AdmissionError("input_missing") from None
    return path


def _expected_digest(value: object) -> str:
    if not isinstance(value, str) or not _SHA.fullmatch(value):
        raise AdmissionError("identity_invalid")
    return value


def _normalized_digest(path: Path, prefix: Path) -> str:
    return hashlib.sha256(
        path.read_bytes().replace(os.fsencode(prefix), b"${PREFIX}")
    ).hexdigest()


def _package_digest(prefix: Path, metadata: dict) -> tuple[str, int]:
    """Hash path names and installed bytes; metadata cannot hide omitted files."""
    files = metadata.get("files")
    if not isinstance(files, list) or any(not isinstance(p, str) for p in files):
        raise AdmissionError("runtime_package_invalid")
    selected = []
    for name in files:
        relative = Path(name)
        if relative.is_absolute() or ".." in relative.parts:
            raise AdmissionError("runtime_package_invalid")
        if name.endswith((".pyc", ".pyo")) or "__pycache__" in relative.parts:
            continue
        selected.append(name)
    if len(set(selected)) != len(selected):
        raise AdmissionError("runtime_package_invalid")
    digest = hashlib.sha256()
    try:
        for name in sorted(selected):
            path = prefix / name
            if not path.resolve().is_relative_to(prefix):
                raise AdmissionError("runtime_link_escape")
            if path.is_symlink():
                # Library aliases are part of the package; external targets fail.
                if not path.resolve().is_relative_to(prefix):
                    raise AdmissionError("runtime_link_escape")
                content = b"L" + os.fsencode(os.readlink(path)).replace(
                    os.fsencode(prefix), b"${PREFIX}"
                )
            elif path.is_file():
                content = b"F" + path.read_bytes().replace(
                    os.fsencode(prefix), b"${PREFIX}"
                )
            else:
                raise AdmissionError("runtime_file_missing")
            digest.update(
                name.encode()
                + b"\0"
                + hashlib.sha256(content).hexdigest().encode()
                + b"\0"
            )
    except OSError:
        raise AdmissionError("runtime_file_unreadable") from None
    return digest.hexdigest(), len(selected)


def load_runtime_binding(path: Path) -> RuntimeBinding:
    """Validate a relocated binding against the repository's fixed build pins."""
    lock_path = (
        Path(__file__).resolve().parents[4] / "config/hitrac_preprocess/tools.lock.json"
    )
    lock = _json(lock_path)
    data = _json(_regular(path))
    try:
        prefix = Path(data["prefix"])
        script = Path(data["script"])
    except (KeyError, TypeError):
        raise AdmissionError("runtime_binding_invalid") from None
    if data.get("schema_version") != "hitrac-runtime-binding-v1":
        raise AdmissionError("runtime_binding_invalid")
    if data.get("lock_sha256") != sha256_file(lock_path):
        raise AdmissionError("runtime_lock_mismatch")
    if (
        not prefix.is_absolute()
        or prefix.resolve() != prefix
        or not _SAFE_PATH.fullmatch(str(prefix))
    ):
        raise AdmissionError("runtime_path_invalid")
    if (
        lock["upstream"]["commit"] != _COMMIT
        or lock["upstream"]["script_sha256"] != SCRIPT_SHA256
    ):
        raise AdmissionError("runtime_lock_mismatch")
    if sha256_file(_regular(script)) != SCRIPT_SHA256:
        raise AdmissionError("upstream_script_changed")
    if not _SAFE_PATH.fullmatch(str(script)):
        raise AdmissionError("runtime_path_invalid")
    expected = {p["name"]: p for p in lock["packages"]}
    seen = set()
    for metadata_path in (prefix / "conda-meta").glob("*.json"):
        metadata = _json(_regular(metadata_path))
        name = metadata.get("name")
        if name not in expected or name in seen:
            raise AdmissionError("runtime_package_set_mismatch")
        seen.add(name)
        pin = expected[name]
        if any(metadata.get(key) != pin[key] for key in ("version", "build", "sha256")):
            raise AdmissionError("runtime_package_identity_mismatch")
        if _package_digest(prefix, metadata) != (
            pin["installed_files_sha256"],
            pin["installed_file_count"],
        ):
            raise AdmissionError("runtime_installed_bytes_changed")
    if seen != set(expected):
        raise AdmissionError("runtime_package_set_mismatch")
    package = prefix / lock["python_site_packages"] / "cLoops2"
    actual_sources = {p.name for p in package.iterdir() if p.name != "__pycache__"}
    if actual_sources != set(lock["upstream"]["package_sources"]):
        raise AdmissionError("upstream_package_changed")
    for name, digest in lock["upstream"]["package_sources"].items():
        if sha256_file(_regular(package / name)) != digest:
            raise AdmissionError("upstream_package_changed")
    tools = {}
    for name, pin in lock["entrypoints"].items():
        entry = Path(pin["path"]) if "path" in pin else prefix / pin["relative_path"]
        # Preserve invocation spelling: bamToBed and cLoops2 are real CLI aliases.
        if "relative_path" in pin and not entry.resolve().is_relative_to(prefix):
            raise AdmissionError("runtime_link_escape")
        if not entry.is_file() or not os.access(entry, os.X_OK):
            raise AdmissionError("runtime_tool_missing")
        if _normalized_digest(entry, prefix) != pin["sha256"]:
            raise AdmissionError("runtime_tool_changed")
        tools[name] = entry
    return RuntimeBinding(
        prefix,
        tools["python"],
        script,
        tools,
        sha256_file(lock_path),
        sha256_file(path),
    )


def _fasta_contigs(path: Path) -> dict[str, int]:
    result: dict[str, int] = {}
    current = None
    try:
        with path.open(encoding="ascii") as handle:
            for raw in handle:
                line = raw.rstrip("\r\n")
                if line.startswith(">"):
                    fields = line[1:].split()
                    if not fields or fields[0] in result:
                        raise AdmissionError("reference_fasta_invalid")
                    current = fields[0]
                    result[current] = 0
                elif current is not None and _SEQUENCE.fullmatch(line):
                    result[current] += len(line)
                else:
                    raise AdmissionError("reference_fasta_invalid")
    except (OSError, UnicodeError):
        raise AdmissionError("reference_fasta_invalid") from None
    if not result or any(length <= 0 for length in result.values()):
        raise AdmissionError("reference_fasta_invalid")
    return result


def load_reference_binding(path: Path, expected_sha256: str) -> ReferenceBinding:
    """Admit six small-index files tied to an externally approved build record."""
    if sha256_file(_regular(path)) != _expected_digest(expected_sha256):
        raise AdmissionError("reference_binding_identity_mismatch")
    data = _json(path)
    try:
        fasta = Path(data["fasta"])
        prefix = Path(data["prefix"])
        files = data["files"]
        contigs = data["contigs"]
    except (KeyError, TypeError):
        raise AdmissionError("reference_binding_invalid") from None
    if data.get("schema_version") != "hitrac-reference-binding-v1" or not isinstance(
        files, dict
    ):
        raise AdmissionError("reference_binding_invalid")
    if set(files) != {"fasta", *INDEX_SUFFIXES} or not prefix.is_absolute():
        raise AdmissionError("reference_index_set_invalid")
    for name, expected in files.items():
        source = fasta if name == "fasta" else Path(f"{prefix}{name}")
        if sha256_file(_regular(source)) != _expected_digest(expected):
            raise AdmissionError("reference_file_identity_mismatch")
        if source.stat().st_size == 0:
            raise AdmissionError("reference_file_empty")
    actual_contigs = _fasta_contigs(fasta)
    if not isinstance(contigs, dict) or contigs != actual_contigs:
        raise AdmissionError("reference_contigs_mismatch")
    return ReferenceBinding(fasta, prefix, files, actual_contigs, expected_sha256)


def _read_id(header: str, mate: int) -> str:
    fields = header[1:].split()
    if not fields:
        raise AdmissionError("fastq_header_invalid")
    name = fields[0]
    if name.endswith(("/1", "/2")):
        if name[-1] != str(mate):
            raise AdmissionError("fastq_mate_role_invalid")
        name = name[:-2]
    if len(fields) > 1 and re.match(r"[0-9]+:", fields[1]):
        if not fields[1].startswith(f"{mate}:"):
            raise AdmissionError("fastq_mate_role_invalid")
    if not name or any(ord(char) < 33 or ord(char) > 126 for char in name):
        raise AdmissionError("fastq_header_invalid")
    return name


def _records(path: Path, mate: int) -> Iterator[str]:
    try:
        _regular(path)
        with path.open("rb") as compressed:
            if compressed.read(2) != b"\x1f\x8b":
                raise AdmissionError("fastq_gzip_invalid")
        with gzip.open(path, "rt", encoding="ascii", newline="") as handle:
            while True:
                header = handle.readline()
                if header == "":
                    return
                sequence, plus, quality = (handle.readline() for _ in range(3))
                if not all((sequence, plus, quality)):
                    raise AdmissionError("fastq_truncated_record")
                lines = [
                    line.rstrip("\r\n") for line in (header, sequence, plus, quality)
                ]
                header, sequence, plus, quality = lines
                if not header.startswith("@") or not plus.startswith("+"):
                    raise AdmissionError("fastq_structure_invalid")
                if not _SEQUENCE.fullmatch(sequence) or len(sequence) != len(quality):
                    raise AdmissionError("fastq_sequence_quality_invalid")
                if any(ord(char) < 33 or ord(char) > 126 for char in quality):
                    raise AdmissionError("fastq_quality_invalid")
                if plus[1:] and plus[1:] != header[1:]:
                    raise AdmissionError("fastq_plus_header_mismatch")
                yield _read_id(header, mate)
    except (OSError, EOFError, UnicodeError, zlib.error):
        raise AdmissionError("fastq_gzip_invalid") from None


def validate_fastq_pair(r1: Path, r2: Path) -> int:
    """Stream strict four-line gzip FASTQ records; allow zero for policy checks."""
    count = 0
    for left, right in zip_longest(_records(r1, 1), _records(r2, 2)):
        if left is None or right is None:
            raise AdmissionError("fastq_pair_count_mismatch")
        if left != right:
            raise AdmissionError("fastq_pair_id_mismatch")
        count += 1
    return count


def _copy_checked(source: Path, destination: Path, expected: str | None = None) -> str:
    _regular(source)
    digest = hashlib.sha256()
    before = source.stat()
    with source.open("rb") as original, destination.open("xb") as staged:
        for block in iter(lambda: original.read(1024 * 1024), b""):
            digest.update(block)
            staged.write(block)
    after = source.stat()
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
        raise AdmissionError("input_changed_during_staging")
    actual = digest.hexdigest()
    if expected is not None and actual != expected:
        raise AdmissionError("reference_changed_during_staging")
    destination.chmod(0o400)
    return actual


def prepare_inputs(
    samples: list[Sample],
    runtime: RuntimeBinding,
    reference: ReferenceBinding,
    attempt: Path,
) -> PreparedInputs:
    """The caller owns a new 0700 attempt; never reuse staging or alter sources."""
    if not samples or len(samples) > 999999:
        raise AdmissionError("sample_set_invalid")
    if (
        not attempt.is_absolute()
        or attempt.resolve() != attempt
        or not _SAFE_PATH.fullmatch(str(attempt))
    ):
        raise AdmissionError("attempt_path_invalid")
    if not attempt.is_dir() or stat.S_IMODE(attempt.stat().st_mode) != 0o700:
        raise AdmissionError("attempt_permissions_invalid")
    ids = set()
    for sample in samples:
        if (
            not isinstance(sample.id, str)
            or not _SAMPLE.fullmatch(sample.id)
            or sample.id in ids
        ):
            raise AdmissionError("sample_id_invalid")
        ids.add(sample.id)
    fastq = attempt / "input"
    reference_dir = attempt / "reference"
    try:
        fastq.mkdir(mode=0o700)
        reference_dir.mkdir(mode=0o700)
    except FileExistsError:
        raise AdmissionError("attempt_staging_exists") from None
    prepared = {}
    for position, sample in enumerate(samples, 1):
        token = f"s{position:06d}"
        digests = {}
        try:
            for mate, source in ((1, sample.r1), (2, sample.r2)):
                digests[f"r{mate}"] = _copy_checked(
                    source, fastq / f"{token}_R{mate}.fastq.gz"
                )
            count = validate_fastq_pair(
                fastq / f"{token}_R1.fastq.gz", fastq / f"{token}_R2.fastq.gz"
            )
        except AdmissionError as exc:
            raise AdmissionError(exc.reason_code, token) from None
        prepared[token] = {
            "display_id": sample.id,
            "raw_pairs": count,
            "sha256": digests,
        }
    prefix = reference_dir / "genome"
    for suffix, digest in reference.files.items():
        source = (
            reference.fasta
            if suffix == "fasta"
            else Path(f"{reference.prefix}{suffix}")
        )
        destination = (
            reference_dir / "reference.fa"
            if suffix == "fasta"
            else Path(f"{prefix}{suffix}")
        )
        _copy_checked(source, destination, digest)
    identity = {
        "runtime_binding_sha256": runtime.binding_sha256,
        "runtime_lock_sha256": runtime.lock_sha256,
        "reference_binding_sha256": reference.binding_sha256,
        "samples": prepared,
    }
    return PreparedInputs(fastq, prefix, prepared, reference.contigs, identity)
