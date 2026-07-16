"""Secure resource closure tests for the bulk RNA-seq runtime boundary."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import pytest

from encode_pipeline.adapters.bulk_rnaseq import resource_closure
from encode_pipeline.adapters.bulk_rnaseq.resource_closure import (
    RIBO_DATABASE_CLOSURE_SCHEME,
    SORTMERNA_INDEX_BINDING_FILENAME,
    SORTMERNA_INDEX_BINDING_SCHEMA_VERSION,
    SORTMERNA_NO_PREBUILT_INDEX_STRATEGY,
    SORTMERNA_VERSION,
    ResourceClosurePolicy,
    safe_regular_file_identity,
    verify_ribo_database_manifest,
    verify_sortmerna_index,
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_database(
    root: Path,
    *,
    filename: str = "db/ref.fa",
    content: bytes = b">rrna\nACGT\n",
) -> tuple[Path, Path]:
    resource = root / filename
    resource.parent.mkdir(parents=True, exist_ok=True)
    resource.write_bytes(content)
    manifest = root / "database-manifest.txt"
    manifest.write_text(f"{filename}\n", encoding="utf-8")
    return manifest, resource


def _database_closure(root: Path, *, content: bytes = b">rrna\nACGT\n"):
    manifest, _ = _write_database(root, content=content)
    result = verify_ribo_database_manifest(
        manifest,
        expected_manifest_sha256=_sha256(manifest),
    )
    assert result.is_success
    assert result.value is not None
    return result.value


def _write_index(root: Path, database_identity: str) -> str:
    root.mkdir(parents=True)
    index_data = root / "idx" / "ref.idx"
    index_data.parent.mkdir()
    index_data.write_bytes(b"deterministic-index\n")
    binding = {
        "schema_version": SORTMERNA_INDEX_BINDING_SCHEMA_VERSION,
        "database_closure_sha256": database_identity,
        "sortmerna_version": "4.3.7",
        "files": [
            {
                "path": "idx/ref.idx",
                "size_bytes": index_data.stat().st_size,
                "sha256": _sha256(index_data),
            }
        ],
    }
    manifest = root / SORTMERNA_INDEX_BINDING_FILENAME
    manifest.write_text(
        json.dumps(binding, sort_keys=True, separators=(",", ":")) + "\n",
        encoding="utf-8",
    )
    return _sha256(manifest)


def _assert_sanitized_failure(result, *, code: str, submitted: Path) -> None:
    assert result.is_failure
    assert result.errors[0].code == code
    public = repr(result.errors[0].to_dict())
    assert str(submitted) not in public
    assert result.errors[0].technical_message is None
    assert result.errors[0].context == {}


def test_database_manifest_verifies_complete_deterministic_closure(tmp_path: Path):
    first_manifest, first_file = _write_database(tmp_path / "first")
    second_manifest, _ = _write_database(tmp_path / "second")

    first = verify_ribo_database_manifest(
        first_manifest,
        expected_manifest_sha256=_sha256(first_manifest),
    )
    second = verify_ribo_database_manifest(
        second_manifest,
        expected_manifest_sha256=_sha256(second_manifest),
    )

    assert first.is_success and second.is_success
    assert first.value is not None and second.value is not None
    assert first.value.manifest_sha256 == _sha256(first_manifest)
    assert first.value.identity_sha256 == second.value.identity_sha256
    assert len(first.value.identity_sha256) == 64
    assert first.value.files == (
        resource_closure.RiboDatabaseFile(
            manifest_entry="db/ref.fa",
            path=first_file,
            size_bytes=first_file.stat().st_size,
            sha256=_sha256(first_file),
        ),
    )
    route = first.value.no_prebuilt_sortmerna_index
    assert route.database_closure_sha256 == first.value.identity_sha256
    assert route.strategy == SORTMERNA_NO_PREBUILT_INDEX_STRATEGY
    assert route.sortmerna_version == SORTMERNA_VERSION == "4.3.7"
    assert route.deterministic_database_inputs is True
    assert route.deterministic_composition_accepted is True
    assert route.requires_runtime_index_build is True
    assert route.scientific_runtime_accepted is False
    assert RIBO_DATABASE_CLOSURE_SCHEME == "sha256-framed-ribo-database-v1"


@pytest.mark.parametrize(
    "manifest_text",
    [
        "https://example.invalid/rrna.fa\n",
        "/absolute/rrna.fa\n",
        "../rrna.fa\n",
        "db/./rrna.fa\n",
        "db//rrna.fa\n",
        "db/ref.fa\ndb/ref.fa\n",
        " db/ref.fa\n",
        "db\\ref.fa\n",
        "db/ref name.fa\n",
        "db/$(touch-PWN).fa\n",
        "db/`touch-PWN`.fa\n",
        "db/ref;touch-PWN.fa\n",
        "db/ref'quote.fa\n",
        'db/ref"quote.fa\n',
    ],
)
def test_database_manifest_rejects_remote_ambiguous_and_duplicate_entries(
    tmp_path: Path,
    manifest_text: str,
):
    manifest, _ = _write_database(tmp_path)
    manifest.write_text(manifest_text, encoding="utf-8")

    result = verify_ribo_database_manifest(
        manifest,
        expected_manifest_sha256=_sha256(manifest),
    )

    _assert_sanitized_failure(
        result,
        code="BULK_RNASEQ_RIBO_DATABASE_INVALID",
        submitted=manifest,
    )


def test_database_manifest_rejects_malformed_and_oversize_content(tmp_path: Path):
    manifest, _ = _write_database(tmp_path)
    manifest.write_bytes(b"db/ref.fa\n\xff")
    malformed = verify_ribo_database_manifest(
        manifest,
        expected_manifest_sha256=_sha256(manifest),
    )
    assert malformed.is_failure

    manifest.write_bytes(b"db/ref.fa\n")
    oversize = verify_ribo_database_manifest(
        manifest,
        expected_manifest_sha256=_sha256(manifest),
        policy=ResourceClosurePolicy(maximum_manifest_bytes=4),
    )
    assert oversize.is_failure


@pytest.mark.parametrize(
    "entries",
    [
        ("a/rrna.fa", "b/rrna.fa"),
        ("a/rrna.fa.gz", "b/rrna.fa"),
        ("a/rrna.fa", "b/rrna.fasta"),
    ],
)
def test_database_manifest_rejects_staging_and_logical_name_collisions(
    tmp_path: Path,
    entries: tuple[str, str],
):
    for entry in entries:
        resource = tmp_path / entry
        resource.parent.mkdir(parents=True, exist_ok=True)
        resource.write_bytes(b">rrna\nACGT\n")
    manifest = tmp_path / "database-manifest.txt"
    manifest.write_text("\n".join(entries) + "\n", encoding="utf-8")

    result = verify_ribo_database_manifest(
        manifest,
        expected_manifest_sha256=_sha256(manifest),
    )

    _assert_sanitized_failure(
        result,
        code="BULK_RNASEQ_RIBO_DATABASE_INVALID",
        submitted=manifest,
    )


def test_database_manifest_rejects_manifest_and_resource_symlinks(tmp_path: Path):
    manifest, resource = _write_database(tmp_path / "real")
    manifest_link = tmp_path / "manifest-link.txt"
    manifest_link.symlink_to(manifest)
    linked_manifest = verify_ribo_database_manifest(
        manifest_link,
        expected_manifest_sha256=_sha256(manifest),
    )
    assert linked_manifest.is_failure

    resource.unlink()
    external = tmp_path / "external.fa"
    external.write_text(">external\n", encoding="utf-8")
    resource.symlink_to(external)
    linked_resource = verify_ribo_database_manifest(
        manifest,
        expected_manifest_sha256=_sha256(manifest),
    )
    assert linked_resource.is_failure


def test_database_manifest_rejects_fifo_without_blocking(tmp_path: Path):
    manifest, resource = _write_database(tmp_path)
    resource.unlink()
    os.mkfifo(resource)

    result = verify_ribo_database_manifest(
        manifest,
        expected_manifest_sha256=_sha256(manifest),
    )

    assert result.is_failure


def test_database_manifest_detects_path_replacement_race(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    manifest, resource = _write_database(tmp_path)
    original = resource_closure._DescriptorRoot.reopen_identity
    replaced = False

    def replace_before_reopen(self, parts, *, is_directory):
        nonlocal replaced
        if parts == ("db", "ref.fa") and not replaced:
            replacement = resource.with_suffix(".replacement")
            replacement.write_bytes(b">replacement\nTGCA\n")
            os.replace(replacement, resource)
            replaced = True
        return original(self, parts, is_directory=is_directory)

    monkeypatch.setattr(
        resource_closure._DescriptorRoot,
        "reopen_identity",
        replace_before_reopen,
    )
    result = verify_ribo_database_manifest(
        manifest,
        expected_manifest_sha256=_sha256(manifest),
    )

    assert replaced is True
    assert result.is_failure


def test_regular_file_identity_is_bounded_and_fail_closed(tmp_path: Path):
    source = tmp_path / "reads.fastq"
    source.write_bytes(b"@r1\nA\n+\n!\n")

    verified = safe_regular_file_identity(
        source,
        expected_sha256=_sha256(source),
        maximum_bytes=source.stat().st_size,
    )

    assert verified.is_success
    assert verified.value is not None
    assert verified.value.path == source
    assert verified.value.sha256 == _sha256(source)
    assert safe_regular_file_identity(source, maximum_bytes=1).is_failure
    assert safe_regular_file_identity(source, expected_sha256="0" * 64).is_failure
    assert safe_regular_file_identity(Path("/dev/null")).is_failure


def test_regular_file_identity_rejects_symlink_and_sanitizes_issue(tmp_path: Path):
    source = tmp_path / "reads.fastq"
    source.write_bytes(b"reads")
    submitted = tmp_path / "submitted.fastq"
    submitted.symlink_to(source)

    result = safe_regular_file_identity(submitted)

    _assert_sanitized_failure(
        result,
        code="BULK_RNASEQ_LOCAL_FILE_INVALID",
        submitted=submitted,
    )


def test_runtime_resource_paths_must_be_absolute(tmp_path: Path):
    source = tmp_path / "reads.fastq"
    source.write_bytes(b"reads")
    manifest, _resource = _write_database(tmp_path / "database")

    with pytest.MonkeyPatch.context() as monkeypatch:
        monkeypatch.chdir(tmp_path)
        assert safe_regular_file_identity("reads.fastq").is_failure
        result = verify_ribo_database_manifest(
            "database/database-manifest.txt",
            expected_manifest_sha256=_sha256(manifest),
        )

    assert result.is_failure
    assert result.errors[0].code == "BULK_RNASEQ_RIBO_DATABASE_INVALID"


def test_sortmerna_index_verifies_tree_identity_and_database_binding(tmp_path: Path):
    database = _database_closure(tmp_path / "database")
    index = tmp_path / "index"
    manifest_sha256 = _write_index(index, database.identity_sha256)

    result = verify_sortmerna_index(
        index,
        expected_index_sha256=manifest_sha256,
        database_closure=database,
    )

    assert result.is_success
    assert result.value is not None
    assert result.value.path == index
    assert result.value.manifest_sha256 == manifest_sha256
    assert result.value.identity_sha256 != manifest_sha256
    assert len(result.value.identity_sha256) == 64
    assert result.value.database_closure_sha256 == database.identity_sha256
    assert result.value.sortmerna_version == "4.3.7"
    assert [entry.relative_path for entry in result.value.entries] == ["idx/ref.idx"]


def test_sortmerna_index_identity_is_location_independent(tmp_path: Path):
    database = _database_closure(tmp_path / "database")
    first = tmp_path / "first"
    second = tmp_path / "second"
    first_manifest = _write_index(first, database.identity_sha256)
    second_manifest = _write_index(second, database.identity_sha256)

    assert first_manifest == second_manifest
    first_result = verify_sortmerna_index(
        first,
        expected_index_sha256=first_manifest,
        database_closure=database,
    )
    second_result = verify_sortmerna_index(
        second,
        expected_index_sha256=second_manifest,
        database_closure=database,
    )
    assert first_result.is_success and second_result.is_success


def test_sortmerna_index_rejects_database_closure_mismatch(tmp_path: Path):
    expected_database = _database_closure(tmp_path / "expected", content=b">a\nA\n")
    other_database = _database_closure(tmp_path / "other", content=b">b\nT\n")
    index = tmp_path / "index"
    manifest_sha256 = _write_index(index, other_database.identity_sha256)

    result = verify_sortmerna_index(
        index,
        expected_index_sha256=manifest_sha256,
        database_closure=expected_database,
    )

    _assert_sanitized_failure(
        result,
        code="BULK_RNASEQ_SORTMERNA_INDEX_INVALID",
        submitted=index,
    )


def test_sortmerna_index_detects_path_replacement_race(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    database = _database_closure(tmp_path / "database")
    index = tmp_path / "index"
    manifest_sha256 = _write_index(index, database.identity_sha256)
    payload = index / "idx" / "ref.idx"
    original = resource_closure._DescriptorRoot.reopen_identity
    replaced = False

    def replace_before_reopen(self, parts, *, is_directory):
        nonlocal replaced
        if parts == ("idx", "ref.idx") and not replaced:
            replacement = payload.with_suffix(".replacement")
            replacement.write_bytes(b"replacement-index\n")
            os.replace(replacement, payload)
            replaced = True
        return original(self, parts, is_directory=is_directory)

    monkeypatch.setattr(
        resource_closure._DescriptorRoot,
        "reopen_identity",
        replace_before_reopen,
    )
    result = verify_sortmerna_index(
        index,
        expected_index_sha256=manifest_sha256,
        database_closure=database,
    )

    assert replaced is True
    assert result.is_failure


@pytest.mark.parametrize("unsafe_kind", ["symlink", "fifo"])
def test_sortmerna_index_rejects_non_regular_entries(
    tmp_path: Path,
    unsafe_kind: str,
):
    database = _database_closure(tmp_path / "database")
    index = tmp_path / "index"
    manifest_sha256 = _write_index(index, database.identity_sha256)
    payload = index / "idx" / "ref.idx"
    payload.unlink()
    if unsafe_kind == "symlink":
        external = tmp_path / "external.idx"
        external.write_bytes(b"unsafe")
        payload.symlink_to(external)
    else:
        os.mkfifo(payload)

    result = verify_sortmerna_index(
        index,
        expected_index_sha256=manifest_sha256,
        database_closure=database,
    )

    assert result.is_failure


def test_sortmerna_index_rejects_malformed_binding_and_unknown_fields(
    tmp_path: Path,
):
    database = _database_closure(tmp_path / "database")
    index = tmp_path / "index"
    _write_index(index, database.identity_sha256)
    descriptor = index / SORTMERNA_INDEX_BINDING_FILENAME
    descriptor.write_text(
        json.dumps(
            {
                "schema_version": SORTMERNA_INDEX_BINDING_SCHEMA_VERSION,
                "database_closure_sha256": database.identity_sha256,
                "sortmerna_version": "4.3.7",
                "files": [],
                "unexpected": True,
            }
        ),
        encoding="utf-8",
    )

    result = verify_sortmerna_index(
        index,
        expected_index_sha256=_sha256(descriptor),
        database_closure=database,
    )

    assert result.is_failure


def test_sortmerna_index_rejects_non_pinned_sortmerna_version(tmp_path: Path):
    database = _database_closure(tmp_path / "database")
    index = tmp_path / "index"
    _write_index(index, database.identity_sha256)
    descriptor = index / SORTMERNA_INDEX_BINDING_FILENAME
    value = json.loads(descriptor.read_bytes())
    value["sortmerna_version"] = "4.3.8"
    descriptor.write_text(
        json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n",
        encoding="utf-8",
    )

    result = verify_sortmerna_index(
        index,
        expected_index_sha256=_sha256(descriptor),
        database_closure=database,
    )

    assert result.is_failure


@pytest.mark.parametrize("mutation", ["tamper", "extra", "missing"])
def test_sortmerna_index_manifest_enforces_exact_file_tree(
    tmp_path: Path,
    mutation: str,
):
    database = _database_closure(tmp_path / "database")
    index = tmp_path / "index"
    manifest_sha256 = _write_index(index, database.identity_sha256)
    payload = index / "idx/ref.idx"
    if mutation == "tamper":
        payload.write_bytes(b"changed-index\n")
    elif mutation == "extra":
        (index / "idx/extra.idx").write_bytes(b"extra\n")
    else:
        payload.unlink()

    result = verify_sortmerna_index(
        index,
        expected_index_sha256=manifest_sha256,
        database_closure=database,
    )

    assert result.is_failure


def test_sortmerna_index_manifest_rejects_duplicate_file_entries(tmp_path: Path):
    database = _database_closure(tmp_path / "database")
    index = tmp_path / "index"
    _write_index(index, database.identity_sha256)
    descriptor = index / SORTMERNA_INDEX_BINDING_FILENAME
    value = json.loads(descriptor.read_bytes())
    value["files"].append(dict(value["files"][0]))
    descriptor.write_text(
        json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n",
        encoding="utf-8",
    )

    result = verify_sortmerna_index(
        index,
        expected_index_sha256=_sha256(descriptor),
        database_closure=database,
    )

    assert result.is_failure
