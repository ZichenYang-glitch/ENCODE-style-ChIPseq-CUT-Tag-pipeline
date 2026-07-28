"""Tests for descriptor-owned local InputFile observation and revalidation."""

from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import os

import pytest

from encode_pipeline.platform.input_registry import build_input_file_revision
from encode_pipeline.services.input_file_access import (
    FileObservation,
    InputFileAccess,
    InputFileAccessError,
)


NOW = datetime(2026, 7, 26, 8, 0, tzinfo=timezone.utc)
PROJECT_ID = "prj_" + "1" * 32
POOL_ID = "stgp_" + "2" * 32
INPUT_FILE_ID = "inpf_" + "3" * 32
REVISION_ID = "inpfr_" + "4" * 32


def test_observe_hashes_nested_regular_file_and_keeps_nlink_as_evidence(
    tmp_path,
) -> None:
    root = tmp_path / "pool"
    nested = root / "cohort"
    nested.mkdir(parents=True)
    source = nested / "reads.fastq"
    source.write_bytes(b"ACGT\n")
    os.link(source, nested / "permitted-hardlink.fastq")

    observation = InputFileAccess(root).observe("cohort/reads.fastq")

    assert observation == FileObservation(
        relative_path="cohort/reads.fastq",
        size_bytes=5,
        content_sha256=hashlib.sha256(b"ACGT\n").hexdigest(),
        path_fingerprint=observation.path_fingerprint,
    )
    assert len(observation.path_fingerprint) == 3
    assert observation.path_fingerprint[-1][3] == 2


@pytest.mark.parametrize("symlink_at", ("directory", "leaf", "root"))
def test_observe_rejects_symlinks_at_every_path_component(
    tmp_path,
    symlink_at: str,
) -> None:
    real_root = tmp_path / "real-pool"
    real_nested = real_root / "cohort"
    real_nested.mkdir(parents=True)
    real_file = real_nested / "reads.fastq"
    real_file.write_bytes(b"ACGT\n")
    root = real_root
    relative_path = "cohort/reads.fastq"
    if symlink_at == "directory":
        alternate = real_root / "alternate"
        alternate.mkdir()
        (alternate / "reads.fastq").write_bytes(b"ACGT\n")
        (real_root / "linked").symlink_to(alternate, target_is_directory=True)
        relative_path = "linked/reads.fastq"
    elif symlink_at == "leaf":
        (real_nested / "linked.fastq").symlink_to(real_file)
        relative_path = "cohort/linked.fastq"
    else:
        linked_root = tmp_path / "linked-pool"
        linked_root.symlink_to(real_root, target_is_directory=True)
        root = linked_root

    with pytest.raises(InputFileAccessError):
        InputFileAccess(root).observe(relative_path)


def test_observe_rejects_traversal_and_non_regular_leaf_without_blocking(
    tmp_path,
) -> None:
    root = tmp_path / "pool"
    root.mkdir()
    fifo = root / "reads.pipe"
    os.mkfifo(fifo)
    access = InputFileAccess(root)

    with pytest.raises(InputFileAccessError):
        access.observe("../escape.fastq")
    with pytest.raises(InputFileAccessError):
        access.observe("reads.pipe")


def test_observe_fails_closed_when_reopen_identity_changes(
    tmp_path,
    monkeypatch,
) -> None:
    root = tmp_path / "pool"
    root.mkdir()
    (root / "reads.fastq").write_bytes(b"ACGT\n")
    access = InputFileAccess(root)
    original = access._reopen_fingerprint

    def changed(relative_path: str):
        fingerprint = original(relative_path)
        return (*fingerprint[:-1], (*fingerprint[-1][:-1], fingerprint[-1][-1] + 1))

    monkeypatch.setattr(access, "_reopen_fingerprint", changed)

    with pytest.raises(InputFileAccessError):
        access.observe("reads.fastq")


def test_reverify_checks_exact_revision_size_and_checksum(tmp_path) -> None:
    root = tmp_path / "pool"
    root.mkdir()
    source = root / "reads.fastq"
    source.write_bytes(b"ACGT\n")
    expected = build_input_file_revision(
        input_file_revision_id=REVISION_ID,
        input_file_id=INPUT_FILE_ID,
        project_id=PROJECT_ID,
        storage_pool_id=POOL_ID,
        revision_number=1,
        relative_path="reads.fastq",
        size_bytes=5,
        content_sha256=hashlib.sha256(b"ACGT\n").hexdigest(),
        created_at=NOW,
    )
    access = InputFileAccess(root)

    assert access.reverify(expected).content_sha256 == expected.content_sha256
    source.write_bytes(b"TGCA\n")
    with pytest.raises(InputFileAccessError):
        access.reverify(expected)


def test_failures_do_not_disclose_private_root_or_relative_path(tmp_path) -> None:
    root = tmp_path / "private-customer-name"
    root.mkdir()
    relative_path = "secret/cohort.fastq"

    with pytest.raises(InputFileAccessError) as captured:
        InputFileAccess(root).observe(relative_path)

    rendered = f"{captured.value!s} {captured.value!r}"
    assert str(root) not in rendered
    assert relative_path not in rendered


@pytest.mark.parametrize("ambiguous_kind", ("double-prefix", "empty-component"))
def test_observe_rejects_noncanonical_private_root(
    tmp_path,
    ambiguous_kind: str,
) -> None:
    root = tmp_path / "pool"
    root.mkdir()
    (root / "reads.fastq").write_bytes(b"ACGT\n")
    raw_root = str(root)
    if ambiguous_kind == "double-prefix":
        ambiguous_root = type(root)(f"/{raw_root}")
        assert os.fspath(ambiguous_root).startswith("//")
    else:
        parent, leaf = raw_root.rsplit("/", 1)
        ambiguous_root = f"{parent}//{leaf}"

    with pytest.raises(InputFileAccessError):
        InputFileAccess(ambiguous_root).observe("reads.fastq")
