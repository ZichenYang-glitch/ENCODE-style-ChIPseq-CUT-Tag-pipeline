"""Tests for descriptor-owned local InputFile observation and revalidation."""

from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import os

import stat

import pytest

from encode_pipeline.platform.input_registry import build_input_file_revision
from encode_pipeline.services.input_file_access import (
    FileObservation,
    InputFileAccess,
    InputFileAccessError,
    _open_absolute_leaf,
    watch_leaf_generation,
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


def test_verified_path_rejects_symlink_replacement_after_reverify(tmp_path) -> None:
    root = tmp_path / "pool"
    nested = root / "reads"
    nested.mkdir(parents=True)
    source = nested / "s1.fastq.gz"
    source.write_bytes(b"ACGT\n")
    replacement = nested / "unregistered.fastq.gz"
    replacement.write_bytes(b"TGCA\n")
    expected = build_input_file_revision(
        input_file_revision_id=REVISION_ID,
        input_file_id=INPUT_FILE_ID,
        project_id=PROJECT_ID,
        storage_pool_id=POOL_ID,
        revision_number=1,
        relative_path="reads/s1.fastq.gz",
        size_bytes=5,
        content_sha256=hashlib.sha256(b"ACGT\n").hexdigest(),
        created_at=NOW,
    )
    access = InputFileAccess(root)
    original = InputFileAccess.reverify

    def reverify_then_replace(self, revision):
        observation = original(self, revision)
        source.unlink()
        source.symlink_to(replacement)
        return observation

    InputFileAccess.reverify = reverify_then_replace  # type: ignore[method-assign]
    try:
        with pytest.raises(InputFileAccessError) as captured:
            access.verified_regular_path(expected)
    finally:
        InputFileAccess.reverify = original  # type: ignore[method-assign]

    rendered = f"{captured.value!s} {captured.value!r}"
    assert str(root) not in rendered
    assert "unregistered.fastq.gz" not in rendered
    assert str(captured.value) == "input file is missing, changed, or unsafe"


def test_verified_path_rejects_replacement_after_second_observation(
    tmp_path,
) -> None:
    root = tmp_path / "pool"
    nested = root / "reads"
    nested.mkdir(parents=True)
    source = nested / "s1.fastq.gz"
    payload = b"ACGT\n"
    source.write_bytes(payload)
    replacement = nested / "unregistered.fastq.gz"
    replacement.write_bytes(b"TGCA\n")
    expected = build_input_file_revision(
        input_file_revision_id=REVISION_ID,
        input_file_id=INPUT_FILE_ID,
        project_id=PROJECT_ID,
        storage_pool_id=POOL_ID,
        revision_number=1,
        relative_path="reads/s1.fastq.gz",
        size_bytes=len(payload),
        content_sha256=hashlib.sha256(payload).hexdigest(),
        created_at=NOW,
    )
    access = InputFileAccess(root)
    original = InputFileAccess.observe
    calls = {"count": 0}

    def observe_then_replace(self, relative_path):
        observation = original(self, relative_path)
        calls["count"] += 1
        if calls["count"] == 2 and source.is_file() and not source.is_symlink():
            source.unlink()
            source.symlink_to(replacement)
        return observation

    InputFileAccess.observe = observe_then_replace  # type: ignore[method-assign]
    try:
        with pytest.raises(InputFileAccessError) as captured:
            access.verified_regular_path(expected)
    finally:
        InputFileAccess.observe = original  # type: ignore[method-assign]

    rendered = f"{captured.value!s} {captured.value!r}"
    assert calls["count"] >= 2
    assert str(root) not in rendered
    assert "unregistered.fastq.gz" not in rendered


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


def _revision(relative_path: str, payload: bytes):
    return build_input_file_revision(
        input_file_revision_id=REVISION_ID,
        input_file_id=INPUT_FILE_ID,
        project_id=PROJECT_ID,
        storage_pool_id=POOL_ID,
        revision_number=1,
        relative_path=relative_path,
        size_bytes=len(payload),
        content_sha256=hashlib.sha256(payload).hexdigest(),
        created_at=NOW,
    )


def _held_descriptor_count() -> int:
    return len(os.listdir("/proc/self/fd"))


def _assert_redacted(exc: InputFileAccessError, *hidden: str) -> None:
    rendered = f"{exc!s} {exc!r}"
    assert str(exc) == "input file is missing, changed, or unsafe"
    assert exc.__cause__ is None
    assert exc.__suppress_context__
    for fragment in hidden:
        assert fragment not in rendered


def test_bind_candidate_opens_regular_file(tmp_path) -> None:
    root = tmp_path / "pool"
    root.mkdir()
    payload = b"ACGT\n"
    source = root / "reads.fastq"
    source.write_bytes(payload)
    before = _held_descriptor_count()
    descriptor = InputFileAccess(root).bind_candidate(
        _revision("reads.fastq", payload),
        str(source),
    )
    try:
        assert stat.S_ISREG(os.fstat(descriptor).st_mode)
        assert os.read(descriptor, len(payload) + 1) == payload
    finally:
        os.close(descriptor)
    assert _held_descriptor_count() == before


def test_open_absolute_leaf_missing_parent_is_redacted_without_descriptor_leak(
    tmp_path,
) -> None:
    root = tmp_path / "pool"
    root.mkdir()
    payload = b"ACGT\n"
    source = root / "reads.fastq"
    source.write_bytes(payload)
    candidate = str(tmp_path / "private-parent-name" / "secret-leaf.fastq")
    access = InputFileAccess(root)
    expected = _revision("reads.fastq", payload)

    with pytest.raises(InputFileAccessError) as captured:
        access.bind_candidate(expected, candidate)
    _assert_redacted(
        captured.value, candidate, "private-parent-name", "secret-leaf.fastq"
    )

    before = _held_descriptor_count()
    for _ in range(128):
        with pytest.raises(InputFileAccessError) as captured:
            _open_absolute_leaf(candidate)
        _assert_redacted(
            captured.value,
            candidate,
            "private-parent-name",
            "secret-leaf.fastq",
        )
    assert _held_descriptor_count() == before


def test_open_absolute_leaf_symlink_parent_is_redacted_without_descriptor_leak(
    tmp_path,
) -> None:
    root = tmp_path / "real-parent"
    root.mkdir()
    payload = b"ACGT\n"
    source = root / "secret-leaf.fastq"
    source.write_bytes(payload)
    linked_parent = tmp_path / "linked-parent"
    linked_parent.symlink_to(root, target_is_directory=True)
    candidate = str(linked_parent / "secret-leaf.fastq")
    access = InputFileAccess(root)
    expected = _revision("secret-leaf.fastq", payload)

    with pytest.raises(InputFileAccessError) as captured:
        access.bind_candidate(expected, candidate)
    _assert_redacted(
        captured.value,
        candidate,
        "linked-parent",
        "real-parent",
        "secret-leaf.fastq",
    )

    before = _held_descriptor_count()
    for _ in range(128):
        with pytest.raises(InputFileAccessError) as captured:
            _open_absolute_leaf(candidate)
        _assert_redacted(
            captured.value,
            candidate,
            "linked-parent",
            "real-parent",
            "secret-leaf.fastq",
        )
    assert _held_descriptor_count() == before


def test_open_absolute_leaf_closes_parent_when_leaf_open_fails(tmp_path) -> None:
    parent = tmp_path / "pool"
    parent.mkdir()
    candidate = str(parent / "missing-leaf.fastq")
    before = _held_descriptor_count()
    for _ in range(128):
        with pytest.raises(InputFileAccessError) as captured:
            _open_absolute_leaf(candidate)
        _assert_redacted(captured.value, candidate, str(parent), "missing-leaf.fastq")
    assert _held_descriptor_count() == before


def test_watch_leaf_generation_untouched_leaf_is_unchanged(tmp_path) -> None:
    parent = tmp_path / "pool"
    parent.mkdir()
    leaf = parent / "reads.fastq"
    leaf.write_bytes(b"ACGT\n")
    sibling = parent / "other.fastq"
    before = _held_descriptor_count()
    with watch_leaf_generation(parent, leaf.name) as watch:
        sibling.write_bytes(b"TGCA\n")
        assert watch.changed() is False
    assert watch.changed() is False
    assert _held_descriptor_count() == before


def test_watch_leaf_generation_equal_length_write_is_a_change(tmp_path) -> None:
    parent = tmp_path / "pool"
    parent.mkdir()
    leaf = parent / "reads.fastq"
    payload = b"ACGT\n"
    leaf.write_bytes(payload)
    with watch_leaf_generation(parent, leaf.name) as watch:
        assert watch.changed() is False
        leaf.write_bytes(payload)
        assert leaf.stat().st_size == len(payload)
        assert watch.changed() is True


def test_watch_leaf_generation_restored_inode_and_bytes_still_changed(
    tmp_path,
) -> None:
    parent = tmp_path / "pool"
    parent.mkdir()
    leaf = parent / "reads.fastq"
    payload = b"ACGT\n"
    leaf.write_bytes(payload)
    identity = (leaf.stat().st_dev, leaf.stat().st_ino)
    with watch_leaf_generation(parent, leaf.name) as watch:
        assert watch.changed() is False
        held = parent / "reads.fastq.held"
        leaf.rename(held)
        leaf.write_bytes(b"TGCA\n")
        held.rename(leaf)
        assert watch.changed() is True
    restored = leaf.stat()
    assert (restored.st_dev, restored.st_ino) == identity
    assert leaf.read_bytes() == payload


def test_watch_leaf_generation_fails_closed_when_parent_cannot_be_armed(
    tmp_path,
) -> None:
    missing = tmp_path / "private-parent-name"
    not_directory = tmp_path / "not-a-directory"
    not_directory.write_bytes(b"x")
    before = _held_descriptor_count()
    for _ in range(16):
        for parent in (missing, not_directory):
            with watch_leaf_generation(parent, "secret-leaf.fastq") as watch:
                assert watch.changed() is True
    assert _held_descriptor_count() == before
