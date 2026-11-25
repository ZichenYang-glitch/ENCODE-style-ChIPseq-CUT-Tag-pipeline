"""An immutable archived member can be read and verified a window at a time."""

from concurrent.futures import ThreadPoolExecutor
import hashlib
import os

import pytest

from encode_pipeline.platform.archived_tracks import ArchiveStoreError, CHUNK_BYTES
from encode_pipeline.services.archive_objects import ArchiveObjects


def test_window_checks_only_intersecting_chunks_and_exact_copy(tmp_path, monkeypatch):
    root = tmp_path / "archive"
    root.mkdir()
    content = bytes(range(256)) * (CHUNK_BYTES // 32)
    source = tmp_path / "signal.bw"
    source.write_bytes(content)
    archive = ArchiveObjects(root)
    digest = hashlib.sha256(content).hexdigest()
    with source.open("rb") as opened:
        member = archive.seal_fd(opened.fileno(), digest)
    source.unlink()
    reads = []
    real_pread = os.pread

    def metered(fd, size, offset):
        reads.append((size, offset))
        return real_pread(fd, size, offset)

    monkeypatch.setattr(os, "pread", metered)
    with archive.open_member(member) as reader:
        assert (
            reader.read(CHUNK_BYTES - 8, CHUNK_BYTES + 11)
            == content[CHUNK_BYTES - 8 : CHUNK_BYTES + 11]
        )
        assert reader.read(2, 2) == b""
    assert reads == [(CHUNK_BYTES, 0), (CHUNK_BYTES, CHUNK_BYTES)]


def test_same_size_replacement_denied_even_outside_changed_window(tmp_path):
    root = tmp_path / "archive"
    root.mkdir()
    source = tmp_path / "signal.bw"
    source.write_bytes(b"a" * CHUNK_BYTES + b"x")
    store = ArchiveObjects(root)
    with source.open("rb") as fd:
        version = store.seal_fd(fd.fileno())
    with store.open_member(version) as reader:
        other = tmp_path / "other"
        other.write_bytes(b"a" * CHUNK_BYTES + b"y")
        other.replace(root / version.sha256)
        with pytest.raises(ArchiveStoreError, match="archive_integrity"):
            reader.read(0, 1)
        with pytest.raises(ArchiveStoreError, match="archive_integrity"):
            store.open_member(version)
    (root / version.sha256).unlink()
    with pytest.raises(ArchiveStoreError, match="archive_offline"):
        store.open_member(version)


def test_retry_and_parallel_publish_reuse_complete_object(tmp_path):
    source = tmp_path / "s"
    source.write_bytes(b"same exact bytes")
    root = tmp_path / "store"
    root.mkdir()
    store = ArchiveObjects(root)

    def publish(_):
        with source.open("rb") as fd:
            return store.seal_fd(fd.fileno())

    with ThreadPoolExecutor(4) as pool:
        items = list(pool.map(publish, range(5)))
    assert len(set(items)) == 1
    assert [p.name for p in root.iterdir()] == [items[0].sha256]
    with (
        source.open("rb") as fd,
        pytest.raises(ArchiveStoreError, match="archive_integrity"),
    ):
        store.seal_fd(fd.fileno(), "0" * 64)
    assert [p.name for p in root.iterdir()] == [items[0].sha256]


def test_partial_failure_and_unsafe_storage_never_return_a_different_member(tmp_path):
    root = tmp_path / "archive"
    root.mkdir()
    source = tmp_path / "fixed"
    source.write_bytes(b"genuine bytes")
    store = ArchiveObjects(root)
    with source.open("rb") as fd:
        with pytest.raises(ArchiveStoreError, match="archive_integrity"):
            store.seal_fd(fd.fileno(), "c" * 64)
        assert not list(root.iterdir())
        pinned = store.seal_fd(fd.fileno())
    with store.open_member(pinned) as opened:
        with pytest.raises(ArchiveStoreError, match="archive_invalid"):
            opened.read(-1, 1)
    payload = root / pinned.sha256
    payload.rename(root / "previous")
    payload.symlink_to(source)
    with pytest.raises(ArchiveStoreError, match="archive_integrity"):
        store.open_member(pinned)
    with pytest.raises(ArchiveStoreError, match="archive_offline"):
        ArchiveObjects(None).open_member(pinned)


def test_detect_inconsistent_content_hash_before_returning_requested_slice(tmp_path):
    from dataclasses import replace

    root = tmp_path / "archive"
    root.mkdir()
    source = tmp_path / "data"
    source.write_bytes(b"abcd")
    store = ArchiveObjects(root)
    with source.open("rb") as fd:
        pinned = store.seal_fd(fd.fileno())
    with store.open_member(replace(pinned, chunks=("0" * 64,))) as handle:
        with pytest.raises(ArchiveStoreError, match="archive_integrity"):
            handle.read(0, 1)
