"""An immutable archived member can be read and verified a window at a time."""

from concurrent.futures import ThreadPoolExecutor
import errno
import hashlib
import os
import stat

import pytest

from encode_pipeline.platform.archived_tracks import ArchiveStoreError, CHUNK_BYTES
from encode_pipeline.services import archive_objects
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


def test_crash_after_atomic_publish_leaves_one_link_and_retry_reuses_it(
    tmp_path, monkeypatch
):
    content = b"atomic publish"
    root = tmp_path / "archive"
    root.mkdir()
    source = tmp_path / "source"
    source.write_bytes(content)
    digest = hashlib.sha256(content).hexdigest()
    store = ArchiveObjects(root)
    real_publish = archive_objects._rename_noreplace
    state = {"armed": True}

    def crash_after_publish(directory, source_name, destination):
        if state["armed"]:
            state["armed"] = False
            real_publish(directory, source_name, destination)
            raise OSError(errno.EIO, "killed at atomic publish")
        return real_publish(directory, source_name, destination)

    def forbid_link(*_args, **_kwargs):
        raise AssertionError("publish must be one rename, not a hard link")

    monkeypatch.setattr(archive_objects, "_rename_noreplace", crash_after_publish)
    monkeypatch.setattr(os, "link", forbid_link)
    with source.open("rb") as fd:
        with pytest.raises(ArchiveStoreError, match="archive_offline"):
            store.seal_fd(fd.fileno())
    published = root / digest
    assert published.is_file()
    assert not published.is_symlink()
    info = published.stat()
    assert info.st_nlink == 1
    assert stat.S_IMODE(info.st_mode) == 0o440
    assert published.read_bytes() == content
    assert [path.name for path in root.iterdir()] == [digest]
    with source.open("rb") as fd:
        member = store.seal_fd(fd.fileno())
    assert member.sha256 == digest
    assert (root / digest).stat().st_nlink == 1
    with store.open_member(member) as reader:
        assert reader.read(0, len(content)) == content


def test_failed_atomic_publish_retries_without_a_second_link(tmp_path, monkeypatch):
    content = b"retry after crash"
    root = tmp_path / "archive"
    root.mkdir()
    orphan = root / ".sealing-orphan"
    orphan.write_bytes(content)
    source = tmp_path / "source"
    source.write_bytes(content)
    digest = hashlib.sha256(content).hexdigest()
    store = ArchiveObjects(root)
    real_publish = archive_objects._rename_noreplace
    state = {"armed": True}

    def crash_before_publish(_directory, _source_name, _destination):
        if state["armed"]:
            state["armed"] = False
            raise OSError(errno.EIO, "killed before atomic publish")
        return real_publish(_directory, _source_name, _destination)

    monkeypatch.setattr(archive_objects, "_rename_noreplace", crash_before_publish)
    with source.open("rb") as fd:
        with pytest.raises(ArchiveStoreError, match="archive_offline"):
            store.seal_fd(fd.fileno())
    assert [path.name for path in root.iterdir()] == [orphan.name]
    assert not (root / digest).exists()
    with source.open("rb") as fd:
        member = store.seal_fd(fd.fileno())
    published = root / digest
    assert published.stat().st_nlink == 1
    assert published.stat().st_ino != orphan.stat().st_ino
    assert stat.S_IMODE(published.stat().st_mode) == 0o440
    with store.open_member(member) as reader:
        assert reader.read(0, len(content)) == content


def test_reuse_rejects_permissions_other_than_0440(tmp_path):
    content = b"permission contract"
    digest = hashlib.sha256(content).hexdigest()
    root = tmp_path / "archive"
    root.mkdir()
    source = tmp_path / "source"
    source.write_bytes(content)
    planted = root / digest
    planted.write_bytes(content)
    store = ArchiveObjects(root)
    for mode in (0o400, 0o444, 0o555):
        os.chmod(planted, mode)
        with source.open("rb") as fd:
            with pytest.raises(ArchiveStoreError, match="archive_integrity"):
                store.seal_fd(fd.fileno())
        assert stat.S_IMODE(planted.stat().st_mode) == mode
        assert planted.read_bytes() == content
        assert [path.name for path in root.iterdir()] == [digest]
    os.chmod(planted, 0o440)
    with source.open("rb") as fd:
        member = store.seal_fd(fd.fileno())
    assert member.sha256 == digest
    with store.open_member(member) as reader:
        assert reader.read(0, len(content)) == content


def test_publish_does_not_replace_symlink_directory_or_different_bytes(tmp_path):
    content = b"sealed-bytes"
    other = b"other-bytes!"
    digest = hashlib.sha256(content).hexdigest()
    root = tmp_path / "archive"
    root.mkdir()
    source = tmp_path / "source"
    source.write_bytes(content)
    outside = tmp_path / "outside"
    outside.write_bytes(content)
    store = ArchiveObjects(root)
    link = root / digest
    link.symlink_to(outside)
    with source.open("rb") as fd:
        with pytest.raises(ArchiveStoreError, match="archive_integrity"):
            store.seal_fd(fd.fileno())
    assert link.is_symlink()
    assert os.readlink(link) == str(outside)
    assert [path.name for path in root.iterdir()] == [digest]

    link.unlink()
    link.mkdir()
    with source.open("rb") as fd:
        with pytest.raises(ArchiveStoreError, match="archive_integrity"):
            store.seal_fd(fd.fileno())
    assert link.is_dir()
    assert not any(path.name.startswith(".sealing-") for path in root.iterdir())

    link.rmdir()
    link.write_bytes(other)
    os.chmod(link, 0o440)
    with source.open("rb") as fd:
        with pytest.raises(ArchiveStoreError, match="archive_integrity"):
            store.seal_fd(fd.fileno())
    assert link.is_file()
    assert link.read_bytes() == other
    assert stat.S_IMODE(link.stat().st_mode) == 0o440
    assert [path.name for path in root.iterdir()] == [digest]


def test_both_read_paths_accumulate_short_preads(tmp_path, monkeypatch):
    real_pread = os.pread

    def short(fd, size, offset):
        return real_pread(fd, min(size, 64 * 1024), offset)

    monkeypatch.setattr(os, "pread", short)
    content = b"q" * CHUNK_BYTES + b"tail"
    root = tmp_path / "archive"
    root.mkdir()
    source = tmp_path / "source"
    source.write_bytes(content)
    store = ArchiveObjects(root)
    with source.open("rb") as fd:
        member = store.seal_fd(fd.fileno())
    published = root / member.sha256
    assert published.stat().st_nlink == 1
    assert stat.S_IMODE(published.stat().st_mode) == 0o440
    with store.open_member(member) as reader:
        assert (
            reader.read(CHUNK_BYTES - 2, CHUNK_BYTES + 4)
            == content[CHUNK_BYTES - 2 : CHUNK_BYTES + 4]
        )
        assert reader.read(0, 0) == b""


def test_source_eof_before_declared_size_is_integrity_failure(tmp_path, monkeypatch):
    real_pread = os.pread
    calls = {"n": 0}

    def early_eof(fd, size, offset):
        calls["n"] += 1
        if calls["n"] == 1:
            return real_pread(fd, min(size, 8), offset)
        return b""

    monkeypatch.setattr(os, "pread", early_eof)
    root = tmp_path / "archive"
    root.mkdir()
    source = tmp_path / "source"
    source.write_bytes(b"x" * 100)
    store = ArchiveObjects(root)
    with source.open("rb") as fd:
        with pytest.raises(ArchiveStoreError, match="archive_integrity"):
            store.seal_fd(fd.fileno())
    assert list(root.iterdir()) == []


def test_reader_eof_before_chunk_end_is_integrity_failure(tmp_path, monkeypatch):
    content = b"y" * 100
    root = tmp_path / "archive"
    root.mkdir()
    source = tmp_path / "source"
    source.write_bytes(content)
    store = ArchiveObjects(root)
    with source.open("rb") as fd:
        member = store.seal_fd(fd.fileno())
    real_pread = os.pread
    calls = {"n": 0}

    def early_eof(fd, size, offset):
        calls["n"] += 1
        if calls["n"] == 1:
            return real_pread(fd, min(size, 4), offset)
        return b""

    monkeypatch.setattr(os, "pread", early_eof)
    with store.open_member(member) as reader:
        with pytest.raises(ArchiveStoreError, match="archive_integrity"):
            reader.read(0, member.size_bytes)
    assert (root / member.sha256).read_bytes() == content


def test_closed_reader_and_illegal_bounds_are_domain_errors(tmp_path):
    root = tmp_path / "archive"
    root.mkdir()
    source = tmp_path / "source"
    source.write_bytes(b"abcdef")
    store = ArchiveObjects(root)
    with source.open("rb") as fd:
        member = store.seal_fd(fd.fileno())
    reader = store.open_member(member)
    reader.close()
    with pytest.raises(ArchiveStoreError, match="archive_invalid"):
        reader.read(0, 1)
    with store.open_member(member) as opened:
        with pytest.raises(ArchiveStoreError, match="archive_invalid"):
            opened.read(0.0, 1)
        with pytest.raises(ArchiveStoreError, match="archive_invalid"):
            opened.read(True, True)
        assert opened.read(1, 4) == b"bcd"
