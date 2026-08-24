"""Contracts for the deterministic two-shard pytest partition."""

from pathlib import Path
import hashlib

import pytest

from scripts.split_deterministic_shards import (
    assign_shards,
    list_test_files,
    main,
    shard_index,
)


REPO_ROOT = Path(__file__).resolve().parents[2]
TEST_ROOT = REPO_ROOT / "test"


def test_repository_partition_is_disjoint_and_exhaustive():
    first, second = assign_shards(TEST_ROOT, 2)
    assert first and second
    assert set(first).isdisjoint(second)
    assert sorted(first + second) == list_test_files(TEST_ROOT)


def test_partition_is_a_pure_function_of_the_tree():
    assert assign_shards(TEST_ROOT, 2) == assign_shards(TEST_ROOT, 2)


def test_assignment_matches_the_documented_sha256_modulo_rule():
    for relative in list_test_files(TEST_ROOT)[:25]:
        expected = (
            int.from_bytes(hashlib.sha256(relative.encode("utf-8")).digest(), "big") % 2
        )
        assert shard_index(relative, 2) == expected


def test_rejects_non_posix_or_empty_paths():
    with pytest.raises(ValueError):
        shard_index("", 2)
    with pytest.raises(ValueError):
        shard_index("nested\\windows.py", 2)


def test_cli_emits_prefixed_sorted_files_for_each_shard(capsys):
    for shard in (1, 2):
        assert main(["--root", "test", "--shard", str(shard), "--of", "2"]) == 0
        out = capsys.readouterr().out.splitlines()
        assert out
        assert out == sorted(out)
        assert all(line.startswith("test/test_") or "/test_" in line for line in out)


def test_cli_rejects_inverted_or_out_of_range_shards(capsys):
    with pytest.raises(SystemExit):
        main(["--root", "test", "--shard", "3", "--of", "2"])
    with pytest.raises(SystemExit):
        main(["--root", "test", "--shard", "0", "--of", "2"])


def test_cli_fails_closed_on_a_tree_without_tests(tmp_path, capsys):
    (tmp_path / "helper.py").write_text("", encoding="utf-8")
    assert main(["--root", str(tmp_path), "--shard", "1", "--of", "2"]) == 2
    assert "no test_*.py" in capsys.readouterr().err


def test_list_test_files_ignores_non_test_and_non_file_entries(tmp_path):
    (tmp_path / "test_alpha.py").write_text("", encoding="utf-8")
    (tmp_path / "notes.py").write_text("", encoding="utf-8")
    nested = tmp_path / "nested"
    nested.mkdir()
    (nested / "test_beta.py").write_text("", encoding="utf-8")
    assert list_test_files(tmp_path) == ["nested/test_beta.py", "test_alpha.py"]
