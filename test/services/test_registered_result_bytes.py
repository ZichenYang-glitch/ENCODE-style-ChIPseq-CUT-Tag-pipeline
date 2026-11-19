"""Failure branches for the private response copy."""

from __future__ import annotations

from hashlib import sha256
from pathlib import Path
from unittest.mock import patch

import pytest

from encode_pipeline.services.registered_result_bytes import (
    RegisteredResultBytesError,
    seal_private_copy,
)


def _seal(tmp_path: Path, *, expected_sha256: str | None = None):
    root = tmp_path / "root"
    root.mkdir()
    payload = b"synthetic-A\t3\n"
    (root / "fixture.tsv").write_bytes(payload)
    response_root = tmp_path / "responses"
    response_root.mkdir()
    return seal_private_copy(
        root=root,
        relative_path="fixture.tsv",
        expected_size=len(payload),
        expected_sha256=expected_sha256 or sha256(payload).hexdigest(),
        response_root=response_root,
    )


def test_short_write_becomes_an_unavailable_copy(tmp_path: Path) -> None:
    with (
        patch(
            "encode_pipeline.services.registered_result_bytes.os.write",
            return_value=0,
        ),
        pytest.raises(RegisteredResultBytesError),
    ):
        _seal(tmp_path)
    assert list((tmp_path / "responses").iterdir()) == []
