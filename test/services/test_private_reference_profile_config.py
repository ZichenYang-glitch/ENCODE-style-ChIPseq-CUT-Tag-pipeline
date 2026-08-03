"""Operator-private Reference Profile configuration tests."""

from __future__ import annotations

import json
import os

import pytest

from encode_pipeline.services.private_reference_profiles import (
    PRIVATE_REFERENCE_PROFILE_SCHEMA_VERSION,
    PrivateReferenceProfileConfigError,
    load_private_reference_profile_config,
)


def _document(secret_path: str) -> dict[str, object]:
    return {
        "schema_version": PRIVATE_REFERENCE_PROFILE_SCHEMA_VERSION,
        "profiles": {
            "grch38-primary": {
                "bindings": {
                    "bulk-rnaseq": {
                        "fasta": secret_path,
                        "fasta_sha256": "a" * 64,
                    }
                }
            }
        },
    }


@pytest.mark.parametrize("suffix", (".json", ".yaml"))
def test_private_reference_config_loads_json_or_yaml_without_repr_leak(
    tmp_path,
    suffix: str,
) -> None:
    secret = str(tmp_path / "operator" / "GRCh38.fa")
    path = tmp_path / f"references{suffix}"
    if suffix == ".json":
        path.write_text(json.dumps(_document(secret)), encoding="utf-8")
    else:
        path.write_text(
            "schema_version: helixweave-reference-profiles-v1\n"
            "profiles:\n"
            "  grch38-primary:\n"
            "    bindings:\n"
            "      bulk-rnaseq:\n"
            f"        fasta: {secret}\n"
            f"        fasta_sha256: {'a' * 64}\n",
            encoding="utf-8",
        )

    config = load_private_reference_profile_config(path)

    assert config.config_keys == ("grch38-primary",)
    assert config.workflow_ids_for("grch38-primary") == ("bulk-rnaseq",)
    assert config.binding_for("grch38-primary", "bulk-rnaseq")["fasta"] == secret
    assert secret not in repr(config)
    assert str(path) not in repr(config)


@pytest.mark.parametrize(
    "document",
    (
        '{"schema_version":"helixweave-reference-profiles-v1",'
        '"schema_version":"helixweave-reference-profiles-v1","profiles":{}}',
        '{"schema_version":"helixweave-reference-profiles-v1","profiles":'
        '{"a":{"bindings":{"bulk-rnaseq":{},"bulk-rnaseq":{}}}}}',
        '{"schema_version":"future","profiles":{}}',
    ),
)
def test_private_reference_config_rejects_duplicate_or_wrong_schema(
    tmp_path,
    document: str,
) -> None:
    path = tmp_path / "references.json"
    path.write_text(document, encoding="utf-8")

    with pytest.raises(PrivateReferenceProfileConfigError) as captured:
        load_private_reference_profile_config(path)

    assert captured.value.reason_code == "REFERENCE_PROFILE_CONFIG_INVALID"
    assert str(path) not in str(captured.value)


def test_private_reference_config_rejects_symlink_and_fifo_without_blocking(
    tmp_path,
) -> None:
    target = tmp_path / "target.json"
    target.write_text(json.dumps(_document("/secret/genome.fa")), encoding="utf-8")
    link = tmp_path / "link.json"
    link.symlink_to(target)
    fifo = tmp_path / "config.fifo"
    os.mkfifo(fifo)

    for path in (link, fifo):
        with pytest.raises(PrivateReferenceProfileConfigError):
            load_private_reference_profile_config(path)


def test_private_reference_config_unknown_key_is_redacted(tmp_path) -> None:
    path = tmp_path / "references.json"
    path.write_text(json.dumps(_document("/secret/genome.fa")), encoding="utf-8")
    config = load_private_reference_profile_config(path)

    with pytest.raises(PrivateReferenceProfileConfigError) as captured:
        config.binding_for("missing", "bulk-rnaseq")

    rendered = f"{captured.value!s} {captured.value!r}"
    assert "grch38-primary" not in rendered
    assert "/secret" not in rendered
