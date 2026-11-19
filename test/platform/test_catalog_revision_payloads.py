"""Catalog validation and stored experiment revision compatibility."""

from base64 import urlsafe_b64decode, urlsafe_b64encode
from datetime import datetime, timezone
import json

import pytest

from encode_pipeline.platform.experiment_catalog import (
    CATALOG_CURSOR_PREFIX,
    CatalogCursorPosition,
    CatalogError,
    TermPin,
    build_experiment_revision_payload_digest,
    decode_catalog_cursor,
    encode_catalog_cursor,
    interpret_sample_catalog,
    interpret_stored_experiment_payload,
    parse_page_filters,
)


SAMPLE_REVISION = "smpr_" + "a" * 32
SPECIES = {"term_id": "term_" + "b" * 32, "term_revision_id": "trmv_" + "c" * 32}
CELL = {"term_id": "term_" + "d" * 32, "term_revision_id": "trmv_" + "e" * 32}
UNKNOWN = {"status": "unknown"}


def _stored_experiment():
    return {
        "schema_version": 1,
        "display_name": "Pilot",
        "assay": UNKNOWN,
        "target": UNKNOWN,
        "batch": UNKNOWN,
        "samples": [
            {
                "ordinal": 0,
                "sample_revision_id": SAMPLE_REVISION,
                "payload_digest": "f" * 64,
                "replicate_role": UNKNOWN,
            }
        ],
    }


def _stored_material():
    return {
        "schema_version": 1,
        "species": {"status": "known", **SPECIES},
        "tissue": UNKNOWN,
        **{
            field: {"status": "unknown", "values": []}
            for field in ("cell_type", "cell_line", "treatment", "tag")
        },
    }


def _lookup(revision_id):
    terms = {
        SPECIES["term_revision_id"]: TermPin(
            **SPECIES, term_type="species", archived=True, display_name="Old species"
        ),
        CELL["term_revision_id"]: TermPin(
            **CELL, term_type="cell_type", archived=False, display_name="T cell"
        ),
    }
    return terms.get(revision_id)


def test_old_experiment_without_group_and_empty_group_have_same_identity():
    older = _stored_experiment()
    old_bytes = json.dumps(older, sort_keys=True, separators=(",", ":"))
    older["samples"][0]["replicate_group"] = ""
    empty_group_bytes = json.dumps(older)
    older["samples"][0]["replicate_group"] = "group-a"
    grouped_bytes = json.dumps(older)

    old_payload = interpret_stored_experiment_payload(old_bytes)
    assert old_payload.samples[0].sample_revision_id == SAMPLE_REVISION
    assert old_payload.samples[0].payload_digest == "f" * 64
    assert old_payload.samples[0].replicate_group is None
    normalized = interpret_stored_experiment_payload(empty_group_bytes)
    assert normalized == old_payload
    grouped = interpret_stored_experiment_payload(grouped_bytes)
    assert grouped.samples[0].replicate_group == "group-a"
    original_digest = build_experiment_revision_payload_digest(old_bytes)
    assert (
        build_experiment_revision_payload_digest(normalized.canonical_json())
        == original_digest
    )
    assert (
        build_experiment_revision_payload_digest(grouped.canonical_json())
        != original_digest
    )


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("ordinal", 1),
        ("sample_revision_id", "latest"),
        ("payload_digest", "not-a-digest"),
        ("replicate_group", "two words"),
        ("replicate_role", {"status": "unknown", **CELL}),
    ],
)
def test_corrupt_stored_sample_pins_cannot_supply_an_experiment_digest(field, value):
    payload = _stored_experiment()
    payload["samples"][0][field] = value
    with pytest.raises(CatalogError) as failure:
        build_experiment_revision_payload_digest(json.dumps(payload))
    assert failure.value.code == "catalog_payload_invalid"


def test_repeated_revision_cannot_supply_an_experiment_digest():
    payload = _stored_experiment()
    payload["samples"].append(payload["samples"][0] | {"ordinal": 1})
    with pytest.raises(CatalogError) as failure:
        build_experiment_revision_payload_digest(json.dumps(payload))
    assert failure.value.code == "catalog_payload_invalid"


@pytest.mark.parametrize(
    ("patch", "code"),
    [
        ({"schema_version": True}, "catalog_payload_invalid"),
        ({"schema_version": 2}, "catalog_payload_invalid"),
        ({"species": "unknown"}, "catalog_payload_invalid"),
        ({"species": {"status": []}}, "catalog_payload_invalid"),
        (
            {"species": {"status": "known", "term_id": SPECIES["term_id"]}},
            "catalog_payload_invalid",
        ),
        ({"species": {"status": "ambiguous", "values": []}}, "catalog_payload_invalid"),
        ({"species": {"status": "known", **CELL}}, "catalog_term_type_mismatch"),
        (
            {
                "species": {
                    "status": "known",
                    **(SPECIES | {"term_revision_id": "trmv_" + "1" * 32}),
                }
            },
            "catalog_payload_invalid",
        ),
        ({"cell_type": None}, "catalog_payload_invalid"),
        ({"cell_type": {"status": {}}}, "catalog_payload_invalid"),
        (
            {"cell_type": {"status": "known", "values": [CELL], "candidates": []}},
            "catalog_payload_invalid",
        ),
        ({"cell_type": {"status": "known", "values": []}}, "catalog_payload_invalid"),
        (
            {"cell_type": {"status": "known", "values": [CELL, CELL]}},
            "catalog_payload_invalid",
        ),
        ({"cell_type": {"status": "known", "values": CELL}}, "catalog_payload_invalid"),
        (
            {
                "cell_type": {
                    "status": "known",
                    "values": [CELL | {"display_name": "T cell"}],
                }
            },
            "catalog_payload_invalid",
        ),
        (
            {"cell_type": {"status": "unknown", "values": [CELL]}},
            "catalog_payload_invalid",
        ),
        (
            {"cell_type": {"status": "unknown", "values": [], "candidates": []}},
            "catalog_payload_invalid",
        ),
        (
            {"cell_type": {"status": "ambiguous", "values": [], "extra": True}},
            "catalog_payload_invalid",
        ),
    ],
)
def test_malformed_or_mistyped_sample_classification_is_rejected(patch, code):
    with pytest.raises(CatalogError) as failure:
        interpret_sample_catalog(
            _stored_material() | {"species": UNKNOWN} | patch, _lookup, for_write=True
        )
    assert failure.value.code == code


def test_archived_term_stays_readable_on_an_existing_sample_revision():
    material = _stored_material()
    stored = interpret_sample_catalog(material, _lookup, for_write=False)
    assert stored.species.term.term_revision_id == SPECIES["term_revision_id"]
    with pytest.raises(CatalogError) as failure:
        interpret_sample_catalog(material, _lookup, for_write=True)
    assert failure.value.code == "catalog_term_archived"


@pytest.mark.parametrize(
    "payload",
    [
        "{",
        "[]",
        "{}",
        json.dumps(_stored_experiment() | {"samples": []}),
        json.dumps(_stored_experiment() | {"samples": ["sample-id"]}),
        json.dumps(_stored_experiment() | {"samples": [{"ordinal": 0}]}),
        json.dumps(_stored_experiment() | {"schema_version": 2}),
    ],
)
def test_invalid_stored_revision_is_unreadable(payload):
    with pytest.raises(CatalogError) as failure:
        interpret_stored_experiment_payload(payload)
    assert failure.value.code == "catalog_payload_invalid"


@pytest.mark.parametrize(
    "change",
    [
        {"v": 2},
        {"kind": "experiments"},
        {"unexpected": True},
        {"position": []},
        {"position": {"created_at": "yesterday", "object_id": "sample"}},
        {"position": {"created_at": "2026-11-18T00:00:00", "object_id": "sample"}},
        {"position": {"created_at": "2026-11-18T00:00:00Z", "object_id": ""}},
        {},
    ],
)
def test_noncanonical_catalog_cursor_is_rejected(change):
    filters = parse_page_filters(
        kind="samples",
        project_id="prj_" + "1" * 32,
        revision_view=None,
        include_archived=False,
        term_type=None,
        created_from=None,
        created_before=None,
        catalog_state=None,
        field_terms={},
        field_statuses={},
    )
    position = CatalogCursorPosition(
        created_at=datetime(2026, 11, 18, tzinfo=timezone.utc), object_id="sample"
    )
    cursor = encode_catalog_cursor(kind="samples", filters=filters, position=position)
    assert decode_catalog_cursor(cursor, kind="samples", filters=filters) == position

    token = cursor.removeprefix(CATALOG_CURSOR_PREFIX)
    body = json.loads(urlsafe_b64decode(token + "=" * (-len(token) % 4)))
    # Default JSON spacing creates an equivalent but noncanonical token for {}.
    changed_body = urlsafe_b64encode(json.dumps(body | change).encode()).decode()
    tampered = CATALOG_CURSOR_PREFIX + changed_body.rstrip("=")
    with pytest.raises(CatalogError) as failure:
        decode_catalog_cursor(tampered, kind="samples", filters=filters)
    assert (failure.value.status_code, failure.value.code) == (
        400,
        "catalog_cursor_invalid",
    )
