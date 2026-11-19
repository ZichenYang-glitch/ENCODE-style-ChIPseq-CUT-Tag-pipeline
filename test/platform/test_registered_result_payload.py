"""Payload branches for registered results that the happy-path HTTP test skips."""

from __future__ import annotations

from datetime import datetime, timezone

import pytest

from encode_pipeline.platform.registered_results import (
    ParentEdge,
    ParentView,
    RegisteredMetadata,
    RegisteredResultPayloadError,
    RevisionRecord,
    build_revision_digest,
    metadata_from_canonical_json,
    parse_metadata,
    public_revision,
    request_sha256,
    validate_idempotency_key,
    validate_project_id,
    validate_result_id,
    validate_result_revision_id,
    validate_stable_key,
)

PROJECT_ID = "prj_" + "a" * 32
RESULT_ID = "rreg_" + "b" * 32
REVISION_ID = "rrev_" + "c" * 32
OTHER_REVISION_ID = "rrev_" + "d" * 32
POOL_ID = "stgp_" + "e" * 32
DIGEST = "ab" * 32
PROFILE_ID = "refpr_" + "f" * 32


def _base(**overrides: object) -> dict:
    payload: dict = {
        "display_name": "Synthetic counts",
        "format": "tsv",
        "origin": {"kind": "external", "source": {"status": "unknown"}},
        "method": {"status": "unknown"},
        "reference": {"status": "unknown"},
    }
    payload.update(overrides)
    return payload


def _record(
    payload: str,
    *,
    parents: tuple[ParentView, ...] = (),
    created_at: datetime | None = None,
    head: str = REVISION_ID,
) -> RevisionRecord:
    return RevisionRecord(
        project_id=PROJECT_ID,
        project_archived_at=None,
        result_id=RESULT_ID,
        stable_key="synthetic",
        result_archived_at=None,
        result_revision_id=REVISION_ID,
        revision_number=1,
        canonical_payload=payload,
        revision_digest=DIGEST,
        created_at=created_at or datetime(2026, 11, 19, tzinfo=timezone.utc),
        storage_pool_id=POOL_ID,
        relative_path="external/fixture39.tsv",
        size_bytes=50,
        content_sha256=DIGEST,
        head_result_revision_id=head,
        parents=parents,
    )


def _raise(call) -> None:
    with pytest.raises(RegisteredResultPayloadError):
        call()


def test_identifiers_and_idempotency_keys_reject_malformed_values() -> None:
    assert validate_project_id(PROJECT_ID) == PROJECT_ID
    assert validate_result_id(RESULT_ID) == RESULT_ID
    assert validate_result_revision_id(REVISION_ID) == REVISION_ID
    assert validate_stable_key("synthetic-A") == "synthetic-A"
    assert validate_idempotency_key("irene-39") == "irene-39"
    for call in (
        lambda: validate_project_id("project"),
        lambda: validate_result_id(None),
        lambda: validate_result_revision_id("rrev_"),
        lambda: validate_stable_key(""),
        lambda: validate_stable_key("../counts"),
        lambda: validate_idempotency_key(""),
        lambda: validate_idempotency_key("a" * 129),
        lambda: validate_idempotency_key(" has space"),
        lambda: validate_idempotency_key(None),
    ):
        _raise(call)


def test_metadata_accepts_known_evidence_and_rejects_the_other_shapes() -> None:
    known = parse_metadata(
        _base(
            origin={
                "kind": "external",
                "source": {
                    "status": "known",
                    "namespace": "synthetic",
                    "accession": "SYNTH1",
                    "version": "1",
                },
            },
            method={"status": "known", "name": "counter", "version": "1"},
            reference={
                "status": "known",
                "kind": "external",
                "namespace": "synthetic",
                "accession": "REF1",
                "version": "1",
            },
        )
    )
    assert known.origin["source"]["accession"] == "SYNTH1"
    assert known.method["name"] == "counter"
    assert known.reference["kind"] == "external"
    profile = parse_metadata(
        _base(
            reference={
                "status": "known",
                "kind": "reference_profile_revision",
                "reference_profile_revision_id": PROFILE_ID,
                "public_identity_sha256": DIGEST,
            }
        )
    )
    assert profile.reference["reference_profile_revision_id"] == PROFILE_ID
    assert (
        parse_metadata(_base(reference={"status": "not_applicable"})).reference[
            "status"
        ]
        == "not_applicable"
    )
    downstream = parse_metadata(
        _base(
            origin={
                "kind": "downstream",
                "parents": [
                    {"result_revision_id": REVISION_ID, "role": "input"},
                    {"result_revision_id": OTHER_REVISION_ID, "role": "control"},
                ],
            }
        )
    )
    assert downstream.parent_refs == (
        (REVISION_ID, "input"),
        (OTHER_REVISION_ID, "control"),
    )
    assert parse_metadata(_base()).parent_refs == ()

    parents = [
        {"result_revision_id": f"rrev_{index:032x}", "role": "input"}
        for index in range(33)
    ]
    invalid = [
        lambda: parse_metadata([]),  # type: ignore[arg-type]
        lambda: parse_metadata({}),
        lambda: parse_metadata(_base(extra=True)),
        lambda: parse_metadata({**_base(), "schema_version": 2}),
        lambda: parse_metadata(_base(display_name="  ")),
        lambda: parse_metadata(_base(display_name="bad\nname")),
        lambda: parse_metadata(_base(display_name="x" * 256)),
        lambda: parse_metadata(_base(display_name=1)),
        lambda: parse_metadata(_base(format="has space")),
        lambda: parse_metadata(_base(origin=[])),
        lambda: parse_metadata(_base(origin={"kind": "external"})),
        lambda: parse_metadata(
            _base(origin={"kind": "external", "source": {}, "parents": []})
        ),
        lambda: parse_metadata(_base(origin={"kind": "internal", "source": {}})),
        lambda: parse_metadata(_base(origin={"kind": "external", "parents": []})),
        lambda: parse_metadata(
            _base(
                origin={
                    "kind": "downstream",
                    "source": {"status": "unknown"},
                }
            )
        ),
        lambda: parse_metadata(
            _base(origin={"kind": "external", "source": {"status": "other"}})
        ),
        lambda: parse_metadata(
            _base(reference={"status": "not_applicable", "kind": "x"})
        ),
        lambda: parse_metadata(
            _base(origin={"kind": "external", "source": {"status": "known"}})
        ),
        lambda: parse_metadata(_base(origin={"kind": "external", "source": "unknown"})),
        lambda: parse_metadata(
            _base(
                origin={
                    "kind": "external",
                    "source": {"status": "unknown", "namespace": "x"},
                }
            )
        ),
        lambda: parse_metadata(
            _base(
                origin={
                    "kind": "external",
                    "source": {
                        "status": "known",
                        "namespace": "a/b",
                        "accession": "1",
                        "version": "1",
                    },
                }
            )
        ),
        lambda: parse_metadata(_base(origin={"kind": "downstream", "parents": []})),
        lambda: parse_metadata(
            _base(origin={"kind": "downstream", "parents": parents})
        ),
        lambda: parse_metadata(
            _base(
                origin={
                    "kind": "downstream",
                    "parents": [{"result_revision_id": REVISION_ID}],
                }
            )
        ),
        lambda: parse_metadata(
            _base(
                origin={
                    "kind": "downstream",
                    "parents": [
                        {"result_revision_id": REVISION_ID, "role": "input"},
                        {"result_revision_id": REVISION_ID, "role": "input"},
                    ],
                }
            )
        ),
        lambda: parse_metadata(
            _base(
                origin={
                    "kind": "downstream",
                    "parents": [
                        {"result_revision_id": REVISION_ID, "role": "bad role"}
                    ],
                }
            )
        ),
        lambda: parse_metadata(_base(method=[])),
        lambda: parse_metadata(_base(method={"status": "known"})),
        lambda: parse_metadata(_base(method={"status": "other"})),
        lambda: parse_metadata(_base(method={"status": "unknown", "name": "x"})),
        lambda: parse_metadata(_base(reference=[])),
        lambda: parse_metadata(_base(reference={"status": "known"})),
        lambda: parse_metadata(_base(reference={"status": "other"})),
        lambda: parse_metadata(_base(reference={"status": "unknown", "kind": "x"})),
        lambda: parse_metadata(
            _base(
                reference={
                    "status": "known",
                    "kind": "reference_profile_revision",
                    "reference_profile_revision_id": "refpr_nope",
                    "public_identity_sha256": DIGEST,
                }
            )
        ),
        lambda: parse_metadata(
            _base(
                reference={
                    "status": "known",
                    "kind": "reference_profile_revision",
                    "reference_profile_revision_id": PROFILE_ID,
                    "public_identity_sha256": "zz",
                }
            )
        ),
        lambda: parse_metadata(
            _base(reference={"status": "known", "kind": "other", "namespace": "x"})
        ),
    ]
    for call in invalid:
        _raise(call)


def test_canonical_payload_round_trip_and_digest_frames() -> None:
    metadata = parse_metadata(_base())
    canonical = metadata.canonical_json()
    assert metadata_from_canonical_json(canonical) == metadata
    decoded = metadata_from_canonical_json(canonical)
    assert isinstance(decoded, RegisteredMetadata)
    for call in (
        lambda: metadata_from_canonical_json(""),
        lambda: metadata_from_canonical_json("{"),
        lambda: metadata_from_canonical_json("[]"),
        lambda: metadata_from_canonical_json(canonical.replace(":", ": ")),
    ):
        _raise(call)

    edge = ParentEdge(
        result_revision_id=REVISION_ID,
        role="input",
        revision_digest=DIGEST,
    )
    digest = build_revision_digest(
        project_id=PROJECT_ID,
        result_id=RESULT_ID,
        result_revision_id=REVISION_ID,
        revision_number=1,
        canonical_payload=canonical,
        storage_pool_id=POOL_ID,
        relative_path="external/fixture39.tsv",
        size_bytes=50,
        content_sha256=DIGEST,
        parents=(edge,),
    )
    assert len(digest) == 64
    hashed = request_sha256("register", PROJECT_ID, None, {"display_name": "A"})
    assert hashed == request_sha256("register", PROJECT_ID, None, {"display_name": "A"})
    assert hashed != request_sha256(
        "register",
        PROJECT_ID,
        RESULT_ID,
        {"display_name": "A"},
    )
    for call in (
        lambda: build_revision_digest(
            project_id=PROJECT_ID,
            result_id=RESULT_ID,
            result_revision_id=REVISION_ID,
            revision_number=True,
            canonical_payload=canonical,
            storage_pool_id=POOL_ID,
            relative_path="external/fixture39.tsv",
            size_bytes=50,
            content_sha256=DIGEST,
            parents=(),
        ),
        lambda: build_revision_digest(
            project_id=PROJECT_ID,
            result_id=RESULT_ID,
            result_revision_id=REVISION_ID,
            revision_number=0,
            canonical_payload=canonical,
            storage_pool_id=POOL_ID,
            relative_path="external/fixture39.tsv",
            size_bytes=50,
            content_sha256=DIGEST,
            parents=(),
        ),
        lambda: build_revision_digest(
            project_id=PROJECT_ID,
            result_id=RESULT_ID,
            result_revision_id=REVISION_ID,
            revision_number=1,
            canonical_payload=canonical,
            storage_pool_id=POOL_ID,
            relative_path="external/fixture39.tsv",
            size_bytes=True,
            content_sha256=DIGEST,
            parents=(),
        ),
        lambda: build_revision_digest(
            project_id=PROJECT_ID,
            result_id=RESULT_ID,
            result_revision_id=REVISION_ID,
            revision_number=1,
            canonical_payload=canonical,
            storage_pool_id=POOL_ID,
            relative_path="external/fixture39.tsv",
            size_bytes=-1,
            content_sha256=DIGEST,
            parents=(),
        ),
        lambda: ParentEdge(
            result_revision_id="nope",
            role="input",
            revision_digest=DIGEST,
        ),
    ):
        _raise(call)


def test_public_revision_hides_storage_and_renders_downstream_parents() -> None:
    metadata = parse_metadata(
        _base(
            origin={
                "kind": "downstream",
                "parents": [{"result_revision_id": OTHER_REVISION_ID, "role": "input"}],
            }
        )
    )
    parent = ParentView(
        result_id=RESULT_ID,
        result_revision_id=OTHER_REVISION_ID,
        revision_number=1,
        revision_status="current",
        current_result_revision_id=OTHER_REVISION_ID,
        display_name="Parent",
        role="input",
        revision_digest=DIGEST,
    )
    published = public_revision(
        _record(
            metadata.canonical_json(),
            parents=(parent,),
            created_at=datetime(2026, 11, 19, 1, 2, 3),
            head=OTHER_REVISION_ID,
        )
    )
    assert published["origin"]["kind"] == "downstream"
    assert published["origin"]["parents"][0]["role"] == "input"
    assert published["revision_status"] == "superseded"
    assert published["created_at"].endswith("Z")
    assert "storage_pool_id" not in published
    external = public_revision(_record(parse_metadata(_base()).canonical_json()))
    assert external["origin"]["kind"] == "external"
    assert external["revision_status"] == "current"
