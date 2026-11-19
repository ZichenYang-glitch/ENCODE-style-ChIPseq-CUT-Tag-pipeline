"""Adoption integrity errors keep revision conflicts distinct from receipts."""

from sqlalchemy.exc import IntegrityError

from encode_pipeline.persistence.catalog_proposal_store import (
    AdoptionConflict,
    _raise_for_adoption_integrity,
)
from encode_pipeline.platform.experiment_catalog import CatalogError


def _integrity(message: str) -> IntegrityError:
    return IntegrityError("INSERT", {}, Exception(message))


def test_experiment_head_collision_stays_a_revision_conflict() -> None:
    exc = _integrity(
        "UNIQUE constraint failed: experiment_revisions.experiment_id, "
        "experiment_revisions.revision_number"
    )
    try:
        _raise_for_adoption_integrity(exc)
    except CatalogError as caught:
        assert caught.code == "catalog_revision_conflict"
    else:
        raise AssertionError("revision collision was not classified")


def test_duplicate_receipt_is_resolved_after_rollback() -> None:
    exc = _integrity("UNIQUE constraint failed: catalog_adoption_receipts.proposal_id")
    try:
        _raise_for_adoption_integrity(exc)
    except AdoptionConflict:
        return
    raise AssertionError("duplicate receipt was classified as a revision conflict")


def test_same_actor_key_is_resolved_as_the_committed_attempt() -> None:
    exc = _integrity(
        "UNIQUE constraint failed: catalog_adoption_receipts.actor_user_id, "
        "catalog_adoption_receipts.idempotency_key"
    )
    try:
        _raise_for_adoption_integrity(exc)
    except AdoptionConflict:
        return
    raise AssertionError("same-request collision was not left for replay")


def test_bound_parameters_do_not_choose_the_constraint() -> None:
    exc = IntegrityError(
        "INSERT INTO experiment_revisions (note) VALUES (?)",
        {"note": "catalog_adoption_receipts catalog_write_receipts UNIQUE"},
        Exception("UNIQUE constraint failed: catalog_adoption_receipts.proposal_id"),
    )
    try:
        _raise_for_adoption_integrity(exc)
    except AdoptionConflict:
        return
    raise AssertionError("receipt constraint was classified from bind text")


def test_revision_columns_ignore_receipt_names_in_the_statement() -> None:
    exc = IntegrityError(
        "INSERT INTO catalog_adoption_receipts (payload) VALUES (?)",
        ("catalog_write_receipts experiment_revision_heads",),
        Exception(
            "UNIQUE constraint failed: experiment_revisions.experiment_id, "
            "experiment_revisions.revision_number"
        ),
    )
    try:
        _raise_for_adoption_integrity(exc)
    except CatalogError as caught:
        assert caught.code == "catalog_revision_conflict"
    else:
        raise AssertionError("revision constraint followed the statement text")
