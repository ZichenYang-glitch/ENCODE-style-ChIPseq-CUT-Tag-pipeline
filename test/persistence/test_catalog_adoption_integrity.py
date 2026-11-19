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
