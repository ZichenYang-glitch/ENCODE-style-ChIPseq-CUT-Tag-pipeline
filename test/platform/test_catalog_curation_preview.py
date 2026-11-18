"""Acceptance and contract tests for M2 offline catalog curation preview."""

from __future__ import annotations

from datetime import timedelta
import json
from pathlib import Path

import pytest

from encode_pipeline.platform.catalog_curation_preview import (
    MAX_CHANGESET_ITEMS,
    PROPOSAL_ADOPTION_WINDOW,
    EXPIRED_PROPOSAL_RETENTION,
    CatalogPreviewError,
    CatalogProposal,
    ProposalCitation,
    ProposalItem,
    SourceReference,
    SYNTHETIC_BASE_TIME,
    SYNTHETIC_EXP_1_ID,
    SYNTHETIC_EXP_2_ID,
    SYNTHETIC_PROJECT_ID,
    SYNTHETIC_PROPOSAL_ID,
    SYNTHETIC_SMP_1_ID,
    build_synthetic_m2_preview_fixture,
    export_m2_preview_contract_examples,
)


REPO_ROOT = Path(__file__).resolve().parents[2]
SPEC_PATH = REPO_ROOT / "docs/architecture/experiment-catalog-m2-preview-contract.md"
EXAMPLES_PATH = (
    REPO_ROOT / "docs/architecture/experiment-catalog-m2-preview-examples.json"
)


def test_m2_preview_spec_and_examples_stay_in_sync() -> None:
    """Published contract examples must match the live synthetic fixture output."""

    spec_text = SPEC_PATH.read_text(encoding="utf-8")
    for case_id in ("M2-P1", "M2-P2", "M2-P3", "M2-P4", "M2-P5", "M2-P6"):
        assert case_id in spec_text

    committed_examples = json.loads(EXAMPLES_PATH.read_text(encoding="utf-8"))
    assert committed_examples == export_m2_preview_contract_examples()


def test_m2_p1_proposal_citations_and_role_separation() -> None:
    """M2-P1: Citations stay verifiable; member sees redacted view and cannot adopt."""

    store = build_synthetic_m2_preview_fixture()
    now = SYNTHETIC_BASE_TIME + timedelta(minutes=1)

    member_view = store.get_proposal_view(
        actor_user_id="usr_member_01",
        proposal_id=SYNTHETIC_PROPOSAL_ID,
        now=now,
    )
    assert member_view["adoptable_by_caller"] is False
    item_3_member = member_view["items"][2]
    assert item_3_member["after"]["cell_type"] == {"status": "ambiguous", "values": []}
    assert "candidates" not in item_3_member["after"]["cell_type"]
    assert "aliases" not in json.dumps(member_view)

    admin_view = store.get_proposal_view(
        actor_user_id="usr_admin_a",
        proposal_id=SYNTHETIC_PROPOSAL_ID,
        now=now,
    )
    assert admin_view["adoptable_by_caller"] is True
    item_3_admin = admin_view["items"][2]
    candidates = item_3_admin["after"]["cell_type"]["candidates"]
    assert [c["display_name"] for c in candidates] == [
        "regulatory T cell",
        "keratinocyte",
    ]
    assert item_3_admin["citations"]["rule_version"] == "m2-catalog-rules/1"
    assert item_3_admin["citations"]["reason_codes"] == [
        "ALIAS_AMBIGUOUS_CANDIDATES",
        "TAG_TERM_MATCH",
    ]
    assert len(item_3_admin["citations"]["source_refs"]) == 4

    with pytest.raises(CatalogPreviewError) as preview_exc:
        store.preview_proposal_subset(
            actor_user_id="usr_member_01",
            proposal_id=SYNTHETIC_PROPOSAL_ID,
            selected_item_ids=["item_1"],
            now=now,
        )
    assert preview_exc.value.status_code == 403
    assert preview_exc.value.code == "ADMINISTRATOR_REQUIRED"

    admin_preview = store.preview_proposal_subset(
        actor_user_id="usr_admin_a",
        proposal_id=SYNTHETIC_PROPOSAL_ID,
        selected_item_ids=["item_1"],
        now=now,
    )
    with pytest.raises(CatalogPreviewError) as adopt_exc:
        store.adopt_proposal_subset(
            actor_user_id="usr_member_01",
            proposal_id=SYNTHETIC_PROPOSAL_ID,
            preview_digest=admin_preview["preview_digest"],
            selected_item_ids=["item_1"],
            expected_revisions=admin_preview["expected_revisions"],
            idempotency_key="idem-member-rejected",
            now=now,
        )
    assert adopt_exc.value.status_code == 403
    assert adopt_exc.value.code == "ADMINISTRATOR_REQUIRED"
    assert store.get_experiment_head(SYNTHETIC_EXP_1_ID)["revision_number"] == 1


def test_m2_p2_subset_selection_and_immutable_preview_digest() -> None:
    """M2-P2: Admin selects a strict subset; excluded items are never adopted."""

    store = build_synthetic_m2_preview_fixture()
    t_preview = SYNTHETIC_BASE_TIME + timedelta(minutes=2)
    t_adopt = SYNTHETIC_BASE_TIME + timedelta(minutes=3)

    preview = store.preview_proposal_subset(
        actor_user_id="usr_admin_a",
        proposal_id=SYNTHETIC_PROPOSAL_ID,
        selected_item_ids=["item_1", "item_3"],
        now=t_preview,
    )
    assert preview["selected_item_ids"] == ["item_1", "item_3"]
    assert preview["excluded_item_ids"] == ["item_2"]
    assert preview["has_conflicts"] is False
    assert preview["expected_revisions"] == {
        SYNTHETIC_EXP_1_ID: 1,
        SYNTHETIC_SMP_1_ID: 1,
    }

    # Tampering with selection while reusing preview_digest fails closed.
    with pytest.raises(CatalogPreviewError) as tamper_exc:
        store.adopt_proposal_subset(
            actor_user_id="usr_admin_a",
            proposal_id=SYNTHETIC_PROPOSAL_ID,
            preview_digest=preview["preview_digest"],
            selected_item_ids=["item_1", "item_2", "item_3"],
            expected_revisions={
                SYNTHETIC_EXP_1_ID: 1,
                SYNTHETIC_EXP_2_ID: 1,
                SYNTHETIC_SMP_1_ID: 1,
            },
            idempotency_key="idem-tampered-subset",
            now=t_adopt,
        )
    assert tamper_exc.value.status_code == 409
    assert tamper_exc.value.code == "catalog_preview_digest_mismatch"

    receipt = store.adopt_proposal_subset(
        actor_user_id="usr_admin_a",
        proposal_id=SYNTHETIC_PROPOSAL_ID,
        preview_digest=preview["preview_digest"],
        selected_item_ids=preview["selected_item_ids"],
        expected_revisions=preview["expected_revisions"],
        idempotency_key="idem-subset-ok",
        now=t_adopt,
    )
    assert receipt["replayed"] is False
    assert receipt["selected_item_ids"] == ["item_1", "item_3"]
    assert receipt["excluded_item_ids"] == ["item_2"]
    assert [t["target_id"] for t in receipt["applied_targets"]] == [
        SYNTHETIC_EXP_1_ID,
        SYNTHETIC_SMP_1_ID,
    ]

    exp1_head = store.get_experiment_head(SYNTHETIC_EXP_1_ID)
    exp2_head = store.get_experiment_head(SYNTHETIC_EXP_2_ID)
    smp1_head = store.get_sample_head(SYNTHETIC_SMP_1_ID)
    assert exp1_head["revision_number"] == 2
    assert exp1_head["payload"]["batch"]["term_id"] == (
        "term_34343434343434343434343434343434"
    )
    assert exp2_head["revision_number"] == 1
    assert exp2_head["payload"]["target"] == {"status": "unknown"}
    assert smp1_head["revision_number"] == 2
    assert smp1_head["catalog"]["cell_type"]["status"] == "ambiguous"
    assert smp1_head["catalog"]["tag"]["status"] == "known"
    # Non-catalog attributes outside catalog block remain untouched.
    assert smp1_head["extra_attributes"] == {
        "legacy_note": "preserve-me-outside-catalog"
    }


def test_m2_p3_concurrent_manual_revision_causes_atomic_conflict() -> None:
    """M2-P3: Competing manual revision on one target aborts the whole batch."""

    store = build_synthetic_m2_preview_fixture()
    t_preview = SYNTHETIC_BASE_TIME + timedelta(minutes=2)
    t_manual = SYNTHETIC_BASE_TIME + timedelta(minutes=3)
    t_adopt = SYNTHETIC_BASE_TIME + timedelta(minutes=4)

    preview = store.preview_proposal_subset(
        actor_user_id="usr_admin_a",
        proposal_id=SYNTHETIC_PROPOSAL_ID,
        selected_item_ids=["item_1", "item_3"],
        now=t_preview,
    )
    manual_result = store.manual_revise_experiment(
        actor_user_id="usr_admin_b",
        experiment_id=SYNTHETIC_EXP_1_ID,
        expected_revision_number=1,
        idempotency_key="idem-admin-b-manual",
        updates={"display_name": "Treg H3K27me3 CUT&Tag (Manual Rev 2)"},
        now=t_manual,
    )
    assert manual_result["revision_number"] == 2

    # Previewing again surfaces conflict_detected before adoption too.
    refreshed_preview = store.preview_proposal_subset(
        actor_user_id="usr_admin_a",
        proposal_id=SYNTHETIC_PROPOSAL_ID,
        selected_item_ids=["item_1", "item_3"],
        now=t_adopt,
    )
    assert refreshed_preview["has_conflicts"] is True
    assert refreshed_preview["conflicts"][0]["target_id"] == SYNTHETIC_EXP_1_ID
    assert refreshed_preview["conflicts"][0]["current_revision_number"] == 2

    with pytest.raises(CatalogPreviewError) as conflict_exc:
        store.adopt_proposal_subset(
            actor_user_id="usr_admin_a",
            proposal_id=SYNTHETIC_PROPOSAL_ID,
            preview_digest=preview["preview_digest"],
            selected_item_ids=preview["selected_item_ids"],
            expected_revisions=preview["expected_revisions"],
            idempotency_key="idem-admin-a-stale",
            now=t_adopt,
        )
    err = conflict_exc.value
    assert err.status_code == 409
    assert err.code == "catalog_revision_conflict"
    assert len(err.conflicts) == 1
    assert err.conflicts[0].target_id == SYNTHETIC_EXP_1_ID
    assert err.conflicts[0].expected_revision_number == 1
    assert err.conflicts[0].current_revision_number == 2

    # Zero partial commit: exp_1 keeps Admin B's rev 2, smp_1 stays at rev 1.
    exp1_head = store.get_experiment_head(SYNTHETIC_EXP_1_ID)
    assert exp1_head["revision_number"] == 2
    assert (
        exp1_head["payload"]["display_name"] == "Treg H3K27me3 CUT&Tag (Manual Rev 2)"
    )
    assert exp1_head["payload"]["batch"] == {"status": "unknown"}
    assert store.get_sample_head(SYNTHETIC_SMP_1_ID)["revision_number"] == 1


def test_m2_p4_two_admins_concurrent_expected_revisions_race() -> None:
    """M2-P4: Two admins concurrently adopting expected_revision=1 produce one winner and one 409."""

    store = build_synthetic_m2_preview_fixture()
    t_preview = SYNTHETIC_BASE_TIME + timedelta(minutes=2)
    t_commit_a = SYNTHETIC_BASE_TIME + timedelta(minutes=3)
    t_commit_b = SYNTHETIC_BASE_TIME + timedelta(minutes=4)

    preview_a = store.preview_proposal_subset(
        actor_user_id="usr_admin_a",
        proposal_id=SYNTHETIC_PROPOSAL_ID,
        selected_item_ids=["item_1", "item_3"],
        now=t_preview,
    )
    preview_b = store.preview_proposal_subset(
        actor_user_id="usr_admin_b",
        proposal_id=SYNTHETIC_PROPOSAL_ID,
        selected_item_ids=["item_1", "item_2"],
        now=t_preview,
    )

    winner_receipt = store.adopt_proposal_subset(
        actor_user_id="usr_admin_a",
        proposal_id=SYNTHETIC_PROPOSAL_ID,
        preview_digest=preview_a["preview_digest"],
        selected_item_ids=preview_a["selected_item_ids"],
        expected_revisions=preview_a["expected_revisions"],
        idempotency_key="idem-race-admin-a",
        now=t_commit_a,
    )
    assert winner_receipt["replayed"] is False

    with pytest.raises(CatalogPreviewError) as loser_exc:
        store.adopt_proposal_subset(
            actor_user_id="usr_admin_b",
            proposal_id=SYNTHETIC_PROPOSAL_ID,
            preview_digest=preview_b["preview_digest"],
            selected_item_ids=preview_b["selected_item_ids"],
            expected_revisions=preview_b["expected_revisions"],
            idempotency_key="idem-race-admin-b",
            now=t_commit_b,
        )
    assert loser_exc.value.status_code == 409
    assert loser_exc.value.code == "catalog_revision_conflict"
    assert loser_exc.value.conflicts[0].target_id == SYNTHETIC_EXP_1_ID
    assert loser_exc.value.conflicts[0].current_revision_number == 2

    # item_2 from losing subset was NOT committed either.
    assert store.get_experiment_head(SYNTHETIC_EXP_1_ID)["revision_number"] == 2
    assert store.get_experiment_head(SYNTHETIC_EXP_2_ID)["revision_number"] == 1
    assert store.get_sample_head(SYNTHETIC_SMP_1_ID)["revision_number"] == 2


def test_m2_p5_dropped_response_idempotent_replay_and_receipt_discovery() -> None:
    """M2-P5: Same Idempotency-Key + same body replays receipt; different body conflicts."""

    store = build_synthetic_m2_preview_fixture()
    t_preview = SYNTHETIC_BASE_TIME + timedelta(minutes=2)
    t_adopt = SYNTHETIC_BASE_TIME + timedelta(minutes=3)
    t_retry = SYNTHETIC_BASE_TIME + timedelta(minutes=4)

    preview = store.preview_proposal_subset(
        actor_user_id="usr_admin_a",
        proposal_id=SYNTHETIC_PROPOSAL_ID,
        selected_item_ids=["item_1", "item_3"],
        now=t_preview,
    )
    first_receipt = store.adopt_proposal_subset(
        actor_user_id="usr_admin_a",
        proposal_id=SYNTHETIC_PROPOSAL_ID,
        preview_digest=preview["preview_digest"],
        selected_item_ids=preview["selected_item_ids"],
        expected_revisions=preview["expected_revisions"],
        idempotency_key="idem-dropped-response",
        now=t_adopt,
    )
    assert first_receipt["replayed"] is False

    # Dropped response retry: even though heads are now rev 2, replay succeeds cleanly.
    replayed_receipt = store.adopt_proposal_subset(
        actor_user_id="usr_admin_a",
        proposal_id=SYNTHETIC_PROPOSAL_ID,
        preview_digest=preview["preview_digest"],
        selected_item_ids=preview["selected_item_ids"],
        expected_revisions=preview["expected_revisions"],
        idempotency_key="idem-dropped-response",
        now=t_retry,
    )
    assert replayed_receipt["replayed"] is True
    assert replayed_receipt["receipt_id"] == first_receipt["receipt_id"]
    assert replayed_receipt["applied_targets"] == first_receipt["applied_targets"]
    assert store.get_experiment_head(SYNTHETIC_EXP_1_ID)["revision_number"] == 2
    assert store.get_sample_head(SYNTHETIC_SMP_1_ID)["revision_number"] == 2

    # Receipt discovery after re-login finds the original receipt.
    discovered = store.list_adoption_receipts(
        actor_user_id="usr_admin_a",
        project_id=SYNTHETIC_PROJECT_ID,
        proposal_id=SYNTHETIC_PROPOSAL_ID,
    )
    assert len(discovered["receipts"]) == 1
    assert discovered["receipts"][0]["receipt_id"] == first_receipt["receipt_id"]

    # Same key with different body returns 409 catalog_idempotency_conflict.
    with pytest.raises(CatalogPreviewError) as conflict_exc:
        store.adopt_proposal_subset(
            actor_user_id="usr_admin_a",
            proposal_id=SYNTHETIC_PROPOSAL_ID,
            preview_digest=preview["preview_digest"],
            selected_item_ids=["item_1"],
            expected_revisions={SYNTHETIC_EXP_1_ID: 1},
            idempotency_key="idem-dropped-response",
            now=t_retry,
        )
    assert conflict_exc.value.status_code == 409
    assert conflict_exc.value.code == "catalog_idempotency_conflict"

    # Missing Idempotency-Key fails with 400 API_REQUEST_INVALID.
    with pytest.raises(CatalogPreviewError) as missing_key_exc:
        store.adopt_proposal_subset(
            actor_user_id="usr_admin_a",
            proposal_id=SYNTHETIC_PROPOSAL_ID,
            preview_digest=preview["preview_digest"],
            selected_item_ids=preview["selected_item_ids"],
            expected_revisions=preview["expected_revisions"],
            idempotency_key=None,
            now=t_retry,
        )
    assert missing_key_exc.value.status_code == 400
    assert missing_key_exc.value.code == "API_REQUEST_INVALID"


def test_m2_p6_expiry_revocation_archived_terms_and_egress_guards() -> None:
    """M2-P6: Expiry, admin revocation, archived terms, and invalid citations fail closed."""

    # 1. Expiry: adopt fails with 410 after 30 minutes; read-only view lasts 7 days.
    expired_store = build_synthetic_m2_preview_fixture()
    preview = expired_store.preview_proposal_subset(
        actor_user_id="usr_admin_a",
        proposal_id=SYNTHETIC_PROPOSAL_ID,
        selected_item_ids=["item_1"],
        now=SYNTHETIC_BASE_TIME + timedelta(minutes=5),
    )
    t_after_30m = SYNTHETIC_BASE_TIME + PROPOSAL_ADOPTION_WINDOW
    view_expired = expired_store.get_proposal_view(
        actor_user_id="usr_admin_a",
        proposal_id=SYNTHETIC_PROPOSAL_ID,
        now=t_after_30m,
    )
    assert view_expired["proposal_state"] == "expired"
    assert view_expired["adoptable_by_caller"] is False

    with pytest.raises(CatalogPreviewError) as exp_exc:
        expired_store.adopt_proposal_subset(
            actor_user_id="usr_admin_a",
            proposal_id=SYNTHETIC_PROPOSAL_ID,
            preview_digest=preview["preview_digest"],
            selected_item_ids=preview["selected_item_ids"],
            expected_revisions=preview["expected_revisions"],
            idempotency_key="idem-expired",
            now=t_after_30m,
        )
    assert exp_exc.value.status_code == 410
    assert exp_exc.value.code == "catalog_proposal_expired"

    t_after_7d = SYNTHETIC_BASE_TIME + EXPIRED_PROPOSAL_RETENTION
    with pytest.raises(CatalogPreviewError) as purged_exc:
        expired_store.get_proposal_view(
            actor_user_id="usr_admin_a",
            proposal_id=SYNTHETIC_PROPOSAL_ID,
            now=t_after_7d,
        )
    assert purged_exc.value.status_code == 404
    assert purged_exc.value.code == "catalog_not_found"

    # 2. Admin downgraded to member before commit -> 403 and zero writes.
    revoked_store = build_synthetic_m2_preview_fixture()
    revoked_preview = revoked_store.preview_proposal_subset(
        actor_user_id="usr_admin_a",
        proposal_id=SYNTHETIC_PROPOSAL_ID,
        selected_item_ids=["item_1"],
        now=SYNTHETIC_BASE_TIME + timedelta(minutes=1),
    )
    revoked_store.set_actor_role("usr_admin_a", role="member", active=True)
    with pytest.raises(CatalogPreviewError) as revoked_exc:
        revoked_store.adopt_proposal_subset(
            actor_user_id="usr_admin_a",
            proposal_id=SYNTHETIC_PROPOSAL_ID,
            preview_digest=revoked_preview["preview_digest"],
            selected_item_ids=revoked_preview["selected_item_ids"],
            expected_revisions=revoked_preview["expected_revisions"],
            idempotency_key="idem-revoked-admin",
            now=SYNTHETIC_BASE_TIME + timedelta(minutes=2),
        )
    assert revoked_exc.value.status_code == 403
    assert revoked_exc.value.code == "ADMINISTRATOR_REQUIRED"
    assert revoked_store.get_experiment_head(SYNTHETIC_EXP_1_ID)["revision_number"] == 1

    # 3. Archived term referenced by proposal fails with 409 catalog_term_archived.
    archived_term_store = build_synthetic_m2_preview_fixture()
    archived_term_store.archive_term_revision("trmv_56565656565656565656565656565656")
    with pytest.raises(CatalogPreviewError) as term_exc:
        archived_term_store.preview_proposal_subset(
            actor_user_id="usr_admin_a",
            proposal_id=SYNTHETIC_PROPOSAL_ID,
            selected_item_ids=["item_1"],
            now=SYNTHETIC_BASE_TIME + timedelta(minutes=1),
        )
    assert term_exc.value.status_code == 409
    assert term_exc.value.code == "catalog_term_archived"

    # 4. Host paths and >100-item changesets are rejected by proposal construction.
    with pytest.raises(CatalogPreviewError) as path_exc:
        SourceReference(
            kind="catalog_projection_row",
            ref_id="/tmp/secret/host/path.tsv",
        )
    assert path_exc.value.status_code == 422
    assert path_exc.value.code == "catalog_payload_invalid"

    too_many_items = [
        ProposalItem(
            item_id=f"item_{idx}",
            target_kind="experiment",
            target_id=f"exp_{idx:032x}",
            base_revision_id=f"expr_{idx:032x}",
            base_revision_number=1,
            allowed_fields=("batch",),
            before={"batch": {"status": "unknown"}},
            after={"batch": {"status": "not_applicable"}},
            citation=ProposalCitation(
                rule_version="m2-catalog-rules/1",
                vocabulary_version="terms-snapshot/2026-11-18",
                reason_codes=("MISSING_EVIDENCE_UNKNOWN",),
                source_refs=(
                    SourceReference(
                        kind="experiment_revision",
                        ref_id=f"expr_{idx:032x}",
                        evidence_token="bulk-bound-check",
                    ),
                ),
            ),
        )
        for idx in range(1, MAX_CHANGESET_ITEMS + 2)
    ]
    with pytest.raises(CatalogPreviewError) as limit_exc:
        CatalogProposal.create(
            proposal_id=SYNTHETIC_PROPOSAL_ID,
            project_id=SYNTHETIC_PROJECT_ID,
            created_by_user_id="usr_member_01",
            source="offline_suggestion",
            items=too_many_items,
            created_at=SYNTHETIC_BASE_TIME,
        )
    assert limit_exc.value.status_code == 422
    assert limit_exc.value.code == "catalog_payload_invalid"


def test_m2_preview_validation_and_error_branches() -> None:
    """Cover term type mismatch, unknown revisions, archived targets, and selection guards."""

    store = build_synthetic_m2_preview_fixture()
    now = SYNTHETIC_BASE_TIME + timedelta(minutes=1)

    # Unknown actor / inactive actor -> 401 AUTHENTICATION_REQUIRED.
    with pytest.raises(CatalogPreviewError) as auth_exc:
        store.get_proposal_view(
            actor_user_id="usr_missing",
            proposal_id=SYNTHETIC_PROPOSAL_ID,
            now=now,
        )
    assert auth_exc.value.status_code == 401
    assert auth_exc.value.code == "AUTHENTICATION_REQUIRED"

    # Empty, duplicate, or unknown selected_item_ids -> 422 catalog_payload_invalid.
    for bad_selection in ([], ["item_1", "item_1"], ["item_unknown"]):
        with pytest.raises(CatalogPreviewError) as sel_exc:
            store.preview_proposal_subset(
                actor_user_id="usr_admin_a",
                proposal_id=SYNTHETIC_PROPOSAL_ID,
                selected_item_ids=bad_selection,
                now=now,
            )
        assert sel_exc.value.status_code == 422
        assert sel_exc.value.code == "catalog_payload_invalid"

    # Non-ordered selected_item_ids and mismatched expected_revisions fail in adopt.
    preview = store.preview_proposal_subset(
        actor_user_id="usr_admin_a",
        proposal_id=SYNTHETIC_PROPOSAL_ID,
        selected_item_ids=["item_1", "item_3"],
        now=now,
    )
    with pytest.raises(CatalogPreviewError) as order_exc:
        store.adopt_proposal_subset(
            actor_user_id="usr_admin_a",
            proposal_id=SYNTHETIC_PROPOSAL_ID,
            preview_digest=preview["preview_digest"],
            selected_item_ids=["item_3", "item_1"],
            expected_revisions=preview["expected_revisions"],
            idempotency_key="idem-bad-order",
            now=now,
        )
    assert order_exc.value.status_code == 409
    assert order_exc.value.code == "catalog_preview_digest_mismatch"

    # Term type mismatch and unknown term revision in a proposal item.
    bad_type_item = ProposalItem(
        item_id="item_bad_type",
        target_kind="experiment",
        target_id=SYNTHETIC_EXP_1_ID,
        base_revision_id="expr_44444444444444444444444444444441",
        base_revision_number=1,
        allowed_fields=("batch",),
        before={"batch": {"status": "unknown"}},
        after={
            "batch": {
                "status": "known",
                "term_id": "term_eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee",
                "term_revision_id": "trmv_ffffffffffffffffffffffffffffffff",
            }
        },
        citation=ProposalCitation(
            rule_version="m2-catalog-rules/1",
            vocabulary_version="terms-snapshot/2026-11-18",
            reason_codes=("BATCH_TERM_MATCH",),
            source_refs=(
                SourceReference(
                    kind="experiment_revision",
                    ref_id="expr_44444444444444444444444444444441",
                    evidence_token="type-mismatch-test",
                ),
            ),
        ),
    )
    bad_type_proposal = CatalogProposal.create(
        proposal_id="cprop_11111111111111111111111111111111",
        project_id=SYNTHETIC_PROJECT_ID,
        created_by_user_id="usr_member_01",
        source="offline_suggestion",
        items=(bad_type_item,),
        created_at=SYNTHETIC_BASE_TIME,
    )
    store._proposals[bad_type_proposal.proposal_id] = bad_type_proposal
    with pytest.raises(CatalogPreviewError) as type_exc:
        store.preview_proposal_subset(
            actor_user_id="usr_admin_a",
            proposal_id=bad_type_proposal.proposal_id,
            selected_item_ids=["item_bad_type"],
            now=now,
        )
    assert type_exc.value.status_code == 409
    assert type_exc.value.code == "catalog_term_type_mismatch"

    # Manual write replay and conflict also enforce idempotency and revision checks.
    manual_1 = store.manual_revise_experiment(
        actor_user_id="usr_admin_b",
        experiment_id=SYNTHETIC_EXP_2_ID,
        expected_revision_number=1,
        idempotency_key="idem-manual-exp2",
        updates={"display_name": "Treg H3K4me3 CUT&Tag v2"},
        now=now,
    )
    manual_replay = store.manual_revise_experiment(
        actor_user_id="usr_admin_b",
        experiment_id=SYNTHETIC_EXP_2_ID,
        expected_revision_number=1,
        idempotency_key="idem-manual-exp2",
        updates={"display_name": "Treg H3K4me3 CUT&Tag v2"},
        now=now,
    )
    assert manual_replay["replayed"] is True
    assert manual_replay["experiment_revision_id"] == manual_1["experiment_revision_id"]

    with pytest.raises(CatalogPreviewError) as manual_idem_conflict:
        store.manual_revise_experiment(
            actor_user_id="usr_admin_b",
            experiment_id=SYNTHETIC_EXP_2_ID,
            expected_revision_number=1,
            idempotency_key="idem-manual-exp2",
            updates={"display_name": "Different body"},
            now=now,
        )
    assert manual_idem_conflict.value.status_code == 409
    assert manual_idem_conflict.value.code == "catalog_idempotency_conflict"

    with pytest.raises(CatalogPreviewError) as manual_rev_conflict:
        store.manual_revise_experiment(
            actor_user_id="usr_admin_b",
            experiment_id=SYNTHETIC_EXP_2_ID,
            expected_revision_number=1,
            idempotency_key="idem-manual-exp2-stale",
            updates={"display_name": "Stale expected rev"},
            now=now,
        )
    assert manual_rev_conflict.value.status_code == 409
    assert manual_rev_conflict.value.code == "catalog_revision_conflict"

    # Archived experiment rejects preview and adoption.
    store._experiments[SYNTHETIC_EXP_1_ID].archived = True
    with pytest.raises(CatalogPreviewError) as arch_exp_exc:
        store.preview_proposal_subset(
            actor_user_id="usr_admin_a",
            proposal_id=SYNTHETIC_PROPOSAL_ID,
            selected_item_ids=["item_1"],
            now=now,
        )
    assert arch_exp_exc.value.status_code == 409
    assert arch_exp_exc.value.code == "catalog_experiment_archived"

    # Disallowed field on sample_catalog is rejected at ProposalItem construction.
    with pytest.raises(CatalogPreviewError) as bad_field_exc:
        ProposalItem(
            item_id="item_bad_field",
            target_kind="sample_catalog",
            target_id=SYNTHETIC_SMP_1_ID,
            base_revision_id="smpr_33333333333333333333333333333331",
            base_revision_number=1,
            allowed_fields=("assay",),
            before={"assay": {"status": "unknown"}},
            after={"assay": {"status": "not_applicable"}},
            citation=bad_type_item.citation,
        )
    assert bad_field_exc.value.status_code == 422
    assert bad_field_exc.value.code == "catalog_payload_invalid"
