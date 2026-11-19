"""Acceptance and contract tests for M2 offline catalog curation preview."""

from __future__ import annotations

from copy import deepcopy
from datetime import timedelta
import json
from pathlib import Path

import pytest

from encode_pipeline.platform.catalog_curation_preview import (
    EXPIRED_PROPOSAL_RETENTION,
    M2_RULE_VERSION,
    M2_VOCABULARY_VERSION,
    MAX_CHANGESET_ITEMS,
    PROPOSAL_ADOPTION_WINDOW,
    CatalogPreviewError,
    CatalogProposal,
    ProposalCitation,
    ProposalItem,
    SourceReference,
    SYNTHETIC_BASE_TIME,
    SYNTHETIC_EXP_1_ID,
    SYNTHETIC_EXP_2_ID,
    SYNTHETIC_OTHER_PROJECT_ID,
    SYNTHETIC_PROJECT_ID,
    SYNTHETIC_PROPOSAL_ID,
    SYNTHETIC_SECOND_PROPOSAL_ID,
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
    for case_id in (
        "M2-P1",
        "M2-P2",
        "M2-P3",
        "M2-P4",
        "M2-P5",
        "M2-P6",
        "M2-P7",
    ):
        assert case_id in spec_text

    committed_examples = json.loads(EXAMPLES_PATH.read_text(encoding="utf-8"))
    assert committed_examples == export_m2_preview_contract_examples()


def test_m2_p1_proposal_citations_and_role_projection() -> None:
    """M2-P1: Member view redacts ambiguous candidates AND their source_refs/tokens."""

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
    # Candidate term_revision_ids and identifying tokens must not leak to members.
    member_refs = item_3_member["citations"]["source_refs"]
    assert [ref["ref_id"] for ref in member_refs] == [
        "smpr_33333333333333333333333333333331",
        "trmv_10101010101010101010101010101010",
    ]
    serialized_member = json.dumps(member_view)
    for secret_token in (
        "trmv_bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb",
        "trmv_dddddddddddddddddddddddddddddddd",
        "alias-candidate-treg",
        "alias-candidate-kc",
        "aliases",
    ):
        assert secret_token not in serialized_member

    # Defense-in-depth: even if visibility is left default on a candidate ref, member view strips it.
    unmarked_ref_citation = ProposalCitation(
        rule_version=M2_RULE_VERSION,
        vocabulary_version=M2_VOCABULARY_VERSION,
        reason_codes=("ALIAS_AMBIGUOUS_CANDIDATES",),
        source_refs=(
            SourceReference(
                kind="sample_revision",
                ref_id="smpr_33333333333333333333333333333331",
                evidence_token="proj-smp-1-catalog-slots",
            ),
            SourceReference(
                kind="term_revision",
                ref_id="trmv_bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb",
                term_revision_id="trmv_bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb",
                evidence_token="unmarked-candidate-token",
            ),
        ),
    )
    projected_member_refs = unmarked_ref_citation.to_object(
        include_candidates=False,
        candidate_term_revision_ids=frozenset(
            {"trmv_bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb"}
        ),
    )["source_refs"]
    assert [ref["ref_id"] for ref in projected_member_refs] == [
        "smpr_33333333333333333333333333333331"
    ]

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
    assert item_3_admin["citations"]["rule_version"] == M2_RULE_VERSION
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

    exp1_head = store.get_experiment_head(SYNTHETIC_EXP_1_ID)
    assert exp1_head["revision_number"] == 2
    assert (
        exp1_head["payload"]["display_name"] == "Treg H3K27me3 CUT&Tag (Manual Rev 2)"
    )
    assert exp1_head["payload"]["batch"] == {"status": "unknown"}
    assert store.get_sample_head(SYNTHETIC_SMP_1_ID)["revision_number"] == 1


def test_m2_p4_single_adoption_outcome_for_disjoint_and_competing_proposals() -> None:
    """M2-P4: One proposal produces at most one adoption outcome; disjoint second confirm fails."""

    store = build_synthetic_m2_preview_fixture()
    t_preview = SYNTHETIC_BASE_TIME + timedelta(minutes=2)
    t_commit_a = SYNTHETIC_BASE_TIME + timedelta(minutes=3)
    t_commit_b = SYNTHETIC_BASE_TIME + timedelta(minutes=4)

    # Admin A previews {item_1, item_3}; Admin B previews disjoint {item_2} on the same proposal.
    preview_a = store.preview_proposal_subset(
        actor_user_id="usr_admin_a",
        proposal_id=SYNTHETIC_PROPOSAL_ID,
        selected_item_ids=["item_1", "item_3"],
        now=t_preview,
    )
    preview_b_disjoint = store.preview_proposal_subset(
        actor_user_id="usr_admin_b",
        proposal_id=SYNTHETIC_PROPOSAL_ID,
        selected_item_ids=["item_2"],
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

    # Reopening the proposal shows adopted state, adoptable_by_caller=False, and original receipt.
    reopened_by_b = store.get_proposal_view(
        actor_user_id="usr_admin_b",
        proposal_id=SYNTHETIC_PROPOSAL_ID,
        now=t_commit_b,
    )
    assert reopened_by_b["proposal_state"] == "adopted"
    assert reopened_by_b["adoptable_by_caller"] is False
    assert reopened_by_b["adoption_receipt"] == winner_receipt

    # Reopening after 7d TTL still preserves the adopted audit receipt.
    reopened_past_retention = store.get_proposal_view(
        actor_user_id="usr_admin_a",
        proposal_id=SYNTHETIC_PROPOSAL_ID,
        now=SYNTHETIC_BASE_TIME + EXPIRED_PROPOSAL_RETENTION + timedelta(days=1),
    )
    assert reopened_past_retention["proposal_state"] == "adopted"
    assert reopened_past_retention["adoption_receipt"] == winner_receipt

    # New preview on an adopted proposal is rejected.
    with pytest.raises(CatalogPreviewError) as re_preview_exc:
        store.preview_proposal_subset(
            actor_user_id="usr_admin_b",
            proposal_id=SYNTHETIC_PROPOSAL_ID,
            selected_item_ids=["item_2"],
            now=t_commit_b,
        )
    assert re_preview_exc.value.status_code == 409
    assert re_preview_exc.value.code == "catalog_proposal_invalidated"

    # Admin B's stale preview for disjoint item_2 cannot write a second receipt.
    with pytest.raises(CatalogPreviewError) as disjoint_exc:
        store.adopt_proposal_subset(
            actor_user_id="usr_admin_b",
            proposal_id=SYNTHETIC_PROPOSAL_ID,
            preview_digest=preview_b_disjoint["preview_digest"],
            selected_item_ids=preview_b_disjoint["selected_item_ids"],
            expected_revisions=preview_b_disjoint["expected_revisions"],
            idempotency_key="idem-race-admin-b-disjoint",
            now=t_commit_b,
        )
    assert disjoint_exc.value.status_code == 409
    assert disjoint_exc.value.code == "catalog_proposal_invalidated"
    assert store.get_experiment_head(SYNTHETIC_EXP_2_ID)["revision_number"] == 1
    assert (
        len(
            store.list_adoption_receipts(
                actor_user_id="usr_admin_a", project_id=SYNTHETIC_PROJECT_ID
            )["receipts"]
        )
        == 1
    )

    # Admin A replaying the original key and body still receives the original receipt.
    replayed = store.adopt_proposal_subset(
        actor_user_id="usr_admin_a",
        proposal_id=SYNTHETIC_PROPOSAL_ID,
        preview_digest=preview_a["preview_digest"],
        selected_item_ids=preview_a["selected_item_ids"],
        expected_revisions=preview_a["expected_revisions"],
        idempotency_key="idem-race-admin-a",
        now=t_commit_b,
    )
    assert replayed["replayed"] is True
    assert replayed["receipt_id"] == winner_receipt["receipt_id"]

    # A separate proposal targeting exp_1 at expected_revision=1 hits 409 catalog_revision_conflict.
    other_prop = CatalogProposal.create(
        proposal_id="cprop_33333333333333333333333333333333",
        project_id=SYNTHETIC_PROJECT_ID,
        created_by_user_id="usr_member_01",
        source="offline_suggestion",
        items=(
            store._proposals[SYNTHETIC_PROPOSAL_ID].items[0],
            store._proposals[SYNTHETIC_PROPOSAL_ID].items[1],
        ),
        created_at=SYNTHETIC_BASE_TIME,
    )
    store._proposals[other_prop.proposal_id] = other_prop
    other_preview = store.preview_proposal_subset(
        actor_user_id="usr_admin_b",
        proposal_id=other_prop.proposal_id,
        selected_item_ids=["item_1", "item_2"],
        now=t_commit_b,
    )
    assert other_preview["has_conflicts"] is True
    with pytest.raises(CatalogPreviewError) as other_conflict:
        store.adopt_proposal_subset(
            actor_user_id="usr_admin_b",
            proposal_id=other_prop.proposal_id,
            preview_digest=other_preview["preview_digest"],
            selected_item_ids=other_preview["selected_item_ids"],
            expected_revisions=other_preview["expected_revisions"],
            idempotency_key="idem-other-prop-stale",
            now=t_commit_b,
        )
    assert other_conflict.value.status_code == 409
    assert other_conflict.value.code == "catalog_revision_conflict"
    assert store.get_experiment_head(SYNTHETIC_EXP_2_ID)["revision_number"] == 1


def test_m2_p5_dropped_response_replay_and_m1a_idempotency_key_scope() -> None:
    """M2-P5: (actor, key) scope replays exact request and rejects reuse across proposals/bodies."""

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
    second_preview = store.preview_proposal_subset(
        actor_user_id="usr_admin_a",
        proposal_id=SYNTHETIC_SECOND_PROPOSAL_ID,
        selected_item_ids=["item_2"],
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

    # Dropped response retry: succeeds cleanly with replayed=True even though proposal is adopted.
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

    # Same admin reusing the exact same Idempotency-Key on a SECOND proposal fails with 409 and zero writes.
    with pytest.raises(CatalogPreviewError) as cross_prop_exc:
        store.adopt_proposal_subset(
            actor_user_id="usr_admin_a",
            proposal_id=SYNTHETIC_SECOND_PROPOSAL_ID,
            preview_digest=second_preview["preview_digest"],
            selected_item_ids=second_preview["selected_item_ids"],
            expected_revisions=second_preview["expected_revisions"],
            idempotency_key="idem-dropped-response",
            now=t_retry,
        )
    assert cross_prop_exc.value.status_code == 409
    assert cross_prop_exc.value.code == "catalog_idempotency_conflict"
    assert store.get_experiment_head(SYNTHETIC_EXP_2_ID)["revision_number"] == 1

    # Same admin reusing same key with different body on same proposal -> 409.
    with pytest.raises(CatalogPreviewError) as diff_body_exc:
        store.adopt_proposal_subset(
            actor_user_id="usr_admin_a",
            proposal_id=SYNTHETIC_PROPOSAL_ID,
            preview_digest=preview["preview_digest"],
            selected_item_ids=["item_1"],
            expected_revisions={SYNTHETIC_EXP_1_ID: 1},
            idempotency_key="idem-dropped-response",
            now=t_retry,
        )
    assert diff_body_exc.value.status_code == 409
    assert diff_body_exc.value.code == "catalog_idempotency_conflict"

    # Receipt discovery after re-login finds only the original receipt.
    discovered = store.list_adoption_receipts(
        actor_user_id="usr_admin_a",
        project_id=SYNTHETIC_PROJECT_ID,
        proposal_id=SYNTHETIC_PROPOSAL_ID,
    )
    assert len(discovered["receipts"]) == 1
    assert discovered["receipts"][0]["receipt_id"] == first_receipt["receipt_id"]

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


def test_m2_p6_baseline_revision_id_provenance_and_version_checks() -> None:
    """M2-P6: Verify base_revision_id + number, source_refs existence/ownership, and rule/vocab versions."""

    now = SYNTHETIC_BASE_TIME + timedelta(minutes=2)

    # Case (a): Same revision number (1), but base_revision_id references a different revision!
    mismatched_id_store = build_synthetic_m2_preview_fixture()
    base_items = mismatched_id_store._proposals[SYNTHETIC_PROPOSAL_ID].items
    wrong_rev_id_item_1 = ProposalItem(
        item_id="item_1",
        target_kind="experiment",
        target_id=SYNTHETIC_EXP_1_ID,
        base_revision_id="expr_99999999999999999999999999999999",
        base_revision_number=1,
        allowed_fields=base_items[0].allowed_fields,
        before=base_items[0].before,
        after=base_items[0].after,
        citation=ProposalCitation(
            rule_version=M2_RULE_VERSION,
            vocabulary_version=M2_VOCABULARY_VERSION,
            reason_codes=("BATCH_TERM_MATCH",),
            source_refs=(
                SourceReference(
                    kind="experiment_revision",
                    ref_id="expr_44444444444444444444444444444441",
                    evidence_token="valid-cited-rev",
                ),
            ),
        ),
    )
    wrong_id_proposal = CatalogProposal.create(
        proposal_id="cprop_44444444444444444444444444444444",
        project_id=SYNTHETIC_PROJECT_ID,
        created_by_user_id="usr_member_01",
        source="offline_suggestion",
        items=(wrong_rev_id_item_1, base_items[1]),
        created_at=SYNTHETIC_BASE_TIME,
    )
    mismatched_id_store._proposals[wrong_id_proposal.proposal_id] = wrong_id_proposal
    preview_wrong_id = mismatched_id_store.preview_proposal_subset(
        actor_user_id="usr_admin_a",
        proposal_id=wrong_id_proposal.proposal_id,
        selected_item_ids=["item_1", "item_2"],
        now=now,
    )
    assert preview_wrong_id["has_conflicts"] is True
    assert preview_wrong_id["conflicts"][0]["expected_revision_number"] == 1
    assert preview_wrong_id["conflicts"][0]["current_revision_number"] == 1
    assert (
        preview_wrong_id["conflicts"][0]["expected_revision_id"]
        == "expr_99999999999999999999999999999999"
    )
    assert (
        preview_wrong_id["conflicts"][0]["current_revision_id"]
        == "expr_44444444444444444444444444444441"
    )

    with pytest.raises(CatalogPreviewError) as wrong_id_adopt_exc:
        mismatched_id_store.adopt_proposal_subset(
            actor_user_id="usr_admin_a",
            proposal_id=wrong_id_proposal.proposal_id,
            preview_digest=preview_wrong_id["preview_digest"],
            selected_item_ids=preview_wrong_id["selected_item_ids"],
            expected_revisions=preview_wrong_id["expected_revisions"],
            idempotency_key="idem-wrong-rev-id",
            now=now,
        )
    assert wrong_id_adopt_exc.value.status_code == 409
    assert wrong_id_adopt_exc.value.code == "catalog_revision_conflict"
    # Zero partial writes across the batch (item_2 also stays at rev 1).
    assert (
        mismatched_id_store.get_experiment_head(SYNTHETIC_EXP_1_ID)["revision_number"]
        == 1
    )
    assert (
        mismatched_id_store.get_experiment_head(SYNTHETIC_EXP_2_ID)["revision_number"]
        == 1
    )

    # Forged caller expected_revisions + digest cannot override server proposal baseline.
    mismatched_id_store.manual_revise_experiment(
        actor_user_id="usr_admin_b",
        experiment_id=SYNTHETIC_EXP_2_ID,
        expected_revision_number=1,
        idempotency_key="idem-advance-exp2",
        updates={"display_name": "Exp 2 Rev 2"},
        now=now,
    )
    forged_expected = {SYNTHETIC_EXP_2_ID: 2}
    from encode_pipeline.platform.catalog_curation_preview import compute_preview_digest

    forged_digest = compute_preview_digest(
        proposal_id=SYNTHETIC_SECOND_PROPOSAL_ID,
        proposal_digest=mismatched_id_store._proposals[
            SYNTHETIC_SECOND_PROPOSAL_ID
        ].proposal_digest,
        selected_item_ids=["item_2"],
        expected_revisions=forged_expected,
    )
    with pytest.raises(CatalogPreviewError) as forged_exc:
        mismatched_id_store.adopt_proposal_subset(
            actor_user_id="usr_admin_a",
            proposal_id=SYNTHETIC_SECOND_PROPOSAL_ID,
            preview_digest=forged_digest,
            selected_item_ids=["item_2"],
            expected_revisions=forged_expected,
            idempotency_key="idem-forged-expected",
            now=now,
        )
    assert forged_exc.value.status_code == 409
    assert forged_exc.value.code == "catalog_proposal_invalidated"
    assert (
        mismatched_id_store.get_experiment_head(SYNTHETIC_EXP_2_ID)["revision_number"]
        == 2
    )

    # Case (b): Invalidated / missing source_ref fails closed with zero writes.
    stale_ref_store = build_synthetic_m2_preview_fixture()
    bad_ref_item_1 = ProposalItem(
        item_id="item_1",
        target_kind="experiment",
        target_id=SYNTHETIC_EXP_1_ID,
        base_revision_id="expr_44444444444444444444444444444441",
        base_revision_number=1,
        allowed_fields=base_items[0].allowed_fields,
        before=base_items[0].before,
        after=base_items[0].after,
        citation=ProposalCitation(
            rule_version=M2_RULE_VERSION,
            vocabulary_version=M2_VOCABULARY_VERSION,
            reason_codes=("BATCH_TERM_MATCH",),
            source_refs=(
                SourceReference(
                    kind="term_revision",
                    ref_id="trmv_00000000000000000000000000000000",
                    term_revision_id="trmv_00000000000000000000000000000000",
                    evidence_token="missing-term-rev",
                ),
            ),
        ),
    )
    bad_ref_proposal = CatalogProposal.create(
        proposal_id="cprop_55555555555555555555555555555555",
        project_id=SYNTHETIC_PROJECT_ID,
        created_by_user_id="usr_member_01",
        source="offline_suggestion",
        items=(bad_ref_item_1, base_items[1]),
        created_at=SYNTHETIC_BASE_TIME,
    )
    stale_ref_store._proposals[bad_ref_proposal.proposal_id] = bad_ref_proposal
    with pytest.raises(CatalogPreviewError) as bad_ref_preview_exc:
        stale_ref_store.preview_proposal_subset(
            actor_user_id="usr_admin_a",
            proposal_id=bad_ref_proposal.proposal_id,
            selected_item_ids=["item_1", "item_2"],
            now=now,
        )
    assert bad_ref_preview_exc.value.status_code == 409
    assert bad_ref_preview_exc.value.code == "catalog_proposal_invalidated"

    bad_ref_digest = compute_preview_digest(
        proposal_id=bad_ref_proposal.proposal_id,
        proposal_digest=bad_ref_proposal.proposal_digest,
        selected_item_ids=["item_1", "item_2"],
        expected_revisions={SYNTHETIC_EXP_1_ID: 1, SYNTHETIC_EXP_2_ID: 1},
    )
    with pytest.raises(CatalogPreviewError) as bad_ref_adopt_exc:
        stale_ref_store.adopt_proposal_subset(
            actor_user_id="usr_admin_a",
            proposal_id=bad_ref_proposal.proposal_id,
            preview_digest=bad_ref_digest,
            selected_item_ids=["item_1", "item_2"],
            expected_revisions={SYNTHETIC_EXP_1_ID: 1, SYNTHETIC_EXP_2_ID: 1},
            idempotency_key="idem-bad-ref-adopt",
            now=now,
        )
    assert bad_ref_adopt_exc.value.status_code == 409
    assert bad_ref_adopt_exc.value.code == "catalog_proposal_invalidated"
    assert (
        stale_ref_store.get_experiment_head(SYNTHETIC_EXP_1_ID)["revision_number"] == 1
    )
    assert (
        stale_ref_store.get_experiment_head(SYNTHETIC_EXP_2_ID)["revision_number"] == 1
    )

    # Case (c): Advancing active rule/vocabulary version invalidates old proposals.
    version_store = build_synthetic_m2_preview_fixture()
    valid_preview = version_store.preview_proposal_subset(
        actor_user_id="usr_admin_a",
        proposal_id=SYNTHETIC_PROPOSAL_ID,
        selected_item_ids=["item_1", "item_2"],
        now=now,
    )
    version_store.set_active_versions(vocabulary_version="terms-snapshot/2026-11-19")
    with pytest.raises(CatalogPreviewError) as vocab_exc:
        version_store.adopt_proposal_subset(
            actor_user_id="usr_admin_a",
            proposal_id=SYNTHETIC_PROPOSAL_ID,
            preview_digest=valid_preview["preview_digest"],
            selected_item_ids=valid_preview["selected_item_ids"],
            expected_revisions=valid_preview["expected_revisions"],
            idempotency_key="idem-stale-vocab",
            now=now,
        )
    assert vocab_exc.value.status_code == 409
    assert vocab_exc.value.code == "catalog_proposal_invalidated"
    assert version_store.get_experiment_head(SYNTHETIC_EXP_1_ID)["revision_number"] == 1

    # Case (d): Cross-project target ownership mismatch fails with 404.
    cross_proj_store = build_synthetic_m2_preview_fixture()
    cross_proj_store._experiments[
        SYNTHETIC_EXP_1_ID
    ].project_id = SYNTHETIC_OTHER_PROJECT_ID
    with pytest.raises(CatalogPreviewError) as cross_proj_exc:
        cross_proj_store.preview_proposal_subset(
            actor_user_id="usr_admin_a",
            proposal_id=SYNTHETIC_PROPOSAL_ID,
            selected_item_ids=["item_1"],
            now=now,
        )
    assert cross_proj_exc.value.status_code == 404
    assert cross_proj_exc.value.code == "catalog_not_found"


def test_m2_p7_expiry_revocation_archived_terms_and_egress_guards() -> None:
    """M2-P7: Expiry, admin revocation, archived terms, and invalid citations fail closed."""

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

    # Admin downgraded to member before commit -> 403 and zero writes.
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

    # Archived term referenced by proposal fails with 409 catalog_term_archived.
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

    # Host paths and >100-item changesets are rejected by proposal construction.
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
                rule_version=M2_RULE_VERSION,
                vocabulary_version=M2_VOCABULARY_VERSION,
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

    with pytest.raises(CatalogPreviewError) as auth_exc:
        store.get_proposal_view(
            actor_user_id="usr_missing",
            proposal_id=SYNTHETIC_PROPOSAL_ID,
            now=now,
        )
    assert auth_exc.value.status_code == 401
    assert auth_exc.value.code == "AUTHENTICATION_REQUIRED"

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
            rule_version=M2_RULE_VERSION,
            vocabulary_version=M2_VOCABULARY_VERSION,
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


def test_m2_projection_row_source_ref_ownership_and_archival_checks() -> None:
    """Cited catalog_projection_row refs must belong to the active project and not be archived."""

    from encode_pipeline.platform.catalog_curation_preview import compute_preview_digest

    now = SYNTHETIC_BASE_TIME + timedelta(minutes=2)

    def build_exp2_proposal_with_projection_ref(
        proposal_id: str, cited_ref_id: str
    ) -> CatalogProposal:
        base_item_2 = (
            build_synthetic_m2_preview_fixture()
            ._proposals[SYNTHETIC_PROPOSAL_ID]
            .items[1]
        )
        item = ProposalItem(
            item_id="item_proj_ref",
            target_kind="experiment",
            target_id=SYNTHETIC_EXP_2_ID,
            base_revision_id=base_item_2.base_revision_id,
            base_revision_number=base_item_2.base_revision_number,
            allowed_fields=base_item_2.allowed_fields,
            before=base_item_2.before,
            after=base_item_2.after,
            citation=ProposalCitation(
                rule_version=M2_RULE_VERSION,
                vocabulary_version=M2_VOCABULARY_VERSION,
                reason_codes=("TARGET_TERM_MATCH",),
                source_refs=(
                    SourceReference(
                        kind="catalog_projection_row",
                        ref_id=cited_ref_id,
                        evidence_token="proj-row-source-check",
                    ),
                ),
            ),
        )
        return CatalogProposal.create(
            proposal_id=proposal_id,
            project_id=SYNTHETIC_PROJECT_ID,
            created_by_user_id="usr_member_01",
            source="offline_suggestion",
            items=(item,),
            created_at=SYNTHETIC_BASE_TIME,
        )

    # 1. Positive case: active same-project catalog_projection_row succeeds in preview and adoption.
    ok_store = build_synthetic_m2_preview_fixture()
    ok_prop = build_exp2_proposal_with_projection_ref(
        "cprop_66666666666666666666666666666666", SYNTHETIC_EXP_1_ID
    )
    ok_store._proposals[ok_prop.proposal_id] = ok_prop
    ok_preview = ok_store.preview_proposal_subset(
        actor_user_id="usr_admin_a",
        proposal_id=ok_prop.proposal_id,
        selected_item_ids=["item_proj_ref"],
        now=now,
    )
    ok_receipt = ok_store.adopt_proposal_subset(
        actor_user_id="usr_admin_a",
        proposal_id=ok_prop.proposal_id,
        preview_digest=ok_preview["preview_digest"],
        selected_item_ids=ok_preview["selected_item_ids"],
        expected_revisions=ok_preview["expected_revisions"],
        idempotency_key="idem-proj-row-ok",
        now=now,
    )
    assert ok_receipt["replayed"] is False
    assert ok_store.get_experiment_head(SYNTHETIC_EXP_2_ID)["revision_number"] == 2

    # 2. Cross-project experiment projection row: target exp_2 is in SYNTHETIC_PROJECT_ID,
    #    but cited exp_1 moved to SYNTHETIC_OTHER_PROJECT_ID -> rejected in preview & adopt.
    cross_exp_store = build_synthetic_m2_preview_fixture()
    cross_exp_prop = build_exp2_proposal_with_projection_ref(
        "cprop_77777777777777777777777777777777", SYNTHETIC_EXP_1_ID
    )
    cross_exp_store._proposals[cross_exp_prop.proposal_id] = cross_exp_prop
    cross_exp_digest = compute_preview_digest(
        proposal_id=cross_exp_prop.proposal_id,
        proposal_digest=cross_exp_prop.proposal_digest,
        selected_item_ids=["item_proj_ref"],
        expected_revisions={SYNTHETIC_EXP_2_ID: 1},
    )
    cross_exp_store._experiments[
        SYNTHETIC_EXP_1_ID
    ].project_id = SYNTHETIC_OTHER_PROJECT_ID
    with pytest.raises(CatalogPreviewError) as cross_exp_preview_exc:
        cross_exp_store.preview_proposal_subset(
            actor_user_id="usr_admin_a",
            proposal_id=cross_exp_prop.proposal_id,
            selected_item_ids=["item_proj_ref"],
            now=now,
        )
    assert cross_exp_preview_exc.value.status_code == 409
    assert cross_exp_preview_exc.value.code == "catalog_proposal_invalidated"

    with pytest.raises(CatalogPreviewError) as cross_exp_adopt_exc:
        cross_exp_store.adopt_proposal_subset(
            actor_user_id="usr_admin_a",
            proposal_id=cross_exp_prop.proposal_id,
            preview_digest=cross_exp_digest,
            selected_item_ids=["item_proj_ref"],
            expected_revisions={SYNTHETIC_EXP_2_ID: 1},
            idempotency_key="idem-proj-row-cross-exp",
            now=now,
        )
    assert cross_exp_adopt_exc.value.status_code == 409
    assert cross_exp_adopt_exc.value.code == "catalog_proposal_invalidated"
    assert (
        cross_exp_store.get_experiment_head(SYNTHETIC_EXP_2_ID)["revision_number"] == 1
    )
    assert (
        cross_exp_store.list_adoption_receipts(
            actor_user_id="usr_admin_a", project_id=SYNTHETIC_PROJECT_ID
        )["receipts"]
        == []
    )

    # 3. Cross-project sample projection row -> rejected in preview & adopt, zero writes/receipts.
    cross_smp_store = build_synthetic_m2_preview_fixture()
    cross_smp_prop = build_exp2_proposal_with_projection_ref(
        "cprop_88888888888888888888888888888888", SYNTHETIC_SMP_1_ID
    )
    cross_smp_store._proposals[cross_smp_prop.proposal_id] = cross_smp_prop
    cross_smp_digest = compute_preview_digest(
        proposal_id=cross_smp_prop.proposal_id,
        proposal_digest=cross_smp_prop.proposal_digest,
        selected_item_ids=["item_proj_ref"],
        expected_revisions={SYNTHETIC_EXP_2_ID: 1},
    )
    cross_smp_store._samples[SYNTHETIC_SMP_1_ID].project_id = SYNTHETIC_OTHER_PROJECT_ID
    for call_fn in (
        lambda: cross_smp_store.preview_proposal_subset(
            actor_user_id="usr_admin_a",
            proposal_id=cross_smp_prop.proposal_id,
            selected_item_ids=["item_proj_ref"],
            now=now,
        ),
        lambda: cross_smp_store.adopt_proposal_subset(
            actor_user_id="usr_admin_a",
            proposal_id=cross_smp_prop.proposal_id,
            preview_digest=cross_smp_digest,
            selected_item_ids=["item_proj_ref"],
            expected_revisions={SYNTHETIC_EXP_2_ID: 1},
            idempotency_key="idem-proj-row-cross-smp",
            now=now,
        ),
    ):
        with pytest.raises(CatalogPreviewError) as exc:
            call_fn()
        assert exc.value.status_code == 409
        assert exc.value.code == "catalog_proposal_invalidated"
    assert (
        cross_smp_store.get_experiment_head(SYNTHETIC_EXP_2_ID)["revision_number"] == 1
    )
    assert (
        cross_smp_store.list_adoption_receipts(
            actor_user_id="usr_admin_a", project_id=SYNTHETIC_PROJECT_ID
        )["receipts"]
        == []
    )

    # 4. Archived same-project experiment projection row -> rejected in preview & adopt, zero writes.
    archived_ref_store = build_synthetic_m2_preview_fixture()
    archived_ref_prop = build_exp2_proposal_with_projection_ref(
        "cprop_99999999999999999999999999999999", SYNTHETIC_EXP_1_ID
    )
    archived_ref_store._proposals[archived_ref_prop.proposal_id] = archived_ref_prop
    archived_ref_digest = compute_preview_digest(
        proposal_id=archived_ref_prop.proposal_id,
        proposal_digest=archived_ref_prop.proposal_digest,
        selected_item_ids=["item_proj_ref"],
        expected_revisions={SYNTHETIC_EXP_2_ID: 1},
    )
    archived_ref_store._experiments[SYNTHETIC_EXP_1_ID].archived = True
    for call_fn in (
        lambda: archived_ref_store.preview_proposal_subset(
            actor_user_id="usr_admin_a",
            proposal_id=archived_ref_prop.proposal_id,
            selected_item_ids=["item_proj_ref"],
            now=now,
        ),
        lambda: archived_ref_store.adopt_proposal_subset(
            actor_user_id="usr_admin_a",
            proposal_id=archived_ref_prop.proposal_id,
            preview_digest=archived_ref_digest,
            selected_item_ids=["item_proj_ref"],
            expected_revisions={SYNTHETIC_EXP_2_ID: 1},
            idempotency_key="idem-proj-row-archived",
            now=now,
        ),
    ):
        with pytest.raises(CatalogPreviewError) as exc:
            call_fn()
        assert exc.value.status_code == 409
        assert exc.value.code == "catalog_experiment_archived"
    assert (
        archived_ref_store.get_experiment_head(SYNTHETIC_EXP_2_ID)["revision_number"]
        == 1
    )
    assert (
        archived_ref_store.list_adoption_receipts(
            actor_user_id="usr_admin_a", project_id=SYNTHETIC_PROJECT_ID
        )["receipts"]
        == []
    )


_SMP_REV = "smpr_33333333333333333333333333333331"
_CONTROL_TERM = "term_78787878787878787878787878787878"
_CONTROL_REV = "trmv_90909090909090909090909090909090"
_TREATMENT_REV = "trmv_89898989898989898989898989898989"


def _experiment_citation(
    *refs: SourceReference, reasons: tuple[str, ...]
) -> ProposalCitation:
    return ProposalCitation(
        rule_version=M2_RULE_VERSION,
        vocabulary_version=M2_VOCABULARY_VERSION,
        reason_codes=reasons,
        source_refs=refs,
    )


def _install_experiment_proposal(
    store,
    *,
    proposal_id: str,
    allowed_fields: tuple[str, ...],
    before: dict,
    after: dict,
    reasons: tuple[str, ...] = ("REPLICATE_ROLE_TERM_MATCH",),
) -> CatalogProposal:
    head = store.get_experiment_head(SYNTHETIC_EXP_2_ID)
    proposal = CatalogProposal.create(
        proposal_id=proposal_id,
        project_id=SYNTHETIC_PROJECT_ID,
        created_by_user_id="usr_member_01",
        source="offline_suggestion",
        created_at=SYNTHETIC_BASE_TIME,
        items=(
            ProposalItem(
                item_id="item_samples",
                target_kind="experiment",
                target_id=SYNTHETIC_EXP_2_ID,
                base_revision_id=head["experiment_revision_id"],
                base_revision_number=head["revision_number"],
                allowed_fields=allowed_fields,
                before=before,
                after=after,
                citation=_experiment_citation(
                    SourceReference(
                        kind="experiment_revision",
                        ref_id=head["experiment_revision_id"],
                        evidence_token="proj-exp-2-sample-pins",
                    ),
                    SourceReference(
                        kind="term_revision",
                        ref_id=_CONTROL_REV,
                        term_revision_id=_CONTROL_REV,
                        evidence_token="vocab-replicate-control",
                    ),
                    reasons=reasons,
                ),
            ),
        ),
    )
    store._proposals[proposal.proposal_id] = proposal
    return proposal


def test_m2_sample_pins_stay_distinct_from_terms_and_adoption_is_normalized() -> None:
    """Sample associations are not term slots, and adoption stores only verified values."""

    from encode_pipeline.platform.catalog_curation_preview import (
        _collect_ambiguous_candidate_revision_ids,
        _validate_and_render_slot,
        compute_preview_digest,
    )

    now = SYNTHETIC_BASE_TIME + timedelta(minutes=2)
    candidate_ids = _collect_ambiguous_candidate_revision_ids(
        {
            "samples": [
                {
                    "replicate_role": {
                        "status": "ambiguous",
                        "candidates": [
                            {
                                "term_revision_id": "trmv_bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb"
                            }
                        ],
                    }
                }
            ]
        }
    )
    assert candidate_ids == frozenset({"trmv_bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb"})

    for bad_samples in (
        [{"ordinal": 0, "sample_revision_id": _SMP_REV}],
        {"status": "unknown"},
    ):
        with pytest.raises(CatalogPreviewError) as slot_exc:
            _validate_and_render_slot(
                "samples",
                bad_samples,
                {},
                include_candidates=True,
                for_write=True,
            )
        assert slot_exc.value.status_code == 422
        assert slot_exc.value.code == "catalog_payload_invalid"
    with pytest.raises(CatalogPreviewError) as unknown_field_exc:
        _validate_and_render_slot(
            "notes",
            {"status": "unknown"},
            {},
            include_candidates=True,
            for_write=True,
        )
    assert unknown_field_exc.value.code == "catalog_payload_invalid"

    with pytest.raises(CatalogPreviewError) as term_shaped_exc:
        ProposalItem(
            item_id="item_term_shaped_samples",
            target_kind="experiment",
            target_id=SYNTHETIC_EXP_2_ID,
            base_revision_id="expr_aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa1",
            base_revision_number=1,
            allowed_fields=("samples",),
            before={"samples": {"status": "unknown"}},
            after={"samples": {"status": "not_applicable"}},
            citation=_experiment_citation(
                SourceReference(
                    kind="experiment_revision",
                    ref_id="expr_aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa1",
                    evidence_token="proj-exp-2-sample-pins",
                ),
                reasons=("REPLICATE_ROLE_TERM_MATCH",),
            ),
        )
    assert term_shaped_exc.value.code == "catalog_payload_invalid"

    manual_store = build_synthetic_m2_preview_fixture()
    original_exp2 = manual_store.get_experiment_head(SYNTHETIC_EXP_2_ID)
    with pytest.raises(CatalogPreviewError) as manual_term_exc:
        manual_store.manual_revise_experiment(
            actor_user_id="usr_admin_a",
            experiment_id=SYNTHETIC_EXP_2_ID,
            expected_revision_number=1,
            idempotency_key="idem-samples-as-term",
            updates={"samples": {"status": "unknown"}},
            now=now,
        )
    assert manual_term_exc.value.code == "catalog_payload_invalid"
    assert manual_store.get_experiment_head(SYNTHETIC_EXP_2_ID) == original_exp2

    # Legal pin change: keep the sample revision and group, retarget replicate_role.
    role_store = build_synthetic_m2_preview_fixture()
    original_samples = deepcopy(
        role_store.get_experiment_head(SYNTHETIC_EXP_2_ID)["payload"]["samples"]
    )
    original_payload = deepcopy(
        role_store.get_experiment_head(SYNTHETIC_EXP_2_ID)["payload"]
    )
    changed_pin = deepcopy(original_samples[0])
    changed_pin["replicate_role"] = {
        "status": "known",
        "term_id": _CONTROL_TERM,
        "term_revision_id": _CONTROL_REV,
    }
    role_proposal = _install_experiment_proposal(
        role_store,
        proposal_id="cprop_aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
        allowed_fields=("samples",),
        before={"samples": original_samples},
        after={"samples": [changed_pin]},
    )
    role_preview = role_store.preview_proposal_subset(
        actor_user_id="usr_admin_a",
        proposal_id=role_proposal.proposal_id,
        selected_item_ids=["item_samples"],
        now=now,
    )
    preview_pin = role_preview["items"][0]["after"]["samples"][0]
    assert preview_pin["sample_revision_id"] == _SMP_REV
    assert preview_pin["replicate_group"] == "biorep-1"
    assert preview_pin["replicate_role"]["display_name"] == "control"
    role_store.adopt_proposal_subset(
        actor_user_id="usr_admin_a",
        proposal_id=role_proposal.proposal_id,
        preview_digest=role_preview["preview_digest"],
        selected_item_ids=role_preview["selected_item_ids"],
        expected_revisions=role_preview["expected_revisions"],
        idempotency_key="idem-sample-pin-role",
        now=now,
    )
    adopted_payload = role_store.get_experiment_head(SYNTHETIC_EXP_2_ID)["payload"]
    adopted_pin = adopted_payload["samples"][0]
    assert adopted_pin["ordinal"] == 0
    assert adopted_pin["sample_revision_id"] == _SMP_REV
    assert adopted_pin["replicate_group"] == "biorep-1"
    assert adopted_pin["replicate_role"] == {
        "status": "known",
        "term_id": _CONTROL_TERM,
        "term_revision_id": _CONTROL_REV,
    }
    assert "display_name" not in adopted_pin["replicate_role"]
    for unchanged in ("assay", "target", "batch", "display_name", "schema_version"):
        assert adopted_payload[unchanged] == original_payload[unchanged]

    # Classification adoption keeps every existing sample pin and drops client extras.
    align_store = build_synthetic_m2_preview_fixture()
    untouched_samples = deepcopy(
        align_store.get_experiment_head(SYNTHETIC_EXP_1_ID)["payload"]["samples"]
    )
    align_preview = align_store.preview_proposal_subset(
        actor_user_id="usr_admin_a",
        proposal_id=SYNTHETIC_PROPOSAL_ID,
        selected_item_ids=["item_1", "item_3"],
        now=now,
    )
    align_store.adopt_proposal_subset(
        actor_user_id="usr_admin_a",
        proposal_id=SYNTHETIC_PROPOSAL_ID,
        preview_digest=align_preview["preview_digest"],
        selected_item_ids=align_preview["selected_item_ids"],
        expected_revisions=align_preview["expected_revisions"],
        idempotency_key="idem-normalized-adopt",
        now=now,
    )
    exp1 = align_store.get_experiment_head(SYNTHETIC_EXP_1_ID)
    sample = align_store.get_sample_head(SYNTHETIC_SMP_1_ID)
    preview_by_id = {item["item_id"]: item for item in align_preview["items"]}
    stored_batch = exp1["payload"]["batch"]
    preview_batch = preview_by_id["item_1"]["after"]["batch"]
    assert stored_batch == {
        "status": "known",
        "term_id": preview_batch["term_id"],
        "term_revision_id": preview_batch["term_revision_id"],
    }
    assert preview_batch["display_name"] == "batch-2026-11"
    assert "display_name" not in stored_batch
    assert exp1["payload"]["samples"] == untouched_samples
    stored_tag = sample["catalog"]["tag"]
    preview_tag = preview_by_id["item_3"]["after"]["tag"]
    assert stored_tag["status"] == "known"
    assert stored_tag["values"] == [
        {
            "term_id": preview_tag["values"][0]["term_id"],
            "term_revision_id": preview_tag["values"][0]["term_revision_id"],
        }
    ]
    assert preview_tag["values"][0]["display_name"] == "m2-curated"
    assert "aliases" not in json.dumps(stored_tag)
    assert sample["extra_attributes"] == {"legacy_note": "preserve-me-outside-catalog"}
    assert sample["catalog"]["species"]["term_revision_id"] == (
        "trmv_67676767676767676767676767676767"
    )

    poison_store = build_synthetic_m2_preview_fixture()
    base_item = poison_store._proposals[SYNTHETIC_PROPOSAL_ID].items[2]
    poisoned_after = deepcopy(base_item.after)
    poisoned_after["tag"]["values"][0] = {
        **poisoned_after["tag"]["values"][0],
        "aliases": ["secret-alias"],
        "note": "donor-17",
    }
    poison = CatalogProposal.create(
        proposal_id="cprop_bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb",
        project_id=SYNTHETIC_PROJECT_ID,
        created_by_user_id="usr_member_01",
        source="offline_suggestion",
        created_at=SYNTHETIC_BASE_TIME,
        items=(
            ProposalItem(
                item_id="item_poison",
                target_kind="sample_catalog",
                target_id=base_item.target_id,
                base_revision_id=base_item.base_revision_id,
                base_revision_number=base_item.base_revision_number,
                allowed_fields=base_item.allowed_fields,
                before=deepcopy(base_item.before),
                after=poisoned_after,
                citation=base_item.citation,
            ),
        ),
    )
    poison_store._proposals[poison.proposal_id] = poison
    poison_digest = compute_preview_digest(
        proposal_id=poison.proposal_id,
        proposal_digest=poison.proposal_digest,
        selected_item_ids=["item_poison"],
        expected_revisions={SYNTHETIC_SMP_1_ID: 1},
    )
    for call in (
        lambda: poison_store.preview_proposal_subset(
            actor_user_id="usr_admin_a",
            proposal_id=poison.proposal_id,
            selected_item_ids=["item_poison"],
            now=now,
        ),
        lambda: poison_store.adopt_proposal_subset(
            actor_user_id="usr_admin_a",
            proposal_id=poison.proposal_id,
            preview_digest=poison_digest,
            selected_item_ids=["item_poison"],
            expected_revisions={SYNTHETIC_SMP_1_ID: 1},
            idempotency_key="idem-poison-tag",
            now=now,
        ),
    ):
        with pytest.raises(CatalogPreviewError) as poison_exc:
            call()
        assert poison_exc.value.status_code == 422
        assert poison_exc.value.code == "catalog_payload_invalid"
    assert poison_store.get_sample_head(SYNTHETIC_SMP_1_ID)["revision_number"] == 1
    assert (
        poison_store.list_adoption_receipts(
            actor_user_id="usr_admin_a", project_id=SYNTHETIC_PROJECT_ID
        )["receipts"]
        == []
    )
    assert "secret-alias" not in json.dumps(
        poison_store.get_sample_head(SYNTHETIC_SMP_1_ID)
    )

    # A sample revision id is not a term, and a pin outside the project is not adopted.
    bad_term_store = build_synthetic_m2_preview_fixture()
    bad_term = _install_experiment_proposal(
        bad_term_store,
        proposal_id="cprop_cccccccccccccccccccccccccccccccc",
        allowed_fields=("batch",),
        before={"batch": {"status": "unknown"}},
        after={
            "batch": {
                "status": "known",
                "term_id": _SMP_REV,
                "term_revision_id": _TREATMENT_REV,
            }
        },
        reasons=("BATCH_TERM_MATCH",),
    )
    with pytest.raises(CatalogPreviewError) as sample_as_term_exc:
        bad_term_store.preview_proposal_subset(
            actor_user_id="usr_admin_a",
            proposal_id=bad_term.proposal_id,
            selected_item_ids=["item_samples"],
            now=now,
        )
    assert sample_as_term_exc.value.code == "catalog_payload_invalid"
    assert (
        bad_term_store.get_experiment_head(SYNTHETIC_EXP_2_ID)["revision_number"] == 1
    )

    extra_store = build_synthetic_m2_preview_fixture()
    extra_pin = deepcopy(
        extra_store.get_experiment_head(SYNTHETIC_EXP_2_ID)["payload"]["samples"][0]
    )
    extra_pin["donor_note"] = "do-not-write"
    extra_proposal = _install_experiment_proposal(
        extra_store,
        proposal_id="cprop_dddddddddddddddddddddddddddddddd",
        allowed_fields=("samples",),
        before={
            "samples": deepcopy(
                extra_store.get_experiment_head(SYNTHETIC_EXP_2_ID)["payload"][
                    "samples"
                ]
            )
        },
        after={"samples": [extra_pin]},
    )
    with pytest.raises(CatalogPreviewError) as extra_exc:
        extra_store.preview_proposal_subset(
            actor_user_id="usr_admin_a",
            proposal_id=extra_proposal.proposal_id,
            selected_item_ids=["item_samples"],
            now=now,
        )
    assert extra_exc.value.code == "catalog_payload_invalid"
    assert "do-not-write" not in json.dumps(
        extra_store.get_experiment_head(SYNTHETIC_EXP_2_ID)
    )

    cross_store = build_synthetic_m2_preview_fixture()
    cross_store._samples[SYNTHETIC_SMP_1_ID].project_id = SYNTHETIC_OTHER_PROJECT_ID
    cross_pin = deepcopy(changed_pin)
    cross_proposal = _install_experiment_proposal(
        cross_store,
        proposal_id="cprop_eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee",
        allowed_fields=("samples",),
        before={"samples": [cross_pin]},
        after={"samples": [cross_pin]},
    )
    with pytest.raises(CatalogPreviewError) as cross_exc:
        cross_store.preview_proposal_subset(
            actor_user_id="usr_admin_a",
            proposal_id=cross_proposal.proposal_id,
            selected_item_ids=["item_samples"],
            now=now,
        )
    assert cross_exc.value.status_code == 409
    assert cross_exc.value.code == "catalog_proposal_invalidated"
    assert cross_store.get_experiment_head(SYNTHETIC_EXP_2_ID)["revision_number"] == 1


@pytest.mark.parametrize(
    "slot",
    [
        "unknown",
        {"status": []},
        {"status": "custom"},
        {"status": "known"},
        {"status": "unknown", "candidates": []},
        {"status": "ambiguous", "values": []},
        {"status": "ambiguous", "candidates": {}},
    ],
)
def test_m2_rejects_malformed_manual_classification_without_consuming_receipt(slot):
    store = build_synthetic_m2_preview_fixture()
    original = store.get_experiment_head(SYNTHETIC_EXP_1_ID)
    kwargs = {
        "actor_user_id": "usr_admin_a",
        "experiment_id": SYNTHETIC_EXP_1_ID,
        "expected_revision_number": 1,
        "idempotency_key": "correctable-slot",
        "now": SYNTHETIC_BASE_TIME + timedelta(minutes=1),
    }
    with pytest.raises(CatalogPreviewError) as failure:
        store.manual_revise_experiment(**kwargs, updates={"batch": slot})
    assert (failure.value.status_code, failure.value.code) == (
        422,
        "catalog_payload_invalid",
    )
    assert store.get_experiment_head(SYNTHETIC_EXP_1_ID) == original

    result = store.manual_revise_experiment(
        **kwargs, updates={"batch": {"status": "not_applicable"}}
    )
    assert result["revision_number"] == 2
    assert store.get_experiment_head(SYNTHETIC_EXP_1_ID)["payload"]["batch"] == {
        "status": "not_applicable"
    }


def test_m2_sample_rename_keeps_material_and_single_candidate_is_redacted():
    store = build_synthetic_m2_preview_fixture()
    now = SYNTHETIC_BASE_TIME + timedelta(minutes=1)
    original = store.get_sample_head(SYNTHETIC_SMP_1_ID)
    species = original["catalog"]["species"]
    pin = {key: species[key] for key in ("term_id", "term_revision_id")}
    renamed = "Renamed sample"
    item = ProposalItem(
        item_id="item_rename",
        target_kind="sample_catalog",
        target_id=SYNTHETIC_SMP_1_ID,
        base_revision_id=original["sample_revision_id"],
        base_revision_number=original["revision_number"],
        allowed_fields=("display_name", "species"),
        before={"display_name": original["display_name"], "species": species},
        after={
            "display_name": renamed,
            "species": {"status": "ambiguous", "candidates": [pin]},
        },
        citation=ProposalCitation(
            rule_version=M2_RULE_VERSION,
            vocabulary_version=M2_VOCABULARY_VERSION,
            reason_codes=("ALIAS_AMBIGUOUS_CANDIDATES",),
            source_refs=(
                SourceReference(
                    kind="sample_revision",
                    ref_id=original["sample_revision_id"],
                    evidence_token="sample-classification",
                ),
            ),
        ),
    )
    proposal = CatalogProposal.create(
        proposal_id="cprop_" + "3" * 32,
        project_id=SYNTHETIC_PROJECT_ID,
        created_by_user_id="usr_member_01",
        source="offline_suggestion",
        items=(item,),
        created_at=SYNTHETIC_BASE_TIME,
    )
    store._proposals[proposal.proposal_id] = proposal
    member = store.get_proposal_view(
        actor_user_id="usr_member_01", proposal_id=proposal.proposal_id, now=now
    )
    assert member["items"][0]["after"]["species"] == {"status": "ambiguous"}

    selection = {
        "actor_user_id": "usr_admin_a",
        "proposal_id": proposal.proposal_id,
        "selected_item_ids": [item.item_id],
        "now": now,
    }
    preview = store.preview_proposal_subset(**selection)
    assert preview["items"][0]["after"]["species"]["candidates"] == [
        pin | {"display_name": "Mus musculus"}
    ]
    assert store.get_sample_head(SYNTHETIC_SMP_1_ID) == original

    store.adopt_proposal_subset(
        **selection,
        preview_digest=preview["preview_digest"],
        expected_revisions=preview["expected_revisions"],
        idempotency_key="rename-sample",
    )
    revised = store.get_sample_head(SYNTHETIC_SMP_1_ID)
    assert revised["display_name"] == renamed
    assert revised["revision_number"] == original["revision_number"] + 1
    assert revised["extra_attributes"] == original["extra_attributes"]
    assert revised["catalog"] == original["catalog"] | {
        "species": {"status": "ambiguous", "candidates": [pin]}
    }
