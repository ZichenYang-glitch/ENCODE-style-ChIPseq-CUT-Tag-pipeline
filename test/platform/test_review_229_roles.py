from copy import deepcopy
from datetime import timedelta

import pytest
from encode_pipeline.platform import catalog_curation_preview as m


@pytest.mark.parametrize("experiment", [m.SYNTHETIC_EXP_2_ID, m.SYNTHETIC_EXP_1_ID])
def test_change_first_role_preserves_every_other_pin(experiment):
    store = m.build_synthetic_m2_preview_fixture()
    head = store.get_experiment_head(experiment)
    pins = head["payload"]["samples"]
    after = deepcopy(pins)
    after[0]["replicate_role"] = {
        "status": "known",
        "term_id": "term_" + "78" * 16,
        "term_revision_id": "trmv_" + "90" * 16,
    }
    item = m.ProposalItem(
        item_id="role",
        target_kind="experiment",
        target_id=experiment,
        base_revision_id=head["experiment_revision_id"],
        base_revision_number=head["revision_number"],
        allowed_fields=("samples",),
        before={"samples": pins},
        after={"samples": after},
        citation=m.ProposalCitation(
            rule_version=m.M2_RULE_VERSION,
            vocabulary_version=m.M2_VOCABULARY_VERSION,
            reason_codes=("REPLICATE_ROLE_TERM_MATCH",),
            source_refs=(
                m.SourceReference(
                    kind="experiment_revision",
                    ref_id=head["experiment_revision_id"],
                    evidence_token="role-review",
                ),
                m.SourceReference(
                    kind="term_revision",
                    ref_id=after[0]["replicate_role"]["term_revision_id"],
                    evidence_token="vocab-control",
                ),
            ),
        ),
    )
    proposal = m.CatalogProposal.create(
        proposal_id="cprop_" + "cf" * 16,
        project_id=m.SYNTHETIC_PROJECT_ID,
        created_by_user_id="usr_member_01",
        source="offline_suggestion",
        items=(item,),
        created_at=m.SYNTHETIC_BASE_TIME,
    )
    store._proposals[proposal.proposal_id] = proposal
    args = dict(
        actor_user_id="usr_admin_a",
        proposal_id=proposal.proposal_id,
        now=m.SYNTHETIC_BASE_TIME + timedelta(minutes=1),
        selected_item_ids=["role"],
    )
    try:
        preview = store.preview_proposal_subset(**args)
    except m.CatalogPreviewError as exc:
        assert store.get_experiment_head(experiment) == head
        assert (
            store.list_adoption_receipts(
                actor_user_id="usr_admin_a", project_id=m.SYNTHETIC_PROJECT_ID
            )["receipts"]
            == []
        )
        pytest.fail(
            f"{len(pins)} pins: {exc.status_code} {exc.code}; revision {head['revision_number']}, zero receipts"
        )
    assert (
        preview["items"][0]["after"]["samples"][0]["replicate_role"]["display_name"]
        == "control"
    )
    store.adopt_proposal_subset(
        **args,
        expected_revisions=preview["expected_revisions"],
        preview_digest=preview["preview_digest"],
        idempotency_key="change-role",
    )
    actual = store.get_experiment_head(experiment)
    assert actual["revision_number"] == head["revision_number"] + 1
    assert actual["payload"]["samples"] == after


@pytest.mark.parametrize(
    "field,collection", [("tag", "values"), ("cell_type", "candidates")]
)
def test_same_term_cannot_be_stored_twice(field, collection):
    from dataclasses import replace

    store = m.build_synthetic_m2_preview_fixture()
    base = store._proposals[m.SYNTHETIC_PROPOSAL_ID]
    item = base.items[2]
    after = deepcopy(item.after)
    pins = after[field][collection]
    pins.append(deepcopy(pins[0]))
    proposal = m.CatalogProposal.create(
        proposal_id="cprop_" + "ce" * 16,
        project_id=base.project_id,
        created_by_user_id=base.created_by_user_id,
        source=base.source,
        items=(replace(item, after=after),),
        created_at=base.created_at,
    )
    store._proposals[proposal.proposal_id] = proposal
    args = dict(
        actor_user_id="usr_admin_a",
        proposal_id=proposal.proposal_id,
        now=m.SYNTHETIC_BASE_TIME + timedelta(minutes=1),
        selected_item_ids=[item.item_id],
    )
    try:
        preview = store.preview_proposal_subset(**args)
    except m.CatalogPreviewError as exc:
        assert (exc.status_code, exc.code) == (422, "catalog_payload_invalid")
        return
    store.adopt_proposal_subset(
        **args,
        expected_revisions=preview["expected_revisions"],
        preview_digest=preview["preview_digest"],
        idempotency_key="duplicate-classification",
    )
    actual = store.get_sample_head(item.target_id)
    occurrences = sum(
        p["term_id"] == pins[0]["term_id"] for p in actual["catalog"][field][collection]
    )
    pytest.fail(
        f"{field}.{collection}: {occurrences} copies stored in sample revision {actual['revision_number']}"
    )
