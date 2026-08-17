from __future__ import annotations

import json
import os
from pathlib import Path

import pytest

from encode_pipeline.deployment.errors import DeploymentError
from encode_pipeline.deployment.backend import _deployment_state_probe
from encode_pipeline.deployment.doctor import FAIL, WARNING
from encode_pipeline.deployment.layout import DeploymentLayout
from encode_pipeline.deployment.operator import OperatorObservation, SERVICE_UNITS
from encode_pipeline.deployment.operator_transaction import (
    OperatorJournalStore,
    OperatorTransaction,
)


IDENTITY = "sha256-" + "a" * 64
PRIOR = "sha256-" + "b" * 64
CANDIDATE = "sha256-" + "c" * 64
BACKUP_SLOT = "sha256-" + "d" * 64
BACKUP_RECEIPT = "sha256-" + "e" * 64
DATABASE = "sha256-" + "f" * 64
TASK = "task-" + "1" * 32


def _store(tmp_path: Path) -> tuple[DeploymentLayout, OperatorJournalStore]:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    (layout.data_root / "operator").mkdir(parents=True, mode=0o700)
    return layout, OperatorJournalStore(
        layout,
        owner_uid=os.getuid(),
        owner_gid=os.getgid(),
    )


def test_terminal_journal_moves_from_single_active_pointer_to_history(
    tmp_path: Path,
) -> None:
    layout, store = _store(tmp_path)

    with store.operation(
        operation="stage",
        task_identity=TASK,
        deployment_identity=IDENTITY,
        component="platform",
        unit=None,
    ) as journal:
        assert store.summary().pending_count == 1
        assert json.loads(layout.operator_transaction_active.read_text())["phase"] == (
            "prepared"
        )
        journal.advance("ingress-pinned")
        journal.complete()

    assert not layout.operator_transaction_active.exists()
    archived = json.loads(
        (layout.operator_transaction_history / f"{TASK}.json").read_text()
    )
    assert archived["phase"] == "complete"
    assert store.summary().pending_count == 0


def test_pre_side_effect_failure_is_durably_aborted_and_does_not_block_retry(
    tmp_path: Path,
) -> None:
    layout, store = _store(tmp_path)

    with pytest.raises(RuntimeError, match="injected"):
        with store.operation(
            operation="activate",
            task_identity=TASK,
            deployment_identity=IDENTITY,
            component="platform",
            unit=None,
        ) as journal:
            journal.advance("runtime-prepared")
            raise RuntimeError("injected")

    archived = json.loads(
        (layout.operator_transaction_history / f"{TASK}.json").read_text()
    )
    assert archived["phase"] == "aborted"
    assert not layout.operator_transaction_active.exists()
    second = "task-" + "2" * 32
    with store.operation(
        operation="stage",
        task_identity=second,
        deployment_identity=IDENTITY,
        component="platform",
        unit=None,
    ) as journal:
        journal.complete()


def test_post_side_effect_failure_leaves_one_global_recovery_record(
    tmp_path: Path,
) -> None:
    layout, store = _store(tmp_path)

    with pytest.raises(RuntimeError, match="injected"):
        with store.operation(
            operation="activate",
            task_identity=TASK,
            deployment_identity=IDENTITY,
            component="platform",
            unit=None,
        ) as journal:
            journal.advance(
                "candidate-selected",
                evidence={
                    "prior_state_identity": PRIOR,
                    "candidate_state_identity": CANDIDATE,
                },
            )
            selected = json.loads(layout.operator_transaction_active.read_text())
            assert selected["prior_state_identity"] == PRIOR
            assert selected["candidate_state_identity"] == CANDIDATE
            journal.advance(
                "service-stopping",
                restart_units=("helixweave-api.service",),
            )
            raise RuntimeError("injected")

    active = json.loads(layout.operator_transaction_active.read_text())
    assert active["phase"] == "recovery-required"
    assert active["prior_state_identity"] == PRIOR
    assert active["candidate_state_identity"] == CANDIDATE
    assert active["restart_units"] == ["helixweave-api.service"]
    assert active["point_of_no_return"] is False
    assert store.summary().recovery_required_count == 1
    with pytest.raises(DeploymentError) as captured:
        with store.operation(
            operation="stage",
            task_identity="task-" + "3" * 32,
            deployment_identity=IDENTITY,
            component="platform",
            unit=None,
        ):
            pass
    assert captured.value.issue.code == "OPERATOR_RECOVERY_REQUIRED"


def test_recovery_fields_and_point_of_no_return_are_fsynced_before_start(
    tmp_path: Path,
) -> None:
    layout, store = _store(tmp_path)

    with pytest.raises(RuntimeError, match="crash-before-start"):
        with store.operation(
            operation="activate",
            task_identity=TASK,
            deployment_identity=IDENTITY,
            component="platform",
            unit=None,
        ) as journal:
            journal.advance(
                "writers-stopped",
                restart_units=(
                    "helixweave-api.service",
                    "helixweave-worker.service",
                ),
                schema_before_heads=("schema-12",),
                target_schema_heads=("schema-13",),
                evidence={
                    "prior_state_identity": PRIOR,
                    "candidate_state_identity": CANDIDATE,
                    "source_database_identity": DATABASE,
                },
            )
            journal.advance(
                "backup-committed",
                evidence={
                    "backup_slot_identity": BACKUP_SLOT,
                    "backup_receipt_identity": BACKUP_RECEIPT,
                },
            )
            journal.advance(
                "service-starting",
                point_of_no_return=True,
            )
            durable = json.loads(layout.operator_transaction_active.read_text())
            assert durable["point_of_no_return"] is True
            assert durable["restart_units"] == [
                "helixweave-api.service",
                "helixweave-worker.service",
            ]
            assert durable["backup_slot_identity"] == BACKUP_SLOT
            assert durable["backup_receipt_identity"] == BACKUP_RECEIPT
            assert durable["source_database_identity"] == DATABASE
            assert durable["schema_before_heads"] == ["schema-12"]
            assert durable["target_schema_heads"] == ["schema-13"]
            raise RuntimeError("crash-before-start")

    assert json.loads(layout.operator_transaction_active.read_text())["phase"] == (
        "recovery-required"
    )


def test_restart_set_is_durable_before_point_of_no_return_and_immutable(
    tmp_path: Path,
) -> None:
    layout, store = _store(tmp_path)
    units = ("helixweave-api.service", "helixweave-worker.service")

    with pytest.raises(RuntimeError, match="crash-before-ponr"):
        with store.operation(
            operation="activate",
            task_identity=TASK,
            deployment_identity=IDENTITY,
            component="platform",
            unit=None,
        ) as journal:
            journal.advance("service-stopping", restart_units=units)
            durable = json.loads(layout.operator_transaction_active.read_text())
            assert durable["restart_units"] == list(units)
            assert durable["point_of_no_return"] is False
            with pytest.raises(DeploymentError) as changed:
                journal.advance(
                    "writers-stopped",
                    restart_units=("helixweave-api.service",),
                )
            assert changed.value.issue.code == "OPERATOR_JOURNAL_INVALID"
            raise RuntimeError("crash-before-ponr")

    recovered = json.loads(layout.operator_transaction_active.read_text())
    assert recovered["phase"] == "recovery-required"
    assert recovered["restart_units"] == list(units)
    assert recovered["point_of_no_return"] is False


def test_hard_crash_phase_is_normalized_before_recovery_controller(
    tmp_path: Path,
) -> None:
    layout, seed_store = _store(tmp_path)
    seed_store._directories(create=True)
    record = OperatorTransaction.create(
        request_identity=IDENTITY,
        operation="activate",
        task_identity=TASK,
        deployment_identity=IDENTITY,
        component="platform",
        unit=None,
        phase="service-stopping",
        restart_units=("helixweave-api.service",),
        prior_active={
            "platform": PRIOR,
            "encode-runtime": None,
            "bulk-rnaseq-runtime": None,
        },
        candidate_active={
            "platform": CANDIDATE,
            "encode-runtime": None,
            "bulk-rnaseq-runtime": None,
        },
        prior_state_identity=PRIOR,
        candidate_state_identity=CANDIDATE,
    )
    seed_store._write_new(layout.operator_transaction_active, record)
    observed: list[str] = []

    class RecoveryController:
        def recover(self, active: OperatorTransaction) -> OperatorTransaction:
            observed.append(active.phase)
            raise RuntimeError("recovery-injected")

    store = OperatorJournalStore(
        layout,
        owner_uid=os.getuid(),
        owner_gid=os.getgid(),
        recovery_controller=RecoveryController(),
    )
    with pytest.raises(RuntimeError, match="recovery-injected"):
        with store.operation(
            operation="stage",
            task_identity="task-" + "2" * 32,
            deployment_identity=IDENTITY,
            component="platform",
            unit=None,
        ):
            pass

    assert observed == ["recovery-required"]
    active = json.loads(layout.operator_transaction_active.read_text())
    assert active["phase"] == "recovery-required"
    assert active["failure_phase"] == "service-stopping"


def test_journal_phase_must_advance_monotonically(tmp_path: Path) -> None:
    _layout, store = _store(tmp_path)
    with store.operation(
        operation="stage",
        task_identity=TASK,
        deployment_identity=IDENTITY,
        component="platform",
        unit=None,
    ) as journal:
        journal.advance("payload-staged")
        for phase in ("payload-staged", "ingress-pinned"):
            with pytest.raises(DeploymentError) as captured:
                journal.advance(phase)
            assert captured.value.issue.code == "OPERATOR_JOURNAL_INVALID"
        journal.complete()


def test_journal_rejects_noncanonical_heads_and_history_replay(tmp_path: Path) -> None:
    _layout, store = _store(tmp_path)
    with store.operation(
        operation="stage",
        task_identity=TASK,
        deployment_identity=IDENTITY,
        component="platform",
        unit=None,
    ) as journal:
        with pytest.raises(DeploymentError) as captured:
            journal.advance(
                "payload-staged",
                target_schema_heads=("schema-2", "schema-1"),
            )
        assert captured.value.issue.code == "OPERATOR_JOURNAL_INVALID"
        journal.complete()

    with pytest.raises(DeploymentError) as replay:
        with store.operation(
            operation="stage",
            task_identity=TASK,
            deployment_identity=IDENTITY,
            component="platform",
            unit=None,
        ):
            pass
    assert replay.value.issue.code == "OPERATOR_TASK_REPLAY"


@pytest.mark.parametrize(
    ("pending", "recovery", "state", "reason"),
    (
        (1, 0, WARNING, "DEPLOYMENT_INTERRUPTED"),
        (0, 1, FAIL, "DEPLOYMENT_RECOVERY_REQUIRED"),
    ),
)
def test_doctor_maps_path_free_operator_journal_counts(
    pending: int,
    recovery: int,
    state: str,
    reason: str,
) -> None:
    observation = OperatorObservation.create(
        state_identity=IDENTITY,
        active={
            component: None
            for component in (
                "platform",
                "encode-runtime",
                "bulk-rnaseq-runtime",
            )
        },
        database_schema_identity=DATABASE,
        database_schema_heads=("schema-12",),
        services={unit: None for unit in SERVICE_UNITS},
        operator_pending_count=pending,
        operator_recovery_required_count=recovery,
    )

    result = _deployment_state_probe(object(), observation)()

    assert result.state == state
    assert result.reason_code == reason
    assert result.evidence_identity is None
