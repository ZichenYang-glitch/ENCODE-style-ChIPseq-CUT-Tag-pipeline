"""Single-active root-owned journal for privileged deployment side effects."""

from __future__ import annotations

from collections.abc import Iterator, Mapping
from contextlib import contextmanager
from dataclasses import dataclass
import fcntl
import json
import os
from pathlib import Path
import re
import stat
from typing import Any, Protocol
import uuid

from encode_pipeline.deployment.canonical import (
    canonical_identity,
    canonical_json_bytes,
    without_key,
)
from encode_pipeline.deployment.errors import DeploymentError, fail
from encode_pipeline.deployment.layout import DeploymentLayout


JOURNAL_SCHEMA = "helixweave-operator-transaction-v3"
JOURNAL_IDENTITY_SCHEME = "helixweave-operator-transaction-identity-v3"
PHASES = (
    "prepared",
    "ingress-pinned",
    "payload-staged",
    "candidate-selected",
    "runtime-prepare-started",
    "runtime-prepared",
    "service-stopping",
    "writers-stopped",
    "bulk-runtime-prepare-started",
    "bulk-runtime-prepared",
    "witness-ready",
    "backup-committed",
    "migration-started",
    "migration-verified",
    "database-published",
    "state-committed",
    "dependency-starting",
    "service-starting",
    "service-observed",
    "complete",
    "aborted",
    "recovery-required",
)
TERMINAL_PHASES = frozenset({"complete", "aborted"})
EVIDENCE_KEYS = frozenset(
    {
        "service_identity",
        "witness_identity",
        "action_receipt_identity",
        "database_prepare_receipt_identity",
        "configuration_identity",
        "runtime_prepare_receipt_identity",
        "runtime_tree_identity",
        "bulk_runtime_prepare_receipt_identity",
        "bulk_runtime_image_set_identity",
        "bulk_runtime_service_identity",
    }
)
RECOVERY_EVIDENCE_KEYS = frozenset(
    {
        "prior_state_identity",
        "candidate_state_identity",
        "backup_slot_identity",
        "backup_receipt_identity",
        "source_database_identity",
        "schema_before_identity",
        "schema_after_identity",
    }
)
_OPERATIONS = frozenset(
    {"stage", "activate", "rollback", "start", "status", "stop", "cleanup", "uninstall"}
)
_COMPONENTS = frozenset({"platform", "encode-runtime", "bulk-rnaseq-runtime"})
_UNITS = (
    "helixweave-redis.service",
    "helixweave-docker-rootless.service",
    "helixweave-api.service",
    "helixweave-worker.service",
)
_UNIT_SET = frozenset(_UNITS)
_WRITER_UNITS = frozenset({"helixweave-api.service", "helixweave-worker.service"})
_DATABASE_MODES = frozenset({"none", "fresh-candidate", "existing-live"})
_IDENTITY_PREFIX = "sha256-"
_TASK_PREFIX = "task-"
_SCHEMA_HEAD = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,127}$")
_MAX_JOURNAL_BYTES = 256 * 1024
_SAFE_ABORT_PHASES = frozenset(
    {
        "prepared",
        "ingress-pinned",
        "payload-staged",
        "candidate-selected",
        "runtime-prepare-started",
        "runtime-prepared",
    }
)
_STAGE_SAFE_ABORT_PHASES = frozenset({"prepared", "ingress-pinned", "payload-staged"})


def _identity(value: object) -> str:
    if (
        not isinstance(value, str)
        or len(value) != 71
        or not value.startswith(_IDENTITY_PREFIX)
        or any(character not in "0123456789abcdef" for character in value[7:])
    ):
        raise fail("OPERATOR_JOURNAL_INVALID", "Operator journal is invalid.")
    return value


def _optional_identity(value: object) -> str | None:
    return None if value is None else _identity(value)


def _task(value: object) -> str:
    if (
        not isinstance(value, str)
        or len(value) != 37
        or not value.startswith(_TASK_PREFIX)
        or any(character not in "0123456789abcdef" for character in value[5:])
    ):
        raise fail("OPERATOR_JOURNAL_INVALID", "Operator journal is invalid.")
    return value


def _restart_units(value: object) -> tuple[str, ...]:
    if not isinstance(value, list):
        raise fail("OPERATOR_JOURNAL_INVALID", "Operator journal is invalid.")
    units = tuple(value)
    if (
        any(not isinstance(unit, str) or unit not in _UNIT_SET for unit in units)
        or len(units) != len(set(units))
        or units != tuple(unit for unit in _UNITS if unit in set(units))
    ):
        raise fail("OPERATOR_JOURNAL_INVALID", "Operator journal is invalid.")
    return units


def _schema_heads(value: object, *, allow_empty: bool) -> tuple[str, ...]:
    if not isinstance(value, list):
        raise fail("OPERATOR_JOURNAL_INVALID", "Operator journal is invalid.")
    heads = tuple(value)
    if (
        (not allow_empty and not heads)
        or len(heads) > 16
        or tuple(sorted(set(heads))) != heads
        or any(
            not isinstance(head, str) or _SCHEMA_HEAD.fullmatch(head) is None
            for head in heads
        )
    ):
        raise fail("OPERATOR_JOURNAL_INVALID", "Operator journal is invalid.")
    return heads


def _component_identities(value: object) -> dict[str, str | None] | None:
    if value is None:
        return None
    if not isinstance(value, dict) or set(value) != _COMPONENTS:
        raise fail("OPERATOR_JOURNAL_INVALID", "Operator journal is invalid.")
    return {
        component: (None if value[component] is None else _identity(value[component]))
        for component in sorted(_COMPONENTS)
    }


@dataclass(frozen=True)
class OperatorTransaction:
    identity: str
    request_identity: str
    operation: str
    task_identity: str
    deployment_identity: str
    component: str | None
    unit: str | None
    phase: str
    failure_phase: str | None
    write_fence: bool
    point_of_no_return: bool
    restart_units: tuple[str, ...]
    prior_running_units: tuple[str, ...]
    prior_active: Mapping[str, str | None] | None
    candidate_active: Mapping[str, str | None] | None
    database_mode: str
    prior_state_identity: str | None
    candidate_state_identity: str | None
    backup_slot_identity: str | None
    backup_receipt_identity: str | None
    source_database_identity: str | None
    schema_before_identity: str | None
    schema_after_identity: str | None
    schema_before_heads: tuple[str, ...]
    target_schema_heads: tuple[str, ...]
    evidence: Mapping[str, str]

    @classmethod
    def create(
        cls,
        *,
        request_identity: str,
        operation: str,
        task_identity: str,
        deployment_identity: str,
        component: str | None,
        unit: str | None,
        phase: str,
        failure_phase: str | None = None,
        write_fence: bool = False,
        point_of_no_return: bool = False,
        restart_units: tuple[str, ...] = (),
        prior_running_units: tuple[str, ...] = (),
        prior_active: Mapping[str, str | None] | None = None,
        candidate_active: Mapping[str, str | None] | None = None,
        database_mode: str = "none",
        prior_state_identity: str | None = None,
        candidate_state_identity: str | None = None,
        backup_slot_identity: str | None = None,
        backup_receipt_identity: str | None = None,
        source_database_identity: str | None = None,
        schema_before_identity: str | None = None,
        schema_after_identity: str | None = None,
        schema_before_heads: tuple[str, ...] = (),
        target_schema_heads: tuple[str, ...] = (),
        evidence: Mapping[str, str] | None = None,
    ) -> "OperatorTransaction":
        value: dict[str, object] = {
            "schema_version": JOURNAL_SCHEMA,
            "request_identity": request_identity,
            "operation": operation,
            "task_identity": task_identity,
            "deployment_identity": deployment_identity,
            "component": component,
            "unit": unit,
            "phase": phase,
            "failure_phase": failure_phase,
            "write_fence": write_fence,
            "point_of_no_return": point_of_no_return,
            "restart_units": list(restart_units),
            "prior_running_units": list(prior_running_units),
            "prior_active": None if prior_active is None else dict(prior_active),
            "candidate_active": (
                None if candidate_active is None else dict(candidate_active)
            ),
            "database_mode": database_mode,
            "prior_state_identity": prior_state_identity,
            "candidate_state_identity": candidate_state_identity,
            "backup_slot_identity": backup_slot_identity,
            "backup_receipt_identity": backup_receipt_identity,
            "source_database_identity": source_database_identity,
            "schema_before_identity": schema_before_identity,
            "schema_after_identity": schema_after_identity,
            "schema_before_heads": list(schema_before_heads),
            "target_schema_heads": list(target_schema_heads),
            "evidence": dict(evidence or {}),
        }
        value["identity"] = canonical_identity(value, scheme=JOURNAL_IDENTITY_SCHEME)
        return cls.from_dict(value)

    @classmethod
    def from_dict(cls, raw: object) -> "OperatorTransaction":
        keys = {
            "schema_version",
            "identity",
            "request_identity",
            "operation",
            "task_identity",
            "deployment_identity",
            "component",
            "unit",
            "phase",
            "failure_phase",
            "write_fence",
            "point_of_no_return",
            "restart_units",
            "prior_running_units",
            "prior_active",
            "candidate_active",
            "database_mode",
            "prior_state_identity",
            "candidate_state_identity",
            "backup_slot_identity",
            "backup_receipt_identity",
            "source_database_identity",
            "schema_before_identity",
            "schema_after_identity",
            "schema_before_heads",
            "target_schema_heads",
            "evidence",
        }
        if not isinstance(raw, dict) or set(raw) != keys:
            raise fail("OPERATOR_JOURNAL_INVALID", "Operator journal is invalid.")
        operation = raw["operation"]
        component = raw["component"]
        unit = raw["unit"]
        phase = raw["phase"]
        evidence = raw["evidence"]
        if (
            raw["schema_version"] != JOURNAL_SCHEMA
            or operation not in _OPERATIONS
            or (component is not None and component not in _COMPONENTS)
            or (unit is not None and unit not in _UNIT_SET)
            or phase not in PHASES
            or (raw["failure_phase"] is not None and raw["failure_phase"] not in PHASES)
            or (phase == "recovery-required" and raw["failure_phase"] is None)
            or (
                raw["failure_phase"] is not None
                and phase not in TERMINAL_PHASES | {"recovery-required"}
            )
            or not isinstance(raw["write_fence"], bool)
            or not isinstance(raw["point_of_no_return"], bool)
            or not isinstance(evidence, dict)
            or not set(evidence).issubset(EVIDENCE_KEYS)
            or any(not isinstance(key, str) for key in evidence)
        ):
            raise fail("OPERATOR_JOURNAL_INVALID", "Operator journal is invalid.")
        units = _restart_units(raw["restart_units"])
        prior_running = _restart_units(raw["prior_running_units"])
        prior_active = _component_identities(raw["prior_active"])
        candidate_active = _component_identities(raw["candidate_active"])
        if (
            (raw["point_of_no_return"] and not _WRITER_UNITS.intersection(units))
            or not set(prior_running).issubset(units)
            or raw["database_mode"] not in _DATABASE_MODES
            or (prior_active is None) != (candidate_active is None)
            or (operation in {"stage", "activate", "rollback"})
            != (component is not None)
            or (operation in {"start", "status", "stop", "cleanup"})
            != (unit is not None)
            or (
                operation == "uninstall" and (component is not None or unit is not None)
            )
        ):
            raise fail("OPERATOR_JOURNAL_INVALID", "Operator journal is invalid.")
        parsed_evidence = {key: _identity(value) for key, value in evidence.items()}
        identity = _identity(raw["identity"])
        if identity != canonical_identity(
            without_key(raw, "identity"), scheme=JOURNAL_IDENTITY_SCHEME
        ):
            raise fail("OPERATOR_JOURNAL_INVALID", "Operator journal is invalid.")
        return cls(
            identity=identity,
            request_identity=_identity(raw["request_identity"]),
            operation=operation,
            task_identity=_task(raw["task_identity"]),
            deployment_identity=_identity(raw["deployment_identity"]),
            component=component,
            unit=unit,
            phase=phase,
            failure_phase=raw["failure_phase"],
            write_fence=raw["write_fence"],
            point_of_no_return=raw["point_of_no_return"],
            restart_units=units,
            prior_running_units=prior_running,
            prior_active=prior_active,
            candidate_active=candidate_active,
            database_mode=raw["database_mode"],
            prior_state_identity=_optional_identity(raw["prior_state_identity"]),
            candidate_state_identity=_optional_identity(
                raw["candidate_state_identity"]
            ),
            backup_slot_identity=_optional_identity(raw["backup_slot_identity"]),
            backup_receipt_identity=_optional_identity(raw["backup_receipt_identity"]),
            source_database_identity=_optional_identity(
                raw["source_database_identity"]
            ),
            schema_before_identity=_optional_identity(raw["schema_before_identity"]),
            schema_after_identity=_optional_identity(raw["schema_after_identity"]),
            schema_before_heads=_schema_heads(
                raw["schema_before_heads"], allow_empty=True
            ),
            target_schema_heads=_schema_heads(
                raw["target_schema_heads"], allow_empty=True
            ),
            evidence=parsed_evidence,
        )

    def to_dict(self) -> dict[str, object]:
        return {
            "schema_version": JOURNAL_SCHEMA,
            "identity": self.identity,
            "request_identity": self.request_identity,
            "operation": self.operation,
            "task_identity": self.task_identity,
            "deployment_identity": self.deployment_identity,
            "component": self.component,
            "unit": self.unit,
            "phase": self.phase,
            "failure_phase": self.failure_phase,
            "write_fence": self.write_fence,
            "point_of_no_return": self.point_of_no_return,
            "restart_units": list(self.restart_units),
            "prior_running_units": list(self.prior_running_units),
            "prior_active": (
                None if self.prior_active is None else dict(self.prior_active)
            ),
            "candidate_active": (
                None if self.candidate_active is None else dict(self.candidate_active)
            ),
            "database_mode": self.database_mode,
            "prior_state_identity": self.prior_state_identity,
            "candidate_state_identity": self.candidate_state_identity,
            "backup_slot_identity": self.backup_slot_identity,
            "backup_receipt_identity": self.backup_receipt_identity,
            "source_database_identity": self.source_database_identity,
            "schema_before_identity": self.schema_before_identity,
            "schema_after_identity": self.schema_after_identity,
            "schema_before_heads": list(self.schema_before_heads),
            "target_schema_heads": list(self.target_schema_heads),
            "evidence": dict(self.evidence),
        }


@dataclass(frozen=True)
class OperatorJournalSummary:
    pending_count: int
    recovery_required_count: int

    def __post_init__(self) -> None:
        if (
            isinstance(self.pending_count, bool)
            or not isinstance(self.pending_count, int)
            or isinstance(self.recovery_required_count, bool)
            or not isinstance(self.recovery_required_count, int)
            or self.pending_count not in {0, 1}
            or self.recovery_required_count not in {0, 1}
            or self.pending_count + self.recovery_required_count > 1
        ):
            raise fail("OPERATOR_JOURNAL_INVALID", "Operator journal is invalid.")


class OperatorRecoveryController(Protocol):
    """Idempotent recovery authority for post-side-effect transactions."""

    def recover(self, record: OperatorTransaction) -> OperatorTransaction: ...


class OperatorJournalHandle:
    def __init__(self, store: "OperatorJournalStore", record: OperatorTransaction):
        self._store = store
        self.record = record

    @property
    def terminal(self) -> bool:
        return self.record.phase in TERMINAL_PHASES

    def advance(
        self,
        phase: str,
        *,
        write_fence: bool | None = None,
        point_of_no_return: bool | None = None,
        restart_units: tuple[str, ...] | None = None,
        prior_running_units: tuple[str, ...] | None = None,
        prior_active: Mapping[str, str | None] | None = None,
        candidate_active: Mapping[str, str | None] | None = None,
        database_mode: str | None = None,
        schema_before_heads: tuple[str, ...] | None = None,
        target_schema_heads: tuple[str, ...] | None = None,
        evidence: Mapping[str, str] | None = None,
    ) -> OperatorTransaction:
        if phase == "complete":
            return self.complete()
        if phase == "aborted":
            return self.abort()
        if phase not in PHASES or self.record.phase in TERMINAL_PHASES | {
            "recovery-required"
        }:
            raise fail(
                "OPERATOR_RECOVERY_REQUIRED",
                "Operator transaction requires recovery.",
                recoverable=True,
            )
        if PHASES.index(phase) <= PHASES.index(self.record.phase):
            raise fail("OPERATOR_JOURNAL_INVALID", "Operator journal is invalid.")
        next_evidence = dict(self.record.evidence)
        recovery = {key: getattr(self.record, key) for key in RECOVERY_EVIDENCE_KEYS}
        for key, value in (evidence or {}).items():
            if key in RECOVERY_EVIDENCE_KEYS:
                parsed = _identity(value)
                if recovery[key] is not None and recovery[key] != parsed:
                    raise fail(
                        "OPERATOR_JOURNAL_INVALID", "Operator journal is invalid."
                    )
                recovery[key] = parsed
            else:
                if key not in EVIDENCE_KEYS:
                    raise fail(
                        "OPERATOR_JOURNAL_INVALID", "Operator journal is invalid."
                    )
                next_evidence[key] = _identity(value)
        selected_write_fence = (
            self.record.write_fence if write_fence is None else write_fence
        )
        selected_point = (
            self.record.point_of_no_return
            if point_of_no_return is None
            else point_of_no_return
        )
        selected_units = (
            self.record.restart_units if restart_units is None else restart_units
        )
        selected_prior_running = (
            self.record.prior_running_units
            if prior_running_units is None
            else prior_running_units
        )
        selected_prior_active = (
            self.record.prior_active if prior_active is None else prior_active
        )
        selected_candidate_active = (
            self.record.candidate_active
            if candidate_active is None
            else candidate_active
        )
        selected_database_mode = (
            self.record.database_mode if database_mode is None else database_mode
        )
        selected_before_heads = (
            self.record.schema_before_heads
            if schema_before_heads is None
            else schema_before_heads
        )
        selected_target_heads = (
            self.record.target_schema_heads
            if target_schema_heads is None
            else target_schema_heads
        )
        if (
            not isinstance(selected_write_fence, bool)
            or not isinstance(selected_point, bool)
            or (self.record.write_fence and not selected_write_fence)
            or (self.record.point_of_no_return and not selected_point)
            or (selected_point and not _WRITER_UNITS.intersection(selected_units))
            or (
                self.record.restart_units
                and selected_units != self.record.restart_units
            )
            or (
                self.record.prior_running_units
                and selected_prior_running != self.record.prior_running_units
            )
            or not set(selected_prior_running).issubset(selected_units)
            or (
                self.record.prior_active is not None
                and selected_prior_active != self.record.prior_active
            )
            or (
                self.record.candidate_active is not None
                and selected_candidate_active != self.record.candidate_active
            )
            or (selected_prior_active is None) != (selected_candidate_active is None)
            or (
                self.record.database_mode != "none"
                and selected_database_mode != self.record.database_mode
            )
            or selected_database_mode not in _DATABASE_MODES
            or (
                self.record.schema_before_heads
                and selected_before_heads != self.record.schema_before_heads
            )
            or (
                self.record.target_schema_heads
                and selected_target_heads != self.record.target_schema_heads
            )
        ):
            raise fail("OPERATOR_JOURNAL_INVALID", "Operator journal is invalid.")
        record = OperatorTransaction.create(
            request_identity=self.record.request_identity,
            operation=self.record.operation,
            task_identity=self.record.task_identity,
            deployment_identity=self.record.deployment_identity,
            component=self.record.component,
            unit=self.record.unit,
            phase=phase,
            failure_phase=self.record.failure_phase,
            write_fence=selected_write_fence,
            point_of_no_return=selected_point,
            restart_units=selected_units,
            prior_running_units=selected_prior_running,
            prior_active=selected_prior_active,
            candidate_active=selected_candidate_active,
            database_mode=selected_database_mode,
            schema_before_heads=selected_before_heads,
            target_schema_heads=selected_target_heads,
            evidence=next_evidence,
            **recovery,
        )
        self._store._replace_active(record, expected_identity=self.record.identity)
        self.record = record
        return record

    def complete(self) -> OperatorTransaction:
        if self.record.phase in TERMINAL_PHASES | {"recovery-required"}:
            raise fail(
                "OPERATOR_RECOVERY_REQUIRED",
                "Operator transaction requires recovery.",
                recoverable=True,
            )
        record = self._copy_with_phase("complete")
        self._store._finish(record, expected_identity=self.record.identity)
        self.record = record
        return record

    def abort(self) -> OperatorTransaction:
        if not self._store.safe_to_abort(self.record):
            raise fail(
                "OPERATOR_RECOVERY_REQUIRED",
                "Operator transaction requires recovery.",
                recoverable=True,
            )
        record = self._copy_with_phase("aborted")
        self._store._finish(record, expected_identity=self.record.identity)
        self.record = record
        return record

    def fail(self) -> OperatorTransaction:
        if self.record.phase in TERMINAL_PHASES | {"recovery-required"}:
            return self.record
        if self._store.safe_to_abort(self.record):
            return self.abort()
        record = self._copy_with_phase("recovery-required")
        self._store._replace_active(record, expected_identity=self.record.identity)
        self.record = record
        return record

    def _copy_with_phase(self, phase: str) -> OperatorTransaction:
        return OperatorTransaction.create(
            request_identity=self.record.request_identity,
            operation=self.record.operation,
            task_identity=self.record.task_identity,
            deployment_identity=self.record.deployment_identity,
            component=self.record.component,
            unit=self.record.unit,
            phase=phase,
            failure_phase=(
                self.record.phase
                if phase == "recovery-required"
                else self.record.failure_phase
            ),
            write_fence=self.record.write_fence,
            point_of_no_return=self.record.point_of_no_return,
            restart_units=self.record.restart_units,
            prior_running_units=self.record.prior_running_units,
            prior_active=self.record.prior_active,
            candidate_active=self.record.candidate_active,
            database_mode=self.record.database_mode,
            prior_state_identity=self.record.prior_state_identity,
            candidate_state_identity=self.record.candidate_state_identity,
            backup_slot_identity=self.record.backup_slot_identity,
            backup_receipt_identity=self.record.backup_receipt_identity,
            source_database_identity=self.record.source_database_identity,
            schema_before_identity=self.record.schema_before_identity,
            schema_after_identity=self.record.schema_after_identity,
            schema_before_heads=self.record.schema_before_heads,
            target_schema_heads=self.record.target_schema_heads,
            evidence=self.record.evidence,
        )


class OperatorJournalStore:
    """One global lock and one authoritative active journal; history is direct-read."""

    def __init__(
        self,
        layout: DeploymentLayout,
        *,
        owner_uid: int = 0,
        owner_gid: int = 0,
        recovery_controller: OperatorRecoveryController | None = None,
    ) -> None:
        self.layout = layout
        self.owner_uid = owner_uid
        self.owner_gid = owner_gid
        self.recovery_controller = recovery_controller

    @contextmanager
    def operation(
        self,
        *,
        operation: str,
        task_identity: str,
        deployment_identity: str,
        component: str | None,
        unit: str | None,
    ) -> Iterator[OperatorJournalHandle]:
        request = {
            "operation": operation,
            "task_identity": task_identity,
            "deployment_identity": deployment_identity,
            "component": component,
            "unit": unit,
        }
        request_identity = canonical_identity(
            request, scheme="helixweave-operator-request-identity-v1"
        )
        directory, history = self._directories(create=True)
        lock = self._open_lock(self.layout.operator_transaction_lock, create=True)
        try:
            try:
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError:
                raise fail(
                    "OPERATOR_BUSY",
                    "Another operator transaction is in progress.",
                    recoverable=True,
                ) from None
            self._reconcile_active(directory, history)
            history_path = history / f"{_task(task_identity)}.json"
            if history_path.exists() or history_path.is_symlink():
                self._read(history_path)
                raise fail(
                    "OPERATOR_TASK_REPLAY",
                    "Operator task identity was already used.",
                )
            record = OperatorTransaction.create(
                request_identity=request_identity,
                operation=operation,
                task_identity=task_identity,
                deployment_identity=deployment_identity,
                component=component,
                unit=unit,
                phase="prepared",
            )
            self._write_new(self.layout.operator_transaction_active, record)
            handle = OperatorJournalHandle(self, record)
            try:
                yield handle
            except BaseException:
                failed = handle.fail()
                if (
                    failed.phase == "recovery-required"
                    and self.recovery_controller is not None
                ):
                    try:
                        self._recover_handle(handle)
                    except BaseException:
                        raise fail(
                            "OPERATOR_RECOVERY_REQUIRED",
                            "Operator transaction requires recovery.",
                            recoverable=True,
                        ) from None
                raise
            if not handle.terminal:
                handle.fail()
                raise fail(
                    "OPERATOR_JOURNAL_INCOMPLETE",
                    "Operator transaction did not reach a terminal state.",
                    recoverable=True,
                )
        finally:
            try:
                fcntl.flock(lock, fcntl.LOCK_UN)
            finally:
                os.close(lock)

    def summary(self) -> OperatorJournalSummary:
        try:
            directory, _history = self._directories(create=False)
        except FileNotFoundError:
            return OperatorJournalSummary(0, 0)
        lock = self._open_lock(self.layout.operator_transaction_lock, create=False)
        try:
            try:
                fcntl.flock(lock, fcntl.LOCK_SH | fcntl.LOCK_NB)
            except BlockingIOError:
                return OperatorJournalSummary(1, 0)
            active = self._read_active(directory)
            if active is None:
                return OperatorJournalSummary(0, 0)
            if active.phase == "recovery-required":
                return OperatorJournalSummary(0, 1)
            return OperatorJournalSummary(1, 0)
        finally:
            try:
                fcntl.flock(lock, fcntl.LOCK_UN)
            finally:
                os.close(lock)

    @staticmethod
    def safe_to_abort(record: OperatorTransaction) -> bool:
        return (
            not record.write_fence
            and not record.point_of_no_return
            and (
                (
                    record.operation == "stage"
                    and record.phase in _STAGE_SAFE_ABORT_PHASES
                )
                or (record.operation != "stage" and record.phase in _SAFE_ABORT_PHASES)
            )
        )

    def _reconcile_active(self, directory: Path, history: Path) -> None:
        active = self._read_active(directory)
        if active is None:
            return
        if active.phase in TERMINAL_PHASES:
            self._archive_and_clear(active, directory, history)
            return
        if self.safe_to_abort(active):
            aborted = OperatorJournalHandle(self, active)._copy_with_phase("aborted")
            self._finish(aborted, expected_identity=active.identity)
            return
        if active.phase != "recovery-required":
            required = OperatorJournalHandle(self, active)._copy_with_phase(
                "recovery-required"
            )
            self._replace_active(required, expected_identity=active.identity)
            active = required
        if self.recovery_controller is not None:
            handle = OperatorJournalHandle(self, active)
            self._recover_handle(handle)
            if handle.terminal:
                return
        raise fail(
            "OPERATOR_RECOVERY_REQUIRED",
            "Operator transaction requires recovery.",
            recoverable=True,
        )

    def _recover_handle(self, handle: OperatorJournalHandle) -> None:
        if self.recovery_controller is None:
            raise fail(
                "OPERATOR_RECOVERY_REQUIRED",
                "Operator transaction requires recovery.",
                recoverable=True,
            )
        prior = handle.record
        recovered = self.recovery_controller.recover(prior)
        self._validate_recovery_result(prior, recovered)
        if recovered.phase in TERMINAL_PHASES:
            self._finish(recovered, expected_identity=prior.identity)
        else:
            self._replace_active(recovered, expected_identity=prior.identity)
        handle.record = recovered

    @staticmethod
    def _validate_recovery_result(
        prior: OperatorTransaction,
        recovered: OperatorTransaction,
    ) -> None:
        immutable = (
            "request_identity",
            "operation",
            "task_identity",
            "deployment_identity",
            "component",
            "unit",
            "failure_phase",
            "restart_units",
            "prior_running_units",
            "prior_active",
            "candidate_active",
            "database_mode",
            "prior_state_identity",
            "candidate_state_identity",
            "backup_slot_identity",
            "backup_receipt_identity",
            "source_database_identity",
            "schema_before_heads",
            "target_schema_heads",
        )
        if (
            not isinstance(recovered, OperatorTransaction)
            or any(getattr(recovered, key) != getattr(prior, key) for key in immutable)
            or recovered.phase not in TERMINAL_PHASES | {"recovery-required"}
            or (prior.point_of_no_return and not recovered.point_of_no_return)
            or (prior.write_fence and not recovered.write_fence)
        ):
            raise fail("OPERATOR_JOURNAL_INVALID", "Operator journal is invalid.")

    def _directories(self, *, create: bool) -> tuple[Path, Path]:
        directory = self.layout.operator_transactions
        history = self.layout.operator_transaction_history
        if not directory.parent.exists() and not directory.parent.is_symlink():
            if not create:
                raise FileNotFoundError
            raise fail(
                "OPERATOR_JOURNAL_UNAVAILABLE", "Operator journal is unavailable."
            )
        self._require_parent_directory(directory.parent)
        if create:
            self._ensure_directory(directory)
            self._ensure_directory(history)
        else:
            if not directory.exists() and not directory.is_symlink():
                raise FileNotFoundError
            self._require_directory(directory)
            self._require_directory(history)
        return directory, history

    def _require_parent_directory(self, path: Path) -> None:
        try:
            observed = path.lstat()
        except OSError:
            raise fail(
                "OPERATOR_JOURNAL_UNAVAILABLE", "Operator journal is unavailable."
            ) from None
        if (
            not stat.S_ISDIR(observed.st_mode)
            or stat.S_ISLNK(observed.st_mode)
            or observed.st_uid != self.owner_uid
            or stat.S_IMODE(observed.st_mode) not in {0o700, 0o710}
        ):
            raise fail("OPERATOR_JOURNAL_UNTRUSTED", "Operator journal is not trusted.")

    def _ensure_directory(self, path: Path) -> None:
        try:
            observed = path.lstat()
        except FileNotFoundError:
            try:
                path.mkdir(mode=0o700)
                os.chown(path, self.owner_uid, self.owner_gid)
                _fsync_directory(path.parent)
                observed = path.lstat()
            except OSError:
                raise fail(
                    "OPERATOR_JOURNAL_UNAVAILABLE",
                    "Operator journal is unavailable.",
                ) from None
        except OSError:
            raise fail(
                "OPERATOR_JOURNAL_UNAVAILABLE", "Operator journal is unavailable."
            ) from None
        self._require_directory(path, observed=observed)

    def _require_directory(
        self,
        path: Path,
        *,
        observed: os.stat_result | None = None,
    ) -> None:
        try:
            value = path.lstat() if observed is None else observed
        except OSError:
            raise fail(
                "OPERATOR_JOURNAL_UNAVAILABLE", "Operator journal is unavailable."
            ) from None
        if (
            not stat.S_ISDIR(value.st_mode)
            or stat.S_ISLNK(value.st_mode)
            or value.st_uid != self.owner_uid
            or value.st_gid != self.owner_gid
            or stat.S_IMODE(value.st_mode) != 0o700
        ):
            raise fail("OPERATOR_JOURNAL_UNTRUSTED", "Operator journal is not trusted.")

    def _open_lock(self, path: Path, *, create: bool) -> int:
        flags = (
            (os.O_RDWR if create else os.O_RDONLY)
            | (os.O_CREAT if create else 0)
            | getattr(os, "O_CLOEXEC", 0)
            | getattr(os, "O_NOFOLLOW", 0)
        )
        try:
            descriptor = os.open(path, flags, 0o600)
            if create:
                os.fchown(descriptor, self.owner_uid, self.owner_gid)
                os.fchmod(descriptor, 0o600)
            observed = os.fstat(descriptor)
        except OSError:
            raise fail(
                "OPERATOR_JOURNAL_UNAVAILABLE", "Operator journal is unavailable."
            ) from None
        if (
            not stat.S_ISREG(observed.st_mode)
            or observed.st_nlink != 1
            or observed.st_uid != self.owner_uid
            or observed.st_gid != self.owner_gid
            or stat.S_IMODE(observed.st_mode) != 0o600
        ):
            os.close(descriptor)
            raise fail("OPERATOR_JOURNAL_UNTRUSTED", "Operator journal is not trusted.")
        return descriptor

    def _read_active(self, directory: Path) -> OperatorTransaction | None:
        path = self.layout.operator_transaction_active
        try:
            path.lstat()
        except FileNotFoundError:
            return None
        except OSError:
            raise fail(
                "OPERATOR_JOURNAL_UNTRUSTED", "Operator journal is not trusted."
            ) from None
        return self._read(path)

    def _replace_active(
        self,
        record: OperatorTransaction,
        *,
        expected_identity: str,
    ) -> None:
        path = self.layout.operator_transaction_active
        current = self._read(path)
        if current.identity != expected_identity:
            raise fail("OPERATOR_JOURNAL_CHANGED", "Operator journal changed.")
        self._atomic_write(path, record, exclusive=False)

    def _finish(
        self,
        record: OperatorTransaction,
        *,
        expected_identity: str,
    ) -> None:
        if record.phase not in TERMINAL_PHASES:
            raise fail("OPERATOR_JOURNAL_INVALID", "Operator journal is invalid.")
        self._replace_active(record, expected_identity=expected_identity)
        directory, history = self._directories(create=False)
        self._archive_and_clear(record, directory, history)

    def _archive_and_clear(
        self,
        record: OperatorTransaction,
        directory: Path,
        history: Path,
    ) -> None:
        if record.phase not in TERMINAL_PHASES:
            raise fail("OPERATOR_JOURNAL_INVALID", "Operator journal is invalid.")
        history_path = history / f"{record.task_identity}.json"
        if history_path.exists() or history_path.is_symlink():
            if self._read(history_path) != record:
                raise fail(
                    "OPERATOR_JOURNAL_UNTRUSTED", "Operator journal is not trusted."
                )
        else:
            self._write_new(history_path, record)
        active_path = self.layout.operator_transaction_active
        current = self._read(active_path)
        if current != record:
            raise fail("OPERATOR_JOURNAL_CHANGED", "Operator journal changed.")
        try:
            active_path.unlink()
            _fsync_directory(directory)
        except OSError:
            raise fail(
                "OPERATOR_JOURNAL_UNAVAILABLE", "Operator journal is unavailable."
            ) from None

    def _write_new(self, path: Path, record: OperatorTransaction) -> None:
        self._atomic_write(path, record, exclusive=True)

    def _atomic_write(
        self,
        path: Path,
        record: OperatorTransaction,
        *,
        exclusive: bool,
    ) -> None:
        temporary = path.parent / f".{path.name}.{uuid.uuid4().hex}.tmp"
        flags = (
            os.O_WRONLY
            | os.O_CREAT
            | os.O_EXCL
            | getattr(os, "O_CLOEXEC", 0)
            | getattr(os, "O_NOFOLLOW", 0)
        )
        descriptor = -1
        try:
            if exclusive and (path.exists() or path.is_symlink()):
                raise OSError
            descriptor = os.open(temporary, flags, 0o600)
            os.fchown(descriptor, self.owner_uid, self.owner_gid)
            os.fchmod(descriptor, 0o600)
            content = canonical_json_bytes(record.to_dict())
            _write_all(descriptor, content)
            os.fsync(descriptor)
            os.close(descriptor)
            descriptor = -1
            if exclusive:
                os.link(temporary, path, follow_symlinks=False)
                temporary.unlink()
            else:
                os.replace(temporary, path)
            _fsync_directory(path.parent)
        except OSError:
            raise fail(
                "OPERATOR_JOURNAL_UNAVAILABLE", "Operator journal is unavailable."
            ) from None
        finally:
            if descriptor >= 0:
                os.close(descriptor)
            try:
                temporary.unlink()
            except FileNotFoundError:
                pass

    def _read(self, path: Path) -> OperatorTransaction:
        flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0)
        descriptor = -1
        try:
            descriptor = os.open(path, flags)
            before = os.fstat(descriptor)
            if (
                not stat.S_ISREG(before.st_mode)
                or before.st_nlink != 1
                or before.st_uid != self.owner_uid
                or before.st_gid != self.owner_gid
                or stat.S_IMODE(before.st_mode) != 0o600
                or not 0 < before.st_size <= _MAX_JOURNAL_BYTES
            ):
                raise OSError
            content = _read_bounded(descriptor, _MAX_JOURNAL_BYTES)
            after = os.fstat(descriptor)
            if len(content) != before.st_size or _witness(before) != _witness(after):
                raise OSError
        except OSError:
            raise fail(
                "OPERATOR_JOURNAL_UNTRUSTED", "Operator journal is not trusted."
            ) from None
        finally:
            if descriptor >= 0:
                os.close(descriptor)
        try:
            raw: Any = json.loads(content, object_pairs_hook=_unique_object)
            record = OperatorTransaction.from_dict(raw)
        except (DeploymentError, UnicodeDecodeError, ValueError, json.JSONDecodeError):
            raise fail(
                "OPERATOR_JOURNAL_UNTRUSTED", "Operator journal is not trusted."
            ) from None
        if canonical_json_bytes(record.to_dict()) != content:
            raise fail("OPERATOR_JOURNAL_UNTRUSTED", "Operator journal is not trusted.")
        return record


def _read_bounded(descriptor: int, maximum: int) -> bytes:
    values = bytearray()
    while len(values) <= maximum:
        chunk = os.read(descriptor, min(65536, maximum + 1 - len(values)))
        if not chunk:
            break
        values.extend(chunk)
    if len(values) > maximum:
        raise OSError
    return bytes(values)


def _write_all(descriptor: int, content: bytes) -> None:
    offset = 0
    while offset < len(content):
        written = os.write(descriptor, content[offset:])
        if written <= 0:
            raise OSError
        offset += written


def _fsync_directory(path: Path) -> None:
    descriptor = os.open(
        path,
        os.O_RDONLY
        | getattr(os, "O_CLOEXEC", 0)
        | getattr(os, "O_DIRECTORY", 0)
        | getattr(os, "O_NOFOLLOW", 0),
    )
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _witness(value: os.stat_result) -> tuple[int, ...]:
    return (
        value.st_dev,
        value.st_ino,
        value.st_mode,
        value.st_uid,
        value.st_gid,
        value.st_nlink,
        value.st_size,
        value.st_mtime_ns,
        value.st_ctime_ns,
    )


def _unique_object(pairs: list[tuple[str, object]]) -> dict[str, object]:
    value: dict[str, object] = {}
    for key, item in pairs:
        if key in value:
            raise ValueError
        value[key] = item
    return value


__all__ = [
    "EVIDENCE_KEYS",
    "JOURNAL_SCHEMA",
    "OperatorJournalHandle",
    "OperatorJournalStore",
    "OperatorJournalSummary",
    "OperatorRecoveryController",
    "OperatorTransaction",
    "PHASES",
    "RECOVERY_EVIDENCE_KEYS",
    "TERMINAL_PHASES",
]
