"""Immutable deployment-state generations with one atomic active pointer."""

from __future__ import annotations

from collections.abc import Callable, Mapping
from contextlib import contextmanager
import fcntl
import json
import os
from pathlib import Path
import stat
from typing import Iterator
import uuid

from encode_pipeline.deployment.canonical import canonical_json_bytes
from encode_pipeline.deployment.errors import DeploymentError, fail
from encode_pipeline.deployment.filesystem import (
    create_directory,
    fsync_directory,
    read_regular_file,
)
from encode_pipeline.deployment.layout import DeploymentLayout
from encode_pipeline.deployment.models import COMPONENTS, DeploymentState


STATE_FILENAME = "deployment-state.json"
MAX_STATE_BYTES = 256 * 1024
FaultInjector = Callable[[str], None]


class StateTransaction:
    """One flock-bound view used for compare-and-swap state transitions."""

    def __init__(
        self,
        store: "StateStore",
        *,
        exclusive: bool,
        expected_owner_uid: int | None,
        expected_owner_gid: int | None,
    ) -> None:
        self._store = store
        self.exclusive = exclusive
        self.expected_owner_uid = expected_owner_uid
        self.expected_owner_gid = expected_owner_gid

    def initialize(self) -> DeploymentState:
        if not self.exclusive:
            raise fail("DEPLOYMENT_STATE_INVALID", "Deployment state is invalid.")
        return self._store._initialize_locked(
            expected_owner_uid=self.expected_owner_uid,
            expected_owner_gid=self.expected_owner_gid,
        )

    def read(self) -> DeploymentState:
        return self._store._read_locked(
            expected_owner_uid=self.expected_owner_uid,
            expected_owner_gid=self.expected_owner_gid,
        )

    def commit(
        self,
        state: DeploymentState,
        *,
        operation: str,
        expected_current_identity: str | None,
        fault: FaultInjector | None = None,
    ) -> None:
        if not self.exclusive:
            raise fail("DEPLOYMENT_STATE_INVALID", "Deployment state is invalid.")
        self._store._commit_locked(
            state,
            operation=operation,
            expected_current_identity=expected_current_identity,
            expected_owner_uid=self.expected_owner_uid,
            expected_owner_gid=self.expected_owner_gid,
            fault=fault,
        )

    def pending_transactions(self) -> tuple[str, ...]:
        return self._store._pending_transactions_locked(
            expected_owner_uid=self.expected_owner_uid,
            expected_owner_gid=self.expected_owner_gid,
        )

    def referenced_identities(self) -> Mapping[str, frozenset[str]]:
        return self._store._referenced_identities_locked(
            expected_owner_uid=self.expected_owner_uid,
            expected_owner_gid=self.expected_owner_gid,
        )


class StateStore:
    """Persist state without mixing release bytes or user data into a generation."""

    def __init__(self, layout: DeploymentLayout) -> None:
        self.layout = layout

    def initialize(
        self,
        *,
        expected_owner_uid: int | None = None,
        expected_owner_gid: int | None = None,
    ) -> DeploymentState:
        with self.transaction(
            exclusive=True,
            expected_owner_uid=expected_owner_uid,
            expected_owner_gid=expected_owner_gid,
        ) as transaction:
            return transaction.initialize()

    def read(
        self,
        *,
        expected_owner_uid: int | None = None,
        expected_owner_gid: int | None = None,
    ) -> DeploymentState:
        with self.transaction(
            exclusive=False,
            expected_owner_uid=expected_owner_uid,
            expected_owner_gid=expected_owner_gid,
        ) as transaction:
            return transaction.read()

    def _read_locked(
        self,
        *,
        expected_owner_uid: int | None,
        expected_owner_gid: int | None,
    ) -> DeploymentState:
        link = self.layout.current_state_link
        try:
            observed = link.lstat()
            target = os.readlink(link)
        except OSError:
            raise fail(
                "DEPLOYMENT_STATE_UNAVAILABLE", "Deployment state is unavailable."
            ) from None
        if (
            not stat.S_ISLNK(observed.st_mode)
            or (
                expected_owner_uid is not None and observed.st_uid != expected_owner_uid
            )
            or (
                expected_owner_gid is not None and observed.st_gid != expected_owner_gid
            )
        ):
            raise fail("DEPLOYMENT_STATE_INVALID", "Deployment state is invalid.")
        target_path = Path(target)
        if (
            target_path.is_absolute()
            or len(target_path.parts) != 2
            or target_path.parts[0] != "generations"
            or target_path.parts[1] in {"", ".", ".."}
        ):
            raise fail("DEPLOYMENT_STATE_INVALID", "Deployment state is invalid.")
        generation = link.parent / target_path
        try:
            directory_stat = generation.lstat()
        except OSError:
            raise fail(
                "DEPLOYMENT_STATE_INVALID", "Deployment state is invalid."
            ) from None
        if (
            not stat.S_ISDIR(directory_stat.st_mode)
            or stat.S_ISLNK(directory_stat.st_mode)
            or stat.S_IMODE(directory_stat.st_mode) != 0o555
            or (
                expected_owner_uid is not None
                and directory_stat.st_uid != expected_owner_uid
            )
            or (
                expected_owner_gid is not None
                and directory_stat.st_gid != expected_owner_gid
            )
        ):
            raise fail("DEPLOYMENT_STATE_INVALID", "Deployment state is invalid.")
        content, file_stat = read_regular_file(
            generation / STATE_FILENAME,
            max_bytes=MAX_STATE_BYTES,
            code="DEPLOYMENT_STATE_INVALID",
        )
        if (
            stat.S_IMODE(file_stat.st_mode) != 0o444
            or (
                expected_owner_uid is not None
                and file_stat.st_uid != expected_owner_uid
            )
            or (
                expected_owner_gid is not None
                and file_stat.st_gid != expected_owner_gid
            )
        ):
            raise fail("DEPLOYMENT_STATE_INVALID", "Deployment state is invalid.")
        try:
            document = json.loads(content, object_pairs_hook=_unique_object)
        except (UnicodeDecodeError, json.JSONDecodeError, ValueError):
            raise fail(
                "DEPLOYMENT_STATE_INVALID", "Deployment state is invalid."
            ) from None
        state = DeploymentState.from_dict(document)
        if generation.name != state.identity or content != canonical_json_bytes(
            state.to_dict()
        ):
            raise fail("DEPLOYMENT_STATE_INVALID", "Deployment state is invalid.")
        return state

    def _initialize_locked(
        self,
        *,
        expected_owner_uid: int | None,
        expected_owner_gid: int | None,
    ) -> DeploymentState:
        current = self._current_state_locked(
            expected_owner_uid=expected_owner_uid,
            expected_owner_gid=expected_owner_gid,
        )
        if current is not None:
            return current
        state = DeploymentState.initial()
        self._commit_locked(
            state,
            operation="initialize",
            expected_current_identity=None,
            expected_owner_uid=expected_owner_uid,
            expected_owner_gid=expected_owner_gid,
        )
        return state

    def commit(
        self,
        state: DeploymentState,
        *,
        operation: str,
        expected_current_identity: str | None,
        expected_owner_uid: int | None = None,
        expected_owner_gid: int | None = None,
        fault: FaultInjector | None = None,
    ) -> None:
        with self.transaction(
            exclusive=True,
            expected_owner_uid=expected_owner_uid,
            expected_owner_gid=expected_owner_gid,
        ) as transaction:
            transaction.commit(
                state,
                operation=operation,
                expected_current_identity=expected_current_identity,
                fault=fault,
            )

    def _commit_locked(
        self,
        state: DeploymentState,
        *,
        operation: str,
        expected_current_identity: str | None,
        expected_owner_uid: int | None,
        expected_owner_gid: int | None,
        fault: FaultInjector | None = None,
    ) -> None:
        if not operation or not operation.replace("-", "").isalnum():
            raise fail(
                "DEPLOYMENT_OPERATION_INVALID", "Deployment operation is invalid."
            )
        current = self._current_state_locked(
            expected_owner_uid=expected_owner_uid,
            expected_owner_gid=expected_owner_gid,
        )
        current_identity = None if current is None else current.identity
        if (
            current_identity != expected_current_identity
            or (current is None and state.generation != 0)
            or (current is not None and state.generation != current.generation + 1)
        ):
            raise fail(
                "DEPLOYMENT_STATE_CHANGED",
                "Deployment state changed before the operation could commit.",
                recoverable=True,
            )
        if self._pending_transactions_locked(
            expected_owner_uid=expected_owner_uid,
            expected_owner_gid=expected_owner_gid,
        ):
            raise fail(
                "DEPLOYMENT_RECOVERY_REQUIRED",
                "An interrupted deployment operation requires recovery.",
                recoverable=True,
            )
        transaction_id = f"{state.generation:08d}-{uuid.uuid4().hex}"
        pending = self.layout.transactions / f"{transaction_id}.pending.json"
        complete = self.layout.transactions / f"{transaction_id}.complete.json"
        plan = {
            "schema_version": "helixweave-deployment-transaction-v1",
            "operation": operation,
            "target_state": state.identity,
        }
        self._write_immutable(pending, canonical_json_bytes(plan), mode=0o444)
        if fault is not None:
            fault("transaction-prepared")

        destination = self.layout.state_generations / state.identity
        if destination.exists() or destination.is_symlink():
            self._verify_generation(
                destination,
                state,
                expected_owner_uid=expected_owner_uid,
                expected_owner_gid=expected_owner_gid,
            )
        else:
            partial = self.layout.state_generations / f".partial-{transaction_id}"
            try:
                partial.mkdir(mode=0o755)
                self._write_immutable(
                    partial / STATE_FILENAME,
                    canonical_json_bytes(state.to_dict()),
                    mode=0o444,
                )
                os.chmod(partial, 0o555)
                fsync_directory(partial)
                os.replace(partial, destination)
                fsync_directory(self.layout.state_generations)
            except DeploymentError:
                raise
            except OSError:
                raise fail(
                    "DEPLOYMENT_STATE_WRITE_FAILED",
                    "Deployment state could not be committed.",
                    recoverable=True,
                ) from None
        if fault is not None:
            fault("generation-committed")

        temporary_link = self.layout.current_state_link.parent / (
            f".current-{transaction_id}"
        )
        try:
            os.symlink(f"generations/{state.identity}", temporary_link)
            os.replace(temporary_link, self.layout.current_state_link)
            fsync_directory(self.layout.current_state_link.parent)
        except OSError:
            raise fail(
                "DEPLOYMENT_ACTIVATION_FAILED",
                "Deployment activation did not change the active generation.",
                recoverable=True,
            ) from None
        if fault is not None:
            fault("pointer-committed")
        try:
            os.replace(pending, complete)
            fsync_directory(self.layout.transactions)
        except OSError:
            raise fail(
                "DEPLOYMENT_RECEIPT_WRITE_FAILED",
                "Deployment activated but its transaction receipt is incomplete.",
                recoverable=True,
            ) from None

    def pending_transactions(
        self,
        *,
        expected_owner_uid: int | None = None,
        expected_owner_gid: int | None = None,
    ) -> tuple[str, ...]:
        with self.transaction(
            exclusive=False,
            expected_owner_uid=expected_owner_uid,
            expected_owner_gid=expected_owner_gid,
        ) as transaction:
            return transaction.pending_transactions()

    def _pending_transactions_locked(
        self,
        *,
        expected_owner_uid: int | None,
        expected_owner_gid: int | None,
    ) -> tuple[str, ...]:
        try:
            observed = self.layout.transactions.lstat()
        except OSError:
            raise fail(
                "DEPLOYMENT_STATE_INVALID", "Deployment state is invalid."
            ) from None
        if (
            not stat.S_ISDIR(observed.st_mode)
            or stat.S_ISLNK(observed.st_mode)
            or (
                expected_owner_uid is not None and observed.st_uid != expected_owner_uid
            )
            or (
                expected_owner_gid is not None and observed.st_gid != expected_owner_gid
            )
        ):
            raise fail("DEPLOYMENT_STATE_INVALID", "Deployment state is invalid.")
        values: list[str] = []
        for path in sorted(self.layout.transactions.iterdir()):
            if not path.name.endswith(".pending.json"):
                continue
            try:
                item = path.lstat()
            except OSError:
                raise fail(
                    "DEPLOYMENT_STATE_INVALID", "Deployment state is invalid."
                ) from None
            if (
                not stat.S_ISREG(item.st_mode)
                or stat.S_ISLNK(item.st_mode)
                or item.st_nlink != 1
                or stat.S_IMODE(item.st_mode) != 0o444
                or (
                    expected_owner_uid is not None and item.st_uid != expected_owner_uid
                )
                or (
                    expected_owner_gid is not None and item.st_gid != expected_owner_gid
                )
            ):
                raise fail("DEPLOYMENT_STATE_INVALID", "Deployment state is invalid.")
            values.append(path.name.removesuffix(".pending.json"))
        return tuple(values)

    def _referenced_identities_locked(
        self,
        *,
        expected_owner_uid: int | None,
        expected_owner_gid: int | None,
    ) -> Mapping[str, frozenset[str]]:
        values: dict[str, set[str]] = {component: set() for component in COMPONENTS}
        try:
            paths = sorted(self.layout.state_generations.iterdir())
        except OSError:
            raise fail(
                "DEPLOYMENT_STATE_INVALID", "Deployment state is invalid."
            ) from None
        if len(paths) > 10_000:
            raise fail("DEPLOYMENT_STATE_INVALID", "Deployment state is invalid.")
        for path in paths:
            if path.name.startswith("."):
                continue
            try:
                directory = path.lstat()
            except OSError:
                raise fail(
                    "DEPLOYMENT_STATE_INVALID", "Deployment state is invalid."
                ) from None
            content, observed = read_regular_file(
                path / STATE_FILENAME,
                max_bytes=MAX_STATE_BYTES,
                code="DEPLOYMENT_STATE_INVALID",
            )
            if (
                not stat.S_ISDIR(directory.st_mode)
                or stat.S_ISLNK(directory.st_mode)
                or stat.S_IMODE(directory.st_mode) != 0o555
                or stat.S_IMODE(observed.st_mode) != 0o444
                or (
                    expected_owner_uid is not None
                    and (
                        directory.st_uid != expected_owner_uid
                        or observed.st_uid != expected_owner_uid
                    )
                )
                or (
                    expected_owner_gid is not None
                    and (
                        directory.st_gid != expected_owner_gid
                        or observed.st_gid != expected_owner_gid
                    )
                )
            ):
                raise fail("DEPLOYMENT_STATE_INVALID", "Deployment state is invalid.")
            try:
                document = json.loads(content, object_pairs_hook=_unique_object)
            except (UnicodeDecodeError, json.JSONDecodeError, ValueError):
                raise fail(
                    "DEPLOYMENT_STATE_INVALID", "Deployment state is invalid."
                ) from None
            state = DeploymentState.from_dict(document)
            if path.name != state.identity or content != canonical_json_bytes(
                state.to_dict()
            ):
                raise fail("DEPLOYMENT_STATE_INVALID", "Deployment state is invalid.")
            for component in COMPONENTS:
                slots = state.components[component]
                values[component].update(
                    identity
                    for identity in (slots.active, slots.previous, slots.staged)
                    if identity is not None
                )
        return {
            component: frozenset(identities) for component, identities in values.items()
        }

    def _verify_generation(
        self,
        path: Path,
        expected: DeploymentState,
        *,
        expected_owner_uid: int | None,
        expected_owner_gid: int | None,
    ) -> None:
        try:
            directory = path.lstat()
        except OSError:
            raise fail(
                "DEPLOYMENT_STATE_INVALID", "Deployment state is invalid."
            ) from None
        content, observed = read_regular_file(
            path / STATE_FILENAME,
            max_bytes=MAX_STATE_BYTES,
            code="DEPLOYMENT_STATE_INVALID",
        )
        if (
            not stat.S_ISDIR(directory.st_mode)
            or stat.S_ISLNK(directory.st_mode)
            or stat.S_IMODE(directory.st_mode) != 0o555
            or stat.S_IMODE(observed.st_mode) != 0o444
            or (
                expected_owner_uid is not None
                and (
                    directory.st_uid != expected_owner_uid
                    or observed.st_uid != expected_owner_uid
                )
            )
            or (
                expected_owner_gid is not None
                and (
                    directory.st_gid != expected_owner_gid
                    or observed.st_gid != expected_owner_gid
                )
            )
            or content != canonical_json_bytes(expected.to_dict())
        ):
            raise fail("DEPLOYMENT_STATE_INVALID", "Deployment state is invalid.")

    @contextmanager
    def transaction(
        self,
        *,
        exclusive: bool,
        expected_owner_uid: int | None = None,
        expected_owner_gid: int | None = None,
    ) -> Iterator[StateTransaction]:
        if exclusive:
            self._ensure_roots(
                expected_owner_uid=expected_owner_uid,
                expected_owner_gid=expected_owner_gid,
            )
            self._ensure_lock_file(
                expected_owner_uid=expected_owner_uid,
                expected_owner_gid=expected_owner_gid,
            )
        else:
            self._validate_roots(
                expected_owner_uid=expected_owner_uid,
                expected_owner_gid=expected_owner_gid,
            )
        descriptor = self._open_lock(
            expected_owner_uid=expected_owner_uid,
            expected_owner_gid=expected_owner_gid,
        )
        operation = fcntl.LOCK_EX if exclusive else fcntl.LOCK_SH
        try:
            try:
                fcntl.flock(descriptor, operation | fcntl.LOCK_NB)
            except BlockingIOError:
                raise fail(
                    "DEPLOYMENT_BUSY",
                    "Another deployment operation is in progress.",
                    recoverable=True,
                ) from None
            yield StateTransaction(
                self,
                exclusive=exclusive,
                expected_owner_uid=expected_owner_uid,
                expected_owner_gid=expected_owner_gid,
            )
        finally:
            try:
                fcntl.flock(descriptor, fcntl.LOCK_UN)
            finally:
                os.close(descriptor)

    def _ensure_roots(
        self,
        *,
        expected_owner_uid: int | None,
        expected_owner_gid: int | None,
    ) -> None:
        data_root = self.layout.data_root
        operator_root = data_root / "operator"
        state_root = self.layout.current_state_link.parent
        if not data_root.exists() and not data_root.is_symlink():
            create_directory(data_root, mode=0o700)
        self._require_trusted_directory(
            data_root,
            expected_owner_uid=expected_owner_uid,
            expected_owner_gid=None,
            mode=None,
        )
        if not operator_root.exists() and not operator_root.is_symlink():
            try:
                operator_root.mkdir(mode=0o700)
            except OSError:
                raise fail(
                    "DEPLOYMENT_STORAGE_UNAVAILABLE",
                    "Deployment storage is unavailable.",
                ) from None
        self._require_trusted_directory(
            operator_root,
            expected_owner_uid=expected_owner_uid,
            expected_owner_gid=None,
            mode=None,
        )
        for path in (
            state_root,
            self.layout.state_generations,
            self.layout.transactions,
        ):
            if not path.exists() and not path.is_symlink():
                try:
                    path.mkdir(mode=0o700)
                except OSError:
                    raise fail(
                        "DEPLOYMENT_STORAGE_UNAVAILABLE",
                        "Deployment storage is unavailable.",
                    ) from None
            self._require_trusted_directory(
                path,
                expected_owner_uid=expected_owner_uid,
                expected_owner_gid=expected_owner_gid,
                mode=0o700,
            )

    def _validate_roots(
        self,
        *,
        expected_owner_uid: int | None,
        expected_owner_gid: int | None,
    ) -> None:
        for path in (
            self.layout.data_root,
            self.layout.data_root / "operator",
        ):
            self._require_trusted_directory(
                path,
                expected_owner_uid=expected_owner_uid,
                expected_owner_gid=None,
                mode=None,
            )
        for path in (
            self.layout.current_state_link.parent,
            self.layout.state_generations,
            self.layout.transactions,
        ):
            self._require_trusted_directory(
                path,
                expected_owner_uid=expected_owner_uid,
                expected_owner_gid=expected_owner_gid,
                mode=0o700,
            )

    def _ensure_lock_file(
        self,
        *,
        expected_owner_uid: int | None,
        expected_owner_gid: int | None,
    ) -> None:
        flags = (
            os.O_RDWR
            | os.O_CREAT
            | os.O_EXCL
            | getattr(os, "O_CLOEXEC", 0)
            | getattr(os, "O_NOFOLLOW", 0)
        )
        try:
            descriptor = os.open(self.layout.state_lock, flags, 0o600)
        except FileExistsError:
            return
        except OSError:
            raise fail(
                "DEPLOYMENT_STATE_WRITE_FAILED",
                "Deployment state could not be committed.",
                recoverable=True,
            ) from None
        try:
            os.fsync(descriptor)
            observed = os.fstat(descriptor)
            if (
                expected_owner_uid is not None and observed.st_uid != expected_owner_uid
            ) or (
                expected_owner_gid is not None and observed.st_gid != expected_owner_gid
            ):
                raise fail(
                    "DEPLOYMENT_STATE_WRITE_FAILED",
                    "Deployment state could not be committed.",
                    recoverable=True,
                )
        finally:
            os.close(descriptor)
        fsync_directory(self.layout.state_lock.parent)

    def _open_lock(
        self,
        *,
        expected_owner_uid: int | None,
        expected_owner_gid: int | None,
    ) -> int:
        flags = os.O_RDWR | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0)
        try:
            descriptor = os.open(self.layout.state_lock, flags)
        except OSError:
            raise fail(
                "DEPLOYMENT_STATE_UNAVAILABLE", "Deployment state is unavailable."
            ) from None
        observed = os.fstat(descriptor)
        if (
            not stat.S_ISREG(observed.st_mode)
            or observed.st_nlink != 1
            or stat.S_IMODE(observed.st_mode) != 0o600
            or (
                expected_owner_uid is not None and observed.st_uid != expected_owner_uid
            )
            or (
                expected_owner_gid is not None and observed.st_gid != expected_owner_gid
            )
        ):
            os.close(descriptor)
            raise fail("DEPLOYMENT_STATE_INVALID", "Deployment state is invalid.")
        return descriptor

    def _current_state_locked(
        self,
        *,
        expected_owner_uid: int | None,
        expected_owner_gid: int | None,
    ) -> DeploymentState | None:
        try:
            self.layout.current_state_link.lstat()
        except FileNotFoundError:
            return None
        except OSError:
            raise fail(
                "DEPLOYMENT_STATE_INVALID", "Deployment state is invalid."
            ) from None
        return self._read_locked(
            expected_owner_uid=expected_owner_uid,
            expected_owner_gid=expected_owner_gid,
        )

    @staticmethod
    def _require_trusted_directory(
        path: Path,
        *,
        expected_owner_uid: int | None,
        expected_owner_gid: int | None,
        mode: int | None,
    ) -> None:
        try:
            observed = path.lstat()
        except OSError:
            raise fail(
                "DEPLOYMENT_STATE_INVALID", "Deployment state is invalid."
            ) from None
        if (
            not stat.S_ISDIR(observed.st_mode)
            or stat.S_ISLNK(observed.st_mode)
            or stat.S_IMODE(observed.st_mode) & 0o022
            or (mode is not None and stat.S_IMODE(observed.st_mode) != mode)
            or (
                expected_owner_uid is not None and observed.st_uid != expected_owner_uid
            )
            or (
                expected_owner_gid is not None and observed.st_gid != expected_owner_gid
            )
        ):
            raise fail("DEPLOYMENT_STATE_INVALID", "Deployment state is invalid.")

    @staticmethod
    def _write_immutable(path: Path, content: bytes, *, mode: int) -> None:
        flags = (
            os.O_WRONLY
            | os.O_CREAT
            | os.O_EXCL
            | getattr(os, "O_CLOEXEC", 0)
            | getattr(os, "O_NOFOLLOW", 0)
        )
        try:
            descriptor = os.open(path, flags, 0o600)
            try:
                written = 0
                while written < len(content):
                    written += os.write(descriptor, content[written:])
                os.fchmod(descriptor, mode)
                os.fsync(descriptor)
            finally:
                os.close(descriptor)
        except OSError:
            raise fail(
                "DEPLOYMENT_STATE_WRITE_FAILED",
                "Deployment state could not be committed.",
                recoverable=True,
            ) from None


def _unique_object(pairs: list[tuple[str, object]]) -> dict[str, object]:
    value: dict[str, object] = {}
    for key, item in pairs:
        if key in value:
            raise ValueError("duplicate JSON key")
        value[key] = item
    return value


__all__ = ["FaultInjector", "StateStore"]
