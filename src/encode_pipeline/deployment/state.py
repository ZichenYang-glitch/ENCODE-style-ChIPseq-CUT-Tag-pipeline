"""Immutable deployment-state generations with one atomic active pointer."""

from __future__ import annotations

from collections.abc import Callable, Mapping
from contextlib import contextmanager
from dataclasses import dataclass
import fcntl
import hashlib
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
from encode_pipeline.deployment.models import (
    BULK_RNASEQ_RUNTIME,
    COMPONENTS,
    ENCODE_RUNTIME,
    PLATFORM,
    DeploymentState,
)


STATE_FILENAME = "deployment-state.json"
PLATFORM_ENV_FILENAME = "platform.env"
MAX_STATE_BYTES = 256 * 1024
MAX_PLATFORM_ENV_BYTES = 64 * 1024
STATE_TRANSACTION_SCHEMA = "helixweave-deployment-transaction-v2"
STATE_RECOVERY_CHOICES = frozenset({"restore-prior", "complete-candidate"})
MAX_STATE_TRANSACTION_ENTRIES = 10_000
FaultInjector = Callable[[str], None]
_NO_COMPLETED_RECOVERY = object()

PLATFORM_ENV_KEYS = (
    "ENCODE_PIPELINE_DATABASE_URL",
    "ENCODE_PIPELINE_WORKSPACE_ROOT",
    "ENCODE_PIPELINE_REDIS_URL",
    "ENCODE_PIPELINE_QUEUE_NAME",
    "ENCODE_PIPELINE_REFERENCE_PROFILE_CONFIG",
    "ENCODE_PIPELINE_REDIS_CONNECT_TIMEOUT_SECONDS",
    "ENCODE_PIPELINE_REDIS_API_READ_TIMEOUT_SECONDS",
    "HELIXWEAVE_ENCODE_RUNTIME_ROOT",
    "HELIXWEAVE_ENCODE_RUNNER_ROOT",
    "HELIXWEAVE_ENCODE_CONDA_PREFIX",
    "HELIXWEAVE_DEPLOYMENT_IDENTITY",
    "HELIXWEAVE_ACTIVE_PLATFORM_IDENTITY",
    "HELIXWEAVE_ACTIVE_API_CONTRACT_SHA256",
    "HELIXWEAVE_BULK_RNASEQ_RUNTIME_ROOT",
    "ENCODE_PIPELINE_MANAGED_DOCKER_EXECUTABLE",
    "ENCODE_PIPELINE_MANAGED_DOCKER_SOCKET",
)


@dataclass(frozen=True)
class PlatformEnvironment:
    content: bytes
    identity: str
    api_contract_sha256: str


@dataclass(frozen=True)
class StateRecoveryPlan:
    transaction_identity: str
    operation: str
    prior_state_identity: str | None
    candidate_state_identity: str

    @classmethod
    def from_bytes(
        cls,
        transaction_identity: str,
        content: bytes,
    ) -> "StateRecoveryPlan":
        try:
            document = json.loads(content, object_pairs_hook=_unique_object)
        except (UnicodeDecodeError, json.JSONDecodeError, ValueError):
            raise fail(
                "DEPLOYMENT_STATE_INVALID", "Deployment state is invalid."
            ) from None
        if (
            not isinstance(document, dict)
            or set(document)
            != {"schema_version", "operation", "prior_state", "target_state"}
            or document["schema_version"] != STATE_TRANSACTION_SCHEMA
            or not isinstance(document["operation"], str)
            or not document["operation"]
            or not document["operation"].replace("-", "").isalnum()
            or not _valid_state_identity(document["target_state"])
            or (
                document["prior_state"] is not None
                and not _valid_state_identity(document["prior_state"])
            )
            or content != canonical_json_bytes(document)
        ):
            raise fail("DEPLOYMENT_STATE_INVALID", "Deployment state is invalid.")
        return cls(
            transaction_identity=transaction_identity,
            operation=document["operation"],
            prior_state_identity=document["prior_state"],
            candidate_state_identity=document["target_state"],
        )


def render_platform_environment(
    layout: DeploymentLayout,
    state: DeploymentState,
    *,
    api_contract_sha256: str,
) -> PlatformEnvironment:
    """Render the only non-secret environment accepted in a state generation."""
    if (
        not isinstance(api_contract_sha256, str)
        or len(api_contract_sha256) != 64
        or any(character not in "0123456789abcdef" for character in api_contract_sha256)
    ):
        raise fail(
            "OPERATOR_CONFIGURATION_INVALID", "Operator configuration is invalid."
        )
    active = {component: state.components[component].active for component in COMPONENTS}
    platform_identity = active[PLATFORM]
    if platform_identity is None:
        raise fail(
            "OPERATOR_CONFIGURATION_INVALID", "Operator configuration is invalid."
        )
    encode_root = ""
    encode_runner_root = ""
    encode_conda_prefix = ""
    if active[ENCODE_RUNTIME] is not None:
        encode_root = str(
            layout.encode_runtimes
            / active[ENCODE_RUNTIME]
            / "payload"
            / "contracts"
            / "encode-runtime"
        )
        materialized = layout.encode_runtime_active_root(active[ENCODE_RUNTIME])
        encode_runner_root = str(materialized / "runner")
        encode_conda_prefix = str(materialized / "conda-envs")
    bulk_root = ""
    if active[BULK_RNASEQ_RUNTIME] is not None:
        bulk_root = str(
            layout.bulk_rnaseq_runtimes
            / active[BULK_RNASEQ_RUNTIME]
            / "payload"
            / "runtime"
        )
    values = {
        "ENCODE_PIPELINE_DATABASE_URL": f"sqlite:///{layout.database}",
        "ENCODE_PIPELINE_WORKSPACE_ROOT": str(layout.workspaces),
        "ENCODE_PIPELINE_REDIS_URL": f"unix://{layout.run_root / 'redis' / 'redis.sock'}",
        "ENCODE_PIPELINE_QUEUE_NAME": "encode-runs",
        "ENCODE_PIPELINE_REFERENCE_PROFILE_CONFIG": str(
            layout.reference_profile_config
        ),
        "ENCODE_PIPELINE_REDIS_CONNECT_TIMEOUT_SECONDS": "2",
        "ENCODE_PIPELINE_REDIS_API_READ_TIMEOUT_SECONDS": "5",
        "HELIXWEAVE_ENCODE_RUNTIME_ROOT": encode_root,
        "HELIXWEAVE_ENCODE_RUNNER_ROOT": encode_runner_root,
        "HELIXWEAVE_ENCODE_CONDA_PREFIX": encode_conda_prefix,
        "HELIXWEAVE_DEPLOYMENT_IDENTITY": state.identity,
        "HELIXWEAVE_ACTIVE_PLATFORM_IDENTITY": platform_identity,
        "HELIXWEAVE_ACTIVE_API_CONTRACT_SHA256": api_contract_sha256,
        "HELIXWEAVE_BULK_RNASEQ_RUNTIME_ROOT": bulk_root,
        "ENCODE_PIPELINE_MANAGED_DOCKER_EXECUTABLE": "/usr/bin/docker",
        "ENCODE_PIPELINE_MANAGED_DOCKER_SOCKET": str(
            layout.run_root / "docker" / "docker.sock"
        ),
    }
    content = "".join(f"{key}={values[key]}\n" for key in PLATFORM_ENV_KEYS).encode(
        "ascii"
    )
    return PlatformEnvironment(
        content=content,
        identity=f"sha256-{hashlib.sha256(content).hexdigest()}",
        api_contract_sha256=api_contract_sha256,
    )


def parse_platform_environment(
    layout: DeploymentLayout,
    state: DeploymentState,
    content: bytes,
) -> PlatformEnvironment:
    """Parse and fully re-derive one generation-bound environment."""
    try:
        if not 0 < len(content) <= MAX_PLATFORM_ENV_BYTES or not content.endswith(
            b"\n"
        ):
            raise ValueError
        lines = content.decode("ascii").splitlines()
        if len(lines) != len(PLATFORM_ENV_KEYS):
            raise ValueError
        values: dict[str, str] = {}
        observed_keys: list[str] = []
        for line in lines:
            key, separator, value = line.partition("=")
            if (
                separator != "="
                or not key
                or key in values
                or any(character in value for character in ("\x00", "\r", "\n"))
            ):
                raise ValueError
            observed_keys.append(key)
            values[key] = value
        if tuple(observed_keys) != PLATFORM_ENV_KEYS:
            raise ValueError
        expected = render_platform_environment(
            layout,
            state,
            api_contract_sha256=values["HELIXWEAVE_ACTIVE_API_CONTRACT_SHA256"],
        )
        if expected.content != content:
            raise ValueError
        return expected
    except (UnicodeError, ValueError, DeploymentError):
        raise fail(
            "DEPLOYMENT_CONFIGURATION_INVALID",
            "Deployment configuration is invalid.",
        ) from None


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
        platform_environment: PlatformEnvironment | None = None,
        fault: FaultInjector | None = None,
    ) -> None:
        if not self.exclusive:
            raise fail("DEPLOYMENT_STATE_INVALID", "Deployment state is invalid.")
        self._store._commit_locked(
            state,
            operation=operation,
            expected_current_identity=expected_current_identity,
            platform_environment=platform_environment,
            expected_owner_uid=self.expected_owner_uid,
            expected_owner_gid=self.expected_owner_gid,
            fault=fault,
        )

    def pending_transactions(self) -> tuple[str, ...]:
        return self._store._pending_transactions_locked(
            expected_owner_uid=self.expected_owner_uid,
            expected_owner_gid=self.expected_owner_gid,
        )

    def recover(
        self,
        *,
        prior_state_identity: str | None,
        candidate_state_identity: str,
        desired: str,
        fault: FaultInjector | None = None,
    ) -> DeploymentState | None:
        if not self.exclusive:
            raise fail("DEPLOYMENT_STATE_INVALID", "Deployment state is invalid.")
        return self._store._recover_locked(
            prior_state_identity=prior_state_identity,
            candidate_state_identity=candidate_state_identity,
            desired=desired,
            expected_owner_uid=self.expected_owner_uid,
            expected_owner_gid=self.expected_owner_gid,
            fault=fault,
        )

    def referenced_identities(self) -> Mapping[str, frozenset[str]]:
        return self._store._referenced_identities_locked(
            expected_owner_uid=self.expected_owner_uid,
            expected_owner_gid=self.expected_owner_gid,
        )


class StateStore:
    """Persist state without mixing release bytes or user data into a generation."""

    def __init__(
        self,
        layout: DeploymentLayout,
        *,
        reader_gid: int | None = None,
        service_gid: int | None = None,
        verify_environment_content: bool = True,
    ) -> None:
        if reader_gid is not None and (
            not isinstance(reader_gid, int)
            or isinstance(reader_gid, bool)
            or reader_gid < 0
        ):
            raise fail("DEPLOYMENT_STATE_INVALID", "Deployment state is invalid.")
        if service_gid is not None and (
            not isinstance(service_gid, int)
            or isinstance(service_gid, bool)
            or service_gid < 0
        ):
            raise fail("DEPLOYMENT_STATE_INVALID", "Deployment state is invalid.")
        if not isinstance(verify_environment_content, bool):
            raise fail("DEPLOYMENT_STATE_INVALID", "Deployment state is invalid.")
        self.layout = layout
        self.reader_gid = reader_gid
        self.service_gid = service_gid
        self.verify_environment_content = verify_environment_content

    def _private_gid(self, owner_gid: int | None) -> int | None:
        return self.reader_gid if self.reader_gid is not None else owner_gid

    def _environment_gid(self, owner_gid: int | None) -> int | None:
        return self.service_gid if self.service_gid is not None else owner_gid

    @property
    def _private_root_mode(self) -> int:
        return 0o750 if self.reader_gid is not None else 0o700

    @property
    def _operator_root_mode(self) -> int:
        return 0o710 if self.reader_gid is not None else 0o700

    @property
    def _public_root_mode(self) -> int:
        return 0o755

    @property
    def _generation_mode(self) -> int:
        return 0o555

    @property
    def _state_file_mode(self) -> int:
        return 0o444

    @property
    def _environment_file_mode(self) -> int:
        return 0o440

    @property
    def _lock_mode(self) -> int:
        return 0o640 if self.reader_gid is not None else 0o600

    @property
    def _transaction_file_mode(self) -> int:
        return 0o440 if self.reader_gid is not None else 0o400

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
            or stat.S_IMODE(directory_stat.st_mode) != self._generation_mode
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
            stat.S_IMODE(file_stat.st_mode) != self._state_file_mode
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
        self._verify_platform_environment(
            generation,
            state,
            expected_owner_uid=expected_owner_uid,
            expected_owner_gid=expected_owner_gid,
        )
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
            platform_environment=None,
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
        platform_environment: PlatformEnvironment | None = None,
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
                platform_environment=platform_environment,
                fault=fault,
            )

    def recover_pending_transaction(
        self,
        *,
        prior_state_identity: str | None,
        candidate_state_identity: str,
        desired: str,
        expected_owner_uid: int | None = None,
        expected_owner_gid: int | None = None,
        fault: FaultInjector | None = None,
    ) -> DeploymentState | None:
        with self.transaction(
            exclusive=True,
            expected_owner_uid=expected_owner_uid,
            expected_owner_gid=expected_owner_gid,
        ) as transaction:
            return transaction.recover(
                prior_state_identity=prior_state_identity,
                candidate_state_identity=candidate_state_identity,
                desired=desired,
                fault=fault,
            )

    def _commit_locked(
        self,
        state: DeploymentState,
        *,
        operation: str,
        expected_current_identity: str | None,
        platform_environment: PlatformEnvironment | None,
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
        selected_environment = self._select_platform_environment(
            current,
            state,
            platform_environment,
            expected_owner_uid=expected_owner_uid,
            expected_owner_gid=expected_owner_gid,
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
        pending = self.layout.state_transactions / f"{transaction_id}.pending.json"
        complete = self.layout.state_transactions / f"{transaction_id}.complete.json"
        plan = {
            "schema_version": STATE_TRANSACTION_SCHEMA,
            "operation": operation,
            "prior_state": current_identity,
            "target_state": state.identity,
        }
        self._write_immutable(
            pending,
            canonical_json_bytes(plan),
            mode=self._transaction_file_mode,
            owner_uid=expected_owner_uid,
            owner_gid=self._private_gid(expected_owner_gid),
        )
        if fault is not None:
            fault("transaction-prepared")

        destination = self.layout.state_generations / state.identity
        if destination.exists() or destination.is_symlink():
            self._verify_generation(
                destination,
                state,
                expected_owner_uid=expected_owner_uid,
                expected_owner_gid=expected_owner_gid,
                expected_environment=selected_environment,
            )
        else:
            partial = self.layout.state_generations / f".partial-{transaction_id}"
            try:
                partial.mkdir(mode=0o755)
                if expected_owner_uid is not None:
                    os.chown(
                        partial,
                        expected_owner_uid,
                        -1 if expected_owner_gid is None else expected_owner_gid,
                    )
                self._write_immutable(
                    partial / STATE_FILENAME,
                    canonical_json_bytes(state.to_dict()),
                    mode=self._state_file_mode,
                    owner_uid=expected_owner_uid,
                    owner_gid=expected_owner_gid,
                )
                if selected_environment is not None:
                    self._write_immutable(
                        partial / PLATFORM_ENV_FILENAME,
                        selected_environment.content,
                        mode=self._environment_file_mode,
                        owner_uid=expected_owner_uid,
                        owner_gid=self._environment_gid(expected_owner_gid),
                    )
                os.chmod(partial, self._generation_mode)
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

        self._switch_pointer_locked(
            state.identity,
            transaction_identity=transaction_id,
            expected_owner_uid=expected_owner_uid,
            expected_owner_gid=expected_owner_gid,
        )
        if fault is not None:
            fault("pointer-committed")
        try:
            os.replace(pending, complete)
            fsync_directory(self.layout.state_transactions)
        except OSError:
            raise fail(
                "DEPLOYMENT_RECEIPT_WRITE_FAILED",
                "Deployment activated but its transaction receipt is incomplete.",
                recoverable=True,
            ) from None

    def _recover_locked(
        self,
        *,
        prior_state_identity: str | None,
        candidate_state_identity: str,
        desired: str,
        expected_owner_uid: int | None,
        expected_owner_gid: int | None,
        fault: FaultInjector | None,
    ) -> DeploymentState | None:
        if (
            desired not in STATE_RECOVERY_CHOICES
            or not _valid_state_identity(candidate_state_identity)
            or (
                prior_state_identity is not None
                and not _valid_state_identity(prior_state_identity)
            )
            or not self.verify_environment_content
        ):
            raise fail("DEPLOYMENT_RECOVERY_INVALID", "Deployment recovery is invalid.")
        pending = self._pending_plans_locked(
            expected_owner_uid=expected_owner_uid,
            expected_owner_gid=expected_owner_gid,
        )
        desired_identity = (
            prior_state_identity
            if desired == "restore-prior"
            else candidate_state_identity
        )
        current = self._current_state_locked(
            expected_owner_uid=expected_owner_uid,
            expected_owner_gid=expected_owner_gid,
        )
        current_identity = None if current is None else current.identity
        if not pending:
            if desired == "restore-prior":
                recovered = self._recover_completed_transaction_locked(
                    prior_state_identity=prior_state_identity,
                    candidate_state_identity=candidate_state_identity,
                    current_identity=current_identity,
                    expected_owner_uid=expected_owner_uid,
                    expected_owner_gid=expected_owner_gid,
                    fault=fault,
                )
                if recovered is not _NO_COMPLETED_RECOVERY:
                    return recovered
            if current_identity != desired_identity:
                raise fail(
                    "DEPLOYMENT_RECOVERY_REQUIRED",
                    "Deployment state requires recovery.",
                    recoverable=True,
                )
            if desired_identity is None:
                return None
            return self._load_generation_locked(
                desired_identity,
                expected_owner_uid=expected_owner_uid,
                expected_owner_gid=expected_owner_gid,
            )
        if len(pending) != 1:
            raise fail(
                "DEPLOYMENT_RECOVERY_REQUIRED",
                "Deployment state requires recovery.",
                recoverable=True,
            )
        pending_path, plan = pending[0]
        if (
            plan.prior_state_identity != prior_state_identity
            or plan.candidate_state_identity != candidate_state_identity
            or current_identity not in {prior_state_identity, candidate_state_identity}
        ):
            raise fail(
                "DEPLOYMENT_RECOVERY_IDENTITY_MISMATCH",
                "Deployment recovery identities do not match.",
            )
        prior = None
        if prior_state_identity is not None:
            prior = self._load_generation_locked(
                prior_state_identity,
                expected_owner_uid=expected_owner_uid,
                expected_owner_gid=expected_owner_gid,
            )
        candidate = None
        candidate_path = self.layout.state_generations / candidate_state_identity
        if candidate_path.exists() or candidate_path.is_symlink():
            candidate = self._load_generation_locked(
                candidate_state_identity,
                expected_owner_uid=expected_owner_uid,
                expected_owner_gid=expected_owner_gid,
            )
        if desired == "complete-candidate" and candidate is None:
            raise fail(
                "DEPLOYMENT_RECOVERY_REQUIRED",
                "Deployment state requires recovery.",
                recoverable=True,
            )
        selected = prior if desired == "restore-prior" else candidate
        if current_identity != desired_identity:
            self._switch_pointer_locked(
                desired_identity,
                transaction_identity=f"recovery-{plan.transaction_identity}",
                expected_owner_uid=expected_owner_uid,
                expected_owner_gid=expected_owner_gid,
            )
        if fault is not None:
            fault("recovery-pointer-committed")
        suffix = "recovered.json" if desired == "restore-prior" else "complete.json"
        receipt = self.layout.state_transactions / (
            f"{plan.transaction_identity}.{suffix}"
        )
        try:
            if receipt.exists() or receipt.is_symlink():
                raise OSError
            os.rename(pending_path, receipt)
            fsync_directory(self.layout.state_transactions)
        except OSError:
            raise fail(
                "DEPLOYMENT_RECEIPT_WRITE_FAILED",
                "Deployment recovery receipt could not be committed.",
                recoverable=True,
            ) from None
        if fault is not None:
            fault("recovery-receipt-committed")
        return selected

    def _recover_completed_transaction_locked(
        self,
        *,
        prior_state_identity: str | None,
        candidate_state_identity: str,
        current_identity: str | None,
        expected_owner_uid: int | None,
        expected_owner_gid: int | None,
        fault: FaultInjector | None,
    ) -> DeploymentState | None | object:
        """Reverse the commit/journal handoff window using its durable receipt."""
        complete = self._matching_transaction_receipt_locked(
            suffix=".complete.json",
            prior_state_identity=prior_state_identity,
            candidate_state_identity=candidate_state_identity,
            expected_owner_uid=expected_owner_uid,
            expected_owner_gid=expected_owner_gid,
        )
        recovered = self._matching_transaction_receipt_locked(
            suffix=".recovered.json",
            prior_state_identity=prior_state_identity,
            candidate_state_identity=candidate_state_identity,
            expected_owner_uid=expected_owner_uid,
            expected_owner_gid=expected_owner_gid,
        )
        if complete is not None and recovered is not None:
            raise fail("DEPLOYMENT_STATE_INVALID", "Deployment state is invalid.")
        if complete is None and recovered is None:
            return _NO_COMPLETED_RECOVERY
        prior = None
        if prior_state_identity is not None:
            prior = self._load_generation_locked(
                prior_state_identity,
                expected_owner_uid=expected_owner_uid,
                expected_owner_gid=expected_owner_gid,
            )
        self._load_generation_locked(
            candidate_state_identity,
            expected_owner_uid=expected_owner_uid,
            expected_owner_gid=expected_owner_gid,
        )
        if current_identity not in {prior_state_identity, candidate_state_identity}:
            raise fail(
                "DEPLOYMENT_RECOVERY_IDENTITY_MISMATCH",
                "Deployment recovery identities do not match.",
            )
        if recovered is not None:
            if current_identity != prior_state_identity:
                raise fail(
                    "DEPLOYMENT_RECOVERY_REQUIRED",
                    "Deployment state requires recovery.",
                    recoverable=True,
                )
            return prior
        assert complete is not None
        complete_path, plan = complete
        if current_identity != prior_state_identity:
            self._switch_pointer_locked(
                prior_state_identity,
                transaction_identity=f"recovery-{plan.transaction_identity}",
                expected_owner_uid=expected_owner_uid,
                expected_owner_gid=expected_owner_gid,
            )
        if fault is not None:
            fault("recovery-pointer-committed")
        receipt = self.layout.state_transactions / (
            f"{plan.transaction_identity}.recovered.json"
        )
        try:
            if receipt.exists() or receipt.is_symlink():
                raise OSError
            os.rename(complete_path, receipt)
            fsync_directory(self.layout.state_transactions)
        except OSError:
            raise fail(
                "DEPLOYMENT_RECEIPT_WRITE_FAILED",
                "Deployment recovery receipt could not be committed.",
                recoverable=True,
            ) from None
        if fault is not None:
            fault("recovery-receipt-committed")
        return prior

    def _matching_transaction_receipt_locked(
        self,
        *,
        suffix: str,
        prior_state_identity: str | None,
        candidate_state_identity: str,
        expected_owner_uid: int | None,
        expected_owner_gid: int | None,
    ) -> tuple[Path, StateRecoveryPlan] | None:
        try:
            entries = sorted(self.layout.state_transactions.iterdir())
        except OSError:
            raise fail(
                "DEPLOYMENT_STATE_INVALID", "Deployment state is invalid."
            ) from None
        if len(entries) > MAX_STATE_TRANSACTION_ENTRIES:
            raise fail("DEPLOYMENT_STATE_INVALID", "Deployment state is invalid.")
        matches: list[tuple[Path, StateRecoveryPlan]] = []
        for path in entries:
            if not path.name.endswith(suffix):
                continue
            transaction_identity = path.name.removesuffix(suffix)
            if not _valid_transaction_identity(transaction_identity):
                raise fail("DEPLOYMENT_STATE_INVALID", "Deployment state is invalid.")
            content, observed = read_regular_file(
                path,
                max_bytes=MAX_STATE_BYTES,
                code="DEPLOYMENT_STATE_INVALID",
            )
            if (
                observed.st_nlink != 1
                or stat.S_IMODE(observed.st_mode) != self._transaction_file_mode
                or (
                    expected_owner_uid is not None
                    and observed.st_uid != expected_owner_uid
                )
                or (
                    self._private_gid(expected_owner_gid) is not None
                    and observed.st_gid != self._private_gid(expected_owner_gid)
                )
            ):
                raise fail("DEPLOYMENT_STATE_INVALID", "Deployment state is invalid.")
            plan = StateRecoveryPlan.from_bytes(transaction_identity, content)
            if (
                plan.prior_state_identity == prior_state_identity
                and plan.candidate_state_identity == candidate_state_identity
            ):
                matches.append((path, plan))
        if len(matches) > 1:
            raise fail("DEPLOYMENT_STATE_INVALID", "Deployment state is invalid.")
        return None if not matches else matches[0]

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
            observed = self.layout.state_transactions.lstat()
        except OSError:
            raise fail(
                "DEPLOYMENT_STATE_INVALID", "Deployment state is invalid."
            ) from None
        if (
            not stat.S_ISDIR(observed.st_mode)
            or stat.S_ISLNK(observed.st_mode)
            or stat.S_IMODE(observed.st_mode) != self._private_root_mode
            or (
                expected_owner_uid is not None and observed.st_uid != expected_owner_uid
            )
            or (
                self._private_gid(expected_owner_gid) is not None
                and observed.st_gid != self._private_gid(expected_owner_gid)
            )
        ):
            raise fail("DEPLOYMENT_STATE_INVALID", "Deployment state is invalid.")
        entries: list[Path] = []
        try:
            for index, path in enumerate(self.layout.state_transactions.iterdir()):
                if index >= MAX_STATE_TRANSACTION_ENTRIES:
                    raise fail(
                        "DEPLOYMENT_STATE_INVALID", "Deployment state is invalid."
                    )
                entries.append(path)
        except DeploymentError:
            raise
        except OSError:
            raise fail(
                "DEPLOYMENT_STATE_INVALID", "Deployment state is invalid."
            ) from None
        values: list[str] = []
        for path in sorted(entries):
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
                or stat.S_IMODE(item.st_mode) != self._transaction_file_mode
                or (
                    expected_owner_uid is not None and item.st_uid != expected_owner_uid
                )
                or (
                    self._private_gid(expected_owner_gid) is not None
                    and item.st_gid != self._private_gid(expected_owner_gid)
                )
            ):
                raise fail("DEPLOYMENT_STATE_INVALID", "Deployment state is invalid.")
            values.append(path.name.removesuffix(".pending.json"))
        return tuple(values)

    def _pending_plans_locked(
        self,
        *,
        expected_owner_uid: int | None,
        expected_owner_gid: int | None,
    ) -> tuple[tuple[Path, StateRecoveryPlan], ...]:
        identities = self._pending_transactions_locked(
            expected_owner_uid=expected_owner_uid,
            expected_owner_gid=expected_owner_gid,
        )
        plans: list[tuple[Path, StateRecoveryPlan]] = []
        for transaction_identity in identities:
            if not _valid_transaction_identity(transaction_identity):
                raise fail("DEPLOYMENT_STATE_INVALID", "Deployment state is invalid.")
            path = self.layout.state_transactions / (
                f"{transaction_identity}.pending.json"
            )
            content, observed = read_regular_file(
                path,
                max_bytes=MAX_STATE_BYTES,
                code="DEPLOYMENT_STATE_INVALID",
            )
            if (
                observed.st_nlink != 1
                or stat.S_IMODE(observed.st_mode) != self._transaction_file_mode
                or (
                    expected_owner_uid is not None
                    and observed.st_uid != expected_owner_uid
                )
                or (
                    self._private_gid(expected_owner_gid) is not None
                    and observed.st_gid != self._private_gid(expected_owner_gid)
                )
            ):
                raise fail("DEPLOYMENT_STATE_INVALID", "Deployment state is invalid.")
            plans.append(
                (
                    path,
                    StateRecoveryPlan.from_bytes(transaction_identity, content),
                )
            )
        return tuple(plans)

    def _load_generation_locked(
        self,
        identity: str,
        *,
        expected_owner_uid: int | None,
        expected_owner_gid: int | None,
    ) -> DeploymentState:
        if not _valid_state_identity(identity):
            raise fail("DEPLOYMENT_STATE_INVALID", "Deployment state is invalid.")
        path = self.layout.state_generations / identity
        content, _observed = read_regular_file(
            path / STATE_FILENAME,
            max_bytes=MAX_STATE_BYTES,
            code="DEPLOYMENT_STATE_INVALID",
        )
        try:
            document = json.loads(content, object_pairs_hook=_unique_object)
        except (UnicodeDecodeError, json.JSONDecodeError, ValueError):
            raise fail(
                "DEPLOYMENT_STATE_INVALID", "Deployment state is invalid."
            ) from None
        state = DeploymentState.from_dict(document)
        if state.identity != identity or content != canonical_json_bytes(
            state.to_dict()
        ):
            raise fail("DEPLOYMENT_STATE_INVALID", "Deployment state is invalid.")
        environment = self._verify_platform_environment(
            path,
            state,
            expected_owner_uid=expected_owner_uid,
            expected_owner_gid=expected_owner_gid,
        )
        self._verify_generation(
            path,
            state,
            expected_owner_uid=expected_owner_uid,
            expected_owner_gid=expected_owner_gid,
            expected_environment=environment,
        )
        return state

    def _switch_pointer_locked(
        self,
        identity: str | None,
        *,
        transaction_identity: str,
        expected_owner_uid: int | None,
        expected_owner_gid: int | None,
    ) -> None:
        if identity is None:
            try:
                self.layout.current_state_link.unlink()
                fsync_directory(self.layout.current_state_link.parent)
            except FileNotFoundError:
                return
            except OSError:
                raise fail(
                    "DEPLOYMENT_ACTIVATION_FAILED",
                    "Deployment activation did not change the active generation.",
                    recoverable=True,
                ) from None
            return
        if not _valid_state_identity(identity):
            raise fail("DEPLOYMENT_STATE_INVALID", "Deployment state is invalid.")
        temporary_link = self.layout.current_state_link.parent / (
            f".current-{transaction_identity}"
        )
        try:
            os.symlink(f"generations/{identity}", temporary_link)
            if expected_owner_uid is not None:
                os.lchown(
                    temporary_link,
                    expected_owner_uid,
                    -1 if expected_owner_gid is None else expected_owner_gid,
                )
            os.replace(temporary_link, self.layout.current_state_link)
            fsync_directory(self.layout.current_state_link.parent)
        except OSError:
            raise fail(
                "DEPLOYMENT_ACTIVATION_FAILED",
                "Deployment activation did not change the active generation.",
                recoverable=True,
            ) from None

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
                or stat.S_IMODE(directory.st_mode) != self._generation_mode
                or stat.S_IMODE(observed.st_mode) != self._state_file_mode
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
            self._verify_platform_environment(
                path,
                state,
                expected_owner_uid=expected_owner_uid,
                expected_owner_gid=expected_owner_gid,
            )
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

    def _select_platform_environment(
        self,
        current: DeploymentState | None,
        candidate: DeploymentState,
        supplied: PlatformEnvironment | None,
        *,
        expected_owner_uid: int | None,
        expected_owner_gid: int | None,
    ) -> PlatformEnvironment | None:
        active = candidate.components[PLATFORM].active
        if active is None:
            if supplied is not None:
                raise fail(
                    "DEPLOYMENT_CONFIGURATION_INVALID",
                    "Deployment configuration is invalid.",
                )
            return None
        if supplied is not None:
            parsed = parse_platform_environment(
                self.layout, candidate, supplied.content
            )
            if parsed != supplied:
                raise fail(
                    "DEPLOYMENT_CONFIGURATION_INVALID",
                    "Deployment configuration is invalid.",
                )
            return parsed
        if current is None or any(
            current.components[component].active
            != candidate.components[component].active
            for component in COMPONENTS
        ):
            raise fail(
                "DEPLOYMENT_CONFIGURATION_REQUIRED",
                "Deployment configuration is required for activation.",
                recoverable=True,
            )
        previous = self._read_platform_environment(
            self.layout.state_generations / current.identity,
            current,
            expected_owner_uid=expected_owner_uid,
            expected_owner_gid=expected_owner_gid,
        )
        if previous is None:
            raise fail(
                "DEPLOYMENT_CONFIGURATION_REQUIRED",
                "Deployment configuration is required for activation.",
                recoverable=True,
            )
        return render_platform_environment(
            self.layout,
            candidate,
            api_contract_sha256=previous.api_contract_sha256,
        )

    def _read_platform_environment(
        self,
        generation: Path,
        state: DeploymentState,
        *,
        expected_owner_uid: int | None,
        expected_owner_gid: int | None,
    ) -> PlatformEnvironment | None:
        path = generation / PLATFORM_ENV_FILENAME
        if state.components[PLATFORM].active is None:
            if path.exists() or path.is_symlink():
                raise fail(
                    "DEPLOYMENT_CONFIGURATION_INVALID",
                    "Deployment configuration is invalid.",
                )
            return None
        if not self.verify_environment_content:
            try:
                observed = path.lstat()
            except OSError:
                raise fail(
                    "DEPLOYMENT_CONFIGURATION_INVALID",
                    "Deployment configuration is invalid.",
                ) from None
            if (
                not stat.S_ISREG(observed.st_mode)
                or stat.S_ISLNK(observed.st_mode)
                or observed.st_nlink != 1
                or not 0 < observed.st_size <= MAX_PLATFORM_ENV_BYTES
                or stat.S_IMODE(observed.st_mode) != self._environment_file_mode
                or (
                    expected_owner_uid is not None
                    and observed.st_uid != expected_owner_uid
                )
                or (
                    self._environment_gid(expected_owner_gid) is not None
                    and observed.st_gid != self._environment_gid(expected_owner_gid)
                )
            ):
                raise fail(
                    "DEPLOYMENT_CONFIGURATION_INVALID",
                    "Deployment configuration is invalid.",
                )
            return None
        content, observed = read_regular_file(
            path,
            max_bytes=MAX_PLATFORM_ENV_BYTES,
            code="DEPLOYMENT_CONFIGURATION_INVALID",
        )
        if (
            stat.S_IMODE(observed.st_mode) != self._environment_file_mode
            or (
                expected_owner_uid is not None and observed.st_uid != expected_owner_uid
            )
            or (
                self._environment_gid(expected_owner_gid) is not None
                and observed.st_gid != self._environment_gid(expected_owner_gid)
            )
        ):
            raise fail(
                "DEPLOYMENT_CONFIGURATION_INVALID",
                "Deployment configuration is invalid.",
            )
        return parse_platform_environment(self.layout, state, content)

    def _verify_platform_environment(
        self,
        generation: Path,
        state: DeploymentState,
        *,
        expected_owner_uid: int | None,
        expected_owner_gid: int | None,
    ) -> PlatformEnvironment | None:
        expected_names = {STATE_FILENAME}
        if state.components[PLATFORM].active is not None:
            expected_names.add(PLATFORM_ENV_FILENAME)
        try:
            observed_names = {entry.name for entry in generation.iterdir()}
        except OSError:
            raise fail(
                "DEPLOYMENT_STATE_INVALID", "Deployment state is invalid."
            ) from None
        if observed_names != expected_names:
            raise fail("DEPLOYMENT_STATE_INVALID", "Deployment state is invalid.")
        return self._read_platform_environment(
            generation,
            state,
            expected_owner_uid=expected_owner_uid,
            expected_owner_gid=expected_owner_gid,
        )

    def _verify_generation(
        self,
        path: Path,
        expected: DeploymentState,
        *,
        expected_owner_uid: int | None,
        expected_owner_gid: int | None,
        expected_environment: PlatformEnvironment | None,
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
            or stat.S_IMODE(directory.st_mode) != self._generation_mode
            or stat.S_IMODE(observed.st_mode) != self._state_file_mode
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
        observed_environment = self._verify_platform_environment(
            path,
            expected,
            expected_owner_uid=expected_owner_uid,
            expected_owner_gid=expected_owner_gid,
        )
        if observed_environment != expected_environment:
            raise fail("DEPLOYMENT_STATE_INVALID", "Deployment state is invalid.")

    @contextmanager
    def transaction(
        self,
        *,
        exclusive: bool,
        expected_owner_uid: int | None = None,
        expected_owner_gid: int | None = None,
    ) -> Iterator[StateTransaction]:
        if exclusive and not self.verify_environment_content:
            raise fail(
                "DEPLOYMENT_STATE_READ_ONLY",
                "Deployment state is read-only for this process.",
            )
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
            exclusive=exclusive,
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
        operator_state_root = self.layout.state_lock.parent
        state_root = self.layout.current_state_link.parent
        if not data_root.exists() and not data_root.is_symlink():
            create_directory(data_root, mode=0o755)
        self._require_trusted_directory(
            data_root,
            expected_owner_uid=expected_owner_uid,
            expected_owner_gid=expected_owner_gid,
            mode=0o755,
        )
        if not operator_root.exists() and not operator_root.is_symlink():
            try:
                operator_root.mkdir(mode=self._operator_root_mode)
                if expected_owner_uid is not None:
                    os.chown(
                        operator_root,
                        expected_owner_uid,
                        -1
                        if self._private_gid(expected_owner_gid) is None
                        else self._private_gid(expected_owner_gid),
                    )
            except OSError:
                raise fail(
                    "DEPLOYMENT_STORAGE_UNAVAILABLE",
                    "Deployment storage is unavailable.",
                ) from None
        self._require_trusted_directory(
            operator_root,
            expected_owner_uid=expected_owner_uid,
            expected_owner_gid=self._private_gid(expected_owner_gid),
            mode=self._operator_root_mode,
        )
        if not operator_state_root.exists() and not operator_state_root.is_symlink():
            try:
                operator_state_root.mkdir(mode=self._private_root_mode)
                if expected_owner_uid is not None:
                    os.chown(
                        operator_state_root,
                        expected_owner_uid,
                        -1
                        if self._private_gid(expected_owner_gid) is None
                        else self._private_gid(expected_owner_gid),
                    )
            except OSError:
                raise fail(
                    "DEPLOYMENT_STORAGE_UNAVAILABLE",
                    "Deployment storage is unavailable.",
                ) from None
        self._require_trusted_directory(
            operator_state_root,
            expected_owner_uid=expected_owner_uid,
            expected_owner_gid=self._private_gid(expected_owner_gid),
            mode=self._private_root_mode,
        )
        for path in (state_root, self.layout.state_generations):
            if not path.exists() and not path.is_symlink():
                try:
                    path.mkdir(mode=self._public_root_mode)
                    if expected_owner_uid is not None:
                        os.chown(
                            path,
                            expected_owner_uid,
                            -1 if expected_owner_gid is None else expected_owner_gid,
                        )
                except OSError:
                    raise fail(
                        "DEPLOYMENT_STORAGE_UNAVAILABLE",
                        "Deployment storage is unavailable.",
                    ) from None
            self._require_trusted_directory(
                path,
                expected_owner_uid=expected_owner_uid,
                expected_owner_gid=expected_owner_gid,
                mode=self._public_root_mode,
            )
        transactions = self.layout.state_transactions
        if not transactions.exists() and not transactions.is_symlink():
            try:
                transactions.mkdir(mode=self._private_root_mode)
                if expected_owner_uid is not None:
                    os.chown(
                        transactions,
                        expected_owner_uid,
                        -1
                        if self._private_gid(expected_owner_gid) is None
                        else self._private_gid(expected_owner_gid),
                    )
            except OSError:
                raise fail(
                    "DEPLOYMENT_STORAGE_UNAVAILABLE",
                    "Deployment storage is unavailable.",
                ) from None
        self._require_trusted_directory(
            transactions,
            expected_owner_uid=expected_owner_uid,
            expected_owner_gid=self._private_gid(expected_owner_gid),
            mode=self._private_root_mode,
        )

    def _validate_roots(
        self,
        *,
        expected_owner_uid: int | None,
        expected_owner_gid: int | None,
    ) -> None:
        for path, mode in (
            (self.layout.data_root, 0o755),
            (self.layout.current_state_link.parent, self._public_root_mode),
            (self.layout.state_generations, self._public_root_mode),
        ):
            self._require_trusted_directory(
                path,
                expected_owner_uid=expected_owner_uid,
                expected_owner_gid=expected_owner_gid,
                mode=mode,
            )
        self._require_trusted_directory(
            self.layout.data_root / "operator",
            expected_owner_uid=expected_owner_uid,
            expected_owner_gid=self._private_gid(expected_owner_gid),
            mode=self._operator_root_mode,
        )
        self._require_trusted_directory(
            self.layout.state_lock.parent,
            expected_owner_uid=expected_owner_uid,
            expected_owner_gid=self._private_gid(expected_owner_gid),
            mode=self._private_root_mode,
        )
        self._require_trusted_directory(
            self.layout.state_transactions,
            expected_owner_uid=expected_owner_uid,
            expected_owner_gid=self._private_gid(expected_owner_gid),
            mode=self._private_root_mode,
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
            descriptor = os.open(self.layout.state_lock, flags, self._lock_mode)
        except FileExistsError:
            return
        except OSError:
            raise fail(
                "DEPLOYMENT_STATE_WRITE_FAILED",
                "Deployment state could not be committed.",
                recoverable=True,
            ) from None
        try:
            if expected_owner_uid is not None:
                os.fchown(
                    descriptor,
                    expected_owner_uid,
                    -1
                    if self._private_gid(expected_owner_gid) is None
                    else self._private_gid(expected_owner_gid),
                )
            os.fchmod(descriptor, self._lock_mode)
            os.fsync(descriptor)
            observed = os.fstat(descriptor)
            if (
                expected_owner_uid is not None and observed.st_uid != expected_owner_uid
            ) or (
                self._private_gid(expected_owner_gid) is not None
                and observed.st_gid != self._private_gid(expected_owner_gid)
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
        exclusive: bool,
        expected_owner_uid: int | None,
        expected_owner_gid: int | None,
    ) -> int:
        access = os.O_RDWR if exclusive else os.O_RDONLY
        flags = access | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0)
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
            or stat.S_IMODE(observed.st_mode) != self._lock_mode
            or (
                expected_owner_uid is not None and observed.st_uid != expected_owner_uid
            )
            or (
                self._private_gid(expected_owner_gid) is not None
                and observed.st_gid != self._private_gid(expected_owner_gid)
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
    def _write_immutable(
        path: Path,
        content: bytes,
        *,
        mode: int,
        owner_uid: int | None,
        owner_gid: int | None,
    ) -> None:
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
                if owner_uid is not None:
                    os.fchown(
                        descriptor,
                        owner_uid,
                        -1 if owner_gid is None else owner_gid,
                    )
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


def _valid_state_identity(value: object) -> bool:
    return (
        isinstance(value, str)
        and len(value) == 71
        and value.startswith("sha256-")
        and all(character in "0123456789abcdef" for character in value[7:])
    )


def _valid_transaction_identity(value: object) -> bool:
    if not isinstance(value, str):
        return False
    generation, separator, nonce = value.partition("-")
    return (
        separator == "-"
        and len(generation) == 8
        and generation.isascii()
        and generation.isdigit()
        and len(nonce) == 32
        and all(character in "0123456789abcdef" for character in nonce)
    )


__all__ = [
    "FaultInjector",
    "PLATFORM_ENV_FILENAME",
    "PLATFORM_ENV_KEYS",
    "PlatformEnvironment",
    "STATE_FILENAME",
    "STATE_RECOVERY_CHOICES",
    "StateRecoveryPlan",
    "StateStore",
    "parse_platform_environment",
    "render_platform_environment",
]
