"""Durable successful-validation snapshots and snapshot-only run creation."""

from __future__ import annotations

from collections.abc import Callable, Mapping
from datetime import datetime, timedelta, timezone
from uuid import uuid4

from encode_pipeline.platform.adapters import VALIDATION_CAPABILITY, WorkflowInputs
from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.platform.results import Issue, Result
from encode_pipeline.platform.snapshots import (
    PAYLOAD_DIGEST_SCHEME,
    VALIDATION_EVIDENCE_OUTCOME,
    ValidatedInputSnapshot,
    ValidatedSnapshotRunCreation,
    build_workflow_inputs_digest,
    canonical_workflow_inputs_json,
)
from encode_pipeline.services.run_repositories import (
    RunRepository,
    ValidatedSnapshotBuildMismatchError as RepositoryBuildMismatchError,
    ValidatedSnapshotExpiredError as RepositoryExpiredError,
    ValidatedSnapshotReplayConflictError as RepositoryReplayConflictError,
)
from encode_pipeline.services.runs import RunService
from encode_pipeline.services.validation import ValidationService
from encode_pipeline.services.workflow_builds import WorkflowBuildIdentityProvider


DEFAULT_SNAPSHOT_TTL = timedelta(minutes=30)


class ValidatedSnapshotNotFoundError(RuntimeError):
    """Snapshot is absent or does not belong to the requested workflow."""


class ValidatedSnapshotExpiredError(RuntimeError):
    """Snapshot expired before its first atomic use."""


class ValidatedSnapshotStaleError(RuntimeError):
    """Snapshot refers to a workflow build that is no longer current."""


class ValidatedSnapshotReplayConflictError(RuntimeError):
    """Snapshot was already consumed with different run metadata."""


class ValidatedSnapshotDataInvalidError(RuntimeError):
    """Persisted snapshot or linked-run evidence failed validation."""


class ValidatedSnapshotBuildUnavailableError(RuntimeError):
    """Current workflow source could not be fingerprinted safely."""


def _now_utc(clock: Callable[[], datetime]) -> datetime:
    value = clock()
    if (
        not isinstance(value, datetime)
        or value.tzinfo is None
        or value.utcoffset() is None
    ):
        raise ValueError("clock must return a timezone-aware datetime")
    return value.astimezone(timezone.utc)


def _build_unavailable_issue() -> Issue:
    return Issue(
        code="VALIDATION_WORKFLOW_BUILD_UNAVAILABLE",
        message="Workflow source identity could not be confirmed for validation.",
        source="workflow_build_identity",
        path="workflow",
    )


def _schema_unavailable_issue() -> Issue:
    return Issue(
        code="VALIDATION_WORKFLOW_SCHEMA_UNAVAILABLE",
        message="Workflow input contract could not be confirmed for validation.",
        source="workflow_schema",
        path="workflow",
    )


class ValidatedInputService:
    """Validate adapter inputs and persist complete immutable success evidence."""

    def __init__(
        self,
        *,
        registry: WorkflowRegistry,
        validation_service: ValidationService,
        build_identity_provider: WorkflowBuildIdentityProvider,
        repository: RunRepository,
        snapshot_id_factory: Callable[[], str] | None = None,
        clock: Callable[[], datetime] | None = None,
        snapshot_ttl: timedelta = DEFAULT_SNAPSHOT_TTL,
    ) -> None:
        if not isinstance(registry, WorkflowRegistry):
            raise ValueError("registry must be a WorkflowRegistry")
        if not isinstance(validation_service, ValidationService):
            raise ValueError("validation_service must be a ValidationService")
        if not isinstance(snapshot_ttl, timedelta) or snapshot_ttl <= timedelta(0):
            raise ValueError("snapshot_ttl must be positive")
        self._registry = registry
        self._validation_service = validation_service
        self._build_identity_provider = build_identity_provider
        self._repository = repository
        self._snapshot_id_factory = snapshot_id_factory or (
            lambda: f"vsnap_{uuid4().hex}"
        )
        self._clock = clock or (lambda: datetime.now(timezone.utc))
        self._snapshot_ttl = snapshot_ttl

    def validate(
        self,
        workflow_id: str,
        inputs: WorkflowInputs,
    ) -> Result[ValidatedInputSnapshot]:
        """Return a durable snapshot only after stable successful validation."""
        try:
            adapter = self._registry.get(workflow_id)
        except KeyError:
            result = self._validation_service.validate(workflow_id, inputs)
            return Result.failure(result.issues)
        if VALIDATION_CAPABILITY not in adapter.capabilities.supports:
            result = self._validation_service.validate(workflow_id, inputs)
            return Result.failure(result.issues)

        before_result = self._build_identity_provider.capture(workflow_id)
        if before_result.is_failure or before_result.value is None:
            return Result.failure([_build_unavailable_issue()])

        try:
            schema = adapter.schema()
        except Exception:
            return Result.failure([_schema_unavailable_issue()])

        validation_result = self._validation_service.validate(workflow_id, inputs)
        if validation_result.is_failure:
            return Result.failure(validation_result.issues)

        after_result = self._build_identity_provider.capture(workflow_id)
        if after_result.is_failure or after_result.value is None:
            return Result.failure([_build_unavailable_issue()])
        if not before_result.value.matches(after_result.value):
            return Result.failure(
                [
                    Issue(
                        code="VALIDATION_WORKFLOW_BUILD_CHANGED",
                        message="Workflow source changed during validation; validate again.",
                        source="workflow_build_identity",
                        path="workflow",
                    )
                ]
            )

        try:
            canonical_payload = canonical_workflow_inputs_json(inputs)
            now = _now_utc(self._clock)
            snapshot = ValidatedInputSnapshot(
                snapshot_id=self._snapshot_id_factory(),
                workflow_id=adapter.metadata.workflow_id,
                adapter_version=adapter.metadata.version,
                schema_version=schema.schema_version,
                schema_dialect=schema.schema_dialect,
                workflow_build_identity=after_result.value,
                canonical_payload=canonical_payload,
                payload_digest_scheme=PAYLOAD_DIGEST_SCHEME,
                payload_digest=build_workflow_inputs_digest(canonical_payload),
                validation_outcome=VALIDATION_EVIDENCE_OUTCOME,
                validation_issue_codes=tuple(
                    issue.code for issue in validation_result.issues
                ),
                validated_at=now,
                expires_at=now + self._snapshot_ttl,
            )
            persisted = self._repository.create_validated_input_snapshot(snapshot)
        except Exception:
            return Result.failure(
                [
                    Issue(
                        code="VALIDATED_SNAPSHOT_PERSISTENCE_FAILED",
                        message="Validated inputs could not be saved. Validate again.",
                        source="persistence",
                        path="inputs",
                    )
                ]
            )
        return Result.success(persisted, issues=validation_result.issues)


class ValidatedRunCreationService:
    """Create one durable run from a server-owned successful validation."""

    def __init__(
        self,
        *,
        run_service: RunService,
        build_identity_provider: WorkflowBuildIdentityProvider,
        clock: Callable[[], datetime] | None = None,
    ) -> None:
        if not isinstance(run_service, RunService):
            raise ValueError("run_service must be a RunService")
        self._run_service = run_service
        self._build_identity_provider = build_identity_provider
        self._clock = clock or (lambda: datetime.now(timezone.utc))

    def create_run(
        self,
        workflow_id: str,
        snapshot_id: str,
        *,
        tags: Mapping[str, str] | None = None,
    ) -> ValidatedSnapshotRunCreation:
        """Consume a snapshot once; identical retries return its canonical run."""
        try:
            snapshot = self._run_service.get_validated_input_snapshot(snapshot_id)
        except KeyError:
            raise ValidatedSnapshotNotFoundError from None
        except (TypeError, ValueError):
            raise ValidatedSnapshotDataInvalidError from None
        if snapshot.workflow_id != workflow_id:
            raise ValidatedSnapshotNotFoundError

        now = _now_utc(self._clock)
        if snapshot.consumed_run_id is not None:
            expected_identity = snapshot.workflow_build_identity
        else:
            if now >= snapshot.expires_at:
                raise ValidatedSnapshotExpiredError
            identity_result = self._build_identity_provider.capture(workflow_id)
            if identity_result.is_failure or identity_result.value is None:
                raise ValidatedSnapshotBuildUnavailableError
            if not snapshot.workflow_build_identity.matches(identity_result.value):
                raise ValidatedSnapshotStaleError
            expected_identity = identity_result.value

        try:
            return self._run_service.create_run_from_validated_snapshot(
                workflow_id,
                snapshot_id,
                expected_build_identity=expected_identity,
                consumed_at=now,
                tags=tags,
            )
        except KeyError:
            raise ValidatedSnapshotNotFoundError from None
        except RepositoryExpiredError:
            raise ValidatedSnapshotExpiredError from None
        except RepositoryBuildMismatchError:
            raise ValidatedSnapshotStaleError from None
        except RepositoryReplayConflictError:
            raise ValidatedSnapshotReplayConflictError from None
        except (TypeError, ValueError):
            raise ValidatedSnapshotDataInvalidError from None
