"""Durable successful-validation snapshots and snapshot-only run creation."""

from __future__ import annotations

from collections.abc import Callable, Mapping
from datetime import datetime, timedelta, timezone
from uuid import uuid4

from encode_pipeline.platform.adapters import (
    VALIDATION_CAPABILITY,
    InputUseDeclaringAdapter,
    WorkflowInputs,
)
from encode_pipeline.platform.data_registry import (
    ProjectSampleBinding,
    ProjectSampleSelection,
)
from encode_pipeline.platform.input_registry import (
    AdapterInputUseContract,
    InputFileRevisionSelection,
    InputProvenanceMode,
    InputUseBindingEnvelope,
    InputUseBindingPlan,
    build_input_use_binding_plan,
)
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
    InputBindingSelectionError,
    ProjectSampleSelectionError,
    RunRepository,
    ValidatedSnapshotBuildMismatchError as RepositoryBuildMismatchError,
    ValidatedSnapshotExpiredError as RepositoryExpiredError,
    ValidatedSnapshotReplayConflictError as RepositoryReplayConflictError,
)
from encode_pipeline.services.runs import RunService
from encode_pipeline.services.validation import ValidationService
from encode_pipeline.services.workflow_builds import WorkflowBuildIdentityProvider
from encode_pipeline.services.workflow_info import resolve_workflow_availability


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


class ValidatedSnapshotExecutionUnavailableError(RuntimeError):
    """Current deployment cannot admit execution for a new run."""


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


def _input_use_contract_unavailable_issue() -> Issue:
    return Issue(
        code="INPUT_USE_CAPABILITY_UNAVAILABLE",
        message="The exact workflow input-use contract could not be confirmed.",
        source="input_registry",
        path="input_selections",
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
        *,
        project_sample_selection: ProjectSampleSelection | None = None,
        input_file_revision_selections: tuple[InputFileRevisionSelection, ...] = (),
    ) -> Result[ValidatedInputSnapshot | None]:
        """Return a durable snapshot only after stable successful validation."""
        input_file_revision_selections = tuple(input_file_revision_selections)
        if any(
            not isinstance(selection, InputFileRevisionSelection)
            for selection in input_file_revision_selections
        ):
            raise ValueError(
                "input_file_revision_selections must contain exact selections"
            )
        try:
            adapter = self._registry.get(workflow_id)
        except KeyError:
            result = self._validation_service.validate(workflow_id, inputs)
            return Result.failure(result.issues)
        if VALIDATION_CAPABILITY not in adapter.capabilities.supports:
            result = self._validation_service.validate(workflow_id, inputs)
            return Result.failure(result.issues)
        if input_file_revision_selections and not isinstance(
            adapter,
            InputUseDeclaringAdapter,
        ):
            return Result.failure(
                [
                    Issue(
                        code="INPUT_USE_CAPABILITY_UNAVAILABLE",
                        message=(
                            "The exact workflow input use does not support "
                            "managed revision selection."
                        ),
                        source="input_registry",
                        path="input_selections",
                    )
                ]
            )
        if input_file_revision_selections and project_sample_selection is None:
            return Result.failure(
                [
                    Issue(
                        code="INPUT_BINDING_SELECTION_INVALID",
                        message=(
                            "Managed input revisions require an exact Project "
                            "and SampleRevision selection."
                        ),
                        source="input_registry",
                        path="project_id",
                    )
                ]
            )

        availability = resolve_workflow_availability(adapter)
        if input_file_revision_selections and availability.execution != "available":
            return Result.failure(
                [
                    Issue(
                        code="INPUT_USE_CAPABILITY_UNAVAILABLE",
                        message=(
                            "Managed provenance is not qualified for this exact "
                            "workflow input use."
                        ),
                        source="input_registry",
                        path="input_selections",
                    )
                ]
            )
        if availability.execution != "available":
            try:
                adapter.schema()
            except Exception:
                return Result.failure([_schema_unavailable_issue()])
            result = self._validation_service.validate(workflow_id, inputs)
            if result.is_failure:
                return Result.failure(result.issues)
            return Result.success(None, issues=result.issues)

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

        input_use_binding_plan: InputUseBindingPlan | None = None
        if project_sample_selection is not None and isinstance(
            adapter,
            InputUseDeclaringAdapter,
        ):
            try:
                declaration_result = adapter.declare_input_uses(
                    inputs,
                    validation_result.value,
                )
            except Exception:
                return Result.failure([_input_use_contract_unavailable_issue()])
            if (
                not isinstance(declaration_result, Result)
                or declaration_result.is_failure
                or not isinstance(
                    declaration_result.value,
                    AdapterInputUseContract,
                )
            ):
                return Result.failure([_input_use_contract_unavailable_issue()])
            try:
                input_use_binding_plan = build_input_use_binding_plan(
                    project_id=project_sample_selection.project_id,
                    workflow_id=workflow_id,
                    contract=declaration_result.value,
                    selections=input_file_revision_selections,
                )
            except (TypeError, ValueError):
                return Result.failure(
                    [
                        Issue(
                            code="INPUT_BINDING_SELECTION_INVALID",
                            message=(
                                "The selected input revisions do not match the "
                                "exact workflow input-use contract."
                            ),
                            source="input_registry",
                            path="input_selections",
                        )
                    ]
                )
            if any(
                input_use.provenance_mode is InputProvenanceMode.MANAGED_REVISION_V1
                for input_use in input_use_binding_plan.input_uses
            ):
                return Result.failure(
                    [
                        Issue(
                            code="INPUT_USE_CAPABILITY_UNAVAILABLE",
                            message=(
                                "Managed provenance is not qualified for this "
                                "exact workflow input use."
                            ),
                            source="input_registry",
                            path="input_selections",
                        )
                    ]
                )

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
            if input_use_binding_plan is None:
                persisted = self._repository.create_validated_input_snapshot(
                    snapshot,
                    project_sample_selection=project_sample_selection,
                )
            else:
                persisted = self._repository.create_validated_input_snapshot(
                    snapshot,
                    project_sample_selection=project_sample_selection,
                    input_use_binding_plan=input_use_binding_plan,
                )
        except ProjectSampleSelectionError:
            return Result.failure(
                [
                    Issue(
                        code="DATA_BINDING_SELECTION_INVALID",
                        message=(
                            "The selected Project or Sample revisions are not "
                            "eligible for a new validation snapshot."
                        ),
                        source="data_registry",
                        path="project_id",
                    )
                ]
            )
        except InputBindingSelectionError:
            return Result.failure(
                [
                    Issue(
                        code="INPUT_BINDING_SELECTION_INVALID",
                        message=(
                            "The selected input revisions are not eligible for "
                            "this validation snapshot."
                        ),
                        source="input_registry",
                        path="input_selections",
                    )
                ]
            )
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

    def get_validated_input_binding(
        self,
        snapshot_id: str,
    ) -> ProjectSampleBinding:
        """Return the safe immutable binding frozen with one snapshot."""
        return self._repository.get_validated_input_binding(snapshot_id)

    def get_validated_input_use_binding(
        self,
        snapshot_id: str,
    ) -> InputUseBindingEnvelope:
        """Return the safe immutable input-use evidence frozen with a snapshot."""
        return self._repository.get_validated_input_use_binding(snapshot_id)


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
            try:
                adapter = self._run_service.registry.get(workflow_id)
            except (KeyError, ValueError):
                raise ValidatedSnapshotExecutionUnavailableError from None
            if resolve_workflow_availability(adapter).execution != "available":
                raise ValidatedSnapshotExecutionUnavailableError
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
