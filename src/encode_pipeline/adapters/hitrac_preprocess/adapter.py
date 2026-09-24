"""Hi-TrAC authoring and private execution adapter; public execution stays closed."""

from __future__ import annotations

from pathlib import Path

from encode_pipeline import __version__
from encode_pipeline.platform.adapters import (
    WorkflowAvailability,
    WorkflowCapabilities,
    WorkflowMetadata,
    WorkflowUpstreamIdentity,
)
from encode_pipeline.platform.reference_profiles import (
    AdapterReferenceBindingIdentity,
    BoundWorkflowReference,
)
from encode_pipeline.platform.results import Result

from .admission import load_reference_binding
from .authoring import build_hitrac_authoring_schema
from .execution import (
    WORKFLOW_ID,
    ExecutionBinding,
    RuntimeAdmission,
    bind_inputs,
    build_command,
    capture_identity,
    failure,
    plan_workspace,
)
from .validation import validate_hitrac_inputs


class HiTracPreprocessAdapter:
    metadata = WorkflowMetadata(
        workflow_id=WORKFLOW_ID,
        name="Hi-TrAC preprocessing",
        version=__version__,
        description="Pinned tracPre2 preprocessing. Public execution is not configured until artifact/QC publication is implemented and qualified.",
        engines=("hitrac-qualification",),
        tags=("hitrac", "preprocessing"),
    )
    upstream_identity = WorkflowUpstreamIdentity(
        name="cLoops2/tracPre2",
        version="git-v0.0.5",
        revision="de6cc732fa00b408551b9f4272933640c08447f1",
    )
    capabilities = WorkflowCapabilities(
        supports=("validation", "input_authoring", "workspace_plan", "command")
    )

    def __init__(
        self,
        *,
        runtime: RuntimeAdmission | None = None,
        binding: ExecutionBinding | None = None,
    ):
        if runtime is not None and not isinstance(runtime, RuntimeAdmission):
            raise ValueError("invalid runtime admission")
        if binding is not None and (
            not isinstance(binding, ExecutionBinding) or binding.runtime != runtime
        ):
            raise ValueError("invalid execution binding")
        self._runtime = runtime
        self._binding = binding

    def schema(self):
        return build_hitrac_authoring_schema()

    def validate(self, inputs):
        return validate_hitrac_inputs(inputs)

    def execution_availability(self):
        # Admission does not grant public execution before H4. No bypass option.
        return WorkflowAvailability(
            execution="not_configured", reason_code="WORKFLOW_EXECUTION_NOT_CONFIGURED"
        )

    def preview_dag(self, inputs):
        return failure("WORKFLOW_CAPABILITY_UNSUPPORTED", "preview_dag")

    def extract_artifacts(self, inputs, workspace):
        # Required WorkflowAdapter protocol method, deliberately no capability.
        return failure("WORKFLOW_CAPABILITY_UNSUPPORTED", "artifact_extract")

    def verify_reference_profile_binding(self, payload):
        try:
            if (
                not isinstance(payload, dict)
                or set(payload) != {"schema_version", "binding", "sha256"}
                or payload["schema_version"] != "hitrac-reference-profile-v1"
            ):
                raise ValueError
            reference = load_reference_binding(
                Path(payload["binding"]), payload["sha256"]
            )
            return Result.success(
                AdapterReferenceBindingIdentity(
                    workflow_id=WORKFLOW_ID,
                    contract_version="hitrac-reference-profile-v1",
                    identity_sha256=reference.binding_sha256,
                )
            )
        except (OSError, ValueError, TypeError, KeyError):
            return failure(
                "HITRAC_REFERENCE_BINDING_INVALID", "reference_profile_revision_id"
            )

    def bind_reference_profile(self, inputs, payload):
        checked = self.validate(inputs)
        if checked.is_failure:
            return checked
        if self._runtime is None:
            return failure(
                "HITRAC_REFERENCE_BINDING_UNAVAILABLE", "reference_profile_revision_id"
            )
        verified = self.verify_reference_profile_binding(payload)
        if verified.is_failure:
            return verified
        try:
            binding = bind_inputs(self._runtime, checked.value, payload)
            adapter = type(self)(runtime=self._runtime, binding=binding)
            return Result.success(
                BoundWorkflowReference(
                    inputs=checked.value, adapter=adapter, identity=verified.value
                )
            )
        except (OSError, ValueError, TypeError, KeyError):
            return failure("HITRAC_INPUT_BINDING_INVALID", "samples")

    def capture_build_identity(self):
        if self._binding is None:
            return failure("HITRAC_EXECUTION_BINDING_REQUIRED")
        try:
            return Result.success(
                capture_identity(self._binding, self.metadata.version)
            )
        except (OSError, ValueError, TypeError, KeyError):
            return failure("HITRAC_EXECUTION_IDENTITY_CHANGED")

    def plan_workspace(self, inputs, workspace):
        if self._binding is None:
            return failure("HITRAC_EXECUTION_BINDING_REQUIRED")
        return plan_workspace(inputs, workspace, self._binding)

    def build_command(self, plan, workspace):
        if self._binding is None:
            return failure("HITRAC_EXECUTION_BINDING_REQUIRED")
        return build_command(plan, workspace, self._binding)
