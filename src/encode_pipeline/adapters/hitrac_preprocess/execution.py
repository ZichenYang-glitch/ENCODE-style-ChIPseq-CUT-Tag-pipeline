"""Server-owned Hi-TrAC binding and private qualification execution plans."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import sys

from encode_pipeline.platform.adapters import CommandSpec, WorkflowInputs, WorkspacePlan
from encode_pipeline.platform.builds import WorkflowBuildIdentity
from encode_pipeline.platform.results import Issue, Result

from .admission import (
    AdmissionError,
    _regular,
    load_reference_binding,
    load_runtime_binding,
    sha256_file,
    validate_fastq_pair,
)
from .qualification import implementation_identity
from .validation import validate_hitrac_inputs

WORKFLOW_ID = "hitrac-preprocess"
REQUEST = "hitrac-request.json"
ATTEMPT = "hitrac-attempt"
ROOT = Path(__file__).resolve().parents[4]


def encoded(value):
    return json.dumps(
        value, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode()


def sha(value):
    return hashlib.sha256(encoded(value)).hexdigest()


def failure(code, path="workflow"):
    return Result.failure(
        [
            Issue(
                code=code,
                message="Hi-TrAC qualification prerequisites could not be confirmed.",
                path=path,
                source=WORKFLOW_ID,
                hint="Check the approved runtime, reference selection and inputs; validate again.",
            )
        ]
    )


def source_identity():
    return implementation_identity()


@dataclass(frozen=True)
class RuntimeAdmission:
    """Operator-created authority, never constructed from WorkflowInputs."""

    binding: Path
    binding_sha256: str
    python: Path
    python_sha256: str
    implementation_sha256: str
    timeout: float = 300

    def verify(self):
        if sha256_file(_regular(self.binding)) != self.binding_sha256:
            raise AdmissionError("runtime_drift")
        runtime = load_runtime_binding(self.binding)
        if (
            sha256_file(self.python) != self.python_sha256
            or self.python != Path(sys.executable).absolute()
        ):
            raise AdmissionError("interpreter_drift")
        if source_identity()["sha256"] != self.implementation_sha256:
            raise AdmissionError("implementation_drift")
        return runtime


def admit_runtime(binding: Path, expected_sha256: str, *, timeout=300):
    """Verify the operator's exact binding in this normally installed checkout."""
    import math

    if re.fullmatch(r"[0-9a-f]{64}", expected_sha256 or "") is None:
        raise AdmissionError("runtime_identity_invalid")
    if (
        isinstance(timeout, bool)
        or not isinstance(timeout, (int, float))
        or not math.isfinite(timeout)
        or timeout <= 0
    ):
        raise AdmissionError("timeout_invalid")
    if sha256_file(_regular(binding)) != expected_sha256:
        raise AdmissionError("runtime_drift")
    load_runtime_binding(binding)
    python = Path(sys.executable).absolute()
    return RuntimeAdmission(
        binding,
        expected_sha256,
        python,
        sha256_file(python),
        source_identity()["sha256"],
        timeout,
    )


def freeze_inputs(inputs: WorkflowInputs):
    result = validate_hitrac_inputs(inputs)
    if result.is_failure:
        raise AdmissionError("inputs_invalid")
    normalized = result.value
    samples = {}
    for index, row in enumerate(normalized.samples, 1):
        token = f"s{index:06d}"
        paths = [_regular(Path(row[key])) for key in ("fastq_1", "fastq_2")]
        before = [p.stat() for p in paths]
        hashes = [sha256_file(p) for p in paths]
        count = validate_fastq_pair(*paths)
        after = [p.stat() for p in paths]
        if [
            (s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns, s.st_ctime_ns)
            for s in before
        ] != [
            (s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns, s.st_ctime_ns) for s in after
        ]:
            raise AdmissionError("input_drift")
        samples[token] = {
            "display_id": row["sample_id"],
            "raw_pairs": count,
            "sha256": dict(zip(("r1", "r2"), hashes)),
        }
    return normalized, samples


@dataclass(frozen=True)
class ExecutionBinding:
    runtime: RuntimeAdmission
    reference: Path
    reference_sha256: str
    inputs_json: bytes
    samples_json: bytes

    def verify(self):
        runtime = self.runtime.verify()
        reference = load_reference_binding(self.reference, self.reference_sha256)
        inputs = WorkflowInputs(**json.loads(self.inputs_json))
        normalized, samples = freeze_inputs(inputs)
        if encoded(samples) != self.samples_json:
            raise AdmissionError("input_drift")
        return runtime, reference, normalized, samples

    def identity(self):
        runtime, reference, inputs, samples = self.verify()
        return {
            "scheme": "sha256-hitrac-execution-binding-v1",
            "implementation": self.runtime.implementation_sha256,
            "python": self.runtime.python_sha256,
            "runtime": runtime.binding_sha256,
            "runtime_lock": runtime.lock_sha256,
            "reference": reference.binding_sha256,
            "inputs": sha({"payload": inputs.to_dict(), "samples": samples}),
            "timeout": self.runtime.timeout,
        }


def bind_inputs(runtime, inputs, payload):
    if (
        not isinstance(payload, dict)
        or set(payload) != {"schema_version", "binding", "sha256"}
        or payload["schema_version"] != "hitrac-reference-profile-v1"
    ):
        raise AdmissionError("reference_binding_invalid")
    runtime.verify()
    path = Path(payload["binding"])
    load_reference_binding(path, payload["sha256"])
    normalized, samples = freeze_inputs(inputs)
    return ExecutionBinding(
        runtime,
        path,
        payload["sha256"],
        encoded(normalized.to_dict()),
        encoded(samples),
    )


def capture_identity(binding: ExecutionBinding, version):
    return WorkflowBuildIdentity(
        workflow_id=WORKFLOW_ID,
        adapter_version=version,
        scheme="sha256-hitrac-execution-binding-v1",
        logical_entrypoint="scripts/run_hitrac_preprocess.py",
        digest=sha(binding.identity()),
        captured_at=datetime.now(timezone.utc),
    )


def workspace_path(workspace):
    path = Path(workspace)
    if (
        not path.is_absolute()
        or path.resolve() != path
        or not re.fullmatch(r"/[A-Za-z0-9_./-]+", str(path))
    ):
        raise AdmissionError("workspace_invalid")
    if (path / ATTEMPT).exists() or (path / ATTEMPT).is_symlink():
        raise AdmissionError("attempt_exists")
    return path


def plan_workspace(inputs, workspace, binding):
    try:
        path = workspace_path(workspace)
        _, _, normalized, samples = binding.verify()
        checked = validate_hitrac_inputs(inputs)
        if (
            checked.is_failure
            or encoded(checked.value.to_dict()) != binding.inputs_json
        ):
            raise AdmissionError("inputs_binding_mismatch")
        payload = {
            "schema_version": "hitrac-private-plan-v1",
            "workspace": str(path),
            "runtime_binding": str(binding.runtime.binding),
            "runtime_sha256": binding.runtime.binding_sha256,
            "reference_binding": str(binding.reference),
            "reference_sha256": binding.reference_sha256,
            "inputs": normalized.to_dict(),
            "samples": samples,
            "identity": binding.identity(),
            "timeout": binding.runtime.timeout,
        }
        return Result.success(WorkspacePlan(files=((REQUEST, encoded(payload)),)))
    except (OSError, ValueError, TypeError, KeyError):
        return failure("HITRAC_PLAN_REJECTED")


def build_command(plan, workspace, binding):
    try:
        path = workspace_path(workspace)
        expected = plan_workspace(
            WorkflowInputs(**json.loads(binding.inputs_json)), path, binding
        )
        if expected.is_failure or plan != expected.value:
            raise AdmissionError("plan_binding_mismatch")
        data = dict(plan.files)[REQUEST]
        runtime = binding.runtime
        script = ROOT / "scripts/run_hitrac_preprocess.py"
        argv = (
            str(runtime.python),
            "-I",
            "-S",
            "-B",
            str(script),
            "--request",
            str(path / REQUEST),
            "--sha256",
            hashlib.sha256(data).hexdigest(),
        )
        return Result.success(
            CommandSpec(
                argv=argv,
                cwd=str(path),
                env={"PYTHONDONTWRITEBYTECODE": "1"},
                redaction_values=(
                    str(path),
                    str(runtime.python),
                    str(ROOT),
                    str(binding.reference),
                    str(runtime.binding),
                ),
            )
        )
    except (OSError, ValueError, TypeError, KeyError):
        return failure("HITRAC_COMMAND_REJECTED")
