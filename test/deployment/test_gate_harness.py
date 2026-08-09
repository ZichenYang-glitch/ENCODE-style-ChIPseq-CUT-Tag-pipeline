from __future__ import annotations

import hashlib
import importlib.util
import inspect
import json
import os
from pathlib import Path
import stat
import sys

import pytest

from encode_pipeline.deployment.canonical import canonical_identity
from encode_pipeline.deployment.errors import DeploymentError, fail
from encode_pipeline.deployment.gate import (
    CLEANUP_PLAN_IDENTITY_SCHEME,
    GATE_OBSERVATION_IDENTITY_SCHEME,
    GATE_STAGES,
    PROCESS_NAMES,
    REQUIRED_CACHE_KINDS,
    RUNTIME_COMPONENTS,
    CacheEvidence,
    CleanupPlan,
    CleanupTarget,
    DeploymentGateRun,
    FilesystemGateObserver,
    GateObservation,
    GatePolicy,
    GateProcessIdentity,
    GateRequest,
    GateStageVerification,
    ResourceEvidence,
    prepare_gate_request,
    verify_cleanup_script,
)
from encode_pipeline.deployment import gate as gate_contract
from encode_pipeline.deployment.models import BULK_RNASEQ_RUNTIME, ENCODE_RUNTIME


ROOT = Path(__file__).resolve().parents[2]
PREPARE_SCRIPT = ROOT / "scripts" / "prepare_deployment_gate.py"
TASK_IDENTITY = f"task-{'1' * 32}"
HEAD_SHA = "2" * 40
RELEASE_IDENTITY = f"sha256-{'3' * 64}"
ENCODE_IDENTITY = f"sha256-{'4' * 64}"
BULK_IDENTITY = f"sha256-{'5' * 64}"
ENVIRONMENT_IDENTITY = f"sha256-{'6' * 64}"
VERIFIER_IDENTITY = f"sha256-{'7' * 64}"
EVIDENCE_IDENTITY = f"sha256-{'8' * 64}"
DISPATCH_IDENTITY = f"sha256-{'9' * 64}"
APPROVAL_IDENTITY = f"sha256-{'a' * 64}"


def _load_prepare_module():
    spec = importlib.util.spec_from_file_location(
        "test_prepare_deployment_gate", PREPARE_SCRIPT
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


prepare_module = _load_prepare_module()


def _policy(tmp_path: Path) -> GatePolicy:
    return GatePolicy.isolated(tmp_path / "fixed", uid=os.getuid(), gid=os.getgid())


def _mkdir(path: Path) -> Path:
    path.mkdir(parents=True, mode=0o700, exist_ok=True)
    path.chmod(0o700)
    return path


def _paths(policy: GatePolicy) -> None:
    _mkdir(policy.gate_root)
    _mkdir(policy.gate_root / "tasks")
    task_root = _mkdir(policy.task_root(TASK_IDENTITY))
    for path in policy.delete_paths(TASK_IDENTITY).values():
        directory = path.parent if path.name == "docker.sock" else path
        _mkdir(directory)
    _mkdir(task_root / "runner" / "prepared")
    for path in policy.retain_paths().values():
        _mkdir(path)
    _mkdir(policy.cleanup_executor.parent)
    policy.cleanup_executor.write_bytes(b"fixed root-owned cleanup executor\n")
    policy.cleanup_executor.chmod(0o555)
    _mkdir(policy.operator_root)
    backend_content = b"fixed root-owned cleanup backend\n"
    closure_document = {
        "schema_version": "helixweave-operator-closure-v1",
        "files": [
            {
                "mode": 0o555,
                "path": "bin/helixweave-gate-cleanup",
                "sha256": hashlib.sha256(backend_content).hexdigest(),
                "size_bytes": len(backend_content),
            }
        ],
    }
    closure_identity = (
        "sha256-"
        + hashlib.sha256(
            json.dumps(
                closure_document,
                allow_nan=False,
                separators=(",", ":"),
                sort_keys=True,
            ).encode("utf-8")
        ).hexdigest()
    )
    closure = _mkdir(policy.operator_root / closure_identity)
    backend_directory = _mkdir(closure / "bin")
    backend = backend_directory / "helixweave-gate-cleanup"
    backend.write_bytes(backend_content)
    backend.chmod(0o555)
    (closure / "closure.json").write_text(
        json.dumps(
            {**closure_document, "identity": closure_identity},
            allow_nan=False,
            separators=(",", ":"),
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    (closure / "closure.json").chmod(0o444)
    backend_directory.chmod(0o555)
    closure.chmod(0o555)
    (policy.operator_root / "current").symlink_to(closure_identity)


def _process(name: str) -> GateProcessIdentity:
    index = PROCESS_NAMES.index(name) + 1
    return GateProcessIdentity.create(
        name=name,
        task_identity=TASK_IDENTITY,
        pid=1000 + index,
        process_start_ticks=9000 + index,
        executable_device=20 + index,
        executable_inode=30 + index,
        cmdline_identity=f"sha256-{40 + index:064x}",
        socket_device=None if name == "runner" else 50 + index,
        socket_inode=None if name == "runner" else 60 + index,
    )


def _observation(policy: GatePolicy) -> GateObservation:
    _paths(policy)
    all_paths = {
        **policy.delete_paths(TASK_IDENTITY),
        **policy.retain_paths(),
    }
    resources = []
    for kind, path in all_paths.items():
        observed_path = path.parent if kind == "dockerd-socket" else path
        resources.append(
            ResourceEvidence.create(
                kind=kind,
                path=path,
                observed=observed_path.stat(),
            )
        )
    caches = tuple(
        CacheEvidence(
            kind,
            f"sha256-{100 + index:064x}",
            policy.gate_root.stat().st_dev,
            policy.gate_root.stat().st_ino + index,
        )
        for index, kind in enumerate(REQUIRED_CACHE_KINDS, start=1)
    )
    executor = gate_contract._observe_cleanup_executor(policy)
    load_limit = min(
        1_000_000,
        (os.cpu_count() or 1) * policy.load_per_cpu_milli,
    )
    return GateObservation.create(
        task_identity=TASK_IDENTITY,
        head_sha=HEAD_SHA,
        release_identity=RELEASE_IDENTITY,
        runtime_identities={
            ENCODE_RUNTIME: ENCODE_IDENTITY,
            BULK_RNASEQ_RUNTIME: BULK_IDENTITY,
        },
        environment_identity=ENVIRONMENT_IDENTITY,
        resources=resources,
        processes=tuple(_process(name) for name in PROCESS_NAMES),
        caches=caches,
        cleanup_executor=executor,
        disk_free_bytes=policy.disk_required_bytes * 2,
        disk_required_bytes=policy.disk_required_bytes,
        load_milli=0,
        load_limit_milli=load_limit,
    )


class StaticObserver:
    def __init__(self, observation: GateObservation) -> None:
        self.observation = observation
        self.calls = 0

    def observe(self, policy: GatePolicy, task_identity: str) -> GateObservation:
        assert policy.task_root(task_identity) == policy.task_root(TASK_IDENTITY)
        self.calls += 1
        return self.observation


def _request(policy: GatePolicy) -> tuple[CleanupPlan, bytes, GateRequest]:
    observation = _observation(policy)
    return prepare_gate_request(
        policy=policy,
        observer=StaticObserver(observation),
        task_identity=TASK_IDENTITY,
        head_sha=HEAD_SHA,
        release_identity=RELEASE_IDENTITY,
        runtime_identities={
            ENCODE_RUNTIME: ENCODE_IDENTITY,
            BULK_RNASEQ_RUNTIME: BULK_IDENTITY,
        },
    )


def _arguments() -> list[str]:
    return [
        "--head-sha",
        HEAD_SHA,
        "--release-identity",
        RELEASE_IDENTITY,
        "--encode-runtime-identity",
        ENCODE_IDENTITY,
        "--bulk-rnaseq-runtime-identity",
        BULK_IDENTITY,
        "--task-identity",
        TASK_IDENTITY,
    ]


def test_public_cli_rejects_caller_reported_process_cache_paths_and_thresholds(
    capsys: pytest.CaptureFixture[str],
) -> None:
    assert prepare_module.main(_arguments()) == 70
    unavailable = json.loads(capsys.readouterr().err)
    assert unavailable["issue"]["code"] == "GATE_OBSERVER_UNAVAILABLE"

    forbidden_arguments = (
        ("--runner-pid", "123"),
        ("--runner-executable-inode", "456"),
        ("--runner-cmdline-identity", ENVIRONMENT_IDENTITY),
        ("--cache", f"actions={ENVIRONMENT_IDENTITY}"),
        ("--required-disk-bytes", "1"),
        ("--max-load-milli", "999999"),
        ("--output-directory", "/tmp/caller-selected"),
        ("--task-root", "/tmp/caller-selected"),
    )
    for name, value in forbidden_arguments:
        assert prepare_module.main([*_arguments(), name, value]) == 64
        error = json.loads(capsys.readouterr().err)
        assert error["issue"]["code"] == "GATE_PREPARATION_ARGUMENT_INVALID"


def test_fixed_policy_observation_binds_processes_resources_cache_and_executor(
    tmp_path: Path,
) -> None:
    policy = _policy(tmp_path)
    plan, script, request = _request(policy)

    assert request.head_sha == HEAD_SHA
    assert request.release_identity == RELEASE_IDENTITY
    assert request.environment_identity == ENVIRONMENT_IDENTITY
    assert request.observation_identity.startswith("sha256-")
    assert {item.component: item.identity for item in request.runtime_identities} == {
        BULK_RNASEQ_RUNTIME: BULK_IDENTITY,
        ENCODE_RUNTIME: ENCODE_IDENTITY,
    }
    assert tuple(item.name for item in plan.processes) == PROCESS_NAMES
    assert all(item.process_start_ticks > 0 for item in plan.processes)
    assert all(item.executable_inode > 0 for item in plan.processes)
    assert plan.executor.file_identity in script.decode("ascii")
    assert plan.executor.closure_identity in script.decode("ascii")
    assert plan.executor.backend_identity in script.decode("ascii")
    assert str(policy.cleanup_executor) in script.decode("ascii")
    assert (
        verify_cleanup_script(plan, script, policy=policy)
        == request.cleanup_script_identity
    )
    assert b"exec /usr/bin/sudo -n " in script
    assert b"rm " not in script

    other_policy = _policy(tmp_path / "other")
    with pytest.raises(DeploymentError) as wrong_executor:
        verify_cleanup_script(plan, script, policy=other_policy)
    assert wrong_executor.value.issue.code == "GATE_CLEANUP_SCRIPT_INVALID"


def test_cleanup_plan_rejects_retained_delete_overlap_in_both_directions(
    tmp_path: Path,
) -> None:
    policy = _policy(tmp_path)
    plan, _script, _request_value = _request(policy)
    raw = plan.to_dict()
    targets = list(raw["targets"])
    historical_index = next(
        index
        for index, item in enumerate(targets)
        if item["kind"] == "historical-evidence"
    )
    original_historical = plan.targets[historical_index]
    overlapping = CleanupTarget.create(
        kind="historical-evidence",
        disposition="retain",
        path=policy.task_root(TASK_IDENTITY),
        device=original_historical.device,
        inode=original_historical.inode,
    )
    targets[historical_index] = overlapping.to_dict()
    raw["targets"] = targets
    raw["identity"] = canonical_identity(
        {key: value for key, value in raw.items() if key != "identity"},
        scheme=CLEANUP_PLAN_IDENTITY_SCHEME,
    )

    with pytest.raises(DeploymentError) as ancestor:
        CleanupPlan.from_dict(raw)
    assert ancestor.value.issue.code == "GATE_CLEANUP_PLAN_INVALID"

    raw = plan.to_dict()
    targets = list(raw["targets"])
    historical_index = next(
        index
        for index, item in enumerate(targets)
        if item["kind"] == "historical-evidence"
    )
    original_historical = plan.targets[historical_index]
    nested = CleanupTarget.create(
        kind="historical-evidence",
        disposition="retain",
        path=policy.task_root(TASK_IDENTITY) / "checkout" / "retained-child",
        device=original_historical.device,
        inode=original_historical.inode,
    )
    targets[historical_index] = nested.to_dict()
    raw["targets"] = targets
    raw["identity"] = canonical_identity(
        {key: value for key, value in raw.items() if key != "identity"},
        scheme=CLEANUP_PLAN_IDENTITY_SCHEME,
    )
    with pytest.raises(DeploymentError) as descendant:
        CleanupPlan.from_dict(raw)
    assert descendant.value.issue.code == "GATE_CLEANUP_PLAN_INVALID"


def test_observation_must_match_requested_head_release_runtime_and_fixed_policy(
    tmp_path: Path,
) -> None:
    policy = _policy(tmp_path)
    observation = _observation(policy)

    with pytest.raises(DeploymentError) as head:
        prepare_gate_request(
            policy=policy,
            observer=StaticObserver(observation),
            task_identity=TASK_IDENTITY,
            head_sha="f" * 40,
            release_identity=RELEASE_IDENTITY,
            runtime_identities={
                ENCODE_RUNTIME: ENCODE_IDENTITY,
                BULK_RNASEQ_RUNTIME: BULK_IDENTITY,
            },
        )
    assert head.value.issue.code == "GATE_OBSERVATION_MISMATCH"

    tampered = observation.to_dict()
    tampered["disk_required_bytes"] += 1
    tampered["disk_free_bytes"] += 1
    tampered["identity"] = canonical_identity(
        {key: value for key, value in tampered.items() if key != "identity"},
        scheme=GATE_OBSERVATION_IDENTITY_SCHEME,
    )
    with pytest.raises(DeploymentError) as policy_mismatch:
        prepare_gate_request(
            policy=policy,
            observer=StaticObserver(GateObservation.from_dict(tampered)),
            task_identity=TASK_IDENTITY,
            head_sha=HEAD_SHA,
            release_identity=RELEASE_IDENTITY,
            runtime_identities={
                ENCODE_RUNTIME: ENCODE_IDENTITY,
                BULK_RNASEQ_RUNTIME: BULK_IDENTITY,
            },
        )
    assert policy_mismatch.value.issue.code == "GATE_OBSERVATION_MISMATCH"


def test_observer_errors_are_redacted(tmp_path: Path) -> None:
    policy = _policy(tmp_path)

    class LeakyObserver:
        def observe(self, policy, task_identity):
            del policy, task_identity
            raise fail("PRIVATE_OBSERVER_ERROR", "/private/path token=secret")

    with pytest.raises(DeploymentError) as captured:
        prepare_gate_request(
            policy=policy,
            observer=LeakyObserver(),
            task_identity=TASK_IDENTITY,
            head_sha=HEAD_SHA,
            release_identity=RELEASE_IDENTITY,
            runtime_identities={
                ENCODE_RUNTIME: ENCODE_IDENTITY,
                BULK_RNASEQ_RUNTIME: BULK_IDENTITY,
            },
        )

    assert captured.value.issue.code == "GATE_OBSERVATION_FAILED"
    assert "private" not in str(captured.value)
    assert "secret" not in str(captured.value)


def test_filesystem_observer_fails_closed_when_fixed_material_is_missing_or_unsafe(
    tmp_path: Path,
) -> None:
    policy = _policy(tmp_path)
    _mkdir(policy.gate_root)
    task_root = _mkdir(policy.task_root(TASK_IDENTITY))

    with pytest.raises(DeploymentError) as missing:
        FilesystemGateObserver().observe(policy, TASK_IDENTITY)
    assert missing.value.issue.code == "GATE_OBSERVATION_INVALID"

    task_root.rmdir()
    outside = _mkdir(tmp_path / "outside")
    task_root.symlink_to(outside, target_is_directory=True)
    with pytest.raises(DeploymentError) as symlink:
        FilesystemGateObserver().observe(policy, TASK_IDENTITY)
    assert symlink.value.issue.code == "GATE_OBSERVATION_INVALID"


@pytest.mark.parametrize(
    ("fault_point", "request_exists"),
    (("cleanup-evidence-synced", False), ("request-synced", True)),
)
def test_prepare_fsyncs_cleanup_before_request_and_retains_interruption_evidence(
    tmp_path: Path,
    fault_point: str,
    request_exists: bool,
) -> None:
    policy = _policy(tmp_path)
    observation = _observation(policy)
    arguments = prepare_module._parser().parse_args(_arguments())

    def fail_at(point: str) -> None:
        if point == fault_point:
            raise RuntimeError("fault injection")

    with pytest.raises(RuntimeError, match="fault injection"):
        prepare_module.prepare(
            arguments,
            observer=StaticObserver(observation),
            policy=policy,
            fault=fail_at,
        )

    output = policy.output_directory(TASK_IDENTITY)
    assert (output / "cleanup-plan.json").exists()
    assert (output / "cleanup.sh").exists()
    assert (output / "gate-request.json").exists() is request_exists
    assert stat.S_IMODE((output / "cleanup-plan.json").stat().st_mode) == 0o600
    assert stat.S_IMODE((output / "cleanup.sh").stat().st_mode) == 0o500


def test_prepare_writes_only_fixed_owned_output_after_observation(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    policy = _policy(tmp_path)
    observation = _observation(policy)
    observer = StaticObserver(observation)

    assert prepare_module.main(_arguments(), observer=observer, policy=policy) == 0
    receipt = json.loads(capsys.readouterr().out)
    output = policy.output_directory(TASK_IDENTITY)
    assert observer.calls == 1
    assert receipt["status"] == "prepared"
    assert str(tmp_path) not in json.dumps(receipt)
    assert stat.S_IMODE((output / "gate-request.json").stat().st_mode) == 0o600
    persisted = GateRequest.from_dict(
        json.loads((output / "gate-request.json").read_text(encoding="utf-8"))
    )
    assert persisted.identity == receipt["request_identity"]


class StaticVerifier:
    def __init__(self, *, failing_stage: str | None = None) -> None:
        self.failing_stage = failing_stage

    def verify(self, request: GateRequest, stage: str) -> GateStageVerification:
        failed = stage == self.failing_stage
        return GateStageVerification.create(
            request_identity=request.identity,
            task_identity=request.task_identity,
            stage=stage,
            status="failed" if failed else "passed",
            reason_code="VERIFIED_FAILURE" if failed else "VERIFIED_OK",
            verifier_identity=VERIFIER_IDENTITY,
            evidence_identity=EVIDENCE_IDENTITY,
        )


def _running_request(tmp_path: Path) -> DeploymentGateRun:
    policy = _policy(tmp_path)
    _plan, _script, request = _request(policy)
    return (
        DeploymentGateRun.prepared(request)
        .dispatch(DISPATCH_IDENTITY)
        .approve(APPROVAL_IDENTITY)
        .start()
    )


def test_stage_cannot_accept_caller_passed_boolean_and_requires_verifier_evidence(
    tmp_path: Path,
) -> None:
    assert "passed" not in inspect.signature(DeploymentGateRun.verify_stage).parameters
    run = _running_request(tmp_path)
    run = run.verify_stage("identity", verifier=StaticVerifier())
    receipt = run.receipts[-1]

    assert receipt.status == "passed"
    assert receipt.verifier_identity == VERIFIER_IDENTITY
    assert receipt.evidence_identity == EVIDENCE_IDENTITY
    rendered = json.dumps(receipt.to_dict())
    assert '"pid"' not in rendered
    assert "cmdline" not in rendered
    assert str(tmp_path) not in rendered


def test_verified_failure_stops_all_gate_stages_except_verified_cleanup(
    tmp_path: Path,
) -> None:
    run = _running_request(tmp_path)
    verifier = StaticVerifier(failing_stage="checkout")
    run = run.verify_stage("identity", verifier=verifier)
    run = run.verify_stage("checkout", verifier=verifier)

    assert run.state == "failed-awaiting-cleanup"
    assert run.failure_stage == "checkout"
    with pytest.raises(DeploymentError):
        run.verify_stage("environment", verifier=verifier)
    terminal = run.verify_cleanup(verifier=StaticVerifier())
    assert terminal.state == "failed"


def test_every_ordered_stage_and_cleanup_requires_verifier(tmp_path: Path) -> None:
    run = _running_request(tmp_path)
    verifier = StaticVerifier()
    for stage in GATE_STAGES:
        run = run.verify_stage(stage, verifier=verifier)
    assert run.state == "awaiting-cleanup"
    assert run.verify_cleanup(verifier=verifier).state == "succeeded"


def test_verifier_exception_is_redacted(tmp_path: Path) -> None:
    class ExplodingVerifier:
        def verify(self, request, stage):
            del request, stage
            raise RuntimeError("/private/path token=secret")

    with pytest.raises(DeploymentError) as captured:
        _running_request(tmp_path).verify_stage(
            "identity", verifier=ExplodingVerifier()
        )
    assert captured.value.issue.code == "GATE_STAGE_VERIFICATION_FAILED"
    assert "private" not in str(captured.value)
    assert "secret" not in str(captured.value)

    class LeakyDeploymentVerifier:
        def verify(self, request, stage):
            del request, stage
            raise fail("PRIVATE_VERIFIER_ERROR", "/private/path token=secret")

    with pytest.raises(DeploymentError) as deployment_error:
        _running_request(tmp_path / "second").verify_stage(
            "identity", verifier=LeakyDeploymentVerifier()
        )
    assert deployment_error.value.issue.code == "GATE_STAGE_VERIFICATION_FAILED"
    assert "private" not in str(deployment_error.value)
    assert "secret" not in str(deployment_error.value)


def test_runtime_component_order_is_fixed() -> None:
    assert RUNTIME_COMPONENTS == (BULK_RNASEQ_RUNTIME, ENCODE_RUNTIME)
