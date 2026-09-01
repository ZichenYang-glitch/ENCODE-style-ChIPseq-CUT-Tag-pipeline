from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path
import stat
import subprocess
from types import SimpleNamespace

import pytest

import encode_pipeline.deployment.bundle_builder as builder_module
import encode_pipeline.deployment.native_contracts as native_module
import encode_pipeline.deployment.operator_action as operator_action_module
from encode_pipeline.deployment.bundle import BundleStore
from encode_pipeline.deployment.bundle_builder import build_encode_runtime_bundle
from encode_pipeline.deployment.canonical import canonical_json_bytes
from encode_pipeline.deployment.encode_runtime_materializer import (
    CONDA_COMPAT_RELATIVE_PATH,
    MICROMAMBA_RELATIVE_PATH,
    RUNTIME_INVENTORY_FILENAME,
    SNAKEMAKE_ACTIVATE_RELATIVE_PATH,
    OfflineEncodeRuntimeMaterializer,
    snakemake_environment_hash,
)
from encode_pipeline.deployment.errors import DeploymentError
from encode_pipeline.deployment.layout import DeploymentLayout
from encode_pipeline.deployment.models import ENCODE_RUNTIME
from encode_pipeline.deployment.native_contracts import (
    ENCODE_MICROMAMBA_PATH,
    ENCODE_PACKAGE_ARCHIVE_ROOT,
)
from encode_pipeline.deployment.operator_action import (
    ENCODE_RUNTIME_PREPARE_UNIT,
    EncodeRuntimeInventory,
    EncodeRuntimePrepareRequest,
    SystemdEncodeRuntimePreparer,
)


IDENTITY_A = f"sha256-{'a' * 64}"
IDENTITY_B = f"sha256-{'b' * 64}"
IDENTITY_C = f"sha256-{'c' * 64}"
TASK_IDENTITY = f"task-{'d' * 32}"
RUNNER_ARCHIVE = b"runner archive"
TOOL_ARCHIVE = b"tool archive"
RUNNER_YAML = b"name: runner\ndependencies:\n  - snakemake=8.30.0\n"
TOOL_YAML = b"name: tool\ndependencies:\n  - samtools=1.20\n"


def _static_elf() -> bytes:
    content = bytearray(120)
    content[:7] = b"\x7fELF\x02\x01\x01"
    content[18:20] = (62).to_bytes(2, "little")
    content[32:40] = (64).to_bytes(8, "little")
    content[54:56] = (56).to_bytes(2, "little")
    content[56:58] = (1).to_bytes(2, "little")
    content[64:68] = (1).to_bytes(4, "little")
    return bytes(content)


def _explicit_lock(filename: str, content: bytes) -> bytes:
    md5 = hashlib.md5(content, usedforsecurity=False).hexdigest()
    return (
        "# platform: linux-64\n"
        "@EXPLICIT\n"
        f"https://conda.anaconda.org/conda-forge/linux-64/{filename}#{md5}\n"
    ).encode()


@dataclass(frozen=True)
class _PreparedRuntime:
    layout: DeploymentLayout
    request: EncodeRuntimePrepareRequest
    source_root: Path
    destination: Path
    source_manifest: tuple[tuple[str, bytes], ...]


def _prepared_runtime(tmp_path: Path, monkeypatch) -> _PreparedRuntime:
    project = tmp_path / "project"
    repository = Path(__file__).resolve().parents[2]
    inventory = (repository / "docs/architecture/artifact-inventory.yaml").read_bytes()
    source_manifest = (
        ("docs/architecture/artifact-inventory.yaml", inventory),
        ("workflow/Snakefile", b'include: "rules/test.smk"\n'),
        (
            "workflow/rules/test.smk",
            b'rule test:\n    conda:\n        "../envs/tool.yml"\n',
        ),
        ("workflow/envs/runner.yml", RUNNER_YAML),
        ("workflow/envs/runner.lock", _explicit_lock("runner-1.conda", RUNNER_ARCHIVE)),
        ("workflow/envs/tool.yml", TOOL_YAML),
        ("workflow/envs/tool.lock", _explicit_lock("tool-1.conda", TOOL_ARCHIVE)),
    )
    for logical_path, content in source_manifest:
        path = project.joinpath(*Path(logical_path).parts)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(content)

    micromamba = tmp_path / "micromamba"
    micromamba.write_bytes(_static_elf())
    micromamba.chmod(0o755)
    cache = tmp_path / "archives"
    cache.mkdir()
    (cache / "runner-1.conda").write_bytes(RUNNER_ARCHIVE)
    (cache / "tool-1.conda").write_bytes(TOOL_ARCHIVE)

    class Provider:
        def __init__(self, _registry, *, project_root):
            assert project_root.is_absolute()

        def source_manifest(self):
            return source_manifest

        def capture(self, _workflow_id):
            return SimpleNamespace(
                is_failure=False,
                value=SimpleNamespace(digest="e" * 64, adapter_version="1.0.0"),
            )

    monkeypatch.setattr(builder_module, "WorkflowBuildIdentityProvider", Provider)
    monkeypatch.setattr(native_module, "WorkflowBuildIdentityProvider", Provider)
    bundle = tmp_path / "encode-runtime.tar"
    manifest = build_encode_runtime_bundle(project, micromamba, cache, bundle)
    layout = DeploymentLayout.isolated(tmp_path / "host")
    BundleStore(layout).stage(
        bundle,
        installed_owner_uid=os.getuid(),
        installed_owner_gid=os.getgid(),
    )
    destination = layout.encode_runtime_active_root(manifest.identity)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.mkdir(mode=0o700)
    destination.chmod(0o700)
    request = EncodeRuntimePrepareRequest.create(
        task_identity=TASK_IDENTITY,
        deployment_identity=manifest.identity,
        authority_platform_identity=IDENTITY_A,
        prior_state_identity=IDENTITY_B,
        candidate_state_identity=IDENTITY_C,
    )
    return _PreparedRuntime(
        layout=layout,
        request=request,
        source_root=layout.component_store(ENCODE_RUNTIME) / manifest.identity,
        destination=destination,
        source_manifest=source_manifest,
    )


class _SyntheticRunner:
    def __init__(
        self,
        *,
        return_code: int = 0,
        exception: BaseException | None = None,
        unsafe_symlink: bool = False,
        hardlinks: bool = False,
        safe_symlink: bool = False,
    ) -> None:
        self.return_code = return_code
        self.exception = exception
        self.unsafe_symlink = unsafe_symlink
        self.hardlinks = hardlinks
        self.safe_symlink = safe_symlink
        self.calls: list[tuple[tuple[str, ...], dict[str, object], bytes]] = []

    def __call__(self, argv, **kwargs):
        arguments = tuple(argv)
        lock = Path(arguments[arguments.index("--file") + 1])
        self.calls.append((arguments, kwargs, lock.read_bytes()))
        prefix = Path(arguments[arguments.index("--prefix") + 1])
        prefix.mkdir(mode=0o700)
        binary = prefix / "bin"
        binary.mkdir(mode=0o700)
        executable = binary / ("snakemake" if prefix.name == "runner" else "samtools")
        executable.write_bytes(b"#!/bin/sh\nexit 0\n")
        executable.chmod(0o755)
        if self.safe_symlink:
            (binary / "tool-link").symlink_to(executable.name)
        if self.unsafe_symlink:
            (binary / "escape").symlink_to("../../../../../outside")
        if self.hardlinks:
            os.link(executable, binary / "hardlink")
        evidence = prefix / "command-evidence"
        evidence.write_text("preserve on failure\n", encoding="utf-8")
        if self.exception is not None:
            raise self.exception
        return SimpleNamespace(returncode=self.return_code)


def _materializer(prepared: _PreparedRuntime, runner: _SyntheticRunner, **kwargs):
    kwargs.setdefault("monotonic_clock", lambda: 0.0)
    return OfflineEncodeRuntimeMaterializer(
        prepared.layout,
        service_uid=os.getuid(),
        service_gid=os.getgid(),
        installed_owner_uid=os.getuid(),
        installed_owner_gid=os.getgid(),
        command_runner=runner,
        **kwargs,
    )


def test_materializer_uses_only_indexed_offline_coordinates_and_fixed_argv(
    tmp_path: Path,
    monkeypatch,
) -> None:
    prepared = _prepared_runtime(tmp_path, monkeypatch)
    runner = _SyntheticRunner(safe_symlink=True)
    monkeypatch.setenv("HELIXWEAVE_PRIVATE_CALLER_VALUE", "must-not-pass")

    receipt = _materializer(prepared, runner, timeout_seconds=11.0).prepare(
        prepared.request
    )

    assert len(runner.calls) == 2
    runner_call, rule_call = runner.calls
    expected_executable = prepared.source_root.joinpath(
        *Path(ENCODE_MICROMAMBA_PATH).parts
    )
    for arguments, options, local_lock in runner.calls:
        assert arguments[0] == str(expected_executable)
        assert arguments[1:5] == (
            "--no-rc",
            "--root-prefix",
            str(prepared.destination / "mamba-root"),
            "create",
        )
        assert "--offline" in arguments
        assert "--always-copy" in arguments
        assert "https://" not in local_lock.decode()
        assert local_lock.startswith(b"# platform: linux-64\n@EXPLICIT\nfile://")
        assert options == {
            "stdin": subprocess.DEVNULL,
            "stdout": subprocess.DEVNULL,
            "stderr": subprocess.DEVNULL,
            "cwd": prepared.destination,
            "env": {"HOME": str(prepared.destination / "mamba-root")},
            "close_fds": True,
            "timeout": 11.0,
            "check": False,
            "shell": False,
            "umask": 0o077,
        }
    assert Path(runner_call[0][runner_call[0].index("--prefix") + 1]) == (
        prepared.destination / "runner"
    )
    expected_hash = snakemake_environment_hash(
        prepared.destination / "conda-envs", TOOL_YAML
    )
    rule_prefix = Path(rule_call[0][rule_call[0].index("--prefix") + 1])
    assert rule_prefix == prepared.destination / "conda-envs" / expected_hash
    assert len(rule_prefix.name) == 32
    assert (
        prepared.destination / "conda-envs" / f"{expected_hash}.env_setup_done"
    ).is_file()
    assert not any(
        path.name.startswith(".partial")
        for path in prepared.destination.parent.iterdir()
    )
    assert (
        prepared.destination.joinpath(
            *Path(MICROMAMBA_RELATIVE_PATH).parts
        ).read_bytes()
        == _static_elf()
    )
    assert (
        prepared.destination.joinpath(*Path(CONDA_COMPAT_RELATIVE_PATH).parts)
        .read_bytes()
        .startswith(b"#!/bin/sh\n")
    )
    assert (
        prepared.destination.joinpath(*Path(SNAKEMAKE_ACTIVATE_RELATIVE_PATH).parts)
        .read_bytes()
        .startswith(b"#!/bin/sh\n")
    )

    inventory_path = prepared.destination / RUNTIME_INVENTORY_FILENAME
    inventory_content = inventory_path.read_bytes()
    inventory = EncodeRuntimeInventory.from_dict(json.loads(inventory_content))
    assert inventory.to_bytes() == inventory_content
    assert receipt.request_identity == prepared.request.identity
    assert receipt.deployment_identity == prepared.request.deployment_identity
    assert receipt.tree_identity == inventory.tree_identity
    assert receipt.inventory_sha256 == hashlib.sha256(inventory_content).hexdigest()
    assert receipt.inventory_size == len(inventory_content)
    assert receipt.entry_count == len(inventory.entries)
    assert RUNTIME_INVENTORY_FILENAME not in {item.path for item in inventory.entries}
    assert any(item.kind == "symlink" for item in inventory.entries)
    execution_entries = {item.path: item for item in inventory.entries}
    assert {
        CONDA_COMPAT_RELATIVE_PATH,
        MICROMAMBA_RELATIVE_PATH,
        SNAKEMAKE_ACTIVATE_RELATIVE_PATH,
    } <= set(execution_entries)
    assert all(
        execution_entries[path].kind == "file" and execution_entries[path].mode == 0o555
        for path in (
            CONDA_COMPAT_RELATIVE_PATH,
            MICROMAMBA_RELATIVE_PATH,
            SNAKEMAKE_ACTIVATE_RELATIVE_PATH,
        )
    )
    assert (
        execution_entries[MICROMAMBA_RELATIVE_PATH].sha256
        == hashlib.sha256(_static_elf()).hexdigest()
    )
    assert stat.S_IMODE(inventory_path.stat().st_mode) == 0o600


def test_materializer_provides_micromamba_a_deterministic_writable_home(
    tmp_path: Path,
    monkeypatch,
) -> None:
    prepared = _prepared_runtime(tmp_path, monkeypatch)

    class HomeRequiringRunner(_SyntheticRunner):
        def __call__(self, argv, **kwargs):
            expected = prepared.destination / "mamba-root"
            if kwargs.get("env") != {"HOME": str(expected)} or not expected.is_dir():
                return SimpleNamespace(returncode=1)
            return super().__call__(argv, **kwargs)

    runner = HomeRequiringRunner()

    _materializer(prepared, runner).prepare(prepared.request)

    assert runner.calls
    assert all(
        options["env"] == {"HOME": str(prepared.destination / "mamba-root")}
        for _arguments, options, _lock in runner.calls
    )


def test_materialized_conda_seam_supports_only_snakemake_830_read_only_probes(
    tmp_path: Path,
    monkeypatch,
) -> None:
    prepared = _prepared_runtime(tmp_path, monkeypatch)
    runner = _SyntheticRunner()

    _materializer(prepared, runner).prepare(prepared.request)

    wrapper = prepared.destination.joinpath(*Path(CONDA_COMPAT_RELATIVE_PATH).parts)
    environment = {
        "LANG": "C.UTF-8",
        "LC_ALL": "C.UTF-8",
        "PATH": "/usr/sbin:/usr/bin:/sbin:/bin",
    }
    info = subprocess.run(
        (str(wrapper), "info", "--json"),
        check=True,
        capture_output=True,
        env=environment,
    )
    version = subprocess.run(
        (str(wrapper), "--version"),
        check=True,
        capture_output=True,
        env=environment,
    )
    config = subprocess.run(
        (str(wrapper), "config", "--get", "channel_priority", "--json"),
        check=True,
        capture_output=True,
        env=environment,
    )
    rejected = subprocess.run(
        (str(wrapper), "env", "create", "--file", "/private/caller.yml"),
        check=False,
        capture_output=True,
        env=environment,
    )

    assert json.loads(info.stdout) == {
        "conda_prefix": str(prepared.destination / "mamba-root"),
        "platform": "linux-64",
    }
    assert version.stdout == b"conda 24.7.1\n"
    assert json.loads(config.stdout) == {"get": {"channel_priority": "strict"}}
    assert rejected.returncode == 64
    assert rejected.stdout == rejected.stderr == b""


def test_materialized_activate_seam_selects_only_a_full_hash_rule_prefix(
    tmp_path: Path,
    monkeypatch,
) -> None:
    prepared = _prepared_runtime(tmp_path, monkeypatch)
    runner = _SyntheticRunner()

    _materializer(prepared, runner).prepare(prepared.request)

    full_hash = snakemake_environment_hash(
        prepared.destination / "conda-envs", TOOL_YAML
    )
    prefix = prepared.destination / "conda-envs" / full_hash
    activate = prepared.destination.joinpath(
        *Path(SNAKEMAKE_ACTIVATE_RELATIVE_PATH).parts
    )
    command = (
        'source "$1" "$2" && '
        '[ "$CONDA_PREFIX" = "$2" ] && '
        '[ "$(command -v samtools)" = "$2/bin/samtools" ]'
    )
    accepted = subprocess.run(
        ("/bin/bash", "-c", command, "helixweave-activate", str(activate), str(prefix)),
        check=False,
        capture_output=True,
        env={"PATH": "/usr/sbin:/usr/bin:/sbin:/bin"},
    )
    rejected = subprocess.run(
        (
            "/bin/bash",
            "-c",
            'source "$1" "$2"',
            "helixweave-activate",
            str(activate),
            str(tmp_path / "caller-selected"),
        ),
        check=False,
        capture_output=True,
        env={"PATH": "/usr/sbin:/usr/bin:/sbin:/bin"},
    )

    assert accepted.returncode == 0
    assert accepted.stdout == accepted.stderr == b""
    assert rejected.returncode == 64
    assert rejected.stdout == rejected.stderr == b""


def test_materializer_shares_one_total_timeout_budget(
    tmp_path: Path,
    monkeypatch,
) -> None:
    prepared = _prepared_runtime(tmp_path, monkeypatch)
    runner = _SyntheticRunner()
    observed_times = iter((0.0, 2.0, 5.0))

    def clock() -> float:
        return next(observed_times, 5.0)

    _materializer(
        prepared,
        runner,
        timeout_seconds=11.0,
        monotonic_clock=clock,
    ).prepare(prepared.request)

    assert [call[1]["timeout"] for call in runner.calls] == [9.0, 6.0]


def test_materializer_deadline_covers_post_create_inventory_work(
    tmp_path: Path,
    monkeypatch,
) -> None:
    prepared = _prepared_runtime(tmp_path, monkeypatch)
    runner = _SyntheticRunner()
    calls = 0

    def clock() -> float:
        nonlocal calls
        calls += 1
        return 0.0 if calls <= 3 else 12.0

    with pytest.raises(DeploymentError) as caught:
        _materializer(
            prepared,
            runner,
            timeout_seconds=11.0,
            monotonic_clock=clock,
        ).prepare(prepared.request)

    assert caught.value.issue.code == "ENCODE_RUNTIME_MATERIALIZATION_FAILED"
    assert len(runner.calls) == 2
    assert prepared.destination.is_dir()
    assert not (prepared.destination / RUNTIME_INVENTORY_FILENAME).exists()


def test_snakemake_environment_hash_matches_version_8_location_formula(
    tmp_path: Path,
) -> None:
    conda_prefix = tmp_path / "scientific" / "conda-envs"
    expected = hashlib.md5(usedforsecurity=False)
    expected.update(os.path.realpath(conda_prefix).encode())
    expected.update(TOOL_YAML)

    observed = snakemake_environment_hash(conda_prefix, TOOL_YAML)

    assert observed == expected.hexdigest()
    assert len(observed) == 32


def test_materializer_breaks_all_output_hardlinks(
    tmp_path: Path,
    monkeypatch,
) -> None:
    prepared = _prepared_runtime(tmp_path, monkeypatch)
    runner = _SyntheticRunner(hardlinks=True)

    _materializer(prepared, runner).prepare(prepared.request)

    regular_files = [
        path
        for path in prepared.destination.rglob("*")
        if path.is_file() and not path.is_symlink()
    ]
    assert regular_files
    assert all(path.stat(follow_symlinks=False).st_nlink == 1 for path in regular_files)
    inventory = EncodeRuntimeInventory.from_dict(
        json.loads((prepared.destination / RUNTIME_INVENTORY_FILENAME).read_bytes())
    )
    executable_entries = {
        item.path: item for item in inventory.entries if item.path.endswith("/hardlink")
    }
    assert executable_entries
    assert all(item.mode == 0o555 for item in executable_entries.values())


@pytest.mark.parametrize(
    ("runner", "expected_code"),
    (
        (_SyntheticRunner(return_code=17), "ENCODE_RUNTIME_MATERIALIZATION_FAILED"),
        (
            _SyntheticRunner(
                exception=subprocess.TimeoutExpired(("micromamba",), timeout=1)
            ),
            "ENCODE_RUNTIME_MATERIALIZATION_FAILED",
        ),
        (_SyntheticRunner(unsafe_symlink=True), "ENCODE_RUNTIME_OUTPUT_INVALID"),
    ),
)
def test_failure_is_redacted_and_preserves_exact_identity_directory(
    tmp_path: Path,
    monkeypatch,
    runner: _SyntheticRunner,
    expected_code: str,
) -> None:
    prepared = _prepared_runtime(tmp_path, monkeypatch)

    with pytest.raises(DeploymentError) as caught:
        _materializer(prepared, runner).prepare(prepared.request)

    assert caught.value.issue.code == expected_code
    assert caught.value.issue.message in {
        "ENCODE runtime materialization failed.",
        "ENCODE runtime output is invalid.",
    }
    assert prepared.destination.is_dir()
    assert list(prepared.destination.rglob("command-evidence"))
    assert not any(
        path.name.startswith(".partial")
        for path in prepared.destination.parent.iterdir()
    )
    assert not (prepared.destination / RUNTIME_INVENTORY_FILENAME).exists()


@pytest.mark.parametrize("fault", ("nonempty", "mode", "owner", "symlink"))
def test_destination_must_be_exact_empty_service_owned_0700(
    tmp_path: Path,
    monkeypatch,
    fault: str,
) -> None:
    prepared = _prepared_runtime(tmp_path, monkeypatch)
    if fault == "nonempty":
        (prepared.destination / "caller-file").write_text("do not remove")
    elif fault == "mode":
        prepared.destination.chmod(0o755)
    elif fault == "symlink":
        prepared.destination.rmdir()
        elsewhere = tmp_path / "elsewhere"
        elsewhere.mkdir(mode=0o700)
        prepared.destination.symlink_to(elsewhere, target_is_directory=True)
    runner = _SyntheticRunner()
    materializer = _materializer(prepared, runner)
    if fault == "owner":
        materializer = OfflineEncodeRuntimeMaterializer(
            prepared.layout,
            service_uid=os.getuid() + 100_000,
            service_gid=os.getgid(),
            installed_owner_uid=os.getuid(),
            installed_owner_gid=os.getgid(),
            command_runner=runner,
        )

    with pytest.raises(DeploymentError) as caught:
        materializer.prepare(prepared.request)

    assert caught.value.issue.code == "ENCODE_RUNTIME_DESTINATION_INVALID"
    assert not runner.calls
    if fault == "nonempty":
        assert (prepared.destination / "caller-file").read_text() == "do not remove"


def test_invalid_installed_bundle_is_rejected_before_any_output(
    tmp_path: Path,
    monkeypatch,
) -> None:
    prepared = _prepared_runtime(tmp_path, monkeypatch)
    source_index = (
        prepared.source_root / "payload/contracts/encode-runtime/package-index.json"
    )
    source_index.chmod(0o644)
    source_index.write_bytes(b"tampered")
    source_index.chmod(0o444)
    runner = _SyntheticRunner()

    with pytest.raises(DeploymentError) as caught:
        _materializer(prepared, runner).prepare(prepared.request)

    assert caught.value.issue.code == "ENCODE_RUNTIME_SOURCE_INVALID"
    assert list(prepared.destination.iterdir()) == []
    assert not runner.calls


def test_materializer_rejects_caller_selected_parameters(
    tmp_path: Path,
    monkeypatch,
) -> None:
    prepared = _prepared_runtime(tmp_path, monkeypatch)
    runner = _SyntheticRunner()

    with pytest.raises(DeploymentError) as caught:
        _materializer(prepared, runner).prepare(  # type: ignore[arg-type]
            {"deployment_identity": prepared.request.deployment_identity}
        )

    assert caught.value.issue.code == "ENCODE_RUNTIME_PREPARE_REQUEST_INVALID"
    assert list(prepared.destination.iterdir()) == []
    assert not runner.calls


def test_encode_unit_timeout_is_outside_the_materializer_budget(
    monkeypatch,
) -> None:
    observed: dict[str, object] = {}

    def run(argv, **kwargs):
        observed["argv"] = tuple(argv)
        observed.update(kwargs)
        return subprocess.CompletedProcess(argv, 0)

    monkeypatch.setattr(operator_action_module.subprocess, "run", run)
    argv = (
        "/usr/bin/systemctl",
        "--no-ask-password",
        "start",
        "--",
        ENCODE_RUNTIME_PREPARE_UNIT,
    )

    result = operator_action_module._run_fixed_systemctl(
        argv, ENCODE_RUNTIME_PREPARE_UNIT
    )

    assert result == 0
    assert observed["argv"] == argv
    assert observed["timeout"] == 1_860.0
    assert observed["timeout"] > operator_action_module._ACTION_TIMEOUT_SECONDS


@pytest.mark.parametrize(
    ("active_state", "expected"),
    ((b"inactive\n", True), (b"failed\n", True), (b"active\n", False)),
)
def test_encode_unit_quiesce_stops_then_requires_a_terminal_state(
    monkeypatch,
    active_state: bytes,
    expected: bool,
) -> None:
    calls: list[tuple[str, ...]] = []

    def run(argv, **kwargs):
        arguments = tuple(argv)
        calls.append(arguments)
        if "show" in arguments:
            kwargs["stdout"].write(active_state)
        return subprocess.CompletedProcess(argv, 0)

    monkeypatch.setattr(operator_action_module.subprocess, "run", run)

    observed = operator_action_module._quiesce_fixed_systemd_unit(
        ENCODE_RUNTIME_PREPARE_UNIT
    )

    assert observed is expected
    assert [call[2] for call in calls] == ["stop", "show"]
    assert all("--" in call for call in calls)


@pytest.mark.parametrize("quiescent", (True, False))
def test_prepare_never_quarantines_before_the_unit_is_terminal(
    tmp_path: Path,
    quiescent: bool,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    events: list[str] = []

    def dispatch(_argv: tuple[str, ...]) -> int:
        events.append("dispatch")
        return 1

    def quiesce() -> bool:
        events.append("quiesce")
        return quiescent

    preparer = SystemdEncodeRuntimePreparer(
        layout,
        service_uid=os.getuid(),
        service_gid=os.getgid(),
        root_uid=os.getuid(),
        root_gid=os.getgid(),
        command_runner=dispatch,
        unit_quiescer=quiesce,
    )
    preparer._require_boundary = lambda: tmp_path  # type: ignore[method-assign]
    preparer._reuse_materialization = (  # type: ignore[method-assign]
        lambda _destination, _request: None
    )
    preparer._prepare_destination = (  # type: ignore[method-assign]
        lambda _destination: events.append("prepare-destination")
    )
    preparer._replace_request = (  # type: ignore[method-assign]
        lambda _boundary, _request: events.append("write-request")
    )
    preparer._prepare_receipt = (  # type: ignore[method-assign]
        lambda _boundary: events.append("prepare-receipt")
    )
    preparer._quarantine_failed = (  # type: ignore[method-assign]
        lambda _destination, _request: events.append("quarantine")
    )
    request = EncodeRuntimePrepareRequest.create(
        task_identity=TASK_IDENTITY,
        deployment_identity=IDENTITY_A,
        authority_platform_identity=IDENTITY_B,
        prior_state_identity=IDENTITY_B,
        candidate_state_identity=IDENTITY_C,
    )

    with pytest.raises(DeploymentError) as caught:
        preparer.prepare(request)

    assert events[:5] == [
        "prepare-destination",
        "write-request",
        "prepare-receipt",
        "dispatch",
        "quiesce",
    ]
    if quiescent:
        assert events[-1] == "quarantine"
        assert caught.value.issue.code == "ENCODE_RUNTIME_PREPARE_FAILED"
    else:
        assert "quarantine" not in events
        assert caught.value.issue.code == "ENCODE_RUNTIME_RECOVERY_REQUIRED"


def test_systemd_materializer_prepares_then_rejects_tampered_reuse_record(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    prepared = _prepared_runtime(tmp_path, monkeypatch)
    prepared.destination.rmdir()
    prepared.destination.parent.chmod(0o555)
    boundary = prepared.layout.encode_runtime_prepare_root
    boundary.mkdir(parents=True)
    boundary.chmod(0o750)
    records = prepared.layout.data_root / "operator" / "runtime-materializations"
    records.mkdir(parents=True)
    records.chmod(0o700)
    runner = _SyntheticRunner(safe_symlink=True)
    commands: list[tuple[str, ...]] = []

    def dispatch(argv: tuple[str, ...]) -> int:
        commands.append(argv)
        request = EncodeRuntimePrepareRequest.from_dict(
            json.loads(prepared.layout.encode_runtime_prepare_request.read_bytes())
        )
        receipt = _materializer(prepared, runner).prepare(request)
        prepared.layout.encode_runtime_prepare_receipt.write_bytes(
            canonical_json_bytes(receipt.to_dict())
        )
        prepared.layout.encode_runtime_prepare_receipt.chmod(0o600)
        return 0

    preparer = SystemdEncodeRuntimePreparer(
        prepared.layout,
        service_uid=os.getuid(),
        service_gid=os.getgid(),
        root_uid=os.getuid(),
        root_gid=os.getgid(),
        command_runner=dispatch,
    )
    original_mkdir = operator_action_module.os.mkdir

    def privileged_destination_mkdir(path, mode=0o777, *, dir_fd=None):
        if Path(path) == prepared.destination and dir_fd is None:
            prepared.destination.parent.chmod(0o755)
            try:
                return original_mkdir(path, mode)
            finally:
                prepared.destination.parent.chmod(0o555)
        if dir_fd is None:
            return original_mkdir(path, mode)
        return original_mkdir(path, mode, dir_fd=dir_fd)

    monkeypatch.setattr(
        operator_action_module.os,
        "mkdir",
        privileged_destination_mkdir,
    )

    receipt = preparer.prepare(prepared.request)

    record = prepared.layout.encode_runtime_materialization_receipt(
        prepared.request.deployment_identity
    )
    inventory_path = prepared.destination / RUNTIME_INVENTORY_FILENAME
    assert receipt.request_identity == prepared.request.identity
    assert commands == [
        (
            "/usr/bin/systemctl",
            "--no-ask-password",
            "start",
            "--",
            ENCODE_RUNTIME_PREPARE_UNIT,
        )
    ]
    assert stat.S_IMODE(prepared.destination.stat().st_mode) == 0o555
    assert stat.S_IMODE(inventory_path.stat().st_mode) == 0o444
    assert stat.S_IMODE(record.stat().st_mode) == 0o400
    assert json.loads(record.read_bytes())["identity"] == receipt.identity

    record.chmod(0o600)
    record.write_bytes(b"{}\n")
    record.chmod(0o400)
    reuse_request = EncodeRuntimePrepareRequest.create(
        task_identity=f"task-{'e' * 32}",
        deployment_identity=prepared.request.deployment_identity,
        authority_platform_identity=prepared.request.authority_platform_identity,
        prior_state_identity=prepared.request.prior_state_identity,
        candidate_state_identity=prepared.request.candidate_state_identity,
    )

    with pytest.raises(DeploymentError) as caught:
        preparer.prepare(reuse_request)

    assert caught.value.issue.code == "ENCODE_RUNTIME_RECOVERY_REQUIRED"
    assert caught.value.issue.recoverable is True
    assert len(commands) == 1
    assert prepared.layout.encode_runtime_prepare_request.read_bytes() == (
        canonical_json_bytes(prepared.request.to_dict())
    )
    assert prepared.destination.is_dir()
    assert inventory_path.is_file()


def test_local_locks_bind_each_indexed_archive_hash(
    tmp_path: Path,
    monkeypatch,
) -> None:
    prepared = _prepared_runtime(tmp_path, monkeypatch)
    runner = _SyntheticRunner()

    _materializer(prepared, runner).prepare(prepared.request)

    package_records = {}
    for path in prepared.source_root.glob(f"{ENCODE_PACKAGE_ARCHIVE_ROOT}/*/*"):
        package_records[path.name] = (
            path.as_uri(),
            hashlib.md5(path.read_bytes(), usedforsecurity=False).hexdigest(),
        )
    observed = b"\n".join(call[2] for call in runner.calls).decode()
    assert set(package_records) == {"runner-1.conda", "tool-1.conda"}
    for uri, md5 in package_records.values():
        assert f"{uri}#{md5}" in observed
