from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, replace
import hashlib
import json
import os
from pathlib import Path
import socket
import stat
import subprocess
from types import SimpleNamespace

import pytest

import encode_pipeline.deployment.bulk_runtime_materializer as materializer_module
from encode_pipeline.adapters.bulk_rnaseq.runtime_assets import (
    RuntimeAssetBinding,
    VerifiedContainerAsset,
    VerifiedRuntimeAssets,
)
from encode_pipeline.deployment.admission import ResolvedContractFacts
from encode_pipeline.deployment.bulk_runtime_materializer import (
    IMAGE_EXACT_READY,
    IMAGE_INVALID,
    IMAGE_MISSING,
    BulkDockerBoundary,
    DOCKER_ENDPOINT,
    DOCKER_EXECUTABLE,
    MaterializedBulkRuntimeVerification,
    OfflineBulkRuntimeMaterializer,
    observe_bulk_docker_boundary,
    probe_bulk_docker_image_state,
    verify_materialized_bulk_runtime,
)
from encode_pipeline.deployment.bundle import BundleStore
from encode_pipeline.deployment.errors import DeploymentError, fail
from encode_pipeline.deployment.layout import DeploymentLayout
from encode_pipeline.deployment.models import (
    BULK_RNASEQ_RUNTIME,
    ContractIdentity,
)
from encode_pipeline.deployment.native_contracts import (
    ProductionNativeContractResolver,
)
from encode_pipeline.deployment.operator_action import (
    BulkRuntimePrepareRequest,
)
from encode_pipeline.platform.results import Issue, Result


PLATFORM_IDENTITY = f"sha256-{'a' * 64}"
BULK_IDENTITY = f"sha256-{'b' * 64}"
PRIOR_STATE_IDENTITY = f"sha256-{'c' * 64}"
CANDIDATE_STATE_IDENTITY = f"sha256-{'d' * 64}"
CONTRACT_IDENTITY = f"sha256-{'e' * 64}"
SERVICE_IDENTITY = f"sha256-{'6' * 64}"
CLIENT_IDENTITY = f"sha256-{'7' * 64}"
ENDPOINT_IDENTITY = f"sha256-{'8' * 64}"
TASK_IDENTITY = f"task-{'f' * 32}"


@dataclass(frozen=True)
class _Prepared:
    layout: DeploymentLayout
    request: BulkRuntimePrepareRequest
    assets: VerifiedRuntimeAssets
    archives: tuple[Path, ...]
    contents: tuple[bytes, ...]
    boundary: BulkDockerBoundary


class _DockerRunner:
    def __init__(self, return_codes: tuple[int, ...] = (0,)) -> None:
        self.return_codes = list(return_codes)
        self.calls: list[tuple[tuple[str, ...], dict[str, object], bytes]] = []

    def __call__(self, argv, **kwargs):
        content = bytearray()
        while chunk := os.read(kwargs["stdin"], 1024):
            content.extend(chunk)
        self.calls.append((tuple(argv), dict(kwargs), bytes(content)))
        return SimpleNamespace(returncode=self.return_codes.pop(0))


class _FullVerifier:
    def __init__(self, *results: Result[VerifiedRuntimeAssets]) -> None:
        self.results = list(results)
        self.calls = 0

    def __call__(self, _binding):
        self.calls += 1
        return self.results.pop(0)


class _StateProbe:
    def __init__(self, states: dict[str, list[str]]) -> None:
        self.states = states
        self.calls: list[str] = []

    def __call__(self, _binding, image, _rootfs):
        self.calls.append(image)
        return self.states[image].pop(0)


class _BoundaryObserver:
    def __init__(self, *boundaries: BulkDockerBoundary) -> None:
        self.boundaries = list(boundaries)
        self.calls = 0

    def __call__(self, uid: int, gid: int) -> BulkDockerBoundary:
        self.calls += 1
        observed = (
            self.boundaries.pop(0) if len(self.boundaries) > 1 else self.boundaries[0]
        )
        assert (uid, gid) == (observed.daemon_uid, observed.daemon_gid)
        return observed


def _install_docker_response(
    monkeypatch: pytest.MonkeyPatch,
    response: bytes,
    *,
    socket_error: bool = False,
) -> list[tuple[object, ...]]:
    remaining = bytearray(response)
    calls: list[tuple[object, ...]] = []

    class FakeSocket:
        def __enter__(self):
            return self

        def __exit__(self, *_args):
            return None

        def settimeout(self, value):
            calls.append(("timeout", value))

        def connect(self, path):
            calls.append(("connect", path))
            if socket_error:
                raise OSError("private socket failure")

        def sendall(self, content):
            calls.append(("request", content))

        def recv(self, _maximum):
            if not remaining:
                return b""
            content = bytes(remaining)
            remaining.clear()
            return content

    def socket_factory(family, kind):
        calls.append(("socket", family, kind))
        return FakeSocket()

    monkeypatch.setattr(materializer_module.socket, "socket", socket_factory)
    return calls


def _not_ready() -> Result[VerifiedRuntimeAssets]:
    return Result.failure((Issue("DOCKER_UNAVAILABLE", "Docker is unavailable."),))


def _prepared(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    *,
    contents: tuple[bytes, ...] = (b"first docker archive",),
) -> _Prepared:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    source = layout.component_store(BULK_RNASEQ_RUNTIME) / BULK_IDENTITY
    runtime = source / "payload" / "runtime"
    archive_root = runtime / "containers" / "assets"
    archive_root.mkdir(parents=True)
    archives: list[Path] = []
    containers: list[VerifiedContainerAsset] = []
    for index, content in enumerate(contents):
        archive = archive_root / f"image-{index}.tar"
        archive.write_bytes(content)
        archive.chmod(0o444)
        archives.append(archive)
        digest = hashlib.sha256(content).hexdigest()
        image_digest = f"sha256:{index + 1:064x}"
        containers.append(
            VerifiedContainerAsset(
                process=f"PROCESS_{index}",
                image=f"registry.invalid/image-{index}@{image_digest}",
                oci_digest=image_digest,
                local_asset=archive,
                local_sha256=digest,
                size_bytes=len(content),
                distribution_manifest=archive_root / f"image-{index}.json",
                distribution_manifest_sha256=f"{index + 2:064x}",
                config_digest=f"sha256:{index + 3:064x}",
                runtime_image=f"sha256:{index + 3:064x}",
                rootfs_diff_ids=(f"sha256:{index + 4:064x}",),
            )
        )
    assets = VerifiedRuntimeAssets(
        root=runtime,
        source_tree=runtime / "source/rnaseq",
        nextflow_executable=runtime / "nextflow/nextflow",
        jdk_archive=runtime / "jdk/corretto.tar.gz",
        jdk_tree=runtime / "jdk/corretto",
        java_executable=runtime / "jdk/corretto/bin/java",
        plugin_root=runtime / "plugins",
        plugin_archive=runtime / "plugins/nf-schema.zip",
        plugin_meta=runtime / "plugins/nf-schema.json",
        plugin_tree=runtime / "plugins/nf-schema",
        container_lock=runtime / "containers/availability-lock.json",
        containers=tuple(containers),
        source_tree_sha256="1" * 64,
        runtime_identity_sha256="2" * 64,
        nextflow_sha256="3" * 64,
        jdk_archive_sha256="4" * 64,
        jdk_tree_sha256="5" * 64,
        java_executable_sha256="6" * 64,
        plugin_archive_sha256="7" * 64,
        plugin_tree_sha256="8" * 64,
        container_inventory_sha256="9" * 64,
        container_lock_sha256="a" * 64,
    )
    contract = ContractIdentity("helixweave.runtime.bulk-rnaseq", CONTRACT_IDENTITY)
    manifest = SimpleNamespace(
        component=BULK_RNASEQ_RUNTIME,
        identity=BULK_IDENTITY,
        provides=(contract,),
    )
    facts = ResolvedContractFacts(
        component=BULK_RNASEQ_RUNTIME,
        deployment_identity=BULK_IDENTITY,
        version="3.26.0",
        contracts=(contract,),
    )

    def verify_installed(_self, component, identity, **kwargs):
        assert component == BULK_RNASEQ_RUNTIME
        assert identity == BULK_IDENTITY
        assert kwargs == {
            "expected_owner_uid": os.getuid(),
            "expected_owner_gid": os.getgid(),
        }
        return manifest

    monkeypatch.setattr(BundleStore, "verify_installed", verify_installed)
    monkeypatch.setattr(
        ProductionNativeContractResolver,
        "resolve",
        lambda _self, observed_root, observed_manifest: (
            facts if observed_root == source and observed_manifest is manifest else None
        ),
    )
    request = BulkRuntimePrepareRequest.create(
        operation="activate",
        task_identity=TASK_IDENTITY,
        candidate_bulk_identity=BULK_IDENTITY,
        authority_platform_identity=PLATFORM_IDENTITY,
        prior_state_identity=PRIOR_STATE_IDENTITY,
        candidate_state_identity=CANDIDATE_STATE_IDENTITY,
        docker_service_identity=SERVICE_IDENTITY,
        docker_client_identity=CLIENT_IDENTITY,
        docker_endpoint_identity=ENDPOINT_IDENTITY,
        docker_daemon_uid=os.getuid(),
        docker_daemon_gid=os.getgid(),
    )
    boundary = BulkDockerBoundary(
        CLIENT_IDENTITY,
        ENDPOINT_IDENTITY,
        os.getuid(),
        os.getgid(),
    )
    return _Prepared(layout, request, assets, tuple(archives), contents, boundary)


def _materializer(
    prepared: _Prepared,
    runner: _DockerRunner,
    full_verifier: _FullVerifier,
    *,
    state_probe: _StateProbe | None = None,
    boundary_observer: _BoundaryObserver | None = None,
    canary_verifier: Callable[[RuntimeAssetBinding], object] | None = None,
) -> OfflineBulkRuntimeMaterializer:
    if state_probe is None:
        state_probe = _StateProbe(
            {
                container.runtime_image: [IMAGE_MISSING, IMAGE_EXACT_READY]
                for container in prepared.assets.containers
            }
        )
    if boundary_observer is None:
        boundary_observer = _BoundaryObserver(prepared.boundary)
    if canary_verifier is None:

        def canary_verifier(_binding: RuntimeAssetBinding) -> object:
            return Result.success(prepared.assets)

    return OfflineBulkRuntimeMaterializer(
        prepared.layout,
        installed_owner_uid=os.getuid(),
        installed_owner_gid=os.getgid(),
        timeout_seconds=11.0,
        command_runner=runner,
        monotonic_clock=lambda: 0.0,
        static_verifier=lambda _binding: Result.success(prepared.assets),
        canary_verifier=canary_verifier,
        full_verifier=full_verifier,
        boundary_observer=boundary_observer,
        image_state_probe=state_probe,
    )


def _request_for_operation(
    prepared: _Prepared,
    operation: str,
) -> BulkRuntimePrepareRequest:
    request = prepared.request
    return BulkRuntimePrepareRequest.create(
        operation=operation,
        task_identity=request.task_identity,
        candidate_bulk_identity=request.candidate_bulk_identity,
        authority_platform_identity=request.authority_platform_identity,
        prior_state_identity=request.prior_state_identity,
        candidate_state_identity=request.candidate_state_identity,
        docker_service_identity=request.docker_service_identity,
        docker_client_identity=request.docker_client_identity,
        docker_endpoint_identity=request.docker_endpoint_identity,
        docker_daemon_uid=request.docker_daemon_uid,
        docker_daemon_gid=request.docker_daemon_gid,
    )


def test_boundary_observer_binds_client_bytes_and_socket_kernel_inode(
    tmp_path: Path,
) -> None:
    client = tmp_path / "docker"
    client.write_bytes(b"fixed docker client")
    client.chmod(0o755)
    socket_path = tmp_path / "docker.sock"
    header = "Num RefCount Protocol Flags Type St Inode Path\n"
    socket_info = SimpleNamespace(
        st_dev=11,
        st_ino=12,
        st_mode=stat.S_IFSOCK | 0o600,
        st_nlink=1,
        st_uid=os.getuid(),
        st_gid=os.getgid(),
        st_size=0,
        st_mtime_ns=1,
        st_ctime_ns=1,
    )

    def proc(kernel_inode: int) -> bytes:
        return (
            header
            + "00000000: 00000002 00000000 00010000 0001 01 "
            + f"{kernel_inode} {socket_path}\n"
        ).encode()

    common = {
        "_client_path": client,
        "_socket_path": socket_path,
        "_socket_lstat": lambda _path: socket_info,
        "_client_owner_uid": os.getuid(),
        "_client_owner_gid": os.getgid(),
    }
    first = observe_bulk_docker_boundary(
        os.getuid(),
        os.getgid(),
        _proc_unix_reader=lambda: proc(91),
        **common,
    )
    changed_kernel_inode = observe_bulk_docker_boundary(
        os.getuid(),
        os.getgid(),
        _proc_unix_reader=lambda: proc(92),
        **common,
    )

    assert first.client_identity.startswith("sha256-")
    assert first.endpoint_identity.startswith("sha256-")
    assert changed_kernel_inode.client_identity == first.client_identity
    assert changed_kernel_inode.endpoint_identity != first.endpoint_identity


@pytest.mark.parametrize(
    ("status", "document", "expected"),
    (
        (200, "exact", IMAGE_EXACT_READY),
        (404, "missing", IMAGE_MISSING),
        (200, "wrong-rootfs", IMAGE_INVALID),
        (500, "failure", IMAGE_INVALID),
    ),
)
def test_production_image_state_probe_is_exact_and_af_unix_only(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    status: int,
    document: str,
    expected: str,
) -> None:
    image = f"sha256:{'1' * 64}"
    layers = (f"sha256:{'2' * 64}",)
    bodies = {
        "exact": {"Id": image, "RootFS": {"Type": "layers", "Layers": list(layers)}},
        "missing": {"message": f"No such image: {image}"},
        "wrong-rootfs": {
            "Id": image,
            "RootFS": {"Type": "layers", "Layers": [f"sha256:{'3' * 64}"]},
        },
        "failure": {"message": "daemon unavailable"},
    }
    body = json.dumps(bodies[document], separators=(",", ":")).encode()
    response = (
        f"HTTP/1.0 {status} status\r\nContent-Length: {len(body)}\r\n\r\n".encode()
        + body
    )
    calls: list[tuple[object, ...]] = []

    class FakeSocket:
        def __enter__(self):
            return self

        def __exit__(self, *_args):
            return None

        def settimeout(self, value):
            calls.append(("timeout", value))

        def connect(self, path):
            calls.append(("connect", path))

        def sendall(self, content):
            calls.append(("request", content))

        def recv(self, _maximum):
            if response:
                content = bytes(response)
                response.clear()
                return content
            return b""

    response = bytearray(response)

    def socket_factory(family, kind):
        calls.append(("socket", family, kind))
        return FakeSocket()

    monkeypatch.setattr(materializer_module.socket, "socket", socket_factory)
    binding = RuntimeAssetBinding(
        root=tmp_path.resolve(),
        docker_executable=materializer_module.DOCKER_EXECUTABLE,
        docker_socket=materializer_module.DOCKER_SOCKET,
    )

    observed = probe_bulk_docker_image_state(binding, image, layers)

    assert observed == expected
    assert calls[0] == ("socket", socket.AF_UNIX, socket.SOCK_STREAM)
    assert ("connect", str(materializer_module.DOCKER_SOCKET)) in calls
    assert not any(b"tcp://" in call[1] for call in calls if call[0] == "request")


@pytest.mark.parametrize(
    ("fault", "response"),
    (
        ("invalid-binding", b""),
        ("socket-error", b""),
        ("empty", b""),
        ("missing-separator", b"HTTP/1.0 200 OK\r\n"),
        ("bad-status", b"FTP/1.0 200 OK\r\n\r\n{}"),
        (
            "duplicate-header",
            b"HTTP/1.0 200 OK\r\nX-Test: one\r\nx-test: two\r\n\r\n{}",
        ),
        (
            "transfer-encoding",
            b"HTTP/1.0 200 OK\r\nTransfer-Encoding: chunked\r\n\r\n{}",
        ),
        (
            "length-mismatch",
            b"HTTP/1.0 200 OK\r\nContent-Length: 9\r\n\r\n{}",
        ),
        (
            "duplicate-json-key",
            b'HTTP/1.0 200 OK\r\n\r\n{"Id":1,"Id":2}',
        ),
        (
            "wrong-missing-message",
            b'HTTP/1.0 404 Missing\r\n\r\n{"message":"different image"}',
        ),
    ),
)
def test_production_image_probe_rejects_malformed_or_ambiguous_daemon_responses(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    fault: str,
    response: bytes,
) -> None:
    image = f"sha256:{'1' * 64}"
    layers = (f"sha256:{'2' * 64}",)
    calls = _install_docker_response(
        monkeypatch,
        response,
        socket_error=fault == "socket-error",
    )
    docker_executable = (
        tmp_path / "docker" if fault == "invalid-binding" else DOCKER_EXECUTABLE
    )
    binding = RuntimeAssetBinding(
        root=tmp_path.resolve(),
        docker_executable=docker_executable,
        docker_socket=materializer_module.DOCKER_SOCKET,
    )

    observed = probe_bulk_docker_image_state(binding, image, layers)

    assert observed == IMAGE_INVALID
    if fault == "invalid-binding":
        assert calls == []
    else:
        assert calls[0] == ("socket", socket.AF_UNIX, socket.SOCK_STREAM)
        assert not any(b"tcp://" in call[1] for call in calls if call[0] == "request")


def test_materializer_uses_archive_fd_and_exact_offline_docker_command(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    prepared = _prepared(tmp_path, monkeypatch)
    runner = _DockerRunner()
    full = _FullVerifier(Result.success(prepared.assets))
    boundary_observer = _BoundaryObserver(prepared.boundary)
    monkeypatch.setenv("DOCKER_HOST", "tcp://attacker.invalid:2375")
    monkeypatch.setenv("DOCKER_CONFIG", str(tmp_path / "private"))

    receipt = _materializer(
        prepared,
        runner,
        full,
        boundary_observer=boundary_observer,
    ).prepare(prepared.request)

    assert full.calls == 1
    assert boundary_observer.calls == 12
    assert len(runner.calls) == 1
    argv, options, consumed = runner.calls[0]
    assert argv == (
        str(DOCKER_EXECUTABLE),
        "--host",
        DOCKER_ENDPOINT,
        "image",
        "load",
    )
    assert str(prepared.archives[0]) not in argv
    assert consumed == prepared.contents[0]
    assert options == {
        "stdin": options["stdin"],
        "stdout": subprocess.DEVNULL,
        "stderr": subprocess.DEVNULL,
        "cwd": prepared.assets.root,
        "env": {
            "HOME": "/nonexistent",
            "LANG": "C.UTF-8",
            "LC_ALL": "C.UTF-8",
            "PATH": "/usr/bin:/bin",
        },
        "close_fds": True,
        "timeout": 11.0,
        "check": False,
        "shell": False,
        "umask": 0o077,
    }
    serialized = json.dumps(receipt.to_dict(), sort_keys=True)
    assert receipt.request_identity == prepared.request.identity
    assert receipt.candidate_bulk_identity == BULK_IDENTITY
    assert receipt.image_count == 1
    assert "loaded_count" not in serialized
    assert str(tmp_path) not in serialized


@pytest.mark.parametrize("mutation", ["missing", "tampered"])
def test_archive_failure_happens_before_any_docker_side_effect(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    mutation: str,
) -> None:
    prepared = _prepared(tmp_path, monkeypatch)
    if mutation == "missing":
        prepared.archives[0].unlink()
    else:
        prepared.archives[0].chmod(0o644)
        prepared.archives[0].write_bytes(b"x" * len(prepared.contents[0]))
        prepared.archives[0].chmod(0o444)
    runner = _DockerRunner()

    with pytest.raises(DeploymentError) as caught:
        _materializer(
            prepared,
            runner,
            _FullVerifier(_not_ready()),
        ).prepare(prepared.request)

    assert caught.value.issue.code == "BULK_RUNTIME_ARCHIVE_INVALID"
    assert runner.calls == []


def test_mixed_retry_loads_only_the_still_missing_image(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    prepared = _prepared(
        tmp_path,
        monkeypatch,
        contents=(b"already ready", b"still missing"),
    )
    first, second = prepared.assets.containers
    state_probe = _StateProbe(
        {
            first.runtime_image: [IMAGE_EXACT_READY],
            second.runtime_image: [IMAGE_MISSING, IMAGE_EXACT_READY],
        }
    )
    runner = _DockerRunner()

    receipt = _materializer(
        prepared,
        runner,
        _FullVerifier(Result.success(prepared.assets)),
        state_probe=state_probe,
    ).prepare(prepared.request)

    assert receipt.image_count == 2
    assert [call[2] for call in runner.calls] == [prepared.contents[1]]


def test_any_invalid_or_unavailable_image_state_causes_zero_load(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    prepared = _prepared(
        tmp_path,
        monkeypatch,
        contents=(b"missing", b"invalid"),
    )
    first, second = prepared.assets.containers
    state_probe = _StateProbe(
        {
            first.runtime_image: [IMAGE_MISSING],
            second.runtime_image: [IMAGE_INVALID],
        }
    )
    runner = _DockerRunner()
    full = _FullVerifier(Result.success(prepared.assets))

    with pytest.raises(DeploymentError) as caught:
        _materializer(
            prepared,
            runner,
            full,
            state_probe=state_probe,
        ).prepare(prepared.request)

    assert caught.value.issue.code == "BULK_RUNTIME_IMAGE_STATE_INVALID"
    assert runner.calls == []
    assert full.calls == 0


def test_canary_failure_with_missing_image_causes_zero_load(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    prepared = _prepared(tmp_path, monkeypatch)
    runner = _DockerRunner()
    full = _FullVerifier(Result.success(prepared.assets))

    with pytest.raises(DeploymentError) as caught:
        _materializer(
            prepared,
            runner,
            full,
            canary_verifier=lambda _binding: _not_ready(),
        ).prepare(prepared.request)

    assert caught.value.issue.code == "BULK_RUNTIME_CANARY_INVALID"
    assert runner.calls == []
    assert full.calls == 0


@pytest.mark.parametrize(
    ("fault", "expected_code"),
    (
        ("boundary-exception", "BULK_RUNTIME_DOCKER_BOUNDARY_INVALID"),
        ("canary-exception", "BULK_RUNTIME_CANARY_INVALID"),
        ("canary-mismatch", "BULK_RUNTIME_CANARY_INVALID"),
        ("state-exception", "BULK_RUNTIME_IMAGE_STATE_INVALID"),
        ("state-unknown", "BULK_RUNTIME_IMAGE_STATE_INVALID"),
        ("full-exception", "BULK_RUNTIME_POST_LOAD_INVALID"),
        ("runner-exception", "BULK_RUNTIME_IMAGE_LOAD_FAILED"),
    ),
)
def test_materializer_maps_boundary_verifier_and_runner_faults(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    fault: str,
    expected_code: str,
) -> None:
    prepared = _prepared(tmp_path, monkeypatch)
    image = prepared.assets.containers[0].runtime_image

    def raise_private(*_args, **_kwargs):
        raise RuntimeError(str(tmp_path))

    def raise_process(*_args, **_kwargs):
        raise OSError(str(tmp_path))

    runner = raise_process if fault == "runner-exception" else _DockerRunner()
    full_verifier = (
        raise_private
        if fault == "full-exception"
        else _FullVerifier(Result.success(prepared.assets))
    )
    boundary_observer = (
        raise_private
        if fault == "boundary-exception"
        else _BoundaryObserver(prepared.boundary)
    )
    if fault == "canary-exception":
        canary_verifier = raise_private
    elif fault == "canary-mismatch":
        changed = replace(prepared.assets, nextflow_sha256="f" * 64)

        def changed_canary(_binding):
            return Result.success(changed)

        canary_verifier = changed_canary
    else:
        canary_verifier = None
    if fault == "state-exception":
        state_probe = raise_private
    elif fault == "state-unknown":

        def unknown_state(*_args):
            return "unexpected"

        state_probe = unknown_state
    elif fault == "full-exception":
        state_probe = _StateProbe({image: [IMAGE_EXACT_READY]})
    else:
        state_probe = None

    with pytest.raises(DeploymentError) as caught:
        _materializer(
            prepared,
            runner,
            full_verifier,
            state_probe=state_probe,
            boundary_observer=boundary_observer,
            canary_verifier=canary_verifier,
        ).prepare(prepared.request)

    assert caught.value.issue.code == expected_code
    assert caught.value.issue.recoverable is True
    assert str(tmp_path) not in str(caught.value)


def test_endpoint_change_after_missing_probe_causes_zero_load(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    prepared = _prepared(tmp_path, monkeypatch)
    changed = replace(prepared.boundary, endpoint_identity=f"sha256-{'9' * 64}")
    observer = _BoundaryObserver(
        prepared.boundary,
        prepared.boundary,
        prepared.boundary,
        prepared.boundary,
        changed,
    )
    runner = _DockerRunner()

    with pytest.raises(DeploymentError) as caught:
        _materializer(
            prepared,
            runner,
            _FullVerifier(Result.success(prepared.assets)),
            boundary_observer=observer,
        ).prepare(prepared.request)

    assert caught.value.issue.code == "BULK_RUNTIME_DOCKER_BOUNDARY_CHANGED"
    assert observer.calls == 5
    assert runner.calls == []


def test_archive_fd_is_rehashed_after_docker_returns(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    prepared = _prepared(tmp_path, monkeypatch)

    class MutatingRunner(_DockerRunner):
        def __call__(self, argv, **kwargs):
            result = super().__call__(argv, **kwargs)
            prepared.archives[0].chmod(0o400)
            return result

    runner = MutatingRunner()

    with pytest.raises(DeploymentError) as caught:
        _materializer(
            prepared,
            runner,
            _FullVerifier(_not_ready()),
        ).prepare(prepared.request)

    assert caught.value.issue.code == "BULK_RUNTIME_ARCHIVE_CHANGED"
    assert prepared.archives[0].is_file()


def test_partial_load_failure_preserves_archives_and_never_deletes_images(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    prepared = _prepared(
        tmp_path,
        monkeypatch,
        contents=(b"first archive", b"second archive"),
    )
    runner = _DockerRunner((0, 7))

    with pytest.raises(DeploymentError) as caught:
        _materializer(
            prepared,
            runner,
            _FullVerifier(_not_ready()),
        ).prepare(prepared.request)

    assert caught.value.issue.code == "BULK_RUNTIME_IMAGE_LOAD_FAILED"
    assert [call[2] for call in runner.calls] == list(prepared.contents)
    assert all(
        path.read_bytes() == content
        for path, content in zip(prepared.archives, prepared.contents, strict=True)
    )
    assert all(
        "rm" not in argv and "prune" not in argv for argv, _options, _ in runner.calls
    )


def test_already_loaded_image_set_is_idempotent_and_has_stable_receipt(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    prepared = _prepared(tmp_path, monkeypatch)
    runner = _DockerRunner()

    ready_probe = _StateProbe(
        {
            container.runtime_image: [IMAGE_EXACT_READY, IMAGE_EXACT_READY]
            for container in prepared.assets.containers
        }
    )
    first = _materializer(
        prepared,
        runner,
        _FullVerifier(Result.success(prepared.assets)),
        state_probe=ready_probe,
    ).prepare(prepared.request)
    second = _materializer(
        prepared,
        runner,
        _FullVerifier(Result.success(prepared.assets)),
        state_probe=ready_probe,
    ).prepare(prepared.request)

    assert first == second
    assert runner.calls == []


def test_post_load_inspection_mismatch_fails_closed(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    prepared = _prepared(tmp_path, monkeypatch)
    changed = replace(prepared.assets, nextflow_sha256="f" * 64)
    runner = _DockerRunner()

    with pytest.raises(DeploymentError) as caught:
        _materializer(
            prepared,
            runner,
            _FullVerifier(Result.success(changed)),
        ).prepare(prepared.request)

    assert caught.value.issue.code == "BULK_RUNTIME_POST_LOAD_INVALID"
    assert len(runner.calls) == 1


def test_request_contract_rejects_paths_environment_images_and_observe() -> None:
    request = BulkRuntimePrepareRequest.create(
        operation="rollback",
        task_identity=TASK_IDENTITY,
        candidate_bulk_identity=BULK_IDENTITY,
        authority_platform_identity=PLATFORM_IDENTITY,
        prior_state_identity=PRIOR_STATE_IDENTITY,
        candidate_state_identity=CANDIDATE_STATE_IDENTITY,
        docker_service_identity=SERVICE_IDENTITY,
        docker_client_identity=CLIENT_IDENTITY,
        docker_endpoint_identity=ENDPOINT_IDENTITY,
        docker_daemon_uid=os.getuid(),
        docker_daemon_gid=os.getgid(),
    )
    assert not any(
        key in request.to_dict()
        for key in ("path", "archive", "environment", "images", "docker_host")
    )
    for injected in (
        {"archive_path": "/private/archive.tar"},
        {"environment": {"DOCKER_HOST": "tcp://attacker.invalid"}},
        {"images": ["attacker/image:latest"]},
    ):
        raw = request.to_dict() | injected
        with pytest.raises(DeploymentError):
            BulkRuntimePrepareRequest.from_dict(raw)
    raw = request.to_dict()
    raw["operation"] = "observe"
    with pytest.raises(DeploymentError):
        BulkRuntimePrepareRequest.from_dict(raw)

    root_daemon = BulkRuntimePrepareRequest.create(
        operation="activate",
        task_identity=TASK_IDENTITY,
        candidate_bulk_identity=BULK_IDENTITY,
        authority_platform_identity=PLATFORM_IDENTITY,
        prior_state_identity=PRIOR_STATE_IDENTITY,
        candidate_state_identity=CANDIDATE_STATE_IDENTITY,
        docker_service_identity=SERVICE_IDENTITY,
        docker_client_identity=CLIENT_IDENTITY,
        docker_endpoint_identity=ENDPOINT_IDENTITY,
        docker_daemon_uid=0,
        docker_daemon_gid=0,
    )
    assert BulkRuntimePrepareRequest.from_dict(root_daemon.to_dict()) == root_daemon


def test_read_only_materialized_verifier_returns_path_free_exact_evidence(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    prepared = _prepared(tmp_path, monkeypatch)
    image = prepared.assets.containers[0].runtime_image
    observer = _BoundaryObserver(prepared.boundary)

    verified = verify_materialized_bulk_runtime(
        prepared.layout,
        BULK_IDENTITY,
        expected_daemon_uid=prepared.request.docker_daemon_uid,
        expected_daemon_gid=prepared.request.docker_daemon_gid,
        expected_client_identity=CLIENT_IDENTITY,
        expected_endpoint_identity=ENDPOINT_IDENTITY,
        installed_owner_uid=os.getuid(),
        installed_owner_gid=os.getgid(),
        _static_verifier=lambda _binding: Result.success(prepared.assets),
        _canary_verifier=lambda _binding: Result.success(prepared.assets),
        _boundary_observer=observer,
        _image_state_probe=_StateProbe({image: [IMAGE_EXACT_READY]}),
    )

    assert isinstance(verified, MaterializedBulkRuntimeVerification)
    assert verified.image_count == 1
    assert verified.docker_client_identity == CLIENT_IDENTITY
    assert verified.docker_endpoint_identity == ENDPOINT_IDENTITY
    assert observer.calls == 5
    serialized = json.dumps(verified.__dict__, sort_keys=True)
    assert str(tmp_path) not in serialized
    assert "path" not in serialized


@pytest.mark.parametrize(
    ("fault", "expected_code"),
    (
        ("static-root", "BULK_RUNTIME_SOURCE_INVALID"),
        ("boundary-exception", "BULK_RUNTIME_DOCKER_BOUNDARY_INVALID"),
        ("canary-exception", "BULK_RUNTIME_CANARY_INVALID"),
        ("canary-mismatch", "BULK_RUNTIME_CANARY_INVALID"),
        ("probe-exception", "BULK_RUNTIME_IMAGE_STATE_INVALID"),
    ),
)
def test_read_only_verifier_maps_dependency_faults_to_stable_contract_errors(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    fault: str,
    expected_code: str,
) -> None:
    prepared = _prepared(tmp_path, monkeypatch)
    image = prepared.assets.containers[0].runtime_image

    def raise_private(*_args):
        raise RuntimeError(str(tmp_path))

    static = (
        replace(prepared.assets, root=tmp_path / "wrong-root")
        if fault == "static-root"
        else prepared.assets
    )
    canary = (
        replace(prepared.assets, nextflow_sha256="f" * 64)
        if fault == "canary-mismatch"
        else prepared.assets
    )
    observer = (
        raise_private
        if fault == "boundary-exception"
        else _BoundaryObserver(prepared.boundary)
    )
    canary_verifier = (
        raise_private
        if fault == "canary-exception"
        else lambda _binding: Result.success(canary)
    )
    image_probe = (
        raise_private
        if fault == "probe-exception"
        else _StateProbe({image: [IMAGE_EXACT_READY]})
    )

    with pytest.raises(DeploymentError) as caught:
        verify_materialized_bulk_runtime(
            prepared.layout,
            BULK_IDENTITY,
            expected_daemon_uid=prepared.request.docker_daemon_uid,
            expected_daemon_gid=prepared.request.docker_daemon_gid,
            expected_client_identity=CLIENT_IDENTITY,
            expected_endpoint_identity=ENDPOINT_IDENTITY,
            installed_owner_uid=os.getuid(),
            installed_owner_gid=os.getgid(),
            _static_verifier=lambda _binding: Result.success(static),
            _canary_verifier=canary_verifier,
            _boundary_observer=observer,
            _image_state_probe=image_probe,
        )

    assert caught.value.issue.code == expected_code
    assert str(tmp_path) not in str(caught.value)


@pytest.mark.parametrize(
    ("state", "code"),
    (
        (IMAGE_MISSING, "BULK_RUNTIME_IMAGE_MISSING"),
        (IMAGE_INVALID, "BULK_RUNTIME_IMAGE_STATE_INVALID"),
    ),
)
def test_read_only_materialized_verifier_rejects_non_exact_images_without_loading(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    state: str,
    code: str,
) -> None:
    prepared = _prepared(tmp_path, monkeypatch)
    image = prepared.assets.containers[0].runtime_image

    def forbidden_prepare(*_args, **_kwargs):
        raise AssertionError("read-only verifier must never invoke the materializer")

    monkeypatch.setattr(OfflineBulkRuntimeMaterializer, "prepare", forbidden_prepare)

    with pytest.raises(DeploymentError) as caught:
        verify_materialized_bulk_runtime(
            prepared.layout,
            BULK_IDENTITY,
            expected_daemon_uid=prepared.request.docker_daemon_uid,
            expected_daemon_gid=prepared.request.docker_daemon_gid,
            expected_client_identity=CLIENT_IDENTITY,
            expected_endpoint_identity=ENDPOINT_IDENTITY,
            installed_owner_uid=os.getuid(),
            installed_owner_gid=os.getgid(),
            _static_verifier=lambda _binding: Result.success(prepared.assets),
            _canary_verifier=lambda _binding: Result.success(prepared.assets),
            _boundary_observer=_BoundaryObserver(prepared.boundary),
            _image_state_probe=_StateProbe({image: [state]}),
        )

    assert caught.value.issue.code == code


def test_read_only_materialized_verifier_rejects_boundary_drift(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    prepared = _prepared(tmp_path, monkeypatch)
    changed = replace(prepared.boundary, endpoint_identity=f"sha256-{'9' * 64}")
    observer = _BoundaryObserver(prepared.boundary, changed)
    image_probe = _StateProbe(
        {prepared.assets.containers[0].runtime_image: [IMAGE_EXACT_READY]}
    )

    with pytest.raises(DeploymentError) as caught:
        verify_materialized_bulk_runtime(
            prepared.layout,
            BULK_IDENTITY,
            expected_daemon_uid=prepared.request.docker_daemon_uid,
            expected_daemon_gid=prepared.request.docker_daemon_gid,
            expected_client_identity=CLIENT_IDENTITY,
            expected_endpoint_identity=ENDPOINT_IDENTITY,
            installed_owner_uid=os.getuid(),
            installed_owner_gid=os.getgid(),
            _static_verifier=lambda _binding: Result.success(prepared.assets),
            _canary_verifier=lambda _binding: Result.success(prepared.assets),
            _boundary_observer=observer,
            _image_state_probe=image_probe,
        )

    assert caught.value.issue.code == "BULK_RUNTIME_DOCKER_BOUNDARY_CHANGED"
    assert image_probe.calls == []


def test_read_only_materialized_verifier_rejects_tampered_installed_bundle_first(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    prepared = _prepared(tmp_path, monkeypatch)
    observer = _BoundaryObserver(prepared.boundary)
    image_probe = _StateProbe(
        {prepared.assets.containers[0].runtime_image: [IMAGE_EXACT_READY]}
    )

    def reject_tampered(_self, _component, _identity, **_kwargs):
        raise fail(
            "DEPLOYMENT_BUNDLE_INVALID",
            "Installed deployment bundle is invalid.",
        )

    monkeypatch.setattr(BundleStore, "verify_installed", reject_tampered)

    with pytest.raises(DeploymentError) as caught:
        verify_materialized_bulk_runtime(
            prepared.layout,
            BULK_IDENTITY,
            expected_daemon_uid=prepared.request.docker_daemon_uid,
            expected_daemon_gid=prepared.request.docker_daemon_gid,
            expected_client_identity=CLIENT_IDENTITY,
            expected_endpoint_identity=ENDPOINT_IDENTITY,
            installed_owner_uid=os.getuid(),
            installed_owner_gid=os.getgid(),
            _static_verifier=lambda _binding: Result.success(prepared.assets),
            _canary_verifier=lambda _binding: Result.success(prepared.assets),
            _boundary_observer=observer,
            _image_state_probe=image_probe,
        )

    assert caught.value.issue.code == "DEPLOYMENT_BUNDLE_INVALID"
    assert observer.calls == 0
    assert image_probe.calls == []


def test_verify_operation_returns_exact_receipt_without_reaching_load_runner(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    prepared = _prepared(tmp_path, monkeypatch)
    request = _request_for_operation(prepared, "verify")
    image = prepared.assets.containers[0].runtime_image
    runner = _DockerRunner()
    full = _FullVerifier(_not_ready())

    receipt = _materializer(
        prepared,
        runner,
        full,
        state_probe=_StateProbe({image: [IMAGE_EXACT_READY]}),
    ).prepare(request)

    assert receipt.request_identity == request.identity
    assert receipt.candidate_bulk_identity == request.candidate_bulk_identity
    assert receipt.image_count == 1
    assert runner.calls == []
    assert full.calls == 0


@pytest.mark.parametrize(
    ("state", "code"),
    (
        (IMAGE_MISSING, "BULK_RUNTIME_IMAGE_MISSING"),
        (IMAGE_INVALID, "BULK_RUNTIME_IMAGE_STATE_INVALID"),
    ),
)
def test_verify_operation_rejects_non_exact_images_without_reaching_load_runner(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    state: str,
    code: str,
) -> None:
    prepared = _prepared(tmp_path, monkeypatch)
    request = _request_for_operation(prepared, "verify")
    image = prepared.assets.containers[0].runtime_image
    runner = _DockerRunner()
    full = _FullVerifier(_not_ready())

    with pytest.raises(DeploymentError) as caught:
        _materializer(
            prepared,
            runner,
            full,
            state_probe=_StateProbe({image: [state]}),
        ).prepare(request)

    assert caught.value.issue.code == code
    assert runner.calls == []
    assert full.calls == 0


def test_verify_operation_is_bound_into_request_identity(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    prepared = _prepared(tmp_path, monkeypatch)
    activate = prepared.request
    verify = _request_for_operation(prepared, "verify")

    assert verify.identity != activate.identity
    assert BulkRuntimePrepareRequest.from_dict(verify.to_dict()) == verify

    substituted = verify.to_dict()
    substituted["operation"] = "activate"
    with pytest.raises(DeploymentError):
        BulkRuntimePrepareRequest.from_dict(substituted)

    forged = replace(verify, identity=activate.identity)
    with pytest.raises(DeploymentError) as caught:
        _materializer(
            prepared,
            _DockerRunner(),
            _FullVerifier(_not_ready()),
        ).prepare(forged)
    assert caught.value.issue.code == "BULK_RUNTIME_PREPARE_REQUEST_INVALID"
