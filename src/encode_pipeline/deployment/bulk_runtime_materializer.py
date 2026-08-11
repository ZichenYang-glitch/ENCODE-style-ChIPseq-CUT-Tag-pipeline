"""Offline materialization of one admitted Bulk RNA-seq Docker image set.

The request selects only content identities.  All filesystem and Docker
coordinates are fixed by the supported deployment layout.  Archive bytes are
fed through already verified file descriptors, and this module never pulls,
removes, prunes, or otherwise mutates Docker state beyond ``image load``.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
import hashlib
import json
import math
import os
from pathlib import Path, PurePosixPath
import re
import socket
import stat
import subprocess
import time
from typing import Any, Literal, Protocol
from urllib.parse import quote

from encode_pipeline.adapters.bulk_rnaseq.runtime_assets import (
    RuntimeAssetBinding,
    VerifiedContainerAsset,
    VerifiedRuntimeAssets,
    verify_runtime_asset_canary,
    verify_runtime_asset_closure,
    verify_runtime_assets,
)
from encode_pipeline.deployment.bundle import BundleStore
from encode_pipeline.deployment.bulk_docker_boundary import (
    BulkDockerBoundary,
    DOCKER_EXECUTABLE,
    DOCKER_SOCKET,
    observe_bulk_docker_boundary,
)
from encode_pipeline.deployment.canonical import canonical_identity
from encode_pipeline.deployment.errors import DeploymentError, fail
from encode_pipeline.deployment.layout import DeploymentLayout
from encode_pipeline.deployment.models import BULK_RNASEQ_RUNTIME
from encode_pipeline.deployment.native_contracts import (
    BULK_RUNTIME_ROOT_PATH,
    ProductionNativeContractResolver,
)
from encode_pipeline.deployment.operator_action import (
    BulkRuntimePrepareReceipt,
    BulkRuntimePrepareRequest,
)


DOCKER_ENDPOINT = f"unix://{DOCKER_SOCKET}"
MATERIALIZATION_TIMEOUT_SECONDS = 14_400.0

_MAX_IMAGES = 4_096
_MAX_DOCKER_CLIENT_BYTES = 256 * 1024 * 1024
_MAX_DOCKER_RESPONSE_BYTES = 8 * 1024 * 1024
_MAX_PROC_NET_UNIX_BYTES = 4 * 1024 * 1024
_DOCKER_PROBE_TIMEOUT_SECONDS = 15.0
_READ_CHUNK_BYTES = 1024 * 1024
_DOCKER_ENVIRONMENT = {
    "HOME": "/nonexistent",
    "LANG": "C.UTF-8",
    "LC_ALL": "C.UTF-8",
    "PATH": "/usr/bin:/bin",
}
_RUNTIME_IDENTITY_SCHEME = "helixweave-bulk-runtime-static-identity-v1"
_IMAGE_SET_IDENTITY_SCHEME = "helixweave-bulk-runtime-image-set-identity-v1"
_PROC_NET_UNIX = Path("/proc/net/unix")
_DIGEST = re.compile(r"^sha256:[0-9a-f]{64}$")
_IDENTITY = re.compile(r"^sha256-[0-9a-f]{64}$")

IMAGE_EXACT_READY = "exact-ready"
IMAGE_MISSING = "missing"
IMAGE_INVALID = "invalid/unavailable"
DockerImageState = Literal["exact-ready", "missing", "invalid/unavailable"]


class DockerCommandRunner(Protocol):
    """The exact ``subprocess.run`` surface used for local image loading."""

    def __call__(
        self,
        argv: Sequence[str],
        *,
        stdin: int,
        stdout: int,
        stderr: int,
        cwd: Path,
        env: Mapping[str, str],
        close_fds: bool,
        timeout: float,
        check: bool,
        shell: bool,
        umask: int,
    ) -> Any: ...


@dataclass(frozen=True)
class _ArchivePlan:
    path: Path
    sha256: str
    size: int


@dataclass(frozen=True)
class _OpenArchive:
    plan: _ArchivePlan
    descriptor: int
    witness: tuple[int, ...]


@dataclass(frozen=True)
class _ImagePlan:
    image: str
    rootfs_diff_ids: tuple[str, ...]
    archive: _ArchivePlan


@dataclass(frozen=True)
class MaterializedBulkRuntimeVerification:
    """Path-free evidence that every admitted image is exactly materialized."""

    runtime_identity: str
    image_set_identity: str
    image_count: int
    docker_client_identity: str
    docker_endpoint_identity: str


def probe_bulk_docker_image_state(
    binding: RuntimeAssetBinding,
    image: str,
    rootfs_diff_ids: tuple[str, ...],
) -> DockerImageState:
    """Return the exact state of one image through the fixed AF_UNIX daemon."""

    if (
        not isinstance(binding, RuntimeAssetBinding)
        or binding.docker_executable != DOCKER_EXECUTABLE
        or binding.docker_socket != DOCKER_SOCKET
        or not _valid_digest(image)
        or not isinstance(rootfs_diff_ids, tuple)
        or not all(_valid_digest(value) for value in rootfs_diff_ids)
    ):
        return IMAGE_INVALID
    target = f"/images/{quote(image, safe='')}/json"
    request = (
        f"GET {target} HTTP/1.0\r\nHost: docker\r\nConnection: close\r\n\r\n"
    ).encode("ascii")
    content = bytearray()
    try:
        with socket.socket(socket.AF_UNIX, socket.SOCK_STREAM) as client:
            client.settimeout(_DOCKER_PROBE_TIMEOUT_SECONDS)
            client.connect(str(DOCKER_SOCKET))
            client.sendall(request)
            while chunk := client.recv(
                min(
                    _READ_CHUNK_BYTES,
                    _MAX_DOCKER_RESPONSE_BYTES + 1 - len(content),
                )
            ):
                content.extend(chunk)
                if len(content) > _MAX_DOCKER_RESPONSE_BYTES:
                    return IMAGE_INVALID
    except (OSError, TimeoutError):
        return IMAGE_INVALID
    try:
        status, body = _parse_docker_http_response(bytes(content))
        document = json.loads(body, object_pairs_hook=_unique_object)
    except (UnicodeDecodeError, ValueError, json.JSONDecodeError):
        return IMAGE_INVALID
    if status == 404:
        if document == {"message": f"No such image: {image}"}:
            return IMAGE_MISSING
        return IMAGE_INVALID
    if status != 200 or not isinstance(document, Mapping):
        return IMAGE_INVALID
    rootfs = document.get("RootFS")
    if (
        document.get("Id") != image
        or not isinstance(rootfs, Mapping)
        or rootfs.get("Type") != "layers"
        or rootfs.get("Layers") != list(rootfs_diff_ids)
    ):
        return IMAGE_INVALID
    return IMAGE_EXACT_READY


def verify_materialized_bulk_runtime(
    layout: DeploymentLayout,
    deployment_identity: str,
    *,
    expected_daemon_uid: int,
    expected_daemon_gid: int,
    expected_client_identity: str,
    expected_endpoint_identity: str,
    installed_owner_uid: int = 0,
    installed_owner_gid: int = 0,
    _bundle_store: BundleStore | None = None,
    _native_resolver: ProductionNativeContractResolver | None = None,
    _static_verifier: Callable[[RuntimeAssetBinding], object] = (
        verify_runtime_asset_closure
    ),
    _canary_verifier: Callable[[RuntimeAssetBinding], object] = (
        verify_runtime_asset_canary
    ),
    _boundary_observer: Callable[[int, int], BulkDockerBoundary] = (
        observe_bulk_docker_boundary
    ),
    _image_state_probe: Callable[
        [RuntimeAssetBinding, str, tuple[str, ...]], DockerImageState
    ] = probe_bulk_docker_image_state,
) -> MaterializedBulkRuntimeVerification:
    """Re-admit and observe one exact image set without mutating Docker state."""

    if (
        not isinstance(layout, DeploymentLayout)
        or not isinstance(deployment_identity, str)
        or _IDENTITY.fullmatch(deployment_identity) is None
        or not isinstance(expected_client_identity, str)
        or _IDENTITY.fullmatch(expected_client_identity) is None
        or not isinstance(expected_endpoint_identity, str)
        or _IDENTITY.fullmatch(expected_endpoint_identity) is None
        or any(
            not isinstance(value, int)
            or isinstance(value, bool)
            or not 0 <= value <= 2**32 - 1
            for value in (
                expected_daemon_uid,
                expected_daemon_gid,
                installed_owner_uid,
                installed_owner_gid,
            )
        )
        or (_bundle_store is not None and not isinstance(_bundle_store, BundleStore))
        or (
            _native_resolver is not None
            and not isinstance(_native_resolver, ProductionNativeContractResolver)
        )
        or not callable(_static_verifier)
        or not callable(_canary_verifier)
        or not callable(_boundary_observer)
        or not callable(_image_state_probe)
    ):
        raise fail(
            "BULK_RUNTIME_VERIFICATION_INVALID",
            "Bulk runtime verification request is invalid.",
            component=BULK_RNASEQ_RUNTIME,
        )
    bundle_store = BundleStore(layout) if _bundle_store is None else _bundle_store
    native_resolver = (
        ProductionNativeContractResolver()
        if _native_resolver is None
        else _native_resolver
    )
    binding, static = _admit_bulk_source(
        layout,
        deployment_identity,
        installed_owner_uid=installed_owner_uid,
        installed_owner_gid=installed_owner_gid,
        bundle_store=bundle_store,
        native_resolver=native_resolver,
        static_verifier=_static_verifier,
    )
    expected_boundary = _observe_expected_boundary(
        _boundary_observer,
        daemon_uid=expected_daemon_uid,
        daemon_gid=expected_daemon_gid,
        client_identity=expected_client_identity,
        endpoint_identity=expected_endpoint_identity,
    )
    try:
        canary_result = _canary_verifier(binding)
    except Exception:
        canary_result = None
    _observe_expected_boundary(
        _boundary_observer,
        daemon_uid=expected_daemon_uid,
        daemon_gid=expected_daemon_gid,
        client_identity=expected_client_identity,
        endpoint_identity=expected_endpoint_identity,
        expected=expected_boundary,
    )
    canary = _verified_value(
        canary_result,
        code="BULK_RUNTIME_CANARY_INVALID",
    )
    if canary != static:
        raise _canary_invalid()

    runtime_identity, image_set_identity, image_count = _runtime_evidence(static)
    for image in _image_plans(static):
        _observe_expected_boundary(
            _boundary_observer,
            daemon_uid=expected_daemon_uid,
            daemon_gid=expected_daemon_gid,
            client_identity=expected_client_identity,
            endpoint_identity=expected_endpoint_identity,
            expected=expected_boundary,
        )
        try:
            state = _image_state_probe(
                binding,
                image.image,
                image.rootfs_diff_ids,
            )
        except Exception:
            state = IMAGE_INVALID
        _observe_expected_boundary(
            _boundary_observer,
            daemon_uid=expected_daemon_uid,
            daemon_gid=expected_daemon_gid,
            client_identity=expected_client_identity,
            endpoint_identity=expected_endpoint_identity,
            expected=expected_boundary,
        )
        if state == IMAGE_MISSING:
            raise _image_missing()
        if state != IMAGE_EXACT_READY:
            raise _image_state_invalid()
    _observe_expected_boundary(
        _boundary_observer,
        daemon_uid=expected_daemon_uid,
        daemon_gid=expected_daemon_gid,
        client_identity=expected_client_identity,
        endpoint_identity=expected_endpoint_identity,
        expected=expected_boundary,
    )
    return MaterializedBulkRuntimeVerification(
        runtime_identity=runtime_identity,
        image_set_identity=image_set_identity,
        image_count=image_count,
        docker_client_identity=expected_boundary.client_identity,
        docker_endpoint_identity=expected_boundary.endpoint_identity,
    )


class OfflineBulkRuntimeMaterializer:
    """Load only verified local archives and attest the resulting image set."""

    def __init__(
        self,
        layout: DeploymentLayout,
        *,
        installed_owner_uid: int = 0,
        installed_owner_gid: int = 0,
        timeout_seconds: float = MATERIALIZATION_TIMEOUT_SECONDS,
        command_runner: DockerCommandRunner = subprocess.run,
        monotonic_clock: Callable[[], float] = time.monotonic,
        bundle_store: BundleStore | None = None,
        native_resolver: ProductionNativeContractResolver | None = None,
        static_verifier: Callable[[RuntimeAssetBinding], object] = (
            verify_runtime_asset_closure
        ),
        canary_verifier: Callable[[RuntimeAssetBinding], object] = (
            verify_runtime_asset_canary
        ),
        full_verifier: Callable[[RuntimeAssetBinding], object] = verify_runtime_assets,
        boundary_observer: Callable[[int, int], BulkDockerBoundary] = (
            observe_bulk_docker_boundary
        ),
        image_state_probe: Callable[
            [RuntimeAssetBinding, str, tuple[str, ...]], DockerImageState
        ] = probe_bulk_docker_image_state,
    ) -> None:
        if (
            not isinstance(layout, DeploymentLayout)
            or any(
                not isinstance(value, int) or isinstance(value, bool) or value < 0
                for value in (installed_owner_uid, installed_owner_gid)
            )
            or not isinstance(timeout_seconds, (int, float))
            or isinstance(timeout_seconds, bool)
            or not math.isfinite(float(timeout_seconds))
            or not 1.0 <= float(timeout_seconds) <= MATERIALIZATION_TIMEOUT_SECONDS
            or not callable(command_runner)
            or not callable(monotonic_clock)
            or not callable(static_verifier)
            or not callable(canary_verifier)
            or not callable(full_verifier)
            or not callable(boundary_observer)
            or not callable(image_state_probe)
            or (bundle_store is not None and not isinstance(bundle_store, BundleStore))
            or (
                native_resolver is not None
                and not isinstance(native_resolver, ProductionNativeContractResolver)
            )
        ):
            raise fail(
                "BULK_RUNTIME_MATERIALIZER_INVALID",
                "Bulk runtime materializer is invalid.",
                component=BULK_RNASEQ_RUNTIME,
            )
        self._layout = layout
        self._installed_owner_uid = installed_owner_uid
        self._installed_owner_gid = installed_owner_gid
        self._timeout_seconds = float(timeout_seconds)
        self._command_runner = command_runner
        self._monotonic_clock = monotonic_clock
        self._bundle_store = (
            BundleStore(layout) if bundle_store is None else bundle_store
        )
        self._native_resolver = (
            ProductionNativeContractResolver()
            if native_resolver is None
            else native_resolver
        )
        self._static_verifier = static_verifier
        self._canary_verifier = canary_verifier
        self._full_verifier = full_verifier
        self._boundary_observer = boundary_observer
        self._image_state_probe = image_state_probe

    def prepare(self, request: BulkRuntimePrepareRequest) -> BulkRuntimePrepareReceipt:
        """Make the fixed candidate image set available without network access."""

        if not isinstance(request, BulkRuntimePrepareRequest):
            raise fail(
                "BULK_RUNTIME_PREPARE_REQUEST_INVALID",
                "Bulk runtime preparation request is invalid.",
                component=BULK_RNASEQ_RUNTIME,
            )
        try:
            if BulkRuntimePrepareRequest.from_dict(request.to_dict()) != request:
                raise ValueError
        except (DeploymentError, TypeError, ValueError):
            raise fail(
                "BULK_RUNTIME_PREPARE_REQUEST_INVALID",
                "Bulk runtime preparation request is invalid.",
                component=BULK_RNASEQ_RUNTIME,
            ) from None
        if request.operation == "verify":
            verified = verify_materialized_bulk_runtime(
                self._layout,
                request.candidate_bulk_identity,
                expected_daemon_uid=request.docker_daemon_uid,
                expected_daemon_gid=request.docker_daemon_gid,
                expected_client_identity=request.docker_client_identity,
                expected_endpoint_identity=request.docker_endpoint_identity,
                installed_owner_uid=self._installed_owner_uid,
                installed_owner_gid=self._installed_owner_gid,
                _bundle_store=self._bundle_store,
                _native_resolver=self._native_resolver,
                _static_verifier=self._static_verifier,
                _canary_verifier=self._canary_verifier,
                _boundary_observer=self._boundary_observer,
                _image_state_probe=self._image_state_probe,
            )
            return _receipt(
                request,
                runtime_identity=verified.runtime_identity,
                image_set_identity=verified.image_set_identity,
                image_count=verified.image_count,
            )
        if request.operation not in {"activate", "rollback"}:
            raise fail(
                "BULK_RUNTIME_PREPARE_REQUEST_INVALID",
                "Bulk runtime preparation request is invalid.",
                component=BULK_RNASEQ_RUNTIME,
            )
        deadline = self._monotonic_clock() + self._timeout_seconds
        binding, static = self._admit_source(request)
        self._require_before_deadline(deadline)
        runtime_identity, image_set_identity, image_count = _runtime_evidence(static)
        expected_boundary = self._require_boundary(request)
        canary = _verified_value(
            self._canary_verify_with_boundary(
                binding,
                request=request,
                expected_boundary=expected_boundary,
            ),
            code="BULK_RUNTIME_CANARY_INVALID",
        )
        if canary != static:
            raise _canary_invalid()
        self._require_before_deadline(deadline)
        image_plans = _image_plans(static)
        missing: list[_ImagePlan] = []
        invalid = False
        for image in image_plans:
            self._require_before_deadline(deadline)
            state = self._probe_image_state(
                binding,
                image,
                request=request,
                expected_boundary=expected_boundary,
            )
            self._require_before_deadline(deadline)
            if state == IMAGE_MISSING:
                missing.append(image)
            elif state != IMAGE_EXACT_READY:
                invalid = True
        if invalid:
            raise _image_state_invalid()

        grouped_missing = _group_missing_images(missing)
        opened: list[_OpenArchive] = []
        try:
            # Open and hash every archive before the first Docker side effect.
            for archive, _images in grouped_missing:
                opened.append(self._open_archive(archive, expected_root=static.root))
                self._require_before_deadline(deadline)
            for opened_archive, (_archive, images) in zip(
                opened, grouped_missing, strict=True
            ):
                self._load_with_boundary(
                    opened_archive,
                    cwd=static.root,
                    deadline=deadline,
                    request=request,
                    expected_boundary=expected_boundary,
                )
                for image in images:
                    self._require_before_deadline(deadline)
                    if (
                        self._probe_image_state(
                            binding,
                            image,
                            request=request,
                            expected_boundary=expected_boundary,
                        )
                        != IMAGE_EXACT_READY
                    ):
                        raise _post_load_invalid()
                    self._require_before_deadline(deadline)
        finally:
            for archive in opened:
                try:
                    os.close(archive.descriptor)
                except OSError:
                    pass

        final = self._full_verify_with_boundary(
            binding,
            request=request,
            expected_boundary=expected_boundary,
        )
        self._require_before_deadline(deadline)
        full = _verified_value(final, code="BULK_RUNTIME_POST_LOAD_INVALID")
        if _runtime_evidence(full) != (
            runtime_identity,
            image_set_identity,
            image_count,
        ):
            raise _post_load_invalid()
        self._require_boundary(request, expected=expected_boundary)
        return _receipt(
            request,
            runtime_identity=runtime_identity,
            image_set_identity=image_set_identity,
            image_count=image_count,
        )

    def _admit_source(
        self, request: BulkRuntimePrepareRequest
    ) -> tuple[RuntimeAssetBinding, VerifiedRuntimeAssets]:
        return _admit_bulk_source(
            self._layout,
            request.candidate_bulk_identity,
            installed_owner_uid=self._installed_owner_uid,
            installed_owner_gid=self._installed_owner_gid,
            bundle_store=self._bundle_store,
            native_resolver=self._native_resolver,
            static_verifier=self._static_verifier,
        )

    def _require_boundary(
        self,
        request: BulkRuntimePrepareRequest,
        *,
        expected: BulkDockerBoundary | None = None,
    ) -> BulkDockerBoundary:
        try:
            observed = self._boundary_observer(
                request.docker_daemon_uid,
                request.docker_daemon_gid,
            )
        except Exception:
            observed = None
        if expected is not None and observed != expected:
            raise fail(
                "BULK_RUNTIME_DOCKER_BOUNDARY_CHANGED",
                "Bulk runtime Docker boundary changed during preparation.",
                component=BULK_RNASEQ_RUNTIME,
                recoverable=True,
            )
        if (
            not isinstance(observed, BulkDockerBoundary)
            or observed.client_identity != request.docker_client_identity
            or observed.endpoint_identity != request.docker_endpoint_identity
            or observed.daemon_uid != request.docker_daemon_uid
            or observed.daemon_gid != request.docker_daemon_gid
        ):
            raise _docker_boundary_invalid()
        return observed

    def _probe_image_state(
        self,
        binding: RuntimeAssetBinding,
        image: _ImagePlan,
        *,
        request: BulkRuntimePrepareRequest,
        expected_boundary: BulkDockerBoundary,
    ) -> DockerImageState:
        self._require_boundary(request, expected=expected_boundary)
        try:
            state = self._image_state_probe(
                binding,
                image.image,
                image.rootfs_diff_ids,
            )
        except Exception:
            state = IMAGE_INVALID
        self._require_boundary(request, expected=expected_boundary)
        if state not in {IMAGE_EXACT_READY, IMAGE_MISSING, IMAGE_INVALID}:
            return IMAGE_INVALID
        return state

    def _load_with_boundary(
        self,
        archive: _OpenArchive,
        *,
        cwd: Path,
        deadline: float,
        request: BulkRuntimePrepareRequest,
        expected_boundary: BulkDockerBoundary,
    ) -> None:
        self._require_boundary(request, expected=expected_boundary)
        try:
            self._load_archive(archive, cwd=cwd, deadline=deadline)
        finally:
            self._require_boundary(request, expected=expected_boundary)

    def _full_verify_with_boundary(
        self,
        binding: RuntimeAssetBinding,
        *,
        request: BulkRuntimePrepareRequest,
        expected_boundary: BulkDockerBoundary,
    ) -> object:
        self._require_boundary(request, expected=expected_boundary)
        try:
            result = self._full_verifier(binding)
        except Exception:
            result = None
        self._require_boundary(request, expected=expected_boundary)
        if result is None:
            raise _post_load_invalid()
        return result

    def _canary_verify_with_boundary(
        self,
        binding: RuntimeAssetBinding,
        *,
        request: BulkRuntimePrepareRequest,
        expected_boundary: BulkDockerBoundary,
    ) -> object:
        self._require_boundary(request, expected=expected_boundary)
        try:
            result = self._canary_verifier(binding)
        except Exception:
            result = None
        self._require_boundary(request, expected=expected_boundary)
        if result is None:
            raise _canary_invalid()
        return result

    def _open_archive(
        self,
        plan: _ArchivePlan,
        *,
        expected_root: Path,
    ) -> _OpenArchive:
        descriptor = -1
        try:
            plan.path.relative_to(expected_root)
            descriptor = os.open(
                plan.path,
                os.O_RDONLY
                | getattr(os, "O_CLOEXEC", 0)
                | getattr(os, "O_NOFOLLOW", 0),
            )
            observed = os.fstat(descriptor)
            if (
                not stat.S_ISREG(observed.st_mode)
                or stat.S_IMODE(observed.st_mode) != 0o444
                or observed.st_nlink != 1
                or observed.st_uid != self._installed_owner_uid
                or observed.st_gid != self._installed_owner_gid
                or observed.st_size != plan.size
            ):
                raise OSError
            witness = _stat_witness(observed)
            if _hash_descriptor(descriptor, expected_size=plan.size) != plan.sha256:
                raise OSError
            after = os.fstat(descriptor)
            if _stat_witness(after) != witness:
                raise OSError
            return _OpenArchive(plan, descriptor, witness)
        except (OSError, ValueError):
            if descriptor >= 0:
                os.close(descriptor)
            raise fail(
                "BULK_RUNTIME_ARCHIVE_INVALID",
                "Bulk runtime archive is invalid.",
                component=BULK_RNASEQ_RUNTIME,
            ) from None

    def _load_archive(
        self,
        archive: _OpenArchive,
        *,
        cwd: Path,
        deadline: float,
    ) -> None:
        self._verify_open_archive(archive)
        remaining = self._remaining(deadline)
        try:
            completed = self._command_runner(
                (
                    str(DOCKER_EXECUTABLE),
                    "--host",
                    DOCKER_ENDPOINT,
                    "image",
                    "load",
                ),
                stdin=archive.descriptor,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                cwd=cwd,
                env=dict(_DOCKER_ENVIRONMENT),
                close_fds=True,
                timeout=remaining,
                check=False,
                shell=False,
                umask=0o077,
            )
        except (OSError, subprocess.SubprocessError):
            raise _image_load_failed() from None
        if getattr(completed, "returncode", None) != 0:
            raise _image_load_failed()
        self._verify_open_archive(archive)
        self._require_before_deadline(deadline)

    def _verify_open_archive(self, archive: _OpenArchive) -> None:
        try:
            before = os.fstat(archive.descriptor)
            if (
                _stat_witness(before) != archive.witness
                or before.st_size != archive.plan.size
                or _hash_descriptor(archive.descriptor, expected_size=archive.plan.size)
                != archive.plan.sha256
                or _stat_witness(os.fstat(archive.descriptor)) != archive.witness
            ):
                raise OSError
        except OSError:
            raise fail(
                "BULK_RUNTIME_ARCHIVE_CHANGED",
                "Bulk runtime archive changed during preparation.",
                component=BULK_RNASEQ_RUNTIME,
                recoverable=True,
            ) from None

    def _remaining(self, deadline: float) -> float:
        remaining = deadline - self._monotonic_clock()
        if not math.isfinite(remaining) or remaining <= 0:
            raise fail(
                "BULK_RUNTIME_PREPARE_TIMEOUT",
                "Bulk runtime preparation timed out.",
                component=BULK_RNASEQ_RUNTIME,
                recoverable=True,
            )
        return min(remaining, self._timeout_seconds)

    def _require_before_deadline(self, deadline: float) -> None:
        self._remaining(deadline)


def _admit_bulk_source(
    layout: DeploymentLayout,
    deployment_identity: str,
    *,
    installed_owner_uid: int,
    installed_owner_gid: int,
    bundle_store: BundleStore,
    native_resolver: ProductionNativeContractResolver,
    static_verifier: Callable[[RuntimeAssetBinding], object],
) -> tuple[RuntimeAssetBinding, VerifiedRuntimeAssets]:
    source_root = layout.component_store(BULK_RNASEQ_RUNTIME) / deployment_identity
    try:
        manifest = bundle_store.verify_installed(
            BULK_RNASEQ_RUNTIME,
            deployment_identity,
            expected_owner_uid=installed_owner_uid,
            expected_owner_gid=installed_owner_gid,
        )
        facts = native_resolver.resolve(source_root, manifest)
        if (
            facts.component != BULK_RNASEQ_RUNTIME
            or facts.deployment_identity != deployment_identity
            or facts.contracts != manifest.provides
        ):
            raise ValueError
        binding = RuntimeAssetBinding(
            root=source_root.joinpath(*PurePosixPath(BULK_RUNTIME_ROOT_PATH).parts),
            docker_executable=DOCKER_EXECUTABLE,
            docker_socket=DOCKER_SOCKET,
        )
        static = _verified_value(
            static_verifier(binding),
            code="BULK_RUNTIME_STATIC_CLOSURE_INVALID",
        )
        if static.root != binding.root:
            raise ValueError
        return binding, static
    except DeploymentError:
        raise
    except Exception:
        raise fail(
            "BULK_RUNTIME_SOURCE_INVALID",
            "Bulk runtime source is invalid.",
            component=BULK_RNASEQ_RUNTIME,
        ) from None


def _observe_expected_boundary(
    observer: Callable[[int, int], BulkDockerBoundary],
    *,
    daemon_uid: int,
    daemon_gid: int,
    client_identity: str,
    endpoint_identity: str,
    expected: BulkDockerBoundary | None = None,
) -> BulkDockerBoundary:
    try:
        observed = observer(daemon_uid, daemon_gid)
    except Exception:
        observed = None
    if expected is not None and observed != expected:
        raise fail(
            "BULK_RUNTIME_DOCKER_BOUNDARY_CHANGED",
            "Bulk runtime Docker boundary changed during verification.",
            component=BULK_RNASEQ_RUNTIME,
            recoverable=True,
        )
    if (
        not isinstance(observed, BulkDockerBoundary)
        or observed.client_identity != client_identity
        or observed.endpoint_identity != endpoint_identity
        or observed.daemon_uid != daemon_uid
        or observed.daemon_gid != daemon_gid
    ):
        raise _docker_boundary_invalid()
    return observed


def _image_plans(assets: VerifiedRuntimeAssets) -> tuple[_ImagePlan, ...]:
    by_image: dict[str, _ImagePlan] = {}
    archive_images: dict[tuple[str, int], str] = {}
    try:
        for container in assets.containers:
            container.local_asset.relative_to(assets.root)
            archive = _ArchivePlan(
                container.local_asset,
                container.local_sha256,
                container.size_bytes,
            )
            plan = _ImagePlan(
                container.runtime_image,
                container.rootfs_diff_ids,
                archive,
            )
            previous = by_image.setdefault(container.runtime_image, plan)
            if previous != plan:
                raise ValueError
            archive_coordinate = (archive.sha256, archive.size)
            prior_image = archive_images.setdefault(
                archive_coordinate, container.runtime_image
            )
            if prior_image != container.runtime_image:
                raise ValueError
    except (AttributeError, TypeError, ValueError):
        raise fail(
            "BULK_RUNTIME_STATIC_CLOSURE_INVALID",
            "Bulk runtime static closure is invalid.",
            component=BULK_RNASEQ_RUNTIME,
        ) from None
    plans = tuple(by_image[image] for image in sorted(by_image))
    if not 0 < len(plans) <= _MAX_IMAGES:
        raise fail(
            "BULK_RUNTIME_STATIC_CLOSURE_INVALID",
            "Bulk runtime static closure is invalid.",
            component=BULK_RNASEQ_RUNTIME,
        )
    return plans


def _group_missing_images(
    images: Sequence[_ImagePlan],
) -> tuple[tuple[_ArchivePlan, tuple[_ImagePlan, ...]], ...]:
    grouped: dict[_ArchivePlan, list[_ImagePlan]] = {}
    for image in images:
        grouped.setdefault(image.archive, []).append(image)
    return tuple(
        (archive, tuple(sorted(grouped[archive], key=lambda item: item.image)))
        for archive in sorted(
            grouped,
            key=lambda item: (item.sha256, item.size, str(item.path)),
        )
    )


def _runtime_evidence(assets: VerifiedRuntimeAssets) -> tuple[str, str, int]:
    if not isinstance(assets, VerifiedRuntimeAssets):
        raise _post_load_invalid()
    images: dict[str, dict[str, object]] = {}
    try:
        for container in assets.containers:
            evidence = _container_evidence(container)
            previous = images.setdefault(container.runtime_image, evidence)
            if previous != evidence:
                raise ValueError
    except (AttributeError, TypeError, ValueError):
        raise _post_load_invalid() from None
    ordered = [images[name] for name in sorted(images)]
    if not 0 < len(ordered) <= _MAX_IMAGES:
        raise _post_load_invalid()
    image_set_identity = canonical_identity(
        {"images": ordered}, scheme=_IMAGE_SET_IDENTITY_SCHEME
    )
    runtime_identity = canonical_identity(
        {
            "runtime_identity_sha256": assets.runtime_identity_sha256,
            "source_tree_sha256": assets.source_tree_sha256,
            "nextflow_sha256": assets.nextflow_sha256,
            "jdk_archive_sha256": assets.jdk_archive_sha256,
            "jdk_tree_sha256": assets.jdk_tree_sha256,
            "java_executable_sha256": assets.java_executable_sha256,
            "plugin_archive_sha256": assets.plugin_archive_sha256,
            "plugin_tree_sha256": assets.plugin_tree_sha256,
            "container_inventory_sha256": assets.container_inventory_sha256,
            "container_lock_sha256": assets.container_lock_sha256,
            "container_process_audit_sha256": (assets.container_process_audit_sha256),
            "network_isolation_executable_sha256": (
                assets.network_isolation_executable_sha256
            ),
            "network_isolation_version_output_sha256": (
                assets.network_isolation_version_output_sha256
            ),
            "image_set_identity": image_set_identity,
        },
        scheme=_RUNTIME_IDENTITY_SCHEME,
    )
    return runtime_identity, image_set_identity, len(ordered)


def _container_evidence(container: VerifiedContainerAsset) -> dict[str, object]:
    return {
        "image": container.image,
        "oci_digest": container.oci_digest,
        "archive_sha256": container.local_sha256,
        "archive_size": container.size_bytes,
        "distribution_manifest_sha256": container.distribution_manifest_sha256,
        "config_digest": container.config_digest,
        "runtime_image": container.runtime_image,
        "rootfs_diff_ids": list(container.rootfs_diff_ids),
    }


def _verified_value(result: object, *, code: str) -> VerifiedRuntimeAssets:
    if not _result_succeeded(result) or not isinstance(
        getattr(result, "value", None), VerifiedRuntimeAssets
    ):
        raise fail(
            code,
            "Bulk runtime verification failed.",
            component=BULK_RNASEQ_RUNTIME,
            recoverable=code
            in {"BULK_RUNTIME_CANARY_INVALID", "BULK_RUNTIME_POST_LOAD_INVALID"},
        )
    return result.value


def _result_succeeded(result: object) -> bool:
    return getattr(result, "is_success", None) is True


def _receipt(
    request: BulkRuntimePrepareRequest,
    *,
    runtime_identity: str,
    image_set_identity: str,
    image_count: int,
) -> BulkRuntimePrepareReceipt:
    return BulkRuntimePrepareReceipt.create(
        request_identity=request.identity,
        candidate_bulk_identity=request.candidate_bulk_identity,
        runtime_identity=runtime_identity,
        image_set_identity=image_set_identity,
        image_count=image_count,
    )


def _hash_descriptor(descriptor: int, *, expected_size: int) -> str:
    os.lseek(descriptor, 0, os.SEEK_SET)
    digest = hashlib.sha256()
    total = 0
    while chunk := os.read(descriptor, _READ_CHUNK_BYTES):
        total += len(chunk)
        if total > expected_size:
            raise OSError
        digest.update(chunk)
    if total != expected_size:
        raise OSError
    os.lseek(descriptor, 0, os.SEEK_SET)
    return digest.hexdigest()


def _read_proc_net_unix() -> bytes:
    descriptor = -1
    try:
        descriptor = os.open(
            _PROC_NET_UNIX,
            os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0),
        )
        before = os.fstat(descriptor)
        if not stat.S_ISREG(before.st_mode) or before.st_uid != 0:
            raise OSError
        content = bytearray()
        while chunk := os.read(
            descriptor,
            min(
                _READ_CHUNK_BYTES,
                _MAX_PROC_NET_UNIX_BYTES + 1 - len(content),
            ),
        ):
            content.extend(chunk)
            if len(content) > _MAX_PROC_NET_UNIX_BYTES:
                raise OSError
        after = os.fstat(descriptor)
        if _stat_witness(before) != _stat_witness(after):
            raise OSError
        return bytes(content)
    except OSError:
        raise _docker_boundary_invalid() from None
    finally:
        if descriptor >= 0:
            os.close(descriptor)


def _socket_kernel_inode(content: bytes, socket_path: Path) -> int:
    try:
        lines = content.decode("ascii", errors="strict").splitlines()
    except UnicodeDecodeError:
        raise ValueError from None
    if not lines or lines[0].split() != [
        "Num",
        "RefCount",
        "Protocol",
        "Flags",
        "Type",
        "St",
        "Inode",
        "Path",
    ]:
        raise ValueError
    matches: list[int] = []
    for line in lines[1:]:
        fields = line.split(maxsplit=7)
        if len(fields) != 8 or fields[7] != str(socket_path):
            continue
        if fields[3] != "00010000" or fields[4] != "0001" or fields[5] != "01":
            raise ValueError
        if not fields[6].isdigit() or int(fields[6]) <= 0:
            raise ValueError
        matches.append(int(fields[6]))
    if len(matches) != 1:
        raise ValueError
    return matches[0]


def _parse_docker_http_response(content: bytes) -> tuple[int, bytes]:
    if not content or len(content) > _MAX_DOCKER_RESPONSE_BYTES:
        raise ValueError
    header, separator, body = content.partition(b"\r\n\r\n")
    if not separator or len(header) > 64 * 1024:
        raise ValueError
    lines = header.split(b"\r\n")
    if not lines:
        raise ValueError
    status_fields = lines[0].split(b" ", 2)
    if (
        len(status_fields) != 3
        or status_fields[0] not in {b"HTTP/1.0", b"HTTP/1.1"}
        or len(status_fields[1]) != 3
        or not status_fields[1].isdigit()
    ):
        raise ValueError
    headers: dict[bytes, bytes] = {}
    for line in lines[1:]:
        name, colon, value = line.partition(b":")
        name = name.strip().lower()
        value = value.strip()
        if (
            not colon
            or not name
            or name in headers
            or any(character < 0x21 or character > 0x7E for character in name)
        ):
            raise ValueError
        headers[name] = value
    if b"transfer-encoding" in headers or b"content-encoding" in headers:
        raise ValueError
    length = headers.get(b"content-length")
    if length is not None and (not length.isdigit() or int(length) != len(body)):
        raise ValueError
    return int(status_fields[1]), body


def _unique_object(pairs: list[tuple[str, object]]) -> dict[str, object]:
    value: dict[str, object] = {}
    for key, item in pairs:
        if key in value:
            raise ValueError
        value[key] = item
    return value


def _valid_digest(value: object) -> bool:
    return isinstance(value, str) and _DIGEST.fullmatch(value) is not None


def _stat_witness(value: os.stat_result) -> tuple[int, ...]:
    return (
        value.st_dev,
        value.st_ino,
        value.st_mode,
        value.st_nlink,
        value.st_uid,
        value.st_gid,
        value.st_size,
        value.st_mtime_ns,
        value.st_ctime_ns,
    )


def _image_load_failed() -> DeploymentError:
    return fail(
        "BULK_RUNTIME_IMAGE_LOAD_FAILED",
        "Bulk runtime image load failed.",
        component=BULK_RNASEQ_RUNTIME,
        recoverable=True,
    )


def _image_state_invalid() -> DeploymentError:
    return fail(
        "BULK_RUNTIME_IMAGE_STATE_INVALID",
        "Bulk runtime image state is invalid or unavailable.",
        component=BULK_RNASEQ_RUNTIME,
        recoverable=True,
    )


def _image_missing() -> DeploymentError:
    return fail(
        "BULK_RUNTIME_IMAGE_MISSING",
        "Bulk runtime image is not materialized.",
        component=BULK_RNASEQ_RUNTIME,
        recoverable=True,
    )


def _docker_boundary_invalid() -> DeploymentError:
    return fail(
        "BULK_RUNTIME_DOCKER_BOUNDARY_INVALID",
        "Bulk runtime Docker boundary is invalid.",
        component=BULK_RNASEQ_RUNTIME,
        recoverable=True,
    )


def _post_load_invalid() -> DeploymentError:
    return fail(
        "BULK_RUNTIME_POST_LOAD_INVALID",
        "Bulk runtime is not ready after image loading.",
        component=BULK_RNASEQ_RUNTIME,
        recoverable=True,
    )


def _canary_invalid() -> DeploymentError:
    return fail(
        "BULK_RUNTIME_CANARY_INVALID",
        "Bulk runtime canary is invalid or unavailable.",
        component=BULK_RNASEQ_RUNTIME,
        recoverable=True,
    )


__all__ = [
    "BulkDockerBoundary",
    "DockerImageState",
    "DOCKER_ENDPOINT",
    "DOCKER_EXECUTABLE",
    "DOCKER_SOCKET",
    "IMAGE_EXACT_READY",
    "IMAGE_INVALID",
    "IMAGE_MISSING",
    "MATERIALIZATION_TIMEOUT_SECONDS",
    "MaterializedBulkRuntimeVerification",
    "OfflineBulkRuntimeMaterializer",
    "observe_bulk_docker_boundary",
    "probe_bulk_docker_image_state",
    "verify_materialized_bulk_runtime",
]
