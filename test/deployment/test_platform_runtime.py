from __future__ import annotations

import base64
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import runpy
import stat
import subprocess
from types import SimpleNamespace
import zipfile

import pytest

from encode_pipeline.deployment.admission import (
    ResolvedContractFacts,
    resolved_facts_compatibility,
)
from encode_pipeline.deployment.canonical import canonical_json_bytes
from encode_pipeline.deployment.database import fresh_database_candidate_path
from encode_pipeline.deployment.errors import DeploymentError, fail
from encode_pipeline.deployment.layout import DeploymentLayout
from encode_pipeline.deployment.models import ContractIdentity, ContractRequirement
from encode_pipeline.deployment.operator_action import (
    BulkRuntimePrepareReceipt,
    BulkRuntimePrepareRequest,
    EncodeRuntimeEntry,
    EncodeRuntimeInventory,
    EncodeRuntimePrepareReceipt,
    EncodeRuntimePrepareRequest,
)
import encode_pipeline.deployment.platform_runtime as runtime_module
from encode_pipeline.deployment.platform_runtime import (
    PLATFORM_RUNTIME_LAUNCHER_MEMBER,
    PLATFORM_RUNTIME_LAUNCHER_PATH,
    PLATFORM_RUNTIME_LOCK_PATH,
    PLATFORM_RUNTIME_SITE_PACKAGES,
    DatabasePrepareReceipt,
    DatabasePrepareRequest,
    DeploymentActionReceipt,
    DeploymentActionRequest,
    LockedWheel,
    PlatformWheelLock,
    ReadinessCheck,
    VERIFICATION_CHECKS,
    build_platform_wheel_lock,
    candidate_service_main,
    collect_platform_runtime_closure,
    inspect_platform_runtime_closure,
    prepare_candidate_database,
    verify_supported_python_runtime,
)


IDENTITY_A = f"sha256-{'a' * 64}"
IDENTITY_B = f"sha256-{'b' * 64}"
IDENTITY_C = f"sha256-{'c' * 64}"
TASK = f"task-{'c' * 32}"


def _record_hash(content: bytes) -> str:
    digest = base64.urlsafe_b64encode(hashlib.sha256(content).digest()).rstrip(b"=")
    return f"sha256={digest.decode('ascii')}"


def _build_wheel(
    wheelhouse: Path,
    *,
    name: str = "helixweave",
    version: str = "1.0.0",
    tag: str = "py3-none-any",
    requires_python: str = ">=3.12,<3.13",
    requirements: tuple[str, ...] = (),
    files: dict[str, bytes] | None = None,
) -> Path:
    normalized = name.replace("-", "_")
    filename = f"{normalized}-{version}-{tag}.whl"
    dist_info = f"{normalized}-{version}.dist-info"
    metadata_lines = [
        "Metadata-Version: 2.3",
        f"Name: {name}",
        f"Version: {version}",
        f"Requires-Python: {requires_python}",
        *[f"Requires-Dist: {item}" for item in requirements],
        "",
        "",
    ]
    members = {
        f"{dist_info}/METADATA": "\n".join(metadata_lines).encode(),
        f"{dist_info}/WHEEL": (
            "Wheel-Version: 1.0\n"
            "Generator: helixweave-test\n"
            "Root-Is-Purelib: true\n"
            f"Tag: {tag}\n\n"
        ).encode(),
        f"{normalized}/__init__.py": b"VALUE = 1\n",
        **(files or {}),
    }
    if name == "helixweave":
        members[PLATFORM_RUNTIME_LAUNCHER_MEMBER] = (
            Path(runtime_module.__file__).parent / "templates" / "helixweave-service"
        ).read_bytes()
    output = io.StringIO()
    writer = csv.writer(output, lineterminator="\n")
    for member, content in sorted(members.items()):
        writer.writerow((member, _record_hash(content), str(len(content))))
    record_name = f"{dist_info}/RECORD"
    writer.writerow((record_name, "", ""))
    members[record_name] = output.getvalue().encode()
    path = wheelhouse / filename
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for member, content in sorted(members.items()):
            info = zipfile.ZipInfo(member)
            info.external_attr = (stat.S_IFREG | 0o644) << 16
            archive.writestr(info, content)
    path.chmod(0o644)
    return path


def _write_lock(path: Path, wheels: tuple[Path, ...]) -> PlatformWheelLock:
    lock = PlatformWheelLock(
        tuple(
            sorted(
                LockedWheel(item.name, hashlib.sha256(item.read_bytes()).hexdigest())
                for item in wheels
            )
        )
    )
    path.write_bytes(lock.to_bytes())
    path.chmod(0o644)
    return lock


def _valid_closure(tmp_path: Path):
    wheelhouse = tmp_path / "wheelhouse"
    wheelhouse.mkdir()
    wheel = _build_wheel(wheelhouse)
    lock_path = tmp_path / "wheel-lock.json"
    lock = _write_lock(lock_path, (wheel,))
    destination = tmp_path / "closure"
    return (
        collect_platform_runtime_closure(wheelhouse, lock_path, destination),
        destination,
        lock,
    )


def test_lock_builder_is_offline_deterministic_and_validates_wheels(
    tmp_path: Path,
) -> None:
    wheelhouse = tmp_path / "wheelhouse"
    wheelhouse.mkdir()
    wheel = _build_wheel(wheelhouse)

    first = build_platform_wheel_lock(wheelhouse)
    second = build_platform_wheel_lock(wheelhouse)

    assert first == second
    assert first == _write_lock(tmp_path / "reviewed-lock.json", (wheel,))
    assert PlatformWheelLock.from_bytes(first.to_bytes()) == first


def test_lock_records_the_patch_neutral_supported_cpython_environment(
    tmp_path: Path,
) -> None:
    wheelhouse = tmp_path / "wheelhouse"
    wheelhouse.mkdir()
    _build_wheel(wheelhouse)

    lock = build_platform_wheel_lock(wheelhouse)

    assert lock.to_dict()["python"] == {
        "abi": "cp312",
        "executable": "/usr/bin/python3.12",
        "implementation_name": "cpython",
        "machine": "x86_64",
        "minimum_python_full_version": "3.12.3",
        "python_version": "3.12",
        "system": "linux",
    }


def test_dependency_markers_use_exact_python_and_reject_kernel_facts(
    tmp_path: Path,
) -> None:
    patch_marker = tmp_path / "patch-marker"
    patch_marker.mkdir()
    _build_wheel(
        patch_marker,
        requirements=('missing; python_full_version >= "3.12.4"',),
    )
    with pytest.raises(DeploymentError) as caught:
        build_platform_wheel_lock(patch_marker)
    assert caught.value.issue.code == "PLATFORM_RUNTIME_CLOSURE_INVALID"

    incompatible_python = tmp_path / "incompatible-python"
    incompatible_python.mkdir()
    _build_wheel(incompatible_python, requires_python=">=3.12.4,<3.13")
    with pytest.raises(DeploymentError) as caught:
        build_platform_wheel_lock(incompatible_python)
    assert caught.value.issue.code == "PLATFORM_RUNTIME_CLOSURE_INVALID"

    kernel_marker = tmp_path / "kernel-marker"
    kernel_marker.mkdir()
    _build_wheel(
        kernel_marker,
        requirements=('missing; platform_release == "unsupported"',),
    )
    with pytest.raises(DeploymentError) as caught:
        build_platform_wheel_lock(kernel_marker)
    assert caught.value.issue.code == "PLATFORM_RUNTIME_CLOSURE_INVALID"


def test_collected_closure_imports_packaging_without_host_site_packages(
    tmp_path: Path,
) -> None:
    import packaging

    wheelhouse = tmp_path / "wheelhouse"
    wheelhouse.mkdir()
    packaging_root = Path(packaging.__file__).parent
    packaging_files = {
        f"packaging/{path.relative_to(packaging_root).as_posix()}": path.read_bytes()
        for path in packaging_root.rglob("*.py")
    }
    platform_wheel = _build_wheel(
        wheelhouse,
        requirements=(f"packaging=={packaging.__version__}",),
        files={
            "runtime_probe.py": (
                b"from packaging.requirements import Requirement\n"
                b"VALUE = Requirement('packaging>=24.2').name\n"
            )
        },
    )
    packaging_wheel = _build_wheel(
        wheelhouse,
        name="packaging",
        version=packaging.__version__,
        files=packaging_files,
    )
    lock_path = tmp_path / "wheel-lock.json"
    _write_lock(lock_path, (platform_wheel, packaging_wheel))
    destination = tmp_path / "closure"
    collect_platform_runtime_closure(wheelhouse, lock_path, destination)
    site_packages = destination.joinpath(*Path(PLATFORM_RUNTIME_SITE_PACKAGES).parts)

    completed = subprocess.run(
        [
            "/usr/bin/python3.12",
            "-I",
            "-S",
            "-c",
            (
                "import pathlib,sys; "
                "site=pathlib.Path(sys.argv[1]).resolve(); "
                "sys.path.insert(0,str(site)); "
                "import packaging,runtime_probe; "
                "assert pathlib.Path(packaging.__file__).resolve().is_relative_to(site); "
                "assert runtime_probe.VALUE == 'packaging'; "
                "print('ok')"
            ),
            str(site_packages),
        ],
        cwd=tmp_path,
        env={"PATH": "/usr/bin:/bin", "PYTHONDONTWRITEBYTECODE": "1"},
        capture_output=True,
        text=True,
        timeout=10,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout == "ok\n"
    assert completed.stderr == ""


def test_locked_closure_constructs_openai_agent_client_offline(
    tmp_path: Path,
) -> None:
    """Prove isolated provider construction with a closure-owned SDK stub.

    The lock, extra selection, transitive dependency graph, extraction, actual
    HelixWeave factory, and actual provider wrapper are production code.  The
    tiny ``openai`` module is deliberately only a constructor-compatible stub:
    this test proves offline import provenance and client wiring, not upstream
    SDK request behavior.
    """
    wheelhouse = tmp_path / "wheelhouse"
    wheelhouse.mkdir()
    package_root = Path(runtime_module.__file__).parents[1]
    provider_files = {
        "encode_pipeline/__init__.py": (package_root / "__init__.py").read_bytes(),
        "encode_pipeline/services/__init__.py": b"",
        "encode_pipeline/services/llm_client.py": (
            package_root / "services" / "llm_client.py"
        ).read_bytes(),
        "encode_pipeline/services/llm_factory.py": (
            package_root / "services" / "llm_factory.py"
        ).read_bytes(),
        "encode_pipeline/services/llm_providers/__init__.py": (
            package_root / "services" / "llm_providers" / "__init__.py"
        ).read_bytes(),
        "encode_pipeline/services/llm_providers/openai_client.py": (
            package_root / "services" / "llm_providers" / "openai_client.py"
        ).read_bytes(),
    }
    wheels = (
        _build_wheel(
            wheelhouse,
            requirements=(
                'api-runtime==1.0.0; extra == "api"',
                'openai==1.0.0; extra == "llm-openai"',
            ),
            files=provider_files,
        ),
        _build_wheel(
            wheelhouse,
            name="api-runtime",
            requirements=("httpx==1.0.0",),
        ),
        _build_wheel(
            wheelhouse,
            name="openai",
            requirements=("httpx==1.0.0", "pydantic==1.0.0"),
            files={
                "openai/__init__.py": (
                    b"from pathlib import Path\n"
                    b"class AsyncOpenAI:\n"
                    b"    def __init__(self, *, api_key):\n"
                    b"        self.api_key = api_key\n"
                    b"        self.origin = Path(__file__).resolve()\n"
                )
            },
        ),
        _build_wheel(
            wheelhouse,
            name="httpx",
            requirements=("httpcore==1.0.0", "certifi==1.0.0"),
        ),
        _build_wheel(
            wheelhouse,
            name="httpcore",
            requirements=("h11==1.0.0",),
        ),
        _build_wheel(
            wheelhouse,
            name="pydantic",
            requirements=("typing-extensions==1.0.0",),
        ),
        _build_wheel(wheelhouse, name="certifi"),
        _build_wheel(wheelhouse, name="h11"),
        _build_wheel(wheelhouse, name="typing-extensions"),
    )
    lock = build_platform_wheel_lock(wheelhouse)
    assert {item.filename for item in lock.wheels} == {item.name for item in wheels}
    lock_path = tmp_path / "wheel-lock.json"
    lock_path.write_bytes(lock.to_bytes())
    lock_path.chmod(0o644)
    destination = tmp_path / "closure"
    collect_platform_runtime_closure(wheelhouse, lock_path, destination)
    site_packages = destination.joinpath(*Path(PLATFORM_RUNTIME_SITE_PACKAGES).parts)
    outside = tmp_path / "outside"
    outside.mkdir()

    completed = subprocess.run(
        [
            "/usr/bin/python3.12",
            "-I",
            "-S",
            "-c",
            (
                "import pathlib,socket,sys; "
                "site=pathlib.Path(sys.argv[1]).resolve(); "
                "assert sys.flags.isolated and sys.flags.no_site; "
                "assert all('site-packages' not in item for item in sys.path); "
                "sys.path.insert(0,str(site)); "
                "network_attempts=[]; "
                "blocked=lambda *args,**kwargs: "
                "(network_attempts.append((args,kwargs)), "
                "(_ for _ in ()).throw(AssertionError('network disabled')))[1]; "
                "socket.socket=blocked; socket.create_connection=blocked; "
                "from encode_pipeline.services.llm_factory import create_llm_client; "
                "import openai; "
                "client=create_llm_client(); "
                "assert type(client).__name__ == 'OpenAILLMClient'; "
                "assert isinstance(client._client, openai.AsyncOpenAI); "
                "assert pathlib.Path(openai.__file__).resolve().is_relative_to(site); "
                "assert client._client.origin.is_relative_to(site); "
                "assert client._client.api_key == 'offline-test-key'; "
                "assert network_attempts == []; "
                "print('ok')"
            ),
            str(site_packages),
        ],
        cwd=outside,
        env={
            "ENCODE_AGENT_LLM_API_KEY": "offline-test-key",
            "ENCODE_AGENT_LLM_PROVIDER": "openai",
            "PATH": "/usr/bin:/bin",
            "PYTHONDONTWRITEBYTECODE": "1",
        },
        capture_output=True,
        text=True,
        timeout=10,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout == "ok\n"
    assert completed.stderr == ""


def test_collector_builds_canonical_wheel_only_application_tree(tmp_path: Path) -> None:
    closure, destination, lock = _valid_closure(tmp_path)

    assert closure.lock_identity == lock.identity
    assert closure == inspect_platform_runtime_closure(destination, lock.identity)
    logical_paths = {item.logical_path for item in closure.files}
    assert PLATFORM_RUNTIME_LOCK_PATH in logical_paths
    assert PLATFORM_RUNTIME_LAUNCHER_PATH in logical_paths
    assert f"{PLATFORM_RUNTIME_SITE_PACKAGES}/helixweave/__init__.py" in logical_paths
    assert any(
        "/wheelhouse/helixweave-1.0.0-py3-none-any.whl" in item
        for item in logical_paths
    )
    assert all(item.source.is_relative_to(destination) for item in closure.files)
    assert all(
        stat.S_IMODE(item.source.stat().st_mode) == item.mode for item in closure.files
    )
    assert all(
        item.mode
        == (0o555 if item.logical_path == PLATFORM_RUNTIME_LAUNCHER_PATH else 0o444)
        for item in closure.files
    )
    assert not any(item.logical_path.endswith(".pyc") for item in closure.files)


@pytest.mark.parametrize("fault", ("hash", "missing", "extra"))
def test_collector_rejects_non_exact_wheelhouse(tmp_path: Path, fault: str) -> None:
    wheelhouse = tmp_path / "wheelhouse"
    wheelhouse.mkdir()
    wheel = _build_wheel(wheelhouse)
    lock_path = tmp_path / "wheel-lock.json"
    lock = _write_lock(lock_path, (wheel,))
    if fault == "hash":
        document = lock.to_dict()
        document["wheels"][0]["sha256"] = "0" * 64
        lock_path.write_bytes(canonical_json_bytes(document))
    elif fault == "missing":
        wheel.unlink()
    else:
        _build_wheel(wheelhouse, name="unexpected")

    with pytest.raises(DeploymentError) as caught:
        collect_platform_runtime_closure(wheelhouse, lock_path, tmp_path / "closure")

    assert caught.value.issue.code in {
        "PLATFORM_RUNTIME_CLOSURE_INVALID",
        "PLATFORM_RUNTIME_LOCK_INVALID",
    }
    assert not (tmp_path / "closure").exists()


@pytest.mark.parametrize(
    ("limit_name", "limit_value"),
    (("MAX_BUNDLE_FILES", 1), ("MAX_BUNDLE_BYTES", 1)),
)
def test_collector_preflights_consumer_aggregate_bounds_before_writing(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    limit_name: str,
    limit_value: int,
) -> None:
    wheelhouse = tmp_path / "wheelhouse"
    wheelhouse.mkdir()
    wheel = _build_wheel(wheelhouse)
    lock_path = tmp_path / "wheel-lock.json"
    _write_lock(lock_path, (wheel,))
    monkeypatch.setattr(runtime_module, limit_name, limit_value)
    destination = tmp_path / "closure"

    with pytest.raises(DeploymentError) as caught:
        collect_platform_runtime_closure(wheelhouse, lock_path, destination)

    assert caught.value.issue.code == "PLATFORM_RUNTIME_CLOSURE_INVALID"
    assert not destination.exists()
    assert not tuple(tmp_path.glob(".closure.*.partial"))


def test_collector_preflights_capacity_before_writing(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    wheelhouse = tmp_path / "wheelhouse"
    wheelhouse.mkdir()
    wheel = _build_wheel(wheelhouse)
    lock_path = tmp_path / "wheel-lock.json"
    _write_lock(lock_path, (wheel,))
    monkeypatch.setattr(
        runtime_module.shutil,
        "disk_usage",
        lambda _path: SimpleNamespace(free=0),
    )
    destination = tmp_path / "closure"

    with pytest.raises(DeploymentError) as caught:
        collect_platform_runtime_closure(wheelhouse, lock_path, destination)

    assert caught.value.issue.code == "DEPLOYMENT_CAPACITY_INSUFFICIENT"
    assert caught.value.issue.recoverable is True
    assert not destination.exists()
    assert not tuple(tmp_path.glob(".closure.*.partial"))


@pytest.mark.parametrize(
    ("tag", "extra_files", "requirements", "dependency"),
    (
        ("cp311-cp311-manylinux_2_17_x86_64", {}, (), False),
        ("cp312-cp312-macosx_13_0_x86_64", {}, (), False),
        ("py3-none-any", {"unsafe.pth": b"import os\n"}, (), False),
        ("py3-none-any", {"../escape.py": b"VALUE = 1\n"}, (), False),
        ("py3-none-any", {"shared.py": b"platform\n"}, ("helper==1.0.0",), True),
    ),
)
def test_collector_rejects_abi_pth_escape_and_collision(
    tmp_path: Path,
    tag: str,
    extra_files: dict[str, bytes],
    requirements: tuple[str, ...],
    dependency: bool,
) -> None:
    wheelhouse = tmp_path / "wheelhouse"
    wheelhouse.mkdir()
    wheels = [
        _build_wheel(
            wheelhouse,
            tag=tag,
            requirements=requirements,
            files=extra_files,
        )
    ]
    if dependency:
        wheels.append(
            _build_wheel(
                wheelhouse,
                name="helper",
                files={"shared.py": b"dependency\n"},
            )
        )
    lock_path = tmp_path / "wheel-lock.json"
    _write_lock(lock_path, tuple(wheels))

    with pytest.raises(DeploymentError) as caught:
        collect_platform_runtime_closure(wheelhouse, lock_path, tmp_path / "closure")

    assert caught.value.issue.code == "PLATFORM_RUNTIME_CLOSURE_INVALID"
    assert not (tmp_path / "closure").exists()


def test_collector_accepts_cp312_and_older_abi3_manylinux_tags(tmp_path: Path) -> None:
    for index, tag in enumerate(
        (
            "cp312-cp312-manylinux_2_17_x86_64",
            "cp39-abi3-manylinux2014_x86_64",
        )
    ):
        case = tmp_path / str(index)
        case.mkdir()
        wheelhouse = case / "wheelhouse"
        wheelhouse.mkdir()
        wheel = _build_wheel(wheelhouse, tag=tag)
        lock_path = case / "wheel-lock.json"
        _write_lock(lock_path, (wheel,))

        closure = collect_platform_runtime_closure(
            wheelhouse, lock_path, case / "closure"
        )

        assert closure.files


@pytest.mark.parametrize("missing", ("openai", "httpcore"))
def test_dependency_graph_requires_complete_selected_provider_closure(
    tmp_path: Path,
    missing: str,
) -> None:
    wheelhouse = tmp_path / "wheelhouse"
    wheelhouse.mkdir()
    wheels = [
        _build_wheel(
            wheelhouse,
            requirements=(
                'api-runtime==1.0.0; extra == "api"',
                'openai==1.0.0; extra == "llm-openai"',
            ),
        ),
        _build_wheel(wheelhouse, name="api-runtime"),
    ]
    if missing != "openai":
        wheels.extend(
            (
                _build_wheel(
                    wheelhouse,
                    name="openai",
                    requirements=("httpx==1.0.0",),
                ),
                _build_wheel(
                    wheelhouse,
                    name="httpx",
                    requirements=("httpcore==1.0.0",),
                ),
            )
        )
    lock_path = tmp_path / "wheel-lock.json"
    _write_lock(lock_path, tuple(wheels))

    with pytest.raises(DeploymentError) as caught:
        collect_platform_runtime_closure(
            wheelhouse,
            lock_path,
            tmp_path / "closure",
        )

    assert caught.value.issue.code == "PLATFORM_RUNTIME_CLOSURE_INVALID"
    assert not (tmp_path / "closure").exists()


def test_dependency_graph_accepts_complete_api_and_openai_root_extras(
    tmp_path: Path,
) -> None:
    wheelhouse = tmp_path / "wheelhouse"
    wheelhouse.mkdir()
    wheels = (
        _build_wheel(
            wheelhouse,
            requirements=(
                'api-runtime==1.0.0; extra == "api"',
                'openai==1.0.0; extra == "llm-openai"',
            ),
        ),
        _build_wheel(wheelhouse, name="api-runtime"),
        _build_wheel(
            wheelhouse,
            name="openai",
            requirements=("httpx==1.0.0",),
        ),
        _build_wheel(
            wheelhouse,
            name="httpx",
            requirements=("httpcore==1.0.0",),
        ),
        _build_wheel(wheelhouse, name="httpcore"),
    )
    lock_path = tmp_path / "wheel-lock.json"
    _write_lock(lock_path, wheels)

    closure = collect_platform_runtime_closure(
        wheelhouse,
        lock_path,
        tmp_path / "closure",
    )

    assert closure.files


def test_dependency_graph_rejects_unreachable_wheels(tmp_path: Path) -> None:
    wheelhouse = tmp_path / "wheelhouse"
    wheelhouse.mkdir()
    platform_wheel = _build_wheel(
        wheelhouse,
        requirements=('needed==1.0.0; extra == "api"',),
    )
    needed = _build_wheel(wheelhouse, name="needed")
    unused = _build_wheel(wheelhouse, name="unused")
    lock_path = tmp_path / "wheel-lock.json"
    _write_lock(lock_path, (platform_wheel, needed, unused))

    with pytest.raises(DeploymentError) as caught:
        collect_platform_runtime_closure(wheelhouse, lock_path, tmp_path / "closure")

    assert caught.value.issue.code == "PLATFORM_RUNTIME_CLOSURE_INVALID"


def test_supported_python_observation_is_path_and_abi_bound() -> None:
    identity = verify_supported_python_runtime(
        executable=Path("/usr/bin/python3.12"),
        version_info=(3, 12, 3),
        implementation_name="cpython",
        implementation_version="3.12.3",
        cache_tag="cpython-312",
        machine="x86_64",
        system="linux",
        libc_version=(2, 17),
    )
    assert identity.startswith("sha256-")

    patched_identity = verify_supported_python_runtime(
        executable=Path("/usr/bin/python3.12"),
        version_info=(3, 12, 13),
        implementation_name="cpython",
        implementation_version="3.12.13",
        cache_tag="cpython-312",
        machine="x86_64",
        system="linux",
        libc_version=(2, 17),
    )
    assert patched_identity.startswith("sha256-")
    assert patched_identity != identity

    with pytest.raises(DeploymentError) as caught:
        verify_supported_python_runtime(
            executable=Path("/usr/bin/python3.12"),
            version_info=(3, 12, 2),
            implementation_name="cpython",
            implementation_version="3.12.2",
            cache_tag="cpython-312",
            machine="x86_64",
            system="linux",
            libc_version=(2, 17),
        )
    assert caught.value.issue.code == "PLATFORM_PYTHON_INCOMPATIBLE"


def _database_request(
    *,
    database_mode: str = "fresh-candidate",
    target_schema_heads: tuple[str, ...] | None = None,
) -> DatabasePrepareRequest:
    from encode_pipeline.persistence.migration_admission import (
        verify_migration_execution_inventory,
    )

    return DatabasePrepareRequest.create(
        operation="activate",
        database_mode=database_mode,
        task_identity=TASK,
        deployment_identity=IDENTITY_A,
        prior_state_identity=IDENTITY_B,
        candidate_state_identity=IDENTITY_C,
        action_receipt_identity=IDENTITY_A,
        backup_receipt_identity=(
            None if database_mode == "fresh-candidate" else IDENTITY_B
        ),
        target_schema_heads=(
            verify_migration_execution_inventory().heads
            if target_schema_heads is None
            else target_schema_heads
        ),
    )


def _encode_prepare_request() -> EncodeRuntimePrepareRequest:
    return EncodeRuntimePrepareRequest.create(
        task_identity=TASK,
        deployment_identity=IDENTITY_B,
        authority_platform_identity=IDENTITY_A,
        prior_state_identity=IDENTITY_C,
        candidate_state_identity=IDENTITY_B,
    )


def _bulk_prepare_request() -> BulkRuntimePrepareRequest:
    return BulkRuntimePrepareRequest.create(
        operation="activate",
        task_identity=TASK,
        candidate_bulk_identity=IDENTITY_B,
        authority_platform_identity=IDENTITY_A,
        prior_state_identity=IDENTITY_C,
        candidate_state_identity=IDENTITY_B,
        docker_service_identity=IDENTITY_A,
        docker_client_identity=IDENTITY_B,
        docker_endpoint_identity=IDENTITY_C,
        docker_daemon_uid=os.getuid(),
        docker_daemon_gid=os.getgid(),
    )


def test_database_prepare_runs_candidate_inventory_and_returns_canonical_receipt(
    tmp_path: Path,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path)
    layout.database.parent.mkdir(parents=True)
    request = _database_request()

    receipt = prepare_candidate_database(
        request,
        layout=layout,
        expected_database_uid=os.getuid(),
        expected_database_gid=os.getgid(),
    )

    candidate = fresh_database_candidate_path(layout, TASK)
    assert candidate.is_file()
    assert not layout.database.exists()
    assert receipt.request_identity == request.identity
    assert receipt.database_before_identity is None
    assert receipt.database_after_identity.startswith("sha256-")
    assert len(receipt.schema_heads) == 1
    content = canonical_json_bytes(receipt.to_dict())
    assert DatabasePrepareReceipt.from_dict(json.loads(content)) == receipt
    assert str(candidate) not in content.decode()


def test_database_prepare_rejects_unknown_head_without_creating_candidate(
    tmp_path: Path,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path)
    layout.database.parent.mkdir(parents=True)
    request = _database_request(target_schema_heads=("unknown-head",))

    with pytest.raises(DeploymentError) as caught:
        prepare_candidate_database(
            request,
            layout=layout,
            expected_database_uid=os.getuid(),
            expected_database_gid=os.getgid(),
        )

    assert caught.value.issue.code == "DB_PREPARE_FAILED"
    assert not fresh_database_candidate_path(layout, TASK).exists()
    assert not layout.database.exists()


def test_database_prepare_mode_is_explicit_and_path_free(tmp_path: Path) -> None:
    layout = DeploymentLayout.isolated(tmp_path)
    layout.database.parent.mkdir(parents=True)
    fresh = _database_request()
    candidate = fresh_database_candidate_path(layout, TASK)

    prepare_candidate_database(
        fresh,
        layout=layout,
        expected_database_uid=os.getuid(),
        expected_database_gid=os.getgid(),
    )
    candidate.rename(layout.database)
    existing = _database_request(database_mode="existing-live")
    receipt = prepare_candidate_database(
        existing,
        layout=layout,
        expected_database_uid=os.getuid(),
        expected_database_gid=os.getgid(),
    )

    assert receipt.database_before_identity is not None
    assert existing.database_mode == "existing-live"
    assert not any("path" in key or "url" in key for key in existing.to_dict())
    injected = existing.to_dict()
    injected["database_path"] = str(tmp_path / "attacker-selected.db")
    with pytest.raises(DeploymentError):
        DatabasePrepareRequest.from_dict(injected)

    with pytest.raises(DeploymentError) as caught:
        prepare_candidate_database(
            fresh,
            layout=layout,
            expected_database_uid=os.getuid(),
            expected_database_gid=os.getgid(),
        )
    assert caught.value.issue.code == "DB_PREPARE_FAILED"
    assert not candidate.exists()


def test_fresh_database_prepare_retry_returns_existing_candidate_receipt(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path)
    layout.database.parent.mkdir(parents=True)
    request = _database_request()
    candidate = fresh_database_candidate_path(layout, TASK)
    first = prepare_candidate_database(
        request,
        layout=layout,
        expected_database_uid=os.getuid(),
        expected_database_gid=os.getgid(),
    )
    before = candidate.stat()
    monkeypatch.setattr(
        runtime_module,
        "_supported_database_ownership",
        lambda: (os.getuid(), os.getgid()),
    )

    resumed = prepare_candidate_database(request, layout=layout)

    after = candidate.stat()
    assert resumed == first
    assert (after.st_dev, after.st_ino, after.st_mtime_ns, after.st_size) == (
        before.st_dev,
        before.st_ino,
        before.st_mtime_ns,
        before.st_size,
    )
    assert not layout.database.exists()


@pytest.mark.parametrize("fault", ("request", "layout", "missing-existing"))
def test_database_prepare_rejects_invalid_public_coordinates(
    tmp_path: Path,
    fault: str,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path)
    layout.database.parent.mkdir(parents=True)
    request: object = _database_request()
    selected_layout: object = layout
    expected = "DB_PREPARE_FAILED"
    if fault == "request":
        request = object()
        expected = "DB_PREPARE_REQUEST_INVALID"
    elif fault == "layout":
        selected_layout = object()
    else:
        request = _database_request(database_mode="existing-live")

    with pytest.raises(DeploymentError) as caught:
        prepare_candidate_database(
            request,
            layout=selected_layout,
            expected_database_uid=os.getuid(),
            expected_database_gid=os.getgid(),
        )

    assert caught.value.issue.code == expected
    assert not layout.database.exists()
    assert str(tmp_path) not in str(caught.value)


def test_database_prepare_preserves_schema_incompatibility_error(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import encode_pipeline.deployment.database as database_module
    import encode_pipeline.persistence.migrations as migrations_module

    layout = DeploymentLayout.isolated(tmp_path)
    layout.database.parent.mkdir(parents=True)
    request = _database_request()
    monkeypatch.setattr(migrations_module, "upgrade_database", lambda _url: None)
    monkeypatch.setattr(
        database_module,
        "inspect_database",
        lambda *_args, **_kwargs: SimpleNamespace(schema_heads=("unexpected",)),
    )

    with pytest.raises(DeploymentError) as caught:
        prepare_candidate_database(
            request,
            layout=layout,
            expected_database_uid=os.getuid(),
            expected_database_gid=os.getgid(),
        )

    assert caught.value.issue.code == "DATABASE_SCHEMA_INCOMPATIBLE"
    assert caught.value.issue.recoverable is False
    assert not fresh_database_candidate_path(layout, TASK).exists()


def test_contract_compatibility_distinguishes_missing_and_present_mismatch() -> None:
    contract = "helixweave.test.runtime"
    platform = ResolvedContractFacts(
        component="platform",
        deployment_identity=IDENTITY_A,
        version="1.0.0",
        contracts=(),
        requirements=(ContractRequirement(contract, (IDENTITY_B,)),),
        database_heads=("head",),
    )
    matching = ResolvedContractFacts(
        component="encode-runtime",
        deployment_identity=IDENTITY_B,
        version="1.0.0",
        contracts=(ContractIdentity(contract, IDENTITY_B),),
    )
    mismatching = ResolvedContractFacts(
        component="encode-runtime",
        deployment_identity=IDENTITY_C,
        version="2.0.0",
        contracts=(ContractIdentity(contract, IDENTITY_C),),
    )
    conflicting_duplicate = ResolvedContractFacts(
        component="bulk-rnaseq-runtime",
        deployment_identity=IDENTITY_C,
        version="1.0.0",
        contracts=(ContractIdentity(contract, IDENTITY_C),),
    )

    assert resolved_facts_compatibility((platform, matching)) == "compatible"
    assert resolved_facts_compatibility((platform,)) == "incomplete"
    assert resolved_facts_compatibility((platform, mismatching)) == "incompatible"
    assert (
        resolved_facts_compatibility((platform, matching, conflicting_duplicate))
        == "incompatible"
    )


@pytest.mark.parametrize(
    ("bulk_active", "bulk_resolves", "expected"),
    (
        (False, False, "incomplete"),
        (True, False, "incompatible"),
        (True, True, "compatible"),
    ),
)
def test_action_probe_never_treats_present_native_failure_as_incomplete(
    monkeypatch: pytest.MonkeyPatch,
    bulk_active: bool,
    bulk_resolves: bool,
    expected: str,
) -> None:
    from encode_pipeline.deployment.bundle import BundleStore
    import encode_pipeline.deployment.database as database_module
    from encode_pipeline.deployment.models import BULK_RNASEQ_RUNTIME, PLATFORM
    from encode_pipeline.deployment.native_contracts import (
        PLATFORM_FRONTEND_CONTRACT,
        PLATFORM_PYTHON_RUNTIME_CONTRACT,
        PLATFORM_REFERENCES_CONTRACT,
        ProductionNativeContractResolver,
    )
    import encode_pipeline.frontend_assets as frontend_assets_module
    from encode_pipeline.persistence.migration_admission import (
        verify_migration_execution_inventory,
    )

    provider_contract = "helixweave.test.bulk-runtime"
    platform_contracts = tuple(
        sorted(
            (
                ContractIdentity(PLATFORM_FRONTEND_CONTRACT, IDENTITY_A),
                ContractIdentity(PLATFORM_PYTHON_RUNTIME_CONTRACT, IDENTITY_B),
                ContractIdentity(PLATFORM_REFERENCES_CONTRACT, IDENTITY_C),
            )
        )
    )
    inventory = verify_migration_execution_inventory()
    platform_facts = ResolvedContractFacts(
        component=PLATFORM,
        deployment_identity=IDENTITY_A,
        version="1.0.0",
        contracts=platform_contracts,
        requirements=(ContractRequirement(provider_contract, (IDENTITY_B,)),),
        database_heads=inventory.heads,
    )
    bulk_contracts = (ContractIdentity(provider_contract, IDENTITY_B),)
    bulk_facts = ResolvedContractFacts(
        component=BULK_RNASEQ_RUNTIME,
        deployment_identity=IDENTITY_B,
        version="1.0.0",
        contracts=bulk_contracts,
    )
    manifests = {
        PLATFORM: SimpleNamespace(
            component=PLATFORM,
            identity=IDENTITY_A,
            provides=platform_contracts,
        ),
        BULK_RNASEQ_RUNTIME: SimpleNamespace(
            component=BULK_RNASEQ_RUNTIME,
            identity=IDENTITY_B,
            provides=bulk_contracts,
        ),
    }

    monkeypatch.setattr(
        BundleStore,
        "verify_installed",
        lambda _self, component, _identity, **_kwargs: manifests[component],
    )

    def resolve(_self, _root, manifest):
        if manifest.component == BULK_RNASEQ_RUNTIME:
            if not bulk_resolves:
                raise ValueError("candidate native closure unavailable")
            return bulk_facts
        return platform_facts

    monkeypatch.setattr(ProductionNativeContractResolver, "resolve", resolve)
    monkeypatch.setattr(
        runtime_module, "verify_supported_python_runtime", lambda: IDENTITY_B
    )
    monkeypatch.setattr(
        frontend_assets_module,
        "load_packaged_frontend_assets",
        lambda: SimpleNamespace(
            manifest=SimpleNamespace(
                identity=IDENTITY_A,
                api_contract_sha256="d" * 64,
            )
        ),
    )
    monkeypatch.setattr(
        runtime_module,
        "_supported_database_ownership",
        lambda: (os.getuid(), os.getgid()),
    )
    monkeypatch.setattr(
        database_module,
        "inspect_database",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(PermissionError()),
    )
    monkeypatch.setattr(
        runtime_module,
        "_permission_readiness",
        lambda *_args: ReadinessCheck("ready", "READY", IDENTITY_A),
    )
    monkeypatch.setattr(
        runtime_module,
        "_reference_readiness",
        lambda *_args: ReadinessCheck("unavailable", "REFERENCE_NOT_READY"),
    )
    monkeypatch.setattr(
        runtime_module,
        "_redis_ping",
        lambda *_args: ReadinessCheck("unavailable", "REDIS_UNAVAILABLE"),
    )
    monkeypatch.setattr(
        runtime_module,
        "_docker_ping",
        lambda *_args: ReadinessCheck("unavailable", "DOCKER_UNAVAILABLE"),
    )
    request = DeploymentActionRequest.create(
        phase="admit",
        operation="activate",
        component=PLATFORM,
        task_identity=TASK,
        deployment_identity=IDENTITY_A,
        authority_platform_identity=IDENTITY_A,
        prior_state_identity=IDENTITY_B,
        candidate_state_identity=IDENTITY_C,
        candidate_active={
            "platform": IDENTITY_A,
            "encode-runtime": None,
            "bulk-rnaseq-runtime": IDENTITY_B if bulk_active else None,
        },
    )
    layout = DeploymentLayout.supported()

    receipt = runtime_module._default_platform_action_probe(
        request, layout.component_store(PLATFORM) / IDENTITY_A
    )

    assert receipt.compatibility == expected
    assert receipt.readiness["configuration"] == ReadinessCheck(
        "not-applicable", "NOT_APPLICABLE"
    )
    assert receipt.readiness["database-schema"].status == "unavailable"


def test_verified_migration_inventory_binds_sorted_known_ancestors() -> None:
    from encode_pipeline.persistence.migration_admission import (
        verify_migration_execution_inventory,
    )

    inventory = verify_migration_execution_inventory()
    identity, target, known = runtime_module._migration_inventory_evidence(inventory)

    assert identity == f"sha256-{inventory.contract_sha256}"
    assert target == inventory.heads
    assert known == tuple(sorted(known))
    assert set(known) == {revision.revision for revision in inventory.revisions} - set(
        inventory.heads
    )


def _action_request() -> DeploymentActionRequest:
    return DeploymentActionRequest.create(
        phase="admit",
        operation="activate",
        component="platform",
        task_identity=TASK,
        deployment_identity=IDENTITY_A,
        authority_platform_identity=IDENTITY_A,
        prior_state_identity=IDENTITY_B,
        candidate_state_identity=IDENTITY_C,
        candidate_active={
            "platform": IDENTITY_A,
            "encode-runtime": IDENTITY_B,
            "bulk-rnaseq-runtime": IDENTITY_C,
        },
    )


def _ready_action_receipt(request: DeploymentActionRequest) -> DeploymentActionReceipt:
    return DeploymentActionReceipt.create(
        request_identity=request.identity,
        status="admitted",
        compatibility="compatible",
        database_before_identity=IDENTITY_A,
        accepted_schema_heads=("head",),
        target_schema_heads=("head",),
        migration_inventory_identity=IDENTITY_C,
        known_schema_revisions=("ancestor",),
        migration_required=False,
        rollback_supported=True,
        api_contract_sha256="d" * 64,
        native_identities={
            "platform": IDENTITY_A,
            "encode-runtime": IDENTITY_B,
            "bulk-rnaseq-runtime": IDENTITY_C,
        },
        frontend_identity=IDENTITY_A,
        reference_compatibility_identity=IDENTITY_B,
        readiness={
            check: ReadinessCheck("ready", "READY", IDENTITY_A)
            for check in VERIFICATION_CHECKS
        },
    )


def test_action_request_and_receipt_are_canonical_path_free_contracts() -> None:
    request = _action_request()
    receipt = _ready_action_receipt(request)
    request_content = canonical_json_bytes(request.to_dict())
    receipt_content = canonical_json_bytes(receipt.to_dict())

    assert DeploymentActionRequest.from_dict(json.loads(request_content)) == request
    assert DeploymentActionReceipt.from_dict(json.loads(receipt_content)) == receipt
    assert receipt.status == "admitted"
    assert b"/opt/" not in receipt_content
    assert b"/var/" not in receipt_content

    noncanonical = json.dumps(request.to_dict(), sort_keys=True, indent=2).encode()
    with pytest.raises(DeploymentError) as caught:
        runtime_module._parse_root_request(
            noncanonical, DeploymentActionRequest.from_dict
        )
    assert caught.value.issue.code == "PLATFORM_SERVICE_REQUEST_UNTRUSTED"

    with_path = request.to_dict()
    with_path["path"] = "/tmp/caller-selected"
    with pytest.raises(DeploymentError):
        DeploymentActionRequest.from_dict(with_path)


def test_action_receipt_rejects_unknown_or_noncanonical_schema_evidence() -> None:
    request = _action_request()
    common = {
        "request_identity": request.identity,
        "status": "admitted",
        "compatibility": "compatible",
        "database_before_identity": IDENTITY_A,
        "target_schema_heads": ("head",),
        "migration_inventory_identity": IDENTITY_C,
        "api_contract_sha256": "d" * 64,
        "native_identities": {
            "platform": IDENTITY_A,
            "encode-runtime": IDENTITY_B,
            "bulk-rnaseq-runtime": IDENTITY_C,
        },
        "frontend_identity": IDENTITY_A,
        "reference_compatibility_identity": IDENTITY_B,
        "readiness": {
            check: ReadinessCheck("ready", "READY", IDENTITY_A)
            for check in VERIFICATION_CHECKS
        },
    }

    with pytest.raises(DeploymentError):
        DeploymentActionReceipt.create(
            **common,
            accepted_schema_heads=("unknown",),
            known_schema_revisions=("ancestor",),
            migration_required=True,
            rollback_supported=False,
        )
    with pytest.raises(DeploymentError):
        DeploymentActionReceipt.create(
            **common,
            accepted_schema_heads=("head",),
            known_schema_revisions=("z-revision", "a-revision"),
            migration_required=False,
            rollback_supported=True,
        )


def test_candidate_service_dispatch_is_a_closed_enumeration(
    capsys: pytest.CaptureFixture[str],
) -> None:
    called: list[str] = []

    assert (
        candidate_service_main(
            ("api",),
            release_root=Path(f"/tmp/{IDENTITY_A}"),
            api_runner=lambda: called.append("api") or 0,
        )
        == 0
    )
    assert (
        candidate_service_main(
            ("worker",),
            release_root=Path(f"/tmp/{IDENTITY_A}"),
            worker_runner=lambda: called.append("worker") or 0,
        )
        == 0
    )
    assert (
        candidate_service_main(
            ("encode-runtime-prepare",),
            release_root=Path(f"/tmp/{IDENTITY_A}"),
            encode_runtime_preparer=lambda: (
                called.append("encode-runtime-prepare") or 0
            ),
        )
        == 0
    )
    assert (
        candidate_service_main(
            ("bulk-runtime-prepare",),
            release_root=Path(f"/tmp/{IDENTITY_A}"),
            bulk_runtime_preparer=lambda: called.append("bulk-runtime-prepare") or 0,
        )
        == 0
    )
    assert called == [
        "api",
        "worker",
        "encode-runtime-prepare",
        "bulk-runtime-prepare",
    ]

    assert (
        candidate_service_main(
            ("encode_pipeline.workers.cli",),
            release_root=Path(f"/tmp/{IDENTITY_A}"),
        )
        == 64
    )
    captured = capsys.readouterr()
    assert "PLATFORM_SERVICE_ACTION_INVALID" in captured.err
    assert "Traceback" not in captured.err


@pytest.mark.parametrize(
    ("mode", "expected_code"),
    (
        ("encode-runtime-prepare", "ENCODE_RUNTIME_PREPARE_REQUEST_MISMATCH"),
        ("bulk-runtime-prepare", "BULK_RUNTIME_PREPARE_REQUEST_MISMATCH"),
        ("db-prepare", "DB_PREPARE_REQUEST_MISMATCH"),
        ("operator-action", "OPERATOR_ACTION_REQUEST_MISMATCH"),
    ),
)
def test_candidate_service_rejects_requests_for_another_release(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    mode: str,
    expected_code: str,
) -> None:
    requests = {
        "encode-runtime-prepare": _encode_prepare_request(),
        "bulk-runtime-prepare": _bulk_prepare_request(),
        "db-prepare": _database_request(),
        "operator-action": _action_request(),
    }
    request = requests[mode]
    monkeypatch.setattr(
        runtime_module,
        "_read_root_request",
        lambda _path: canonical_json_bytes(request.to_dict()),
    )

    exit_code = candidate_service_main(
        (mode,),
        release_root=Path(f"/tmp/{IDENTITY_C}"),
    )

    captured = capsys.readouterr()
    error = json.loads(captured.err)
    assert exit_code == 65
    assert captured.out == ""
    assert error["issue"]["code"] == expected_code
    assert error["issue"]["recoverable"] is False
    assert "Traceback" not in captured.err


@pytest.mark.parametrize(
    ("mode", "expected_code"),
    (
        ("encode-runtime-prepare", "ENCODE_RUNTIME_PREPARE_RECEIPT_INVALID"),
        ("bulk-runtime-prepare", "BULK_RUNTIME_PREPARE_RECEIPT_INVALID"),
        ("operator-action", "OPERATOR_ACTION_RECEIPT_INVALID"),
    ),
)
def test_candidate_service_rejects_invalid_child_receipts(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    mode: str,
    expected_code: str,
) -> None:
    action_probe = None
    if mode == "encode-runtime-prepare":
        from encode_pipeline.deployment import encode_runtime_materializer

        request = _encode_prepare_request()
        monkeypatch.setattr(
            encode_runtime_materializer.OfflineEncodeRuntimeMaterializer,
            "prepare",
            lambda *_args: object(),
        )
    elif mode == "bulk-runtime-prepare":
        from encode_pipeline.deployment import bulk_runtime_materializer

        request = _bulk_prepare_request()
        monkeypatch.setattr(
            bulk_runtime_materializer.OfflineBulkRuntimeMaterializer,
            "prepare",
            lambda *_args: object(),
        )
    else:
        request = _action_request()

        def invalid_action_probe(*_args):
            return object()

        action_probe = invalid_action_probe
    monkeypatch.setattr(
        runtime_module,
        "_read_root_request",
        lambda _path: canonical_json_bytes(request.to_dict()),
    )

    exit_code = candidate_service_main(
        (mode,),
        release_root=Path(f"/tmp/{IDENTITY_A}"),
        action_probe=action_probe,
    )

    captured = capsys.readouterr()
    error = json.loads(captured.err)
    assert exit_code == 65
    assert captured.out == ""
    assert error["issue"]["code"] == expected_code
    assert error["issue"]["recoverable"] is False


def test_candidate_service_emits_database_prepare_receipt(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    request = _database_request()
    receipt = DatabasePrepareReceipt.create(
        request_identity=request.identity,
        database_before_identity=None,
        database_after_identity=IDENTITY_B,
        schema_heads=request.target_schema_heads,
    )
    monkeypatch.setattr(
        runtime_module,
        "_read_root_request",
        lambda _path: canonical_json_bytes(request.to_dict()),
    )
    monkeypatch.setattr(
        runtime_module,
        "prepare_candidate_database",
        lambda observed: receipt if observed == request else None,
    )

    assert (
        candidate_service_main(
            ("db-prepare",),
            release_root=Path(f"/tmp/{IDENTITY_A}"),
        )
        == 0
    )
    captured = capsys.readouterr()
    assert captured.err == ""
    assert DatabasePrepareReceipt.from_dict(json.loads(captured.out)) == receipt


@pytest.mark.parametrize(
    ("failure", "exit_code", "expected_code", "recoverable"),
    (
        ("deployment", 69, "CANDIDATE_DEPENDENCY_UNAVAILABLE", True),
        ("unexpected", 70, "PLATFORM_SERVICE_FAILED", False),
    ),
)
def test_candidate_service_maps_failures_to_stable_public_errors(
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
    failure: str,
    exit_code: int,
    expected_code: str,
    recoverable: bool,
) -> None:
    def runner() -> int:
        if failure == "deployment":
            raise fail(
                "CANDIDATE_DEPENDENCY_UNAVAILABLE",
                "Candidate dependency is unavailable.",
                recoverable=True,
            )
        raise RuntimeError(str(tmp_path))

    assert (
        candidate_service_main(
            ("api",),
            release_root=Path(f"/tmp/{IDENTITY_A}"),
            api_runner=runner,
        )
        == exit_code
    )
    captured = capsys.readouterr()
    error = json.loads(captured.err)
    assert captured.out == ""
    assert error["issue"] == {
        "code": expected_code,
        "message": error["issue"]["message"],
        "recoverable": recoverable,
    }
    assert str(tmp_path) not in captured.err
    assert "Traceback" not in captured.err


@pytest.mark.parametrize("fault", ("invalid-document", "relative-request-path"))
def test_candidate_database_service_rejects_untrusted_request_documents(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
    fault: str,
) -> None:
    kwargs = {}
    if fault == "invalid-document":
        monkeypatch.setattr(
            runtime_module,
            "_read_root_request",
            lambda _path: canonical_json_bytes({"unexpected": True}),
        )
        expected = "DATABASE_PREPARE_REQUEST_INVALID"
    else:
        kwargs["database_request_path"] = Path("relative/request.json")
        expected = "PLATFORM_SERVICE_REQUEST_UNTRUSTED"

    assert (
        candidate_service_main(
            ("db-prepare",),
            release_root=Path(f"/tmp/{IDENTITY_A}"),
            **kwargs,
        )
        == 65
    )
    captured = capsys.readouterr()
    error = json.loads(captured.err)
    assert captured.out == ""
    assert error["issue"]["code"] == expected


def test_candidate_service_runs_the_fixed_offline_encode_materializer(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    from encode_pipeline.deployment import encode_runtime_materializer

    request = EncodeRuntimePrepareRequest.create(
        task_identity=TASK,
        deployment_identity=IDENTITY_B,
        authority_platform_identity=IDENTITY_A,
        prior_state_identity=IDENTITY_C,
        candidate_state_identity=IDENTITY_B,
    )
    entries = (
        EncodeRuntimeEntry.from_dict(
            {
                "path": "conda-envs/0123456789abcdef.env_setup_done",
                "kind": "file",
                "sha256": hashlib.sha256(b"").hexdigest(),
                "size": 0,
                "mode": 0o444,
                "target": None,
            }
        ),
        EncodeRuntimeEntry.from_dict(
            {
                "path": "runner/bin/conda",
                "kind": "file",
                "sha256": hashlib.sha256(b"conda").hexdigest(),
                "size": len(b"conda"),
                "mode": 0o555,
                "target": None,
            }
        ),
        EncodeRuntimeEntry.from_dict(
            {
                "path": "runner/bin/snakemake",
                "kind": "file",
                "sha256": hashlib.sha256(b"runner").hexdigest(),
                "size": len(b"runner"),
                "mode": 0o555,
                "target": None,
            }
        ),
    )
    receipt = EncodeRuntimePrepareReceipt.create(
        request_identity=request.identity,
        deployment_identity=request.deployment_identity,
        inventory=EncodeRuntimeInventory.create(entries),
    )
    monkeypatch.setattr(
        runtime_module,
        "_read_root_request",
        lambda _path: canonical_json_bytes(request.to_dict()),
    )
    monkeypatch.setattr(
        encode_runtime_materializer.OfflineEncodeRuntimeMaterializer,
        "prepare",
        lambda _self, observed: receipt if observed == request else None,
    )

    assert (
        candidate_service_main(
            ("encode-runtime-prepare",),
            release_root=Path(f"/tmp/{IDENTITY_A}"),
        )
        == 0
    )

    captured = capsys.readouterr()
    assert captured.err == ""
    assert EncodeRuntimePrepareReceipt.from_dict(json.loads(captured.out)) == receipt


def test_candidate_service_runs_the_fixed_offline_bulk_materializer(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    from encode_pipeline.deployment import bulk_runtime_materializer

    request = BulkRuntimePrepareRequest.create(
        operation="activate",
        task_identity=TASK,
        candidate_bulk_identity=IDENTITY_B,
        authority_platform_identity=IDENTITY_A,
        prior_state_identity=IDENTITY_C,
        candidate_state_identity=IDENTITY_B,
        docker_service_identity=IDENTITY_A,
        docker_client_identity=IDENTITY_B,
        docker_endpoint_identity=IDENTITY_C,
        docker_daemon_uid=os.getuid(),
        docker_daemon_gid=os.getgid(),
    )
    receipt = BulkRuntimePrepareReceipt.create(
        request_identity=request.identity,
        candidate_bulk_identity=request.candidate_bulk_identity,
        runtime_identity=IDENTITY_A,
        image_set_identity=IDENTITY_C,
        image_count=2,
    )
    monkeypatch.setattr(
        runtime_module,
        "_read_root_request",
        lambda path: (
            canonical_json_bytes(request.to_dict())
            if path == runtime_module.BULK_RUNTIME_PREPARE_REQUEST_PATH
            else b""
        ),
    )
    monkeypatch.setattr(
        bulk_runtime_materializer.OfflineBulkRuntimeMaterializer,
        "prepare",
        lambda _self, observed: receipt if observed == request else None,
    )

    assert (
        candidate_service_main(
            ("bulk-runtime-prepare",),
            release_root=Path(f"/tmp/{IDENTITY_A}"),
        )
        == 0
    )

    captured = capsys.readouterr()
    assert captured.err == ""
    assert BulkRuntimePrepareReceipt.from_dict(json.loads(captured.out)) == receipt


def test_candidate_launcher_only_forwards_api_secrets() -> None:
    content = (
        Path(runtime_module.__file__).parent / "templates" / "helixweave-service"
    ).read_text()

    assert '"ENCODE_AGENT_LLM_PROVIDER"' in content
    assert '"ENCODE_AGENT_LLM_API_KEY"' in content
    assert '"ENCODE_AGENT_LLM_MODEL"' in content
    assert 'API_SERVICE_ENVIRONMENT if command == "api" else set()' in content
    assert '"PYTHONDONTWRITEBYTECODE": "1"' in content


def test_candidate_action_emits_the_single_canonical_receipt(
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    request = _action_request()
    receipt = _ready_action_receipt(request)
    monkeypatch.setattr(
        runtime_module,
        "_read_root_request",
        lambda _path: canonical_json_bytes(request.to_dict()),
    )

    assert (
        candidate_service_main(
            ("operator-action",),
            release_root=Path(f"/tmp/{IDENTITY_A}"),
            action_probe=lambda observed, _root: (
                receipt if observed == request else None
            ),
        )
        == 0
    )

    captured = capsys.readouterr()
    assert captured.err == ""
    assert DeploymentActionReceipt.from_dict(json.loads(captured.out)) == receipt


@pytest.mark.parametrize(
    ("probe", "response", "request_prefix"),
    (
        (runtime_module._redis_ping, b"+PONG\r\n", b"*1\r\n$4\r\nPING\r\n"),
        (
            runtime_module._docker_ping,
            b"HTTP/1.0 200 OK\r\nContent-Length: 2\r\n\r\nOK",
            b"GET /_ping HTTP/1.0\r\n",
        ),
    ),
)
def test_runtime_readiness_uses_real_bounded_socket_protocols(
    monkeypatch: pytest.MonkeyPatch,
    probe,
    response: bytes,
    request_prefix: bytes,
) -> None:
    class FakePath:
        def lstat(self):
            return SimpleNamespace(
                st_mode=stat.S_IFSOCK | 0o660,
                st_dev=11,
                st_ino=12,
            )

        def __str__(self) -> str:
            return "/fixed/service.sock"

    class FakeSocket:
        def __init__(self) -> None:
            self.sent = b""
            self.responses = [response, b""]

        def __enter__(self):
            return self

        def __exit__(self, *_args) -> None:
            return None

        def settimeout(self, value: float) -> None:
            assert value == 2.0

        def connect(self, path: str) -> None:
            assert path == "/fixed/service.sock"

        def sendall(self, content: bytes) -> None:
            self.sent = content

        def recv(self, _maximum: int) -> bytes:
            return self.responses.pop(0)

    client = FakeSocket()
    monkeypatch.setattr(runtime_module.socket, "socket", lambda *_args: client)

    readiness = probe(FakePath())

    assert readiness.status == "ready"
    assert readiness.reason_code == "READY"
    assert readiness.identity is not None
    assert client.sent.startswith(request_prefix)


def test_launcher_uses_fixed_candidate_coordinate_and_isolated_python() -> None:
    content = (
        Path(runtime_module.__file__).parent / "templates" / "helixweave-service"
    ).read_text(encoding="utf-8")

    assert 'PYTHON = "/usr/bin/python3.12"' in content
    assert '"operator-action"' in content
    assert '"db-prepare"' in content
    assert '"encode-runtime-prepare"' in content
    assert '"bulk-runtime-prepare"' in content
    assert '"python3.12"' in content
    assert '"site-packages"' in content
    assert '"-I", "-S", "-c"' in content
    assert '"HELIXWEAVE_ENCODE_RUNNER_ROOT"' in content
    assert '"HELIXWEAVE_ENCODE_CONDA_PREFIX"' in content
    assert "dict(os.environ)" not in content
    assert "shell=True" not in content
    assert "subprocess" not in content
    assert "PYTHONPATH" not in content


def test_candidate_launcher_dispatches_bulk_runtime_prepare(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    template = Path(runtime_module.__file__).parent / "templates" / "helixweave-service"
    namespace = runpy.run_path(str(template), run_name="helixweave_service_test")
    main = namespace["main"]
    release = Path(f"/opt/helixweave/releases/platform/{IDENTITY_A}")
    site_packages = release / "payload/platform/lib/python3.12/site-packages"
    main.__globals__["_fixed_release"] = lambda: (release, site_packages)
    main.__globals__["Path"] = lambda _value: SimpleNamespace(
        stat=lambda *, follow_symlinks: SimpleNamespace(
            st_mode=stat.S_IFREG | 0o555,
            st_uid=0,
        )
    )
    observed: dict[str, object] = {}

    def fake_execve(executable, argv, environment) -> None:
        observed.update(
            executable=executable,
            argv=argv,
            environment=environment,
        )

    monkeypatch.setattr(namespace["os"], "execve", fake_execve)

    assert main(("bulk-runtime-prepare",)) == 70
    assert observed["executable"] == "/usr/bin/python3.12"
    assert "bulk-runtime-prepare" in observed["argv"][-1]
    assert "PYTHONPATH" not in observed["environment"]


def test_lock_rejects_duplicate_keys_and_noncanonical_bytes() -> None:
    with pytest.raises(DeploymentError):
        PlatformWheelLock.from_bytes(
            b'{"python":{},"python":{},"schema_version":"x","wheels":[]}\n'
        )

    lock = PlatformWheelLock((LockedWheel("x-1.0-py3-none-any.whl", "0" * 64),))
    with pytest.raises(DeploymentError):
        PlatformWheelLock.from_bytes(
            json.dumps(lock.to_dict(), indent=2, sort_keys=True).encode()
        )


def test_collected_closure_detects_mode_and_content_tampering(tmp_path: Path) -> None:
    closure, destination, _lock = _valid_closure(tmp_path)
    target = closure.source_for(PLATFORM_RUNTIME_LAUNCHER_PATH)
    target.chmod(0o755)

    with pytest.raises(DeploymentError) as caught:
        inspect_platform_runtime_closure(destination)

    assert caught.value.issue.code == "PLATFORM_RUNTIME_CLOSURE_INVALID"
    target.chmod(0o555)
    wheel = next(
        item.source for item in closure.files if "/wheelhouse/" in item.logical_path
    )
    wheel.chmod(0o644)
    wheel.write_bytes(wheel.read_bytes() + b"tamper")
    wheel.chmod(0o444)
    with pytest.raises(DeploymentError):
        inspect_platform_runtime_closure(destination)


def test_inspector_rederives_site_and_launcher_but_ignores_other_bundle_contracts(
    tmp_path: Path,
) -> None:
    closure, destination, _lock = _valid_closure(tmp_path)
    destination.chmod(0o755)
    (destination / "manifest.json").write_bytes(b"unrelated bundle manifest\n")
    (destination / "manifest.json").chmod(0o444)
    contract_root = destination / "payload" / "contracts" / "platform"
    contract_root.chmod(0o755)
    (contract_root / "unrelated-contract.json").write_bytes(b"{}\n")
    (contract_root / "unrelated-contract.json").chmod(0o444)
    contract_root.chmod(0o555)
    destination.chmod(0o555)

    assert inspect_platform_runtime_closure(destination) == closure

    launcher = closure.source_for(PLATFORM_RUNTIME_LAUNCHER_PATH)
    original_launcher = launcher.read_bytes()
    launcher.chmod(0o755)
    launcher.write_bytes(original_launcher + b"# tampered\n")
    launcher.chmod(0o555)
    with pytest.raises(DeploymentError):
        inspect_platform_runtime_closure(destination)
    launcher.chmod(0o755)
    launcher.write_bytes(original_launcher)
    launcher.chmod(0o555)

    site_file = closure.source_for(
        f"{PLATFORM_RUNTIME_SITE_PACKAGES}/helixweave/__init__.py"
    )
    site_file.chmod(0o644)
    site_file.write_bytes(b"VALUE = 2\n")
    site_file.chmod(0o444)
    with pytest.raises(DeploymentError):
        inspect_platform_runtime_closure(destination)


@pytest.mark.parametrize("kind", ("file", "directory"))
def test_inspector_rejects_extra_members_inside_fixed_closure_roots(
    tmp_path: Path,
    kind: str,
) -> None:
    _closure, destination, _lock = _valid_closure(tmp_path)
    bin_root = destination / "payload" / "platform" / "bin"
    bin_root.chmod(0o755)
    unexpected = bin_root / "unexpected"
    if kind == "file":
        unexpected.write_bytes(b"unexpected\n")
        unexpected.chmod(0o444)
    else:
        unexpected.mkdir(mode=0o555)
    bin_root.chmod(0o555)

    with pytest.raises(DeploymentError):
        inspect_platform_runtime_closure(destination)


def test_wheel_input_must_be_regular_non_writable_file(tmp_path: Path) -> None:
    wheelhouse = tmp_path / "wheelhouse"
    wheelhouse.mkdir()
    wheel = _build_wheel(wheelhouse)
    lock_path = tmp_path / "wheel-lock.json"
    _write_lock(lock_path, (wheel,))
    wheel.chmod(0o666)

    with pytest.raises(DeploymentError):
        collect_platform_runtime_closure(wheelhouse, lock_path, tmp_path / "closure")

    assert not (tmp_path / "closure").exists()


def test_runtime_tree_contains_no_symlinks(tmp_path: Path) -> None:
    closure, destination, _lock = _valid_closure(tmp_path)
    assert closure.files
    assert not any(path.is_symlink() for path in destination.rglob("*"))
    assert os.access(destination, os.R_OK)
