from __future__ import annotations

import json
import os
from pathlib import Path
import signal
import subprocess
import sys

import pytest

import scripts.run_local_platform as local_platform
from scripts.run_local_platform import (
    EnvironmentCheck,
    RUNTIME_IMPORTS,
    RuntimeConfig,
    WorkflowCheck,
    _http_ready,
    _session_process_groups,
    build_parser,
    build_shared_environment,
    cleanup_service_sessions,
    config_from_args,
    prepare_runtime,
    run_environment_doctor,
    run_workflow_doctor,
)
from scripts.results_visibility_fixture import prepare_results_visibility_fixture
from platform_runtime import (
    cleanup_owned_runtime_root,
    prepare_bulk_authoring_fixture,
    prepare_owned_runtime_root,
    write_manifest,
)


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
OWNERSHIP_SENTINEL = ".encode-platform-playwright-owned"


def test_environment_doctor_accepts_the_locked_project_toolchain(tmp_path, monkeypatch):
    frontend_root = tmp_path / "frontend"
    (frontend_root / "node_modules" / ".bin").mkdir(parents=True)
    (frontend_root / "node_modules" / ".bin" / "vite").touch()
    (frontend_root / "node_modules" / ".bin" / "playwright").touch()
    (frontend_root / "package-lock.json").touch()
    monkeypatch.setattr(local_platform.sys, "version_info", (3, 12, 13))
    monkeypatch.setattr(
        local_platform,
        "_import_runtime_dependency",
        lambda module: None,
    )
    monkeypatch.setattr(
        local_platform.shutil,
        "which",
        lambda executable: f"/tools/{executable}",
    )
    versions = {
        "snakemake": "8.30.0\n",
        "redis-server": "Redis server v=7.4.7 sha=00000000 malloc=libc bits=64\n",
        "node": "v22.18.0\n",
        "npm": "10.9.3\n",
    }
    monkeypatch.setattr(
        local_platform,
        "_read_tool_version",
        lambda executable, _arguments: versions[executable],
    )
    monkeypatch.setattr(
        local_platform,
        "_read_reachable_redis_version",
        lambda _redis_url: (_ for _ in ()).throw(
            RuntimeError("the configured Redis server version is unavailable")
        ),
    )
    monkeypatch.setattr(local_platform, "_port_available", lambda _host, _port: True)

    checks = run_environment_doctor(frontend_root)

    assert checks == (
        EnvironmentCheck("Python", "3.12.13"),
        EnvironmentCheck("Python packages", "available"),
        EnvironmentCheck("Snakemake", "8.30.0"),
        EnvironmentCheck("Redis server", "7.4.7"),
        EnvironmentCheck("Node.js", "22.18.0"),
        EnvironmentCheck("npm", "10.9.3"),
        EnvironmentCheck("Frontend dependencies", "available"),
    )


def test_environment_doctor_requires_product_runtime_not_test_tools() -> None:
    assert "uvicorn" in RUNTIME_IMPORTS
    assert "pytest" not in RUNTIME_IMPORTS


def test_environment_doctor_does_not_require_playwright(tmp_path, monkeypatch) -> None:
    frontend_root = tmp_path / "frontend"
    (frontend_root / "node_modules" / ".bin").mkdir(parents=True)
    (frontend_root / "node_modules" / ".bin" / "vite").touch()
    (frontend_root / "package-lock.json").touch()
    monkeypatch.setattr(local_platform.sys, "version_info", (3, 12, 13))
    monkeypatch.setattr(
        local_platform,
        "_import_runtime_dependency",
        lambda module: None,
    )
    monkeypatch.setattr(
        local_platform.shutil,
        "which",
        lambda executable: f"/tools/{executable}",
    )
    versions = {
        "snakemake": "8.30.0\n",
        "redis-server": "Redis server v=7.4.7 sha=00000000 malloc=libc bits=64\n",
        "node": "v22.18.0\n",
        "npm": "10.9.3\n",
    }
    monkeypatch.setattr(
        local_platform,
        "_read_tool_version",
        lambda executable, _arguments: versions[executable],
    )
    monkeypatch.setattr(
        local_platform,
        "_read_reachable_redis_version",
        lambda _redis_url: (_ for _ in ()).throw(
            RuntimeError("the configured Redis server version is unavailable")
        ),
    )
    monkeypatch.setattr(local_platform, "_port_available", lambda _host, _port: True)

    checks = run_environment_doctor(frontend_root)

    assert EnvironmentCheck("Frontend dependencies", "available") in checks


def test_workflow_doctor_reports_both_products_without_optional_bulk_runtime() -> None:
    checks = run_workflow_doctor(REPOSITORY_ROOT, environ={})

    assert checks == (
        WorkflowCheck(
            workflow_id="encode-style-chipseq-cuttag-atac-mnase",
            name="ENCODE-style ChIP-seq/CUT&Tag/ATAC/MNase",
            authoring="available",
            execution="available",
            reason_code="WORKFLOW_EXECUTION_READY",
        ),
        WorkflowCheck(
            workflow_id="bulk-rnaseq",
            name="Bulk RNA-seq",
            authoring="available",
            execution="not_configured",
            reason_code="WORKFLOW_EXECUTION_NOT_CONFIGURED",
        ),
    )


def test_workflow_doctor_fails_closed_for_incomplete_trial_tree(tmp_path) -> None:
    with pytest.raises(RuntimeError, match="product information is unavailable"):
        run_workflow_doctor(tmp_path, environ={})


def test_workflow_doctor_redacts_partial_bulk_runtime_configuration() -> None:
    private_value = "/private/operator/runtime"
    checks = run_workflow_doctor(
        REPOSITORY_ROOT,
        environ={"HELIXWEAVE_BULK_RNASEQ_RUNTIME_ROOT": private_value},
    )

    bulk = checks[1]
    assert bulk.execution == "unavailable"
    assert bulk.reason_code == "WORKFLOW_EXECUTION_UNAVAILABLE"
    assert private_value not in repr(checks)


def test_workflow_doctor_redacts_registry_composition_failure(
    monkeypatch,
) -> None:
    private_detail = "/private/operator/runtime"

    def fail_registry(**_kwargs):
        raise RuntimeError(private_detail)

    monkeypatch.setattr(
        "encode_pipeline.services.defaults.create_default_workflow_registry",
        fail_registry,
    )

    with pytest.raises(
        RuntimeError, match="workflow product information is unavailable"
    ) as caught:
        run_workflow_doctor(REPOSITORY_ROOT, environ={})

    assert private_detail not in str(caught.value)


def test_workflow_doctor_rejects_unexpected_registry_identity(monkeypatch) -> None:
    class UnexpectedDescriptor:
        class metadata:
            workflow_id = "unexpected-workflow"

    class UnexpectedWorkflowInfo:
        def __init__(self, _registry):
            pass

        def list_descriptors(self):
            return (UnexpectedDescriptor(),)

    monkeypatch.setattr(
        "encode_pipeline.services.workflow_info.WorkflowInfoService",
        UnexpectedWorkflowInfo,
    )

    with pytest.raises(
        RuntimeError, match="workflow product information is unavailable"
    ):
        run_workflow_doctor(REPOSITORY_ROOT, environ={})


@pytest.mark.parametrize(
    ("status", "expected"), [(200, True), (204, True), (404, False)]
)
def test_http_readiness_requires_a_success_response(
    monkeypatch,
    status,
    expected,
) -> None:
    class Response:
        def __enter__(self):
            return self

        def __exit__(self, *_args):
            return None

        def __init__(self, response_status):
            self.status = response_status

    class Opener:
        def open(self, _url, *, timeout):
            assert timeout == 1
            return Response(status)

    monkeypatch.setattr(local_platform, "build_opener", lambda *_args: Opener())

    assert _http_ready("http://127.0.0.1:8000/ready") is expected


@pytest.mark.parametrize(
    "redis_url",
    [
        "rediss://127.0.0.1:6379/0",
        "redis://private-host:6379/0",
        "redis://user@localhost:6379/0",
        "redis://:secret@localhost:6379/0",
        "redis://localhost:6379/1",
        "redis://localhost:6379/0?private=true",
        "redis://localhost:6379/0#private",
        "redis://localhost:0/0",
        "redis://localhost:70000/0",
    ],
)
def test_automatic_redis_start_target_rejects_unsafe_coordinates(redis_url) -> None:
    assert local_platform._automatic_redis_start_target(redis_url) is None


def test_automatic_redis_start_target_defaults_only_the_missing_port() -> None:
    assert local_platform._automatic_redis_start_target("redis://localhost/0") == (
        "localhost",
        6379,
    )


def test_environment_doctor_accepts_a_versioned_reachable_redis_when_binary_is_absent(
    tmp_path, monkeypatch
):
    frontend_root = tmp_path / "frontend"
    (frontend_root / "node_modules" / ".bin").mkdir(parents=True)
    (frontend_root / "node_modules" / ".bin" / "vite").touch()
    (frontend_root / "node_modules" / ".bin" / "playwright").touch()
    (frontend_root / "package-lock.json").touch()
    monkeypatch.setattr(local_platform.sys, "version_info", (3, 12, 13))
    monkeypatch.setattr(
        local_platform,
        "_import_runtime_dependency",
        lambda module: None,
    )
    versions = {
        "snakemake": "8.30.0\n",
        "node": "v22.18.0\n",
        "npm": "10.9.3\n",
    }

    def tool_version(executable, _arguments):
        if executable == "redis-server":
            raise RuntimeError("redis-server is unavailable on PATH")
        return versions[executable]

    monkeypatch.setattr(local_platform, "_read_tool_version", tool_version)
    monkeypatch.setattr(
        local_platform.shutil,
        "which",
        lambda executable: (
            None if executable == "redis-server" else f"/tools/{executable}"
        ),
    )
    monkeypatch.setattr(
        local_platform,
        "_read_reachable_redis_version",
        lambda _redis_url: (7, 4, 9),
        raising=False,
    )

    checks = run_environment_doctor(
        frontend_root,
        redis_url="redis://127.0.0.1:6379/14",
    )

    assert EnvironmentCheck("Redis server", "7.4.9") in checks


def test_environment_doctor_external_redis_failure_does_not_expose_connection(
    tmp_path, monkeypatch
):
    frontend_root = tmp_path / "frontend"
    (frontend_root / "node_modules" / ".bin").mkdir(parents=True)
    (frontend_root / "node_modules" / ".bin" / "vite").touch()
    (frontend_root / "node_modules" / ".bin" / "playwright").touch()
    (frontend_root / "package-lock.json").touch()
    monkeypatch.setattr(local_platform.sys, "version_info", (3, 12, 13))
    monkeypatch.setattr(
        local_platform,
        "_import_runtime_dependency",
        lambda module: None,
    )
    versions = {
        "snakemake": "8.30.0\n",
        "node": "v22.18.0\n",
        "npm": "10.9.3\n",
    }

    def tool_version(executable, _arguments):
        if executable == "redis-server":
            raise RuntimeError("redis-server is unavailable on PATH")
        return versions[executable]

    def unavailable_redis(_redis_url):
        raise RuntimeError("redis://user:secret@private-host:6379/0")

    monkeypatch.setattr(local_platform, "_read_tool_version", tool_version)
    monkeypatch.setattr(
        local_platform.shutil,
        "which",
        lambda executable: (
            None if executable == "redis-server" else f"/tools/{executable}"
        ),
    )
    monkeypatch.setattr(
        local_platform,
        "_read_reachable_redis_version",
        unavailable_redis,
    )

    with pytest.raises(RuntimeError) as caught:
        run_environment_doctor(
            frontend_root,
            redis_url="redis://user:secret@private-host:6379/0",
        )

    assert str(caught.value) == "the configured Redis server version is unavailable"
    assert "secret" not in str(caught.value)
    assert "private-host" not in str(caught.value)


def test_environment_doctor_prefers_the_configured_reachable_redis(
    tmp_path, monkeypatch
):
    frontend_root = tmp_path / "frontend"
    (frontend_root / "node_modules" / ".bin").mkdir(parents=True)
    (frontend_root / "node_modules" / ".bin" / "vite").touch()
    (frontend_root / "node_modules" / ".bin" / "playwright").touch()
    (frontend_root / "package-lock.json").touch()
    monkeypatch.setattr(local_platform.sys, "version_info", (3, 12, 13))
    monkeypatch.setattr(
        local_platform,
        "_import_runtime_dependency",
        lambda module: None,
    )
    monkeypatch.setattr(
        local_platform.shutil,
        "which",
        lambda executable: f"/tools/{executable}",
    )
    versions = {
        "snakemake": "8.30.0\n",
        "node": "v22.18.0\n",
        "npm": "10.9.3\n",
    }

    def tool_version(executable, _arguments):
        if executable == "redis-server":
            raise RuntimeError("redis-server could not report a version")
        return versions[executable]

    monkeypatch.setattr(local_platform, "_read_tool_version", tool_version)
    monkeypatch.setattr(
        local_platform,
        "_read_reachable_redis_version",
        lambda _redis_url: (7, 4, 9),
    )

    checks = run_environment_doctor(
        frontend_root,
        redis_url="redis://127.0.0.1:6379/14",
    )

    assert EnvironmentCheck("Redis server", "7.4.9") in checks


def test_environment_doctor_rejects_unreachable_custom_redis_even_with_binary(
    tmp_path, monkeypatch
):
    frontend_root = tmp_path / "frontend"
    (frontend_root / "node_modules" / ".bin").mkdir(parents=True)
    (frontend_root / "node_modules" / ".bin" / "vite").touch()
    (frontend_root / "package-lock.json").touch()
    monkeypatch.setattr(local_platform.sys, "version_info", (3, 12, 13))
    monkeypatch.setattr(
        local_platform, "_import_runtime_dependency", lambda _module: None
    )
    monkeypatch.setattr(
        local_platform.shutil,
        "which",
        lambda executable: f"/tools/{executable}",
    )
    versions = {
        "snakemake": "8.30.0\n",
        "redis-server": "Redis server v=7.4.7 sha=00000000 malloc=libc bits=64\n",
        "node": "v22.18.0\n",
        "npm": "10.9.3\n",
    }
    binary_calls: list[str] = []

    def tool_version(executable, _arguments):
        binary_calls.append(executable)
        return versions[executable]

    monkeypatch.setattr(local_platform, "_read_tool_version", tool_version)
    monkeypatch.setattr(
        local_platform,
        "_read_reachable_redis_version",
        lambda _redis_url: (_ for _ in ()).throw(
            RuntimeError("redis://private-host:6379/14")
        ),
    )

    with pytest.raises(
        RuntimeError, match="configured Redis server version is unavailable"
    ) as caught:
        run_environment_doctor(
            frontend_root,
            redis_url="redis://private-host:6379/14",
        )

    assert "redis-server" not in binary_calls
    assert "private-host" not in str(caught.value)


def test_environment_doctor_does_not_replace_an_old_reachable_redis_with_binary(
    tmp_path, monkeypatch
):
    frontend_root = tmp_path / "frontend"
    (frontend_root / "node_modules" / ".bin").mkdir(parents=True)
    (frontend_root / "node_modules" / ".bin" / "vite").touch()
    (frontend_root / "package-lock.json").touch()
    monkeypatch.setattr(local_platform.sys, "version_info", (3, 12, 13))
    monkeypatch.setattr(
        local_platform, "_import_runtime_dependency", lambda _module: None
    )
    monkeypatch.setattr(
        local_platform.shutil,
        "which",
        lambda executable: f"/tools/{executable}",
    )
    versions = {
        "snakemake": "8.30.0\n",
        "redis-server": "Redis server v=7.4.7 sha=00000000 malloc=libc bits=64\n",
        "node": "v22.18.0\n",
        "npm": "10.9.3\n",
    }
    binary_calls: list[str] = []

    def tool_version(executable, _arguments):
        binary_calls.append(executable)
        return versions[executable]

    monkeypatch.setattr(local_platform, "_read_tool_version", tool_version)
    monkeypatch.setattr(
        local_platform,
        "_read_reachable_redis_version",
        lambda _redis_url: (6, 2, 14),
    )

    with pytest.raises(RuntimeError, match="Redis server 7.x is required"):
        run_environment_doctor(frontend_root)

    assert "redis-server" not in binary_calls


def test_reachable_redis_version_wraps_malformed_url_without_disclosure(
    monkeypatch,
):
    def malformed_url(*_args, **_kwargs):
        raise ValueError("redis://user:secret@private-host:not-a-port/0")

    monkeypatch.setattr("redis.Redis.from_url", malformed_url)

    with pytest.raises(RuntimeError) as caught:
        local_platform._read_reachable_redis_version(
            "redis://user:secret@private-host:not-a-port/0"
        )

    assert str(caught.value) == "the configured Redis server version is unavailable"
    assert "secret" not in str(caught.value)
    assert "private-host" not in str(caught.value)


def test_reachable_redis_version_rejects_malformed_info_and_closes_connection(
    monkeypatch,
):
    class MalformedInfoConnection:
        closed = False

        def info(self, *, section):
            assert section == "server"
            return []

        def close(self):
            self.closed = True

    connection = MalformedInfoConnection()
    monkeypatch.setattr(
        "redis.Redis.from_url",
        lambda *_args, **_kwargs: connection,
    )

    with pytest.raises(RuntimeError) as caught:
        local_platform._read_reachable_redis_version("redis://private-host:6379/0")

    assert str(caught.value) == "the configured Redis server version is unavailable"
    assert connection.closed is True
    assert "private-host" not in str(caught.value)


@pytest.mark.parametrize(
    ("overrides", "message"),
    [
        ({"python": (3, 11, 9)}, "Python 3.12"),
        ({"snakemake": "8.29.0\n"}, "Snakemake 8.30.0"),
        (
            {
                "redis-server": (
                    "Redis server v=6.2.14 sha=00000000 malloc=libc bits=64\n"
                )
            },
            "Redis server 7.x",
        ),
        ({"node": "v18.20.0\n"}, "Node.js 20 or newer"),
    ],
)
def test_environment_doctor_rejects_incompatible_tools(
    tmp_path, monkeypatch, overrides, message
):
    frontend_root = tmp_path / "frontend"
    (frontend_root / "node_modules" / ".bin").mkdir(parents=True)
    (frontend_root / "node_modules" / ".bin" / "vite").touch()
    (frontend_root / "node_modules" / ".bin" / "playwright").touch()
    (frontend_root / "package-lock.json").touch()
    monkeypatch.setattr(
        local_platform.sys,
        "version_info",
        overrides.get("python", (3, 12, 13)),
    )
    monkeypatch.setattr(
        local_platform,
        "_import_runtime_dependency",
        lambda module: None,
    )
    monkeypatch.setattr(
        local_platform.shutil,
        "which",
        lambda executable: f"/tools/{executable}",
    )
    versions = {
        "snakemake": "8.30.0\n",
        "redis-server": "Redis server v=7.4.7 sha=00000000 malloc=libc bits=64\n",
        "node": "v22.18.0\n",
        "npm": "10.9.3\n",
    }
    versions.update({key: value for key, value in overrides.items() if key != "python"})
    monkeypatch.setattr(
        local_platform,
        "_read_tool_version",
        lambda executable, _arguments: versions[executable],
    )
    monkeypatch.setattr(
        local_platform,
        "_read_reachable_redis_version",
        lambda _redis_url: (_ for _ in ()).throw(
            RuntimeError("the configured Redis server version is unavailable")
        ),
    )
    monkeypatch.setattr(local_platform, "_port_available", lambda _host, _port: True)

    with pytest.raises(RuntimeError, match=message):
        run_environment_doctor(frontend_root)


def test_environment_doctor_reports_missing_packages_without_exception_details(
    tmp_path, monkeypatch
):
    frontend_root = tmp_path / "frontend"
    private_detail = "private-package-location"
    monkeypatch.setattr(local_platform.sys, "version_info", (3, 12, 13))

    def missing_dependency(module):
        raise ModuleNotFoundError(private_detail)

    monkeypatch.setattr(
        local_platform,
        "_import_runtime_dependency",
        missing_dependency,
    )

    with pytest.raises(RuntimeError) as caught:
        run_environment_doctor(frontend_root)

    message = str(caught.value)
    assert "Python package fastapi is unavailable" in message
    assert private_detail not in message


def test_doctor_mode_has_no_port_runtime_or_process_side_effects(
    tmp_path, monkeypatch, capsys
):
    frontend_root = tmp_path / "frontend"
    runtime_root = tmp_path / "runtime"
    calls: list[str] = []
    monkeypatch.setattr(
        local_platform,
        "run_environment_doctor",
        lambda _root, *, redis_url: (
            calls.append("doctor") or (EnvironmentCheck("Python", "3.12.13"),)
        ),
    )
    monkeypatch.setattr(
        local_platform,
        "run_workflow_doctor",
        lambda _root, *, environ: (
            calls.append("workflows")
            or (
                WorkflowCheck(
                    "encode-style-chipseq-cuttag-atac-mnase",
                    "ENCODE-style ChIP-seq/CUT&Tag/ATAC/MNase",
                    "available",
                    "available",
                    "WORKFLOW_EXECUTION_READY",
                ),
                WorkflowCheck(
                    "bulk-rnaseq",
                    "Bulk RNA-seq",
                    "available",
                    "not_configured",
                    "WORKFLOW_EXECUTION_NOT_CONFIGURED",
                ),
            )
        ),
    )
    monkeypatch.setattr(
        local_platform,
        "_port_available",
        lambda _host, _port: calls.append("port") or True,
    )

    assert (
        local_platform.main(
            [
                "--doctor",
                "--frontend-root",
                str(frontend_root),
                "--runtime-root",
                str(runtime_root),
            ]
        )
        == 0
    )

    assert calls == ["doctor", "workflows"]
    assert not runtime_root.exists()
    output = capsys.readouterr().out
    assert "HelixWeave environment check" in output
    assert "Python: 3.12.13" in output
    assert "bulk-rnaseq: authoring=available execution=not_configured" in output
    assert "Bulk RNA-seq execution is optional" in output
    assert "private" not in output.lower()


def test_launcher_help_uses_product_brand_and_keeps_compatibility_flags():
    help_text = build_parser().format_help()

    assert "HelixWeave local stack" in help_text
    assert "--doctor" in help_text
    assert "--input-authoring-demo" in help_text
    assert "--results-visibility-demo" in help_text


def test_shared_environment_uses_one_absolute_lifecycle_configuration(tmp_path):
    project_root = tmp_path / "project"
    config = RuntimeConfig(
        project_root=project_root,
        frontend_root=tmp_path / "frontend",
        runtime_root=tmp_path / "runtime",
        redis_url="redis://127.0.0.1:6380/0",
        queue_name="browser-e2e",
        api_host="127.0.0.1",
        api_port=8010,
        frontend_host="127.0.0.1",
        frontend_port=4173,
        readiness_timeout=30,
    )

    environment = build_shared_environment(
        config, {"PATH": "/tools", "NO_PROXY": "127.0.0.1"}
    )

    assert environment["ENCODE_PIPELINE_DATABASE_URL"] == (
        f"sqlite:///{tmp_path / 'runtime' / 'platform.db'}"
    )
    assert environment["ENCODE_PIPELINE_WORKSPACE_ROOT"] == str(
        tmp_path / "runtime" / "workspaces"
    )
    assert environment["ENCODE_PIPELINE_REDIS_URL"] == config.redis_url
    assert environment["ENCODE_PIPELINE_QUEUE_NAME"] == "browser-e2e"
    assert environment["PYTHONPATH"] == str(project_root / "src")
    assert environment["VITE_API_PROXY_TARGET"] == "http://127.0.0.1:8010"
    assert "127.0.0.1" in environment["NO_PROXY"]
    assert environment["NO_PROXY"].split(",").count("127.0.0.1") == 1
    assert environment["NO_PROXY"] == environment["no_proxy"]


def test_supervisor_exposes_frontend_only_after_api_and_worker_are_ready(
    tmp_path, monkeypatch
):
    config = RuntimeConfig(
        project_root=tmp_path / "project",
        frontend_root=tmp_path / "frontend",
        runtime_root=tmp_path / "runtime",
        redis_url="redis://127.0.0.1:6380/0",
        queue_name="browser-e2e",
        api_host="127.0.0.1",
        api_port=8010,
        frontend_host="127.0.0.1",
        frontend_port=4173,
        readiness_timeout=30,
    )
    supervisor = local_platform.PlatformSupervisor(config)
    order: list[str] = []
    monkeypatch.setattr(
        supervisor,
        "_start",
        lambda name, _argv, *, cwd: order.append(name),
    )
    monkeypatch.setattr(
        supervisor,
        "_wait_backend_ready",
        lambda: order.append("backend-ready"),
        raising=False,
    )

    supervisor._start_services()

    assert order == ["api", "worker", "backend-ready", "frontend"]


def test_supervisor_starts_only_the_bounded_configured_local_redis(
    tmp_path, monkeypatch
):
    config = RuntimeConfig(
        project_root=tmp_path / "project",
        frontend_root=tmp_path / "frontend",
        runtime_root=tmp_path / "runtime",
        redis_url="redis://localhost:6380/0",
        queue_name="browser-e2e",
        api_host="127.0.0.1",
        api_port=8010,
        frontend_host="127.0.0.1",
        frontend_port=4173,
        readiness_timeout=30,
    )
    supervisor = local_platform.PlatformSupervisor(config)
    redis_ready = iter((False, True))
    monkeypatch.setattr(
        local_platform,
        "_redis_ready",
        lambda _redis_url: next(redis_ready),
    )
    starts: list[tuple[str, list[str], Path]] = []
    monkeypatch.setattr(
        supervisor,
        "_start",
        lambda name, argv, *, cwd: starts.append((name, argv, cwd)),
    )

    def wait_until(condition, timeout, description, check_processes):
        assert timeout == 30
        assert description == "Redis"
        assert check_processes == supervisor._assert_processes_alive
        assert condition() is True

    monkeypatch.setattr(local_platform, "_wait_until", wait_until)

    supervisor._ensure_redis()

    assert starts == [
        (
            "redis",
            [
                "redis-server",
                "--bind",
                "localhost",
                "--port",
                "6380",
                "--dir",
                str(config.runtime_root / "redis"),
                "--save",
                "",
                "--appendonly",
                "no",
            ],
            config.runtime_root,
        )
    ]
    assert (config.runtime_root / "redis").is_dir()


def test_supervisor_refuses_to_start_redis_for_an_unsafe_coordinate(
    tmp_path, monkeypatch
):
    config = RuntimeConfig(
        project_root=tmp_path / "project",
        frontend_root=tmp_path / "frontend",
        runtime_root=tmp_path / "runtime",
        redis_url="redis://private-host:6380/0",
        queue_name="browser-e2e",
        api_host="127.0.0.1",
        api_port=8010,
        frontend_host="127.0.0.1",
        frontend_port=4173,
        readiness_timeout=30,
    )
    supervisor = local_platform.PlatformSupervisor(config)
    monkeypatch.setattr(local_platform, "_redis_ready", lambda _redis_url: False)

    with pytest.raises(RuntimeError, match="unavailable and cannot be started locally"):
        supervisor._ensure_redis()

    assert not config.runtime_root.exists()


def test_session_process_groups_never_include_the_test_process_group(monkeypatch):
    own_group = 88
    monkeypatch.setattr("scripts.run_local_platform.Path.iterdir", lambda _path: ())
    monkeypatch.setattr("scripts.run_local_platform.os.getpgrp", lambda: own_group)

    groups = _session_process_groups(99)

    assert groups == ()
    assert own_group not in groups
    assert signal.SIGKILL.value > 0


def test_cleanup_service_sessions_signals_every_process_group(tmp_path, monkeypatch):
    session_id = 10001
    horse_group = 10002
    pid_file = tmp_path / "service-pids.json"
    pid_file.write_text(json.dumps({"worker": session_id}), encoding="utf-8")
    scans = iter(((horse_group, session_id), ()))
    monkeypatch.setattr(
        "scripts.run_local_platform._session_process_groups",
        lambda _session_id: next(scans, ()),
    )
    signals: list[tuple[int, signal.Signals]] = []
    monkeypatch.setattr(
        "scripts.run_local_platform.os.killpg",
        lambda group, signum: signals.append((group, signum)),
    )

    assert cleanup_service_sessions(pid_file) == 0
    assert (horse_group, signal.SIGTERM) in signals
    assert (session_id, signal.SIGTERM) in signals


def test_playwright_runtime_reset_rejects_a_mismatched_owner(tmp_path):
    runtime_root = tmp_path / "runtime"
    runtime_root.mkdir()
    unrelated = runtime_root / "unrelated.txt"
    unrelated.write_text("keep", encoding="utf-8")
    (runtime_root / OWNERSHIP_SENTINEL).write_text("owner-a", encoding="utf-8")
    environment = dict(os.environ)
    environment.update(
        {
            "ENCODE_PIPELINE_E2E_ROOT": str(runtime_root),
            "ENCODE_PIPELINE_E2E_OWNER": "owner-b",
        }
    )

    completed = subprocess.run(
        [sys.executable, str(REPOSITORY_ROOT / "test/browser/platform_runtime.py")],
        cwd=REPOSITORY_ROOT,
        env=environment,
        capture_output=True,
        text=True,
        timeout=10,
        check=False,
    )

    assert completed.returncode != 0
    assert "not owned" in completed.stderr
    assert unrelated.read_text(encoding="utf-8") == "keep"


def test_playwright_runtime_reset_repairs_read_only_owned_directories(tmp_path):
    runtime_root = tmp_path / "runtime"
    read_only_directory = runtime_root / "workspaces" / "attempt" / "sortmerna"
    read_only_directory.mkdir(parents=True)
    stale_index = read_only_directory / "index.bursttrie_0.dat"
    stale_index.write_bytes(b"sealed index")
    (runtime_root / OWNERSHIP_SENTINEL).write_text("owner-a", encoding="utf-8")
    read_only_directory.chmod(0o550)
    read_only_directory.parent.chmod(0o550)

    prepared = prepare_owned_runtime_root(str(runtime_root), "owner-a")

    assert prepared == runtime_root
    assert prepared.is_dir()
    assert (prepared / OWNERSHIP_SENTINEL).read_text(encoding="utf-8") == "owner-a"
    assert not stale_index.exists()


def test_playwright_runtime_cleanup_repairs_owned_tree_without_following_symlinks(
    tmp_path,
):
    runtime_root = tmp_path / "runtime"
    read_only_directory = runtime_root / "workspaces" / "attempt" / "sortmerna"
    read_only_directory.mkdir(parents=True)
    (read_only_directory / "index.bursttrie_0.dat").write_bytes(b"sealed index")
    outside = tmp_path / "outside"
    outside.mkdir()
    outside_sentinel = outside / "keep.txt"
    outside_sentinel.write_text("keep", encoding="utf-8")
    (runtime_root / "outside-link").symlink_to(outside, target_is_directory=True)
    (runtime_root / OWNERSHIP_SENTINEL).write_text("owner-a", encoding="utf-8")
    read_only_directory.chmod(0o550)
    read_only_directory.parent.chmod(0o550)

    cleanup_owned_runtime_root(str(runtime_root), "owner-a")

    assert not runtime_root.exists()
    assert outside_sentinel.read_text(encoding="utf-8") == "keep"


def test_playwright_runtime_cleanup_rejects_a_symlinked_root(tmp_path):
    owned_runtime = tmp_path / "owned-runtime"
    owned_runtime.mkdir()
    sentinel = owned_runtime / OWNERSHIP_SENTINEL
    sentinel.write_text("owner-a", encoding="utf-8")
    runtime_link = tmp_path / "runtime-link"
    runtime_link.symlink_to(owned_runtime, target_is_directory=True)

    with pytest.raises(ValueError, match="symbolic links"):
        cleanup_owned_runtime_root(str(runtime_link), "owner-a")

    assert sentinel.read_text(encoding="utf-8") == "owner-a"


def test_playwright_runtime_cleanup_cli_repairs_read_only_owned_directories(
    tmp_path,
):
    runtime_root = tmp_path / "runtime"
    read_only_directory = runtime_root / "workspaces" / "attempt" / "sortmerna"
    read_only_directory.mkdir(parents=True)
    (read_only_directory / "index.bursttrie_0.dat").write_bytes(b"sealed index")
    (runtime_root / OWNERSHIP_SENTINEL).write_text("owner-a", encoding="utf-8")
    read_only_directory.chmod(0o550)
    environment = dict(os.environ)
    environment.update(
        {
            "ENCODE_PIPELINE_E2E_ROOT": str(runtime_root),
            "ENCODE_PIPELINE_E2E_OWNER": "owner-a",
        }
    )

    completed = subprocess.run(
        [
            sys.executable,
            str(REPOSITORY_ROOT / "test/browser/platform_runtime.py"),
            "--cleanup-owned-runtime",
        ],
        cwd=REPOSITORY_ROOT,
        env=environment,
        capture_output=True,
        text=True,
        timeout=10,
        check=False,
    )

    assert completed.returncode == 0
    assert completed.stdout == ""
    assert completed.stderr == ""
    assert not runtime_root.exists()


def test_playwright_runtime_manifest_contains_only_controlled_fixture_inputs(
    tmp_path,
):
    runtime_root = tmp_path / "runtime"
    runtime_root.mkdir()
    inputs = prepare_results_visibility_fixture(runtime_root / "project")
    bulk_authoring = prepare_bulk_authoring_fixture(runtime_root)

    write_manifest(runtime_root, "browser-queue", inputs, bulk_authoring)

    raw = json.loads((runtime_root / "runtime.json").read_text(encoding="utf-8"))
    assert set(raw) == {
        "workflowId",
        "samplesPath",
        "resultsConfig",
        "cancelConfig",
        "emptyConfig",
        "malformedConfig",
        "expectedQcSummary",
        "runtimeRoot",
        "workspaceRoot",
        "markerRoot",
        "queueName",
        "bulkWorkflowId",
        "bulkSamplesPath",
        "bulkConfig",
        "bulkOptions",
        "bulkExpectedExecution",
        "bulkRequiredArtifactOutputTypes",
        "bulkRequiredQcMetricKeys",
    }
    assert raw["samplesPath"] == str(inputs.samples_path)
    assert raw["resultsConfig"] == inputs.results_config
    assert raw["cancelConfig"] == inputs.cancel_config
    assert raw["emptyConfig"] == inputs.empty_config
    assert raw["malformedConfig"] == inputs.malformed_config
    assert raw["expectedQcSummary"] == inputs.expected_qc_summary
    assert all(raw[name] == value for name, value in bulk_authoring.items())
    assert raw["bulkWorkflowId"] == "bulk-rnaseq"
    assert raw["bulkOptions"] == {}
    assert raw["bulkExpectedExecution"] == "not_configured"
    assert raw["bulkRequiredArtifactOutputTypes"] == []
    assert raw["bulkRequiredQcMetricKeys"] == []
    bulk_samples_path = Path(raw["bulkSamplesPath"])
    assert bulk_samples_path.is_relative_to(runtime_root)
    bulk_rows = bulk_samples_path.read_text(encoding="utf-8").splitlines()
    assert len(bulk_rows) == 4
    assert [row.split("\t")[0] for row in bulk_rows[1:]] == [
        "pairedA",
        "pairedA",
        "singleB",
    ]
    assert [row.split("\t")[3] for row in bulk_rows[1:]] == ["PE", "PE", "SE"]
    assert [row.split("\t")[6] for row in bulk_rows[1:]] == [
        "reverse",
        "reverse",
        "auto",
    ]
    assert all(row.split("\t")[7] == "ILLUMINA" for row in bulk_rows[1:])
    bulk_standard = raw["bulkConfig"]["standard"]
    assert bulk_standard["analysis"] == {
        "alignment": "star",
        "quantification": "salmon",
    }
    assert bulk_standard["trimming"] == {"enabled": True, "tool": "fastp"}
    assert bulk_standard["umi"]["mode"] == "read_name"
    assert bulk_standard["ribosomal_rna_removal"]["tool"] == "sortmerna"
    assert all(
        "samples" not in raw[name]
        for name in (
            "resultsConfig",
            "cancelConfig",
            "emptyConfig",
            "malformedConfig",
        )
    )
    serialized = json.dumps(raw)
    assert "sqlite:" not in serialized
    assert "redis:" not in serialized
    assert "ENCODE_PIPELINE_" not in serialized


def test_launcher_defaults_remain_ordinary_platform_runtime():
    args = build_parser().parse_args([])

    config = config_from_args(args)

    assert config.project_root == REPOSITORY_ROOT
    assert config.runtime_root == REPOSITORY_ROOT / ".local/platform-demo"
    assert config.queue_name == "encode-pipeline-demo"


def test_results_demo_prepares_isolated_project_and_bounded_inputs(tmp_path):
    runtime_root = tmp_path / "runtime" / "results-demo"
    args = build_parser().parse_args(
        [
            "--results-visibility-demo",
            "--runtime-root",
            str(runtime_root),
        ]
    )

    config, inputs_path = prepare_runtime(args)

    assert config.runtime_root == runtime_root
    assert config.project_root == runtime_root / "results-visibility-project"
    assert inputs_path == runtime_root / "results-visibility-inputs.json"
    raw = json.loads(inputs_path.read_text(encoding="utf-8"))
    assert set(raw) == {"workflowId", "samplesPath", "resultsConfig"}
    assert raw["workflowId"] == "encode-style-chipseq-cuttag-atac-mnase"
    assert raw["resultsConfig"]["threads"] == 1
    assert raw["resultsConfig"]["outdir"] == "results"
    assert "samples" not in raw["resultsConfig"]
    assert Path(raw["samplesPath"]) == config.project_root / "samples.tsv"
    serialized = json.dumps(raw)
    assert "sqlite:" not in serialized
    assert "redis:" not in serialized
    assert "ENCODE_PIPELINE_" not in serialized


def test_results_demo_has_dedicated_default_runtime_without_preparing_it():
    args = build_parser().parse_args(["--results-visibility-demo"])

    config = config_from_args(args)

    assert config.runtime_root == REPOSITORY_ROOT / ".local/results-visibility-demo"
    assert config.queue_name == "encode-pipeline-results-demo"


def test_input_authoring_demo_alias_uses_the_controlled_results_runtime():
    args = build_parser().parse_args(["--input-authoring-demo"])

    config = config_from_args(args)

    assert args.results_visibility_demo is True
    assert config.runtime_root == REPOSITORY_ROOT / ".local/results-visibility-demo"
    assert config.queue_name == "encode-pipeline-results-demo"


def test_results_demo_rejects_an_explicit_project_root(tmp_path):
    args = build_parser().parse_args(
        [
            "--results-visibility-demo",
            "--project-root",
            str(tmp_path / "project"),
            "--runtime-root",
            str(tmp_path / "runtime"),
        ]
    )

    with pytest.raises(ValueError, match="project-root"):
        prepare_runtime(args)


def test_main_checks_ports_then_prepares_demo_before_starting_processes(
    tmp_path,
    monkeypatch,
):
    order: list[str] = []
    config = RuntimeConfig(
        project_root=REPOSITORY_ROOT,
        frontend_root=REPOSITORY_ROOT / "frontend",
        runtime_root=tmp_path / "runtime",
        redis_url="redis://127.0.0.1:6380/0",
        queue_name="results-demo",
        api_host="127.0.0.1",
        api_port=8010,
        frontend_host="127.0.0.1",
        frontend_port=4173,
        readiness_timeout=30,
    )

    def configured(_args):
        order.append("configure")
        return config

    def prepared(_args):
        order.append("prepare")
        return config, None

    def available(_host, _port):
        order.append("port")
        return True

    def doctor(_frontend_root, *, redis_url):
        order.append("doctor")
        return (EnvironmentCheck("Python", "3.12.13"),)

    class Supervisor:
        def __init__(self, _config):
            order.append("process")

        def run(self):
            return 0

    monkeypatch.setattr(local_platform, "config_from_args", configured)
    monkeypatch.setattr(local_platform, "prepare_runtime", prepared)
    monkeypatch.setattr(local_platform, "run_environment_doctor", doctor)
    monkeypatch.setattr(
        local_platform,
        "run_workflow_doctor",
        lambda _root, *, environ: order.append("workflows") or (),
    )
    monkeypatch.setattr(local_platform, "_port_available", available)
    monkeypatch.setattr(local_platform, "PlatformSupervisor", Supervisor)

    assert (
        local_platform.main(
            [
                "--results-visibility-demo",
                "--runtime-root",
                str(tmp_path / "runtime"),
            ]
        )
        == 0
    )

    assert order == [
        "doctor",
        "workflows",
        "configure",
        "port",
        "port",
        "prepare",
        "process",
    ]
