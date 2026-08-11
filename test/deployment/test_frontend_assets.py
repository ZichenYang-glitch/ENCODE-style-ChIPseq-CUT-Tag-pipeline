"""Precompiled frontend identity, verification, and same-origin serving tests."""

from __future__ import annotations

import asyncio
import hashlib
import json
from pathlib import Path
import runpy

import httpx
import pytest

from encode_pipeline.deployment.errors import DeploymentError
from encode_pipeline.deployment.frontend import (
    API_CONTRACT_ENVIRONMENT,
    ManifestStaticApplication,
    create_app as create_deployed_app,
)
from encode_pipeline.frontend_assets import (
    AssetRecord,
    FrontendAssetManifest,
    VerifiedFrontendAssets,
    build_frontend_assets,
    canonical_json_sha256,
    load_packaged_frontend_assets,
    parse_manifest_bytes,
    verify_asset_directory,
    verify_frontend_api_contract,
)


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
PACKAGING_SCRIPT = REPOSITORY_ROOT / "scripts" / "package_frontend_assets.py"


def _write_build(root: Path) -> Path:
    root.mkdir(parents=True)
    (root / "assets").mkdir()
    (root / "index.html").write_bytes(
        b'<!doctype html><div id="root"></div><script src="/assets/app-a1.js"></script>\n'
    )
    (root / "favicon.svg").write_bytes(b"<svg></svg>\n")
    (root / "assets" / "app-a1.js").write_bytes(b"window.HELIXWEAVE = true;\n")
    (root / "assets" / "app-a1.css").write_bytes(b"body { color: black; }\n")
    return root


def _verified(root: Path, *, schema: object | None = None) -> VerifiedFrontendAssets:
    return build_frontend_assets(
        _write_build(root),
        frontend_version="0.3.0",
        package_lock_bytes=b'{"lockfileVersion":3}\n',
        openapi_bytes=json.dumps(
            schema
            if schema is not None
            else {"openapi": "3.1.0", "info": {"version": "0.3.0"}},
            indent=2,
        ).encode("utf-8"),
    )


def _request(app, method: str, path: str, **kwargs) -> httpx.Response:
    async def run() -> httpx.Response:
        transport = httpx.ASGITransport(app=app)
        async with httpx.AsyncClient(
            transport=transport,
            base_url="http://testserver",
            follow_redirects=False,
        ) as client:
            response = await client.request(method, path, **kwargs)
            await response.aread()
            return response

    return asyncio.run(run())


def test_manifest_is_canonical_complete_and_content_addressed(tmp_path: Path) -> None:
    first = _verified(tmp_path / "first")
    second = _verified(tmp_path / "second")

    assert first.manifest == second.manifest
    manifest_bytes = first.manifest.to_bytes()
    assert manifest_bytes.endswith(b"\n")
    assert b"\n" not in manifest_bytes[:-1]
    assert parse_manifest_bytes(manifest_bytes) == first.manifest
    assert [item.path for item in first.manifest.files] == [
        "assets/app-a1.css",
        "assets/app-a1.js",
        "favicon.svg",
        "index.html",
    ]
    assert set(first.content) == {item.path for item in first.manifest.files}
    assert first.manifest.entrypoint == "index.html"
    assert (
        first.manifest.package_lock_sha256
        == hashlib.sha256(b'{"lockfileVersion":3}\n').hexdigest()
    )

    changed_root = _write_build(tmp_path / "changed")
    (changed_root / "assets" / "app-a1.js").write_bytes(b"window.HELIXWEAVE = false;\n")
    changed = build_frontend_assets(
        changed_root,
        frontend_version="0.3.0",
        package_lock_bytes=b'{"lockfileVersion":3}\n',
        openapi_bytes=b'{"openapi":"3.1.0","info":{"version":"0.3.0"}}',
    )
    assert changed.manifest.identity != first.manifest.identity


def test_openapi_identity_is_semantic_and_exact_compatibility_is_fail_closed(
    tmp_path: Path,
) -> None:
    compact = b'{"info":{"version":"0.3.0"},"openapi":"3.1.0"}'
    pretty = b'{\n  "openapi": "3.1.0",\n  "info": {"version": "0.3.0"}\n}\n'
    assert canonical_json_sha256(compact) == canonical_json_sha256(pretty)

    assets = _verified(
        tmp_path / "dist",
        schema={"openapi": "3.1.0", "info": {"version": "0.3.0"}},
    )
    verify_frontend_api_contract(
        assets.manifest,
        {"info": {"version": "0.3.0"}, "openapi": "3.1.0"},
    )
    with pytest.raises(DeploymentError) as caught:
        verify_frontend_api_contract(
            assets.manifest,
            {"openapi": "3.1.0", "info": {"version": "0.3.1"}},
        )
    assert caught.value.issue.code == "FRONTEND_API_INCOMPATIBLE"
    assert "0.3.1" not in caught.value.issue.message


@pytest.mark.parametrize(
    "mutation",
    ("same-length", "missing", "extra"),
)
def test_asset_directory_verification_rejects_every_closure_drift(
    tmp_path: Path,
    mutation: str,
) -> None:
    assets = _verified(tmp_path / "dist")
    if mutation == "same-length":
        target = tmp_path / "dist" / "assets" / "app-a1.js"
        target.write_bytes(b"X" * len(target.read_bytes()))
    elif mutation == "missing":
        (tmp_path / "dist" / "favicon.svg").unlink()
    else:
        (tmp_path / "dist" / "unexpected.txt").write_bytes(b"unexpected")

    with pytest.raises(DeploymentError) as caught:
        verify_asset_directory(assets.manifest, tmp_path / "dist")
    assert caught.value.issue.code == "FRONTEND_ASSET_INTEGRITY_FAILED"
    assert str(tmp_path) not in caught.value.issue.message


def test_asset_scanner_rejects_symlinks_and_unsafe_manifest_paths(
    tmp_path: Path,
) -> None:
    build = _write_build(tmp_path / "dist")
    private_target = tmp_path / "private-secret.js"
    private_target.write_bytes(b"secret")
    (build / "assets" / "linked.js").symlink_to(private_target)

    with pytest.raises(DeploymentError) as caught:
        build_frontend_assets(
            build,
            frontend_version="0.3.0",
            package_lock_bytes=b"lock",
            openapi_bytes=b"{}",
        )
    assert caught.value.issue.code == "FRONTEND_ASSET_PACKAGE_INVALID"
    assert str(private_target) not in caught.value.issue.message

    unsafe = AssetRecord(
        path="../private.js",
        size_bytes=1,
        sha256=hashlib.sha256(b"x").hexdigest(),
    )
    with pytest.raises(DeploymentError) as unsafe_error:
        FrontendAssetManifest.create(
            frontend_version="0.3.0",
            package_lock_sha256=hashlib.sha256(b"lock").hexdigest(),
            api_contract_sha256=hashlib.sha256(b"api").hexdigest(),
            entrypoint="index.html",
            files=(unsafe,),
        )
    assert unsafe_error.value.issue.code == "FRONTEND_ASSET_MANIFEST_INVALID"


def test_manifest_parser_rejects_noncanonical_duplicate_and_unknown_fields(
    tmp_path: Path,
) -> None:
    manifest = _verified(tmp_path / "dist").manifest
    document = manifest.to_dict()
    noncanonical = json.dumps(document, indent=2).encode("utf-8")
    with pytest.raises(DeploymentError, match="FRONTEND_ASSET_MANIFEST_INVALID"):
        parse_manifest_bytes(noncanonical)

    duplicate = b'{"schema_version":"a","schema_version":"b"}\n'
    with pytest.raises(DeploymentError, match="FRONTEND_ASSET_MANIFEST_INVALID"):
        parse_manifest_bytes(duplicate)

    document["unexpected"] = True
    from encode_pipeline.deployment.canonical import canonical_json_bytes

    with pytest.raises(DeploymentError, match="FRONTEND_ASSET_MANIFEST_INVALID"):
        parse_manifest_bytes(canonical_json_bytes(document))

    entrypoint_record = next(
        item for item in manifest.files if item.path == manifest.entrypoint
    )
    duplicate_path = AssetRecord(
        path="index.html",
        size_bytes=entrypoint_record.size_bytes + 1,
        sha256=hashlib.sha256(b"duplicate").hexdigest(),
    )
    with pytest.raises(DeploymentError, match="FRONTEND_ASSET_MANIFEST_INVALID"):
        FrontendAssetManifest.create(
            frontend_version="0.3.0",
            package_lock_sha256=hashlib.sha256(b"lock").hexdigest(),
            api_contract_sha256=hashlib.sha256(b"api").hexdigest(),
            entrypoint="index.html",
            files=(*manifest.files, duplicate_path),
        )

    with pytest.raises(DeploymentError, match="FRONTEND_ASSET_MANIFEST_INVALID"):
        build_frontend_assets(
            tmp_path / "dist",
            frontend_version="0.3.0",
            package_lock_bytes=b"",
            openapi_bytes=b"{}",
        )


def test_packaging_script_replaces_only_generated_payload_and_is_offline(
    tmp_path: Path,
) -> None:
    module = runpy.run_path(str(PACKAGING_SCRIPT))
    source = PACKAGING_SCRIPT.read_text(encoding="utf-8")
    assert "subprocess" not in source
    assert "urllib" not in source
    assert "socket" not in source

    dist = _write_build(tmp_path / "dist")
    package_json = tmp_path / "package.json"
    package_json.write_text(
        '{"name":"helixweave-frontend","version":"0.3.0"}\n',
        encoding="utf-8",
    )
    package_lock = tmp_path / "package-lock.json"
    package_lock.write_text('{"lockfileVersion":3}\n', encoding="utf-8")
    openapi = tmp_path / "openapi.json"
    openapi.write_text('{"openapi":"3.1.0"}\n', encoding="utf-8")
    output = tmp_path / "package" / "frontend_assets"
    output.mkdir(parents=True)
    sentinel = output / "__init__.py"
    sentinel.write_text("# retained\n", encoding="utf-8")
    (output / "static").mkdir()
    (output / "static" / "stale.js").write_bytes(b"stale")

    identity = module["package_frontend_assets"](
        dist_root=dist,
        package_json=package_json,
        package_lock=package_lock,
        openapi=openapi,
        output_root=output,
    )

    assert identity.startswith("sha256-")
    assert sentinel.read_text(encoding="utf-8") == "# retained\n"
    assert not (output / "static" / "stale.js").exists()
    parsed = parse_manifest_bytes((output / "asset-manifest.json").read_bytes())
    verified = verify_asset_directory(parsed, output / "static")
    assert verified.manifest.identity == identity


def test_manifest_static_application_serves_only_verified_bytes_and_spa_routes(
    tmp_path: Path,
) -> None:
    assets = _verified(tmp_path / "dist")
    app = ManifestStaticApplication(assets)

    index = _request(app, "GET", "/")
    assert index.status_code == 200
    assert index.content == assets.content["index.html"]
    assert index.headers["cache-control"] == "no-cache"
    assert index.headers["x-content-type-options"] == "nosniff"
    assert index.headers["x-frame-options"] == "DENY"
    assert index.headers["referrer-policy"] == "no-referrer"
    assert "content-security-policy" not in index.headers

    script = _request(app, "GET", "/assets/app-a1.js")
    assert script.status_code == 200
    assert script.content == assets.content["assets/app-a1.js"]
    assert script.headers["cache-control"].endswith("immutable")
    assert script.headers["etag"].startswith('"')

    head = _request(app, "HEAD", "/assets/app-a1.js")
    assert head.status_code == 200
    assert head.content == b""
    assert int(head.headers["content-length"]) == len(script.content)
    unchanged = _request(
        app,
        "GET",
        "/assets/app-a1.js",
        headers={"if-none-match": script.headers["etag"]},
    )
    assert unchanged.status_code == 304
    assert unchanged.content == b""
    assert unchanged.headers["x-frame-options"] == "DENY"

    deep_link = _request(
        app,
        "GET",
        "/runs/run-123",
        headers={"accept": "text/html"},
    )
    assert deep_link.status_code == 200
    assert deep_link.content == index.content
    assert (
        _request(app, "GET", "/settings", headers={"accept": "text/html"}).status_code
        == 404
    )
    assert _request(app, "GET", "/assets/missing.js").status_code == 404
    assert (
        _request(
            app, "GET", "/runs/run-123", headers={"accept": "application/json"}
        ).status_code
        == 404
    )
    assert _request(app, "POST", "/").status_code == 405

    # Serving remains bound to the verified in-memory bytes, not mutable paths.
    (tmp_path / "dist" / "index.html").write_bytes(b"tampered")
    assert _request(app, "GET", "/").content == index.content


def test_deployed_factory_preserves_api_routes_and_serves_frontend_same_origin(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from fastapi import FastAPI
    import encode_pipeline.api.main as api_main

    class Closeable:
        def __init__(self) -> None:
            self.closed = False

        def close(self) -> None:
            self.closed = True

    captured: dict[str, object] = {}

    def create_test_api(**kwargs) -> FastAPI:
        captured.update(kwargs)
        api = FastAPI(title="HelixWeave API", version="0.3.0")

        @api.get("/api/v1/ping")
        async def ping() -> dict[str, bool]:
            return {"ok": True}

        api.state.run_queue = Closeable()
        api.state.persistence = Closeable()
        return api

    schema = create_test_api().openapi()
    assets = _verified(tmp_path / "dist", schema=schema)
    monkeypatch.setattr(api_main, "create_app", create_test_api)

    app = create_deployed_app(
        database_url=f"sqlite:///{tmp_path / 'runtime.db'}",
        workspace_root=tmp_path / "workspaces",
        project_root=REPOSITORY_ROOT,
        verified_assets=assets,
        api_contract_sha256=assets.manifest.api_contract_sha256,
    )
    try:
        ping = _request(app, "GET", "/api/v1/ping")
        assert ping.status_code == 200
        assert ping.json() == {"ok": True}
        assert _request(app, "GET", "/openapi.json").json() == schema
        assert _request(app, "GET", "/").content == assets.content["index.html"]
        assert (
            _request(
                app,
                "GET",
                "/workflows/bulk-rnaseq",
                headers={"accept": "text/html"},
            ).status_code
            == 200
        )
        assert app.state.frontend_asset_identity == assets.manifest.identity
        assert (
            app.state.frontend_api_contract_sha256
            == assets.manifest.api_contract_sha256
        )
        from encode_pipeline.persistence import open_existing_run_persistence
        from encode_pipeline.platform.registry import WorkflowRegistry

        assert captured["persistence_opener"] is open_existing_run_persistence
        assert isinstance(captured["registry"], WorkflowRegistry)
        assert captured["project_root"] == REPOSITORY_ROOT
    finally:
        app.state.run_queue.close()
        app.state.persistence.close()


def test_repository_package_assets_are_verified() -> None:
    assets = load_packaged_frontend_assets()

    assert assets.manifest.frontend_version == "0.3.0"
    assert assets.content["index.html"].startswith(b"<!doctype html>")
    assert any(path.startswith("assets/") for path in assets.content)
    assert assets.manifest.api_contract_sha256 == canonical_json_sha256(
        (REPOSITORY_ROOT / "frontend" / "openapi.json").read_bytes()
    )
    assert (
        assets.manifest.package_lock_sha256
        == hashlib.sha256(
            (REPOSITORY_ROOT / "frontend" / "package-lock.json").read_bytes()
        ).hexdigest()
    )


def test_deployed_factory_rejects_old_schema_without_running_migrations(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import sqlite3

    from encode_pipeline.persistence import runtime as persistence_runtime
    from encode_pipeline.persistence.runtime import (
        DatabaseSchemaNotReadyError,
        open_run_persistence,
    )

    database = tmp_path / "platform.db"
    prepared = open_run_persistence(f"sqlite:///{database}")
    prepared.close()
    with sqlite3.connect(database) as connection:
        connection.execute("UPDATE alembic_version SET version_num = ?", ("old",))
    monkeypatch.setattr(
        persistence_runtime,
        "upgrade_database",
        lambda _database_url: (_ for _ in ()).throw(
            AssertionError("supported service must not run Alembic")
        ),
    )
    assets = load_packaged_frontend_assets()

    with pytest.raises(DatabaseSchemaNotReadyError) as raised:
        create_deployed_app(
            database_url=f"sqlite:///{database}",
            workspace_root=tmp_path / "workspaces",
            project_root=REPOSITORY_ROOT,
            verified_assets=assets,
            api_contract_sha256=assets.manifest.api_contract_sha256,
        )

    assert raised.value.reason_code == "DATABASE_SCHEMA_NOT_CURRENT"
    with sqlite3.connect(database) as connection:
        assert connection.execute(
            "SELECT version_num FROM alembic_version"
        ).fetchone() == ("old",)


@pytest.mark.parametrize("configured", (None, "0" * 64))
def test_deployed_factory_rejects_api_contract_before_api_factory_side_effects(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    configured: str | None,
) -> None:
    import encode_pipeline.api.main as api_main

    assets = _verified(tmp_path / "dist")
    invoked = False

    def create_test_api(**_kwargs):
        nonlocal invoked
        invoked = True
        raise AssertionError("API factory must not be invoked")

    monkeypatch.setattr(api_main, "create_app", create_test_api)
    monkeypatch.delenv(API_CONTRACT_ENVIRONMENT, raising=False)

    with pytest.raises(DeploymentError) as caught:
        create_deployed_app(
            verified_assets=assets,
            api_contract_sha256=configured,
        )

    assert caught.value.issue.code == "FRONTEND_API_INCOMPATIBLE"
    assert invoked is False
