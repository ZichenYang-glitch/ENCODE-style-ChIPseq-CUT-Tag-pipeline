"""Verified, package-owned production frontend assets.

The frontend is built once on a controlled build host.  Runtime hosts consume
only the resulting immutable bytes and this canonical manifest; they never run
Node.js, npm, Vite, or source checkout tooling.
"""

from encode_pipeline.frontend_assets.manifest import (
    FRONTEND_ASSET_IDENTITY_SCHEME,
    FRONTEND_ASSET_MANIFEST_SCHEMA,
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

__all__ = [
    "FRONTEND_ASSET_IDENTITY_SCHEME",
    "FRONTEND_ASSET_MANIFEST_SCHEMA",
    "AssetRecord",
    "FrontendAssetManifest",
    "VerifiedFrontendAssets",
    "build_frontend_assets",
    "canonical_json_sha256",
    "load_packaged_frontend_assets",
    "parse_manifest_bytes",
    "verify_asset_directory",
    "verify_frontend_api_contract",
]
