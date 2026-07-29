#!/usr/bin/env python3
"""Fail-closed HelixWeave source provenance checks.

This module is deliberately standard-library-only and never imports product
code.  The clean checkout bootstrap loads it under ``python -I -S`` before
making source or third-party package paths importable.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from email.parser import Parser
from importlib import util
import json
import os
from pathlib import Path
import re
import sys
import sysconfig
import tomllib
from urllib.parse import urlparse
from urllib.request import url2pathname


DISTRIBUTION_NAME = "helixweave"
DISTRIBUTION_VERSION = "0.3.0"
IMPORT_NAMESPACE = "encode_pipeline"
_KNOWN_DISTRIBUTION_NAMES = (
    "helixweave",
    "encode-pipeline",
)
_METADATA_SUFFIXES = (".dist-info", ".egg-info")


@dataclass(frozen=True)
class SourceProvenanceError(Exception):
    """A public-safe, stable provenance failure."""

    reason_code: str
    guidance: str

    def __str__(self) -> str:
        return f"source provenance check failed [{self.reason_code}]: {self.guidance}"


@dataclass(frozen=True)
class DistributionClaim:
    """One physical distribution metadata claimant."""

    metadata_path: Path
    site_root: Path
    name: str
    version: str
    top_levels: frozenset[str]
    record_entries: tuple[str, ...] | None
    direct_url: dict[str, object] | None


@dataclass(frozen=True)
class CheckoutLayout:
    """Resolved checkout identity."""

    repository_root: Path
    source_root: Path
    package_root: Path
    initializer: Path
    version: str


def _fail(reason_code: str, guidance: str) -> None:
    raise SourceProvenanceError(reason_code, guidance)


def _normalized_distribution_name(value: str) -> str:
    return re.sub(r"[-_.]+", "-", value).lower()


def _resolved_existing_directory(
    path: Path,
    *,
    reason_code: str,
    guidance: str,
) -> Path:
    try:
        resolved = path.resolve(strict=True)
    except (OSError, RuntimeError):
        _fail(reason_code, guidance)
    if not resolved.is_dir():
        _fail(reason_code, guidance)
    return resolved


def _resolved_existing_file(
    path: Path,
    *,
    reason_code: str,
    guidance: str,
) -> Path:
    try:
        resolved = path.resolve(strict=True)
    except (OSError, RuntimeError):
        _fail(reason_code, guidance)
    if not resolved.is_file():
        _fail(reason_code, guidance)
    return resolved


def _is_within(path: Path, root: Path) -> bool:
    return path == root or path.is_relative_to(root)


def _checkout_layout(repository_root: Path) -> CheckoutLayout:
    root = _resolved_existing_directory(
        repository_root,
        reason_code="repository_root_invalid",
        guidance="provide the expected existing repository root explicitly",
    )
    source_root = _resolved_existing_directory(
        root / "src",
        reason_code="source_root_invalid",
        guidance="the expected repository must contain its canonical src directory",
    )
    if source_root.parent != root:
        _fail(
            "source_root_invalid",
            "the expected repository src directory must not escape its root",
        )
    package_root = _resolved_existing_directory(
        source_root / IMPORT_NAMESPACE,
        reason_code="source_root_invalid",
        guidance="the expected src directory must contain encode_pipeline",
    )
    if not _is_within(package_root, source_root):
        _fail(
            "source_root_invalid",
            "the expected package directory must not escape its src directory",
        )
    initializer = _resolved_existing_file(
        package_root / "__init__.py",
        reason_code="source_root_invalid",
        guidance="the expected package must have its canonical initializer",
    )
    try:
        project = tomllib.loads((root / "pyproject.toml").read_text(encoding="utf-8"))
        project_metadata = project["project"]
        project_name = project_metadata["name"]
        version = project_metadata["version"]
    except (KeyError, OSError, TypeError, UnicodeError, tomllib.TOMLDecodeError):
        _fail(
            "repository_metadata_invalid",
            "repair the checkout project metadata before retrying",
        )
    if (
        not isinstance(project_name, str)
        or _normalized_distribution_name(project_name) != DISTRIBUTION_NAME
        or not isinstance(version, str)
        or not version
    ):
        _fail(
            "repository_metadata_invalid",
            "repair the checkout project identity before retrying",
        )
    return CheckoutLayout(root, source_root, package_root, initializer, version)


def _venv_root_from_executable() -> Path | None:
    executable = Path(os.path.abspath(sys.executable))
    candidates = (
        executable.parent.parent,
        executable.parent,
    )
    for candidate in candidates:
        if (candidate / "pyvenv.cfg").is_file():
            return candidate.resolve()
    return None


def environment_site_roots() -> tuple[Path, ...]:
    """Return the selected interpreter environment's unprocessed site roots."""

    venv_root = _venv_root_from_executable()
    candidates: list[Path] = []
    if venv_root is not None:
        if os.name == "nt":
            candidates.append(venv_root / "Lib" / "site-packages")
        else:
            candidates.append(
                venv_root
                / "lib"
                / f"python{sys.version_info.major}.{sys.version_info.minor}"
                / "site-packages"
            )
    else:
        for key in ("purelib", "platlib"):
            raw = sysconfig.get_path(key)
            if raw:
                candidates.append(Path(raw))

    roots: list[Path] = []
    for candidate in candidates:
        try:
            resolved = candidate.resolve(strict=True)
        except (OSError, RuntimeError):
            continue
        if resolved.is_dir() and resolved not in roots:
            roots.append(resolved)
    if not roots:
        _fail(
            "environment_site_root_invalid",
            "use an interpreter with one discoverable environment site-packages root",
        )
    return tuple(roots)


def _active_site_roots() -> tuple[Path, ...]:
    roots = list(environment_site_roots())
    for entry in sys.path:
        if not isinstance(entry, str) or not entry:
            continue
        candidate = Path(entry)
        if candidate.name not in {"site-packages", "dist-packages"}:
            continue
        try:
            resolved = candidate.resolve(strict=True)
        except (OSError, RuntimeError):
            continue
        if resolved.is_dir() and resolved not in roots:
            roots.append(resolved)
    return tuple(roots)


def _site_roots(explicit: tuple[Path, ...] | None) -> tuple[Path, ...]:
    if explicit is None:
        return _active_site_roots()
    roots: list[Path] = []
    for raw in explicit:
        root = _resolved_existing_directory(
            raw,
            reason_code="environment_site_root_invalid",
            guidance="use the selected interpreter environment site-packages",
        )
        if root in roots:
            _fail(
                "namespace_mapping_conflict",
                "remove duplicate or aliased environment package roots",
            )
        roots.append(root)
    if not roots:
        _fail(
            "environment_site_root_invalid",
            "use the selected interpreter environment site-packages",
        )
    return tuple(roots)


def _metadata_path_looks_known(path: Path) -> bool:
    normalized = _normalized_distribution_name(path.name)
    return any(
        normalized.startswith(f"{known}-") for known in _KNOWN_DISTRIBUTION_NAMES
    )


def _read_optional_text(path: Path, *, known: bool) -> str | None:
    if not path.exists():
        return None
    try:
        if not path.is_file():
            raise OSError
        return path.read_text(encoding="utf-8")
    except (OSError, UnicodeError):
        if known:
            _fail(
                "distribution_metadata_invalid",
                "repair the HelixWeave distribution metadata before retrying",
            )
        return None


def _parse_record(raw: str | None, *, known: bool) -> tuple[str, ...] | None:
    if raw is None:
        return None
    try:
        rows = tuple(csv.reader(raw.splitlines()))
    except (csv.Error, TypeError, UnicodeError):
        if known:
            _fail(
                "distribution_metadata_invalid",
                "repair the HelixWeave distribution inventory before retrying",
            )
        return None
    entries: list[str] = []
    for row in rows:
        if not row or not row[0]:
            if known:
                _fail(
                    "distribution_metadata_invalid",
                    "repair the HelixWeave distribution inventory before retrying",
                )
            return None
        entry = row[0].replace("\\", "/")
        path = Path(entry)
        if path.is_absolute():
            if known:
                _fail(
                    "distribution_metadata_invalid",
                    "repair the HelixWeave distribution inventory before retrying",
                )
            return None
        entries.append(entry)
    return tuple(entries)


def _parse_direct_url(raw: str | None, *, known: bool) -> dict[str, object] | None:
    if raw is None:
        return None
    try:
        value = json.loads(raw)
    except (TypeError, ValueError):
        value = None
    if (
        not isinstance(value, dict)
        or not isinstance(value.get("url"), str)
        or not value["url"]
    ):
        if known:
            _fail(
                "distribution_metadata_invalid",
                "repair the HelixWeave direct URL metadata before retrying",
            )
        return None
    return value


def _metadata_claim(path: Path, site_root: Path) -> DistributionClaim | None:
    path_hint = _metadata_path_looks_known(path)
    metadata_name = "METADATA" if path.name.endswith(".dist-info") else "PKG-INFO"
    raw_metadata = _read_optional_text(path / metadata_name, known=path_hint)
    if raw_metadata is None:
        if path_hint:
            _fail(
                "distribution_metadata_invalid",
                "repair the HelixWeave distribution metadata before retrying",
            )
        return None
    try:
        parsed = Parser().parsestr(raw_metadata)
        name = parsed["Name"]
        version = parsed["Version"]
    except (TypeError, UnicodeError):
        name = version = None
    normalized_name = (
        _normalized_distribution_name(name) if isinstance(name, str) else None
    )
    known = path_hint or normalized_name in _KNOWN_DISTRIBUTION_NAMES

    raw_top_level = _read_optional_text(path / "top_level.txt", known=known)
    top_levels = frozenset(
        line.strip() for line in (raw_top_level or "").splitlines() if line.strip()
    )
    raw_record = _read_optional_text(path / "RECORD", known=known)
    record_entries = _parse_record(raw_record, known=known)
    record_claims = record_entries is not None and any(
        entry.split("/", 1)[0] == IMPORT_NAMESPACE for entry in record_entries
    )
    claims_namespace = IMPORT_NAMESPACE in top_levels or record_claims
    if not known and not claims_namespace:
        return None
    if (
        not isinstance(name, str)
        or not name.strip()
        or not isinstance(version, str)
        or not version.strip()
    ):
        _fail(
            "distribution_metadata_invalid",
            "repair the HelixWeave distribution identity metadata before retrying",
        )
    direct_url = _parse_direct_url(
        _read_optional_text(path / "direct_url.json", known=known),
        known=known,
    )
    return DistributionClaim(
        metadata_path=path,
        site_root=site_root,
        name=name,
        version=version,
        top_levels=top_levels,
        record_entries=record_entries,
        direct_url=direct_url,
    )


def _physical_claimants(site_roots: tuple[Path, ...]) -> tuple[DistributionClaim, ...]:
    claims: list[DistributionClaim] = []
    for site_root in site_roots:
        try:
            children = tuple(sorted(site_root.iterdir(), key=lambda path: path.name))
        except OSError:
            _fail(
                "distribution_metadata_invalid",
                "repair the selected environment metadata before retrying",
            )
        for child in children:
            if not child.is_dir() or not child.name.endswith(_METADATA_SUFFIXES):
                continue
            claim = _metadata_claim(child, site_root)
            if claim is not None:
                claims.append(claim)
    if not claims:
        _fail(
            "distribution_missing",
            "install HelixWeave in the selected source mode before retrying",
        )
    if len(claims) != 1:
        _fail(
            "distribution_claimant_conflict",
            "remove every additional encode_pipeline distribution claimant",
        )
    return tuple(claims)


def _verify_claim_identity(claim: DistributionClaim, expected_version: str) -> None:
    if _normalized_distribution_name(claim.name) != DISTRIBUTION_NAME:
        _fail(
            "distribution_identity_mismatch",
            "remove conflicting namespace owners from the selected environment",
        )
    if claim.version != expected_version:
        _fail(
            "distribution_version_mismatch",
            "install the exact HelixWeave version declared by the selected source",
        )


def _editable_root(direct_url: dict[str, object] | None) -> Path | None:
    if direct_url is None:
        return None
    directory_info = direct_url.get("dir_info")
    if (
        not isinstance(directory_info, dict)
        or directory_info.get("editable") is not True
    ):
        return None
    raw_url = direct_url.get("url")
    if not isinstance(raw_url, str):
        _fail(
            "distribution_metadata_invalid",
            "repair the HelixWeave editable source metadata before retrying",
        )
    parsed = urlparse(raw_url)
    if (
        parsed.scheme != "file"
        or parsed.netloc not in {"", "localhost"}
        or parsed.query
        or parsed.fragment
    ):
        _fail(
            "distribution_metadata_invalid",
            "use canonical local editable metadata for checkout mode",
        )
    try:
        return Path(url2pathname(parsed.path)).resolve(strict=True)
    except (OSError, RuntimeError, ValueError):
        _fail(
            "distribution_source_mismatch",
            "reinstall the expected checkout into the selected environment",
        )


def _record_owned_pth_names(claim: DistributionClaim) -> tuple[str, ...]:
    if claim.record_entries is None:
        return ()
    return tuple(
        entry
        for entry in claim.record_entries
        if "/" not in entry and entry.endswith(".pth")
    )


def _product_relevant_executable_pth(path: Path, line: str) -> bool:
    text = f"{path.name}\n{line}".lower()
    return any(
        token in text
        for token in (
            IMPORT_NAMESPACE,
            DISTRIBUTION_NAME,
            "encode-pipeline",
            "__editable__",
        )
    )


def _pth_mappings(
    site_roots: tuple[Path, ...],
) -> tuple[tuple[Path, Path], ...]:
    namespace_mappings: list[tuple[Path, Path]] = []
    for site_root in site_roots:
        for hook_name in ("sitecustomize.py", "usercustomize.py"):
            if (site_root / hook_name).exists():
                _fail(
                    "startup_hook_unsafe",
                    "remove environment startup customization before retrying",
                )
        try:
            finder_files = tuple(site_root.glob("__editable__*finder.py"))
            pth_files = tuple(sorted(site_root.glob("*.pth")))
        except OSError:
            _fail(
                "pth_mapping_unsafe",
                "repair the selected environment path metadata before retrying",
            )
        if finder_files:
            _fail(
                "pth_mapping_unsafe",
                "replace executable editable finders with one plain source mapping",
            )
        for pth_path in pth_files:
            try:
                raw = pth_path.read_text(encoding="utf-8")
            except (OSError, UnicodeError):
                _fail(
                    "pth_mapping_unsafe",
                    "repair unreadable environment path metadata before retrying",
                )
            if "\x00" in raw:
                _fail(
                    "pth_mapping_unsafe",
                    "repair unsafe environment path metadata before retrying",
                )
            for raw_line in raw.splitlines():
                line = raw_line.strip()
                if not line or line.startswith("#"):
                    continue
                if re.match(r"^import(?:[ \t]|$)", line):
                    if _product_relevant_executable_pth(pth_path, line):
                        _fail(
                            "pth_mapping_unsafe",
                            "replace executable product .pth hooks with a plain mapping",
                        )
                    continue
                candidate = Path(line)
                if not candidate.is_absolute():
                    candidate = site_root / candidate
                try:
                    mapping = candidate.resolve(strict=True)
                except (OSError, RuntimeError):
                    if "__editable__" in pth_path.name.lower():
                        _fail(
                            "editable_mapping_invalid",
                            "repair the recorded editable source mapping before retrying",
                        )
                    continue
                if not mapping.is_dir():
                    if "__editable__" in pth_path.name.lower():
                        _fail(
                            "editable_mapping_invalid",
                            "repair the recorded editable source mapping before retrying",
                        )
                    continue
                if (mapping / IMPORT_NAMESPACE).is_dir():
                    namespace_mappings.append((pth_path.resolve(), mapping))
    return tuple(namespace_mappings)


def _audit_checkout(
    layout: CheckoutLayout,
    site_roots: tuple[Path, ...],
) -> DistributionClaim:
    try:
        source_metadata = tuple(layout.source_root.glob("*.egg-info"))
    except OSError:
        _fail(
            "distribution_metadata_invalid",
            "repair source-tree build metadata before retrying",
        )
    for metadata_path in source_metadata:
        if not metadata_path.is_dir() or not _metadata_path_looks_known(metadata_path):
            continue
        raw_direct_url = _read_optional_text(
            metadata_path / "direct_url.json",
            known=False,
        )
        if raw_direct_url is None:
            continue
        direct_url = _parse_direct_url(raw_direct_url, known=True)
        if _editable_root(direct_url) != layout.repository_root:
            _fail(
                "distribution_source_mismatch",
                "repair source metadata that refers to another checkout",
            )

    (claim,) = _physical_claimants(site_roots)
    _verify_claim_identity(claim, layout.version)
    if IMPORT_NAMESPACE not in claim.top_levels or claim.record_entries is None:
        _fail(
            "namespace_ownership_unproven",
            "install metadata that explicitly inventories encode_pipeline ownership",
        )
    if _editable_root(claim.direct_url) != layout.repository_root:
        _fail(
            "distribution_source_mismatch",
            "reinstall the expected checkout into the selected environment",
        )
    owned_pth_names = _record_owned_pth_names(claim)
    if len(owned_pth_names) != 1:
        _fail(
            "editable_mapping_invalid",
            "record exactly one plain editable source mapping",
        )
    mappings = _pth_mappings(site_roots)
    if len(mappings) != 1:
        _fail(
            "namespace_mapping_conflict",
            "remove orphan, stale, duplicate, or aliased namespace mappings",
        )
    mapping_path, mapping_root = mappings[0]
    if (
        mapping_path.parent != claim.site_root
        or mapping_path.name != owned_pth_names[0]
        or mapping_root != layout.source_root
    ):
        _fail(
            "editable_mapping_invalid",
            "bind the sole recorded editable mapping to the expected checkout src",
        )
    return claim


def _audit_installed(
    site_roots: tuple[Path, ...],
) -> DistributionClaim:
    (claim,) = _physical_claimants(site_roots)
    _verify_claim_identity(claim, DISTRIBUTION_VERSION)
    if claim.direct_url is not None:
        if _editable_root(claim.direct_url) is not None:
            _fail(
                "installed_editable_forbidden",
                "install a clean wheel or sdist into the isolated environment",
            )
        if "dir_info" in claim.direct_url or not isinstance(
            claim.direct_url.get("archive_info"),
            dict,
        ):
            _fail(
                "installed_source_mismatch",
                "install a clean wheel or rebuilt sdist artifact",
            )
    if claim.record_entries is None or not any(
        entry == f"{IMPORT_NAMESPACE}/__init__.py" for entry in claim.record_entries
    ):
        _fail(
            "namespace_ownership_unproven",
            "install an artifact whose RECORD explicitly owns encode_pipeline",
        )
    if _pth_mappings(site_roots):
        _fail(
            "namespace_mapping_conflict",
            "remove source mappings from the installed artifact environment",
        )
    return claim


def _package_spec() -> object:
    try:
        spec = util.find_spec(IMPORT_NAMESPACE)
    except (AttributeError, ImportError, ModuleNotFoundError, ValueError):
        spec = None
    if spec is None:
        _fail(
            "module_not_found",
            "install the requested artifact or editable checkout before retrying",
        )
    return spec


def _resolved_search_locations(spec: object) -> tuple[Path, ...]:
    raw_locations = getattr(spec, "submodule_search_locations", None)
    if raw_locations is None:
        _fail(
            "module_search_location_mismatch",
            "ensure the package has one canonical search location",
        )
    try:
        return tuple(Path(location).resolve(strict=True) for location in raw_locations)
    except (OSError, RuntimeError):
        _fail(
            "module_search_location_mismatch",
            "ensure every package search location exists and is canonical",
        )


def _resolved_origin(spec: object) -> Path:
    raw_origin = getattr(spec, "origin", None)
    if not isinstance(raw_origin, str) or raw_origin in {"built-in", "frozen"}:
        _fail(
            "module_origin_mismatch",
            "ensure the package resolves from the selected source mode",
        )
    return _resolved_existing_file(
        Path(raw_origin),
        reason_code="module_origin_mismatch",
        guidance="ensure the package resolves from the selected source mode",
    )


def _verify_checkout_spec(layout: CheckoutLayout) -> None:
    spec = _package_spec()
    locations = _resolved_search_locations(spec)
    if len(locations) != 1 or locations[0] != layout.package_root:
        _fail(
            "module_search_location_mismatch",
            "run from one canonical package location in the expected checkout",
        )
    if _resolved_origin(spec) != layout.initializer:
        _fail(
            "module_origin_mismatch",
            "reinstall the expected checkout into the selected environment",
        )


def _verify_installed_spec(
    claim: DistributionClaim,
    site_roots: tuple[Path, ...],
) -> None:
    spec = _package_spec()
    locations = _resolved_search_locations(spec)
    origin = _resolved_origin(spec)
    if (
        len(locations) != 1
        or not any(_is_within(locations[0], root) for root in site_roots)
        or not any(_is_within(origin, root) for root in site_roots)
        or not _is_within(origin, claim.site_root)
    ):
        _fail(
            "installed_source_mismatch",
            "remove source checkouts and external editable paths before retrying",
        )
    try:
        relative_origin = origin.relative_to(claim.site_root).as_posix()
    except ValueError:
        _fail(
            "installed_source_mismatch",
            "install the artifact into the selected environment site-packages",
        )
    if claim.record_entries is None or relative_origin not in claim.record_entries:
        _fail(
            "namespace_ownership_unproven",
            "install an artifact whose RECORD owns the resolved package origin",
        )


def _prepend_audited_paths(source_root: Path | None, roots: tuple[Path, ...]) -> None:
    additions = [str(root) for root in roots]
    if source_root is not None:
        additions.insert(0, str(source_root))
    for addition in reversed(additions):
        while addition in sys.path:
            sys.path.remove(addition)
        sys.path.insert(0, addition)


def verify_checkout(
    repository_root: Path,
    *,
    site_roots: tuple[Path, ...] | None = None,
) -> None:
    """Prove source, metadata, and mapping identity for one checkout."""

    layout = _checkout_layout(repository_root)
    roots = _site_roots(site_roots)
    _audit_checkout(layout, roots)
    _verify_checkout_spec(layout)


def bootstrap_checkout(repository_root: Path) -> tuple[Path, ...]:
    """Audit, then activate, the selected checkout without processing site hooks."""

    layout = _checkout_layout(repository_root)
    roots = _site_roots(environment_site_roots())
    _audit_checkout(layout, roots)
    _prepend_audited_paths(layout.source_root, roots)
    _verify_checkout_spec(layout)
    return roots


def _installed_environment_is_isolated() -> bool:
    return _venv_root_from_executable() is not None or sys.prefix != sys.base_prefix


def verify_installed_artifact(
    *,
    site_roots: tuple[Path, ...] | None = None,
) -> None:
    """Prove one installed wheel or rebuilt sdist owns the namespace."""

    if not _installed_environment_is_isolated():
        _fail(
            "installed_environment_invalid",
            "select installed-artifact mode only in an isolated environment",
        )
    roots = _site_roots(site_roots)
    claim = _audit_installed(roots)
    _verify_installed_spec(claim, roots)


def bootstrap_installed_artifact() -> tuple[Path, ...]:
    """Audit, then activate, one installed artifact without site startup hooks."""

    if not _installed_environment_is_isolated():
        _fail(
            "installed_environment_invalid",
            "select installed-artifact mode only in an isolated environment",
        )
    roots = _site_roots(environment_site_roots())
    claim = _audit_installed(roots)
    _prepend_audited_paths(None, roots)
    _verify_installed_spec(claim, roots)
    return roots


def require_checkout(repository_root: Path) -> None:
    """Exit safely when a compatibility wrapper has invalid provenance."""

    try:
        verify_checkout(repository_root)
    except SourceProvenanceError as error:
        print(error, file=sys.stderr)
        raise SystemExit(2) from None


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Verify HelixWeave Python source provenance."
    )
    subparsers = parser.add_subparsers(dest="mode", required=True)
    checkout = subparsers.add_parser(
        "checkout",
        help="verify one explicitly requested source checkout",
    )
    checkout.add_argument("--repository-root", type=Path, required=True)
    subparsers.add_parser(
        "installed-artifact",
        help="verify a clean wheel or rebuilt sdist installation",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if args.mode == "checkout":
            verify_checkout(args.repository_root)
        else:
            verify_installed_artifact()
    except SourceProvenanceError as error:
        print(error, file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
