#!/usr/bin/env python3
"""Fail-closed HelixWeave source provenance checks.

This module is deliberately standard-library-only and never imports product
code.  The clean checkout bootstrap loads it under ``python -I -S`` before
making source or third-party package paths importable.
"""

from __future__ import annotations

import argparse
import ast
import base64
import csv
from dataclasses import dataclass
from email.parser import Parser
import hashlib
from importlib import machinery
from importlib import util
import json
import os
from pathlib import Path
import re
import stat
import sys
import sysconfig
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
_PTH_INVENTORY_GUIDANCE = "remove unrecognized executable .pth startup hooks"


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


@dataclass(frozen=True)
class AuditedExecutablePth:
    """One exact non-product startup hook and its locked owner."""

    filename: str
    sha256: str
    distribution: str
    versions: frozenset[str]
    conda_builds: frozenset[tuple[str, str]]


_AUDITED_EXECUTABLE_PTH = {
    spec.filename: spec
    for spec in (
        AuditedExecutablePth(
            filename="a1_coverage.pth",
            sha256="ef2ed06d19867ec669c09a804060666a9cd5e383af0a9d11aa2de79b77d448e8",
            distribution="coverage",
            versions=frozenset({"7.15.1"}),
            conda_builds=frozenset({("7.15.1", "py312h8a5da7c_0")}),
        ),
        AuditedExecutablePth(
            filename="coloredlogs.pth",
            sha256="dda83a855986efa5cd87f0248b0199c0086eb0e8e7fece7d6741959c5ce39536",
            distribution="coloredlogs",
            versions=frozenset({"15.0.1"}),
            conda_builds=frozenset({("15.0.1", "pyhd8ed1ab_4")}),
        ),
        AuditedExecutablePth(
            filename="distutils-precedence.pth",
            sha256="2638ce9e2500e572a5e0de7faed6661eb569d1b696fcba07b0dd223da5f5d224",
            distribution="setuptools",
            versions=frozenset({"82.0.1", "83.0.0"}),
            conda_builds=frozenset(
                {
                    ("82.0.1", "pyh332efcf_0"),
                    ("83.0.0", "pyh332efcf_0"),
                }
            ),
        ),
        AuditedExecutablePth(
            filename="sphinxcontrib_jsmath-1.0.1-py3.9-nspkg.pth",
            sha256="a328acccc2310e241f406ef5fd39a60ce5ffdd678a9454be561d28af01a81e54",
            distribution="sphinxcontrib-jsmath",
            versions=frozenset({"1.0.1"}),
            conda_builds=frozenset({("1.0.1", "pyhd8ed1ab_1")}),
        ),
    )
}


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
    if initializer.parent != package_root:
        _fail(
            "source_root_invalid",
            "the expected package initializer must not escape its package directory",
        )
    project_file = _resolved_existing_file(
        root / "pyproject.toml",
        reason_code="repository_metadata_invalid",
        guidance="repair the checkout project metadata before retrying",
    )
    if project_file.parent != root:
        _fail(
            "repository_metadata_invalid",
            "the checkout project metadata must not escape its repository root",
        )
    try:
        project_name, version = _bounded_project_identity(
            project_file.read_text(encoding="utf-8")
        )
    except (OSError, TypeError, UnicodeError, ValueError):
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


def _bounded_project_identity(raw: str) -> tuple[str, str]:
    """Read the two literal identity fields without requiring TOML in Python 3.10."""

    in_project = False
    values: dict[str, str] = {}
    for raw_line in raw.splitlines():
        line = raw_line.strip()
        if line.startswith("["):
            if line == "[project]":
                if in_project:
                    raise ValueError
                in_project = True
                continue
            if in_project:
                break
            continue
        if not in_project or not line or line.startswith("#"):
            continue
        match = re.fullmatch(r"(name|version)\s*=\s*(.+)", line)
        if match is None:
            continue
        key, encoded = match.groups()
        if key in values:
            raise ValueError
        try:
            value = ast.literal_eval(encoded)
        except (SyntaxError, ValueError):
            raise ValueError from None
        if not isinstance(value, str) or not value:
            raise ValueError
        values[key] = value
    if set(values) != {"name", "version"}:
        raise ValueError
    return values["name"], values["version"]


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


def _read_optional_text(
    path: Path,
    *,
    known: bool,
    container: Path | None = None,
) -> str | None:
    try:
        is_symlink = path.is_symlink()
        exists = path.exists()
    except OSError:
        is_symlink = exists = False
    if is_symlink:
        _fail(
            "distribution_metadata_invalid",
            "repair metadata whose file identity escapes the selected environment",
        )
    if not exists:
        return None
    try:
        resolved = path.resolve(strict=True)
        if not resolved.is_file() or (
            container is not None and (resolved.parent != container or resolved != path)
        ):
            raise OSError
        return resolved.read_text(encoding="utf-8")
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


def _metadata_identity(
    raw: str,
    *,
    known: bool,
) -> tuple[str, str] | None:
    try:
        parsed = Parser().parsestr(raw)
        names = parsed.get_all("Name", [])
        versions = parsed.get_all("Version", [])
        valid = (
            not parsed.defects
            and len(names) == 1
            and isinstance(names[0], str)
            and bool(names[0].strip())
            and len(versions) == 1
            and isinstance(versions[0], str)
            and bool(versions[0].strip())
        )
    except (TypeError, UnicodeError):
        valid = False
    if not valid:
        if known:
            _fail(
                "distribution_metadata_invalid",
                "repair the HelixWeave distribution identity metadata before retrying",
            )
        return None
    return names[0].strip(), versions[0].strip()


def _mapping_inventory_claims_namespace(
    record_entries: tuple[str, ...] | None,
    site_root: Path,
    path_mappings: tuple[tuple[Path, Path], ...],
) -> bool:
    if record_entries is None:
        return False
    recorded_pth_names = {
        entry for entry in record_entries if "/" not in entry and entry.endswith(".pth")
    }
    for pth_path, mapping_root in path_mappings:
        if pth_path.parent != site_root or pth_path.name not in recorded_pth_names:
            continue
        try:
            package_root = (mapping_root / IMPORT_NAMESPACE).resolve(strict=True)
        except (OSError, RuntimeError):
            continue
        if package_root.is_dir() and _is_within(package_root, mapping_root):
            return True
    return False


def _metadata_claim(
    path: Path,
    site_root: Path,
    path_mappings: tuple[tuple[Path, Path], ...],
) -> DistributionClaim | None:
    path_hint = _metadata_path_looks_known(path)
    raw_top_level = _read_optional_text(
        path / "top_level.txt",
        known=path_hint,
        container=path,
    )
    top_levels = frozenset(
        line.strip() for line in (raw_top_level or "").splitlines() if line.strip()
    )
    raw_record = _read_optional_text(
        path / "RECORD",
        known=path_hint,
        container=path,
    )
    record_entries = _parse_record(raw_record, known=path_hint)
    record_claims = record_entries is not None and any(
        entry.split("/", 1)[0] == IMPORT_NAMESPACE for entry in record_entries
    )
    claims_namespace = (
        IMPORT_NAMESPACE in top_levels
        or record_claims
        or _mapping_inventory_claims_namespace(
            record_entries,
            site_root,
            path_mappings,
        )
    )

    metadata_name = "METADATA" if path.name.endswith(".dist-info") else "PKG-INFO"
    raw_metadata = _read_optional_text(
        path / metadata_name,
        known=path_hint or claims_namespace,
        container=path,
    )
    if raw_metadata is None:
        if path_hint or claims_namespace:
            _fail(
                "distribution_metadata_invalid",
                "repair the HelixWeave distribution metadata before retrying",
            )
        return None
    identity = _metadata_identity(
        raw_metadata,
        known=path_hint or claims_namespace,
    )
    if identity is None:
        return None
    name, version = identity
    known = (
        path_hint
        or claims_namespace
        or _normalized_distribution_name(name) in _KNOWN_DISTRIBUTION_NAMES
    )
    if not known:
        return None

    raw_top_level = _read_optional_text(
        path / "top_level.txt",
        known=True,
        container=path,
    )
    top_levels = frozenset(
        line.strip() for line in (raw_top_level or "").splitlines() if line.strip()
    )
    raw_record = _read_optional_text(
        path / "RECORD",
        known=True,
        container=path,
    )
    record_entries = _parse_record(raw_record, known=True)
    direct_url = _parse_direct_url(
        _read_optional_text(
            path / "direct_url.json",
            known=True,
            container=path,
        ),
        known=True,
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


def _physical_claimants(
    site_roots: tuple[Path, ...],
    path_mappings: tuple[tuple[Path, Path], ...],
) -> tuple[DistributionClaim, ...]:
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
            if not child.name.endswith(_METADATA_SUFFIXES):
                continue
            try:
                resolved_child = child.resolve(strict=True)
            except (OSError, RuntimeError):
                _fail(
                    "distribution_metadata_invalid",
                    "repair metadata whose identity escapes the selected environment",
                )
            if (
                child.is_symlink()
                or not resolved_child.is_dir()
                or resolved_child.parent != site_root
                or resolved_child != child
            ):
                _fail(
                    "distribution_metadata_invalid",
                    "repair metadata whose identity escapes the selected environment",
                )
            claim = _metadata_claim(resolved_child, site_root, path_mappings)
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


def _record_inventory_owners(
    site_root: Path,
    pth_path: Path,
    spec: AuditedExecutablePth,
    raw: bytes,
) -> tuple[tuple[str, str], ...] | None:
    expected_hash = "sha256=" + base64.urlsafe_b64encode(
        hashlib.sha256(raw).digest()
    ).decode("ascii").rstrip("=")
    owners: list[tuple[str, str]] = []
    try:
        metadata_paths = tuple(
            sorted(
                (
                    path
                    for path in site_root.iterdir()
                    if path.name.endswith(_METADATA_SUFFIXES)
                ),
                key=lambda path: path.name,
            )
        )
    except OSError:
        return None
    for metadata_path in metadata_paths:
        try:
            resolved = metadata_path.resolve(strict=True)
            if (
                metadata_path.is_symlink()
                or resolved != metadata_path
                or resolved.parent != site_root
                or not resolved.is_dir()
            ):
                return None
            record_path = resolved / "RECORD"
            if record_path.is_symlink():
                return None
            if not record_path.exists():
                continue
            if record_path.resolve(strict=True).parent != resolved:
                return None
            rows = tuple(
                csv.reader(record_path.read_text(encoding="utf-8").splitlines())
            )
        except (OSError, RuntimeError, UnicodeError, csv.Error):
            return None
        matching_rows = tuple(row for row in rows if row and row[0] == pth_path.name)
        if not matching_rows:
            continue
        if len(matching_rows) != 1 or len(matching_rows[0]) != 3:
            return None
        _, recorded_hash, recorded_size = matching_rows[0]
        if recorded_hash != expected_hash or recorded_size != str(len(raw)):
            return None
        metadata_name = (
            "METADATA" if resolved.name.endswith(".dist-info") else "PKG-INFO"
        )
        try:
            identity_path = resolved / metadata_name
            if (
                identity_path.is_symlink()
                or identity_path.resolve(strict=True).parent != resolved
            ):
                return None
            identity = _metadata_identity(
                identity_path.read_text(encoding="utf-8"),
                known=False,
            )
        except (OSError, RuntimeError, UnicodeError):
            return None
        if identity is None:
            return None
        name, version = identity
        if (
            _normalized_distribution_name(name)
            != _normalized_distribution_name(spec.distribution)
            or version not in spec.versions
        ):
            return None
        owners.append((_normalized_distribution_name(name), version))
        if len(owners) > 1:
            return None
    return tuple(owners)


def _conda_prefix_for_site_root(site_root: Path) -> tuple[bool, Path | None]:
    if site_root.name != "site-packages":
        return True, None
    python_root = site_root.parent
    if python_root.name == "Lib":
        prefix = python_root.parent
    elif python_root.parent.name == "lib" and re.fullmatch(
        r"python\d+\.\d+", python_root.name
    ):
        prefix = python_root.parent.parent
    else:
        return True, None
    candidate = prefix / "conda-meta"
    try:
        mode = candidate.lstat().st_mode
    except FileNotFoundError:
        for parent in prefix.parents:
            try:
                (parent / "conda-meta").lstat()
            except FileNotFoundError:
                continue
            except OSError:
                return False, None
            return False, None
        return True, None
    except OSError:
        return False, None
    if not stat.S_ISDIR(mode):
        return False, None
    try:
        metadata = candidate.resolve(strict=True)
    except (OSError, RuntimeError):
        return False, None
    if metadata != candidate or metadata.parent != prefix or not metadata.is_dir():
        return False, None
    return True, prefix


def _conda_inventory_owners(
    site_root: Path,
    pth_path: Path,
    spec: AuditedExecutablePth,
    raw: bytes,
) -> tuple[tuple[str, str], ...] | None:
    prefix_valid, prefix = _conda_prefix_for_site_root(site_root)
    if not prefix_valid:
        return None
    if prefix is None:
        return ()
    try:
        relative_path = pth_path.relative_to(prefix).as_posix()
        inventory_paths = {
            relative_path,
            f"{site_root.name}/{pth_path.name}",
        }
        metadata_root = (prefix / "conda-meta").resolve(strict=True)
        metadata_paths = tuple(sorted(metadata_root.glob("*.json")))
    except (OSError, RuntimeError, ValueError):
        return None
    owners: list[tuple[str, str]] = []
    expected_digest = hashlib.sha256(raw).hexdigest()
    for metadata_path in metadata_paths:
        try:
            resolved = metadata_path.resolve(strict=True)
            if (
                metadata_path.is_symlink()
                or resolved != metadata_path
                or resolved.parent != metadata_root
                or not resolved.is_file()
            ):
                return None
            metadata = json.loads(resolved.read_text(encoding="utf-8"))
        except (OSError, RuntimeError, UnicodeError, ValueError):
            return None
        if not isinstance(metadata, dict):
            return None
        paths_data = metadata.get("paths_data")
        if not isinstance(paths_data, dict):
            continue
        paths = paths_data.get("paths")
        if not isinstance(paths, list):
            return None
        matching_paths = tuple(
            entry
            for entry in paths
            if isinstance(entry, dict) and entry.get("_path") in inventory_paths
        )
        if not matching_paths:
            continue
        if len(matching_paths) != 1:
            return None
        entry = matching_paths[0]
        name = metadata.get("name")
        version = metadata.get("version")
        build = metadata.get("build")
        if (
            not isinstance(name, str)
            or not isinstance(version, str)
            or not isinstance(build, str)
            or _normalized_distribution_name(name)
            != _normalized_distribution_name(spec.distribution)
            or (version, build) not in spec.conda_builds
            or entry.get("path_type") != "hardlink"
            or entry.get("sha256") != expected_digest
            or entry.get("sha256_in_prefix") != expected_digest
            or entry.get("size_in_bytes") != len(raw)
        ):
            return None
        owners.append((_normalized_distribution_name(name), version))
        if len(owners) > 1:
            return None
    if not owners:
        return None
    return tuple(owners)


def _audited_nonproduct_executable_pth(
    path: Path,
    raw: bytes,
    site_root: Path,
) -> bool:
    spec = _AUDITED_EXECUTABLE_PTH.get(path.name)
    if spec is None or hashlib.sha256(raw).hexdigest() != spec.sha256:
        return False
    record_owners = _record_inventory_owners(site_root, path, spec, raw)
    conda_owners = _conda_inventory_owners(site_root, path, spec, raw)
    if record_owners is None or conda_owners is None:
        return False
    owners = (*record_owners, *conda_owners)
    if not owners or len(record_owners) > 1 or len(conda_owners) > 1:
        return False
    expected_name = _normalized_distribution_name(spec.distribution)
    return len(set(owners)) == 1 and all(
        name == expected_name and version in spec.versions for name, version in owners
    )


def _reject_startup_customizations(site_root: Path) -> None:
    hook_stems = ("sitecustomize", "usercustomize")
    import_suffixes = tuple(
        suffix.casefold()
        for suffix in (
            *machinery.SOURCE_SUFFIXES,
            *machinery.BYTECODE_SUFFIXES,
            *machinery.EXTENSION_SUFFIXES,
        )
    )
    try:
        children = tuple(site_root.iterdir())
    except OSError:
        _fail(
            "startup_hook_unsafe",
            "repair the selected environment startup metadata before retrying",
        )
    for child in children:
        name = child.name.casefold()
        if any(
            name == stem or any(name == f"{stem}{suffix}" for suffix in import_suffixes)
            for stem in hook_stems
        ):
            _fail(
                "startup_hook_unsafe",
                "remove environment startup customization before retrying",
            )


def _pth_mappings(
    site_roots: tuple[Path, ...],
) -> tuple[tuple[Path, Path], ...]:
    path_mappings: list[tuple[Path, Path]] = []
    for site_root in site_roots:
        _reject_startup_customizations(site_root)
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
                resolved_pth = pth_path.resolve(strict=True)
                if (
                    pth_path.is_symlink()
                    or resolved_pth.parent != site_root
                    or resolved_pth != pth_path
                    or not resolved_pth.is_file()
                ):
                    raise OSError
                raw_bytes = resolved_pth.read_bytes()
                raw = raw_bytes.decode("utf-8")
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
            executable_lines = tuple(
                line
                for line in (raw_line.strip() for raw_line in raw.splitlines())
                if line and re.match(r"^import(?:[ \t]|$)", line)
            )
            if executable_lines and not _audited_nonproduct_executable_pth(
                resolved_pth,
                raw_bytes,
                site_root,
            ):
                _fail(
                    "pth_mapping_unsafe",
                    _PTH_INVENTORY_GUIDANCE,
                )
            for raw_line in raw.splitlines():
                line = raw_line.strip()
                if not line or line.startswith("#"):
                    continue
                if re.match(r"^import(?:[ \t]|$)", line):
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
                if not mapping.exists():
                    if "__editable__" in pth_path.name.lower():
                        _fail(
                            "editable_mapping_invalid",
                            "repair the recorded editable source mapping before retrying",
                        )
                    continue
                path_mappings.append((pth_path.resolve(), mapping))
    return tuple(path_mappings)


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
        if not _metadata_path_looks_known(metadata_path):
            continue
        try:
            resolved_metadata = metadata_path.resolve(strict=True)
        except (OSError, RuntimeError):
            _fail(
                "distribution_metadata_invalid",
                "repair source-tree build metadata before retrying",
            )
        if (
            metadata_path.is_symlink()
            or not resolved_metadata.is_dir()
            or resolved_metadata.parent != layout.source_root
            or resolved_metadata != metadata_path
        ):
            _fail(
                "distribution_metadata_invalid",
                "repair source-tree build metadata before retrying",
            )
        raw_direct_url = _read_optional_text(
            resolved_metadata / "direct_url.json",
            known=False,
            container=resolved_metadata,
        )
        if raw_direct_url is None:
            continue
        direct_url = _parse_direct_url(raw_direct_url, known=True)
        if _editable_root(direct_url) != layout.repository_root:
            _fail(
                "distribution_source_mismatch",
                "repair source metadata that refers to another checkout",
            )

    mappings = _pth_mappings(site_roots)
    (claim,) = _physical_claimants(site_roots, mappings)
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
    mappings = _pth_mappings(site_roots)
    (claim,) = _physical_claimants(site_roots, mappings)
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
    if mappings:
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
