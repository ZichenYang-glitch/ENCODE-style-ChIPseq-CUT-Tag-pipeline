#!/usr/bin/env python3
"""Fail-closed source provenance checks without importing product modules."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from importlib import metadata, util
import json
from pathlib import Path
import re
import sys
import sysconfig
from urllib.parse import urlparse
from urllib.request import url2pathname


DISTRIBUTION_NAME = "helixweave"
IMPORT_NAMESPACE = "encode_pipeline"
_KNOWN_DISTRIBUTION_NAMES = (
    "helixweave",
    "encode-pipeline",
)


@dataclass(frozen=True)
class SourceProvenanceError(Exception):
    """A public-safe, stable provenance failure."""

    reason_code: str
    guidance: str

    def __str__(self) -> str:
        return f"source provenance check failed [{self.reason_code}]: {self.guidance}"


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


def _distribution_claims_namespace(distribution: metadata.Distribution) -> bool:
    try:
        top_level = distribution.read_text("top_level.txt")
    except (OSError, UnicodeError):
        _fail(
            "distribution_metadata_invalid",
            "repair the Python package inventory metadata before retrying",
        )
    if top_level is not None:
        names = {line.strip() for line in top_level.splitlines() if line.strip()}
        if IMPORT_NAMESPACE in names:
            return True
    try:
        files = distribution.files or ()
    except (csv.Error, OSError, TypeError, UnicodeError, ValueError):
        _fail(
            "distribution_metadata_invalid",
            "repair the Python package inventory metadata before retrying",
        )
    return any(
        tuple(str(file).replace("\\", "/").split("/", 1))[0] == IMPORT_NAMESPACE
        for file in files
    )


def _is_known_distribution(distribution: metadata.Distribution) -> bool:
    try:
        name = distribution.metadata["Name"]
    except (KeyError, OSError, UnicodeError):
        return False
    return (
        isinstance(name, str)
        and _normalized_distribution_name(name) in _KNOWN_DISTRIBUTION_NAMES
    )


def _claiming_distributions() -> tuple[metadata.Distribution, ...]:
    try:
        distributions = tuple(metadata.distributions())
        known_distributions = tuple(
            distribution
            for name in _KNOWN_DISTRIBUTION_NAMES
            for distribution in metadata.distributions(name=name)
        )
    except (OSError, TypeError, ValueError):
        _fail(
            "distribution_metadata_invalid",
            "repair the isolated environment metadata before retrying",
        )
    claiming = [
        distribution
        for distribution in distributions
        if _distribution_claims_namespace(distribution)
        or _is_known_distribution(distribution)
    ]
    claiming.extend(known_distributions)
    if not claiming:
        _fail(
            "distribution_missing",
            "install HelixWeave in the selected source mode before retrying",
        )
    return tuple(claiming)


def _distribution_name(distribution: metadata.Distribution) -> str:
    try:
        name = distribution.metadata["Name"]
    except (KeyError, OSError, UnicodeError):
        name = None
    if not isinstance(name, str) or not name.strip():
        _fail(
            "distribution_metadata_invalid",
            "repair the HelixWeave distribution metadata before retrying",
        )
    return name


def _distribution_root(distribution: metadata.Distribution) -> Path:
    try:
        return Path(distribution.locate_file("")).resolve(strict=True)
    except (OSError, RuntimeError, TypeError, ValueError):
        _fail(
            "distribution_metadata_invalid",
            "repair the HelixWeave distribution metadata before retrying",
        )


def _direct_url(
    distribution: metadata.Distribution,
) -> dict[str, object] | None:
    try:
        raw = distribution.read_text("direct_url.json")
    except (OSError, UnicodeError):
        _fail(
            "distribution_metadata_invalid",
            "repair the HelixWeave direct URL metadata before retrying",
        )
    if raw is None:
        return None
    try:
        parsed = json.loads(raw)
    except (TypeError, ValueError):
        _fail(
            "distribution_metadata_invalid",
            "repair the HelixWeave direct URL metadata before retrying",
        )
    if not isinstance(parsed, dict):
        _fail(
            "distribution_metadata_invalid",
            "repair the HelixWeave direct URL metadata before retrying",
        )
    if not isinstance(parsed.get("url"), str) or not parsed["url"]:
        _fail(
            "distribution_metadata_invalid",
            "repair the HelixWeave direct URL metadata before retrying",
        )
    return parsed


def _editable_root(direct_url: dict[str, object]) -> Path | None:
    directory_info = direct_url.get("dir_info")
    if not isinstance(directory_info, dict):
        return None
    if directory_info.get("editable") is not True:
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
            "reinstall the requested checkout into an isolated environment",
        )


def _verify_distribution_identity(
    distributions: tuple[metadata.Distribution, ...],
) -> None:
    for distribution in distributions:
        if (
            _normalized_distribution_name(_distribution_name(distribution))
            != DISTRIBUTION_NAME
        ):
            _fail(
                "distribution_identity_mismatch",
                "remove conflicting namespace owners from the isolated environment",
            )


def verify_checkout(repository_root: Path) -> None:
    """Prove the import namespace and editable metadata belong to one checkout."""

    root = _resolved_existing_directory(
        repository_root,
        reason_code="repository_root_invalid",
        guidance="provide an existing repository root explicitly",
    )
    source_root = _resolved_existing_directory(
        root / "src",
        reason_code="source_root_invalid",
        guidance="the requested repository must contain its canonical src directory",
    )
    if source_root.parent != root:
        _fail(
            "source_root_invalid",
            "the requested repository src directory must not escape its root",
        )
    package_root = _resolved_existing_directory(
        source_root / IMPORT_NAMESPACE,
        reason_code="source_root_invalid",
        guidance="the requested src directory must contain encode_pipeline",
    )
    if not _is_within(package_root, source_root):
        _fail(
            "source_root_invalid",
            "the requested package directory must not escape its src directory",
        )
    expected_init = _resolved_existing_file(
        package_root / "__init__.py",
        reason_code="source_root_invalid",
        guidance="the requested package must have its canonical initializer",
    )

    spec = _package_spec()
    locations = _resolved_search_locations(spec)
    if len(locations) != 1 or locations[0] != package_root:
        _fail(
            "module_search_location_mismatch",
            "run from one canonical package location in the requested checkout",
        )
    if _resolved_origin(spec) != expected_init:
        _fail(
            "module_origin_mismatch",
            "reinstall the requested checkout into an isolated environment",
        )

    distributions = _claiming_distributions()
    _verify_distribution_identity(distributions)
    for distribution in distributions:
        distribution_root = _distribution_root(distribution)
        direct_url = _direct_url(distribution)
        editable_root = _editable_root(direct_url) if direct_url is not None else None
        if _is_within(distribution_root, source_root):
            if direct_url is not None and editable_root != root:
                _fail(
                    "distribution_source_mismatch",
                    "reinstall the requested checkout into an isolated environment",
                )
            continue
        if editable_root != root:
            _fail(
                "distribution_source_mismatch",
                "reinstall the requested checkout into an isolated environment",
            )


def _installed_site_roots() -> tuple[Path, ...]:
    if sys.prefix == sys.base_prefix:
        _fail(
            "installed_environment_invalid",
            "select installed-artifact mode only in an isolated environment",
        )
    roots: set[Path] = set()
    for name in ("purelib", "platlib"):
        raw_path = sysconfig.get_path(name)
        if not raw_path:
            continue
        root = _resolved_existing_directory(
            Path(raw_path),
            reason_code="installed_environment_invalid",
            guidance="use the current isolated environment site-packages",
        )
        if not _is_within(root, Path(sys.prefix).resolve()):
            _fail(
                "installed_environment_invalid",
                "use the current isolated environment site-packages",
            )
        roots.add(root)
    if not roots:
        _fail(
            "installed_environment_invalid",
            "use an isolated environment with canonical site-packages",
        )
    return tuple(sorted(roots))


def _within_any(path: Path, roots: tuple[Path, ...]) -> bool:
    return any(_is_within(path, root) for root in roots)


def verify_installed_artifact() -> None:
    """Prove an installed wheel or sdist owns the namespace in this environment."""

    site_roots = _installed_site_roots()
    distributions = _claiming_distributions()
    _verify_distribution_identity(distributions)
    for distribution in distributions:
        direct_url = _direct_url(distribution)
        if direct_url is not None:
            if _editable_root(direct_url) is not None:
                _fail(
                    "installed_editable_forbidden",
                    "install a clean wheel or sdist into the isolated environment",
                )
            if "dir_info" in direct_url or not isinstance(
                direct_url.get("archive_info"),
                dict,
            ):
                _fail(
                    "installed_source_mismatch",
                    "install a clean wheel or sdist into the isolated environment",
                )
        if not _within_any(_distribution_root(distribution), site_roots):
            _fail(
                "installed_source_mismatch",
                "remove source checkouts and external metadata from the environment",
            )
        try:
            installed_metadata = distribution.read_text("METADATA")
        except (OSError, UnicodeError):
            installed_metadata = None
        if installed_metadata is None:
            _fail(
                "installed_source_mismatch",
                "install the artifact as site-packages distribution metadata",
            )

    spec = _package_spec()
    locations = _resolved_search_locations(spec)
    if (
        len(locations) != 1
        or not _within_any(locations[0], site_roots)
        or not _within_any(_resolved_origin(spec), site_roots)
    ):
        _fail(
            "installed_source_mismatch",
            "remove source checkouts and external editable paths before retrying",
        )


def require_checkout(repository_root: Path) -> None:
    """Exit safely when a source-checkout entry point has invalid provenance."""

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
        help="verify a clean wheel or sdist installation",
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
