#!/usr/bin/env python3
"""Audit the pinned nf-core/rnaseq process-to-container universe read-only."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import stat
import sys
from typing import Any


PROJECT_ROOT = Path(__file__).resolve().parents[1]
CONTRACT_ROOT = PROJECT_ROOT / "src/encode_pipeline/contracts/nfcore_rnaseq"
SOURCE_MANIFEST = CONTRACT_ROOT / "source-manifest-3.26.0.json"
CONTAINER_INVENTORY = CONTRACT_ROOT / "container-inventory-3.26.0.json"
COMMITTED_AUDIT = CONTRACT_ROOT / "container-process-audit-3.26.0.json"

PROJECT = "nf-core/rnaseq"
RELEASE = "3.26.0"
COMMIT = "e7ca46272c8f9d5ceee3f71759f4ba551d3217a4"
SOURCE_TREE_SHA256 = "4f779bd8934d41896fbf137ff31158b02daacd7efbb40ed8cee55b9c8f757722"
RESERVED_DEFAULT_DENY_LABEL = "helixweave_verified_container_v1"
PROCESS_UNIVERSE_COUNT = 78
INCLUDED_PROCESS_COUNT = 56
EXCLUDED_PROCESS_COUNT = 22
CONFIG_CONTAINER_ASSIGNMENT_COUNT = 78
DEFAULT_CONFIG_CONTAINER_ASSIGNMENT_COUNT = 1

EXCLUSIONS = {
    "BBMAP_BBSPLIT": "unsupported_contaminant_screening",
    "BRACKEN_BRACKEN": "unsupported_contaminant_screening",
    "CUSTOM_CATADDITIONALFASTA": "unsupported_additional_fasta",
    "CUSTOM_RSEMMERGECOUNTS": "unsupported_rsem_route",
    "HISAT2_ALIGN": "unsupported_hisat2_route",
    "HISAT2_BUILD": "unsupported_hisat2_route",
    "HISAT2_EXTRACTSPLICESITES": "unsupported_hisat2_route",
    "KALLISTO_INDEX": "unsupported_kallisto_route",
    "KALLISTO_QUANT": "unsupported_kallisto_route",
    "KRAKEN2_KRAKEN2": "unsupported_contaminant_screening",
    "PARABRICKS_RNAFQ2BAM": "unsupported_accelerated_aligner",
    "RIBODETECTOR": "unsupported_ribodetector",
    "RSEM_CALCULATEEXPRESSION": "unsupported_rsem_route",
    "RSEM_PREPAREREFERENCE": "platform_owned_transcript_fasta",
    "RUSTQC": "unsupported_rustqc",
    "SENTIEON_RSEMCALCULATEEXPRESSION": "unsupported_commercial_aligner",
    "SENTIEON_RSEMPREPAREREFERENCE": "unsupported_commercial_aligner",
    "SENTIEON_STARALIGN": "unsupported_commercial_aligner",
    "SEQKIT_STATS": "unsupported_ribodetector_support",
    "SYLPHTAX_TAXPROF": "unsupported_contaminant_screening",
    "SYLPH_PROFILE": "unsupported_contaminant_screening",
    "UMICOLLAPSE": "unsupported_umi_dedup_tool",
}

_PROCESS_RE = re.compile(r"(?m)^\s*process\s+([A-Z][A-Z0-9_]*)\s*\{")
_INPUT_RE = re.compile(r"(?m)^\s*input\s*:")
_SINGLE_QUOTED_RE = re.compile(r"'([^']+)'")
_FIXED_CONTAINER_RE = re.compile(r'(?m)^\s*container\s+"([^"$\n]+)"\s*$')
_CONFIG_CONTAINER_RE = re.compile(r"(?m)^\s*container\s*=")
_WITH_NAME_RE = re.compile(r"(?m)^\s*withName:\s*(['\"])(.*?)\1\s*\{")
_ALIASED_INCLUDE_RE = re.compile(
    r"include\s*\{\s*([A-Z][A-Z0-9_]*)\s+as\s+([A-Z][A-Z0-9_]*)\s*\}",
    re.MULTILINE,
)

_DEFAULT_CONTAINER_CONFIG = "conf/modules/prepare_genome.config"
_PROFILE_CONTAINER_CONFIGS = frozenset({"conf/arm.config"})
_EXPECTED_DEFAULT_CONTAINER_OVERRIDES = {
    ".*PARABRICKS_STARGENOMEGENERATE": "unsupported_parabricks_alias",
}


class AuditError(ValueError):
    """The staged source does not match the committed immutable contract."""


def _sha256(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def _canonical_sha256(value: object) -> str:
    return _sha256(
        json.dumps(
            value,
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
        ).encode("utf-8")
    )


def _load_object(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_bytes())
    if not isinstance(value, dict):
        raise AuditError(f"{path.name} must contain a JSON object")
    return value


def _normalise_image_coordinate(value: str) -> str:
    if "@" in value or any(char.isspace() for char in value):
        raise AuditError(f"invalid image coordinate: {value!r}")
    final = value.rsplit("/", 1)[-1]
    if ":" not in final:
        raise AuditError(f"unversioned image coordinate: {value!r}")
    first = value.split("/", 1)[0]
    if "." not in first and ":" not in first and first != "localhost":
        value = f"quay.io/{value}"
    if "/" not in value:
        raise AuditError(f"unqualified image coordinate: {value!r}")
    return value


def _declared_image_coordinates(content: str) -> tuple[str, ...]:
    input_match = _INPUT_RE.search(content)
    header = content[: input_match.start()] if input_match else content
    fixed = _FIXED_CONTAINER_RE.search(header)
    if fixed:
        return (_normalise_image_coordinate(fixed.group(1)),)
    container_offset = header.find("container")
    if container_offset < 0:
        raise AuditError("process has no container directive")
    candidates: list[str] = []
    for value in _SINGLE_QUOTED_RE.findall(header[container_offset:]):
        final = value.rsplit("/", 1)[-1]
        if value.startswith(("http://", "https://")) or ":" not in final:
            continue
        candidates.append(_normalise_image_coordinate(value))
    coordinates = tuple(sorted(set(candidates)))
    if not coordinates:
        raise AuditError("process has no Docker image coordinate")
    return coordinates


def _verify_source_tree(
    source_root: Path,
    source_manifest: dict[str, Any],
) -> dict[str, str]:
    expected_entries = source_manifest.get("files")
    if not isinstance(expected_entries, list):
        raise AuditError("source manifest files must be a list")
    expected = {entry["path"]: entry for entry in expected_entries}
    if len(expected) != len(expected_entries):
        raise AuditError("source manifest contains duplicate paths")

    actual: dict[str, dict[str, object]] = {}
    source_hashes: dict[str, str] = {}
    for path in sorted(source_root.rglob("*")):
        relative = path.relative_to(source_root).as_posix()
        info = path.lstat()
        if stat.S_ISDIR(info.st_mode):
            continue
        if not stat.S_ISREG(info.st_mode):
            raise AuditError(f"source entry is not a regular file: {relative}")
        content = path.read_bytes()
        if RESERVED_DEFAULT_DENY_LABEL.encode("ascii") in content:
            raise AuditError("reserved default-deny label already exists upstream")
        digest = _sha256(content)
        actual[relative] = {
            "path": relative,
            "size_bytes": len(content),
            "sha256": digest,
            "executable": bool(info.st_mode & 0o111),
        }
        source_hashes[relative] = digest
    if set(actual) != set(expected):
        raise AuditError("source file set differs from the pinned manifest")
    if any(actual[path] != expected[path] for path in actual):
        raise AuditError("source file identity differs from the pinned manifest")
    ordered_paths = sorted(actual, key=PurePosixPath)
    if (
        _canonical_sha256([actual[path] for path in ordered_paths])
        != SOURCE_TREE_SHA256
    ):
        raise AuditError("source tree digest differs from the pinned identity")
    return source_hashes


def _container_config_assignments(
    source_root: Path,
    source_hashes: dict[str, str],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Audit every config-level container assignment in the pinned source.

    Only assignments in the pipeline's default-loaded config graph can affect
    the fixed ``docker`` profile.  Profile-only and upstream-test configs are
    nevertheless included in the all-assignment digest so a newly introduced
    assignment cannot escape review.
    """
    aliases: list[tuple[str, str, str]] = []
    for source_file in sorted(source_root.rglob("*.nf")):
        relative = source_file.relative_to(source_root).as_posix()
        for source_process, alias in _ALIASED_INCLUDE_RE.findall(
            source_file.read_text(encoding="utf-8")
        ):
            aliases.append((alias, source_process, relative))

    assignments: list[dict[str, object]] = []
    defaults: list[dict[str, object]] = []
    assignment_files: set[str] = set()
    for config_file in sorted(source_root.rglob("*.config")):
        relative = config_file.relative_to(source_root).as_posix()
        content = config_file.read_text(encoding="utf-8")
        container_matches = list(_CONFIG_CONTAINER_RE.finditer(content))
        if not container_matches:
            continue
        assignment_files.add(relative)
        selector_matches = list(_WITH_NAME_RE.finditer(content))
        for container_match in container_matches:
            preceding = [
                match
                for match in selector_matches
                if match.start() < container_match.start()
            ]
            if not preceding:
                raise AuditError(
                    f"config container assignment has no withName selector: {relative}"
                )
            selector_match = preceding[-1]
            next_selector = next(
                (
                    match
                    for match in selector_matches
                    if match.start() > container_match.start()
                ),
                None,
            )
            block_end = next_selector.start() if next_selector else len(content)
            block = content[container_match.start() : block_end]
            image_coordinates = []
            for value in _SINGLE_QUOTED_RE.findall(block):
                final = value.rsplit("/", 1)[-1]
                if value.startswith(("http://", "https://")) or ":" not in final:
                    continue
                image_coordinates.append(_normalise_image_coordinate(value))
            coordinates = sorted(set(image_coordinates))
            if not coordinates:
                raise AuditError(
                    f"config container assignment has no Docker image: {relative}"
                )
            selector = selector_match.group(2)
            matching_aliases = []
            try:
                selector_re = re.compile(selector)
            except re.error as exc:
                raise AuditError(f"invalid withName selector: {selector!r}") from exc
            for alias, source_process, alias_file in sorted(aliases):
                if selector_re.fullmatch(alias) is None:
                    continue
                matching_aliases.append(
                    {
                        "alias": alias,
                        "source_process": source_process,
                        "source_file": alias_file,
                        "source_file_sha256": source_hashes[alias_file],
                    }
                )
            if relative == _DEFAULT_CONTAINER_CONFIG:
                load_scope = "default"
            elif relative in _PROFILE_CONTAINER_CONFIGS:
                load_scope = "unselected_arm64_profile"
            elif "tests" in PurePosixPath(relative).parts:
                load_scope = "upstream_test_only"
            else:
                raise AuditError(
                    f"unclassified config container assignment: {relative}"
                )
            entry: dict[str, object] = {
                "source_file": relative,
                "source_file_sha256": source_hashes[relative],
                "selector": selector,
                "image_coordinates": coordinates,
                "matched_aliases": matching_aliases,
                "load_scope": load_scope,
            }
            assignments.append(entry)
            if load_scope == "default":
                reason = _EXPECTED_DEFAULT_CONTAINER_OVERRIDES.get(selector)
                if reason is None:
                    raise AuditError(
                        f"default config container selector is not owned: {selector}"
                    )
                if not matching_aliases:
                    raise AuditError(
                        f"default config container selector has no audited alias: {selector}"
                    )
                default_entry = {
                    **entry,
                    "platform_override": "deny",
                    "override_reason": reason,
                }
                defaults.append(default_entry)

    expected_assignment_files = {
        _DEFAULT_CONTAINER_CONFIG,
        *_PROFILE_CONTAINER_CONFIGS,
        "modules/local/star_genomeparams_upgrade/tests/nextflow.legacy_index.config",
        "modules/nf-core/parabricks/rnafq2bam/tests/nextflow.config",
        "subworkflows/local/align_star/tests/nextflow.upgraded_legacy.config",
    }
    if assignment_files != expected_assignment_files:
        raise AuditError("config container assignment file set differs")
    if len(assignments) != CONFIG_CONTAINER_ASSIGNMENT_COUNT:
        raise AuditError("config container assignment count differs")
    if len(defaults) != DEFAULT_CONFIG_CONTAINER_ASSIGNMENT_COUNT:
        raise AuditError("default config container assignment count differs")
    if {entry["selector"] for entry in defaults} != set(
        _EXPECTED_DEFAULT_CONTAINER_OVERRIDES
    ):
        raise AuditError("default config container selector set differs")
    return assignments, defaults


def build_audit(source_root: Path) -> dict[str, Any]:
    """Build the deterministic audit value without modifying either tree."""
    source_manifest = _load_object(SOURCE_MANIFEST)
    inventory = _load_object(CONTAINER_INVENTORY)
    source_hashes = _verify_source_tree(source_root, source_manifest)
    config_assignments, default_config_assignments = _container_config_assignments(
        source_root,
        source_hashes,
    )
    inventory_entries = inventory.get("entries")
    if not isinstance(inventory_entries, list):
        raise AuditError("container inventory entries must be a list")
    included = {entry["process"]: entry for entry in inventory_entries}
    if len(included) != INCLUDED_PROCESS_COUNT:
        raise AuditError("included process count differs")
    if set(included) & set(EXCLUSIONS):
        raise AuditError("included and excluded process sets overlap")

    discovered: dict[str, dict[str, object]] = {}
    for source_file in sorted(source_root.rglob("*.nf")):
        relative = source_file.relative_to(source_root)
        if "tests" in relative.parts:
            continue
        content = source_file.read_text(encoding="utf-8")
        processes = _PROCESS_RE.findall(content)
        if not processes:
            continue
        if len(processes) != 1:
            raise AuditError(f"expected one process in {relative.as_posix()}")
        process = processes[0]
        if process in discovered:
            raise AuditError(f"duplicate process declaration: {process}")
        discovered[process] = {
            "process": process,
            "source_file": relative.as_posix(),
            "source_file_sha256": source_hashes[relative.as_posix()],
            "image_coordinates": list(_declared_image_coordinates(content)),
        }

    expected_processes = set(included) | set(EXCLUSIONS)
    if set(discovered) != expected_processes:
        missing = sorted(expected_processes - set(discovered))
        extra = sorted(set(discovered) - expected_processes)
        raise AuditError(f"process universe differs; missing={missing}, extra={extra}")
    if len(discovered) != PROCESS_UNIVERSE_COUNT:
        raise AuditError("process universe count differs")

    processes: list[dict[str, object]] = []
    for process in sorted(discovered):
        entry = discovered[process]
        if process in included:
            inventory_entry = included[process]
            expected = {
                "process": process,
                "source_file": inventory_entry["source_file"],
                "source_file_sha256": inventory_entry["source_file_sha256"],
                "image_coordinates": [inventory_entry["image_coordinate"]],
            }
            if entry != expected:
                raise AuditError(f"inventory differs from source for {process}")
            entry["disposition"] = "included"
        else:
            entry["disposition"] = "excluded"
            entry["exclusion_reason"] = EXCLUSIONS[process]
        processes.append(entry)

    return {
        "schema_version": "1.1.0",
        "project": PROJECT,
        "release": RELEASE,
        "commit": COMMIT,
        "source_tree_sha256": SOURCE_TREE_SHA256,
        "reserved_default_deny_label": RESERVED_DEFAULT_DENY_LABEL,
        "process_universe_count": PROCESS_UNIVERSE_COUNT,
        "included_process_count": INCLUDED_PROCESS_COUNT,
        "excluded_process_count": EXCLUDED_PROCESS_COUNT,
        "config_container_assignment_count": CONFIG_CONTAINER_ASSIGNMENT_COUNT,
        "config_container_assignments_sha256": _canonical_sha256(config_assignments),
        "default_config_container_assignment_count": (
            DEFAULT_CONFIG_CONTAINER_ASSIGNMENT_COUNT
        ),
        "default_config_container_assignments_sha256": _canonical_sha256(
            default_config_assignments
        ),
        "default_config_container_assignments": default_config_assignments,
        "processes_sha256": _canonical_sha256(processes),
        "processes": processes,
    }


def _canonical_bytes(value: object) -> bytes:
    return (json.dumps(value, indent=2, ensure_ascii=False) + "\n").encode("utf-8")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_root", type=Path)
    parser.add_argument(
        "--check",
        action="store_true",
        help="compare with the committed audit instead of emitting JSON",
    )
    args = parser.parse_args(argv)
    source_root = args.source_root.resolve()
    if not source_root.is_dir() or source_root.is_symlink():
        parser.error("source_root must be a real directory")
    try:
        content = _canonical_bytes(build_audit(source_root))
    except (AuditError, OSError, UnicodeError, json.JSONDecodeError) as exc:
        print(f"audit failed: {exc}", file=sys.stderr)
        return 1
    if args.check:
        try:
            committed = COMMITTED_AUDIT.read_bytes()
        except OSError as exc:
            print(f"audit failed: {exc}", file=sys.stderr)
            return 1
        if committed != content:
            print("committed container process audit differs", file=sys.stderr)
            return 1
        audit = json.loads(content)
        print(f"processes={PROCESS_UNIVERSE_COUNT}")
        print(f"config_container_assignments={CONFIG_CONTAINER_ASSIGNMENT_COUNT}")
        print(
            "config_container_assignments_sha256="
            f"{audit['config_container_assignments_sha256']}"
        )
        print(
            "default_config_container_assignments="
            f"{DEFAULT_CONFIG_CONTAINER_ASSIGNMENT_COUNT}"
        )
        print(
            "default_config_container_assignments_sha256="
            f"{audit['default_config_container_assignments_sha256']}"
        )
        print(f"audit_sha256={_sha256(content)}")
        return 0
    sys.stdout.buffer.write(content)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
