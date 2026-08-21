"""Source-owned, fail-closed admission for Alembic migration execution."""

from __future__ import annotations

import ast
from contextlib import contextmanager
from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import stat
import tempfile
from typing import Any, Iterator


MIGRATION_EXECUTION_INVENTORY_FILE = "migration-execution-inventory-1.0.0.json"
MIGRATION_EXECUTION_INVENTORY_ID = "helixweave-platform-migrations"
MIGRATION_EXECUTION_INVENTORY_SCHEMA_VERSION = "1.0.0"
MIGRATION_EXECUTION_INVENTORY_SCHEME = (
    "sha256-canonical-migration-execution-inventory-v1"
)
# This source constant is the trust anchor for the committed inventory bytes.
# Updating the inventory is a two-file review action: regenerate the JSON, then
# deliberately synchronize these exact coordinates.
MIGRATION_EXECUTION_INVENTORY_SIZE_BYTES = 3935
MIGRATION_EXECUTION_INVENTORY_SHA256 = (
    "4bbaaa55bc5b35d00e34432771215effce72b40bfcd7a1a0041da165aec7cef9"
)

_MAXIMUM_INVENTORY_BYTES = 512 * 1024
_MAXIMUM_MIGRATION_BYTES = 2 * 1024 * 1024
_MAXIMUM_REVISIONS = 1024
_SHA256 = re.compile(r"^[0-9a-f]{64}$")
_REVISION = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,127}$")
_REVISION_PATH = re.compile(r"^alembic/versions/[A-Za-z0-9][A-Za-z0-9_.-]*\.py$")
_METADATA_NAMES = (
    "revision",
    "down_revision",
    "branch_labels",
    "depends_on",
)


class MigrationAdmissionError(RuntimeError):
    """A stable, public-safe migration execution rejection."""

    def __init__(self, reason_code: str) -> None:
        self.reason_code = reason_code
        super().__init__(f"migration execution admission failed [{self.reason_code}]")

    def __str__(self) -> str:
        return f"migration execution admission failed [{self.reason_code}]"


@dataclass(frozen=True)
class VerifiedMigrationSource:
    """One exact source file captured before Alembic is importable."""

    path: str
    size_bytes: int
    sha256: str
    content: bytes


@dataclass(frozen=True)
class VerifiedMigrationRevision:
    """One exact admitted Alembic revision and its static graph metadata."""

    source: VerifiedMigrationSource
    revision: str
    down_revision: tuple[str, ...]
    branch_labels: tuple[str, ...]
    depends_on: tuple[str, ...]


@dataclass(frozen=True)
class VerifiedMigrationExecutionInventory:
    """Captured migration bytes that may be materialized for Alembic."""

    inventory_sha256: str
    contract_sha256: str
    environment: VerifiedMigrationSource
    revisions: tuple[VerifiedMigrationRevision, ...]
    bases: tuple[str, ...]
    heads: tuple[str, ...]


@dataclass(frozen=True)
class MigrationInventoryTrustAnchor:
    """Literal size/hash coordinates carried by one platform release."""

    size_bytes: int
    sha256: str

    def __post_init__(self) -> None:
        if (
            not isinstance(self.size_bytes, int)
            or isinstance(self.size_bytes, bool)
            or not 0 < self.size_bytes <= _MAXIMUM_INVENTORY_BYTES
            or not isinstance(self.sha256, str)
            or _SHA256.fullmatch(self.sha256) is None
        ):
            raise MigrationAdmissionError("MIGRATION_EXECUTION_TRUST_ANCHOR_INVALID")


def migration_inventory_trust_anchor_from_source(
    source_bytes: bytes,
) -> MigrationInventoryTrustAnchor:
    """Extract only the two reviewed literal anchors from candidate source.

    Deployment admission uses this to validate a newer wheel without importing
    or executing that wheel.  The outer wheel/bundle identity remains the
    release trust boundary; this parser only preserves the inventory's existing
    two-file review contract across platform versions.
    """
    try:
        content = _bounded_bytes(
            source_bytes,
            maximum_bytes=_MAXIMUM_MIGRATION_BYTES,
            reason_code="MIGRATION_EXECUTION_TRUST_ANCHOR_INVALID",
        )
        tree = ast.parse(content, filename="migration_admission.py")
        values: dict[str, object] = {}
        expected = {
            "MIGRATION_EXECUTION_INVENTORY_SIZE_BYTES",
            "MIGRATION_EXECUTION_INVENTORY_SHA256",
        }
        for statement in tree.body:
            if not isinstance(statement, (ast.Assign, ast.AnnAssign)):
                continue
            targets = (
                statement.targets
                if isinstance(statement, ast.Assign)
                else [statement.target]
            )
            names = [target.id for target in targets if isinstance(target, ast.Name)]
            selected = [name for name in names if name in expected]
            if not selected:
                continue
            if len(selected) != 1 or selected[0] in values:
                raise _AdmissionFailure("MIGRATION_EXECUTION_TRUST_ANCHOR_INVALID")
            value_node = statement.value
            if value_node is None:
                raise _AdmissionFailure("MIGRATION_EXECUTION_TRUST_ANCHOR_INVALID")
            value = ast.literal_eval(value_node)
            values[selected[0]] = value
        if set(values) != expected:
            raise _AdmissionFailure("MIGRATION_EXECUTION_TRUST_ANCHOR_INVALID")
        return MigrationInventoryTrustAnchor(
            size_bytes=values["MIGRATION_EXECUTION_INVENTORY_SIZE_BYTES"],
            sha256=values["MIGRATION_EXECUTION_INVENTORY_SHA256"],
        )
    except MigrationAdmissionError:
        raise
    except _AdmissionFailure as error:
        raise MigrationAdmissionError(error.reason_code) from None
    except (SyntaxError, TypeError, UnicodeError, ValueError):
        raise MigrationAdmissionError(
            "MIGRATION_EXECUTION_TRUST_ANCHOR_INVALID"
        ) from None


class _DuplicateJsonKey(ValueError):
    pass


class _AdmissionFailure(Exception):
    def __init__(self, reason_code: str) -> None:
        self.reason_code = reason_code
        super().__init__(reason_code)


def verify_migration_execution_inventory(
    *,
    persistence_root: Path | None = None,
    inventory_bytes: bytes | None = None,
    trust_anchor: MigrationInventoryTrustAnchor | None = None,
) -> VerifiedMigrationExecutionInventory:
    """Verify and capture every source-owned Alembic env/revision script.

    The default root is the physical package containing this module, which is
    the same location used by :mod:`encode_pipeline.persistence.migrations` in
    both a source checkout and an installed distribution.  No revision module
    is imported while this admission decision is made.
    """

    try:
        anchor = (
            MigrationInventoryTrustAnchor(
                MIGRATION_EXECUTION_INVENTORY_SIZE_BYTES,
                MIGRATION_EXECUTION_INVENTORY_SHA256,
            )
            if trust_anchor is None
            else trust_anchor
        )
        if not isinstance(anchor, MigrationInventoryTrustAnchor):
            raise _AdmissionFailure("MIGRATION_EXECUTION_TRUST_ANCHOR_INVALID")
        root = (
            Path(os.path.abspath(__file__)).parent
            if persistence_root is None
            else _require_absolute_path(persistence_root)
        )
        with _open_absolute_directory(root) as root_descriptor:
            raw_inventory = (
                _read_regular_at(
                    root_descriptor,
                    MIGRATION_EXECUTION_INVENTORY_FILE,
                    maximum_bytes=_MAXIMUM_INVENTORY_BYTES,
                    invalid_reason="MIGRATION_EXECUTION_INVENTORY_INVALID",
                )
                if inventory_bytes is None
                else _bounded_bytes(
                    inventory_bytes,
                    maximum_bytes=_MAXIMUM_INVENTORY_BYTES,
                    reason_code="MIGRATION_EXECUTION_INVENTORY_INVALID",
                )
            )
            if (
                len(raw_inventory) != anchor.size_bytes
                or hashlib.sha256(raw_inventory).hexdigest() != anchor.sha256
            ):
                raise _AdmissionFailure("MIGRATION_EXECUTION_INVENTORY_INVALID")
            parsed = _parse_inventory(raw_inventory)
            environment, revisions = _capture_declared_sources(
                root_descriptor,
                parsed,
            )
        return VerifiedMigrationExecutionInventory(
            inventory_sha256=parsed["inventory_sha256"],
            contract_sha256=hashlib.sha256(raw_inventory).hexdigest(),
            environment=environment,
            revisions=revisions,
            bases=tuple(parsed["bases"]),
            heads=tuple(parsed["heads"]),
        )
    except MigrationAdmissionError:
        raise
    except _AdmissionFailure as error:
        raise MigrationAdmissionError(error.reason_code) from None
    except (OSError, TypeError, UnicodeError, ValueError):
        raise MigrationAdmissionError("MIGRATION_EXECUTION_INVENTORY_INVALID") from None


@contextmanager
def admitted_migration_script_location(
    *,
    persistence_root: Path | None = None,
    inventory_bytes: bytes | None = None,
) -> Iterator[tuple[Path, VerifiedMigrationExecutionInventory]]:
    """Yield a private snapshot containing only inventory-admitted bytes."""

    verified = verify_migration_execution_inventory(
        persistence_root=persistence_root,
        inventory_bytes=inventory_bytes,
    )
    temporary: tempfile.TemporaryDirectory[str] | None = None
    try:
        temporary = tempfile.TemporaryDirectory(
            prefix="helixweave-migration-snapshot-",
            ignore_cleanup_errors=True,
        )
        snapshot_root = Path(temporary.name)
        os.chmod(snapshot_root, 0o700)
        alembic_root = snapshot_root / "alembic"
        versions_root = alembic_root / "versions"
        alembic_root.mkdir(mode=0o700)
        versions_root.mkdir(mode=0o700)
        _write_snapshot_file(
            alembic_root / "env.py",
            verified.environment.content,
        )
        for revision in verified.revisions:
            _write_snapshot_file(
                versions_root / PurePosixPath(revision.source.path).name,
                revision.source.content,
            )
    except OSError:
        if temporary is not None:
            temporary.cleanup()
        raise MigrationAdmissionError("MIGRATION_EXECUTION_SNAPSHOT_FAILED") from None
    try:
        yield alembic_root, verified
    finally:
        temporary.cleanup()


def build_migration_execution_inventory(
    persistence_root: Path,
) -> dict[str, Any]:
    """Build canonical inventory data from a reviewed persistence source tree."""

    try:
        root = _require_absolute_path(persistence_root)
        with _open_absolute_directory(root) as root_descriptor:
            alembic_descriptor = _open_directory_at(
                root_descriptor,
                "alembic",
                reason_code="MIGRATION_EXECUTION_INVENTORY_INVALID",
            )
            try:
                environment_content = _read_regular_at(
                    alembic_descriptor,
                    "env.py",
                    maximum_bytes=_MAXIMUM_MIGRATION_BYTES,
                    invalid_reason="MIGRATION_EXECUTION_INVENTORY_INVALID",
                )
                versions_descriptor = _open_directory_at(
                    alembic_descriptor,
                    "versions",
                    reason_code="MIGRATION_EXECUTION_INVENTORY_INVALID",
                )
                try:
                    names = _revision_source_names(versions_descriptor)
                    revisions: list[dict[str, Any]] = []
                    for name in names:
                        content = _read_regular_at(
                            versions_descriptor,
                            name,
                            maximum_bytes=_MAXIMUM_MIGRATION_BYTES,
                            invalid_reason="MIGRATION_REVISION_UNSAFE_TYPE",
                        )
                        metadata = _parse_revision_metadata(content)
                        revisions.append(
                            {
                                "path": f"alembic/versions/{name}",
                                "size_bytes": len(content),
                                "sha256": hashlib.sha256(content).hexdigest(),
                                "revision": metadata["revision"],
                                "down_revision": list(metadata["down_revision"]),
                                "branch_labels": list(metadata["branch_labels"]),
                                "depends_on": list(metadata["depends_on"]),
                            }
                        )
                finally:
                    os.close(versions_descriptor)
            finally:
                os.close(alembic_descriptor)
        bases, heads = _validate_graph_entries(revisions)
        payload: dict[str, Any] = {
            "schema_version": MIGRATION_EXECUTION_INVENTORY_SCHEMA_VERSION,
            "inventory_id": MIGRATION_EXECUTION_INVENTORY_ID,
            "scheme": MIGRATION_EXECUTION_INVENTORY_SCHEME,
            "environment": {
                "path": "alembic/env.py",
                "size_bytes": len(environment_content),
                "sha256": hashlib.sha256(environment_content).hexdigest(),
            },
            "revision_count": len(revisions),
            "revisions": revisions,
            "bases": list(bases),
            "heads": list(heads),
        }
        return {
            **payload,
            "inventory_sha256": hashlib.sha256(_canonical_json(payload)).hexdigest(),
        }
    except _AdmissionFailure as error:
        raise ValueError("migration execution inventory source is invalid") from error
    except (OSError, TypeError, UnicodeError, ValueError) as error:
        raise ValueError("migration execution inventory source is invalid") from error


def canonical_migration_execution_inventory_bytes(
    document: object,
) -> bytes:
    """Serialize one inventory document deterministically."""

    return _canonical_json(document) + b"\n"


def _capture_declared_sources(
    root_descriptor: int,
    inventory: dict[str, Any],
) -> tuple[
    VerifiedMigrationSource,
    tuple[VerifiedMigrationRevision, ...],
]:
    alembic_descriptor = _open_directory_at(
        root_descriptor,
        "alembic",
        reason_code="MIGRATION_EXECUTION_INVENTORY_INVALID",
    )
    try:
        environment_spec = inventory["environment"]
        environment_content = _read_regular_at(
            alembic_descriptor,
            "env.py",
            maximum_bytes=_MAXIMUM_MIGRATION_BYTES,
            invalid_reason="MIGRATION_EXECUTION_INVENTORY_INVALID",
        )
        environment = _verified_source(
            environment_spec,
            environment_content,
            mismatch_reason="MIGRATION_EXECUTION_INVENTORY_INVALID",
        )
        versions_descriptor = _open_directory_at(
            alembic_descriptor,
            "versions",
            reason_code="MIGRATION_EXECUTION_INVENTORY_INVALID",
        )
        try:
            observed_names = _revision_source_names(versions_descriptor)
            expected_by_name = {
                PurePosixPath(item["path"]).name: item
                for item in inventory["revisions"]
            }
            observed = set(observed_names)
            expected = set(expected_by_name)
            if observed - expected:
                raise _AdmissionFailure("MIGRATION_REVISION_UNKNOWN")
            if expected - observed:
                raise _AdmissionFailure("MIGRATION_REVISION_MISSING")

            revisions: list[VerifiedMigrationRevision] = []
            for name in observed_names:
                item = expected_by_name[name]
                try:
                    content = _read_regular_at(
                        versions_descriptor,
                        name,
                        maximum_bytes=_MAXIMUM_MIGRATION_BYTES,
                        invalid_reason="MIGRATION_REVISION_UNSAFE_TYPE",
                    )
                except FileNotFoundError:
                    raise _AdmissionFailure("MIGRATION_REVISION_MISSING") from None
                source = _verified_source(
                    item,
                    content,
                    mismatch_reason="MIGRATION_REVISION_DIGEST_MISMATCH",
                )
                metadata = _parse_revision_metadata(content)
                if (
                    metadata["revision"] != item["revision"]
                    or metadata["down_revision"] != tuple(item["down_revision"])
                    or metadata["branch_labels"] != tuple(item["branch_labels"])
                    or metadata["depends_on"] != tuple(item["depends_on"])
                ):
                    raise _AdmissionFailure("MIGRATION_REVISION_GRAPH_INVALID")
                revisions.append(
                    VerifiedMigrationRevision(
                        source=source,
                        revision=metadata["revision"],
                        down_revision=metadata["down_revision"],
                        branch_labels=metadata["branch_labels"],
                        depends_on=metadata["depends_on"],
                    )
                )
        finally:
            os.close(versions_descriptor)
    finally:
        os.close(alembic_descriptor)
    return environment, tuple(revisions)


def _verified_source(
    specification: dict[str, Any],
    content: bytes,
    *,
    mismatch_reason: str,
) -> VerifiedMigrationSource:
    observed_sha256 = hashlib.sha256(content).hexdigest()
    if (
        len(content) != specification["size_bytes"]
        or observed_sha256 != specification["sha256"]
    ):
        raise _AdmissionFailure(mismatch_reason)
    return VerifiedMigrationSource(
        path=specification["path"],
        size_bytes=len(content),
        sha256=observed_sha256,
        content=content,
    )


def _parse_inventory(content: bytes) -> dict[str, Any]:
    try:
        document = json.loads(
            content.decode("utf-8"),
            object_pairs_hook=_unique_object,
        )
    except (
        UnicodeDecodeError,
        json.JSONDecodeError,
        _DuplicateJsonKey,
    ) as error:
        raise _AdmissionFailure("MIGRATION_EXECUTION_INVENTORY_INVALID") from error
    expected_fields = {
        "schema_version",
        "inventory_id",
        "scheme",
        "environment",
        "revision_count",
        "revisions",
        "bases",
        "heads",
        "inventory_sha256",
    }
    if not isinstance(document, dict) or set(document) != expected_fields:
        raise _AdmissionFailure("MIGRATION_EXECUTION_INVENTORY_INVALID")
    if content != canonical_migration_execution_inventory_bytes(document):
        raise _AdmissionFailure("MIGRATION_EXECUTION_INVENTORY_INVALID")
    supplied_identity = document["inventory_sha256"]
    payload = {
        key: value for key, value in document.items() if key != "inventory_sha256"
    }
    if (
        document["schema_version"] != MIGRATION_EXECUTION_INVENTORY_SCHEMA_VERSION
        or document["inventory_id"] != MIGRATION_EXECUTION_INVENTORY_ID
        or document["scheme"] != MIGRATION_EXECUTION_INVENTORY_SCHEME
        or not isinstance(supplied_identity, str)
        or _SHA256.fullmatch(supplied_identity) is None
        or supplied_identity != hashlib.sha256(_canonical_json(payload)).hexdigest()
    ):
        raise _AdmissionFailure("MIGRATION_EXECUTION_INVENTORY_INVALID")

    environment = document["environment"]
    _validate_source_spec(
        environment,
        expected_path="alembic/env.py",
        reason_code="MIGRATION_EXECUTION_INVENTORY_INVALID",
    )
    revisions = document["revisions"]
    revision_count = document["revision_count"]
    if (
        isinstance(revision_count, bool)
        or not isinstance(revision_count, int)
        or revision_count <= 0
        or revision_count > _MAXIMUM_REVISIONS
        or not isinstance(revisions, list)
        or len(revisions) != revision_count
    ):
        raise _AdmissionFailure("MIGRATION_EXECUTION_INVENTORY_INVALID")

    normalized_paths: set[str] = set()
    revision_ids: set[str] = set()
    previous_path: str | None = None
    for item in revisions:
        if not isinstance(item, dict) or set(item) != {
            "path",
            "size_bytes",
            "sha256",
            "revision",
            "down_revision",
            "branch_labels",
            "depends_on",
        }:
            raise _AdmissionFailure("MIGRATION_EXECUTION_INVENTORY_INVALID")
        _validate_source_spec(
            item,
            expected_path=None,
            reason_code="MIGRATION_EXECUTION_INVENTORY_INVALID",
        )
        path = item["path"]
        revision = item["revision"]
        normalized_path = path.casefold()
        if (
            _REVISION_PATH.fullmatch(path) is None
            or PurePosixPath(path).name == "__init__.py"
            or previous_path is not None
            and path <= previous_path
        ):
            raise _AdmissionFailure("MIGRATION_EXECUTION_INVENTORY_INVALID")
        if normalized_path in normalized_paths or revision in revision_ids:
            raise _AdmissionFailure("MIGRATION_REVISION_DUPLICATE")
        normalized_paths.add(normalized_path)
        revision_ids.add(revision)
        previous_path = path
        if not isinstance(revision, str) or _REVISION.fullmatch(revision) is None:
            raise _AdmissionFailure("MIGRATION_EXECUTION_INVENTORY_INVALID")
        for field in ("down_revision", "branch_labels", "depends_on"):
            _validate_revision_list(item[field])

    computed_bases, computed_heads = _validate_graph_entries(revisions)
    bases = document["bases"]
    heads = document["heads"]
    _validate_revision_list(bases, require_nonempty=True)
    _validate_revision_list(heads, require_nonempty=True)
    if (
        tuple(bases) != computed_bases
        or tuple(heads) != computed_heads
        or len(bases) != 1
        or len(heads) != 1
    ):
        raise _AdmissionFailure("MIGRATION_REVISION_GRAPH_INVALID")
    return document


def _validate_source_spec(
    value: object,
    *,
    expected_path: str | None,
    reason_code: str,
) -> None:
    if not isinstance(value, dict):
        raise _AdmissionFailure(reason_code)
    if expected_path is not None and set(value) != {
        "path",
        "size_bytes",
        "sha256",
    }:
        raise _AdmissionFailure(reason_code)
    path = value.get("path")
    size_bytes = value.get("size_bytes")
    file_sha256 = value.get("sha256")
    if (
        not isinstance(path, str)
        or expected_path is not None
        and path != expected_path
        or isinstance(size_bytes, bool)
        or not isinstance(size_bytes, int)
        or size_bytes <= 0
        or size_bytes > _MAXIMUM_MIGRATION_BYTES
        or not isinstance(file_sha256, str)
        or _SHA256.fullmatch(file_sha256) is None
    ):
        raise _AdmissionFailure(reason_code)


def _validate_revision_list(
    value: object,
    *,
    require_nonempty: bool = False,
) -> None:
    if (
        not isinstance(value, list)
        or require_nonempty
        and not value
        or any(
            not isinstance(item, str) or _REVISION.fullmatch(item) is None
            for item in value
        )
        or value != sorted(set(value))
    ):
        raise _AdmissionFailure("MIGRATION_REVISION_GRAPH_INVALID")


def _validate_graph_entries(
    revisions: list[dict[str, Any]],
) -> tuple[tuple[str, ...], tuple[str, ...]]:
    if not revisions or len(revisions) > _MAXIMUM_REVISIONS:
        raise _AdmissionFailure("MIGRATION_REVISION_GRAPH_INVALID")
    by_revision: dict[str, dict[str, Any]] = {}
    for item in revisions:
        revision = item.get("revision")
        if not isinstance(revision, str) or _REVISION.fullmatch(revision) is None:
            raise _AdmissionFailure("MIGRATION_REVISION_GRAPH_INVALID")
        if revision in by_revision:
            raise _AdmissionFailure("MIGRATION_REVISION_DUPLICATE")
        by_revision[revision] = item

    branch_labels: set[str] = set()
    for item in revisions:
        for label in tuple(item.get("branch_labels", ())):
            if label in by_revision or label in branch_labels:
                raise _AdmissionFailure("MIGRATION_REVISION_GRAPH_INVALID")
            branch_labels.add(label)

    referenced_as_parent: set[str] = set()
    for revision, item in by_revision.items():
        down_revision = tuple(item.get("down_revision", ()))
        depends_on = tuple(item.get("depends_on", ()))
        for dependency in (*down_revision, *depends_on):
            if dependency == revision or dependency not in by_revision:
                raise _AdmissionFailure("MIGRATION_REVISION_GRAPH_INVALID")
        referenced_as_parent.update(down_revision)

    dependency_counts: dict[str, int] = {}
    dependents: dict[str, list[str]] = {revision: [] for revision in by_revision}
    for revision, item in by_revision.items():
        dependencies = (
            *tuple(item.get("down_revision", ())),
            *tuple(item.get("depends_on", ())),
        )
        dependency_counts[revision] = len(dependencies)
        for dependency in dependencies:
            dependents[dependency].append(revision)
    ready = sorted(
        revision
        for revision, dependency_count in dependency_counts.items()
        if dependency_count == 0
    )
    visited_count = 0
    while ready:
        revision = ready.pop()
        visited_count += 1
        for dependent in dependents[revision]:
            dependency_counts[dependent] -= 1
            if dependency_counts[dependent] == 0:
                ready.append(dependent)
    if visited_count != len(by_revision):
        raise _AdmissionFailure("MIGRATION_REVISION_GRAPH_INVALID")

    bases = tuple(
        sorted(
            revision
            for revision, item in by_revision.items()
            if not tuple(item.get("down_revision", ()))
        )
    )
    heads = tuple(sorted(set(by_revision) - referenced_as_parent))
    if len(bases) != 1 or len(heads) != 1:
        raise _AdmissionFailure("MIGRATION_REVISION_GRAPH_INVALID")

    reachable: set[str] = set()
    pending = [heads[0]]
    while pending:
        revision = pending.pop()
        if revision in reachable:
            continue
        reachable.add(revision)
        pending.extend(tuple(by_revision[revision].get("down_revision", ())))
    if reachable != set(by_revision):
        raise _AdmissionFailure("MIGRATION_REVISION_GRAPH_INVALID")
    return bases, heads


def _parse_revision_metadata(content: bytes) -> dict[str, object]:
    try:
        source = content.decode("utf-8")
        tree = ast.parse(source, filename="<migration-revision>")
    except (SyntaxError, UnicodeDecodeError) as error:
        raise _AdmissionFailure("MIGRATION_REVISION_GRAPH_INVALID") from error
    assignments: dict[str, object] = {}
    for node in tree.body:
        name: str | None = None
        value_node: ast.expr | None = None
        if isinstance(node, ast.Assign) and len(node.targets) == 1:
            target = node.targets[0]
            if isinstance(target, ast.Name):
                name = target.id
                value_node = node.value
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
            name = node.target.id
            value_node = node.value
        if name not in _METADATA_NAMES:
            continue
        if name in assignments or value_node is None:
            raise _AdmissionFailure("MIGRATION_REVISION_DUPLICATE")
        try:
            assignments[name] = ast.literal_eval(value_node)
        except (SyntaxError, ValueError) as error:
            raise _AdmissionFailure("MIGRATION_REVISION_GRAPH_INVALID") from error
    if set(assignments) != set(_METADATA_NAMES):
        raise _AdmissionFailure("MIGRATION_REVISION_GRAPH_INVALID")
    revision = assignments["revision"]
    if not isinstance(revision, str) or _REVISION.fullmatch(revision) is None:
        raise _AdmissionFailure("MIGRATION_REVISION_GRAPH_INVALID")
    return {
        "revision": revision,
        "down_revision": _normalize_revision_literal(assignments["down_revision"]),
        "branch_labels": _normalize_revision_literal(assignments["branch_labels"]),
        "depends_on": _normalize_revision_literal(assignments["depends_on"]),
    }


def _normalize_revision_literal(value: object) -> tuple[str, ...]:
    if value is None:
        return ()
    if isinstance(value, str):
        values = (value,)
    elif isinstance(value, (tuple, list)):
        values = tuple(value)
    else:
        raise _AdmissionFailure("MIGRATION_REVISION_GRAPH_INVALID")
    if (
        any(
            not isinstance(item, str) or _REVISION.fullmatch(item) is None
            for item in values
        )
        or tuple(sorted(set(values))) != values
    ):
        raise _AdmissionFailure("MIGRATION_REVISION_GRAPH_INVALID")
    return values


def _revision_source_names(versions_descriptor: int) -> tuple[str, ...]:
    try:
        names = tuple(
            sorted(
                name
                for name in os.listdir(versions_descriptor)
                if name.endswith(".py") and name != "__init__.py"
            )
        )
    except OSError as error:
        raise _AdmissionFailure("MIGRATION_EXECUTION_INVENTORY_INVALID") from error
    if not names or len(names) > _MAXIMUM_REVISIONS:
        raise _AdmissionFailure("MIGRATION_REVISION_GRAPH_INVALID")
    normalized: set[str] = set()
    for name in names:
        if "/" in name or "\\" in name or name.casefold() in normalized:
            raise _AdmissionFailure("MIGRATION_REVISION_DUPLICATE")
        normalized.add(name.casefold())
        try:
            info = os.stat(
                name,
                dir_fd=versions_descriptor,
                follow_symlinks=False,
            )
        except OSError as error:
            raise _AdmissionFailure("MIGRATION_REVISION_UNSAFE_TYPE") from error
        if not stat.S_ISREG(info.st_mode) or info.st_nlink != 1:
            raise _AdmissionFailure("MIGRATION_REVISION_UNSAFE_TYPE")
    return names


def _read_regular_at(
    directory_descriptor: int,
    name: str,
    *,
    maximum_bytes: int,
    invalid_reason: str,
) -> bytes:
    if (
        not isinstance(name, str)
        or not name
        or "/" in name
        or "\\" in name
        or name in {".", ".."}
    ):
        raise _AdmissionFailure(invalid_reason)
    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0)
    nofollow = getattr(os, "O_NOFOLLOW", None)
    nonblock = getattr(os, "O_NONBLOCK", None)
    if nofollow is None or nonblock is None:
        raise _AdmissionFailure(invalid_reason)
    flags |= nofollow | nonblock
    try:
        descriptor = os.open(
            name,
            flags,
            dir_fd=directory_descriptor,
        )
    except FileNotFoundError:
        raise
    except OSError as error:
        raise _AdmissionFailure(invalid_reason) from error
    try:
        before = os.fstat(descriptor)
        if (
            not stat.S_ISREG(before.st_mode)
            or before.st_nlink != 1
            or before.st_size <= 0
            or before.st_size > maximum_bytes
        ):
            raise _AdmissionFailure(invalid_reason)
        chunks: list[bytes] = []
        remaining = before.st_size
        while remaining:
            chunk = os.read(descriptor, min(remaining, 1024 * 1024))
            if not chunk:
                raise _AdmissionFailure(invalid_reason)
            chunks.append(chunk)
            remaining -= len(chunk)
        if os.read(descriptor, 1):
            raise _AdmissionFailure(invalid_reason)
        after = os.fstat(descriptor)
        if _file_identity(before) != _file_identity(after):
            raise _AdmissionFailure(invalid_reason)
        return b"".join(chunks)
    finally:
        os.close(descriptor)


@contextmanager
def _open_absolute_directory(path: Path) -> Iterator[int]:
    if not isinstance(path, Path) or not path.is_absolute():
        raise _AdmissionFailure("MIGRATION_EXECUTION_INVENTORY_INVALID")
    directory_flags = _directory_flags()
    descriptors: list[int] = []
    try:
        current = os.open("/", directory_flags)
        descriptors.append(current)
        for component in path.parts[1:]:
            descriptor = os.open(
                component,
                directory_flags,
                dir_fd=current,
            )
            descriptors.append(descriptor)
            info = os.fstat(descriptor)
            if not stat.S_ISDIR(info.st_mode):
                raise _AdmissionFailure("MIGRATION_EXECUTION_INVENTORY_INVALID")
            current = descriptor
        yield current
    except _AdmissionFailure:
        raise
    except OSError as error:
        raise _AdmissionFailure("MIGRATION_EXECUTION_INVENTORY_INVALID") from error
    finally:
        for descriptor in reversed(descriptors):
            try:
                os.close(descriptor)
            except OSError:
                pass


def _open_directory_at(
    directory_descriptor: int,
    name: str,
    *,
    reason_code: str,
) -> int:
    if not isinstance(name, str) or "/" in name or "\\" in name:
        raise _AdmissionFailure(reason_code)
    descriptor: int | None = None
    transferred = False
    try:
        descriptor = os.open(
            name,
            _directory_flags(),
            dir_fd=directory_descriptor,
        )
        if not stat.S_ISDIR(os.fstat(descriptor).st_mode):
            raise _AdmissionFailure(reason_code)
        transferred = True
        return descriptor
    except _AdmissionFailure:
        raise
    except OSError as error:
        raise _AdmissionFailure(reason_code) from error
    finally:
        if descriptor is not None and not transferred:
            try:
                os.close(descriptor)
            except OSError:
                pass


def _directory_flags() -> int:
    nofollow = getattr(os, "O_NOFOLLOW", None)
    directory = getattr(os, "O_DIRECTORY", None)
    if nofollow is None or directory is None:
        raise _AdmissionFailure("MIGRATION_EXECUTION_INVENTORY_INVALID")
    return os.O_RDONLY | nofollow | directory | getattr(os, "O_CLOEXEC", 0)


def _write_snapshot_file(path: Path, content: bytes) -> None:
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_CLOEXEC", 0)
    descriptor = os.open(path, flags, 0o400)
    try:
        view = memoryview(content)
        while view:
            written = os.write(descriptor, view)
            if written <= 0:
                raise OSError("snapshot write failed")
            view = view[written:]
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _file_identity(info: os.stat_result) -> tuple[int, int, int, int, int]:
    return (
        info.st_dev,
        info.st_ino,
        info.st_mode,
        info.st_nlink,
        info.st_size,
    )


def _require_absolute_path(value: object) -> Path:
    if not isinstance(value, Path) or not value.is_absolute():
        raise _AdmissionFailure("MIGRATION_EXECUTION_INVENTORY_INVALID")
    return value


def _bounded_bytes(
    value: object,
    *,
    maximum_bytes: int,
    reason_code: str,
) -> bytes:
    if not isinstance(value, bytes) or not value or len(value) > maximum_bytes:
        raise _AdmissionFailure(reason_code)
    return value


def _unique_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise _DuplicateJsonKey
        result[key] = value
    return result


def _canonical_json(value: object) -> bytes:
    return json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
