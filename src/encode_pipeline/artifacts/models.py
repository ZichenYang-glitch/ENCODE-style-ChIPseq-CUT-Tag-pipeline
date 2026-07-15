"""Artifact catalog data model and validation."""

from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True)
class Artifact:
    """Immutable representation of one pipeline artifact entry."""

    id: str
    description: str
    scope: str
    level: str
    assay_gate: str
    path_template: str
    producing_rule: str
    tool: str
    manifest_output_type: Optional[str]
    pipeline_done: bool
    rule_all: bool
    config_gate: Optional[str]
    notes: str


VALID_SCOPES = {"sample", "experiment", "project", "reference"}
VALID_LEVELS = {"per_sample", "pooled_experiment", "project", "reference"}
VALID_ASSAY_GATES = {
    "all",
    "peak_centric",
    "mnase",
    "cuttag",
    "chipseq",
    "atac",
    "idr",
}

REQUIRED_FIELDS = {
    "id",
    "description",
    "scope",
    "level",
    "assay_gate",
    "path_template",
    "producing_rule",
    "tool",
    "manifest_output_type",
    "pipeline_done",
    "rule_all",
    "config_gate",
    "notes",
}

REQUIRED_STR_FIELDS = {
    "id",
    "description",
    "scope",
    "level",
    "assay_gate",
    "path_template",
    "producing_rule",
    "tool",
    "notes",
}


def validate_artifact(entry: dict) -> list[str]:
    """Validate a raw dict before constructing Artifact.

    Returns a list of error message strings. Empty list means valid.
    Never raises — returns errors for all malformed input.
    """
    if not isinstance(entry, dict):
        return [f"entry must be dict, got {type(entry).__name__}"]

    errors = []
    eid = entry.get("id", "?")
    if not isinstance(eid, str):
        eid = repr(eid)

    try:
        actual = set(entry.keys())
        missing = REQUIRED_FIELDS - actual
        extra = actual - REQUIRED_FIELDS
    except AttributeError:
        return [f"entry must be dict, got {type(entry).__name__}"]
    if missing:
        errors.append(f"id={eid}: missing fields {sorted(missing)}")
    if extra:
        errors.append(f"id={eid}: extra fields {sorted(extra)}")

    for field in sorted(REQUIRED_STR_FIELDS):
        val = entry.get(field)
        if not isinstance(val, str):
            errors.append(f"id={eid}: {field} must be str, got {type(val).__name__}")

    if isinstance(entry.get("scope"), str) and entry["scope"] not in VALID_SCOPES:
        errors.append(f"id={eid}: invalid scope '{entry['scope']}'")
    if isinstance(entry.get("level"), str) and entry["level"] not in VALID_LEVELS:
        errors.append(f"id={eid}: invalid level '{entry['level']}'")
    if (
        isinstance(entry.get("assay_gate"), str)
        and entry["assay_gate"] not in VALID_ASSAY_GATES
    ):
        errors.append(f"id={eid}: invalid assay_gate '{entry['assay_gate']}'")

    pt = entry.get("path_template")
    if isinstance(pt, str) and not pt.startswith("results/"):
        errors.append(f"id={eid}: path_template must start with 'results/'")

    mot = entry.get("manifest_output_type")
    if mot is not None and not isinstance(mot, str):
        errors.append(f"id={eid}: manifest_output_type must be str or None")

    cg = entry.get("config_gate")
    if cg is not None and not isinstance(cg, str):
        errors.append(f"id={eid}: config_gate must be str or None")

    if not isinstance(entry.get("pipeline_done"), bool):
        errors.append(
            f"id={eid}: pipeline_done must be bool, got "
            f"{type(entry.get('pipeline_done')).__name__}"
        )
    if not isinstance(entry.get("rule_all"), bool):
        errors.append(
            f"id={eid}: rule_all must be bool, got "
            f"{type(entry.get('rule_all')).__name__}"
        )

    return errors


def artifacts_by_id(artifacts):
    """Index artifacts by id. Raises ValueError on duplicate id."""
    result = {}
    for a in artifacts:
        if a.id in result:
            raise ValueError(f"duplicate artifact id: {a.id}")
        result[a.id] = a
    return result


def artifacts_by_manifest_output_type(artifacts):
    """Index artifacts by manifest_output_type (non-null only).

    Raises ValueError on duplicate manifest_output_type.
    """
    result = {}
    for a in artifacts:
        mt = a.manifest_output_type
        if mt is None:
            continue
        if mt in result:
            raise ValueError(f"duplicate manifest_output_type: {mt}")
        result[mt] = a
    return result


def filter_artifacts(
    artifacts,
    *,
    assay_gate=None,
    scope=None,
    level=None,
    pipeline_done=None,
    rule_all=None,
):
    """Filter artifacts by field values. None means no filter."""
    results = list(artifacts)
    if assay_gate is not None:
        results = [a for a in results if a.assay_gate == assay_gate]
    if scope is not None:
        results = [a for a in results if a.scope == scope]
    if level is not None:
        results = [a for a in results if a.level == level]
    if pipeline_done is not None:
        results = [a for a in results if a.pipeline_done == pipeline_done]
    if rule_all is not None:
        results = [a for a in results if a.rule_all == rule_all]
    return results


def load_artifacts(path: str) -> list[Artifact]:
    """Load and validate artifact-inventory.yaml.

    Returns list of validated Artifact instances.
    Raises ValueError if the file is malformed or any entry fails validation.
    """
    import yaml

    with open(path) as fh:
        data = yaml.safe_load(fh)

    if not isinstance(data, dict) or "artifacts" not in data:
        raise ValueError(
            f"Invalid artifact inventory: missing 'artifacts' key in {path}"
        )
    if not isinstance(data["artifacts"], list):
        raise ValueError(
            f"Invalid artifact inventory: 'artifacts' must be a list in {path}"
        )

    all_errors = []
    artifacts = []
    for entry in data["artifacts"]:
        errors = validate_artifact(entry)
        if errors:
            all_errors.extend(errors)
        else:
            artifacts.append(Artifact(**entry))

    if all_errors:
        raise ValueError("\n".join(all_errors))

    return artifacts
