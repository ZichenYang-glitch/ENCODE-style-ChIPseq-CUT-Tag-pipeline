# Stage 45: Artifact Dataclass Spike — Design Spec

**Date:** 2026-06-06
**Branch:** stage45-artifact-dataclass-spike
**Status:** spec written, awaiting plan and approval

## Purpose

Introduce a frozen Python `Artifact` dataclass matching `docs/architecture/artifact-inventory.yaml` for docs/tests validation only. This is a spike — not runtime artifactization.

## Scope

### In scope

1. `workflow/lib/__init__.py` — minimal package marker (new file).
2. `workflow/lib/artifact.py` — frozen `Artifact` dataclass (13 fields), `validate_artifact()` function, `load_artifacts()` loader.
3. `test/test_stage45_artifact_model.py` — 22 checks covering dataclass construction, field types, enumeration validation, immutability, loader round-trip, `validate_artifact()` error cases, and malformed-input hardening.
4. `docs/architecture/artifact-roadmap.md` — minimal update marking Stage 45 as implemented.

### Out of scope

- No `workflow/Snakefile` changes.
- No `workflow/rules/*.smk` changes.
- No `rule all` changes.
- No `artifact_path()` function.
- No Snakemake rules consuming `Artifact`.
- No `make_manifest.py` changes.
- No `output-contract.md` generation from artifacts.
- No `paths.smk` changes.

## Design

### `workflow/lib/__init__.py`

Minimal package marker:

```python
# workflow/lib — pipeline support library (stdlib-only)
```

### `workflow/lib/artifact.py`

Three components: the `Artifact` dataclass, `validate_artifact()`, and `load_artifacts()`.

#### Artifact dataclass

```python
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
```

Design decisions:

- `frozen=True` — immutable after construction. Tests verify `FrozenInstanceError` on mutation attempt.
- Enum-like fields (`scope`, `level`, `assay_gate`) are plain `str`, not `Enum` subclasses. For a validation-only spike, string whitelists are simpler and avoid updating two places when adding values.
- `Optional[str]` for nullable fields (`manifest_output_type`, `config_gate`).
- `bool` for `pipeline_done` and `rule_all`.

#### Validation helper

```python
VALID_SCOPES = {"sample", "experiment", "project", "reference"}
VALID_LEVELS = {"per_sample", "pooled_experiment", "project", "reference"}
VALID_ASSAY_GATES = {
    "all", "peak_centric", "mnase", "cuttag", "chipseq", "atac", "idr",
}


def validate_artifact(entry: dict) -> list[str]:
    """Validate a raw dict before constructing Artifact.

    Returns a list of error message strings. Empty list means valid.
    """
    errors = []

    required = {
        "id", "description", "scope", "level", "assay_gate",
        "path_template", "producing_rule", "tool",
        "manifest_output_type", "pipeline_done", "rule_all",
        "config_gate", "notes",
    }
    actual = set(entry.keys())
    missing = required - actual
    extra = actual - required
    eid = entry.get("id", "?")

    if missing:
        errors.append(f"id={eid}: missing fields {sorted(missing)}")
    if extra:
        errors.append(f"id={eid}: extra fields {sorted(extra)}")

    if entry.get("scope") not in VALID_SCOPES:
        errors.append(f"id={eid}: invalid scope '{entry.get('scope')}'")
    if entry.get("level") not in VALID_LEVELS:
        errors.append(f"id={eid}: invalid level '{entry.get('level')}'")
    if entry.get("assay_gate") not in VALID_ASSAY_GATES:
        errors.append(f"id={eid}: invalid assay_gate '{entry.get('assay_gate')}'")

    pt = entry.get("path_template", "")
    if not pt.startswith("results/"):
        errors.append(f"id={eid}: path_template must start with 'results/'")

    mot = entry.get("manifest_output_type")
    if mot is not None and not isinstance(mot, str):
        errors.append(f"id={eid}: manifest_output_type must be str or None")

    if not isinstance(entry.get("pipeline_done"), bool):
        errors.append(f"id={eid}: pipeline_done must be bool")
    if not isinstance(entry.get("rule_all"), bool):
        errors.append(f"id={eid}: rule_all must be bool")

    return errors
```

`validate_artifact()` accepts a raw dict, not an `Artifact` instance. It returns `list[str]` errors. Artifact construction happens only after validation passes, keeping the dataclass itself free of validation logic.

#### Loader function

```python
def load_artifacts(path: str) -> list[Artifact]:
    """Load and validate artifact-inventory.yaml.

    Returns list of validated Artifact instances.
    Raises ValueError if the file is malformed or any entry fails validation.
    """
    import yaml

    with open(path) as fh:
        data = yaml.safe_load(fh)

    if not isinstance(data, dict) or "artifacts" not in data:
        raise ValueError(f"Invalid artifact inventory: missing 'artifacts' key in {path}")

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
```

### `test/test_stage45_artifact_model.py`

Import mechanism:

```python
import os
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO_ROOT, "workflow"))

from lib.artifact import Artifact, load_artifacts, validate_artifact
```

22 checks following the established test pattern (`_check()` framework, `main()`, exit code):

| # | Check name | What it verifies |
|:---|:---|:---|
| 1 | yaml_parses | `load_artifacts()` returns a non-empty list |
| 2 | entry_count | `len(load_artifacts(...)) == 62` |
| 3 | all_frozen_artifacts | Every entry is an `Artifact` instance |
| 4 | unique_ids | All artifact IDs are unique |
| 5 | frozen_no_mutation | Assigning `artifact.id = "x"` raises `FrozenInstanceError` |
| 6 | schema_13_fields | Every `Artifact` has exactly 13 fields (via `dataclasses.fields()`) |
| 7 | valid_scope | All `scope` in `{sample, experiment, project, reference}` |
| 8 | valid_level | All `level` in `{per_sample, pooled_experiment, project, reference}` |
| 9 | valid_assay_gate | All `assay_gate` in `{all, peak_centric, mnase, cuttag, chipseq, atac, idr}` |
| 10 | path_starts_with_results | All `path_template` start with `results/` |
| 11 | manifest_output_type_str_or_none | All `manifest_output_type` are `str` or `None` |
| 12 | pipeline_done_bool | All `pipeline_done` are `bool` |
| 13 | rule_all_bool | All `rule_all` are `bool` |
| 14 | from_dict_roundtrip | Construct `Artifact` from a minimal dict, verify field access |
| 15 | validate_rejects_bad_scope | `validate_artifact({"id":"x","scope":"bad",...})` returns errors |
| 16 | validate_rejects_missing_fields | `validate_artifact({"id":"x"})` returns errors for missing fields |

### Roadmap update

In `docs/architecture/artifact-roadmap.md`, add one line under Stage 45:

```
- **Stage 45** (implemented 2026-06-06): Artifact Dataclass Spike —
  `workflow/lib/artifact.py` (frozen dataclass, loader, validation helpers)
  with `test/test_stage45_artifact_model.py` (22 checks).
```

## Verification

```bash
# 1. Package structure
test -f workflow/lib/__init__.py && test -f workflow/lib/artifact.py
python3 -c "import sys; sys.path.insert(0, 'workflow'); from lib.artifact import Artifact, load_artifacts, validate_artifact; print('import OK')"

# 2. Dataclass instantiation
python3 -c "
import sys; sys.path.insert(0, 'workflow')
from lib.artifact import Artifact
a = Artifact(id='test', description='', scope='sample', level='per_sample',
    assay_gate='all', path_template='results/test.bam', producing_rule='test',
    tool='test', manifest_output_type=None, pipeline_done=False,
    rule_all=True, config_gate=None, notes='')
assert a.id == 'test'
assert a.frozen == True  # via AttributeError on mutation attempt
print('dataclass OK')
"

# 3. Loader round-trip
python3 -c "
import os, sys
REPO_ROOT = os.path.dirname(os.path.abspath('.'))
sys.path.insert(0, os.path.join(REPO_ROOT, 'workflow'))
from lib.artifact import load_artifacts
arts = load_artifacts('docs/architecture/artifact-inventory.yaml')
print(f'Loaded {len(arts)} artifacts')
assert len(arts) == 62
print('loader OK')
"

# 4. Full test suite
python3 test/test_stage45_artifact_model.py

# 5. Existing tests unchanged
python3 test/test_stage43_artifact_inventory.py
python3 test/test_validation_stress.py
python3 test/test_stage28_release_readiness.py
python3 test/test_no_hardcoded_paths.py
git diff --check
```

## Non-goals (reaffirmed)

- No `workflow/Snakefile` changes.
- No `workflow/rules/*.smk` changes.
- No `rule all` changes.
- No `artifact_path()`.
- No Snakemake rule consumption of `Artifact`.
- No `make_manifest.py` changes.
- No `output-contract.md` generation.
- No `paths.smk` changes.
