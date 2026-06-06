# Stage 45: Artifact Dataclass Spike — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Introduce a frozen Python `Artifact` dataclass matching `artifact-inventory.yaml` for docs/tests validation only — no runtime artifactization.

**Architecture:** `workflow/lib/artifact.py` with `Artifact` dataclass (frozen, 13 fields), `validate_artifact()` (dict → list[str]), and `load_artifacts()` (YAML → list[Artifact]). Tests import via `sys.path.insert(0, os.path.join(REPO_ROOT, "workflow"))`.

**Tech Stack:** Python 3 stdlib (`dataclasses`, `typing`), PyYAML 6.0.3 (loader only).

---

### Task 1: Create `workflow/lib/` package

**Files:**
- Create: `workflow/lib/__init__.py`

- [ ] **Step 1: Create the package directory and marker**

```bash
mkdir -p workflow/lib
```

- [ ] **Step 2: Write `workflow/lib/__init__.py`**

```python
# workflow/lib — pipeline support library (stdlib-only)
```

---

### Task 2: Create `workflow/lib/artifact.py`

**Files:**
- Create: `workflow/lib/artifact.py`

- [ ] **Step 1: Write the complete module**

```python
"""Stage 45: Artifact dataclass spike — docs/tests validation only.

Provides:
  - Artifact: frozen dataclass matching artifact-inventory.yaml
  - validate_artifact(entry: dict) -> list[str]
  - load_artifacts(path: str) -> list[Artifact]
"""

from dataclasses import dataclass
from typing import Optional


# ---------------------------------------------------------------------------
# Dataclass
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

VALID_SCOPES = {"sample", "experiment", "project", "reference"}
VALID_LEVELS = {"per_sample", "pooled_experiment", "project", "reference"}
VALID_ASSAY_GATES = {
    "all", "peak_centric", "mnase", "cuttag", "chipseq", "atac", "idr",
}

REQUIRED_FIELDS = {
    "id", "description", "scope", "level", "assay_gate",
    "path_template", "producing_rule", "tool",
    "manifest_output_type", "pipeline_done", "rule_all",
    "config_gate", "notes",
}


def validate_artifact(entry: dict) -> list[str]:
    """Validate a raw dict before constructing Artifact.

    Returns a list of error message strings. Empty list means valid.
    """
    errors = []
    eid = entry.get("id", "?")

    # 13 fields exactly
    actual = set(entry.keys())
    missing = REQUIRED_FIELDS - actual
    extra = actual - REQUIRED_FIELDS
    if missing:
        errors.append(f"id={eid}: missing fields {sorted(missing)}")
    if extra:
        errors.append(f"id={eid}: extra fields {sorted(extra)}")

    # Enum-like fields
    if entry.get("scope") not in VALID_SCOPES:
        errors.append(
            f"id={eid}: invalid scope '{entry.get('scope')}'"
        )
    if entry.get("level") not in VALID_LEVELS:
        errors.append(
            f"id={eid}: invalid level '{entry.get('level')}'"
        )
    if entry.get("assay_gate") not in VALID_ASSAY_GATES:
        errors.append(
            f"id={eid}: invalid assay_gate '{entry.get('assay_gate')}'"
        )

    # path_template must start with results/
    pt = entry.get("path_template", "")
    if not pt.startswith("results/"):
        errors.append(
            f"id={eid}: path_template must start with 'results/'"
        )

    # manifest_output_type is str or None
    mot = entry.get("manifest_output_type")
    if mot is not None and not isinstance(mot, str):
        errors.append(
            f"id={eid}: manifest_output_type must be str or None"
        )

    # pipeline_done and rule_all must be bool
    if not isinstance(entry.get("pipeline_done"), bool):
        errors.append(f"id={eid}: pipeline_done must be bool")
    if not isinstance(entry.get("rule_all"), bool):
        errors.append(f"id={eid}: rule_all must be bool")

    return errors


# ---------------------------------------------------------------------------
# Loader
# ---------------------------------------------------------------------------

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

- [ ] **Step 2: Verify import succeeds**

```bash
python3 -c "
import sys; sys.path.insert(0, 'workflow')
from lib.artifact import Artifact, validate_artifact, load_artifacts
print('import OK')
"
```

- [ ] **Step 3: Verify basic construction**

```bash
python3 -c "
import sys; sys.path.insert(0, 'workflow')
from lib.artifact import Artifact
a = Artifact(
    id='test', description='desc', scope='sample', level='per_sample',
    assay_gate='all', path_template='results/test.bam', producing_rule='test_rule',
    tool='test_tool', manifest_output_type=None, pipeline_done=False,
    rule_all=True, config_gate=None, notes=''
)
print(f'id={a.id}, frozen={a}')
# Verify frozen
try:
    a.id = 'mutated'
    print('FAIL: should have raised FrozenInstanceError')
except Exception as e:
    print(f'frozen OK: {type(e).__name__}')
"
```

- [ ] **Step 4: Verify loader round-trip**

```bash
python3 -c "
import os, sys
sys.path.insert(0, 'workflow')
from lib.artifact import load_artifacts
arts = load_artifacts('docs/architecture/artifact-inventory.yaml')
print(f'Loaded {len(arts)} artifacts')
assert len(arts) == 62, f'expected 62, got {len(arts)}'
print('loader OK')
"
```

---

### Task 3: Create `test/test_stage45_artifact_model.py`

**Files:**
- Create: `test/test_stage45_artifact_model.py`

- [ ] **Step 1: Write the test file with all 22 checks**

```python
#!/usr/bin/env python3
"""Stage 45 artifact model tests.

Validates the Artifact dataclass, loader, and validate_artifact() helper.
"""

import os
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO_ROOT, "workflow"))

from lib.artifact import (
    Artifact,
    load_artifacts,
    validate_artifact,
    VALID_SCOPES,
    VALID_LEVELS,
    VALID_ASSAY_GATES,
)

_PASS = 0
_FAIL = 0


def _check(name, condition, detail=""):
    global _PASS, _FAIL
    if condition:
        print(f"PASS: {name}")
        _PASS += 1
    else:
        msg = f"FAIL: {name}"
        if detail:
            msg += f"  -- {detail}"
        print(msg)
        _FAIL += 1


def main():
    global _PASS, _FAIL

    inv_path = os.path.join(REPO_ROOT, "docs", "architecture",
                            "artifact-inventory.yaml")

    # 1. Load returns non-empty list
    artifacts = load_artifacts(inv_path)
    _check("1-load_artifacts_returns_list",
           isinstance(artifacts, list) and len(artifacts) > 0,
           f"got {type(artifacts).__name__} with {len(artifacts) if isinstance(artifacts, list) else 'N/A'} entries")

    # 2. Entry count
    _check("2-entry_count", len(artifacts) == 62,
           f"expected 62, got {len(artifacts)}")

    # 3. All entries are Artifact instances
    all_art = all(isinstance(a, Artifact) for a in artifacts)
    bad_types = [a for a in artifacts if not isinstance(a, Artifact)]
    _check("3-all_are_artifacts", all_art,
           f"{len(bad_types)} entries not Artifact instances")

    # 4. Unique IDs
    ids = [a.id for a in artifacts]
    dupes = sorted(set(i for i in ids if ids.count(i) > 1))
    _check("4-unique_ids", len(dupes) == 0,
           f"duplicates: {dupes}")

    # 5. Frozen — mutation raises FrozenInstanceError
    try:
        artifacts[0].id = "mutated"
        _check("5-frozen_no_mutation", False,
               "expected FrozenInstanceError but no exception raised")
    except Exception as e:
        _check("5-frozen_no_mutation",
               type(e).__name__ == "FrozenInstanceError",
               f"got {type(e).__name__}: {e}")

    # 6. 13 fields via dataclasses.fields()
    from dataclasses import fields
    field_names = {f.name for f in fields(Artifact)}
    _check("6-schema_13_fields", field_names == {
        "id", "description", "scope", "level", "assay_gate",
        "path_template", "producing_rule", "tool",
        "manifest_output_type", "pipeline_done", "rule_all",
        "config_gate", "notes",
    })

    # 7. Valid scope
    bad_scope = [(a.id, a.scope) for a in artifacts
                 if a.scope not in VALID_SCOPES]
    _check("7-valid_scope", len(bad_scope) == 0,
           f"invalid scopes: {bad_scope}")

    # 8. Valid level
    bad_level = [(a.id, a.level) for a in artifacts
                 if a.level not in VALID_LEVELS]
    _check("8-valid_level", len(bad_level) == 0,
           f"invalid levels: {bad_level}")

    # 9. Valid assay_gate
    bad_gate = [(a.id, a.assay_gate) for a in artifacts
                if a.assay_gate not in VALID_ASSAY_GATES]
    _check("9-valid_assay_gate", len(bad_gate) == 0,
           f"invalid assay_gates: {bad_gate}")

    # 10. path_template starts with results/
    bad_path = [(a.id, a.path_template) for a in artifacts
                if not a.path_template.startswith("results/")]
    _check("10-path_starts_with_results", len(bad_path) == 0,
           f"bad paths: {bad_path}")

    # 11. manifest_output_type is str or None
    bad_mot = [(a.id, type(a.manifest_output_type).__name__)
               for a in artifacts
               if a.manifest_output_type is not None
               and not isinstance(a.manifest_output_type, str)]
    _check("11-manifest_output_type_str_or_none", len(bad_mot) == 0,
           f"bad types: {bad_mot}")

    # 12. pipeline_done is bool
    bad_pd = [(a.id, type(a.pipeline_done).__name__) for a in artifacts
              if not isinstance(a.pipeline_done, bool)]
    _check("12-pipeline_done_bool", len(bad_pd) == 0,
           f"not bool: {bad_pd}")

    # 13. rule_all is bool
    bad_ra = [(a.id, type(a.rule_all).__name__) for a in artifacts
              if not isinstance(a.rule_all, bool)]
    _check("13-rule_all_bool", len(bad_ra) == 0,
           f"not bool: {bad_ra}")

    # 14. Round-trip from dict
    d = {
        "id": "test_roundtrip",
        "description": "round-trip test",
        "scope": "project",
        "level": "project",
        "assay_gate": "all",
        "path_template": "results/test/test.txt",
        "producing_rule": "test_rule",
        "tool": "test_tool",
        "manifest_output_type": None,
        "pipeline_done": True,
        "rule_all": False,
        "config_gate": None,
        "notes": "",
    }
    a = Artifact(**d)
    _check("14-from_dict_roundtrip",
           a.id == "test_roundtrip" and a.pipeline_done is True
           and a.rule_all is False and a.manifest_output_type is None)

    # 15. validate_artifact rejects bad scope
    bad = dict(d)
    bad["scope"] = "invalid_scope"
    errors = validate_artifact(bad)
    _check("15-validate_rejects_bad_scope", len(errors) > 0,
           f"expected errors, got {errors}")

    # 16. validate_artifact rejects missing fields
    minimal = {"id": "minimal"}
    errors = validate_artifact(minimal)
    _check("16-validate_rejects_missing_fields", len(errors) > 0,
           f"expected errors, got {errors}")

    total = _PASS + _FAIL
    print(f"\n{_PASS} passed, {_FAIL} failed, {total} total")
    sys.exit(0 if _FAIL == 0 else 1)


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run the tests**

```bash
python3 test/test_stage45_artifact_model.py
```
Expected: 16/16 PASS.

---

### Task 4: Update artifact roadmap

**Files:**
- Modify: `docs/architecture/artifact-roadmap.md`

- [ ] **Step 1: Mark Stage 45 as implemented**

Find the line in the roadmap mentioning Stage 45 and update it to show implemented status.

The current text under the "Staged Path" section has:

```
- **Stage 45:** Artifact Dataclass Spike — `workflow/lib/artifact.py` for tests/docs only.
```

Replace with:

```
- **Stage 45** (implemented 2026-06-06): Artifact Dataclass Spike —
  `workflow/lib/artifact.py` (frozen dataclass, loader, validation helpers)
  with `test/test_stage45_artifact_model.py` (22 checks).
```

---

### Task 5: Run existing test suite (regression check)

- [ ] **Step 1: Run all existing tests**

```bash
python3 test/test_stage43_artifact_inventory.py
```
Expected: 10/10 PASS.

```bash
python3 test/test_validation_stress.py
```
Expected: 40/40 PASS.

```bash
python3 test/test_stage28_release_readiness.py
```
Expected: 11/11 PASS.

```bash
python3 test/test_no_hardcoded_paths.py
```
Expected: PASS.

```bash
git diff --check
```
Expected: no output (clean).

---

### Task 6: Commit

- [ ] **Step 1: Stage and commit**

```bash
git add workflow/lib/__init__.py workflow/lib/artifact.py \
  test/test_stage45_artifact_model.py \
  docs/architecture/artifact-roadmap.md \
  docs/superpowers/specs/2026-06-06-stage45-artifact-dataclass-spike-design.md \
  docs/superpowers/plans/2026-06-06-stage45-artifact-dataclass-spike.md

git commit -m "$(cat <<'EOF'
feat: add Artifact dataclass spike (Stage 45)

Introduce frozen Artifact dataclass matching artifact-inventory.yaml
(13 fields), validate_artifact() helper, and load_artifacts() loader
for docs/tests validation only. Not runtime artifactization.
EOF
)"
```

- [ ] **Step 2: Verify commit**

```bash
git log --oneline -1
git status
```
