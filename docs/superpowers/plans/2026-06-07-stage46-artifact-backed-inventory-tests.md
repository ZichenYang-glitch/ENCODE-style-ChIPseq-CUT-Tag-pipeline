# Stage 46: Artifact-Backed Inventory Tests — Implementation Plan

> **For agentic workers:** Use superpowers:executing-plans inline. Steps use checkbox (`- [ ]`) syntax.

**Goal:** Refactor `test/test_stage43_artifact_inventory.py` to use `load_artifacts()` and `Artifact` from Stage 45 as its canonical inventory read path.

**Architecture:** Replace manual `yaml.safe_load()` + raw dict iteration with `load_artifacts()` + Artifact attribute access. Derive schema from `dataclasses.fields(Artifact)` instead of a duplicate `REQUIRED_FIELDS` set.

**Tech Stack:** Python 3 stdlib (`ast`, `dataclasses`, `os`, `re`, `sys`), PyYAML (indirect, via `load_artifacts()`).

---

### Task 1: Refactor `test/test_stage43_artifact_inventory.py`

**Files:**
- Modify: `test/test_stage43_artifact_inventory.py`

- [ ] **Step 1: Replace YAML loading with load_artifacts()**

Current (lines 93-99):
```python
    inv_path = os.path.join(REPO_ROOT, "docs", "architecture",
                            "artifact-inventory.yaml")
    import yaml

    with open(inv_path) as fh:
        data = yaml.safe_load(fh)
```

Replacement:
```python
    sys.path.insert(0, os.path.join(REPO_ROOT, "workflow"))
    from lib.artifact import Artifact, load_artifacts
    from dataclasses import fields

    inv_path = os.path.join(REPO_ROOT, "docs", "architecture",
                            "artifact-inventory.yaml")
```

The `import yaml` block moves inside `load_artifacts()`. The test no longer directly imports `yaml`.

- [ ] **Step 2: Replace check 1 (yaml_parses → load_artifacts_succeeds)**

Current (lines 101-108):
```python
    _check("1-yaml_parses", data is not None and "artifacts" in data)
    if data is None or "artifacts" not in data:
        print(f"\n{_PASS} passed, {_FAIL} failed, {_PASS + _FAIL} total")
        sys.exit(1 if _FAIL else 0)

    artifacts = data["artifacts"]
    _check("1b-has_entries", len(artifacts) >= 50,
           f"expected >=50, got {len(artifacts)}")
```

Replacement:
```python
    try:
        artifacts = load_artifacts(inv_path)
    except Exception as e:
        _check("1-load_artifacts_succeeds", False,
               f"load_artifacts raised {type(e).__name__}: {e}")
        print(f"\n{_PASS} passed, {_FAIL} failed, {_PASS + _FAIL} total")
        sys.exit(1)

    _check("1-load_artifacts_succeeds", True)
    _check("1b-has_entries", len(artifacts) >= 50,
           f"expected >=50, got {len(artifacts)}")
```

- [ ] **Step 3: Replace check 2 (unique_ids) — dict access → attribute**

Current (lines 110-119):
```python
    seen = set()
    dupes = []
    for a in artifacts:
        aid = a["id"]
        if aid in seen:
            dupes.append(aid)
        seen.add(aid)
```

Replacement:
```python
    seen = set()
    dupes = []
    for a in artifacts:
        if a.id in seen:
            dupes.append(a.id)
        seen.add(a.id)
```

- [ ] **Step 4: Replace check 3 (no_workstation_paths) — dict access → attribute**

Current (lines 121-129):
```python
    BAD_PREFIXES = ("/home/", "/data/", "/mnt/")
    bad_paths = []
    for a in artifacts:
        pt = a.get("path_template", "")
        if any(pt.startswith(p) for p in BAD_PREFIXES):
            bad_paths.append(f"{a['id']}: {pt}")
```

Replacement:
```python
    BAD_PREFIXES = ("/home/", "/data/", "/mnt/")
    bad_paths = []
    for a in artifacts:
        if any(a.path_template.startswith(p) for p in BAD_PREFIXES):
            bad_paths.append(f"{a.id}: {a.path_template}")
```

- [ ] **Step 5: Replace check 4 (no_future_keywords) — dict access → attribute**

Current (lines 131-142):
```python
    BAD_KEYWORDS = ("Artifact dataclass", "artifact_path()",
                    "paths.smk", "AssayPolicy YAML",
                    "global target resolver")
    bad_notes = []
    for a in artifacts:
        notes = a.get("notes", "") or ""
        for kw in BAD_KEYWORDS:
            if kw.lower() in notes.lower():
                bad_notes.append(f"{a['id']}: mentions '{kw}'")
```

Replacement:
```python
    BAD_KEYWORDS = ("Artifact dataclass", "artifact_path()",
                    "paths.smk", "AssayPolicy YAML",
                    "global target resolver")
    bad_notes = []
    for a in artifacts:
        for kw in BAD_KEYWORDS:
            if kw.lower() in a.notes.lower():
                bad_notes.append(f"{a.id}: mentions '{kw}'")
```

- [ ] **Step 6: Replace checks 5-6 (manifest coverage) — dict access → attribute**

Current lines 151-166 (lines 153, 162-163 change):
```python
    inventory_manifest_types = set()
    for a in artifacts:
        mt = a.get("manifest_output_type")
        if mt is not None:
            inventory_manifest_types.add(mt)

    missing = manifest_types - inventory_manifest_types
    _check("5-manifest_coverage", len(missing) == 0,
           f"manifest types not in inventory: {sorted(missing)}")

    # 6. ...
    non_null_mt = {a["manifest_output_type"] for a in artifacts
                   if a.get("manifest_output_type") is not None}
```

Replacement:
```python
    inventory_manifest_types = {a.manifest_output_type for a in artifacts
                                if a.manifest_output_type is not None}

    missing = manifest_types - inventory_manifest_types
    _check("5-manifest_coverage", len(missing) == 0,
           f"manifest types not in inventory: {sorted(missing)}")

    non_null_mt = {a.manifest_output_type for a in artifacts
                   if a.manifest_output_type is not None}
```

- [ ] **Step 7: Replace checks 7-8 (output-contract cross-ref) — dict access → attribute**

Current lines 181-187 change:
```python
    inventory_ids_all = {a["id"] for a in artifacts}
```

Replacement:
```python
    inventory_ids_all = {a.id for a in artifacts}
```

- [ ] **Step 8: Replace check 9 (schema completeness) — use dataclasses.fields(Artifact)**

Current (lines 191-212):
```python
    REQUIRED_FIELDS = {
        "id", "description", "scope", "level", "assay_gate",
        "path_template", "producing_rule", "tool", "manifest_output_type",
        "pipeline_done", "rule_all", "config_gate", "notes",
    }
    schema_errors = []
    for a in artifacts:
        keys = set(a.keys())
        missing_fields = REQUIRED_FIELDS - keys
        extra_fields = keys - REQUIRED_FIELDS
        if missing_fields:
            schema_errors.append(f"{a['id']}: missing {sorted(missing_fields)}")
        if extra_fields:
            schema_errors.append(f"{a['id']}: extra {sorted(extra_fields)}")
        if not isinstance(a.get("pipeline_done"), bool):
            schema_errors.append(
                f"{a['id']}: pipeline_done not bool ({type(a.get('pipeline_done')).__name__})")
        if not isinstance(a.get("rule_all"), bool):
            schema_errors.append(
                f"{a['id']}: rule_all not bool ({type(a.get('rule_all')).__name__})")
    _check("9-schema_completeness", len(schema_errors) == 0,
           f"schema errors: {schema_errors[:5]}{'...' if len(schema_errors) > 5 else ''}")
```

Replacement:
```python
    artifact_field_names = {f.name for f in fields(Artifact)}
    schema_errors = []
    for a in artifacts:
        if not isinstance(a, Artifact):
            schema_errors.append(f"not an Artifact: {type(a).__name__}")
            continue
        actual_fields = set(vars(a).keys())
        missing = artifact_field_names - actual_fields
        extra = actual_fields - artifact_field_names
        if missing:
            schema_errors.append(f"{a.id}: missing fields {sorted(missing)}")
        if extra:
            schema_errors.append(f"{a.id}: extra fields {sorted(extra)}")
    _check("9-schema_completeness",
           len(artifact_field_names) == 13 and len(schema_errors) == 0,
           f"field count={len(artifact_field_names)}; errors: {schema_errors[:5]}{'...' if len(schema_errors) > 5 else ''}")
```

- [ ] **Step 9: Run updated test**

```bash
python3 test/test_stage43_artifact_inventory.py
```
Expected: 10/10 PASS (checks 1, 1b, 2-9).

---

### Task 2: Update artifact roadmap

**Files:**
- Modify: `docs/architecture/artifact-roadmap.md`

- [ ] **Step 1: Add Stage 46 completion note**

Find the section after Stage 45 and add:

```
### Stage 46: Artifact-Backed Inventory Tests (implemented 2026-06-07)

Test-layer adoption: `test_stage43_artifact_inventory.py` now uses
`load_artifacts()` and `Artifact` from `workflow/lib/artifact.py`
as its canonical inventory read path. No DAG/runtime changes.
```

---

### Task 3: Run full test suite (regression + new)

- [ ] **Step 1: Run all acceptance tests**

```bash
python3 test/test_stage43_artifact_inventory.py
```
Expected: 10/10 PASS.

```bash
python3 test/test_stage45_artifact_model.py
```
Expected: 22/22 PASS.

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
Expected: clean.

```bash
rm -rf workflow/lib/__pycache__
git status --short
```
Expected: only `?? reconstruction-deep-research-report.md:Zone.Identifier` (untracked, not committed).

---

### Task 4: Commit

- [ ] **Step 1: Stage and commit**

```bash
git add test/test_stage43_artifact_inventory.py \
  docs/architecture/artifact-roadmap.md \
  docs/superpowers/specs/2026-06-07-stage46-artifact-backed-inventory-tests-design.md \
  docs/superpowers/plans/2026-06-07-stage46-artifact-backed-inventory-tests.md

git commit -m "$(cat <<'EOF'
test: use Artifact model for inventory checks (Stage 46)

Refactor test_stage43_artifact_inventory.py to use load_artifacts()
and Artifact from workflow/lib/artifact.py as the canonical inventory
read path. Schema completeness derived from dataclasses.fields(Artifact).
No DAG/runtime changes.
EOF
)"
```
