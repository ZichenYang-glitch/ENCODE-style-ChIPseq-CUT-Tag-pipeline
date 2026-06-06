# Stage 48: Artifact Query Helpers — Implementation Plan

> **For agentic workers:** Use superpowers:executing-plans inline.

**Goal:** Add 3 query helpers to `workflow/lib/artifact.py` and refactor 3 test files to use them.

---

### Task 1: Add query helpers to `workflow/lib/artifact.py`

**Files:**
- Modify: `workflow/lib/artifact.py`

- [ ] **Step 1: Add `artifacts_by_id()` after the `REQUIRED_STR_FIELDS` block**

```python
def artifacts_by_id(artifacts):
    """Index artifacts by id. Raises ValueError on duplicate id.

    Args:
        artifacts: iterable of Artifact objects.

    Returns:
        dict[str, Artifact]: artifact id -> Artifact mapping.
    """
    result = {}
    for a in artifacts:
        if a.id in result:
            raise ValueError(f"duplicate artifact id: {a.id}")
        result[a.id] = a
    return result


def artifacts_by_manifest_output_type(artifacts):
    """Index artifacts by manifest_output_type (non-null only).

    Args:
        artifacts: iterable of Artifact objects.

    Returns:
        dict[str, Artifact]: manifest_output_type -> Artifact mapping.
        Entries with manifest_output_type=None are excluded.

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


def filter_artifacts(artifacts, *, assay_gate=None, scope=None,
                     level=None, pipeline_done=None, rule_all=None):
    """Filter artifacts by field values. None means no filter.

    Args:
        artifacts: iterable of Artifact objects.
        assay_gate: filter by assay_gate value.
        scope: filter by scope value.
        level: filter by level value.
        pipeline_done: filter by pipeline_done bool.
        rule_all: filter by rule_all bool.

    Returns:
        list[Artifact]: filtered artifacts.
    """
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
```

---

### Task 2: Add query-helper tests to `test/test_stage45_artifact_model.py`

- [ ] **Step 1: Add imports for new helpers**

```python
from lib.artifact import (
    Artifact, load_artifacts, validate_artifact,
    VALID_SCOPES, VALID_LEVELS, VALID_ASSAY_GATES,
    artifacts_by_id, artifacts_by_manifest_output_type, filter_artifacts,
)
```

- [ ] **Step 2: Add checks 23-28 after existing check 22**

```python
    # 23. artifacts_by_id returns 62 entries
    by_id = artifacts_by_id(artifacts)
    _check("23-artifacts_by_id_62", len(by_id) == 62,
           f"expected 62, got {len(by_id)}")

    # 24. duplicate id raises ValueError
    dup_test = [Artifact(**d), Artifact(**d)]
    try:
        artifacts_by_id(dup_test)
        _check("24-by_id_duplicate_raises", False,
               "expected ValueError but no exception raised")
    except ValueError as e:
        _check("24-by_id_duplicate_raises",
               "duplicate" in str(e).lower(),
               f"got ValueError: {e}")

    # 25. artifacts_by_manifest_output_type — keys + values match
    by_mot = artifacts_by_manifest_output_type(artifacts)
    non_null_mot = {a.manifest_output_type for a in artifacts
                    if a.manifest_output_type is not None}
    _check("25-by_mot_keys_match", set(by_mot.keys()) == non_null_mot,
           f"key mismatch")
    _check("26-by_mot_values_match",
           all(v.manifest_output_type == k for k, v in by_mot.items()),
           "value.manifest_output_type != key")

    # 27. duplicate manifest_output_type raises ValueError
    dup_mot_arts = [
        Artifact(id="a", description="", scope="sample", level="per_sample",
                 assay_gate="all", path_template="results/a.bam",
                 producing_rule="r", tool="t",
                 manifest_output_type="final_bam",
                 pipeline_done=False, rule_all=True, config_gate=None, notes=""),
        Artifact(id="b", description="", scope="sample", level="per_sample",
                 assay_gate="all", path_template="results/b.bam",
                 producing_rule="r", tool="t",
                 manifest_output_type="final_bam",
                 pipeline_done=False, rule_all=True, config_gate=None, notes=""),
    ]
    try:
        artifacts_by_manifest_output_type(dup_mot_arts)
        _check("27-by_mot_duplicate_raises", False,
               "expected ValueError but no exception raised")
    except ValueError as e:
        _check("27-by_mot_duplicate_raises",
               "duplicate" in str(e).lower(),
               f"got ValueError: {e}")

    # 28. filter_artifacts(assay_gate="mnase") returns 13 entries
    mnase_filtered = filter_artifacts(artifacts, assay_gate="mnase")
    _check("28-filter_mnase_13", len(mnase_filtered) == 13
           and all(a.assay_gate == "mnase" for a in mnase_filtered),
           f"expected 13 mnase, got {len(mnase_filtered)}")

    # 29. filter_artifacts(scope="sample", pipeline_done=True)
    sample_pd = filter_artifacts(artifacts, scope="sample",
                                 pipeline_done=True)
    _check("29-filter_sample_pipeline_done",
           len(sample_pd) > 0
           and all(a.scope == "sample" for a in sample_pd)
           and all(a.pipeline_done is True for a in sample_pd),
           f"got {len(sample_pd)} entries, invariant violated")
```

- [ ] **Step 1: Run updated test**

```bash
python3 test/test_stage45_artifact_model.py
```
Expected: 29/29 PASS (22 original + 7 new = 29 checks).

---

### Task 3: Refactor `test/test_stage43_artifact_inventory.py`

- [ ] **Step 1: Add import**

```python
from lib.artifact import Artifact, load_artifacts, artifacts_by_manifest_output_type
```

- [ ] **Step 2: Replace check 6 local dict construction**

Current:
```python
    non_null_mt = {a.manifest_output_type for a in artifacts
                   if a.manifest_output_type is not None}
```

Replacement:
```python
    by_mot = artifacts_by_manifest_output_type(artifacts)
    non_null_mt = set(by_mot.keys())
```

---

### Task 4: Refactor `test/test_stage47_mnase_path_contract.py`

- [ ] **Step 1: Add imports**

```python
from lib.artifact import load_artifacts, filter_artifacts, artifacts_by_id
```

- [ ] **Step 2: Replace manual filtering and indexing**

Current:
```python
    mnase_arts = [a for a in artifacts if a.assay_gate == "mnase"]
```
Replacement:
```python
    mnase_arts = filter_artifacts(artifacts, assay_gate="mnase")
```

Current:
```python
    inv_by_id = {a.id: a for a in mnase_arts}
```
Replacement:
```python
    inv_by_id = artifacts_by_id(mnase_arts)
```

---

### Task 5: Update artifact roadmap

Add Stage 48 section between Stage 47 and Stage 49+.

---

### Task 6: Run full test suite

```bash
python3 test/test_stage45_artifact_model.py       # 29/29
python3 test/test_stage43_artifact_inventory.py   # 10/10
python3 test/test_stage47_mnase_path_contract.py  # 17/17
python3 test/test_validation_stress.py            # 40/40
python3 test/test_stage28_release_readiness.py    # 11/11
python3 test/test_no_hardcoded_paths.py           # PASS
git diff --check
```

### Task 7: Commit

```bash
git add workflow/lib/artifact.py \
  test/test_stage45_artifact_model.py \
  test/test_stage43_artifact_inventory.py \
  test/test_stage47_mnase_path_contract.py \
  docs/architecture/artifact-roadmap.md \
  docs/superpowers/specs/2026-06-07-stage48-artifact-query-helpers-design.md \
  docs/superpowers/plans/2026-06-07-stage48-artifact-query-helpers.md

git commit -m "$(cat <<'EOF'
test: add Artifact query helpers (Stage 48)

Add artifacts_by_id(), artifacts_by_manifest_output_type(), and
filter_artifacts() to workflow/lib/artifact.py. Refactor Stage 43
and Stage 47 contract tests to use them. Add 7 query-helper checks
to Stage 45 model test (22->29).
EOF
)"
```
