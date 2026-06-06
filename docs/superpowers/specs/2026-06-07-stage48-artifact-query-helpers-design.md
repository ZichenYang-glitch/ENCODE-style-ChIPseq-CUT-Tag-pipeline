# Stage 48: Artifact Query Helpers — Design Spec

**Date:** 2026-06-07
**Branch:** stage48-artifact-query-helpers
**Status:** implemented 2026-06-07

## Purpose

Add 3 tiny query/index helpers to `workflow/lib/artifact.py` and refactor existing contract tests to use them. Pure docs/tests/model — no DAG/runtime changes.

## Scope

### In scope

1. `workflow/lib/artifact.py` — add `artifacts_by_id()`, `artifacts_by_manifest_output_type()`, `filter_artifacts()`.
2. `test/test_stage45_artifact_model.py` — add 6 checks (23-28) for the new helpers.
3. `test/test_stage43_artifact_inventory.py` — replace local set/dict construction with helpers where appropriate.
4. `test/test_stage47_mnase_path_contract.py` — use `filter_artifacts()` and `artifacts_by_id()`.
5. `docs/architecture/artifact-roadmap.md` — Stage 48 section.

### Out of scope

- No `.smk`, `workflow/Snakefile`, `rule all` changes.
- No `make_manifest.py` changes.
- No `artifact_path()` or path rendering.
- No DAG/runtime behavior changes.

## Design

### Three helpers

```python
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

def filter_artifacts(artifacts, *, assay_gate=None, scope=None,
                     level=None, pipeline_done=None, rule_all=None):
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
```

### Test additions (Stage 45, checks 23-28)

| # | Check | What |
|:---|:---|:---|
| 23 | `artifacts_by_id_62` | `artifacts_by_id(artifacts)` returns dict with 62 keys |
| 24 | `by_id_duplicate_raises` | Duplicate id raises ValueError |
| 25 | `by_manifest_output_type` | Keys = all non-null manifest_output_type values; each `v.manifest_output_type == key` |
| 26 | `by_mot_duplicate_raises` | Duplicate manifest_output_type raises ValueError |
| 27 | `filter_mnase_13` | `filter_artifacts(artifacts, assay_gate="mnase")` returns 13 entries |
| 28 | `filter_sample_pipeline_done` | `filter_artifacts(artifacts, scope="sample", pipeline_done=True)` returns non-empty list where all have scope=="sample" and pipeline_done==True |

### Test refactors

**Stage 43** — Replace local dict construction in check 6:
```python
non_null_mt = {a.manifest_output_type for a in artifacts
               if a.manifest_output_type is not None}
```
→
```python
by_mot = artifacts_by_manifest_output_type(artifacts)
non_null_mt = set(by_mot.keys())
```

**Stage 47** — Replace manual filtering and indexing:
```python
mnase_arts = [a for a in artifacts if a.assay_gate == "mnase"]
inv_by_id = {a.id: a for a in mnase_arts}
```
→
```python
mnase_arts = filter_artifacts(artifacts, assay_gate="mnase")
inv_by_id = artifacts_by_id(mnase_arts)
```

### Roadmap

Add Stage 48 section, push target-helper adoption to Stage 49+.

## Verification

```bash
python3 test/test_stage45_artifact_model.py       # 29/29 PASS
python3 test/test_stage43_artifact_inventory.py   # 10/10 PASS
python3 test/test_stage47_mnase_path_contract.py  # 17/17 PASS
python3 test/test_validation_stress.py            # 40/40 PASS
python3 test/test_stage28_release_readiness.py    # 11/11 PASS
python3 test/test_no_hardcoded_paths.py           # PASS
git diff --check                                    # clean
```
