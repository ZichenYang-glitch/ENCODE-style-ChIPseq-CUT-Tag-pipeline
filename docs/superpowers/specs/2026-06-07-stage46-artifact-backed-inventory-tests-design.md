# Stage 46: Artifact-Backed Inventory Tests — Design Spec

**Date:** 2026-06-07
**Branch:** stage46-artifact-backed-inventory-tests
**Status:** implemented 2026-06-07

## Purpose

Refactor `test/test_stage43_artifact_inventory.py` to consume Stage 45's Artifact model (`load_artifacts`, `Artifact`, `dataclasses.fields`). This is a test-layer adoption step — no new checks, no new logic, just a different data source.

## Scope

### In scope

- Replace manual `yaml.safe_load()` with `load_artifacts()` from `workflow/lib/artifact.py`.
- Replace raw dict access (`a["id"]`, `a.get(...)`) with Artifact attribute access (`a.id`, `a.path_template`).
- Replace the duplicate `REQUIRED_FIELDS` set with `dataclasses.fields(Artifact)`.
- Rename check 1 from "yaml_parses" to "load_artifacts_succeeds".
- Update `docs/architecture/artifact-roadmap.md` with Stage 46 completion note.

### Out of scope

- No `.smk`, `workflow/Snakefile`, `rule all` changes.
- No `scripts/make_manifest.py` changes.
- No `workflow/lib/artifact.py` changes (consumed, not modified).
- No `output-contract.md` changes.
- No new checks. No check semantics changed.
- No DAG/runtime behavior changes.

## Design

### Import addition

```python
import sys
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO_ROOT, "workflow"))
from lib.artifact import Artifact, load_artifacts
from dataclasses import fields
```

The standalone `import yaml` is removed — `load_artifacts()` handles YAML internally.

### Check-by-check changes

| # | Check name (new) | Data source | Change |
|:---|:---|:---|:---|
| 1 | `load_artifacts_succeeds` | `artifacts = load_artifacts(inv_path)` | Replaces `yaml.safe_load()` + dict check. If loader raises, test fails. |
| 1b | `has_entries` | `len(artifacts) >= 50` | Unchanged logic. |
| 2 | `unique_ids` | `a.id` | Dict access → attribute. |
| 3 | `no_workstation_paths` | `a.path_template` | Dict access → attribute. `path_template` guaranteed `str`. |
| 4 | `no_future_keywords_in_notes` | `a.notes` | Dict access → attribute. `notes` guaranteed `str`. |
| 5 | `manifest_coverage` | `a.manifest_output_type` | Dict access → attribute. Already validated as `str\|None`. |
| 6 | `no_extra_manifest_types` | Same | Same pattern. |
| 7 | `oc_types_in_inventory` | `a.id` | Dict access → attribute. |
| 8 | `inventory_ids_in_oc` | Same | Same pattern. |
| 9 | `schema_completeness` | `{f.name for f in fields(Artifact)}` | No duplicate `REQUIRED_FIELDS`. Check verifies 13 fields via `dataclasses.fields()`, then confirms all entries are `Artifact` instances. Bool/string validation stays in Stage 45 test. |

### What stays unchanged

- `_extract_manifest_output_types()` — AST extraction from `make_manifest.py`.
- `_check()` framework, `main()`, exit code.
- `BAD_KEYWORDS` list.
- Output-contract regex extraction from `docs/output-contract.md`.
- Check 1b through 8 logic — only the data access path changes.

### Roadmap update

Add one section to `docs/architecture/artifact-roadmap.md`:

```
### Stage 46: Artifact-Backed Inventory Tests (implemented 2026-06-07)

Test-layer adoption: `test_stage43_artifact_inventory.py` now uses
`load_artifacts()` and `Artifact` from `workflow/lib/artifact.py`
as its canonical inventory read path. No DAG/runtime changes.
```

## Verification

```bash
python3 test/test_stage43_artifact_inventory.py      # 10/10 PASS
python3 test/test_stage45_artifact_model.py           # 22/22 PASS
python3 test/test_validation_stress.py               # 40/40 PASS
python3 test/test_stage28_release_readiness.py        # 11/11 PASS
python3 test/test_no_hardcoded_paths.py              # PASS
git diff --check                                       # clean
```
