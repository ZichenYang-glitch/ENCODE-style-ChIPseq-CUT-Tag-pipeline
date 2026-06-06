# Stage 50: Output Contract Generation Dry-Run — Design Spec

**Date:** 2026-06-07
**Branch:** stage50-output-contract-dry-run
**Status:** implemented 2026-06-07

## Purpose

Add a dry-run test that renders output-contract rows from the artifact inventory and verifies bidirectional coverage against `docs/output-contract.md`. Proves equivalence without adopting generation.

## Scope

### In scope

- `test/test_stage50_output_contract_dry_run.py` — 8 checks.
- Parse output-contract.md "Current Output Types" markdown tables.
- Render comparable rows from Artifact inventory via `load_artifacts()` and `artifacts_by_id()`.
- Bidirectional comparison: output_type/id set equivalence AND path_template/path equivalence.
- Normalize `<exp>` → `<experiment>` placeholder in OC paths (the only documented mismatch pattern).
- Update `docs/architecture/artifact-roadmap.md`.

### Out of scope

- No writing to `docs/output-contract.md`.
- No generator CLI or Makefile/CI wiring.
- No `.smk`, `workflow/Snakefile`, `rule all` changes.
- No `make_manifest.py`, `artifact.py` changes.
- No `artifact_path()`, target helper, DAG/runtime, manifest generation changes.

## Design

### Path normalization

Inspection shows 21 experiment-level OC entries use `<exp>` as placeholder while inventory uses `<experiment>`. All 41 sample-level entries match exactly. Single normalization:

```python
def _normalize_oc_path(path):
    return path.replace("<exp>", "<experiment>")
```

### 8 checks

| # | Check | What |
|:---|:---|:---|
| 1 | `load_artifacts_succeeds` | `load_artifacts()` loads inventory |
| 2 | `oc_section_markers` | output-contract.md has required section markers |
| 3 | `oc_rows_non_empty` | Parsed output type rows > 0 from OC |
| 4 | `no_duplicate_oc_types` | No duplicate `output_type` in OC tables |
| 5 | `oc_types_in_inventory` | Every OC output_type has a matching inventory id |
| 6 | `inventory_ids_in_oc` | Every inventory artifact id has a matching OC output_type |
| 7 | `path_templates_match` | Normalized OC path == inventory path_template for all matching pairs |
| 8 | `dry_run_output_type_set_matches` | `set(inventory_ids) == set(oc_types)` |

### What stays unchanged

- `docs/output-contract.md` — read only.
- No `.smk`, `Snakefile`, `make_manifest.py`, `artifact.py`.
- Stage 43 checks 7-8 remain — this test adds path-level comparison.

### Roadmap

Add Stage 50 section. Target helpers → Stage 51+.

## Verification

```bash
python3 test/test_stage50_output_contract_dry_run.py        # 8/8 PASS
python3 test/test_stage49_manifest_artifact_contract.py     # 8/8 PASS
python3 test/test_stage43_artifact_inventory.py             # 10/10 PASS
python3 test/test_stage45_artifact_model.py                 # 29/29 PASS
python3 test/test_stage47_mnase_path_contract.py            # 17/17 PASS
python3 test/test_stage28_release_readiness.py              # 11/11 PASS
python3 test/test_no_hardcoded_paths.py                     # PASS
git diff --check                                               # clean
```
