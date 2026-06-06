# Stage 49: Manifest Artifact Contract Tests — Design Spec

**Date:** 2026-06-07
**Branch:** stage49-manifest-artifact-coverage
**Status:** implemented 2026-06-07

## Purpose

Add a dedicated test that verifies manifest output_type vocabulary matches the artifact inventory, using Stage 48 query helpers and Stage 45 Artifact model. Pure test-layer addition.

## Scope

### In scope

- Create `test/test_stage49_manifest_artifact_contract.py` (8 checks).
- Inline AST extraction from `scripts/make_manifest.py` (read-only).
- Compare manifest vocabulary against `artifacts_by_manifest_output_type()` keys.
- Validate tracked Artifacts: producing_rule non-empty, path_template starts with results/.
- Update `docs/architecture/artifact-roadmap.md`.

### Out of scope

- No `.smk`, `workflow/Snakefile`, `rule all` changes.
- No `scripts/make_manifest.py` edits.
- No `workflow/lib/artifact.py` changes.
- No `artifact_path()`, target helper, DAG/runtime, or manifest generation changes.
- Stage 43 checks 5-6 are preserved — this test adds orthogonal coverage.

## Design

### 8 checks

| # | Check | What |
|:---|:---|:---|
| 1 | `load_artifacts_succeeds` | `load_artifacts()` loads inventory |
| 2 | `manifest_types_from_ast` | AST extracts non-empty set from `make_manifest.py` |
| 3 | `no_duplicate_manifest_output_type` | `artifacts_by_manifest_output_type()` doesn't raise ValueError |
| 4 | `manifest_vs_inventory_bidirectional` | `manifest_types == set(by_mot.keys())` |
| 5 | `tracked_mapping_values_match_keys` | `all(a.manifest_output_type == key for key, a in by_mot.items())` |
| 6 | `producing_rule_non_empty` | Every tracked Artifact has non-empty `producing_rule` |
| 7 | `path_starts_results` | Every tracked Artifact's `path_template` starts with `results/` |
| 8 | `tracked_artifacts_non_empty` | `len(by_mot) > 0` |

### Key design decisions

- Check 3 (`no_duplicate_manifest_output_type`) is `try/except ValueError` — if duplicate exists, record failure and exit early because later checks depend on the manifest-output-type mapping.
- No standalone `mot_not_null` — `artifacts_by_manifest_output_type()` already excludes null values.
- Stage 43 checks 5-6 remain unchanged.
- AST extraction is inlined in the new test.

### Roadmap

Add Stage 49 section, pull Stage 50+ for target helpers.

## Verification

```bash
python3 test/test_stage49_manifest_artifact_contract.py   # 8/8 PASS
python3 test/test_stage43_artifact_inventory.py            # 10/10 PASS
python3 test/test_stage45_artifact_model.py                # 29/29 PASS
python3 test/test_stage47_mnase_path_contract.py           # 17/17 PASS
python3 test/test_stage25_manifest_stress.py               # 18/18 PASS
python3 test/test_stage28_release_readiness.py             # 11/11 PASS
python3 test/test_no_hardcoded_paths.py                    # PASS
git diff --check                                              # clean
```
