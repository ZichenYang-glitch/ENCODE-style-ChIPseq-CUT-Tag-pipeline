# Stage 52: Release Readiness Sweep — Design Spec

**Date:** 2026-06-07
**Branch:** stage52-release-readiness-sweep
**Status:** implemented 2026-06-07

## Purpose

Sweep documentation and tests for stale statements after Stages 40-51. Confirm internal consistency. Not a feature stage — not artifact runtime adoption.

## Findings Summary

### High severity (7 items)

| # | File | Issue |
|:---|:---|:---|
| H1 | `README.md` L18-20 | Pre-release version description stale ("12 stages", v0.2.0-rc1) |
| H2 | `README.md` L467-484 | CI test suite list missing Stages 39-50 tests |
| H3 | `CHANGELOG.md` | Missing [Unreleased] entries for Stages 41-51 |
| H4 | `KNOWN_ISSUES.md` L8-17 | Current Status stops at Stage 27c |
| H5 | `artifact-roadmap.md` L11-50 | "Current Baseline (v0.2 / Stage 40)" section describes pre-extraction state |
| H6 | `artifact-roadmap.md` | Stages 41-44 missing "implemented" status markers |
| H7 | `configuration.md` L121-128 | MNase section documents only legacy `mono_range` |

### Medium severity (4 items)

| # | File | Issue |
|:---|:---|:---|
| M1 | `README.md` L32-34 | MNase feature list omits sub/di BAMs |
| M2 | `README.md` L426-432 | Repository layout omits metadata/targets/paths.smk + workflow/lib/ |
| M3 | `assay-policy.md` L152,247 | Two references to Stage 40 as "deferred" |
| M4 | `KNOWN_ISSUES.md` L382 | Stage 37 "In progress" → should be "Completed" |

### Low severity (3 items)

| # | File | Issue |
|:---|:---|:---|
| L1 | `README.md` L288-333 | Output tree missing `03_fragments/` |
| L2 | `output-contract.md` L129-143 | Gating Conditions missing pooled MNase entries |
| L3 | `artifact-roadmap.md` L251,296-302 | Stage gate numbers stale (48/50+ → 53+) |

## Fixes Applied

### README.md
- Update pre-release version description
- Expand MNase feature list
- Update output tree with `03_fragments/`
- Update repository layout with new files
- Update CI test suite list

### CHANGELOG.md
- Add [Unreleased] entries for Stages 41-51

### KNOWN_ISSUES.md
- Update Current Status
- Update deferred list
- Fix Stage 37 status

### artifact-roadmap.md
- Update Current Baseline section
- Add "implemented" markers to Stages 41-44
- Fix stage gate numbers

### configuration.md
- Expand MNase config section with fragments/dyad_range/callers

### assay-policy.md
- Fix "deferred to Stage 40" references

### output-contract.md
- Add pooled MNase gating conditions

## Verification

```bash
python3 test/test_stage28_release_readiness.py
python3 test/test_no_hardcoded_paths.py
python3 test/test_stage43_artifact_inventory.py
python3 test/test_stage45_artifact_model.py
python3 test/test_stage47_mnase_path_contract.py
python3 test/test_stage49_manifest_artifact_contract.py
python3 test/test_stage50_output_contract_dry_run.py
git diff --check
```
