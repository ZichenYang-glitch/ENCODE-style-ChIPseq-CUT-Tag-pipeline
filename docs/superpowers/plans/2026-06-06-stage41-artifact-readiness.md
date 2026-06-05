# Stage 41: Artifact Readiness — Implementation Plan

**Date:** 2026-06-06
**Status:** implemented after Stage 40 merge
**Spec:** `docs/superpowers/specs/2026-06-06-stage41-artifact-readiness-design.md`

---

## Pre-requisite: Branch setup

```bash
git switch main
git pull --ff-only origin main
git switch -c stage41-artifact-readiness
```

---

## Tasks

### Task 1: Archive reconstruction report

- [ ] Create directory: `mkdir -p docs/archive`
- [ ] Move: `mv reconstruction-deep-research-report.md docs/archive/reconstruction-deep-research-report-shell-era.md`
- [ ] Prepend archive header (see spec Section 2A)
- [ ] Do NOT touch `reconstruction-deep-research-report.md:Zone.Identifier`

### Task 2: Create artifact roadmap

- [ ] Create directory: `mkdir -p docs/architecture`
- [ ] Write `docs/architecture/artifact-roadmap.md` with:
  - "Current baseline" section (Snakemake + target helpers + assay rule files)
  - Staged path (Stage 41 through Stage 50+)
  - Each stage: goal, scope, non-goals, trigger condition
  - Explicit hard gates (what must NOT happen yet)

### Task 3: Create developer checklist

- [ ] Create directory: `mkdir -p docs/developer`
- [ ] Write `docs/developer/add-output-or-assay-checklist.md` with 11 numbered items covering:
  1. Sample sheet schema / validation
  2. Config schema
  3. Workflow rule file
  4. Target helper wiring
  5. pipeline_done wiring
  6. Manifest rows
  7. MultiQC custom content
  8. output-contract update
  9. Smoke/dry-run tests
  10. Negative tests for assay gating
  11. Docs (assay-policy, configuration, README, CHANGELOG)

### Task 4: Output contract consistency pass

- [ ] Read `docs/output-contract.md` — extract all `output_type` values
- [ ] Read `scripts/make_manifest.py` — extract all `output_type` string literals
- [ ] Read `workflow/Snakefile` `_*_targets()` functions — identify all expanded paths
- [ ] Read `workflow/rules/report.smk` `pipeline_done` inputs — identify all named inputs
- [ ] Compare vocabularies. Report any mismatches.
- [ ] Fix `docs/output-contract.md` only if a concrete mismatch exists.
- [ ] Do NOT modify `.smk`, `Snakefile`, or manifest code.

### Task 5: Run tests

- [ ] `python3 test/test_stage28_release_readiness.py`
- [ ] `python3 test/test_no_hardcoded_paths.py`
- [ ] `git diff --check`

### Task 6: Commit

- [ ] `git add` all new/modified docs
- [ ] Commit with message: `docs: add Stage 41 artifact readiness documentation`
