# Stage 56: Artifact Runtime Boundary ADR — Design Spec

**Date:** 2026-06-19
**Status:** Implemented
**Author:** YangZiChen-glitch (Kaslana)
**Scope:** Documentation-only — formalize the Artifact/Snakemake boundary

## 1. Purpose

Formally adopt the Artifact Runtime Boundary ADR and update the artifact
roadmap to reflect the corrected direction. No code changes, no Snakemake
changes, no test changes beyond documentation consistency.

## 2. Deliverables

### A. Artifact Runtime Boundary ADR

**File:** `docs/superpowers/specs/2026-06-19-artifact-runtime-boundary-adr.md`

Already written and locally staged. Defines:
- Artifact = product/result language (output catalog, manifest schema, API schema seed)
- Snakemake = workflow execution engine (rules, DAG, wildcards, dependencies)
- 5-level taxonomy (Level 1 catalog through Level 4 avoidance)
- CLI-only: maintain catalog + contract tests, no runtime adoption
- FastAPI: schema seed + thin `resolve_artifact_path()` outside Snakemake

### B. Artifact Roadmap Update

**File:** `docs/architecture/artifact-roadmap.md`

Update the "Staged Path" section to reference the ADR and clarify the boundary.
Specifically:
- Add a reference to the ADR in the introduction
- Add explicit boundary statement: Artifact = result semantics / platform schema; Snakemake = workflow execution engine
- Clarify Level 2.5 resolver is allowed only for future backend, not for .smk runtime
- No runtime adoption now

## 3. Files

| File | Action |
|------|--------|
| `docs/superpowers/specs/2026-06-19-artifact-runtime-boundary-adr.md` | Already exists — commit as-is |
| `docs/superpowers/specs/2026-06-19-stage56-artifact-boundary-adr.md` | Create (this file) |
| `docs/superpowers/plans/2026-06-19-stage56-artifact-boundary-adr.md` | Create |
| `docs/architecture/artifact-roadmap.md` | Update — reference ADR, clarify boundary |

## 4. Non-goals

- No Snakemake, rule, config, or validate_samples.py changes
- No new tests (existing release/doc tests are sufficient)
- No Co-Authored-By

## 5. Verification

```bash
python3 test/test_stage28_release_readiness.py
python3 test/test_no_hardcoded_paths.py
git diff --check
```
