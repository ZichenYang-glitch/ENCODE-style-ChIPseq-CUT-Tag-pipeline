# Stage 51: Artifact Adoption Decision Record — Design Spec

**Date:** 2026-06-07
**Branch:** stage51-artifact-adoption-decision
**Status:** implemented 2026-06-07

## Purpose

Create a decision record evaluating whether the project should proceed from artifact contract tests (Stages 41-50) into runtime-adjacent artifact adoption. Documentation-only — no code changes.

## Scope

### In scope

- `docs/architecture/artifact-adoption-decision.md` — full decision record.
- `docs/architecture/artifact-roadmap.md` — Stage 51 section.
- Spec/plan docs under `docs/superpowers/`.

### Out of scope

- No code, test, `.smk`, `Snakefile`, `make_manifest.py`, `artifact.py` changes.
- No DAG/runtime behavior changes.

## Design

### Decision record structure

1. Context: summary of Stages 41-50 achievements
2. Decision: pause runtime artifact adoption
3. Candidate options A-F with value/risk/preconditions/blast-radius/stage
4. Recommendation: pause, return to release hardening/science

### Candidate options

| Option | Risk | Status |
|:---|:---|:---|
| A. Output-contract generation adoption | Low | Possible next if generation is worth it |
| B. Manifest generation from inventory | Medium | Requires import-path decision |
| C. Broader paths.smk migration | Medium | Assay-by-assay only |
| D. Artifact-assisted target helpers | High | Gated on demonstrated maintenance pain |
| E. Global target resolver | Very high | Explicit team decision gate |
| F. Pause — return to release hardening/science | Low | **Recommended** |

### Recommendation

Pause runtime artifact adoption. Return to release hardening and scientific features. The artifact contract test infrastructure remains in place and self-sustaining.

## Verification

```bash
python3 test/test_stage28_release_readiness.py   # 11/11 PASS
python3 test/test_no_hardcoded_paths.py          # PASS
git diff --check                                    # clean
```
