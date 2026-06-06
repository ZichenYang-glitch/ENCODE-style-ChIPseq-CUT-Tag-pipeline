# Stage 52: Release Readiness Sweep — Implementation Plan

> **For agentic workers:** Documentation-only. No code changes.

**Goal:** Sweep documentation for stale statements after Stages 40-51. Confirm internal consistency.

---

### Task 1: Fix README.md

- Update pre-release version description to reflect post-rc1 changes
- Expand repo layout with `workflow/lib/` and new `.smk` files
- Update test suite description

### Task 2: Fix CHANGELOG.md

- Add [Unreleased] entries for Stages 41-51 and Stage 52

### Task 3: Fix KNOWN_ISSUES.md

- Update Current Status through Stage 51
- Fix Stage 37 status from "In progress" to completed
- Fix Stage 27a wording (planned, not implemented)

### Task 4: Fix artifact-roadmap.md

- Update Current Baseline section
- Add "implemented" markers to Stages 41-44
- Fix stage gate numbers (48/50+ → 53+)
- Add Stage 52 completion note

### Task 5: Fix configuration.md

- Expand MNase config section with fragments/dyad_range/callers

### Task 6: Fix assay-policy.md

- Fix "deferred to Stage 40" → "implemented in Stage 40"

### Task 7: Fix output-contract.md

- Add pooled MNase gating conditions

### Task 8: Run acceptance checks

```bash
python3 test/test_stage28_release_readiness.py   # 11/11
python3 test/test_no_hardcoded_paths.py          # PASS
git diff --check
```

### Task 9: Commit

```bash
git commit -m "docs: sweep release readiness after artifact stages (Stage 52)"
```
