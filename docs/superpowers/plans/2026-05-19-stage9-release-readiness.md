# Stage 9: Release Readiness — Implementation Plan

> **Goal:** Polish README Limitations and KNOWN_ISSUES for release readiness. Documentation only — no code changes.

**Architecture:** Two file edits. README: restructure Limitations + add sample-sheet note + release checklist. KNOWN_ISSUES: replace Current Status + refresh Scope Gap + fix High Priority #3.

**Files:** `README.md`, `KNOWN_ISSUES.md` — no workflow/config/schema/script/test changes.

---

### Task 1: README.md — Limitations restructure

- Replace the flat Limitations section (~lines 335-352) with the grouped version from the spec (Section 2.1).
- The replacement text is in `docs/superpowers/specs/2026-05-19-stage9-release-readiness-design.md` Section 2.1.
- Verify: `grep "### Assay support" README.md` returns a match.
- Verify: `grep "exactly 2 treatment" README.md` returns a match.

### Task 2: README.md — Sample-sheet note + Release checklist

- After the sample sheet example table (~line 142), append the compact-replicate note from spec Section 2.2.
- Under Developer Notes, after the Local execution section, append the Release checklist from spec Section 2.3.
- Verify: `grep "compact replicate" README.md` returns a match.
- Verify: `grep "Release checklist" README.md` returns a match.

### Task 3: KNOWN_ISSUES.md — Current Status

- Replace the existing `## Current Status` section (lines 6-29) with the concise snapshot from spec Section 3.2.
- Verify the old verbose status lines are gone.

### Task 4: KNOWN_ISSUES.md — Scope Gap

- Replace lines 31-48 with the refreshed Scope Gap from spec Section 3.1.
- Verify: `grep "ENCODE-inspired" KNOWN_ISSUES.md` returns a match.
- Verify: `grep "Replicate-aware" KNOWN_ISSUES.md` returns a match (completed, not missing).

### Task 5: KNOWN_ISSUES.md — High Priority #3

- Replace lines 245-248 with the fixed version from spec Section 3.3.
- Verify: `grep "no-control real execution" KNOWN_ISSUES.md` returns a match.

### Task 6: Verification

```bash
python3 test/test_validation_stress.py
python3 test/test_stage8_smoke_profiles.py
git status --short --untracked-files=all
```

### Task 7: Commit
