# Stage 57 Implementation Plan: Shell Safety Hardening

**Date:** 2026-06-19
**Status:** Implemented
**Scope:** P0 safety — add `set -e -o pipefail` to all Snakemake shell blocks

## Files changed

| File | Change |
|------|--------|
| `workflow/rules/common.smk` | +10 guards |
| `workflow/rules/qc.smk` | +6 guards |
| `workflow/rules/report.smk` | +2 guards |
| `test/test_stage57_shell_safety.py` | Create — enforcement test |
| `docs/superpowers/specs/2026-06-19-stage57-shell-safety-design.md` | Create |
| `docs/superpowers/plans/2026-06-19-stage57-shell-safety.md` | Create |

## Steps

1. Add `set -e -o pipefail` to 18 shell blocks
2. Write enforcement test parsing all .smk shell blocks
3. Verify: test passes, release tests pass, diff clean

## Non-goals

- No rule logic, path, target, config, or manifest changes
- No `set -euo pipefail`
- No Co-Authored-By
