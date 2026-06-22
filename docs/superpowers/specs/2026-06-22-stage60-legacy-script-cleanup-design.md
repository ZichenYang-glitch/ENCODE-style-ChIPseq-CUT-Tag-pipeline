# Stage 60: Legacy Script Cleanup — Design Spec

**Date:** 2026-06-22
**Status:** Implemented
**Author:** YangZiChen-glitch (Kaslana)
**Scope:** P0 hardening — remove/deprecate legacy monolithic shell entrypoint

## 1. Purpose

The legacy `scripts/chipseq.sh` (~400-line shell script) was the pre-Snakemake
entrypoint. It is superseded by `workflow/Snakefile` + `config/config.yaml`.
Users with old notes or muscle memory may accidentally run the old script.
Stage 60 archives the historical content and replaces the script with a clear
deprecation shim.

## 2. Audit

`chipseq.sh` references found in 5 files (excluding itself):
- README.md: 1 reference (legacy section)
- KNOWN_ISSUES.md: 2 references (medium-priority tech debt)
- docs/archive/: 8 references (historical design docs)
- .claude/settings.local.json: 1 reference (permissions config)
- Zero test files execute or reference the script

## 3. Approach

**Option C: archive + shim.** Move historical full script to
`docs/archive/scripts/chipseq-legacy.sh`, replace `scripts/chipseq.sh`
with a deprecation shim that exits 1.

## 4. Files

| File | Action |
|------|--------|
| `docs/archive/scripts/chipseq-legacy.sh` | Create — historical full script |
| `scripts/chipseq.sh` | Replace with deprecation shim |
| `README.md` | Update legacy section |
| `KNOWN_ISSUES.md` | Update issues #5 and #7 |
| `test/test_stage60_legacy_script.py` | 8 tests |
| `docs/superpowers/plans/2026-06-22-stage60-legacy-script-cleanup.md` | Implementation plan |

## 5. Non-goals

- No `.smk`, `Snakefile`, config, or scientific changes
- No Co-Authored-By
