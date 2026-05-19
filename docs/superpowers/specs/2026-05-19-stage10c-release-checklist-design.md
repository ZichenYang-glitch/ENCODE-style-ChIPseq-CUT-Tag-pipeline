# Stage 10c: Release Checklist Execution for v0.1.0-beta — Design Spec

**Date:** 2026-05-19
**Status:** design review → immediate execution
**Scope:** one new file — run local checks and record results in `docs/release-checks/v0.1.0-beta.md`
**Excluded:** workflow/script/schema/config/test/env changes, README/CHANGELOG updates, release tagging

---

## 1. Goals

1. Create `docs/release-checks/v0.1.0-beta.md` with actual check results (not an empty checklist).
2. Run the six local checks from the README release checklist and record PASS/FAIL/SKIP, exit code, notes, date/time.
3. Check CI status and record it.
4. No pipeline, config, or doc changes beyond this one file.

---

## 2. File to create

`docs/release-checks/v0.1.0-beta.md`

### Structure

```markdown
# v0.1.0-beta Release Verification

**Date:** 2026-05-19
**Time:** <actual execution time>
**Checklist reference:** README > Developer Notes > Release checklist
**Overall result:** <PASS if all local checks pass>

## Local checks

### 1. Validate default config and stress-test validation

- **`python3 scripts/validate_samples.py --config config/config.yaml`**
  - Result: PASS/FAIL
  - Exit code: N

- **`python3 test/test_validation_stress.py`**
  - Result: PASS/FAIL
  - Exit code: N
  - Notes: <e.g. 15/15>

### 2. Default DAG check

- **`snakemake -s workflow/Snakefile --configfile config/config.yaml -n --quiet`**
  - Result: PASS/FAIL
  - Exit code: N

### 3. Dry-run smoke profiles

- **`python3 test/test_stage8_smoke_profiles.py`**
  - Result: PASS/FAIL
  - Exit code: N
  - Notes: <e.g. 7/7>

### 4. Tiny real execution

- **`python3 test/test_stage8b_tiny_execution.py`**
  - Result: PASS / FAIL / SKIP
  - Exit code: <0/1/2>
  - Notes: <SKIP if exit 2 with missing-tool notes>

### 5. Repo hygiene

- **`git status --short --untracked-files=all`**
  - Result: PASS/FAIL
  - Notes: no `results/`, `.snakemake/`, `*.fq`, `*.fq.gz`, `*.bam`, `*.bai`, `*.bw`

### 6. CI status

- **GitHub Actions fast-checks on main:**
  - Status: <passing / failing / not checked>
  - Link: <run URL or "pending manual check">
  - Date checked: <date>

- **Manual real-execution (workflow_dispatch):**
  - Status: Not run for v0.1.0-beta local verification; manual workflow_dispatch remains optional.

## Results summary

| # | Check | Result | Exit code | Notes |
|---|-------|--------|-----------|-------|
| 1a | validate_samples config | | | |
| 1b | validation_stress | | | |
| 2 | DAG check | | | |
| 3 | smoke profiles | | | |
| 4 | tiny execution | | | |
| 5 | repo hygiene | | | |
| 6a | CI fast-checks | | | |
| 6b | CI real-execution | Not run | — | optional manual dispatch |
```

---

## 3. Execution rules

### Exit code handling

- **validate_samples:** 0 = PASS, non-zero = FAIL with stderr notes.
- **validation_stress:** 0 = PASS, non-zero = FAIL with failure count.
- **DAG check:** 0 = PASS (DAG resolved), non-zero = FAIL with Snakemake error.
- **smoke profiles:** 0 = PASS (7/7), non-zero = FAIL with profile names.
- **tiny execution:** 0 = PASS, 1 = FAIL (tool/runtime error), 2 = SKIP (missing tools, environment not available).
- **repo hygiene:** 0 = PASS (no unexpected artifacts), non-zero lines = FAIL with listed files.

### Recording format

Each check records:
- **Result:** one of `PASS`, `FAIL`, `SKIP`, or `Not run`.
- **Exit code:** the numeric exit code (or `—` for non-shell checks like CI).
- **Notes:** one-line context — test counts, SKIP reason, CI run link, unexpected artifact names.

### CI checks

- **fast-checks:** Check via `gh run list --workflow=ci.yml --branch main --limit 1 --json status,conclusion,databaseId,url` or GitHub Actions UI. Record status, conclusion, and URL.
- **real-execution:** Not run locally. Record as "Not run for v0.1.0-beta local verification; manual workflow_dispatch remains optional."

---

## 4. Out-of-scope

- No workflow, rule, script, schema, config, test, or env changes.
- No README, CHANGELOG, or KNOWN_ISSUES changes.
- No release tagging.
- No new pipeline features.

---

## 5. Self-review

- **No feature creep:** Exactly one new file under `docs/release-checks/`. All checks mirror the README release checklist.
- **Commands match README:** All six commands are copied verbatim from the README release checklist section.
- **No workflow changes:** None of the files in workflow/, scripts/, config/, test/ are touched.
- **Clear status labels:** PASS/FAIL/SKIP/Not run — each check has an unambiguous result label.
- **Tiny execution exit 2 handling:** Explicitly mapped to SKIP with missing-tool notes, not PASS.
- **CI real-execution not deferred to Stage 10d:** Recorded as "Not run" with the approved wording, not deferred.
