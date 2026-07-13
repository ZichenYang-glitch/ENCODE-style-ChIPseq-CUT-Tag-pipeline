# PR132 Read-only QC Workbench Implementation Plan

> **Implementation discipline:** Follow test-driven development task by task;
> keep `RunProgressPanel` as the only canonical run poller and use only generated
> QC operations for QC server state.

**Goal:** Add a URL-driven, responsive QC workbench which truthfully presents
persisted PR131 metrics without changing backend or scientific behavior.

**Architecture:** The route owns URL state, `RunProgressPanel` owns the one run
snapshot and polling decision, and a focused `QcWorkbench` uses Orval plus
TanStack Query infinite pagination. Pure helpers validate outcomes and cursors;
a separate list component owns responsive presentation.

**Tech stack:** React 18, React Router 6, TanStack Query v5, Orval, TypeScript,
Tailwind CSS, lucide-react, Vitest, Testing Library, Playwright.

---

### Task 1: Pure QC outcome and pagination invariants

**Files:**
- Create: `frontend/src/features/run-qc/qcState.ts`
- Create: `frontend/src/features/run-qc/qcState.test.ts`

- [ ] Write failing tests for latest indexed/failed/invalidated outcome,
  malformed count/reason, truncated history, durable metric IDs, duplicate page
  rows, valid cursor advancement, and invalid/repeated cursor detection.
- [ ] Run the focused test and confirm RED because the module is absent.
- [ ] Implement pure helpers without I/O or React state.
- [ ] Run the focused test and confirm GREEN.

### Task 2: Responsive, exact-value metric list

**Files:**
- Create: `frontend/src/features/run-qc/QcMetricList.tsx`
- Create: `frontend/src/features/run-qc/QcMetricList.test.tsx`

- [ ] Write failing tests for exact decimal text, every required field,
  null-safe labels, backend `qc_flag`, source-artifact control, Load more
  pending state, table semantics, mobile entry semantics, and long accessible
  text.
- [ ] Implement the compact semantic table and mobile vertical entries using
  generated `QcMetricResponse` only.
- [ ] Add lucide source navigation with tooltip and accessible label; make all
  identifiers and values overflow-safe.
- [ ] Run list tests and confirm GREEN.

### Task 3: Query orchestration and honest states

**Files:**
- Create: `frontend/src/features/run-qc/QcWorkbench.tsx`
- Create: `frontend/src/features/run-qc/QcWorkbench.test.tsx`

- [ ] Mock only generated `listRunQcMetrics`. Add failing tests for active run,
  pending/invalidated, unconfirmed, failed, loading skeleton, confirmed empty,
  nonzero-empty inconsistency, two-page load, invalid/repeated cursor,
  cursor-not-found reset, initial transport error, stale cached data, and retry.
- [ ] Implement a run-scoped `useInfiniteQuery`, page size 50, explicit Load
  more, strict success-envelope validation, safe local errors, and first-page
  reset using the existing QueryClient.
- [ ] Preserve loaded rows on refetch/next-page errors and never render raw
  exception or technical issue data.
- [ ] Run focused workbench tests and confirm GREEN.

### Task 4: URL-driven QC tab and one canonical poller

**Files:**
- Modify: `frontend/src/routes/runs/detail.tsx`
- Modify: `frontend/src/features/run-progress/RunProgressPanel.tsx`
- Modify: `frontend/src/features/run-progress/RunProgressPanel.test.tsx`
- Create: `frontend/src/routes/__tests__/qc-workbench-route.test.tsx`

- [ ] Add failing route tests for `?view=qc`, refresh/history restoration,
  unknown-view fallback, switching views, and source artifact navigation to
  `?view=artifacts&artifact=<id>`.
- [ ] Add failing polling tests proving succeeded QC pending uses the existing
  run query and indexed/failed/unconfirmed outcomes stop it.
- [ ] Extend the view union, tabs, keyboard behavior, and route callbacks. Clear
  stale artifact selection when entering QC/Activity.
- [ ] Derive the QC outcome from the existing event snapshot and render
  `QcWorkbench`; do not add another run query or interval.
- [ ] Run route, panel, artifact, and QC tests and confirm GREEN.

### Task 5: Real desktop/mobile Playwright gate

**Files:**
- Modify: `frontend/e2e/durable-run.spec.ts`

- [ ] Extend the existing real success flow after worker completion to open QC,
  wait for the confirmed-zero state, assert `?view=qc`, reload, and restore.
- [ ] Check 1440x900, 1024x768, 390x844, and 360x800 for visible shared
  controls, selected QC tab, no horizontal overflow, and usable status text.
- [ ] Attach full-page screenshots for desktop and mobile inspection.
- [ ] Run the real Redis/RQ/Snakemake Playwright launcher with zero skips and
  inspect screenshots at full resolution.
- [ ] Confirm the test did not add a QC source or modify the deterministic or
  production scientific workflow; non-empty real QC remains PR134.

### Task 6: Full regression and publication gate

**Files:**
- Review: complete `origin/main...HEAD` diff

- [ ] Run focused and all frontend Vitest tests, TypeScript typecheck, and the
  production build.
- [ ] Export OpenAPI and regenerate Orval; require zero generated drift and no
  package-lock changes.
- [ ] Run `PYTHONPATH=src python3 -m pytest`, relevant real durable worker tests,
  Ruff check, Ruff format check, and `git diff --check`.
- [ ] Confirm no Redis, worker, horse, Snakemake, Vite, or helper remains.
- [ ] Commit and push without attribution, create Draft PR132, wait for every
  required GitHub check, and compare local/remote/PR HEAD.
- [ ] Request one independent read-only merge-gate review. Fix all blocking or
  important findings, rerun affected/full gates, and request re-review.
- [ ] If and only if CI is green, review is clear, scope/worktree are clean, and
  base is latest main, mark Ready and squash merge with
  `--match-head-commit`; fetch and verify the merge commit before PR133.
