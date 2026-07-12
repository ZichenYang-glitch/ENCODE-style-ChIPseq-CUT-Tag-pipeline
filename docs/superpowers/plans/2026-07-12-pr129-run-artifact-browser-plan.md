# PR129 Run Artifact Browser Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a URL-addressable, responsive Artifact Browser to the existing run workbench using persisted PR128 artifact references.

**Architecture:** Keep `RunProgressPanel` as the single canonical run polling and action controller. Add focused artifact state, list, inspector, and query components that call Orval generated artifact operations directly; the run route owns URL state.

**Tech Stack:** React 18, React Router 6, TanStack Query v5, Orval, TypeScript, Tailwind CSS, lucide-react, Vitest, Testing Library, Playwright.

---

### Task 1: URL-driven workbench tabs

**Files:**
- Modify: `frontend/src/routes/runs/detail.tsx`
- Modify: `frontend/src/features/run-progress/RunProgressPanel.tsx`
- Modify: `frontend/src/routes/__tests__/router.test.tsx`

- [ ] Add route tests starting at `/runs/stub-run-1?view=artifacts`, switching
  tabs, selecting a synthetic artifact callback, and navigating Back/Forward.
  Assert `view` and `artifact` are restored from the router location.
- [ ] Run `npm test -- --run src/routes/__tests__/router.test.tsx` and confirm
  the new tests fail because the tab controls do not exist.
- [ ] Parse URL state with `useSearchParams`. Treat `artifact` as implying the
  Artifacts view, remove both parameters for Activity, and push history for tab
  and selection changes.
- [ ] Validate decoded artifact deep links against
  `^[A-Za-z][A-Za-z0-9_.-]{0,127}$`. Add `/`, `..`, encoded separator, and
  overlong tests proving `getRunArtifact` is never called for invalid IDs.
- [ ] Add an accessible tablist to `RunProgressPanel`, with the canonical run
  header and actions above both tab panels. Keep current creation mode behavior
  unchanged when no persisted `runId` is supplied.
- [ ] Run the route and existing RunProgressPanel tests and confirm both old and
  new behavior pass.

### Task 2: Pure artifact state and formatting helpers

**Files:**
- Create: `frontend/src/features/run-artifacts/artifactState.ts`
- Create: `frontend/src/features/run-artifacts/artifactState.test.ts`

- [ ] Write failing tests for latest-outcome selection, missing/invalid indexed
  counts, duplicate artifact IDs across pages, repeated/empty cursors, IEC byte
  formatting, and invalid timestamp fallback.
- [ ] Run the focused test and confirm the helpers are missing.
- [ ] Implement `artifactExtractionOutcome(events)`, `flattenArtifactPages`,
  `safeNextArtifactCursor`, `formatBytes`, and `formatProducedTime` as pure
  functions over generated model types and existing run event types.
- [ ] Run the focused test and confirm all helper cases pass.

### Task 3: Responsive artifact list

**Files:**
- Create: `frontend/src/features/run-artifacts/ArtifactList.tsx`
- Create: `frontend/src/features/run-artifacts/ArtifactList.test.tsx`

- [ ] Add failing tests for all required columns/fields, selected state,
  accessible row selection, full long-path text, Load more pending/disabled
  behavior, and mobile entry markup.
- [ ] Implement a compact semantic table for `md` and wider and bordered mobile
  entries below `md`. Render generated `ArtifactReferenceResponse` fields
  directly, using pure byte/time formatters.
- [ ] Add break/`min-w-0` classes and `title` attributes so long IDs and paths
  remain complete and never force horizontal scrolling.
- [ ] Run list tests and verify they pass.

### Task 4: Artifact detail inspector and copy states

**Files:**
- Create: `frontend/src/features/run-artifacts/ArtifactInspector.tsx`
- Create: `frontend/src/features/run-artifacts/ArtifactInspector.test.tsx`

- [ ] Add failing tests for metadata whitelist rendering, opaque URI and
  relative path, missing optional values, stable loading/error placeholders,
  Retry, clipboard success, clipboard rejection, labels, and absence of a
  Download button.
- [ ] Implement a border-separated inspector with generated detail response
  types. Use `navigator.clipboard.writeText`, lucide `Copy`, `title`, and
  `aria-label`; announce success/failure through a polite live region.
- [ ] Run inspector tests and verify all states pass.

### Task 5: Artifact query orchestration and honest states

**Files:**
- Create: `frontend/src/features/run-artifacts/ArtifactBrowser.tsx`
- Create: `frontend/src/features/run-artifacts/ArtifactBrowser.test.tsx`
- Modify: `frontend/src/features/run-progress/RunProgressPanel.tsx`
- Modify: `frontend/src/features/run-progress/RunProgressPanel.test.tsx`

- [ ] Mock only the Orval generated `listRunArtifacts` and `getRunArtifact`
  operations. Add failing tests for pending-success, indexing, extraction
  failure, indexed empty, first-load error, retained pages after refetch error,
  two-page Load more, repeated cursor stop, deep-link detail, and Retry.
- [ ] Implement `useInfiniteQuery` with query key `['run-artifacts', runId]`,
  initial cursor `undefined`, page size 50, explicit Load more, deduplicated
  flattening, and defensive cursor termination.
- [ ] Enable artifact list/detail only after the latest outcome is indexed.
  Render redacted local messages for `ApiError` and unknown failures without
  printing exception strings or `technical_message`.
- [ ] Fail closed unless successful envelopes have `ok=true`, the requested
  `run_id`, and matching artifact identity. Require indexed `artifact_count` to
  be a non-negative safe integer. A truncated event snapshot with no visible
  outcome shows “status cannot be confirmed” and Refresh, never true empty.
- [ ] Derive one `shouldPollRunSnapshot` for lifecycle polling or an active
  succeeded/indexing Artifact view. Use it for `refetchInterval`, the existing
  15-minute timeout timer, the paused notice, and Refresh window reset; do not
  create a second run query or unbounded indexing poll.
- [ ] Label extraction failure recovery “Refresh status”; reserve “Retry” for
  list/detail API refetches because no extraction retry route exists.
- [ ] Render the list and inspector in a single responsive grid, inspector
  second in DOM order. Run focused artifact and RunProgressPanel tests.

### Task 6: Real browser artifact success and responsive evidence

**Files:**
- Modify: `test/browser/platform_runtime.py`
- Modify: `frontend/e2e/durable-run.spec.ts`

- [ ] Change only the deterministic browser test project so the success task
  writes both `result/complete.txt` and
  `results/multiqc/result_manifest.tsv`; cancellation writes neither terminal
  output before acknowledgement.
- [ ] Extend the desktop success test to open Artifacts, wait for the
  `result_manifest` row, open the inspector, assert no Download action, reload
  `?view=artifacts&artifact=...`, and verify the detail restores.
- [ ] At 1440×900, 1024×768, 390×844, and 360×800, assert no document horizontal
  overflow and that tabs, Refresh, row selection, inspector, and Copy are
  operable. Attach desktop and mobile screenshots for inspection.
- [ ] Run the real Playwright gate with Redis and the independent DurableWorker;
  confirm success/cancellation tests pass and cleanup leaves no owned process.
- [ ] Inspect attached screenshots at full resolution for overlap, clipping,
  wrapping, hierarchy, and truthful status text.

### Task 7: Regression, drift, and publication gate

**Files:**
- Review: complete `origin/main...HEAD` diff

- [ ] Run all frontend tests, TypeScript typecheck, and production build.
- [ ] Run OpenAPI export and Orval regeneration with `PYTHONPATH=../src`; require
  zero generated diff and no package-lock churn.
- [ ] Run `PYTHONPATH=src python3 -m pytest -q`, real durable Redis/RQ/Snakemake
  gates, changed-file Ruff, Ruff format check, and `git diff --check`.
- [ ] Confirm the diff contains no backend contract, artifact extraction,
  production Snakemake, QC, download, authentication, or deferred backlog work.
- [ ] Request a final independent read-only UX/code review, fix only blocking
  findings, and rerun affected plus full gates.
- [ ] Commit without AI/vendor attribution, push
  `agent/pr129-artifact-browser-ui`, create Draft PR129 against main, and wait
  for GitHub fast-checks, browser-e2e, frontend, lint, and lock-check to pass.
  Keep the PR Draft and do not merge.
