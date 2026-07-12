# PR126 Real Run Frontend Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Connect the durable local execution API to the existing React UI, prove success and truthful cancellation in a real browser, and provide one reliable local launcher.

**Architecture:** Keep generated Orval operations behind the existing adapter and use TanStack Query as the only browser server-state layer. A run receives a durable URL immediately after creation; the detail route owns preflight, start, cancel, bounded polling, structured issues, events, and logs. A Python supervisor starts Redis, FastAPI, one RQ worker, and Vite with one SQLite/workspace configuration; Playwright drives that same real stack with a deterministic tiny workflow.

**Tech Stack:** React 18, React Router 6, TanStack Query 5, Orval, Vitest, Playwright, FastAPI, RQ/Redis, SQLite, Snakemake.

---

### Task 1: Generated-client adapter and safe API failures

**Files:**
- Modify: `frontend/src/api/runClient.ts`
- Modify: `frontend/src/api/generated-client-adapters.ts`
- Modify: `frontend/src/api/fetcher.ts`
- Test: `frontend/src/api/runClient.test.ts`
- Test: `frontend/src/api/generated-client-adapters.test.ts`
- Test: `frontend/src/api/fetcher.test.ts`

- [ ] Add failing tests that require `startRun(runId)` to call the generated `/start` operation and require non-2xx responses to retain only public Issue fields (`code`, `message`, `severity`, `path`, `source`, `hint`) while discarding `technical_message` and arbitrary `context`.
- [ ] Run `npm test -- --run src/api/runClient.test.ts src/api/generated-client-adapters.test.ts src/api/fetcher.test.ts`; expect failures for the missing method and issue envelope.
- [ ] Add `startRun` to `RunApiClient`, its HTTP/stub implementations, and the generated adapter. Extend `ApiError` with a typed public issue array populated by the shared fetcher without retaining unsafe fields.
- [ ] Re-run the focused command; expect all tests to pass.

### Task 2: Truthful run detail state and actions

**Files:**
- Modify: `frontend/src/features/run-progress/RunProgressPanel.tsx`
- Modify: `frontend/src/features/run-progress/RunStatusBadge.tsx`
- Modify: `frontend/src/features/run-progress/RunEventFeed.tsx`
- Create: `frontend/src/features/run-progress/RunIssuePanel.tsx`
- Test: `frontend/src/features/run-progress/RunProgressPanel.test.tsx`
- Create: `frontend/src/features/run-progress/RunIssuePanel.test.tsx`

- [ ] Add failing component tests for PLANNED start, RUNNING cancellation remaining RUNNING after 202, cancellation-request restoration from events, retryable cancellation errors, terminal-only polling stop, unknown-status safe fallback, structured failure rendering, and bounded cursor collection.
- [ ] Run `npm test -- --run src/features/run-progress`; expect the new behavioral assertions to fail.
- [ ] Replace local async action bookkeeping with TanStack mutations and invalidation. Poll only known active states for a bounded window; stop on PLANNED, terminal, unknown, or the bound. Fetch at most five pages per events/log stream, rejecting repeated cursors.
- [ ] Render only public Issue fields. Replace arbitrary event-context JSON with a lifecycle-field allowlist. Make unknown status neutral and disable start/cancel for unknown states.
- [ ] Re-run `npm test -- --run src/features/run-progress`; expect all tests to pass.

### Task 3: Durable create-to-detail navigation

**Files:**
- Modify: `frontend/src/features/workflow-detail/WorkflowDetail.tsx`
- Modify: `frontend/src/routes/workflows/detail.tsx`
- Modify: `frontend/src/routes/runs/detail.tsx`
- Test: `frontend/src/routes/__tests__/real-preflight-path.test.tsx`
- Test: `frontend/src/routes/__tests__/router.test.tsx`

- [ ] Add a failing route test proving create success immediately navigates to `/runs/:runId`, preflight is invoked from the destination route, a preflight failure leaves the durable URL visible, and browser reload can recover from GET without repeating preflight.
- [ ] Run `npm test -- --run src/routes/__tests__/real-preflight-path.test.tsx src/routes/__tests__/router.test.tsx`; expect ordering/navigation failures.
- [ ] Seed the run-progress query after create and navigate with one-use route state requesting preflight. Consume that state using `replace`; offer an explicit `Run preflight` action for a reloaded CREATED run.
- [ ] Re-run the focused route tests; expect all tests to pass.

### Task 4: Real local runtime supervisor

**Files:**
- Create: `scripts/run_local_platform.py`
- Create: `test/browser/platform_runtime.py`
- Create: `test/browser/test_platform_runtime.py`
- Modify: `frontend/vite.config.ts`
- Modify: `.gitignore`
- Modify: `docs/development/local-platform-runtime.md`

- [ ] Add failing Python tests for command construction, shared environment values, readiness timeout, and reverse-order process-group cleanup.
- [ ] Run `PYTHONPATH=src python3 -m pytest test/browser/test_platform_runtime.py -q`; expect import/behavior failures.
- [ ] Implement a foreground supervisor that optionally starts `redis-server`, then starts FastAPI, one RQ worker, and Vite in separate process groups; share explicit DB, queue, workspace, Redis URL, API proxy target, and runtime directory; terminate then kill every group on exit.
- [ ] Document `python3 scripts/run_local_platform.py` with prerequisites, ports, data path, readiness, and Ctrl-C cleanup. Ignore only local runtime and Playwright artifacts.
- [ ] Re-run the focused Python tests; expect all tests to pass.

### Task 5: Real Playwright success and cancellation gates

**Files:**
- Modify: `frontend/package.json`
- Modify mechanically: `frontend/package-lock.json`
- Create: `frontend/playwright.config.ts`
- Create: `frontend/e2e/durable-run.spec.ts`
- Modify: `test/browser/platform_runtime.py`

- [ ] Install official `@playwright/test` and add `test:e2e`; configure one desktop success project and one mobile cancellation project with one worker and retained failure artifacts.
- [ ] Build a deterministic isolated workflow fixture with a quick success config and a long-running cancellation config. Persist only browser-safe fixture inputs in a runtime manifest.
- [ ] Write the desktop test: validate, create, verify durable URL and PLANNED, start, observe QUEUED/RUNNING, SUCCEEDED, stdout marker, reload, and re-check persisted state/logs.
- [ ] Write the mobile test: start the long run, capture a real cancel response with HTTP 202 and `running`, verify “Cancellation requested” without a fake state, await backend-confirmed CANCELLED, and assert controls remain operable without horizontal overflow.
- [ ] Run `npm run test:e2e`; expect both projects to pass without skips and cleanup to leave no Redis/worker/Snakemake descendants.

### Task 6: CI, full validation, and Draft PR

**Files:**
- Modify: `.github/workflows/ci.yml`
- Modify if generated drift exists: `frontend/openapi.json`
- Modify mechanically if generated drift exists: `frontend/src/api/generated/**`

- [ ] Add a required browser E2E job with a Redis service, the existing `ci-fast` environment, Node 20, Chromium installation, `npm ci`, and `npm run test:e2e`; upload Playwright artifacts only on failure.
- [ ] Run focused backend/frontend tests, then `PYTHONPATH=src python3 -m pytest`, real durable worker tests, `npm test -- --run`, `npm run typecheck`, `npm run build`, OpenAPI export/generation drift, Ruff check/format check, and `git diff --check`; fix every failure in scope.
- [ ] Review the complete diff for lifecycle truthfulness, cancellation races, issue redaction, polling bounds, child cleanup, responsive behavior, generated drift, and unrelated files.
- [ ] Commit without attribution, push `agent/pr126-real-run-frontend`, create Draft PR #126, and wait for every GitHub check to finish green.
- [ ] Ask two read-only reviewers to perform the final durable-state/cancellation and frontend/E2E/demo merge-gate reviews; fix every blocking finding, re-run validation, and leave the PR unmerged.

## Self-review

- Spec coverage: Tasks 1–3 cover the real UI lifecycle, truthful cancellation, reload recovery, structured issues, generated adapter boundary, bounded polling, and responsive controls. Tasks 4–5 cover the one-command runtime and non-mocked success/cancel browser gates. Task 6 covers generated drift, required CI, full validation, independent review, and Draft-only delivery.
- Scope: No artifact/QC, auth, multi-tenant, HPC, second adapter, reconciler, immutable bundle, or Agent write capability is introduced.
- Type consistency: `startRun`, public Issue fields, route-state preflight, `RunSnapshot`, and runtime manifest names are used consistently across their defining and consuming tasks.
