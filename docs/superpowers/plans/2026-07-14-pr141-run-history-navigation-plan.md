# PR141 Run History and Navigation Implementation Plan

> Implement one persistent run-history read path and one global navigation
> hierarchy. Preserve all lifecycle, QC, artifact, and workflow semantics.

## 1. Domain and repository contract (TDD)

- Add workflow-neutral `RunSummary` and run-history cursor/filter/page
  primitives adjacent to run lifecycle domain code.
- Add failing InMemory parity tests for descending tie order, filters, cursor
  identity/filter checks, deleted cursor, bounded pages, and corrupt mapping,
  selected-row, and sentinel identities.
- Make the shared `RunSummary` projection reject status/timestamp evidence
  drift, including terminal rows without `ended_at`, active rows with
  `ended_at`, and running/succeeded rows without `started_at`.
- Add a dedicated `RunRepository.list_run_summaries` contract while preserving
  `list_runs()` and its full-record creation-order behavior unchanged.
- Add equivalent SQLAlchemy tests proving keyset pages, tie ordering, filters,
  same-time greater/smaller-ID inserts, strict summary-only boundary lookup,
  exact selected columns, and an SQL-level limit.

## 2. SQLite indexes and migration

- Add explicit global, workflow-only, status-only, and workflow+status
  composite indexes for the four supported query shapes.
- Add the next Alembic revision without transforming existing rows.
- Test upgrade from the current main head, index presence, row preservation,
  downgrade, re-upgrade, and SQLite plans without a temporary order-by B-tree.

## 3. Service cursor and public page semantics (TDD)

- Implement canonical bounded cursor encoding/decoding and filter binding.
- Add typed cursor-invalid, cursor-not-found, and data-invalid failures.
- Extend `RunService` with a bounded `list_run_history` method that requests
  `limit + 1` and constructs the next cursor from the last public row.
- Test malformed/cross-filter/deleted cursors, exact page boundaries, insertion
  stability, invalid limits/filters, and persisted-data corruption.

## 4. FastAPI contract and generated client (TDD)

- Add strict `RunSummaryResponse` and `RunHistoryResponse` Pydantic models.
- Add synchronous `GET /api/v1/runs` before the run-detail route with explicit
  `operation_id=listRuns` and dedicated safe issue envelopes.
- Add API tests for success, empty, filters/pages, request bounds, error mapping,
  projection rejection, operation ID, threadpool execution, and leakage.
- Add dedicated `listRuns` request-validation and uncaught-exception handler
  tests so both 400 and 500 remain disclosure-safe `RunHistoryResponse` shapes.
- Export OpenAPI and regenerate Orval mechanically; verify generated-only drift.

## 5. Global navigation and `/runs` workbench (TDD)

- Add failing AppShell tests for the single Runs/Workflows/contextual authoring
  hierarchy and active rules.
- Register lazy `/runs` before `/runs/:runId` and build a focused RunHistoryPage.
- Use generated `listRuns`, TanStack `useInfiniteQuery`, strict client response
  guards, URL workflow/status filters, duplicate-cursor defense, and restart
  recovery.
- Add pure URL-filter parsing tests for duplicate, empty, unknown-status,
  oversized, and control-character values; invalid filters never issue a query
  or render their raw value. Prevent placeholder data from crossing exact
  filter query keys.
- Distinguish cursor protocol/replacement errors (reload first page, retaining
  loaded rows) from transport failures (retry, retaining same-key cached rows).
- Build compact desktop table/mobile records, truthful loading/error/empty/
  cached-data states, refresh/load-more actions, and succeeded QC/Artifact links.
- Add Vitest route/component coverage for all required states and browser
  back/forward behavior.

## 6. Real browser gate and responsive evidence

- Extend the existing desktop real-success Playwright path after canonical
  success, and add an independent successful run in the actual mobile project.
- Navigate through global Runs, find the real run ID, enter Activity, then visit
  QC and Artifacts deep links.
- Capture and inspect 1440×900, 1024×768, 390×844, and 360×800 screenshots;
  assert no page-level horizontal overflow or inaccessible controls.
- Preserve the existing cancellation/process-cleanup gate.

## 7. Final verification and delivery

- Run focused tests followed by full Python, durable Redis/RQ/Snakemake,
  frontend Vitest, typecheck, build, Playwright, OpenAPI/Orval drift, Ruff,
  formatting, and `git diff --check`.
- Confirm no task-owned Redis, worker, horse, Snakemake, Vite, or helper remains.
- Request two independent read-only merge-gate reviews (persistence/security and
  frontend UX/accessibility), fix every blocking/important finding, and re-run
  affected gates.
- Review `origin/main...HEAD`, commit without attribution trailers, push
  `codex/pr141-run-history-navigation`, create Draft PR141, and wait for exact
  HEAD CI to be fully green. Do not merge.
