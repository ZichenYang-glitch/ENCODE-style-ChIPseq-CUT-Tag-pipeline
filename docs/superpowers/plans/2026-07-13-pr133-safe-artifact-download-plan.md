# PR133 Safe Artifact Download Implementation Plan

> **Implementation discipline:** Use test-driven development for every security
> invariant. Keep paths out of the route and use the generated operation in the
> frontend; do not implement Range, checksums, or a second client.

**Goal:** Add a safe run-scoped artifact download endpoint and a truthful
artifact-inspector download action.

**Architecture:** `ArtifactDownloadService` resolves SQLite metadata and returns
an already-open descriptor-owned stream plan. A thin synchronous FastAPI route
projects stable JSON failures or a `StreamingResponse`. OpenAPI/Orval generate
the operation, while the existing transport boundary adds a shared blob
decoder. TanStack Query owns frontend mutation state.

**Tech stack:** Python 3, FastAPI, Pydantic, SQLAlchemy-backed `RunService`,
POSIX descriptors, Starlette `StreamingResponse`/`BackgroundTask`, OpenAPI,
Orval, React, TanStack Query v5, lucide-react, Vitest, Playwright.

---

### Task 1: Descriptor-owned download service

**Files:**
- Create: `src/encode_pipeline/services/artifact_downloads.py`
- Create: `test/services/test_artifact_downloads.py`

- [ ] Write failing tests for known-run lookup, missing/cross-run artifacts,
  persisted identity/path/name/MIME/size corruption, and a successful bounded
  stream.
- [ ] Add failing absolute/traversal/backslash/NUL/noncanonical path tests and
  symlink tests for workspace root, run, parent, and final components.
- [ ] Add failing directory, FIFO, device, missing, larger/smaller, inode swap,
  in-place mutation, early EOF, growth, read-error, generator-close, and
  idempotent-close tests.
- [ ] Implement `ArtifactDownloadService.prepare` and
  `ArtifactDownloadPlan.iter_bytes` with descriptor-chain fingerprints,
  64-KiB reads, expected-size bounds, `Result/Issue`, and sanitized reason codes.
- [ ] Run the focused service suite and confirm GREEN with no leaked descriptor.

### Task 2: Thin binary FastAPI contract

**Files:**
- Modify: `src/encode_pipeline/api/models.py`
- Modify: `src/encode_pipeline/api/dependencies.py`
- Modify: `src/encode_pipeline/api/main.py`
- Modify: `src/encode_pipeline/api/routes/artifacts.py`
- Create: `test/api/test_routes_artifact_download.py`
- Modify: `test/api/test_openapi_export.py`

- [ ] Write failing API tests for exact bytes, Content-Length, safe controlled
  filename, MIME/nosniff/cache headers, full response to Range, and unchanged
  run/events.
- [ ] Add missing/cross-run 404 parity, symlink/missing/size 409, corrupted row
  and unexpected 500, path/error disclosure, SQLite reopen, and descriptor-close
  tests.
- [ ] Add an AST/import boundary test proving the route does not use `os`,
  `Path`, `open`, or `FileResponse` for workspace access.
- [ ] Add stable `downloadRunArtifact` operation and binary/error OpenAPI schema
  assertions.
- [ ] Compose `ArtifactDownloadService`, add the dependency, synchronous route,
  `StreamingResponse`, and idempotent background cleanup. Run focused API tests.

### Task 3: Mechanical binary Orval operation

**Files:**
- Modify: `frontend/orval.config.js`
- Modify: `frontend/src/api/fetcher.ts`
- Modify: `frontend/src/api/fetcher.test.ts`
- Mechanically regenerate: `frontend/openapi.json`
- Mechanically regenerate: `frontend/src/api/generated/`

- [ ] Add failing fetcher tests showing JSON operations still call `json()`,
  binary operations call `blob()` even with `application/json` MIME, and JSON
  error envelopes remain redacted.
- [ ] Refactor the existing transport into one shared request/error path with
  exported JSON `fetcher` and blob `blobFetcher` decoders.
- [ ] Configure only `downloadRunArtifact` to use `blobFetcher`; export OpenAPI
  and regenerate Orval mechanically.
- [ ] Assert the generated operation URL, method, return type, mutator, and zero
  hand edits through contract/drift tests.

### Task 4: Artifact inspector download interaction

**Files:**
- Modify: `frontend/src/features/run-artifacts/ArtifactBrowser.tsx`
- Modify: `frontend/src/features/run-artifacts/ArtifactBrowser.test.tsx`
- Modify: `frontend/src/features/run-artifacts/ArtifactInspector.tsx`
- Modify: `frontend/src/features/run-artifacts/ArtifactInspector.test.tsx`
- Modify: `frontend/src/features/run-artifacts/artifactState.ts`
- Modify: `frontend/src/features/run-artifacts/artifactState.test.ts`

- [ ] Write failing tests for generated mutation use, pending disabled label,
  safe success download, URL revoke, redacted failure, Retry, selection identity,
  and controlled filename fallback.
- [ ] Add a pure filename helper and use a TanStack mutation in
  `ArtifactBrowser`; never write an endpoint path or DTO.
- [ ] Add the compact lucide Download action and polite status feedback to the
  existing inspector without changing list/detail truthfulness or layout.
- [ ] Run artifact, route, fetcher, and full frontend tests, typecheck, and build.

### Task 5: Real browser download gate

**Files:**
- Modify: `frontend/e2e/durable-run.spec.ts`

- [ ] Extend the real success run to request the selected persisted
  `result_manifest` through the generated UI action.
- [ ] Assert a real Playwright download, safe suggested filename, exact
  manifest bytes, and no page navigation or local path disclosure.
- [ ] Recheck desktop/mobile overflow and retain the real cancellation test.
- [ ] Run real Redis/RQ/Snakemake Playwright with zero skips and verify all
  Redis, worker, horse, Snakemake, Vite, and helper processes are reaped.

### Task 6: Full security and publication gate

**Files:**
- Review: complete `origin/main...HEAD` diff

- [ ] Run focused service/API/frontend tests, full Python, real
  Redis/RQ/SIGALRM/tiny/cancellation, full frontend, typecheck, production build,
  Playwright, OpenAPI/Orval drift, changed-file Ruff/format, and diff check.
- [ ] Audit descriptor ownership, midstream error behavior, identity races,
  header injection, cross-run isolation, error disclosure, and non-goal scope.
- [ ] Commit and push without attribution, create Draft PR133, and wait for all
  required CI jobs.
- [ ] Request one final independent read-only security/merge-gate review. Fix
  all blocking or important findings and request re-review of the new HEAD.
- [ ] If and only if CI/review/scope/head/base/worktree conditions all pass,
  mark Ready and squash merge with the exact expected head SHA; fetch and verify
  the merge before creating PR134.
