# PR137 Validated Snapshot Submission Implementation Plan

**Goal:** Persist a successful backend validation as an opaque immutable
snapshot and require that snapshot for HTTP run creation, then connect the
PR136 workbench to validate, create, preflight, and explicit start.

**Architecture:** A workflow-neutral snapshot domain value is persisted through
InMemory/SQLAlchemy repository parity. `ValidatedInputService` owns validation
and identity capture; `ValidatedRunCreationService` owns first-use build checks;
`RunService` and the repository atomically consume the snapshot and create the
run. FastAPI routes are thin and the browser uses generated Orval operations
through TanStack Query.

---

### Task 0: Reproducible local environment doctor

- [x] Create `.local/envs/ci-fast` from the existing lock and install only
  project API/development dependencies.
- [x] Provide persistent Redis 7 in that environment without a system service.
- [x] Run frontend install, typecheck, build, and tests.
- [x] Start the real results demo and verify SQLite, worker, Snakemake, QC,
  artifact, and byte-identical download, then prove process cleanup.
- [x] Add a side-effect-free launcher doctor with controlled failures and tests.

### Task 1: Snapshot domain and repository TDD

- [x] Add canonical payload/digest, validation evidence, snapshot, and atomic
  creation-result domain values with strict invariants.
- [x] Add failing InMemory parity tests for create/read, tamper, expiry,
  cross-workflow, replay, replay conflict, and concurrent consumption.
- [x] Extend the repository contract and implement atomic InMemory behavior.
- [x] Add the SQLAlchemy model and Alembic migration from the current head.
- [x] Implement SQLAlchemy create/read/consume with `BEGIN IMMEDIATE`, strict row
  projection, rollback, restart, and concurrency tests.

### Task 2: Validation and creation services

- [x] Add failing tests for successful snapshot persistence, adapter validation
  failure, build capture failure/change, schema binding, and safe integer
  rejection.
- [x] Implement `ValidatedInputService` with before/after identity capture.
- [x] Add `RunService` snapshot read/atomic-create methods and a
  `ValidatedRunCreationService` for current identity/expiry/replay semantics.
- [x] Test no-run-on-failure and one-run-under-concurrency guarantees.

### Task 3: FastAPI, OpenAPI, and generated client

- [x] Add strict snapshot response and snapshot-only create request models.
- [x] Update validate/create routes, dependencies, composition, status/error
  envelopes, request limits, persistent restart tests, and explicit operation
  IDs.
- [x] Prove raw mutable inputs and client validation claims are rejected.
- [x] Export OpenAPI and regenerate Orval mechanically; update only adapters
  that bridge existing UI abstractions.

### Task 4: Workbench submission flow

- [x] Extend the draft reducer with a semantic revision/serialization identity
  and immediate snapshot invalidation on config/sample/options edits.
- [x] Add Review validation, backend Issue display, expiry state, and
  snapshot-only Create action using generated operations and TanStack mutations.
- [x] Navigate to the durable run URL and request existing preflight after a
  first or replay create response; leave Start explicit on the run page.
- [x] Update the legacy workflow detail path so it cannot bypass snapshots.
- [x] Cover successful validation, stale edit, backend rejection, lost response,
  replay, 409/unconfirmed state, URL/history, and accessibility.

### Task 5: Real browser path and full gates

- [x] Add real FastAPI browser coverage for author → validate → snapshot →
  create → preflight and refresh-safe run URL, without claiming PR138 results.
- [x] Run Python focused/full and durable Redis/RQ/Snakemake gates.
- [x] Run frontend Vitest, typecheck/build, desktop/mobile Playwright, OpenAPI
  export/Orval drift, Ruff/format/diff checks, and process cleanup checks.
- [ ] Request one independent adversarial read-only review, fix all blocking or
  important findings, and rerun affected/full gates.
- [ ] Commit, push, open Draft PR137, verify exact reviewed/local/remote/CI HEAD,
  wait for all required checks, mark ready, squash merge with
  `--match-head-commit`, fetch, and verify the merge on `origin/main`.
