# PR135 Adapter-Owned Authoring Contract Implementation Plan

**Goal:** Publish a versioned renderable adapter input contract and make inline
sample rows validate, persist, and materialize through the existing scientific
and platform boundaries.

**Architecture:** Add workflow-neutral, top-level-frozen schema metadata with
defensive copies at construction and serialization, project the uniform
platform input ceilings, publish an ENCODE Draft 2020-12 partial/complete
contract, bridge inline rows to the existing scientific validator with a
private deterministic TSV, and enforce actual authoring request bytes before
Pydantic. Keep lifecycle, persistence schema, worker scheduling, Snakemake
targets, and scientific rules unchanged.

**Tech stack:** Python 3, dataclasses, FastAPI/Pydantic, SQLAlchemy/SQLite,
JSON Schema Draft 2020-12, OpenAPI, Orval, pytest, React/Vite regression tests.

---

### Task 1: Workflow-neutral schema and input-limit contract

**Files:**

- Modify: `src/encode_pipeline/platform/adapters.py`
- Modify: `test/platform/test_adapters.py`

- [x] Write failing tests for version/dialect, per-surface coverage, modes,
  exact platform ceilings, top-level freezing, boundary copy isolation,
  JSON-safe schemas, and shared row bounds.
- [x] Add top-level-frozen `WorkflowSchemaCoverage`, `WorkflowAuthoringModes`,
  `WorkflowInputModes`, and `WorkflowInputLimits` values.
- [x] Require contract `1.0.0` limit values to equal the transport/domain
  ceilings; defer adapter-specific narrower values until both boundaries can
  enforce them in a future contract version.
- [x] Upgrade `WorkflowSchema.to_dict()` to the complete versioned envelope.
- [x] Export the authoring primitives, dialect, and hard ceilings from the
  public `encode_pipeline.platform` façade.
- [x] Apply the same structural ceilings in `WorkflowInputs` without adding
  ENCODE-specific field knowledge.
- [x] Run the focused platform contract tests.

### Task 2: ENCODE schema and inline scientific validation

**Files:**

- Create: `src/encode_pipeline/adapters/encode_authoring.py`
- Modify: `src/encode_pipeline/adapters/encode.py`
- Modify: `test/adapters/test_encode_adapter.py`
- Modify: workspace planning/materialization tests as needed

- [x] Write failing schema validity and tiny-profile coverage tests.
- [x] Write failing inline precedence, empty/unknown/control-character,
  cleanup, concurrency, path-policy, and no-temp-leak tests.
- [x] Publish stable config/sample/options JSON Schemas and add the
  `input_authoring` capability.
- [x] Implement fixed-header temporary TSV rendering and one private
  validation path that calls the existing config/sample validators.
- [x] Remove temporary paths from retained results and use the existing
  external-input policy during validation.
- [x] Prove two plans and equivalent inline/path inputs materialize identical
  canonical workspace bytes.

### Task 3: Typed API projection and pre-parser request boundary

**Files:**

- Modify: `src/encode_pipeline/api/models.py`
- Create: `src/encode_pipeline/api/request_limits.py`
- Modify: `src/encode_pipeline/api/main.py`
- Modify: `src/encode_pipeline/api/routes/workflows.py`
- Modify: `src/encode_pipeline/api/routes/runs.py`
- Modify: API route tests

- [x] Write failing tests for typed `schema`, null 404 schema, inline validate
  and create, bounded rows/columns/cells, and safe null validation value.
- [x] Write pure-ASGI tests for exact/+1, multiple chunks, absent and lying
  `Content-Length`, correct envelopes, zero downstream calls, and unaffected
  non-authoring routes.
- [x] Add strict Pydantic input aliases and strongly typed schema response
  models using JSON-value projections.
- [x] Keep `schema` required-but-nullable and project JSON values as a
  recursive union that Orval can generate without hand-written DTOs.
- [x] Add the route-scoped actual-byte middleware and explicit 413 contracts.
- [x] Convert validate/create endpoints to synchronous functions.
- [x] Add create -> SQLite reopen -> plan -> materialize coverage with inline
  snapshots and prove no temporary path is persisted.

### Task 4: OpenAPI and frontend contract synchronization

**Files:**

- Modify mechanically: `frontend/openapi.json`
- Regenerate mechanically: `frontend/src/api/generated/`
- Modify minimally: existing frontend compatibility adapter/types/tests
- Modify: `test/api/test_openapi_export.py`

- [x] Add OpenAPI assertions for stable operation IDs, typed schema envelope,
  request bounds, and 413 validation/create responses.
- [x] Export OpenAPI and run Orval regeneration; never hand-edit generated
  files.
- [x] Map the new response into the existing workflow detail compatibility
  surface without implementing the PR136 workbench.
- [x] Run frontend tests, typecheck, production build, and drift checks.

### Task 5: Full verification, review, and Draft PR

- [x] Run focused and full Python tests with `PYTHONPATH=src`.
- [x] Run the real Redis/RQ/SIGALRM/cancellation/tiny Snakemake gate.
- [x] Run Playwright regression, OpenAPI/Orval zero-drift, Ruff, format, and
  `git diff --check`.
- [x] Inspect the complete `origin/main...HEAD` diff for schema truthfulness,
  local-path disclosure, request bypasses, deterministic workspace bytes, and
  scope drift.
- [ ] Request an independent read-only merge-gate review and fix all blocking
  or important findings.
- [ ] Commit without attribution, push, create Draft PR135, wait for all
  required GitHub checks to pass, and keep the PR unmerged.
