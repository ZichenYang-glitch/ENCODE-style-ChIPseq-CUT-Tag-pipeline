# PR131 Read-only QC Metrics API Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Expose a bounded, lossless, run-isolated, read-only API for persisted workflow-neutral QC metrics.

**Architecture:** Extend the existing run repository and `RunService` with metric-ID keyset pagination, validate persisted values at the repository boundary, and project them through a dedicated strict Pydantic/FastAPI contract. OpenAPI remains the only source of generated TypeScript operations and DTOs; the request path has no workspace, adapter, worker, or lifecycle dependency.

**Tech Stack:** Python 3, FastAPI, Pydantic v2, SQLAlchemy 2, SQLite, pytest, OpenAPI, Orval, TypeScript.

---

## File structure

- `src/encode_pipeline/services/run_repositories.py`: paginated repository protocol, InMemory parity, and shared durable read validation.
- `src/encode_pipeline/persistence/repositories.py`: bounded run-scoped SQL keyset query and cursor-row validation.
- `src/encode_pipeline/services/runs.py`: positive-limit validation and repository delegation.
- `src/encode_pipeline/api/models.py`: strict lossless public metric and list-envelope models.
- `src/encode_pipeline/api/routes/qc_metrics.py`: thin read-only list route and stable QC issues.
- `src/encode_pipeline/api/routes/__init__.py`: router registration.
- `src/encode_pipeline/api/main.py`: QC-specific request-validation envelope.
- `test/persistence/test_qc_metrics_repository.py`: InMemory/SQLAlchemy query parity and corruption gates.
- `test/services/test_run_service.py`: service pagination boundary.
- `test/api/test_routes_qc_metrics.py`: API contract, redaction, reopen, and precision tests.
- `test/api/test_openapi_export.py`: stable operation and response-schema assertions.
- `frontend/openapi.json`, `frontend/src/api/generated/`: mechanical contract output only.

### Task 1: Add bounded repository and service queries

**Files:**
- Modify: `src/encode_pipeline/services/run_repositories.py`
- Modify: `src/encode_pipeline/persistence/repositories.py`
- Modify: `src/encode_pipeline/services/runs.py`
- Modify: `test/persistence/test_qc_metrics_repository.py`
- Modify: `test/services/test_run_service.py`

- [ ] Write failing parity tests that seed three deterministic metrics plus a
  second run, then assert sorted full reads, `limit=2`, `after=<metric_id>`, an
  unknown/cross-run cursor `KeyError`, and no cross-run rows.
- [ ] Add failing corruption tests proving both a selected row and the cursor
  row itself are rejected by shared durable validation rather than skipped.
- [ ] Run `PYTHONPATH=src python3 -m pytest -q test/persistence/test_qc_metrics_repository.py test/services/test_run_service.py` and verify the new keyword arguments fail.
- [ ] Change the protocol and service signature to
  `list_qc_metrics(run_id, *, after=None, limit=None)`; retain no-limit behavior
  for internal indexing tests and reject non-positive service limits.
- [ ] Implement InMemory sorted slicing after validating the cursor object and
  each selected object with `_validate_qc_metric_fields`.
- [ ] Implement SQLAlchemy cursor lookup with both keys, convert and validate
  that full row, then query `metric_id > after ORDER BY metric_id LIMIT limit`.
  Make `_qc_metric_from_row` call the shared validator before returning.
- [ ] Re-run the focused tests and require all to pass.

### Task 2: Define the lossless disclosure-safe API contract

**Files:**
- Modify: `src/encode_pipeline/api/models.py`
- Create: `test/api/test_routes_qc_metrics.py`

- [ ] Write failing model/route tests for exact `9007199254740993.123456789012`
  string output, nullable fields, UTC time, deterministic ID mismatch, unsafe
  display text, bad scope/token/unit/flag/source fields, and requested-run
  mismatch. Assert no private value or `technical_message` is returned.
- [ ] Add `QcMetricResponse` with `value: str`, controlled literal fields,
  canonical decimal projection, semantic-ID recomputation, safe public text,
  scope consistency, and UTC timestamp normalization.
- [ ] Add `RunQcMetricsResponse` with only the stable envelope fields.
- [ ] Run `PYTHONPATH=src python3 -m pytest -q test/api/test_routes_qc_metrics.py`
  and verify model-focused cases pass while the missing route cases remain red.

### Task 3: Add the thin read-only FastAPI route

**Files:**
- Create: `src/encode_pipeline/api/routes/qc_metrics.py`
- Modify: `src/encode_pipeline/api/routes/__init__.py`
- Modify: `src/encode_pipeline/api/main.py`
- Modify: `test/api/test_routes_qc_metrics.py`

- [ ] Add failing route cases for an empty run, stable two-page traversal,
  default/1/100 limits, rejected 0/101 limits, malformed cursor, unknown run,
  unknown cursor, cross-run cursor, damaged selected row, damaged cursor row,
  strict field projection, and SQLite reopen.
- [ ] Implement a synchronous FastAPI route so Starlette runs bounded SQLite
  I/O in its worker threadpool. Call only `RunService.get_run` and
  `RunService.list_qc_metrics(limit=limit + 1)`.
- [ ] Return stable sanitized issues: `RUN_NOT_FOUND` (404),
  `RUN_QC_METRIC_CURSOR_NOT_FOUND` (400), and
  `RUN_QC_METRIC_DATA_INVALID` (500).
- [ ] Register the router and add a `listRunQcMetrics` branch to the request
  validation handler that returns `RunQcMetricsResponse` and never includes
  Pydantic error text.
- [ ] Run the route test module and require all cases to pass.

### Task 4: Regenerate and verify OpenAPI/Orval

**Files:**
- Modify: `test/api/test_openapi_export.py`
- Regenerate: `frontend/openapi.json`
- Regenerate: `frontend/src/api/generated/`

- [ ] Add a failing operation-gate assertion for
  `GET /api/v1/runs/{run_id}/qc-metrics -> listRunQcMetrics`, declared QC
  envelopes, no default 422, and `QcMetricResponse.value` as OpenAPI string.
- [ ] Run `PYTHONPATH=src python3 -m pytest -q test/api/test_openapi_export.py`
  and verify the committed contract drift is red.
- [ ] Run `npm --prefix frontend run openapi:regenerate`; do not edit generated
  files manually.
- [ ] Run the OpenAPI/export and generated-client drift tests, frontend tests,
  typecheck, and production build.

### Task 5: Full verification, review, and publication

**Files:**
- Review: complete `origin/main...HEAD` diff

- [ ] Run focused persistence/service/API tests and
  `PYTHONPATH=src python3 -m pytest -q`.
- [ ] Run Ruff check, Ruff format check, and `git diff --check`.
- [ ] Confirm no migration, indexing/parser, worker, Snakemake, manifest,
  inventory, React product source, or package-lock file changed.
- [ ] Request one independent read-only merge-gate review focused on run
  isolation, cursor-row validation, Decimal precision, disclosure, and generated
  contract drift; fix every blocking finding and re-run affected gates.
- [ ] Commit without attribution trailers, push
  `agent/pr131-readonly-qc-api`, create Draft PR131, and wait for every required
  GitHub check to succeed. Do not merge.
