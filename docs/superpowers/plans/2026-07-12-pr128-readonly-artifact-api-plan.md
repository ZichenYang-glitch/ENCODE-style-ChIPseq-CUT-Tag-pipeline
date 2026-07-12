# PR128 Read-only Artifact API Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Expose bounded, run-isolated, disclosure-safe list and detail APIs for persisted `RunArtifactRef` records.

**Architecture:** Add run-scoped keyset queries to the repository and `RunService`, then project those domain records through strict Pydantic response models in a dedicated FastAPI artifacts router. OpenAPI remains the only source for Orval DTOs and operations; no filesystem or extraction dependency enters the request path.

**Tech Stack:** Python 3, FastAPI, Pydantic v2, SQLAlchemy 2, SQLite, pytest, OpenAPI, Orval, TypeScript.

---

### Task 1: Run-scoped repository queries

**Files:**
- Modify: `src/encode_pipeline/services/run_repositories.py`
- Modify: `src/encode_pipeline/persistence/repositories.py`
- Modify: `src/encode_pipeline/services/runs.py`
- Test: `test/services/test_run_repositories.py`
- Test: `test/persistence/test_sqlalchemy_run_repository.py`
- Test: `test/services/test_run_service.py`

- [ ] Add failing InMemory and SQLAlchemy tests for deterministic row order,
  `after=<artifact_id>`, bounded `limit`, invalid same-run cursor, cross-run
  cursor rejection, and `(run_id, artifact_id)` detail isolation.
- [ ] Change the repository protocol to:

```python
def list_artifacts(
    self,
    run_id: str,
    *,
    after: str | None = None,
    limit: int | None = None,
) -> tuple[RunArtifactRef, ...]: ...

def get_artifact(self, run_id: str, artifact_id: str) -> RunArtifactRef: ...
```

- [ ] Implement InMemory cursor lookup against values sorted by `artifact_id`;
  raise `KeyError` when the cursor or detail does not exist in that run.
- [ ] Implement SQLAlchemy keyset pagination by first resolving the cursor row
  with both `run_id` and `artifact_id`, then selecting rows for that run with a
  greater artifact ID, ordered by artifact ID and limited in SQL.
- [ ] Add `RunService.list_artifacts(..., after, limit)` validation and
  `RunService.get_artifact(run_id, artifact_id)` without exposing ORM rows.
- [ ] Run:

```bash
PYTHONPATH=src python3 -m pytest -q \
  test/services/test_run_repositories.py \
  test/persistence/test_sqlalchemy_run_repository.py \
  test/services/test_run_service.py
```

### Task 2: Strict Pydantic artifact projection

**Files:**
- Modify: `src/encode_pipeline/api/models.py`
- Test: `test/api/test_routes_artifacts.py`

- [ ] Add failing tests that insert malicious legacy records containing an
  absolute relative path, `file://` URI, separator-bearing name, hostile known
  metadata string, and arbitrary extra metadata. Assert unsafe known fields
  produce a path-free error and unknown keys are never returned.
- [ ] Add controlled metadata, `ArtifactReferenceResponse`,
  `RunArtifactsResponse`, and `RunArtifactDetailResponse` models.
- [ ] Require canonical `results/...` POSIX paths, non-negative sizes, logical
  identifiers, safe MIME, safe leaf names, and an exact opaque URI:

```text
run://runs/{quoted_run_id}/artifacts/{artifact_id}
```

- [ ] Expose `relative_path`, `output_type`, and `size_bytes` as explicit fields;
  whitelist PR127 catalog/logical metadata fields and ignore all other
  persisted keys so arbitrary adapter metadata cannot become a new disclosure
  channel without an explicit API contract update.

### Task 3: Read-only FastAPI artifact routes

**Files:**
- Create: `src/encode_pipeline/api/routes/artifacts.py`
- Modify: `src/encode_pipeline/api/routes/__init__.py`
- Test: `test/api/test_routes_artifacts.py`

- [ ] Add failing route tests for list, empty list, detail, stable order, two-page
  traversal, limit 1 and 100, invalid 0/101 limits, unknown run, unknown
  artifact, cross-run artifact, invalid/cross-run cursor, and SQLite reopen.
- [ ] Implement list endpoint with `Query(default=50, ge=1, le=100)`, fetch
  `limit + 1`, and return the last visible artifact ID as `next_cursor`.
- [ ] Implement detail using only `RunService.get_artifact(run_id,
  artifact_id)`.
- [ ] Return stable envelopes:

```text
404 RUN_NOT_FOUND
404 RUN_ARTIFACT_NOT_FOUND
400 RUN_ARTIFACT_CURSOR_NOT_FOUND
500 RUN_ARTIFACT_DATA_INVALID
```

- [ ] Set explicit operation IDs `listRunArtifacts` and `getRunArtifact` and
  document 400/404/500 response models.
- [ ] Run `PYTHONPATH=src python3 -m pytest -q test/api/test_routes_artifacts.py`.

### Task 4: OpenAPI and Orval contract

**Files:**
- Modify: `test/api/test_openapi_export.py`
- Regenerate: `frontend/openapi.json`
- Regenerate: `frontend/src/api/generated/`

- [ ] Extend `EXPECTED_OPERATIONS` with both artifact endpoints and assert
  unique stable operation IDs.
- [ ] Run `npm --prefix frontend run openapi:regenerate`; do not edit generated
  files manually.
- [ ] Run:

```bash
PYTHONPATH=src python3 -m pytest -q \
  test/api/test_openapi_export.py \
  test/api/test_generated_client_drift.py
npm --prefix frontend run typecheck
npm --prefix frontend test -- --run
npm --prefix frontend run build
```

### Task 5: Full verification, review, and publication

**Files:**
- Review: complete `origin/main...HEAD` diff

- [ ] Run focused persistence/API tests, then
  `PYTHONPATH=src python3 -m pytest -q`.
- [ ] Run the existing real Redis/RQ/Snakemake durable gates to prove the
  read-only query change did not alter execution behavior.
- [ ] Run changed-file Ruff, `ruff format --check`, and `git diff --check`.
- [ ] Confirm no migration, frontend product code, package-lock, workflow,
  manifest, inventory, or Stage 51 file changed.
- [ ] Request one independent read-only merge-gate review for run isolation,
  disclosure, pagination, persistence, and OpenAPI; fix every blocker.
- [ ] Commit without attribution trailers, push
  `agent/pr128-readonly-artifact-api`, open Draft PR128, and wait until CI,
  browser E2E, frontend, lint, and lock-check are green. Do not merge.
