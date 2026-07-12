# PR127 Adapter Artifact Extraction Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Safely and atomically index ENCODE manifest-visible regular-file outputs after durable local execution succeeds.

**Architecture:** The ENCODE adapter produces neutral candidates from the existing manifest builder and artifact inventory. A platform service validates lifecycle, build identity, paths, filesystem types, and metadata, then asks `RunService` for one-transaction replacement plus event. Worker composition invokes the service only after canonical success and isolates extraction failures from RQ execution success.

**Tech Stack:** Python 3 dataclasses and pathlib/os stat, existing Result/Issue domain model, SQLAlchemy 2, SQLite, RQ/Redis, pytest.

---

### Task 1: Adapter-neutral extraction contract

**Files:**
- Modify: `src/encode_pipeline/platform/adapters.py`
- Modify: `src/encode_pipeline/platform/__init__.py`
- Test: `test/platform/test_adapters.py`

- [ ] Add failing tests for immutable candidate metadata, protocol method presence, and capability discovery.
- [ ] Run `PYTHONPATH=src python3 -m pytest test/platform/test_adapters.py -q` and confirm the new imports fail.
- [ ] Add `ExtractedArtifactCandidate` and `WorkflowAdapter.extract_artifacts(workspace)` without platform or persistence imports.
- [ ] Update test adapters with an explicit unsupported method and rerun the focused tests.

### Task 2: ENCODE manifest/catalog mapping and build identity

**Files:**
- Modify: `src/encode_pipeline/adapters/encode.py`
- Modify: `src/encode_pipeline/services/workflow_builds.py`
- Modify: `test/browser/platform_runtime.py`
- Test: `test/adapters/test_encode_adapter.py`
- Test: `test/services/test_workflow_builds.py`

- [ ] Add failing tests for present/missing/not-applicable rows, mixed assays, dynamic `biorep<N>` full matches, directory aggregate skips, unknown types, and result-manifest inclusion.
- [ ] Add a failing identity test that changes only `artifact-inventory.yaml` and expects a different digest.
- [ ] Implement direct `build_manifest_rows` and catalog mapping with controlled metadata and no TSV parsing, shell, glob, or content hashing.
- [ ] Add the inventory file to the build source manifest and copy it into the controlled browser project.
- [ ] Run the adapter, identity, Stage 49, Stage 50, and target-selection tests.

### Task 3: Atomic artifact replacement

**Files:**
- Modify: `src/encode_pipeline/services/run_repositories.py`
- Modify: `src/encode_pipeline/persistence/repositories.py`
- Modify: `src/encode_pipeline/services/runs.py`
- Test: `test/services/test_run_repositories.py`
- Test: `test/services/test_run_service.py`
- Test: `test/persistence/test_sqlalchemy_run_repository.py`

- [ ] Add failing shared contract tests for full replacement, SUCCEEDED-only enforcement, stable ordering, identical no-op, changed-set event, and rollback on invalid insert.
- [ ] Add `replace_artifacts(..., expected_status, event) -> RunEvent | None` to the repository contract.
- [ ] Implement InMemory swap plus event under one lock.
- [ ] Implement SQLAlchemy `BEGIN IMMEDIATE`, canonical status check, full delete/insert, equality/no-op check, and success event in one transaction.
- [ ] Add `RunService.replace_artifacts` validation and the safe `artifacts_indexed` draft.
- [ ] Reopen SQLite in a second service and prove the complete set and event survive restart.

### Task 4: Artifact extraction safety service

**Files:**
- Create: `src/encode_pipeline/services/artifact_extraction.py`
- Modify: `src/encode_pipeline/services/defaults.py`
- Test: `test/services/test_artifact_extraction.py`

- [ ] Add failing tests for lifecycle/capability/identity checks and explicit retry.
- [ ] Add table-driven failing tests for absolute, traversal, backslash, empty, NUL, duplicate, symlink-component, directory, FIFO, and metadata path attacks.
- [ ] Implement deterministic framed IDs, opaque URIs, ended-at timestamps, safe `lstat`, containment, regular-file checks, and bounded scalar metadata.
- [ ] Validate the complete candidate set before calling replace.
- [ ] Convert failures to a generic Result and a best-effort `artifact_extraction_failed` event containing only a stable reason code.
- [ ] Prove failed validation preserves a previously indexed complete set.

### Task 5: Worker success integration and failure isolation

**Files:**
- Modify: `src/encode_pipeline/workers/runtime.py`
- Modify: `src/encode_pipeline/workers/jobs.py`
- Test: `test/workers/test_runtime.py`
- Test: `test/workers/test_jobs.py`
- Test: `test/workers/test_tiny_execution_e2e.py`

- [ ] Add failing composition and job tests showing extraction happens only after canonical SUCCEEDED.
- [ ] Add a SQLite-backed worker test that reopens `RunService` and reads indexed refs/events.
- [ ] Inject the service into `WorkerRuntime` and call it from the LocalExecutionService success branch.
- [ ] Ensure Result failure and unexpected exception both leave the RQ job successful and the run SUCCEEDED with a safe failure event.
- [ ] Prove FAILED/CANCELLED/timeout paths never invoke extraction.

### Task 6: Verification, review, and publication

**Files:**
- Modify only implementation/test documentation if a real failure requires it.

- [ ] Run all focused tests and `PYTHONPATH=src python3 -m pytest`.
- [ ] Run real Redis/RQ/SIGALRM/cancellation/tiny Snakemake tests with zero integration skips.
- [ ] Run OpenAPI export drift without regenerating frontend files, changed-file Ruff/format, and `git diff --check`.
- [ ] Review the complete `origin/main...HEAD` diff for path leakage, transaction gaps, event duplication, worker/RQ truthfulness, and Stage 51 scope.
- [ ] Commit without attribution trailers, push `agent/pr127-adapter-artifact-extraction`, and create Draft PR127.
- [ ] Wait for every required GitHub Action to pass; fix real failures.
- [ ] Run final independent read-only domain/persistence and manifest/scientific merge-gate reviews; fix every blocking finding and rerun CI.
