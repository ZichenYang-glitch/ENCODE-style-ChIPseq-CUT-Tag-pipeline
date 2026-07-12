# PR130 QC Summary Indexing Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Persist a small workflow-neutral set of trusted ENCODE QC summary metrics after successful durable execution without changing scientific targets or run terminal state.

**Architecture:** Existing artifact indexing supplies a complete persisted source set. The platform safely reads only adapter-declared summary artifacts through descriptor-relative `O_NOFOLLOW` traversal, the ENCODE adapter maps fixed TSV contracts into neutral Decimal candidates, and the repository atomically validates the expected artifact set while replacing QC rows and recording redacted outcomes.

**Tech Stack:** Python 3.11, dataclasses, `Decimal`, stdlib `csv`, SQLAlchemy 2, Alembic, SQLite, RQ, pytest, Ruff.

---

## File structure

- `src/encode_pipeline/platform/adapters.py`: optional QC adapter protocol and neutral source/candidate values.
- `src/encode_pipeline/platform/runs.py`: immutable durable `RunQcMetric` domain value.
- `src/encode_pipeline/adapters/encode_qc.py`: exact ENCODE TSV headers, metric catalog, parser, and semantic bounds.
- `src/encode_pipeline/adapters/encode.py`: advertise `qc_summary_extract` and delegate source types/parsing.
- `src/encode_pipeline/services/qc_summary_indexing.py`: lifecycle/build checks, descriptor-safe source reads, candidate validation, deterministic IDs, and redacted outcomes.
- `src/encode_pipeline/services/run_repositories.py`: repository QC contracts and InMemory parity.
- `src/encode_pipeline/services/runs.py`: locked QC persistence/query/failure methods.
- `src/encode_pipeline/persistence/models.py`: private SQLAlchemy QC row.
- `src/encode_pipeline/persistence/repositories.py`: transactional SQL implementation and artifact-driven QC invalidation.
- `src/encode_pipeline/persistence/alembic/versions/20260712_05_run_qc_metrics.py`: additive QC table migration.
- `src/encode_pipeline/services/defaults.py`, `src/encode_pipeline/workers/runtime.py`, `src/encode_pipeline/workers/jobs.py`: worker-only composition and post-success ordering.
- `test/adapters/test_encode_qc.py`, `test/services/test_qc_summary_indexing.py`: parser and platform boundary tests.
- Existing platform, persistence, migration, runtime, job, and workflow-build tests: parity and integration coverage.

### Task 1: Add neutral QC contracts and domain model

**Files:**
- Modify: `src/encode_pipeline/platform/adapters.py`
- Modify: `src/encode_pipeline/platform/runs.py`
- Modify: `test/platform/test_adapters.py`
- Create: `test/platform/test_qc_metrics.py`

- [ ] **Step 1: Write failing contract and domain tests**

Cover immutable source metadata, bytes-only document content, Decimal-only
candidates, optional runtime protocol recognition, safe scope combinations,
and `RunQcMetric.to_dict()` copy behavior. Use an adapter that implements:

```python
class QcAdapter:
    def qc_source_output_types(self):
        return ("qc_summary",)

    def extract_qc_metrics(self, inputs, sources):
        return Result.success(())
```

- [ ] **Step 2: Run tests and verify the new imports fail**

Run: `PYTHONPATH=src python3 -m pytest test/platform/test_adapters.py test/platform/test_qc_metrics.py -q`

Expected: collection failure for missing QC contract/domain names.

- [ ] **Step 3: Implement the minimal neutral values**

Add frozen `QcSourceArtifact`, `QcSourceDocument`, and
`ExtractedQcMetricCandidate` dataclasses plus runtime-checkable
`QcSummaryExtractingAdapter`. Add frozen `RunQcMetric` with explicit fields and
no metadata mapping. Keep semantic validation in the platform service.

- [ ] **Step 4: Run the focused tests**

Run: `PYTHONPATH=src python3 -m pytest test/platform/test_adapters.py test/platform/test_qc_metrics.py -q`

Expected: pass.

- [ ] **Step 5: Commit**

```bash
git add src/encode_pipeline/platform/adapters.py src/encode_pipeline/platform/runs.py test/platform/test_adapters.py test/platform/test_qc_metrics.py
git commit -m "feat: define QC summary contracts"
```

### Task 2: Implement the ENCODE fixed-contract parser

**Files:**
- Create: `src/encode_pipeline/adapters/encode_qc.py`
- Modify: `src/encode_pipeline/adapters/encode.py`
- Create: `test/adapters/test_encode_qc.py`
- Modify: `test/adapters/test_encode_adapter.py`

- [ ] **Step 1: Write failing parser tests**

Build in-memory `QcSourceDocument` fixtures for all three exact writer headers.
Assert the 13 locked metric keys, scopes, IDs, assays, units, `NA` skipping,
pooled `ok/mismatch/unknown` flag mapping, and mixed-assay sorting. Reject wrong
headers, missing/extra rows, identity mismatches, duplicate sources, malformed
UTF-8, illegal units, exponential or over-precision decimals, negative counts,
fractions outside 0..1, NaN, and Infinity.

- [ ] **Step 2: Verify RED**

Run: `PYTHONPATH=src python3 -m pytest test/adapters/test_encode_qc.py test/adapters/test_encode_adapter.py -q`

Expected: missing parser/delegates.

- [ ] **Step 3: Implement strict stdlib parsing**

Use `csv.DictReader(io.StringIO(content.decode("utf-8", errors="strict")), delimiter="\t")` with exact header equality and exactly one non-empty row. A
code-owned catalog maps each source field to metric key, label, unit, and
bounds. Parse via `Decimal`; accept literal `NA` as omitted, reject exponent
syntax and values with more than 12 fractional digits.

- [ ] **Step 4: Advertise and delegate the optional capability**

Add `qc_summary_extract` to ENCODE capabilities. Return only the three exact
source output types and delegate parsing to `encode_qc` without reading paths.

- [ ] **Step 5: Run parser/adapter tests**

Run: `PYTHONPATH=src python3 -m pytest test/adapters/test_encode_qc.py test/adapters/test_encode.py test/adapters/test_encode_adapter.py test/platform/test_adapters.py -q`

Expected: pass.

- [ ] **Step 6: Commit**

```bash
git add src/encode_pipeline/adapters/encode.py src/encode_pipeline/adapters/encode_qc.py test/adapters/test_encode_qc.py test/adapters/test_encode_adapter.py
git commit -m "feat: map ENCODE QC summaries"
```

### Task 3: Add lossless QC persistence and migration

**Files:**
- Modify: `src/encode_pipeline/persistence/models.py`
- Create: `src/encode_pipeline/persistence/alembic/versions/20260712_05_run_qc_metrics.py`
- Modify: `test/persistence/test_migrations.py`

- [ ] **Step 1: Write migration tests first**

Require `run_qc_metrics`, head `20260712_05`, run cascade, unique
`(run_id, metric_id)`, run index, scope/flag/scope-ID checks, and downgrade to
`20260712_04` without changing existing tables. Seed a current-head database,
upgrade, and verify no lifecycle/artifact rows change.

- [ ] **Step 2: Verify RED**

Run: `PYTHONPATH=src python3 -m pytest test/persistence/test_migrations.py -q`

Expected: missing table/head revision assertions fail.

- [ ] **Step 3: Add the private row and additive Alembic revision**

Store `value_text` as bounded text, not SQL `NUMERIC`. Include checks for
scope/flag and required IDs. Do not add an artifact foreign key.

- [ ] **Step 4: Run migration tests**

Run: `PYTHONPATH=src python3 -m pytest test/persistence/test_migrations.py -q`

Expected: pass including independent downgrade/re-upgrade.

- [ ] **Step 5: Commit**

```bash
git add src/encode_pipeline/persistence/models.py src/encode_pipeline/persistence/alembic/versions/20260712_05_run_qc_metrics.py test/persistence/test_migrations.py
git commit -m "feat: add durable QC metric schema"
```

### Task 4: Add atomic repository parity and artifact invalidation

**Files:**
- Modify: `src/encode_pipeline/services/run_repositories.py`
- Modify: `src/encode_pipeline/services/runs.py`
- Modify: `src/encode_pipeline/persistence/repositories.py`
- Modify: `test/services/test_run_service.py`
- Modify: `test/persistence/test_sqlalchemy_run_repository.py`

- [ ] **Step 1: Write failing InMemory and SQLAlchemy tests**

Cover complete replace, stable sorting, restart round-trip, canonical decimal
text above `2^53` and at the 26+12 digit boundary, exact retry no-op, changed
set replacement, insert/event rollback, duplicate IDs, expected artifact set
order independence, stale expected source rejection, concurrent artifact/QC
replacement, and repeated identical failure no-op.

Use the domain-owned framed SHA-256 helper for every fixture ID. Prove both
repositories atomically reject a well-formed digest from different semantic
coordinates and cannot persist identical coordinates under two legal hashes.
Use the same domain token validator in service and repository parity tests,
including successful `-S1`, `.S1`, `_S1` and rejected slash, exact `..`,
control, empty, and overlong identifiers.

Also prove a non-equivalent artifact replacement atomically clears QC and adds
one `qc_metrics_invalidated`, while equivalent reordered artifacts preserve QC
and first-ever artifact indexing writes no invalidation.

- [ ] **Step 2: Verify RED**

Run: `PYTHONPATH=src python3 -m pytest test/services/test_run_service.py test/persistence/test_sqlalchemy_run_repository.py -q`

Expected: missing repository/service operations.

- [ ] **Step 3: Implement repository operations**

Add `replace_qc_metrics`, `list_qc_metrics`, and atomic
`record_qc_metrics_failure`. SQLAlchemy must call `_begin_write`, re-read the
run and all current `RunArtifactRow` values, compare complete sets after sorting
by artifact ID, validate every source, then replace rows and add the event in
the same transaction. Serialize Decimal with a strict canonical helper and
parse it back without float conversion.

- [ ] **Step 4: Implement artifact-driven invalidation**

Before changing a non-equivalent artifact set, clear QC rows. After the
`artifacts_indexed` event, add `qc_metrics_invalidated` only if prior QC rows or
an outcome exist and latest outcome is not already invalidated. Mirror this
under the InMemory lock.

- [ ] **Step 5: Run persistence tests**

Run: `PYTHONPATH=src python3 -m pytest test/services/test_run_service.py test/persistence/test_sqlalchemy_run_repository.py test/persistence/test_runtime.py -q`

Expected: pass.

- [ ] **Step 6: Commit**

```bash
git add src/encode_pipeline/services/run_repositories.py src/encode_pipeline/services/runs.py src/encode_pipeline/persistence/repositories.py test/services/test_run_service.py test/persistence/test_sqlalchemy_run_repository.py
git commit -m "feat: persist complete QC metric indexes"
```

### Task 5: Implement the platform-safe QC indexing service

**Files:**
- Create: `src/encode_pipeline/services/qc_summary_indexing.py`
- Create: `test/services/test_qc_summary_indexing.py`

- [ ] **Step 1: Write failing lifecycle and safety tests**

Cover deterministic IDs, safe source references, build/capability checks,
successful empty sets, idempotent retry, prior-set preservation on failure,
FAILED/CANCELLED exclusion, redacted failure events, and exact artifact-set
validation. Attack absolute, traversal, backslash, NUL, duplicate path,
symlinked workspace component, final symlink, directory, FIFO, oversized file,
persisted `size_bytes` larger/smaller/type/bound violations, read-time mutation,
and a race that swaps a source for an outside symlink before final open. Assert
outside sentinel bytes and size-mismatched content never reach the adapter and
that a prior complete QC set survives failure.

- [ ] **Step 2: Verify RED**

Run: `PYTHONPATH=src python3 -m pytest test/services/test_qc_summary_indexing.py -q`

Expected: service import failure.

- [ ] **Step 3: Implement descriptor-only source reads**

Traverse from `os.open("/", O_DIRECTORY | O_CLOEXEC)` using `dir_fd`,
`O_NOFOLLOW`, and directory descriptors. Open the final source with
`O_RDONLY | O_NOFOLLOW | O_NONBLOCK | O_CLOEXEC`, require regular `fstat` and
exact agreement with validated persisted `size_bytes`, read at most expected
size + 1 through `os.read`, require exact final content length, compare pre/post
fd identity, and close all descriptors in `finally` blocks. Document that
same-size content replacement remains undetectable without a checksum.

- [ ] **Step 4: Implement candidate validation and persistence**

Validate bounded tokens, scope-ID consistency, finite Decimal grammar and
precision, controlled units/flags, source membership, duplicate source paths,
and duplicate semantic metric IDs before calling the atomic RunService method.
Use the domain-owned identifier validator and framed metric-ID builder in both
service and repository; the repository must recompute each ID before writing.
Map all exceptions to stable reason codes with no path or exception text.

- [ ] **Step 5: Run service tests**

Run: `PYTHONPATH=src python3 -m pytest test/services/test_qc_summary_indexing.py test/services/test_artifact_extraction.py -q`

Expected: pass.

- [ ] **Step 6: Commit**

```bash
git add src/encode_pipeline/services/qc_summary_indexing.py test/services/test_qc_summary_indexing.py
git commit -m "feat: safely index QC summaries"
```

### Task 6: Compose the worker post-success path and identity coverage

**Files:**
- Modify: `src/encode_pipeline/services/defaults.py`
- Modify: `src/encode_pipeline/workers/runtime.py`
- Modify: `src/encode_pipeline/workers/jobs.py`
- Modify: `test/services/test_workflow_builds.py`
- Modify: `test/workers/test_runtime.py`
- Modify: `test/workers/test_jobs.py`
- Modify: `test/workers/test_tiny_execution_e2e.py`

- [ ] **Step 1: Write failing composition and worker tests**

Require runtime reconstruction of `QcSummaryIndexingService`; success ordering
of execution → artifacts → QC; persisted metrics readable through a newly
opened SQLite service; empty source success count 0; parser/safety failure
leaving run `SUCCEEDED` with one redacted failure event; and no QC call on
FAILED/CANCELLED/timeout paths. Prove modifying `encode_qc.py` changes the build
digest.

- [ ] **Step 2: Verify RED**

Run: `PYTHONPATH=src python3 -m pytest test/services/test_workflow_builds.py test/workers/test_runtime.py test/workers/test_jobs.py test/workers/test_tiny_execution_e2e.py -q`

Expected: missing composition and worker assertions fail.

- [ ] **Step 3: Add default/runtime composition**

Construct the QC service from the same registry, RunService, identity provider,
and workspace root as artifact extraction. Store it on `WorkerRuntime` and
preserve close/failure cleanup.

- [ ] **Step 4: Add contained post-success invocation**

Pass only a successful complete artifact tuple to QC indexing. If artifact
extraction fails, write a stable QC source-unavailable failure. Contain expected
and defensive QC failures without raising execution failure or changing
`SUCCEEDED`.

- [ ] **Step 5: Run worker tests**

Run: `PYTHONPATH=src python3 -m pytest test/services/test_workflow_builds.py test/workers/test_runtime.py test/workers/test_jobs.py test/workers/test_tiny_execution_e2e.py -q`

Expected: pass.

- [ ] **Step 6: Commit**

```bash
git add src/encode_pipeline/services/defaults.py src/encode_pipeline/workers/runtime.py src/encode_pipeline/workers/jobs.py test/services/test_workflow_builds.py test/workers/test_runtime.py test/workers/test_jobs.py test/workers/test_tiny_execution_e2e.py
git commit -m "feat: index QC after durable execution"
```

### Task 7: Full verification and scope audit

**Files:**
- Modify only files required by discovered PR130 failures.

- [ ] **Step 1: Run focused PR130 tests**

```bash
PYTHONPATH=src python3 -m pytest \
  test/platform/test_qc_metrics.py \
  test/adapters/test_encode_qc.py \
  test/services/test_qc_summary_indexing.py \
  test/persistence/test_migrations.py \
  test/persistence/test_sqlalchemy_run_repository.py \
  test/workers/test_runtime.py test/workers/test_jobs.py -q
```

- [ ] **Step 2: Run the full Python suite against this worktree**

Run: `PYTHONPATH=src python3 -m pytest`

Expected: all tests pass; only established environment/profile skips remain.

- [ ] **Step 3: Run real durable worker gates**

Run the repository CI command covering real Redis/RQ/SIGALRM/cancellation/tiny
execution. Confirm zero skips in the selected real tests.

- [ ] **Step 4: Prove excluded surfaces have no drift**

Run OpenAPI export and Orval regeneration drift checks without accepting any
generated changes. Assert no diff in `frontend/`, OpenAPI, package-lock,
Snakemake targets, rules, manifest code, or artifact inventory.

- [ ] **Step 5: Run lint and diff checks**

```bash
PYTHONPATH=src python3 -m ruff check <changed-python-files>
PYTHONPATH=src python3 -m ruff format --check <changed-python-files>
git diff --check origin/main...HEAD
```

- [ ] **Step 6: Review the complete diff**

Audit lifecycle isolation, descriptor closure, concurrency, decimal round-trip,
idempotency, failure redaction, build identity, migration downgrade, and scope.
Then request the single final read-only merge-gate review required by the user.

- [ ] **Step 7: Publish without merging**

Push `agent/pr130-qc-summary-indexing`, create Draft PR130 targeting `main`, wait
for every required GitHub check, fix real failures minimally, and retain Draft.
