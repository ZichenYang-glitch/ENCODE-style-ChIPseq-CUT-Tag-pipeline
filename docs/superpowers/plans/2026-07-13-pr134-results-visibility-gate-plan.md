# PR134 Results Visibility Gate Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Prove a real run can produce persisted QC metrics, navigate to their source artifact, and download exact bytes in desktop/mobile browsers, while giving the local demo the same deterministic fixture.

**Architecture:** Add one sentinel-owned controlled-project builder under `scripts/` and reuse it from Playwright and an opt-in local demo flag. The fixture emits valid, absent, or malformed machine-readable QC sources without touching the scientific workflow; existing worker indexing, SQLite APIs, React views, and the PR133 download boundary remain unchanged.

**Tech Stack:** Python 3, Snakemake, Redis/RQ, SQLite/SQLAlchemy, FastAPI, React Router, TanStack Query, Orval, Playwright, Vitest, pytest.

---

### Task 1: Sentinel-owned deterministic results fixture

**Files:**
- Create: `scripts/results_visibility_fixture.py`
- Create: `test/browser/test_results_visibility_fixture.py`

- [x] **Step 1: Write failing ownership and deterministic-input tests**

Add tests that call `prepare_results_visibility_fixture(project_root)` and
require a `ResultsVisibilityInputs` value with absolute `samples_path`, four
configs whose `threads` values are 1–4, and `outdir == "results"`. Require the
generated project to contain package source, artifact inventory, default
profile, Snakefile, task helper, samples, and an exact ownership sentinel.

```python
inputs = prepare_results_visibility_fixture(tmp_path / "project")
assert inputs.results_config["threads"] == 1
assert inputs.cancel_config["threads"] == 2
assert inputs.empty_config["threads"] == 3
assert inputs.malformed_config["threads"] == 4
assert {config["outdir"] for config in inputs.configs} == {"results"}
assert inputs.samples_path.is_absolute()
```

Add a second test that places an unrelated file in an existing project root
without the sentinel and requires `ValueError` without deleting that file.
Add a third test that prepares twice and requires only the sentinel-owned
project to be replaced while a sibling SQLite marker remains intact.

- [x] **Step 2: Run the new tests and confirm RED**

Run:

```bash
PYTHONPATH=src python3 -m pytest -q \
  test/browser/test_results_visibility_fixture.py
```

Expected: collection fails because `scripts.results_visibility_fixture` does
not exist.

- [x] **Step 3: Implement the fixture builder**

Create a frozen input dataclass and one public preparation function:

```python
@dataclass(frozen=True)
class ResultsVisibilityInputs:
    samples_path: Path
    results_config: dict[str, object]
    cancel_config: dict[str, object]
    empty_config: dict[str, object]
    malformed_config: dict[str, object]
    expected_qc_summary: str

def prepare_results_visibility_fixture(
    project_root: Path,
    *,
    repository_root: Path = REPOSITORY_ROOT,
) -> ResultsVisibilityInputs:
    ...
```

Resolve and reject `/`, `/tmp`, home, the repository root, and paths outside a
dedicated child. If the directory exists, require the exact
`.encode-results-visibility-fixture` sentinel before `shutil.rmtree`. Write the
sentinel immediately after creating the new root; on later failure, remove
only that owned root.

Copy current `src/encode_pipeline`, `pyproject.toml`, artifact inventory, the
minimal default profile, and the existing platform worker sample sheet. Build
the QC content from the fixed `_QC_HEADER`: all columns `NA` except `sample=C1`,
`assay=chipseq`, `total_reads=1000`, `frip=0.125`, `peak_count=50`,
`percent_duplication=0.2`, `estimated_library_size=9007199254740993`,
`nrf=0.8`, `pbc1=0.75`, and `pbc2=3.0`.

Write a controlled Snakefile with explicit mode selection and separate result
and QC outputs. The cancellation mode must not schedule the QC output. The
helper writes a truthful manifest for results/empty/malformed modes, a valid or
malformed QC TSV as requested, and the existing cancellation PID evidence.

- [x] **Step 4: Run focused fixture and build-identity tests; confirm GREEN**

Run:

```bash
PYTHONPATH=src python3 -m pytest -q \
  test/browser/test_results_visibility_fixture.py \
  test/services/test_workflow_builds.py
```

Expected: all pass and no path outside the test temporary root is changed.

- [x] **Step 5: Commit the fixture boundary**

```bash
git add scripts/results_visibility_fixture.py \
  test/browser/test_results_visibility_fixture.py
git commit -m "test: add deterministic results visibility fixture"
```

### Task 2: Reuse the fixture in the real platform harness

**Files:**
- Modify: `test/browser/platform_runtime.py`
- Modify: `test/browser/test_platform_runtime.py`

- [ ] **Step 1: Write failing runtime-manifest assertions**

Extend the browser runtime test to require `runtime.json` to contain
`resultsConfig`, `cancelConfig`, `emptyConfig`, `malformedConfig`, and the
expected QC summary text, with no database URL, Redis URL, environment value,
or workspace absolute path beyond the already controlled runtime fields.

- [ ] **Step 2: Run the focused tests and confirm RED**

Run:

```bash
PYTHONPATH=src python3 -m pytest -q \
  test/browser/test_platform_runtime.py
```

Expected: runtime keys and non-empty QC assertions fail against the old inline
project.

- [ ] **Step 3: Replace the inline browser project with the shared builder**

Delete `create_controlled_project` and its embedded scripts from
`test/browser/platform_runtime.py`. Call
`prepare_results_visibility_fixture(project_root)` and pass the returned input
object into `write_manifest`. Serialize only the four input configs, sample
path, expected QC summary, runtime roots, workflow ID, and queue name.

Do not change the existing bundled-workflow tiny execution or cancellation
tests. The real Playwright stack supplies the independent API/worker/Snakemake
process boundary for this controlled fixture.

- [ ] **Step 4: Run focused runtime tests; confirm GREEN**

```bash
PYTHONPATH=src python3 -m pytest -q test/browser/test_platform_runtime.py
```

Expected: all pass; the serialized runtime input contains only the controlled
fixture contract. Task 3 proves durable artifacts and metrics through the real
independent-process stack.

- [ ] **Step 5: Commit real harness reuse**

```bash
git add test/browser/platform_runtime.py \
  test/browser/test_platform_runtime.py
git commit -m "test: index real QC results in durable gates"
```

### Task 3: Browser run-to-QC-to-download product gate

**Files:**
- Modify: `frontend/e2e/durable-run.spec.ts`

- [ ] **Step 1: Write the non-empty QC and source-download assertions**

Extend `RuntimeManifest` with the shared fixture values. After the results run
is SUCCEEDED, navigate to QC before Artifacts and assert:

```typescript
await expect(page.getByRole('heading', { name: 'Indexed QC metrics' })).toBeVisible();
await expect(page.getByText('Total reads', { exact: true })).toBeVisible();
await expect(page.getByText('1000', { exact: true })).toBeVisible();
await page.getByRole('button', {
  name: 'Open source artifact for Total reads',
}).click();
await expect(page.getByText('results/C1/01_qc/C1.qc_summary.tsv', {
  exact: true,
})).toBeVisible();
```

Download `C1.qc_summary.tsv`, require the suggested filename and exact bytes,
reload the artifact deep link, then return to QC and reload that deep link.

Add helpers that run the empty and malformed configs through the same real
create/preflight/start path. Require `No indexed QC metrics` only for empty;
require `QC metric indexing failed` and canonical `succeeded` for malformed.

- [ ] **Step 2: Update screenshot assertions for real QC content**

Change `captureQcViewport` to require `qc-metric-list`, the total-reads row,
the source action, and no horizontal overflow. Keep 1440x900, 1024x768,
390x844, and 360x800 captures. Change artifact captures to the selected QC
summary source and preserve desktop inspector-right/mobile inspector-below
bounding-box assertions.

- [ ] **Step 3: Run the real Playwright acceptance gate**

Run:

```bash
ENCODE_PIPELINE_E2E_REDIS_URL=redis://127.0.0.1:6392/14 \
PYTHONPATH=../src npm --prefix frontend run test:e2e
```

The pre-PR134 baseline is the already-recorded confirmed-empty QC path. With
Tasks 1–2 and the new assertions present, both Playwright projects must pass
with zero skips and the results run must expose eight metrics.

- [ ] **Step 4: Verify browser downloads and cleanup**

Require the direct Range canary to stay full `200` with `Accept-Ranges: none`.
After Playwright exits, assert its owner-marked runtime root is gone and no
Redis/worker/horse/Snakemake/Vite/helper process from the unique queue or
runtime remains.

- [ ] **Step 5: Commit the browser gate**

```bash
git add frontend/e2e/durable-run.spec.ts
git commit -m "test: prove run QC artifact download visibility"
```

### Task 4: One-command results demo and runtime documentation

**Files:**
- Modify: `scripts/run_local_platform.py`
- Modify: `test/browser/test_platform_runtime.py`
- Modify: `docs/development/local-platform-runtime.md`

- [ ] **Step 1: Write failing demo-mode tests**

Parse `--results-visibility-demo` and require it to use
`.local/results-visibility-demo` when `--runtime-root` is not provided, create
`results-visibility-project` through the shared builder, and write
`results-visibility-inputs.json` containing only workflow ID, success config,
and samples path. Require ordinary launcher defaults to remain
`.local/platform-demo` and the repository project root.

Test that explicit `--project-root` plus demo mode is rejected and that fixture
preparation happens before any port/process mutation.

- [ ] **Step 2: Run runtime tests and confirm RED**

Run:

```bash
PYTHONPATH=src python3 -m pytest -q \
  test/browser/test_platform_runtime.py
```

Expected: parser/config assertions fail because the demo flag is absent.

- [ ] **Step 3: Implement the opt-in demo preparation**

Add a boolean parser flag and retain explicit knowledge of whether project and
runtime roots were supplied. In demo mode, call the shared fixture builder,
select the dedicated default runtime root, replace only the sentinel-owned
project, and atomically write the bounded JSON input file. Point
`RuntimeConfig.project_root` at that project before constructing the
supervisor.

Print:

```text
Results demo inputs: <runtime>/results-visibility-inputs.json
Platform ready: http://127.0.0.1:5173
```

Do not auto-create a run or write SQLite outside normal API/worker paths.

- [ ] **Step 4: Document one-command use and exact stop/data semantics**

Add the command, input-file purpose, URL, `.local/results-visibility-demo/`
layout, prerequisites, Ctrl-C cleanup, and explicit statement that the fixture
does not alter the scientific workflow. Update deterministic worker E2E text
to include persisted artifact/QC/download visibility and remove the stale
statement that artifact extraction and QC UI are out of scope.

- [ ] **Step 5: Run focused tests and commit**

```bash
PYTHONPATH=src python3 -m pytest -q \
  test/browser/test_platform_runtime.py \
  test/browser/test_results_visibility_fixture.py
git add scripts/run_local_platform.py \
  test/browser/test_platform_runtime.py \
  docs/development/local-platform-runtime.md
git commit -m "feat: add results visibility demo mode"
```

### Task 5: Full gate, Draft PR, independent review, and conditional merge

**Files:**
- Review: complete `origin/main...HEAD`
- Modify only if a real in-scope failure requires it

- [ ] **Step 1: Run complete local verification**

Run focused fixture/runtime/worker tests, then:

```bash
PYTHONPATH=src python3 -m pytest -q
ENCODE_PIPELINE_TEST_REDIS_URL=redis://127.0.0.1:6391/15 \
  PYTHONPATH=src python3 -m pytest -q \
  test/workers/test_redis_process_integration.py \
  test/workers/test_cancellation_e2e.py \
  test/workers/test_tiny_execution_e2e.py
npm --prefix frontend test -- --run
npm --prefix frontend run typecheck
npm --prefix frontend run build
```

Expected: full suites pass; four real worker tests and both Playwright projects
execute with zero skips.

- [ ] **Step 2: Prove contract and formatting drift are zero**

Run OpenAPI export/Orval regeneration twice and require no additional diff.
Run changed-file Ruff check, Ruff format check, `git diff --check`, and the
no-hardcoded-paths guard. Confirm there is no OpenAPI/generated/package-lock
change because PR134 adds no API or dependency.

- [ ] **Step 3: Inspect screenshots and process cleanup**

Write local screenshots to a temporary directory, inspect all four viewport
captures, and require no overlap, clipping, empty panel, unusable action, or
horizontal document overflow. Confirm no Redis, worker, horse, Snakemake,
Vite, or helper process remains and every successful temporary runtime root is
deleted.

- [ ] **Step 4: Publish Draft PR134 and wait for CI**

Push `agent/pr134-results-visibility-gate`, create Draft PR134 against current
`origin/main`, and wait for fast-checks, browser-e2e, frontend, lint, and
lock-check to finish. Treat skips in the new browser/worker gates or flaky
failures as merge blockers.

- [ ] **Step 5: Request final independent read-only review**

Ask one Agency reviewer to inspect the exact PR HEAD for fixture containment,
scientific-workflow non-interference, canonical persistence, QC truthfulness,
download bytes, URL/mobile layout, cleanup, docs, and scope. Fix every
blocking/important finding, rerun affected gates, push, wait for green CI, and
request re-review of the new exact HEAD.

- [ ] **Step 6: Start the tryable demo and conditionally merge**

Start `python3 scripts/run_local_platform.py --results-visibility-demo` on
available local ports, verify readiness and the printed input file, and report
the live frontend URL while it remains usable. Stop it only after the requested
demo handoff or before ending if the environment cannot retain background
services safely.

Only when PR HEAD/local/remote match, base is latest `origin/main`, required CI
is green, real gates did not skip, independent review has no
blocking/important findings, worktree is clean, and diff is in scope: mark the
PR Ready and squash merge with `--match-head-commit <exact-head>`. Fetch and
verify the merge commit on `origin/main`, then mark the long-term Goal complete.

## Self-review

- Spec coverage: Tasks 1–2 cover one safe reusable fixture and runtime
  integration; Task 3 covers real durable indexing plus the exact browser
  chain, states, viewports, reload,
  bytes, and cleanup; Task 4 covers the one-command demo; Task 5 covers every
  required local/CI/review/merge condition.
- Scope: no API, generated client, dependency, persistence schema, parser,
  catalog, scientific target, lifecycle, input authoring, or Agent capability
  changes are planned.
- Type consistency: `ResultsVisibilityInputs`, four `*Config` runtime keys,
  `prepare_results_visibility_fixture`, fixture modes, and demo filenames are
  defined once and used consistently in later tasks.
