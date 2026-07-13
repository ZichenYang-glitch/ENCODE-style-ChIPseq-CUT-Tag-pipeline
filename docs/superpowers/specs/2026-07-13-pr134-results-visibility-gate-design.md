# PR134 Results Visibility Gate Design

## Objective and boundary

PR134 closes the results-visibility sub-milestone with one deterministic,
non-mocked browser path:

```text
validate -> create -> preflight -> start -> RQ worker -> Snakemake
  -> SUCCEEDED -> QC metrics -> source artifact -> safe download
```

The run, events, artifact references, and QC metrics remain canonical SQLite
state. Redis/RQ still carries only scheduling data. The browser reads the
existing PR128/PR131 APIs and downloads through the PR133 endpoint. This PR
does not change an API, generated client, scientific Snakefile, manifest
generator, artifact catalog, QC parser, Stage 51, or lifecycle semantics.

## Considered approaches

### Extend only the inline Playwright project

The current browser runtime could write one QC TSV in its embedded helper.
That is the smallest diff, but the one-command demo would exercise a different
project and the controlled output contract would remain hidden in a large test
harness. This is rejected.

### Reuse one controlled results fixture in E2E and the demo

A focused Python fixture builder creates a sentinel-owned project containing
the current adapter package, artifact inventory, default profile, a tiny
Snakefile, and deterministic helper scripts. Playwright and an explicit local
demo flag use the same builder. This is selected because it gives tests and
humans the same workflow source while leaving production scientific targets
untouched.

### Modify the bundled scientific tiny profile or rule targets

Enabling the real QC rules would require bioinformatics inputs and tools, or a
test-only branch inside scientific rules. Either weakens the scientific/test
boundary and risks changing target selection. This is rejected.

## Controlled fixture contract

`scripts/results_visibility_fixture.py` owns creation of the controlled
project. It accepts a dedicated project directory and refuses to replace an
existing directory unless an exact fixture sentinel proves ownership. It
never removes the surrounding runtime root, SQLite database, workspace tree,
or an unrelated path.

The generated project has only the files already required by the worker build
identity boundary:

- current `src/encode_pipeline` source;
- `pyproject.toml` and the artifact inventory;
- a minimal default Snakemake profile;
- the existing deterministic sample sheet;
- a tiny `workflow/Snakefile` and task helper.

The helper has four configuration-selected modes. `threads` is used only as a
fixture selector inside this controlled Snakefile:

| Mode | threads | terminal state | artifact/QC outcome |
| --- | ---: | --- | --- |
| results | 1 | SUCCEEDED | manifest + valid `C1.qc_summary.tsv`, eight metrics |
| cancel | 2 | CANCELLED | no completion, manifest, or QC output |
| empty | 3 | SUCCEEDED | manifest only, confirmed zero QC metrics |
| malformed | 4 | SUCCEEDED | manifest + malformed QC source, indexing failed |

The valid source uses the existing exact ENCODE `qc_summary` header and one
`C1`/`chipseq` row. It supplies deterministic values for total reads, FRiP,
peak count, percent duplication, estimated library size, NRF, PBC1, and PBC2.
All other columns are `NA`. These produce eight persisted metrics without
inventing thresholds or QC flags.

The manifest is a truthful machine-readable fixture output. Results and
malformed modes record the QC summary and manifest as present; empty records
only the manifest. The cancellation helper records its process group and
normal child, then sleeps without writing any result. No large scientific
input is read and no bioinformatics executable is invoked.

## Real browser gate

The existing official Playwright harness continues to start real FastAPI,
Redis, one independent `DurableWorker`, Vite, and Snakemake against an isolated
SQLite/workspace root and unique queue. API calls are not mocked.

The main results test performs the product sequence in this order:

1. Validate inputs, create the durable run URL, complete preflight, and start.
2. Wait for canonical `SUCCEEDED` and persisted execution logs.
3. Open `?view=qc`; require the non-empty indexed state, exact decimal values,
   sample/assay data, and the persisted source-artifact action.
4. Reload the QC deep link and verify the same database-backed metrics.
5. Follow the source action to
   `?view=artifacts&artifact=<artifact_id>` and require the QC summary relative
   path and metadata.
6. Download that artifact through the real binary endpoint and compare the
   filename and complete bytes with the fixture source.
7. Reload the artifact deep link and recheck persisted detail.

The same browser test creates one empty and one malformed run. It requires a
confirmed empty message only after `qc_metrics_indexed(metric_count=0)`, and a
redacted indexing-failed message while the canonical run remains SUCCEEDED.
These paths do not read raw source content through an API.

The existing cancellation test remains a real mobile path. A 202 cancel keeps
the visible state RUNNING until worker acknowledgement; final CANCELLED is
required and the helper and child process must be gone.

## Layout, screenshots, and URL recovery

The E2E does not redesign the workbench. It exercises the existing Activity,
QC, and Artifacts tabs at 1440x900, 1024x768, 390x844, and 360x800. QC rows,
artifact list/inspector, source action, download button, and shared run actions
must remain visible, non-empty, and operable without document horizontal
overflow.

Screenshots are attached to Playwright results and can also be written to an
explicit local directory. Bounding-box assertions verify desktop side-by-side
artifact layout and mobile stacked layout. Reload assertions verify both QC
and selected-artifact URLs rather than relying on in-memory React state.

## One-command local demo

`python3 scripts/run_local_platform.py --results-visibility-demo` prepares the
same controlled project under a dedicated runtime root and then starts the
existing supervised Redis/API/worker/Vite stack. It writes a bounded public
input manifest containing the success configuration and sample-sheet path, and
prints that path together with the normal frontend URL. The user still performs
validate/create/preflight/start in the UI, so the demo does not seed or forge
lifecycle state.

Unless the caller supplies another runtime root, this mode uses
`.local/results-visibility-demo/`; it does not mix with the ordinary
`.local/platform-demo/` database. Repeated starts may replace only the
sentinel-owned controlled project. SQLite and workspaces remain durable until
the user removes that demo root. One Ctrl-C retains the existing bounded
process-session cleanup.

The ordinary `python3 scripts/run_local_platform.py` behavior remains
unchanged. Prerequisites remain the API extra, locked frontend packages,
Snakemake, and local Redis or `redis-server`.

## Failure and cleanup behavior

Fixture preparation fails before services start if the target is unsafe,
already exists without the sentinel, or cannot be built completely. Partial
fixture directories are removed only when ownership has already been
established.

Every Playwright invocation keeps its random runtime root, owner token, unique
queue, and process-session evidence. Normal completion and failure both stop
API, worker, work horse, Snakemake, helpers, Vite, and launcher-owned Redis.
The runtime wrapper deletes its temporary root only after successful teardown;
failure artifacts remain available for diagnosis. Cleanup never flushes a
shared Redis database or deletes an unowned path.

## Tests and acceptance

Focused Python tests cover fixture ownership, deterministic files/inputs,
valid/empty/malformed modes, and demo argument/runtime isolation. The real
worker gate verifies that the controlled success run yields durable artifacts
and eight SQLite QC metrics after reopening persistence.

Playwright proves non-empty QC, source navigation, artifact metadata, exact
download bytes, deep-link reload, confirmed empty, redacted indexing failure,
mobile cancellation, four viewport sizes, and runtime/process cleanup. The
required browser CI job executes these tests without skip.

Acceptance also requires the full Python suite, the real
Redis/RQ/SIGALRM/cancellation/tiny gate, frontend tests, typecheck, production
build, OpenAPI/Orval zero drift, changed-file Ruff/format checks,
`git diff --check`, a clean worktree, green GitHub checks, and an independent
read-only merge-gate review.

## Residual risks and non-goals

The fixture proves platform integration, not scientific correctness or
real-data performance; existing scientific tests retain that responsibility.
Artifact downloads still have the documented no-checksum same-metadata
replacement risk and browser Blob memory cost. Authentication, multi-tenancy,
HPC, object storage, immutable workflow bundles, reconciler/leases, thresholds,
automatic scientific conclusions, input authoring, a second adapter, and Agent
write actions remain outside this PR.
