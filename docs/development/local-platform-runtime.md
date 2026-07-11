# Durable Local Platform Runtime

The local workflow-platform runtime stores canonical run lifecycle state in
file-backed SQLite. Redis and RQ provide only an out-of-process scheduling
boundary: queued payloads contain a stable run ID, while each worker reopens the
database and rebuilds the registry and execution dependencies for that run.
SQLite stores the canonical assignment, dispatch marker, and atomic worker claim.
RQ's unique job ID prevents duplicate scheduling while a reusable Redis job
record exists, but it is not the lifecycle truth or the final idempotency guard.
A job ID reused for another run is rejected, and a stale or unsuccessful RQ job
is reported instead of being treated as a successful enqueue.

## Run locally

Install the API extra, make Redis available, and start the FastAPI factory:

```bash
python -m pip install -e ".[api]"
uvicorn encode_pipeline.api.main:create_app --factory --reload
```

The API extra installs the platform runtime, Redis client, and RQ worker. It
does not install the scientific environment. Run the worker in an environment
where the `snakemake` executable is on `PATH`; the bundled workflow continues
to manage its rule-specific tools through its existing Snakemake configuration.

By default, the application uses:

```text
~/.encode-pipeline/platform.db
~/.encode-pipeline/workspaces/
redis://localhost:6379/0
RQ queue: encode-pipeline
```

The app factory applies the bundled Alembic migrations before accepting
requests. Run records, events, stdout/stderr chunks, and artifact references
survive a normal API process restart.

The API and worker must receive the same absolute database and workspace paths.
The Redis URL may include credentials, so it must not be copied into public run
events or issues:

```bash
export ENCODE_PIPELINE_DATABASE_URL="sqlite:////absolute/path/platform.db"
export ENCODE_PIPELINE_WORKSPACE_ROOT="/absolute/path/workspaces"
export ENCODE_PIPELINE_REDIS_URL="redis://localhost:6379/0"
export ENCODE_PIPELINE_QUEUE_NAME="encode-pipeline"
export ENCODE_PIPELINE_JOB_TIMEOUT_SECONDS="604800"

# Terminal 1
uvicorn encode_pipeline.api.main:create_app --factory

# Terminal 2, with the same environment
encode-worker
```

Only file-backed SQLite URLs are supported in this phase. An in-memory SQLite
database is deliberately rejected because it cannot preserve state across the
migration, API, and worker connections. Relative SQLite and workspace paths are
also rejected so processes with different working directories cannot silently
open different state.

`encode-worker --burst` processes the jobs currently in the configured queue
and exits. It is intended for integration tests and operational checks; the
normal local worker stays running.

`ENCODE_PIPELINE_JOB_TIMEOUT_SECONDS` is the workflow deadline enforced by the
worker's `ProcessRunner` (seven days by default); it is not the outer RQ timeout.
RQ sets its job timeout to that deadline plus a fixed 30-second cleanup grace.
The production `encode-worker` uses a phase-aware `DurableWorker`. Its
death-penalty adapter leaves RQ's configured exception unchanged and converts a
`JobTimeoutException` to `WorkerHardTimeout` only when the timeout signal is
actually delivered while the worker is inside the wrapped `job._execute` main
phase. `WorkerHardTimeout` derives directly from `BaseException`, not
`Exception`, so application `except Exception` handlers cannot accidentally
normalize or swallow that outer deadline. Job preparation, surrounding
`Job.perform` bookkeeping, and success, failure, or stopped callbacks remain
outside the main phase and retain RQ's native timeout behavior. Other RQ timeout
exception types are also left unchanged.

When the workflow deadline expires, that outer window lets the worker terminate
the directly spawned Snakemake process, perform a best-effort drain of output
that is immediately available, and persist the durable failure before RQ
applies its own timeout. The post-termination drain stops at the first of
256 KiB, eight non-blocking selector iterations, or 50 milliseconds; a
descendant that keeps a pipe open cannot extend that bounded drain.

If the outer RQ timeout is delivered, the timeout exception passes through the
process and execution layers. The worker makes a best-effort SQLite transition
to `failed` with reason code `WORKER_JOB_TIMEOUT`, then re-raises the original
`WorkerHardTimeout` so RQ also records the job failure. The database mapping can
itself fail, so this is not a substitute for later operational reconciliation.
The grace does not extend the workflow's execution budget, and it does not
guarantee termination of Snakemake descendants or a process group; that cleanup
boundary remains PR125.

RQ result metadata is retained for one day and failure metadata for seven days.
These are scheduling-retention settings only; durable run state and events
remain in SQLite.

## Submission and execution lifecycle

Preflight remains an explicit prerequisite. A successful preflight validates
the adapter inputs, materializes the per-run workspace, completes a Snakemake
dry-run, and leaves the run in `planned`. Execution is then submitted explicitly:

```bash
curl -X POST "http://127.0.0.1:8000/api/v1/runs/${RUN_ID}/start"
```

`POST /api/v1/runs/{run_id}/start` returns `202` only after the queue has
accepted the run's stable job identity and SQLite records a dispatch marker and
the `queued` transition. Runs that have not completed preflight return `409`.
An unavailable queue returns `503` with `RUN_QUEUE_UNAVAILABLE`; an identity or
state conflict returns `409` with `RUN_START_CONFLICT`. Public errors do not
include Redis URLs, credentials, command lines, or workspace paths.

The normal lifecycle is:

```text
PLANNED --start--> QUEUED --worker claim--> RUNNING --> SUCCEEDED | FAILED
```

The RQ payload contains only `run_id`, serialized with RQ's JSON serializer.
The worker reopens SQLite, rebuilds the default registry and service graph,
reconstructs the execution plan from the persisted inputs, and verifies that
the materialized workspace still matches that plan. It does not receive an
adapter, service, configuration object, or other process-local instance from
the API. After atomically claiming the SQLite assignment, it executes the
rebuilt Snakemake command and records the terminal state. A missing or modified
workspace, command-construction error, process startup error, timeout, non-zero
exit, or handled worker exception becomes a durable `failed` run with a stable,
public-safe reason code.

SQLite remains canonical throughout this sequence. RQ status can help an
operator inspect scheduling, but API reads, restart recovery, duplicate-worker
suppression, events, logs, and terminal outcomes all use SQLite state.

## Worker initialization failure repair

Runtime bootstrap, dependency composition, and adapter reconstruction can fail
before the normal execution service is usable. The worker then makes a
best-effort attempt to reopen the canonical SQLite database without depending
on the adapter registry. It mutates the run only when the persisted assignment
strictly matches the dequeued run ID, job ID, `rq` backend, and queue name. With
that ownership proof, the fallback can repair a missing dispatch marker,
converge a still-`planned` run through `queued`, and durably transition the run
to `failed` with reason code `WORKER_INITIALIZATION_FAILED`. Identity drift or
an unrelated lifecycle state is left unchanged.

This repair still depends on reaching the canonical database. If the configured
database URL is invalid, unavailable, or points somewhere that does not contain
the matching assignment, the worker cannot truthfully write canonical state.
RQ records the job failure, but SQLite may remain at its previous state and
requires operator investigation or a later reconciliation path. The worker does
not fabricate a recovery from Redis metadata or mark a different database as
authoritative.

## Submission failure windows and idempotency

Submission deliberately separates a durable reservation from confirmed
dispatch so each crash window has an explicit meaning:

- Before enqueue, SQLite may contain a `planned` run and an assignment without
  `dispatched_at`. This is a reservation, not worker ownership. It survives an
  API restart and a retry reuses the same job ID.
- If Redis cannot confirm enqueue, the API does not synthesize `queued`. The
  reservation normally remains `planned` and retryable. A concurrent worker may
  already have advanced SQLite, so clients should treat the returned run and a
  subsequent `GET` as authoritative rather than infer state from the HTTP error
  alone.
- If Redis accepted the job but the API stopped before persisting the marker or
  transition, the worker can persist `dispatched_at`, atomically move
  `planned` to `queued`, and then claim the assignment. Retrying start converges
  on the same RQ job identity.
- Once `queued` and `dispatched_at` are durable, an API restart preserves the
  run without asking Redis. Once the worker claim is durable, a later API
  restart also preserves `running`.
- Repeated start requests reuse the stable assignment and cannot append a
  second `queued` event. A duplicate worker invocation sees `claimed_at` and
  returns without executing Snakemake again. A backend job whose run ID,
  function, queue, or stable job ID conflicts with SQLite is rejected.

These rules cover API failure and retry races. A hard worker crash after an
atomic claim still leaves durable evidence that the worker owned the run, but
there is no heartbeat or lease reconciler yet; the API does not guess that the
claimed job is dead.

## Restart semantics

Restart recovery follows durable ownership instead of treating every active run
the same:

- `validating` is currently API-owned preflight work. An API restart marks it
  `failed` with `RUN_INTERRUPTED_BY_API_RESTART`.
- `queued` is worker-owned only when its RQ assignment has a persisted
  `dispatched_at` marker. `running` requires the stronger `claimed_at` marker,
  which is written atomically with the worker's first durable event. An API
  restart preserves those records without querying Redis.
- An assignment row alone is only a reservation. A `queued` or `running` run
  without the marker required for its state is an orphaned legacy or invalid
  record. It is marked `failed` with
  `RUN_ORPHANED_AFTER_API_RESTART`.
- `created`, `planned`, and terminal runs remain unchanged.

Recovery checks ownership and writes the failed state plus one
`run_recovered_after_restart` event in a single SQLite write transaction.
Re-running recovery is idempotent. A claim marker proves the hand-off occurred;
it does not prove that a worker is still alive. Worker crash/heartbeat
reconciliation is deliberately not guessed by the API process.

## Persistent process logs and environment isolation

The process runner drains stdout and stderr concurrently and appends bounded
line chunks to the run log tables while Snakemake is executing. Logs written
before a later process failure remain readable through
`GET /api/v1/runs/{run_id}/logs`; they are not stored in Redis or deferred until
the terminal transition. The default capture and persistence limit is
10,000,000 bytes per stream for one invocation. Additional output is drained to
avoid blocking the child, but is truncated with a structured warning. A log
persistence callback failure terminates the direct child and is mapped to a
durable execution failure.

Snakemake receives a reduced environment containing the runtime values needed
to locate the toolchain and trust store, such as `PATH`, selected Conda values,
locale, home, temporary-directory, and TLS certificate settings. Platform,
agent, cloud, source-control, and model-provider credentials are not inherited.
A `CommandSpec` is also rejected if it attempts to override a protected secret
name. In particular, the SQLite and Redis connection settings used by the API
and worker are not copied into the scientific subprocess environment.

## Deterministic worker E2E

`test/profiles/platform_worker_tiny/` is the platform execution profile. It has
one control-only sample while `use_control` is disabled, so the bundled
workflow resolves to a deterministic `all` job without invoking bioinformatics
tools or reading large scientific inputs. This is intentionally an engine and
lifecycle check, not a substitute for the scientific workflow tests.

`test/workers/test_tiny_execution_e2e.py` uses a real Redis queue, file-backed
SQLite, a separate `encode-worker --burst` process, and the real Snakemake
executable. It proves create → preflight → start, closes and reopens the API
before worker execution, observes `queued` survive restart, and then verifies
`queued` → `running` → `succeeded`, the durable dispatch and claim markers, RQ
completion, and persisted Snakemake log markers. Run it with:

```bash
ENCODE_PIPELINE_TEST_REDIS_URL=redis://localhost:6379/0 \
  PYTHONPATH=src python -m pytest test/workers/test_tiny_execution_e2e.py -v
```

The test skips when either the dedicated Redis test URL or `snakemake` on
`PATH` is unavailable. Real scientific data execution remains covered by the
existing Snakemake profiles and execution harness.

## Current cancellation boundary

The existing cancel route can durably move a non-terminal run to `cancelled`,
and worker completion checks will not overwrite that terminal state. In PR124,
however, cancellation does not cancel the RQ job, signal a running Snakemake
process, or terminate descendants. A process already running may continue and
append logs until it exits even though SQLite reports `cancelled`.

RQ job cancellation, real process cancellation, POSIX process-group creation
and termination, cancellation races, and orphan-process tests are explicitly
deferred to PR125. Authentication, multi-tenancy, HPC scheduling, artifact
extraction, and QC UI are also outside this local-execution milestone.
