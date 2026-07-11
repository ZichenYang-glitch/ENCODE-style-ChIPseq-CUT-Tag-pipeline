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

Execution jobs have a configurable timeout (seven days by default). RQ result
metadata is retained for one day and failure metadata for seven days. These are
scheduling-retention settings only; durable run state and events remain in
SQLite.

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
reconciliation is not implemented in PR123 and is deliberately not guessed by
the API process.

## Current boundary

The API's preflight remains the existing local background dry-run path. The RQ
worker currently validates its durable job identity, rebuilds dependencies from
SQLite and the default registry, and atomically claims the assignment while
recording one handshake event. A duplicate worker invocation observes the
existing claim and returns without recording a second handshake. No public
route submits execution yet. Real Snakemake execution, persistent streaming
logs, worker failure reconciliation, process cancellation, authentication, HPC
scheduling, and artifact extraction remain later milestones.
