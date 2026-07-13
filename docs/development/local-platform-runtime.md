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

### One-command demo stack

Create the repository-local locked environment once. The environment lives
under the ignored `.local/` directory and is never part of a commit:

```bash
micromamba create -p .local/envs/ci-fast --file workflow/envs/ci-fast.lock
./.local/envs/ci-fast/bin/python -m pip install -e ".[api,dev]"
npm --prefix frontend ci
export PATH="$PWD/.local/envs/ci-fast/bin:$PATH"
python scripts/run_local_platform.py --doctor
```

The doctor checks Python 3.12, the API/test imports, Snakemake 8.30.0, Redis
server 7 or newer, Node.js 20 or newer, npm, and the locked frontend install.
It exits before opening ports or creating runtime data. Failures name the
missing prerequisite without printing environment variables or private paths.

After the doctor succeeds, one foreground command starts the complete local
stack:

```bash
python scripts/run_local_platform.py
```

Open `http://127.0.0.1:5173`. The supervisor starts or reuses local Redis,
FastAPI on port 8000, one independent RQ worker, and Vite on port 5173. It
uses one shared file-backed database and workspace tree:

```text
.local/platform-demo/platform.db
.local/platform-demo/workspaces/
.local/platform-demo/logs/{redis,api,worker,frontend}.log
```

Press Ctrl-C once to terminate the API, worker, any RQ work horse and normal
descendants, Vite, and a Redis process that the launcher started. An already
running Redis at `redis://127.0.0.1:6379/0` is reused and is not stopped or
flushed. The default queue is `encode-pipeline-demo`; Redis remains scheduling
metadata only, while the SQLite file remains canonical.

The launcher deliberately requires the same locked toolchain as the doctor.
It normally verifies the local `redis-server` binary before startup. When that
binary is absent, an explicitly configured, already-running Redis is accepted
only after a bounded server-version check confirms Redis 7 or newer; otherwise
startup fails before creating runtime data. The launcher performs bounded
readiness checks for Redis, the API, the registered worker, and the frontend;
startup failures name the failed service and point to its local log without
printing Redis credentials.

Use explicit flags for an isolated demo or non-default ports:

```bash
python scripts/run_local_platform.py \
  --runtime-root /tmp/encode-platform-demo \
  --redis-url redis://127.0.0.1:6380/0 \
  --queue-name encode-platform-demo \
  --api-port 8010 \
  --frontend-port 4173
```

Vite proxies `/api` to the configured API port. The launcher sets the narrow
`VITE_API_PROXY_TARGET` value for that process; it does not enable broad CORS.

### One-command results visibility demo

An opt-in deterministic project exercises the full results path without
scientific data or bioinformatics tools:

```bash
python scripts/run_local_platform.py --results-visibility-demo
```

Open the printed `http://127.0.0.1:5173` URL. The launcher also prints the path
to `.local/results-visibility-demo/results-visibility-inputs.json`. In the
workflow page, choose **Author inputs**, paste the JSON `resultsConfig` into
advanced YAML mode (JSON is valid YAML), and import the TSV at `samplesPath`.
In Review, use Validate inputs and Create run; on the durable run URL, wait for
preflight and then use Start run. The real worker and Snakemake process produce
eight persisted QC metrics, a sample QC summary artifact, and a result
manifest. The QC source action opens the persisted artifact detail, whose
Download action returns the exact TSV.

This mode uses a separate durable root:

```text
.local/results-visibility-demo/platform.db
.local/results-visibility-demo/workspaces/
.local/results-visibility-demo/results-visibility-project/
.local/results-visibility-demo/results-visibility-inputs.json
.local/results-visibility-demo/logs/{redis,api,worker,frontend}.log
```

Its default RQ queue is `encode-pipeline-results-demo`, distinct from the
ordinary `encode-pipeline-demo` queue even when both use the same local Redis.

The controlled project is replaced on restart only when its ownership sentinel
matches. SQLite and workspaces are retained. An unowned directory is never
deleted. The ordinary launcher still uses the repository's scientific workflow
and `.local/platform-demo/`; the demo flag does not add a test branch to the
scientific Snakefile, manifest generator, artifact catalog, or QC parser.

Press Ctrl-C once to stop the demo stack and all worker/Snakemake descendants.
The prerequisites and configurable Redis/API/frontend ports are the same as
for the ordinary launcher.

### Manual processes

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
export ENCODE_PIPELINE_REDIS_CONNECT_TIMEOUT_SECONDS="2"
export ENCODE_PIPELINE_REDIS_API_READ_TIMEOUT_SECONDS="5"
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

The Redis timeout values must be positive finite numbers. The API queue client
uses the two-second connection timeout and five-second command read timeout by
default. The synchronous start and running-cancel routes run in FastAPI's worker
threadpool, so a slow or unavailable Redis endpoint cannot block the server
event loop. A start timeout returns the public-safe
`503 RUN_QUEUE_UNAVAILABLE` response; a stop timeout returns the equally
sanitized, retryable `503 RUN_CANCELLATION_UNAVAILABLE`. The worker connection
uses the same finite
connection timeout but deliberately does not inherit the API read timeout. RQ
sets the worker socket read timeout from its longer blocking dequeue interval,
so an idle worker can wait for work without being disconnected every five
seconds.

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

The runtime dependency is restricted to the RQ 2.10 minor series
(`rq>=2.10,<2.11`) because the phase boundary currently relies on the private
`Job._execute` hook. A compatibility canary verifies that exact dispatch
contract plus RQ's public stop-command signature, stopped-callback signature,
work-horse process-group setup, and `killpg` behavior. CI runs real Redis
independent-process SIGALRM and process-cancellation tests.
Any RQ minor upgrade must deliberately update that guard and pass the real
signal test before the supported range is widened.

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
The grace does not extend the workflow's execution budget. User cancellation is
a separate RQ stop-command path described below; it terminates the RQ horse
process group rather than relying on `ProcessRunner`'s direct-child cleanup.

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
An unavailable queue returns `503` with `RUN_QUEUE_UNAVAILABLE`; a missing
workflow build identity returns `409` with
`RUN_WORKFLOW_BUILD_IDENTITY_MISSING`; and a queue/job identity or lifecycle
conflict returns `409` with `RUN_START_CONFLICT`. Public errors do not include
Redis URLs, credentials, command lines, or workspace paths.

The normal lifecycle is:

```text
PLANNED --start--> QUEUED --worker claim--> RUNNING --> SUCCEEDED | FAILED
                                                   \--stop acknowledged--> CANCELLED
```

### Durable workflow build identity

A run planned by the current preflight path has both a successful dry-run and a
one-to-one workflow build identity in SQLite. The identity records the workflow
ID and adapter version, the logical `workflow/Snakefile` entrypoint, a versioned
hashing scheme, a SHA-256 digest, and its capture time. The digest is computed
from logical relative paths and file contents, never absolute paths,
modification times, or workspace state. Its controlled manifest covers the
Python adapter/platform package, workflow source and lock/config files, the
default Snakemake profile, top-level execution scripts, and `pyproject.toml`.
Missing files, non-regular files, and symlinks fail closed.

Preflight captures the identity before planning and again after the Snakemake
dry-run. Source drift between those captures fails the run with
`PREFLIGHT_WORKFLOW_BUILD_CHANGED`. On success, the `validating -> planned`
transition, identity row, and final preflight event are committed in one SQLite
transaction. Existing databases are deliberately not backfilled: a legacy
`planned` run without an identity is rejected by the start endpoint with
`409 RUN_WORKFLOW_BUILD_IDENTITY_MISSING`.

After a worker has validated the durable run/job assignment and converged the
run to `queued`, but before it claims the assignment, enters `running`, or
starts a process, it fingerprints its reconstructed local source tree. A
missing identity, unreadable local source, or mismatch becomes a durable
`failed` run with `RUN_WORKFLOW_BUILD_IDENTITY_MISSING`,
`RUN_WORKFLOW_BUILD_IDENTITY_UNAVAILABLE`, or
`RUN_WORKFLOW_BUILD_IDENTITY_MISMATCH`, respectively. The assignment remains
unclaimed and Snakemake is not started. The API and worker compose the command
builder and identity provider from the same explicit project source root, so
the verified logical entrypoint is the one used to build the command.

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

`test/workers/test_cancellation_e2e.py` creates a deterministic, long-running
tiny Snakefile under a temporary controlled project root. It uses the same real
Redis, real Snakemake executable, file-backed SQLite, and independent
`DurableWorker` boundary. The test records the Snakemake/helper/child PIDs and
shared RQ horse process group, requests cancellation through HTTP, then proves
that the group is gone, no completion marker exists, logs remain persisted,
SQLite is `cancelled`, and RQ reports `STOPPED` in its failed-job registry. Both
tests run in the pull-request `fast-checks` job, where Redis and Snakemake are
required rather than optional.

Ordinary local full-suite runs skip the real-process tests when the dedicated
Redis test URL is absent. Once that URL enables the cancellation gate,
`snakemake` and POSIX process groups are mandatory and a missing prerequisite is
a test failure. Real scientific data execution remains covered by the existing
Snakemake profiles and execution harness.

The required browser E2E uses the same controlled results project as the demo,
but an invocation-owned temporary database, workspace, queue, and service
sessions. It proves real validate → create → preflight → start → SUCCEEDED,
non-empty QC persistence, QC-source artifact navigation, exact safe download,
confirmed empty and redacted indexing-failure states, deep-link reload, mobile
cancellation, and bounded cleanup. This remains a platform integration gate;
scientific correctness remains in the existing workflow tests.

## Running cancellation semantics

The cancel route still moves `created`, `validating`, `planned`, or `queued`
runs immediately to `cancelled`. A queued RQ job may still dequeue, but its
worker rechecks SQLite and exits without starting `ProcessRunner` or Snakemake.

A running run has no synthetic `cancelling` status. Cancellation proceeds in
two durable phases:

1. SQLite atomically stores `cancellation_requested_at` and the first
   `cancellation_reason` on the claimed execution assignment, together with one
   `cancellation_requested` event. The run remains `running`; `ended_at` and the
   public run record's `cancellation_reason` remain unset.
2. The synchronous API route strictly matches the run/job/backend/queue and RQ
   job function, arguments, origin, started status, worker identity, and
   registered module-level stopped callback. It then calls RQ 2.10's public
   `send_stop_job_command`. A successful publish returns HTTP `202` while
   SQLite remains active; if acknowledgement or natural completion commits
   before the response is assembled, the route instead returns the canonical
   terminal snapshot with HTTP `200`. Publication itself is never treated as
   process-termination acknowledgement.
3. RQ kills and waits for the work horse process group. Snakemake and its normal
   descendants stay in that group because `ProcessRunner` does not create a new
   session. Only the parent's `on_stopped(job, connection)` callback can
   atomically store `cancellation_acknowledged_at`, move the expected
   `running` row to `cancelled`, set `ended_at` and the persisted reason, and
   append the single terminal event.

Repeated HTTP requests preserve the first intent and event but resend the stop
command while the run remains active, because Redis Pub/Sub delivery is not a
durable acknowledgement. Repeated callbacks are no-ops. Natural success or
failure and stop acknowledgement race through SQLite's expected-status write;
whichever terminal commit wins is preserved, so no second terminal event can be
created. A valid stopped callback without user cancellation intent atomically
moves an exactly matched `planned`, `queued`, or `running` execution to `failed`
with reason `WORKER_STOP_WITHOUT_CANCELLATION`, rather than impersonating a user
cancellation or leaving an RQ-stopped job active in SQLite.

An assignment/configuration conflict returns `409 RUN_CANCELLATION_CONFLICT`
without sending a stop. Redis, RQ, OS, missing-job, stale-job, or callback
identity failures return sanitized, retryable
`503 RUN_CANCELLATION_UNAVAILABLE`; an intent may remain durable for a later
retry, but SQLite is never changed to `cancelled` on that error. A pre-PR125 job
without the stopped callback is deliberately not killed.

After a successful RQ stop, the scheduling metadata is `STOPPED` and the job is
present in RQ's failed-job registry; API lifecycle reads still use SQLite. If
the stopped callback cannot reach SQLite after the process group is gone, the
run can remain `running` with a durable intent and requires operator diagnosis.
Heartbeat/lease reconciliation is intentionally not invented here.

Authentication, multi-tenancy, HPC scheduling, object storage, immutable
workflow bundles, automatic QC thresholds/conclusions, and a second adapter
remain outside this local runtime.
