# Run Recovery and Administrator Operations

Status: accepted for PR #172

## Context

SQLite is the canonical run lifecycle and result-metadata store. RQ records
prove that an execution was scheduled and may provide current operational
evidence, but they are neither a second lifecycle store nor sufficient by
themselves to manufacture a terminal outcome. The existing durable dispatch
and claim markers prove historical handoff, not continuing worker ownership.

An API restart therefore preserves a queued run with a dispatch marker and a
running run with a claim marker. That policy is intentionally safe, but it can
leave an operator-visible run active after an RQ job disappears, a worker exits
before its callback closes SQLite state, result indexing is interrupted, or
managed cleanup cannot be confirmed.

## Decision

HelixWeave provides an explicit, workflow-neutral administrator recovery
service and local `helixweave admin run` commands. It does not add an HTTP
mutation route or browser control before authentication and roles exist.

Diagnosis combines immutable SQLite lifecycle, assignment, workflow-build,
result-attempt, and event evidence with a bounded RQ observation. Output is a
stable, path-free projection. Redis connection details, worker names, process
identifiers, host names, workspace paths, inputs, environment values, command
lines, and raw exceptions are never returned.

Queue evidence distinguishes:

- exact schedulable ownership;
- an exact started job with a locally verifiable worker owner;
- an exact started job whose owner cannot be proven;
- an exact terminal job;
- a missing job;
- identity drift; and
- an unavailable queue.

No elapsed-time threshold is used to declare a worker dead. A missing job
after a durable claim is ambiguous and cannot authorize mutation.

Administrator mutations require the caller's exact run status and full
durable assignment identity. Repositories compare those values and all
monotonic markers under the SQLite write transaction before changing lifecycle
or appending the public-safe audit event.

An administrator fail is allowed only when no live owner remains:

- a dispatched, unclaimed queued assignment may be failed when its exact job
  is missing or terminal; and
- a claimed assignment may be failed only when its exact RQ job is terminal
  and configured managed-container cleanup succeeds.

Queue unavailability, identity drift, a live owner, an unproven started owner,
a missing claimed job, or pending cancellation intent refuses the operation.
The fail transaction changes only lifecycle and its audit event. It cannot
create success, artifact, QC, publication, or scientific evidence.

A safe requeue republishes only the same dispatched assignment that has never
been claimed. It does not create a new execution attempt, clear a marker, or
move a running or terminal run back to queued. Two nullable, monotonic fields
on `run_execution_assignments` record the one-time administrator requeue
request and queue confirmation. The request and its event commit atomically
before Redis is mutated. An exact terminal RQ record may be removed only after
the request is durable, then the same identity is enqueued. Confirmation is
recorded after RQ accepts or already contains that identity. A crash can retry
the durable request and stable job ID while the assignment remains unclaimed.
The replacement job carries a path-free marker bound to the request timestamp;
that marker is queue evidence, not lifecycle truth. At every exact worker entry
or terminal callback, the worker confirms a still-pending request in SQLite
before build admission, cleanup, claim, or terminal failure. This preserves
delivery proof after the RQ record expires. An explicit administrator failure
event instead abandons a pending request without falsely recording queue
acceptance.
The SQLite claim is the final execution guard: any worker that won a race has
already made the assignment ineligible, and a duplicate delivery exits before
scientific execution. A second recovery request after confirmation is a stable
conflict; the operator must fail the run and create a new run instead of cycling
lifecycle state.

Doctor visibility aggregates diagnostic gaps for the database, queue,
terminal callback, result indexing, and cleanup. It is read-only and
machine-readable. Missing configuration is distinct from unavailable or
inconsistent state.

## Consequences

- SQLite remains the only source of lifecycle truth.
- Dispatch, claim, cancellation, and requeue markers remain monotonic.
- RQ and process observations only narrow what an administrator may do.
- Requeue cannot duplicate a claimed execution.
- Result indexing and cleanup gaps are visible without being reclassified as
  scientific success or failure.
- The assignment schema, repositories, recovery service, RQ adapter, and
  migration enter the Bulk RNA-seq execution implementation closure and make
  the prior protected qualification stale.
- Automatic leases, heartbeat reconciliation, background repair, multi-node
  scheduling, deployment management, authentication, and a browser admin
  surface remain out of scope.
