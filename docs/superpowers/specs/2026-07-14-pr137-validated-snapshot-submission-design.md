# PR137 Validated Snapshot Submission Design

## Objective and boundary

PR137 turns successful backend validation into a durable, immutable input
authorization that can be consumed to create exactly one canonical run. The
browser workbench sends its current Config/Samples/Options draft to the real
adapter validator, receives an opaque snapshot only on success, and creates a
run by sending that snapshot ID rather than resending mutable inputs.

The PR does not move validation into RQ, add a lifecycle state, change
preflight/execution, persist browser drafts, or alter the scientific validator.
PR138 owns the complete real browser input-to-results gate.

## Local development prerequisite

The repository's `workflow/envs/ci-fast.lock` is the reproducible local Python
and Snakemake base. A project-local environment may live under the ignored
`.local/envs/ci-fast` directory; no interpreter path is committed. The local
platform launcher performs a side-effect-free doctor check before probing
ports, creating runtime directories, or starting processes. It requires Python
3.12, Snakemake 8.30.0, Redis server 7 or newer, Node 20 or newer, npm, the API
and development Python packages, and locked frontend dependencies. Failures
name the missing requirement and a remediation action without printing PATH,
environment variables, private filesystem locations, or subprocess output.
An existing Redis is accepted only when `redis-server` is absent and a bounded,
sanitized version probe confirms Redis 7 or newer. API and worker readiness
also precede Vite startup so the browser readiness URL cannot race the backend.

## Durable snapshot model

`ValidatedInputSnapshot` is a workflow-neutral domain value stored behind the
existing `RunRepository` boundary. It contains:

- an opaque random `vsnap_<uuid>` identifier;
- workflow ID, adapter version, schema version, and schema dialect;
- the complete workflow build identity captured around validation;
- one canonical JSON text containing Config/Samples/Options;
- a versioned framed SHA-256 payload digest;
- successful adapter-validation evidence: outcome, bounded issue codes, and
  validation time;
- a fixed first-use expiry; and
- optional atomic consumption evidence: run ID and consumption time.

Canonical JSON uses sorted keys, compact UTF-8 JSON, finite numbers, and only
JavaScript-safe integers because the public OpenAPI client represents JSON
numbers as TypeScript numbers. The canonical text is re-parsed, reserialized,
and rehashed at every repository write/read boundary. This catches malformed or
tampered durable data before it can become a run. Domain values expose fresh
`WorkflowInputs` copies rather than mutable ORM rows.

The initial expiry is 30 minutes. Expiry is checked for first consumption only.
Once a snapshot has produced a run, an identical retry returns that same run
even after expiry or a later deployment; this is required to recover from a
lost create response. A consumed snapshot can never create another run.

## Validation and build identity

`ValidatedInputService` composes the existing `ValidationService`, registry,
build identity provider, and repository. It captures build identity immediately
before reading the adapter schema contract and again after adapter validation,
so both the contract version and validation evidence are inside the same stable
build window. A missing identity, unavailable schema contract, or identity
mismatch fails closed and creates no snapshot. On success, the service persists
the snapshot only after the adapter result, schema contract, canonical payload,
digest, and validation evidence are complete.

Run creation checks the path workflow ID against the snapshot and recaptures
the current build immediately before first consumption. A missing or changed
identity fails closed before a run row exists. The later preflight build
identity remains independently authoritative for execution, so a change after
creation still cannot cross the existing PLANNED execution gate.

The snapshot stores source identity but does not package an immutable wheel or
workflow bundle. That existing deployment limitation remains explicit.

## Atomic create, replay, and concurrency semantics

The repository exposes one atomic consume-and-create operation. Under the
InMemory lock or SQLite `BEGIN IMMEDIATE`, it:

1. reloads and validates the snapshot;
2. verifies workflow and expected build identity;
3. if already consumed, validates the linked run, requires its creation time to
   equal the pre-expiry snapshot consumption time, checks requested tags, and
   returns that canonical run as an idempotent replay;
4. otherwise checks expiry, creates the CREATED run and its single initial
   event, and writes consumption evidence in the same transaction.

Concurrent identical requests converge on one run. A replay with different
tags is a conflict; inputs cannot differ because create accepts no inputs.
Snapshot IDs are single-use authorization tokens, not deterministic payload
IDs: separately successful validations produce separate snapshots.

Public failure semantics are stable and disclosure-safe:

- missing or cross-workflow snapshot: 404 `VALIDATED_SNAPSHOT_NOT_FOUND`;
- expired before first use: 409 `VALIDATED_SNAPSHOT_EXPIRED`;
- stale workflow build: 409 `VALIDATED_SNAPSHOT_STALE`;
- different replay metadata: 409 `VALIDATED_SNAPSHOT_REPLAY_CONFLICT`;
- corrupted durable snapshot or linked run: 409
  `VALIDATED_SNAPSHOT_DATA_INVALID`; and
- unavailable identity/schema/persistence: controlled 503/500 envelopes
  without raw exception text.

No error distinguishes an unknown ID from a snapshot belonging to another
workflow. No snapshot payload, build source path, database address, or internal
exception is returned.

## HTTP and OpenAPI contract

`POST /api/v1/workflows/{workflow_id}/validate` retains operation ID
`validateWorkflow`. A successful response now includes a safe snapshot
projection: ID, workflow ID, schema/adapter versions, payload digest,
validated time, and expiry. Failed validation has `snapshot: null` and never
persists an authorization.

`POST /api/v1/workflows/{workflow_id}/runs` retains operation ID `createRun`
but its request is now only `{snapshot_id, tags}`. A first create returns 201;
an idempotent replay returns 200 with the same `RunResponse`. Raw config,
samples, options, or a client-provided `validated` flag are forbidden by
Pydantic `extra="forbid"`.

Routes remain thin. They call the snapshot services and project domain values
through Pydantic response models. ORM rows and the workspace are never accessed
by routes. OpenAPI and Orval files are generated mechanically.

## Frontend flow and stale-draft rule

The PR136 workbench remains the canonical in-memory draft owner. Review adds
generated-operation-backed TanStack mutations:

1. Validate the deterministic Review payload.
2. Retain the returned snapshot together with the draft's deterministic
   serialization.
3. Any semantic config, sample, or option mutation immediately discards that
   snapshot. Tab changes, Form/YAML view changes, and a no-op serialization do
   not.
4. Create sends only the retained snapshot ID and navigates to
   `/runs/{run_id}` on the first or replay response.
5. The run page performs existing preflight and retains its explicit Start
   action.

Validation issues remain backend facts and are rendered through controlled
Issue fields. A lost create response leaves the snapshot available for retry;
server idempotency returns the canonical run. A conflict or unknown response is
described as unconfirmed and never fabricates a run. Each create attempt also
carries the immutable workflow/snapshot or draft-revision authority that
started it. If inputs change before its response arrives, neither workbench
navigates nor starts preflight; it warns that the earlier request may have
created a canonical run. The legacy workflow-detail validator is updated to the
same snapshot contract and guards both validation and creation responses with
the current workflow/input authority, so a late response cannot restore or act
on a snapshot after an edit or navigation.

## Tests and residual limits

Domain, InMemory, SQLAlchemy, migration, service, API, OpenAPI, generated-client,
workbench, and run-progress tests cover tamper detection, cross-workflow use,
stale identity, expiry, replay, conflicting replay, restart recovery, and
concurrent create. Browser coverage in this PR proves real validate → snapshot
→ create → preflight navigation; PR138 owns the full execution/results path.

Residual limits are deliberate: snapshots are not user-scoped because auth and
multi-tenancy are out of scope; no cleanup/reconciler deletes expired unused
snapshots; and a workflow source change after create but before preflight is
rejected by the existing preflight build gate rather than an immutable bundle.
