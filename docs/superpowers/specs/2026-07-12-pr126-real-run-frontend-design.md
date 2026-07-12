# PR126 Real Run Frontend Design

## Objective

Complete Milestone B's browser-facing durable local execution path without
changing the PR123–125 lifecycle contract. The product path is validation,
durable run creation, preflight, explicit execution, progress/log observation,
truthful cancellation, and terminal recovery from the canonical SQLite record.

This PR does not add authentication, multi-tenancy, HPC, artifact extraction,
QC dashboards, another adapter, workflow-source bundles, leases/reconciliation,
or Agent write operations.

## Considered approaches

### Extend the existing run-progress boundary

Keep React Router, the existing run route, the generated-client adapter, and the
run-progress components. Add TanStack Query mutations, safe issue rendering,
bounded pagination/polling, official Playwright, and a small process supervisor.

This is the selected approach because it preserves established boundaries and
keeps the PR focused on the missing product and integration paths.

### Replace the run page with a new frontend state machine

A new route-level controller could separate every server resource into its own
query and model transitions explicitly. It would also duplicate working UI and
increase the risk of a broad frontend rewrite. This is rejected for PR126.

### Add streaming APIs or a full Compose stack

SSE/WebSockets would require a new backend contract. Full Compose would require
new API, worker, frontend, and scientific images; the existing runner image only
starts Snakemake. Both expand the architecture beyond the final Milestone B PR
and are rejected.

## Frontend architecture

`RunApiClient` gains `startRun(runId)`. The generated-client adapter invokes the
existing Orval-generated `startRun` function. Components continue to consume
the stable DTO adapter and never import SQLAlchemy, RQ, or backend internals.

The run URL and SQLite-backed GET response are the recovery boundary. TanStack
Query owns the run snapshot, events, and stdout/stderr. Start, preflight, and
cancel use `useMutation`; responses seed canonical server records and targeted
query invalidation refreshes related resources. The UI never synthesizes a
QUEUED, CANCELLED, or CANCELLING state.

Run creation seeds the returned CREATED record and navigates immediately to
`/runs/:runId` with one-use route state requesting preflight. The detail page
consumes that state once and triggers the mutation. A hard reload does not infer
permission to mutate a CREATED run; it presents a retryable Run preflight action.
Thus a failed or lost preflight request never loses the durable run URL.

PLANNED alone enables Start run. Cancel is enabled only for known active states:
CREATED, VALIDATING, PLANNED, QUEUED, and RUNNING. Unknown statuses receive a
neutral badge, stop automatic polling, disable lifecycle mutations, and retain
manual Refresh.

## Cancellation truthfulness

A successful RUNNING cancel may return `202` with a RUNNING record. The cache
keeps that record and the page displays “Cancellation requested” while polling.
The same hint is reconstructed from a `cancellation_requested` event after
reload. Because persisted intent does not prove Redis delivery, RUNNING retains
a Retry cancellation action. Only a backend record whose status is CANCELLED is
rendered as cancelled. A `200 + CANCELLED` callback race is accepted as already
acknowledged. A `409` or `503` preserves the current run and renders safe,
structured issues without fabricating a terminal state.

## Query bounds and errors

Automatic polling runs at a fixed interval only for CREATED, VALIDATING, QUEUED,
and RUNNING, stops for terminal, PLANNED, and unknown states, and pauses after a
fixed observation window. Manual Refresh remains available.

Events and each log stream follow cursors for at most five fixed-size pages.
Repeated cursors terminate collection and the UI reports truncation instead of
looping or loading without bounds.

The fetcher preserves a shape-validated public issue envelope on non-2xx
responses. The adapter maps it into the existing UI Issue type. The run page
shows code, message, source, and hint only. It never displays technical_message,
raw JavaScript exception strings, arbitrary event context, Redis details,
commands, or filesystem paths. RunRecord.error and RunEvent.issue use the same
safe view. Execution stdout/stderr remain visible as explicit run logs.

## Responsive behavior

Existing Button and lucide components are reused. Lifecycle actions have
accessible names, pending labels, and a shared mutation gate. Metadata collapses
to one column on narrow screens, run IDs can wrap, and log containers constrain
horizontal overflow. Desktop and 390×844 mobile Playwright flows assert viewport
fit and operable controls.

## Real browser E2E

Official `@playwright/test` runs Chromium with one worker. A Python runtime
harness creates a controlled project, temporary file-backed SQLite database,
unique RQ queue, workspace, and marker directory. It starts real FastAPI, a
continuous independent DurableWorker, and Vite against a real Redis service.

The controlled workflow copies the current package source and uses a tiny
Snakefile. A valid `threads=1` config completes quickly with deterministic log
markers. A valid `threads=2` config starts a long helper and normal child process
for cancellation. No API is mocked and no scientific dataset is used.

The desktop flow proves validate → create and URL retention → preflight →
PLANNED → start → durable terminal success → logs → reload recovery. The mobile
flow waits for RUNNING and the long-task log, captures the real cancel response,
proves `202 + RUNNING` is displayed as requested rather than cancelled, then
waits for canonical CANCELLED. It rejects any CANCELLING label.

The harness checks Redis, API, worker, and frontend readiness. Each child starts
in an isolated session. SIGINT, SIGTERM, normal completion, and test failure all
perform bounded TERM/KILL cleanup across worker horse and descendant process
groups. The Redis queue/job is unique and removed without flushing unrelated
data.

## Local demo

`python3 scripts/run_local_platform.py` is the documented one-command launcher.
It starts a local redis-server, FastAPI, a continuous worker, and Vite using the
same absolute SQLite/workspace environment. Defaults are Redis 6379, API 8000,
frontend 5173, and `.local/platform-demo/` for durable data and logs. Ports and
data root are configurable. Readiness failures terminate the full process tree;
Ctrl-C is the stop command.

The launcher requires the API package, `snakemake`, `redis-server`, Node/npm,
and installed frontend dependencies. Vite continues to proxy `/api`; no broad
CORS policy is introduced.

## CI and acceptance

A required browser-e2e job uses Redis 7.4.9, the existing ci-fast micromamba
environment, Node 20, the official Chromium install, and the real runtime
harness. It cannot silently skip. Failure artifacts contain Playwright output,
not platform credentials.

Acceptance also requires focused and full Python tests, existing real
Redis/RQ/Snakemake tests, frontend unit tests, typecheck, production build,
OpenAPI/Orval drift, changed-file Ruff and format checks, `git diff --check`,
process-residue checks, and independent read-only merge-gate review.
