# Run Lifecycle and Progress Boundary Design

> **Status:** design review
> **Scope:** docs-only boundary and API contract for workflow run lifecycle and progress tracking
> **Date:** 2026-07-03

This document defines the next platform phase after the validation MVP and the read-only agent layer: a workflow run lifecycle and progress-tracking boundary. It is a specification document only. It does not implement FastAPI routes, service code, frontend code, Snakemake execution, databases, or storage backends.

The design answers the following questions:

1. What is the minimal run lifecycle model?
2. What are the run primitives?
3. Where should the boundary live?
4. What API contract should future FastAPI expose?
5. How should the frontend progress UI consume it later?
6. How should the agent interact with it later?
7. What is explicitly out of scope for the first implementation?

## 1. Goal and non-goals

### 1.1 Goal

Add a stable, workflow-neutral run lifecycle and progress-tracking boundary to the workflow platform. The boundary must support:

- Creating a run for a registered workflow.
- Tracking a finite, observable run status.
- Emitting ordered run events that describe status and stage changes.
- Appending log chunks for later display.
- Recording artifact references produced by a run.
- Cancelling a run before or during execution.
- Exposing lifecycle operations and progress inspection through API endpoints, while keeping later agent tools read-only.

### 1.2 Non-goals

This phase explicitly does **not** implement or specify:

| Out of scope | Why |
|---|---|
| Real Snakemake execution, process spawning, subprocess management, signal handling | Execution runtime is a future phase. PR99 defines only the boundary where a runtime will plug in. |
| HPC profiles, scheduler integration, cloud batch, SLURM/LSF/PBS | Compute environment abstraction is future work. |
| Authentication and authorization | No user model exists yet. |
| Persistence DB schema, migrations, ORM models | RunService may start as an in-memory stub. Storage is an implementation detail. |
| Artifact browser, artifact download endpoints, artifact content serving | Artifact references are opaque pointers only in this phase. |
| QC dashboard, QC metric aggregation, QC interpretation | QC data is not surfaced beyond possible artifact references. |
| Provenance graph, RO-Crate, lineage tracking | Provenance is a future platform phase. |
| Multi-user collaboration, workspace sharing, concurrency policies | Single-user platform semantics for now. |
| DAG visualization | DAG previews are adapter capability future work, not run lifecycle. |
| Server-Sent Events or WebSocket streaming | Polling is the first transport; streaming is a future optimization. |
| Human-in-the-loop approval gates | No state-changing agent actions are permitted. |
| Real-time filesystem log tailing | Log chunks are platform-owned event data, not live file streams. |
| Extending `WorkflowAdapter` with `preview_dag`, `plan_workspace`, or `build_command` | Those methods remain unsupported in the ENCODE adapter until execution planning is designed. |

## 2. Design approach

PR99 adopts a **platform-service boundary spec** approach:

- Define concrete run primitives in the platform layer.
- Define a thin `RunService` / `RunManager` service boundary that owns the lifecycle state machine, event ordering, and artifact references.
- Define HTTP endpoints that expose the service boundary without leaking execution internals.
- Keep `WorkflowAdapter` as a workflow capability boundary. Adapters expose capabilities such as `validation`; they do not store or manage runs.
- Reserve an execution-adapter boundary for a future `RunDriver` or `ExecutionAdapter` interface that plugs into `RunService`.

This approach was reviewed by Backend/API Architect, Workflow Execution Architect, Frontend UX Architect, Security/Safety Reviewer, Test/QA Architect, and Scope Guard/Code Reviewer perspectives. All reviewers preferred the platform-service boundary approach over a minimal API-only spec or an adapter-extended lifecycle spec.

## 3. Run lifecycle model

### 3.1 Status enum

```text
created
validating
planned
queued
running
succeeded
failed
cancelled
```

`created`, `validating`, `planned`, `queued`, and `running` are **active** states. `succeeded`, `failed`, and `cancelled` are **terminal** states.

### 3.2 Valid transitions

```text
                    +------------+
                    |  created   |
                    +------+-----+
                           |
                           v
                    +------------+
                    | validating |
                    +------+-----+
                           |
              +------------+------------+
              |                         |
              v                         v
       +------------+            +------------+
       |   failed   |            |   planned  |
       +------------+            +------+-----+
                                        |
                                        v
                                 +------------+
                                 |   queued   |
                                 +------+-----+
                                        |
                                        v
                                 +------------+
                                 |  running   |
                                 +------+-----+
                                        |
                    +-------------------+-------------------+
                    |                   |                   |
                    v                   v                   v
             +------------+      +------------+      +------------+
             | succeeded  |      |   failed   |      | cancelled  |
             +------------+      +------------+      +------------+
```

A run may transition from any active state to `cancelled`. A run may transition from `created` directly to `cancelled`. Terminal states are absorbing: no further transitions occur after reaching `succeeded`, `failed`, or `cancelled`.

The `validating` state is entered immediately after creation because the platform runs `ValidationService.validate()` before any execution planning. If validation fails, the run transitions to `failed`. If validation succeeds, the run transitions to `planned`.

The `planned` and `queued` states are preparation states for future execution. In the first implementation they may be simulated or stubbed, but they are retained in the model so that later execution integration does not require a new status vocabulary.

### 3.3 Cancellation semantics

Cancellation is a platform-level state transition. In PR99 there is no external process, scheduler job, or HPC task to kill. The cancel endpoint records intent, transitions the run to `cancelled`, and halts future event and log emission. The endpoint is idempotent: cancelling an already-terminal run returns the current terminal state without error.

A future implementation may extend cancellation to signal an execution driver, but the platform status transition remains the canonical record.

## 4. Run primitives

The following primitives are **non-executable reference sketches**. They illustrate the shapes that future platform models, Pydantic request/response models, and TypeScript types should follow. They are not committed as code files.

### 4.1 RunId

Opaque, unique run identifier. A UUID string is recommended.

### 4.2 RunStatus

```python
# Python-ish
from enum import Enum

class RunStatus(str, Enum):
    CREATED = "created"
    VALIDATING = "validating"
    PLANNED = "planned"
    QUEUED = "queued"
    RUNNING = "running"
    SUCCEEDED = "succeeded"
    FAILED = "failed"
    CANCELLED = "cancelled"

    @property
    def is_terminal(self) -> bool:
        return self in {
            RunStatus.SUCCEEDED,
            RunStatus.FAILED,
            RunStatus.CANCELLED,
        }
```

```typescript
// TypeScript-ish
type RunStatus =
  | 'created'
  | 'validating'
  | 'planned'
  | 'queued'
  | 'running'
  | 'succeeded'
  | 'failed'
  | 'cancelled';

function isTerminal(status: RunStatus): boolean {
  return ['succeeded', 'failed', 'cancelled'].includes(status);
}
```

### 4.3 RunRecord

```python
# Python-ish
from dataclasses import dataclass
from datetime import datetime
from typing import Any

@dataclass(frozen=True)
class RunRecord:
    run_id: str
    workflow_id: str
    inputs: dict[str, Any]  # WorkflowInputs shape
    status: RunStatus
    created_at: datetime
    updated_at: datetime
    started_at: datetime | None
    ended_at: datetime | None
    current_stage: str | None
    cancellation_reason: str | None
    error: Issue | None
    tags: dict[str, str]
```

```typescript
// TypeScript-ish
interface RunRecord {
  run_id: string;
  workflow_id: string;
  inputs: WorkflowInputs;
  status: RunStatus;
  created_at: string; // ISO 8601
  updated_at: string;
  started_at: string | null;
  ended_at: string | null;
  current_stage: string | null;
  cancellation_reason: string | null;
  error: Issue | null;
  tags: Record<string, string>;
}
```

`started_at` is set on transition into `running`. `ended_at` is set on transition into a terminal state. `error` is non-null only when `status` is `failed`. `cancellation_reason` is non-null only when `status` is `cancelled`.

### 4.4 RunEvent

Run events are immutable, ordered records of state changes.

```python
# Python-ish
from dataclasses import dataclass
from datetime import datetime
from typing import Any

@dataclass(frozen=True)
class RunEvent:
    event_id: str
    run_id: str
    sequence: int
    event_type: str  # e.g. "status_changed", "stage_started", "stage_completed", "issue_added", "cancelled"
    timestamp: datetime
    status: RunStatus | None
    stage: str | None
    message: str
    context: dict[str, Any]
    issue: Issue | None
```

```typescript
// TypeScript-ish
interface RunEvent {
  event_id: string;
  run_id: string;
  sequence: number;
  event_type:
    | 'status_changed'
    | 'stage_started'
    | 'stage_completed'
    | 'issue_added'
    | 'cancelled';
  timestamp: string;
  status: RunStatus | null;
  stage: string | null;
  message: string;
  context: Record<string, unknown>;
  issue: Issue | null;
}
```

Events within a run are strictly ordered by `sequence`. The platform guarantees monotonic `sequence` values and a strict total order. The `event_id` is an opaque cursor for pagination.

### 4.5 RunLogChunk

Log chunks are append-only, human-readable text fragments.

```python
# Python-ish
from dataclasses import dataclass
from datetime import datetime

@dataclass(frozen=True)
class RunLogChunk:
    chunk_id: str
    run_id: str
    stream_name: str  # e.g. "stdout", "stderr", "engine"
    sequence: int
    timestamp: datetime
    lines: list[str]
```

```typescript
// TypeScript-ish
interface RunLogChunk {
  chunk_id: string;
  run_id: string;
  stream_name: string;
  sequence: number;
  timestamp: string;
  lines: string[];
}
```

Like events, log chunks are ordered by `sequence` within a stream. The frontend should treat log chunks as display-only text; it must not parse log lines for UI state.

### 4.6 RunArtifactRef

Artifact references are opaque pointers to outputs produced by a run. They do not contain file contents.

```python
# Python-ish
from dataclasses import dataclass
from datetime import datetime
from typing import Any

@dataclass(frozen=True)
class RunArtifactRef:
    artifact_id: str
    run_id: str
    artifact_type: str  # e.g. "file", "directory", "qc_report", "log"
    name: str
    uri: str
    mime_type: str | None
    produced_at: datetime
    metadata: dict[str, Any]
```

```typescript
// TypeScript-ish
interface RunArtifactRef {
  artifact_id: string;
  run_id: string;
  artifact_type: 'file' | 'directory' | 'qc_report' | 'log';
  name: string;
  uri: string;
  mime_type: string | null;
  produced_at: string;
  metadata: Record<string, unknown>;
}
```

The `uri` is intentionally opaque. It may use a platform scheme such as `run://runs/{run_id}/artifacts/{name}` rather than exposing absolute filesystem paths. Artifact content serving and download endpoints are out of scope.

## 5. Boundary location

### 5.1 Platform primitives

Run primitives (`RunId`, `RunStatus`, `RunRecord`, `RunEvent`, `RunLogChunk`, `RunArtifactRef`) belong in the platform layer, alongside `Result`, `Issue`, `WorkflowMetadata`, `WorkflowInputs`, and the `WorkflowAdapter` protocol. They are workflow-neutral and JSON-serializable.

### 5.2 Service boundary

A `RunService` or `RunManager` service owns:

- Creating runs.
- Enforcing valid status transitions.
- Emitting ordered run events.
- Appending log chunks.
- Recording artifact references.
- Handling cancellation.
- Coordinating with `ValidationService` for the `validating` state.

The service is an internal platform boundary, not a deployed microservice. Future implementation may place it in `src/encode_pipeline/services/run.py`.

### 5.3 API routes

FastAPI routes remain thin translators:

- Parse HTTP request bodies into platform inputs.
- Call `RunService` methods.
- Serialize platform primitives into the JSON envelopes documented below.
- Return appropriate HTTP status codes.

Routes do not contain lifecycle logic, execution logic, or workflow-specific behavior.

### 5.4 Execution adapter boundary

A future execution adapter or `RunDriver` interface will plug into `RunService`. It will be responsible for taking a `RunRecord` and its planned workspace/command and turning it into actual execution progress. PR99 does not define this interface.

The existing `WorkflowAdapter` protocol remains unchanged. Its `preview_dag`, `plan_workspace`, and `build_command` methods remain unsupported in the ENCODE adapter until execution planning is designed.

## 6. API contract

All endpoints use endpoint-specific envelopes that include `ok` and `issues`. Issue objects follow the existing `IssueResponse` shape from the validation MVP contract.

### 6.1 `POST /api/v1/workflows/{workflow_id}/runs`

Creates a new run for the specified workflow.

Request body:

```json
{
  "inputs": {
    "config": { "use_control": false },
    "samples": "samples.tsv",
    "options": { "strict_inputs": false }
  },
  "idempotency_key": "550e8400-e29b-41d4-a716-446655440000"
}
```

`inputs` follows the `WorkflowInputs` shape. `idempotency_key` is optional.

Success response (HTTP 201):

```json
{
  "ok": true,
  "run_id": "6ba7b810-9dad-11d1-80b4-00c04fd430c8",
  "workflow_id": "encode-style-chipseq-cuttag-atac-mnase",
  "status": "created",
  "created_at": "2026-07-03T12:00:00Z",
  "updated_at": "2026-07-03T12:00:00Z",
  "started_at": null,
  "ended_at": null,
  "current_stage": null,
  "cancellation_reason": null,
  "error": null,
  "tags": {},
  "issues": []
}
```

Idempotent replay response (HTTP 200): same body with `ok: true` when the same `idempotency_key` is reused.

Error responses:

- `404 Not Found` if `workflow_id` is not registered. Issue code `WORKFLOW_NOT_FOUND`.
- `409 Conflict` if the workflow does not support the `run` capability. Issue code `WORKFLOW_CAPABILITY_UNSUPPORTED`.
- `400 Bad Request` if the request body is malformed. Issue code `API_REQUEST_INVALID`.

### 6.2 `GET /api/v1/runs/{run_id}`

Returns the current run record.

Success response (HTTP 200):

```json
{
  "ok": true,
  "run": {
    "run_id": "6ba7b810-9dad-11d1-80b4-00c04fd430c8",
    "workflow_id": "encode-style-chipseq-cuttag-atac-mnase",
    "inputs": {
      "config": { "use_control": false },
      "samples": "samples.tsv",
      "options": { "strict_inputs": false }
    },
    "status": "validating",
    "created_at": "2026-07-03T12:00:00Z",
    "updated_at": "2026-07-03T12:00:05Z",
    "started_at": null,
    "ended_at": null,
    "current_stage": "validation",
    "cancellation_reason": null,
    "error": null,
    "tags": {}
  },
  "issues": []
}
```

Error response:

- `404 Not Found` if `run_id` does not exist. Issue code `RUN_NOT_FOUND`.

### 6.3 `GET /api/v1/runs/{run_id}/events`

Returns ordered run events. Supports cursor-based pagination.

Query parameters:

- `after` (optional): opaque event cursor. Returns events after this cursor.
- `limit` (optional): integer, default 50, maximum 1000.

Success response (HTTP 200):

```json
{
  "ok": true,
  "run_id": "6ba7b810-9dad-11d1-80b4-00c04fd430c8",
  "events": [
    {
      "event_id": "evt-1",
      "run_id": "6ba7b810-9dad-11d1-80b4-00c04fd430c8",
      "sequence": 1,
      "event_type": "status_changed",
      "timestamp": "2026-07-03T12:00:00Z",
      "status": "created",
      "stage": null,
      "message": "Run created.",
      "context": { "previous_status": null, "new_status": "created" },
      "issue": null
    },
    {
      "event_id": "evt-2",
      "run_id": "6ba7b810-9dad-11d1-80b4-00c04fd430c8",
      "sequence": 2,
      "event_type": "status_changed",
      "timestamp": "2026-07-03T12:00:01Z",
      "status": "validating",
      "stage": "validation",
      "message": "Validation started.",
      "context": { "previous_status": "created", "new_status": "validating" },
      "issue": null
    }
  ],
  "next_cursor": "evt-2",
  "issues": []
}
```

`next_cursor` is null when no further events exist. Clients should stop polling when the run is terminal and `next_cursor` is null.

### 6.4 `GET /api/v1/runs/{run_id}/logs`

Returns log chunks for one run.

Query parameters:

- `stream` (optional): stream name, default `"stdout"`.
- `after` (optional): opaque chunk cursor.
- `limit` (optional): integer, default 50, maximum 1000.

Success response (HTTP 200):

```json
{
  "ok": true,
  "run_id": "6ba7b810-9dad-11d1-80b4-00c04fd430c8",
  "stream": "stdout",
  "chunks": [
    {
      "chunk_id": "log-1",
      "run_id": "6ba7b810-9dad-11d1-80b4-00c04fd430c8",
      "stream_name": "stdout",
      "sequence": 1,
      "timestamp": "2026-07-03T12:00:10Z",
      "lines": ["Initializing run...", "Loaded sample sheet."]
    }
  ],
  "next_cursor": "log-1",
  "issues": []
}
```

### 6.5 `POST /api/v1/runs/{run_id}/cancel`

Requests cancellation of a run.

Request body:

```json
{
  "reason": "User requested cancellation."
}
```

Success response (HTTP 200):

```json
{
  "ok": true,
  "run_id": "6ba7b810-9dad-11d1-80b4-00c04fd430c8",
  "workflow_id": "encode-style-chipseq-cuttag-atac-mnase",
  "status": "cancelled",
  "created_at": "2026-07-03T12:00:00Z",
  "updated_at": "2026-07-03T12:00:15Z",
  "started_at": null,
  "ended_at": "2026-07-03T12:00:15Z",
  "current_stage": null,
  "cancellation_reason": "User requested cancellation.",
  "error": null,
  "tags": {},
  "issues": []
}
```

Cancellation is idempotent. If the run is already terminal, the current terminal state is returned with `ok: true`.

Error response:

- `404 Not Found` if `run_id` does not exist. Issue code `RUN_NOT_FOUND`.

## 7. Frontend progress UI consumption

The frontend should implement progress tracking with **polling first**. Server-Sent Events and WebSocket streaming are deferred until latency or scale requirements justify the added connection-management complexity.

### 7.1 Polling strategy

- Poll `GET /api/v1/runs/{run_id}/events` with the `after` cursor.
- Start with a 2-second interval while the run is active.
- Back off to 5–10 seconds in `queued`.
- Stop or slow to 30 seconds in terminal states.
- Handle empty responses by keeping the current cursor and continuing to poll.
- On receiving a terminal event with `next_cursor: null`, perform one final poll, then stop.

### 7.2 Progress UI components

- **Status badge**: displays `RunStatus` with color and label.
- **Stage timeline**: driven by `stage_started` / `stage_completed` events.
- **Event feed**: structured list of `RunEvent` messages.
- **Log panel**: append-only display of `RunLogChunk.lines`, grouped by stream.
- **Cancel button**: enabled for active states; disabled for terminal states.

### 7.3 Events vs logs

- **Events** are structured, machine-readable state changes. The UI derives status, stage, and timeline from events.
- **Logs** are human-readable text streams. The UI displays them but does not parse them for state.

### 7.4 Error handling

- Network failures show a retry banner without losing already-fetched events.
- A failed run displays the `error` Issue from the run record prominently.
- Unknown run IDs show a 404 message with a link back to the workflow catalog.

## 8. Agent interaction

The agent layer remains **read-only** in this phase. It may inspect runs, events, logs, and artifact references, but it may not create, cancel, delete, or mutate runs.

### 8.1 Allowed future tools

Examples of safe read-only tool names for a future implementation PR:

- `list_runs`
- `get_run_status`
- `list_run_events`
- `explain_run_status`
- `summarize_run_issues`

### 8.2 Forbidden tool prefixes

Any new tool must pass the existing `ReadOnlyToolRegistry` deny list. Names matching `run_*`, `submit_*`, `execute_*`, `cancel_*`, `delete_*`, `update_*`, `apply_*`, `write_*`, or `modify_*` are rejected.

### 8.3 Agent context

When explaining a run, the agent receives:

- The run record (status, stage, timestamps).
- Recent run events.
- Structured issues from validation or failure.
- Schema hints for the workflow.

It does not receive:

- Absolute filesystem paths.
- Environment variables or secrets.
- Raw command lines.
- Unredacted log lines.

### 8.4 Safety framing

The agent sidebar should continue to display a read-only label such as **Run Assistant — Read Only**. The agent may answer questions like "Why is my run stuck in validating?" or "What does the `planned` status mean?" It must not offer to cancel, retry, or submit runs.

## 9. State-machine invariants and event guarantees

The spec requires the following invariants for any future implementation:

1. **Terminal absorption**: once a run reaches `succeeded`, `failed`, or `cancelled`, no further status transitions occur.
2. **Monotonic sequences**: `RunEvent.sequence` and `RunLogChunk.sequence` are strictly increasing within a run and stream.
3. **Status/event consistency**: the `RunRecord.status` always matches the latest `RunEvent.status` when the latest event carries one.
4. **Timestamp ordering**: `created_at` <= `updated_at` <= any event timestamp. `started_at` is non-null only after entering `running`. `ended_at` is non-null only after entering a terminal state.
5. **Idempotent cancel**: repeated cancel requests for the same run return success with the current terminal state.
6. **One event per transition**: each status transition emits exactly one `status_changed` event.
7. **Event cursor stability**: once an event has been assigned an `event_id`, it is never modified or reordered.

## 10. Testability notes

A future implementation should be testable without real execution:

- `RunService` should accept a test double for the execution driver.
- Tests can inject synthetic adapter outcomes to drive status transitions.
- Contract tests should verify that repeated event polls reconstruct the same ordered sequence.
- Transition-matrix tests should assert illegal transitions are rejected.
- Cancellation tests should assert idempotency and terminal absorption.
- Redaction tests should verify that log/event content does not expose absolute paths or secrets before API or agent exposure.

## 11. Out-of-scope recap

This PR99 spec does not include and must not be interpreted to require:

- Real Snakemake execution, process spawning, subprocess management, or signal handling.
- HPC profiles, scheduler integration, cloud batch, or container orchestration.
- Authentication, authorization, multi-user isolation, or workspace sharing.
- Persistence DB schema, migrations, or ORM models.
- Artifact browser, artifact content serving, or download endpoints.
- QC dashboard, QC metric interpretation, or QC visualization.
- Provenance graph, RO-Crate, or lineage tracking.
- DAG visualization or interactive graph rendering.
- Server-Sent Events, WebSocket streaming, or real-time file tailing.
- Human-in-the-loop approval gates or agent-driven mutations.
- Extension of `WorkflowAdapter` with `preview_dag`, `plan_workspace`, or `build_command`.

## 12. File touch guardrail

PR99 changes only this markdown file under `docs/development/`. It does not change:

- Python source files.
- FastAPI routes or models.
- Frontend TypeScript/React files.
- `pyproject.toml` or lockfiles.
- CI configuration.
- Snakemake files or workflow rules.
- `docs/superpowers/`, `research/`, `encode-pipeline-architecture-report.md`, or `agent-layer-architecture-report.md`.

## 13. Verification checklist

Before PR99 is submitted:

- `git diff --name-only origin/main..HEAD` shows only `docs/development/2026-07-03-run-lifecycle-progress-boundary-design.md`.
- `git diff --check origin/main..HEAD` is clean.
- Placeholder scan over this file returns no matches for common unfinished-marker patterns.
- Commit messages contain no attribution, tool, or vendor trailers.
- The spec contains no executable code files, no Pydantic/TypeScript source files, and no implementation plan for PR100.
