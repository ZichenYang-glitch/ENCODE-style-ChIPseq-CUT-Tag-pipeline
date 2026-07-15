# HelixWeave Architecture Overview

HelixWeave is a workflow-neutral, local/small-team omics platform. This page is
the maintained entry point for its durable architecture. It intentionally
describes ownership and contracts rather than implementation history.

## System boundary

The platform turns reviewed inputs into durable evidence:

```text
author inputs
  -> adapter validation
  -> immutable validated snapshot
  -> SQLite run lifecycle
  -> Redis/RQ execution handoff
  -> workflow execution
  -> adapter artifact/QC extraction
  -> safe API projection
  -> browser evidence review
```

The supported deployment is a workstation or a small trusted team sharing a
filesystem and runtime. SQLite and local workspaces are product choices for
that scope, not placeholders for an implied hosted service.

## Ownership

Workflow-neutral code owns concepts that must work for every adapter:

- `platform/` defines domain contracts such as results, issues, workflows,
  runs, events, artifacts, QC, and lifecycle states.
- `services/` orchestrates validation, planning, submission, execution,
  cancellation, indexing, and read-only Agent context.
- `api/` translates HTTP requests and responses over service contracts.
- `persistence/` implements SQLAlchemy repositories and Alembic migrations.
- `workers/` implements Redis/RQ process and queue mechanics.
- `frontend/` consumes the generated API client and presents workflow/run UX.

Adapter code owns workflow-specific schemas, validation, workspace planning,
commands, and artifact/QC extraction. The Snakemake files under `workflow/`
own scientific rules, targets, config keys, environments, and output behavior.
Generic platform layers must not depend on those private details.

See [Adapter conformance](../development/adapter-conformance.md) for the
executable adapter boundary and [assay policy](../assay-policy.md) for bundled
workflow behavior.

## Canonical state and execution

SQLite is the canonical store for run identity, lifecycle, events, logs, and
indexed result metadata. Lifecycle transitions use repository contracts and
compare-and-set ownership so API and worker processes cannot both claim truth.

Redis/RQ is an execution handoff. Queue state may help operate workers, but it
does not replace SQLite lifecycle state. Jobs carry durable run identity and
reconstruct dependencies inside the worker process.

Long-running Snakemake processes execute outside HTTP handlers. Cancellation
is complete only after process-group termination is acknowledged and the
canonical lifecycle transition succeeds. Recovery must never claim that an
orphaned process is still running.

Lifecycle, repository, worker, cancellation, and recovery tests are the
executable detail behind these invariants.

## Validated snapshot contract

Run submission is based on an immutable server-validated snapshot, not on
untrusted browser validation or a mutable request body. The snapshot binds the
workflow identity, adapter schema/validation result, normalized inputs, and
workspace plan used for execution.

The worker executes that durable identity. Revalidation or reconstruction must
not silently change what the user submitted. Public platform code may use
workflow-neutral snapshot fields but must not expose adapter-private validation
payloads.

## HTTP and generated-client contract

FastAPI routes remain thin over service objects. Pydantic models define public
request, response, and error envelopes. The exported OpenAPI document is the
source for the Orval-generated frontend client and mocks.

Generated files under `frontend/src/api/generated/` are never edited by hand.
Backend contract changes require regeneration and a zero-drift check. Stable
routes remain under `/api/v1`; transport errors and structured workflow issues
are distinct concepts.

## Artifacts, QC, and downloads

Adapters map successful workflow outputs into platform-owned artifact and QC
records. The platform indexes metadata in SQLite while files remain in the
run workspace. Artifact URIs are stable identifiers, not arbitrary filesystem
paths.

Downloads resolve an indexed artifact through fail-closed workspace policy.
The API never accepts an arbitrary absolute path and does not reveal workspace
locations, symlink targets, environment values, or raw exception text.

The bundled workflow's scientific vocabulary is documented in the
[output contract](../output-contract.md), [QC guide](../qc-interpretation.md),
and [artifact inventory](artifact-inventory.yaml). Runtime artifact generation
remains owned by Snakemake; see the
[artifact adoption decision](artifact-adoption-decision.md).

## Read-only Agent boundary

The Agent may interpret structured schemas, issues, run state, logs, artifacts,
and QC through platform services. It may explain or propose next steps, but it
must not submit, start, cancel, mutate, or delete runs.

Agent explanations are not provenance. Context construction applies the same
redaction and path rules as public APIs, and unsupported operations fail closed.
Service import-boundary, redaction, route, and frontend tests enforce this
contract.

## Compatibility anchors

Maintenance work preserves the `encode_pipeline` import namespace, the
`encode-pipeline` distribution, existing `encode-*` CLIs, `/api/v1` routes,
workflow/adapter identity, Alembic history, SQLite fields, Redis/RQ job
identity, artifact URIs, environment variables, and Snakemake behavior unless
a separately reviewed migration explicitly changes one of those contracts.

## Operational references

- [Local platform runtime](../development/local-platform-runtime.md)
- [Configuration reference](../configuration.md)
- [Sample sheet reference](../sample-sheet.md)
- [Reproducibility policy](../reproducibility-policy.md)
- [Current product roadmap](../development/workflow-platform-agent-roadmap.md)
