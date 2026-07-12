# PR131 Read-only QC Metrics API Design

## Scope

PR131 exposes the workflow-neutral `RunQcMetric` rows already indexed by PR130.
The request path reads only the run repository through `RunService`; it never
opens a workspace path, parses a source summary, invokes an adapter, changes a
run, or retries indexing. SQLite remains the sole source of QC metadata.

The API adds one operation:

- `GET /api/v1/runs/{run_id}/qc-metrics`, with the explicit operation ID
  `listRunQcMetrics`.

QC charts, thresholds, reindexing, mutation, deletion, downloads, source-file
content, and frontend product work remain out of scope.

## Alternatives considered

1. Load every metric and paginate in the route. This preserves the current
   repository signature but permits an unbounded database read and process
   allocation.
2. Add a separate QC query repository or read-model service. This could support
   future filters, but one bounded list operation does not justify a second
   persistence abstraction.
3. Extend the existing repository and `RunService` with run-scoped metric-ID
   keyset pagination. This is selected because it matches the artifact API,
   keeps SQL pagination bounded, and preserves the established dependency
   direction.

## Query and pagination contract

The endpoint accepts `after` and `limit`. `limit` defaults to 50 and is bounded
to 1 through 100. `after`, when present, must use the durable
`qcmetric-<64 lowercase hexadecimal characters>` grammar. The repository orders
rows by `metric_id` ascending and the route requests `limit + 1` rows. When an
extra row exists, `next_cursor` is the last returned metric ID.

Before applying `metric_id > after`, the repository resolves the cursor with
both `run_id` and `metric_id` and validates the complete cursor row. A missing
or cross-run cursor returns HTTP 400 with `RUN_QC_METRIC_CURSOR_NOT_FOUND`. A
damaged cursor row fails closed like any damaged result row and returns the
path-free `RUN_QC_METRIC_DATA_INVALID` response. Cursor values are never echoed.

An unknown run returns HTTP 404 with the existing `RUN_NOT_FOUND` issue. An
existing run with no metrics returns a successful empty page. Atomic QC
replacement can occur between page requests; ordering is stable within one
generation, while a replaced cursor may become invalid. Clients must restart
from the first page after that response. The API does not claim snapshot
pagination across generations.

## Public response projection

`QcMetricResponse` exposes only:

- `metric_id`, `metric_key`, and `display_name`;
- `value` as canonical non-exponent decimal text;
- `unit` and `scope`;
- nullable `sample_id`, `experiment_id`, `assay`, and `qc_flag`;
- `source_artifact_id` and UTC `produced_at`.

`value` is a Pydantic `str`, not `Decimal`, so OpenAPI and Orval generate a
lossless TypeScript string. The projection uses PR130's canonical decimal
serializer and validates the exact persisted grammar. It also requires the
domain record's `run_id` to equal the requested run, recomputes the deterministic
metric ID from its semantic coordinates, validates scope/identifier consistency,
restricts unit and flag vocabularies, and requires an aware timestamp.

The envelope contains `ok`, `run_id`, `qc_metrics`, `next_cursor`, and `issues`.
No workspace path, database URL, source TSV bytes, command, environment value,
ORM row, arbitrary metadata, exception string, or `technical_message` is part
of the model. Unsafe display text, malformed identifiers, invalid decimal
values, or mismatched IDs fail the whole page with HTTP 500 and
`RUN_QC_METRIC_DATA_INVALID`; the response contains no offending value.

## Service and persistence boundary

`RunRepository.list_qc_metrics(run_id, after, limit)` has matching InMemory and
SQLAlchemy semantics. Both implementations validate the cursor and every
returned domain value with PR130's shared durable validator before returning.
SQLAlchemy uses the existing `(run_id, metric_id)` uniqueness index and applies
ordering and limit in SQL. No migration is required.

`RunService.list_qc_metrics` validates positive internal limits and delegates.
The route is a synchronous FastAPI endpoint, so Starlette runs its bounded
`RunService.get_run` plus paginated SQLite read in the maintained worker
threadpool instead of blocking the event loop. It performs no Redis, worker,
workspace, or file I/O and never uses a SQLAlchemy session directly.
Database/infrastructure
exceptions continue through the application's generic redacted 500 boundary;
known persisted-row and projection validation failures use the QC-specific
invalid-data envelope.

FastAPI request validation for `listRunQcMetrics` has a dedicated handler branch
that returns `RunQcMetricsResponse` with `API_REQUEST_INVALID`, an allowlisted
message/path/hint, and no technical detail. This prevents the general validation
envelope from drifting from the operation's declared contract.

## OpenAPI and generated client

The new router is tagged `qc-metrics`, registered under `/api/v1`, and declares
stable 400, 404, and 500 response models. The OpenAPI operation gate asserts the
explicit unique ID and the absence of an unreachable default 422 response.
`frontend/openapi.json` and the Orval output are regenerated mechanically; no
TypeScript DTO or client is written by hand and no frontend product code or
package lock changes.

## Verification and residual risks

Tests cover repository parity, stable order, bounded pages, cursor isolation,
cursor-row corruption, empty results, unknown runs, exact large/fractional
decimals, nullable fields, SQLite reopen, strict projection, disclosure
fail-closed behavior, request bounds, operation IDs, OpenAPI export, and Orval
drift. Full Python and frontend regression gates remain required.

Residual risks are limited to the existing local deployment boundary: keyset
pagination is not a multi-request snapshot, and a database outage uses the
generic redacted infrastructure error rather than a QC-specific availability
code. PR131 deliberately does not add authentication; run isolation is enforced
by query predicates, while authorization remains a later platform concern.
