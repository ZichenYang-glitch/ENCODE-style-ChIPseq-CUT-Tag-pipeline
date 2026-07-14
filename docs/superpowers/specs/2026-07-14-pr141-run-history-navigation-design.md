# PR141 Run History and Navigation Design

**Status:** locked for implementation after bounded design review
**Scope:** read-only run history, global navigation, and run-result discoverability

## Outcome

PR141 adds one canonical, persistent run-history view without changing run
lifecycle or result semantics. `GET /api/v1/runs` reads SQLite through
`RunService` and `RunRepository`; `/runs` renders that contract through the
generated Orval operation. The global primary navigation becomes exactly
`Runs`, `Workflows`, and (only when a workflow is in the URL) `New analysis`.
QC and Artifacts remain run-local tabs.

## API contract

`GET /api/v1/runs` has explicit `operation_id=listRuns` and query parameters:

- `limit`: integer, default 50, inclusive range 1–100.
- `after`: optional opaque cursor, bounded to 1,024 ASCII characters.
- `workflow_id`: optional trimmed, non-control string, maximum 255 characters.
- `status`: optional exact `RunStatus` value.

The success envelope is:

```json
{
  "ok": true,
  "runs": [
    {
      "run_id": "...",
      "workflow_id": "...",
      "status": "succeeded",
      "created_at": "2026-07-14T01:02:03Z",
      "updated_at": "2026-07-14T01:03:03Z",
      "started_at": "2026-07-14T01:02:10Z",
      "ended_at": "2026-07-14T01:03:03Z",
      "current_stage": "qc_summary_indexing"
    }
  ],
  "next_cursor": null,
  "issues": []
}
```

The repository returns a dedicated workflow-neutral `RunSummary` read model.
The SQL query selects only its nine public columns; it does not load inputs,
tags, cancellation reasons, persisted Issues, technical messages, commands,
workspace/database paths, or environment values. The existing full-record
`list_runs()` remains unchanged for internal recovery code. Summary validation
rejects unsafe identifiers, control characters, path-like text, invalid
lifecycle values, naive timestamps, inconsistent timestamp ordering, and
status/evidence drift. Terminal status is equivalent to a non-null `ended_at`;
`RUNNING` and `SUCCEEDED` require `started_at`; pre-execution non-terminal
statuses forbid it. A corrupt selected row fails the whole request closed with the stable,
disclosure-safe `RUN_HISTORY_DATA_INVALID` issue.

Malformed or filter-mismatched cursors return `RUN_HISTORY_CURSOR_INVALID`.
A structurally valid cursor whose boundary run has been deleted returns
`RUN_HISTORY_CURSOR_NOT_FOUND`. FastAPI query-shape failures use a dedicated
`RunHistoryResponse` `API_REQUEST_INVALID` envelope. None of these issues echo
the cursor or persisted values.

## Cursor and keyset semantics

Ordering is always `created_at DESC, run_id DESC`. The cursor is a versioned,
canonical base64url JSON token with exact keys for:

- boundary `created_at` in normalized UTC form;
- boundary `run_id`;
- the normalized `workflow_id` filter or null;
- the status filter or null.

The token is opaque protocol state, not a credential and not a confidentiality
mechanism. Decoding is strict: bounded/canonical base64url, exact JSON fields
and types, supported version, normalized timestamp, safe IDs, and an exact
match with the current filters. The repository additionally requires the
boundary run still to exist with the same `created_at` and to satisfy the
current workflow/status filters. That makes deleted boundaries, lifecycle
changes out of a status-filtered set, and cross-filter reuse deterministic.

The service asks the repository for `limit + 1`, returns at most `limit`, and
builds `next_cursor` from the last returned row only when another row exists.
An inserted row with a strictly greater `(created_at, run_id)` sort key is ahead
of the boundary and is skipped by the current traversal; a strictly smaller
key, including a same-time lexically smaller ID or a clock-regressed timestamp,
may appear on a later page. Insertion ahead of the boundary does not repeat
previously returned rows. Rows after a deleted non-boundary item may naturally
disappear; history is a live read model, not a snapshot.

`RunRepository.list_runs()` remains unchanged for existing internal
recovery/test code. A separate `list_run_summaries` contract applies filters,
keyset predicates, descending order, and `LIMIT` in SQL; it never loads all
runs before pagination. Its boundary lookup selects only summary coordinates.
InMemory implements the same validation and ordering.

Every selected `limit + 1` row, including the sentinel, is validated before
the service truncates the public page. Boundary rows, selected rows, and the
row used to construct `next_cursor` all pass the same `RunSummary` validator.
InMemory additionally requires each mapping key to equal `record.run_id`.
Identity or safety failures uniformly become `RUN_HISTORY_DATA_INVALID`, so a
bad sentinel cannot be deferred to a later request.

A small Alembic migration adds four explicit ordering indexes: global
`(created_at, run_id)`, workflow-only `(workflow_id, created_at, run_id)`,
status-only `(status, created_at, run_id)`, and combined
`(workflow_id, status, created_at, run_id)`. SQLite query-plan tests require the
matching index and reject a temporary order-by B-tree for each filter shape.
The migration does not change stored lifecycle data.

The route is a synchronous FastAPI endpoint, so SQLAlchemy work runs in
FastAPI's threadpool rather than on the event loop.
Both global exception branches recognize `listRuns`: query validation returns
a `RunHistoryResponse` with `API_REQUEST_INVALID`, while an otherwise uncaught
driver/runtime failure returns the same envelope shape with a generic
`INTERNAL_SERVER_ERROR`. Neither branch includes cursor, filter, SQL, or
exception text.

## Frontend information architecture

`AppShell` owns the only primary navigation:

1. `Runs` — current for `/runs` and `/runs/:runId`.
2. `Workflows` — current for the catalog and exact workflow detail.
3. `New analysis` — derived from a URL workflow ID and current only on that
   workflow's authoring route.

There is no global QC/Artifacts navigation and no invented overview or data
page. The brand continues to return to Workflows. On a run detail URL, Runs is
current and returns to history; Activity, QC, and Artifacts retain their
existing query-string deep links.

`/runs` calls only the generated `listRuns` operation inside TanStack Query.
The query key contains normalized URL filters. `useInfiniteQuery` loads 50 rows
per page and accepts only a non-repeated next cursor. A malformed/replaced
cursor shows a recovery action that resets to the first page. Manual refresh
refetches canonical pages while already rendered data remains visible; a
transport failure with cached rows is an error banner, never a false empty
state.

The workflow filter options come from the existing workflow catalog query, but
an otherwise safe current URL value remains visible even if catalog loading
fails. A pure URL parser requires at most one `workflow_id` and one `status`.
Workflow values must be trimmed, 1–255 characters, and control-free; status
must be an exact backend enum. Duplicate, empty, oversized, control-containing,
or unknown-status parameters produce a controlled “filters could not be used”
state with Reset; they are not sent to the API and their raw content is not
rendered. Browser back/forward passes through this same parser.

Status options are the exact backend enum values. Changing either filter
updates the URL and starts at page one. Query data is retained only while
refetching the exact same filter key. A different filter key may use only its
own existing cache; it never shows another filter's rows under a new label.

Cursor values remain opaque to the UI. The client checks only their public
ASCII/base64url shape, length, and non-repetition. A response containing an
unsafe or repeated next cursor, and API errors
`RUN_HISTORY_CURSOR_INVALID`/`RUN_HISTORY_CURSOR_NOT_FOUND`, enter a dedicated
protocol-recovery state: already loaded rows remain visible and “Reload from
first page” clears that exact query. Ordinary transport failures instead offer
Retry and retain same-key cached rows. Neither class becomes an empty state.

Desktop uses a compact semantic table. Mobile uses the same data in separated
vertical records. Long IDs wrap with full accessible title text. Every row has
an Activity/run-detail link. Succeeded rows additionally link to
`?view=qc` and `?view=artifacts` without implying either collection is nonempty.
Loading, unfiltered/filtered empty, initial error, cached-data error, loading
more, and refresh states are explicit and use status/alert semantics.

## Testing and gate

- Domain/service/repository parity: stable tie ordering, filters, multi-page
  traversal, same-time greater/smaller-ID insertion, deleted boundary,
  cross-filter cursor, invalid token, invalid limits, corrupt boundary/selected/
  sentinel rows, SQL `LIMIT`, and index query plans.
- API: exact operation ID, threadpool execution, safe summary projection,
  success/empty/filter pages, every stable error envelope, and no sensitive
  fields in serialized responses/OpenAPI.
- Frontend: direct `/runs`, nav active rules, URL filters/history, responsive
  data/empty/loading/error/cached-error/cursor recovery, pagination, refresh,
  and succeeded QC/Artifact links.
- Real Playwright extends the existing deterministic run gate: author and run a
  real analysis, find its canonical ID in `/runs`, enter Activity, and verify
  QC/Artifacts deep links at 1440×900 and 1024×768. A separate, independently
  created successful run exercises the same history-to-results path in the
  actual mobile Playwright project at 390×844 and 360×800; viewport resizing in
  the desktop project is supplementary evidence only. Every viewport must have
  no document-level horizontal overflow.

## Explicit non-goals and residual limits

This PR does not add QC sources, global QC, run deletion, arbitrary search,
authentication, tenancy, draft persistence, another adapter, or workflow
changes. Pagination reflects live SQLite state rather than a historical
snapshot. Cursors are opaque and strictly/canonically validated as protocol
data, but are not signed and are not authorization tokens. A client can
reconstruct a cursor for an existing boundary; that only chooses a read
position. Future authorization must be applied independently to both
cursor-boundary validation and row queries.
