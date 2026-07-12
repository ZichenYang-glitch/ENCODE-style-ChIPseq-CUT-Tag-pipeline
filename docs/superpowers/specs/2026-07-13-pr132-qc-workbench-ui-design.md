# PR132 Read-only QC Workbench Design

## Scope and intent

PR132 adds a quiet, URL-addressable QC view to the existing run workbench. It
renders only the persisted, workflow-neutral metrics exposed by PR131 through
the generated `listRunQcMetrics` operation. It does not change the QC API,
indexer, adapter, scientific workflow, thresholds, or lifecycle state.

The existing run progress query remains the sole owner of canonical run,
event, and log polling. QC pages are independent read-only server state and do
not create a second run poller.

## Alternatives considered

1. Add a standalone QC route which reloads the run and events. This separates
   files but creates competing run polling and can make status/actions disagree.
2. Fetch every QC page inside `RunProgressPanel`. This keeps one component but
   couples lifecycle mutations, logs, pagination, and responsive presentation.
3. Keep `RunProgressPanel` as the workbench controller and pass its canonical
   status/outcome into a focused `QcWorkbench`. The workbench calls the
   generated QC operation through TanStack Query. This is selected because it
   preserves one lifecycle truth with a bounded component boundary.

## URL and navigation contract

- `/runs/:runId` and `?view=activity` show Activity.
- `?view=artifacts` shows Artifacts and preserves the existing optional
  `artifact=<artifact_id>` detail selection.
- `?view=qc` shows QC.
- Selecting QC or Activity removes a stale `artifact` parameter. Selecting a
  metric's source artifact sets both `view=artifacts` and the persisted opaque
  artifact ID.
- An `artifact` parameter continues to imply the Artifacts view. Unknown view
  values fall back to Activity.
- Search-parameter changes push history, so refresh, Back, Forward, and copied
  deep links restore the same view.

The three tabs use a single `tablist`, stable tab/panel IDs, roving tab index,
and Left/Right/Home/End keyboard navigation. Run ID, status, timestamps,
Start/Cancel/Refresh, cancellation acknowledgement, and bounded-polling notice
remain above all three views.

## QC outcome and honest state map

The latest event by sequence among `qc_metrics_indexed`,
`qc_metrics_indexing_failed`, and `qc_metrics_invalidated` defines the visible
QC outcome. An indexed event is valid only when `metric_count` is a
non-negative safe integer. A failure may display only an allowlisted stable
reason code. A missing outcome in a truncated event snapshot is unconfirmed,
not pending or empty.

```text
canonical run snapshot
  ├─ status != succeeded ───────────── QC available after successful completion
  └─ status == succeeded
       ├─ no QC outcome ────────────── indexing pending + bounded canonical poll
       ├─ latest = invalidated ─────── indexing pending + bounded canonical poll
       ├─ latest = failed ──────────── redacted failure + Refresh status
       ├─ malformed/truncated ──────── unconfirmed + Refresh status
       └─ latest = indexed
            ├─ first page loading ──── fixed-height skeleton rows
            ├─ initial API failure ─── redacted failure + Retry
            ├─ API failure + cache ─── preserve rows + Retry status
            ├─ count=0 + empty page ─ confirmed empty
            ├─ count>0 + empty page ─ inconsistent/not-ready + Retry
            ├─ cursor changed/missing ─ preserve rows + reload from first page
            └─ rows ───────────────── compact table or mobile entries
```

The run query polls at its existing interval only when lifecycle status is
active, or when the selected Artifacts/QC view is waiting for its corresponding
post-success outcome. The existing 15-minute pause and manual Refresh govern
that same query. The QC list itself is not periodically polled.

`qc_flag` is rendered exactly as the persisted backend fact (`pass`, `warning`,
`fail`, or not reported). The UI does not calculate thresholds, combine flags,
or infer scientific quality.

## Pagination and query safety

`QcWorkbench` uses TanStack Query v5 `useInfiniteQuery` and calls the Orval
generated `listRunQcMetrics` function directly. The query key is scoped by
`runId`, the initial cursor is absent, page size is 50, and more data is loaded
only through an explicit button.

Successful envelopes must have `ok=true`, the requested `run_id`, a metric
array, valid `qcmetric-<64 lowercase hex>` IDs, and safe scalar fields before
they render. Exception strings and issue technical details are never shown.
Pages are flattened in response order and duplicate metric IDs are ignored.

A next cursor is accepted only when it has the durable metric-ID grammar and
has not appeared as the current or an earlier page parameter. Invalid, empty,
or repeated cursors stop pagination. The UI explicitly says pagination could
not continue safely and offers a first-page reload. The same recovery is used
when the API returns `RUN_QC_METRIC_CURSOR_NOT_FOUND`, which is expected if an
atomic QC replacement occurs between page requests. It does not claim a
multi-page snapshot.

Refetch or next-page transport failures preserve confirmed cached rows. Retry
messages are local and disclosure-safe; raw URLs, paths, exception text, and
`technical_message` never enter the DOM.

## Component boundaries

### `RunDetailPage`

Owns URL decoding and callbacks for Activity, Artifacts, QC, artifact
selection, and source-artifact navigation. It performs no data fetching.

### `RunProgressPanel`

Remains the only run controller. It derives artifact and QC outcomes from the
one event snapshot, decides the one canonical polling condition, renders the
shared header/actions and tabs, and passes read-only values into each view.

### `qcState.ts`

Contains pure helpers for outcome derivation, durable metric/cursor validation,
defensive cursor continuation, page flattening, and produced-time formatting.

### `QcWorkbench`

Owns the generated infinite query, first-page reload, honest async states, and
the loaded-count summary. It does not read workspace files or mutate the run.

### `QcMetricList`

Renders a semantic fixed-layout table at medium widths and above and
border-separated vertical entries below that breakpoint. It shows exact
decimal text, display/key, unit, scope, sample/experiment, assay, backend flag,
source artifact, and produced time. Long values remain complete accessible
text and wrap without document overflow. The source action uses a lucide icon,
tooltip, and accessible name to navigate to the existing artifact inspector.

## Visual and accessibility rules

The view reuses the current gray/white workbench, teal accent, status colors,
`Button`, CSS variables, and lucide-react. It uses border separation rather
than nested cards, gradients, oversized headings, or a new UI framework.

Loading rows have stable minimum dimensions. Loading, indexing, retry, cursor,
and empty messages use polite status semantics. Table headers describe the
combined fields; mobile entries use `dl` labels. Controls remain reachable by
keyboard, and long metric IDs/keys/values cannot create horizontal scrolling at
1440x900, 1024x768, 390x844, or 360x800.

## Test and E2E boundary

Vitest and Testing Library cover outcome precedence, invalid context,
truncation, pagination, duplicate/invalid cursors, exact decimal rendering,
responsive markup, source navigation, URL history/deep links, active/pending,
confirmed empty, failure, initial error, cached-data error, retry, and
first-page reload.

The existing real Redis/RQ/Snakemake browser profile currently has no supported
QC source and therefore truthfully produces `qc_metrics_indexed` with zero
metrics. PR132's real Playwright coverage opens `?view=qc`, confirms the empty
outcome, reloads the deep link, and checks desktop/mobile layout and overflow.
PR134, by explicit milestone boundary, will extend the deterministic fixture to
produce a real non-empty QC summary and prove the complete results chain.

OpenAPI and Orval are regeneration/drift checks only. No generated file,
backend contract, package lock, test workflow output, or production scientific
target changes in PR132.
