# PR129 Run Artifact Browser Design

## Scope and intent

PR129 turns the persisted artifact references exposed by PR128 into a quiet,
repeatable run-workbench experience. It adds URL-addressable Activity and
Artifacts views, bounded artifact pagination, and a metadata inspector. It does
not add downloads, file reads, QC interpretation, backend routes, or another
source of run lifecycle truth.

The existing run progress query remains the only owner of canonical run,
event, and log polling. Artifact list and detail queries read only the generated
PR128 operations.

## Alternatives considered

1. Add a standalone artifact route with its own run query. This is visually
   simple but creates competing run polling and makes header actions diverge.
2. Move every run mutation and query into the route before adding tabs. This
   yields a pure container/presenter split but is a broad refactor of proven
   execution UI.
3. Keep `RunProgressPanel` as the run controller and workbench shell, while
   delegating artifact data and presentation to focused components. This is the
   selected approach because both views share one canonical run snapshot and
   one action surface with a bounded diff.

## URL and navigation contract

- `/runs/:runId` and `?view=activity` show Activity.
- `?view=artifacts` shows the Artifact Browser.
- `?view=artifacts&artifact=<artifact_id>` opens the inspector.
- An `artifact` parameter implies the Artifacts view even if `view` is absent,
  so copied deep links recover safely.
- Switching to Activity removes `view` and `artifact`. Selecting an artifact
  sets both parameters. Each user action pushes browser history, so refresh,
  Back, Forward, and direct navigation reproduce the same state.
- Unknown `view` values fall back to Activity without crashing.
- The decoded `artifact` value must match
  `^[A-Za-z][A-Za-z0-9_.-]{0,127}$` before the generated detail operation is
  called. Invalid, traversal-like, separator-bearing, or overlong deep links
  show a safe unavailable-link state and never issue a detail request.

Tabs use `role=tablist`, `role=tab`, `aria-selected`, stable panel IDs, and
left/right arrow keyboard movement. Run ID, status, timestamps, cancellation
notice, polling notice, failure, and Start/Cancel/Refresh controls stay above
the tab panels and remain available in both views.

## UX state map

```text
canonical run snapshot
  ├─ loading/error/missing ───────────── existing run shell state
  ├─ status != succeeded ─────────────── artifacts pending-success message
  └─ status == succeeded
       ├─ no extraction outcome ──────── indexing state + bounded canonical
       │                                  run snapshot polling
       ├─ latest outcome = failed ────── redacted failure + Refresh
       └─ latest outcome = indexed
            ├─ artifact API loading ──── stable list skeleton
            ├─ API error + cached data ─ preserve rows + Retry status
            ├─ API error, no data ────── failure state + Retry
            ├─ count == 0 + empty API ─ true empty state
            ├─ count > 0 + empty API ── inconsistent/index updating + Retry
            └─ rows present
                 ├─ no selection ────── list + selection prompt
                 ├─ detail loading ──── list + stable inspector skeleton
                 ├─ detail failure ──── list + redacted Retry state
                 └─ detail success ──── list + inspector + copy feedback
```

The latest `artifacts_indexed` or `artifact_extraction_failed` event by event
sequence determines the outcome, so an explicit retry can supersede an earlier
failure. An indexed event must report a non-negative `artifact_count`. A zero
count is considered truly empty only after the API also returns an empty first
page.

The existing run query continues at its current interval while the Artifacts
view is waiting for an extraction outcome. A single derived
`shouldPollRunSnapshot` value covers lifecycle polling and succeeded/indexing;
it controls the query interval, the 15-minute timeout timer, the paused notice,
and Refresh window reset together. This uses the same query key and observer,
not a second polling loop, and cannot poll indexing indefinitely.

An `artifacts_indexed` count must be a non-negative safe integer. Invalid event
context cannot produce an empty state. If the event snapshot is truncated and
contains no visible extraction outcome, the UI says the outcome cannot be
confirmed and offers Refresh rather than claiming indexing or emptiness.

## Component boundaries

### `RunDetailPage`

Owns only URL decoding and navigation callbacks. It passes `activeView`, the
selected artifact ID, and tab/selection callbacks into the workbench.

### `RunProgressPanel`

Remains the sole run controller: run/events/logs query, mutations, bounded
polling, header facts, action buttons, and tab shell. Activity rendering remains
unchanged in behavior. It derives the extraction outcome from its existing
event snapshot and passes that value to `ArtifactBrowser`.

### `artifactState.ts`

Contains pure, unit-tested helpers for outcome derivation, cursor validation,
page flattening/deduplication, byte formatting, and timestamp formatting. It
does not fetch or render.

### `ArtifactBrowser`

Owns TanStack Query v5 `useInfiniteQuery` for `listRunArtifacts` and `useQuery`
for `getRunArtifact`. Both call Orval generated operations directly. The list
uses a page size of 50 and an explicit Load more button. Repeated, empty, or
already-seen cursors stop pagination safely. Refetch errors retain cached pages
and show Retry instead of replacing data with an empty state.

Successful list/detail envelopes are checked for `ok`, requested `run_id`, and
selected artifact identity before rendering. A malformed success response is a
safe unavailable state, never an empty result. A failed extraction offers
“Refresh status” because PR129 has no extraction retry route; only artifact API
failures use “Retry”.

### `ArtifactList`

Renders one semantic table at `md` and above and bordered vertical entries
below `md`. Both present output type, filename, relative path, size, scope,
sample/experiment, assay, and produced time. IDs and paths stay in the DOM as
complete accessible text, use `title`, and wrap without horizontal overflow.

### `ArtifactInspector`

Renders after the list in DOM order. A responsive grid places it on the right
at desktop widths and below the list on narrower viewports. Selecting or deep
linking an artifact on a narrow viewport scrolls the inspector into view, and a
visible Back to artifacts control returns to the list without changing the URL
selection. It displays only generated whitelisted metadata, opaque URI, and
relative path. Copy buttons use lucide icons, visible tooltips through `title`,
accessible labels, and live success or failure feedback. There is no Download
control.

## Visual system and responsiveness

The browser reuses the existing white/gray workbench, teal accent, status
colors, `Button`, `Panel`, CSS variables, and lucide-react. Layout uses border
separation instead of nested cards, no gradients or marketing surfaces, and no
radius larger than the existing 8 px.

Desktop uses a compact table plus a right inspector. Mobile hides the table and
shows vertical records with the inspector underneath. The content container,
cells, IDs, paths, and URI all use `min-width: 0`, wrapping, and overflow-safe
rules. Playwright checks 1440×900, 1024×768, 390×844, and 360×800.

Loading rows and inspector placeholders have fixed minimum heights. Indexing,
copy, retry, pagination, and loading updates use `role=status` or `aria-live`
without interruptive alerts.

## Testing and E2E

Vitest and Testing Library cover URL tabs/history/deep links, extraction outcome
states, infinite pagination and repeated cursors, detail recovery, copy success
and failure, retained results on refetch failure, skeletons, empty results, and
retry actions.

The deterministic browser workflow writes its existing completion marker plus
`results/multiqc/result_manifest.tsv`. The real desktop success E2E opens the
Artifacts tab, waits for `result_manifest`, selects it, verifies the inspector,
reloads the selected deep link, and checks restoration. The same completed run
is inspected at tablet and mobile viewport sizes, with screenshots and explicit
horizontal-overflow/control checks. Cancellation coverage remains unchanged.

OpenAPI and Orval are regenerated only as a drift check. No backend, generated
contract, package dependency, workflow production target, or scientific data
fixture is changed.
