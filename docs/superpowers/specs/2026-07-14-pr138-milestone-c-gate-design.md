# PR138 Milestone C Input-to-Results Gate Design

## Objective and boundary

PR138 closes Milestone C with one non-mocked browser path:

```text
workflow -> author config and inline samples -> validate -> immutable snapshot
  -> create -> preflight -> explicit start -> RUNNING -> SUCCEEDED
  -> persisted QC -> persisted artifact -> exact safe download
```

The path uses the existing FastAPI, SQLite, Redis/RQ worker, Snakemake,
OpenAPI/Orval client, React workbench, and results views. It does not add an API,
lifecycle state, scientific rule, parser, artifact/QC contract, or frontend
state engine. Scientific correctness remains covered by the scientific workflow
tests; this gate proves platform integration and state truthfulness.

## Baseline findings and selected approach

PR134 already provides a sentinel-owned deterministic workflow that produces
QC, artifacts, exact download bytes, empty QC, malformed QC, and a cancellable
process tree. PR136 proves schema-driven editing and URL history. PR137 proves
real validation, snapshot creation, and durable run creation, but stops at
PLANNED. The PR134 results test still creates runs through the legacy workflow
detail form, so no single gate currently spans authoring through results.

PR138 will reuse the one controlled project and real service supervisor. The
results test will create every run through `/workflows/{id}/new-run`, the real
adapter schema, inline TSV rows, the validation snapshot, snapshot-only create,
preflight, and explicit start. This is selected over creating another fixture,
another Playwright stack, or changing the scientific tiny workflow.

## Input contract and path privacy

The controlled config documents become true authoring configs: they omit the
`samples` key. The browser reads the fixture TSV and submits its content as
inline rows through the existing workbench. The runtime manifest exposes the
controlled sample file only to Playwright's file chooser; the HTTP payload does
not contain that browser-side path.

This matches the adapter contract in which config, samples, and options are
sibling inputs. Workspace materialization remains responsible for writing
`config/samples.tsv` and adding the workspace-relative config reference. Tests
will inspect validate requests and require no `samples` config key or fixture
sample path while still requiring the complete row content.

## Truthful lifecycle and results assertions

The controlled success task emits its entered log marker and waits briefly
before writing outputs. The browser must observe canonical RUNNING while the
worker owns the run, then canonical SUCCEEDED. It may not accept an immediate
terminal transition as proof of the intermediate state.

The same run must then expose eight persisted QC metrics, navigate from a
metric to its persisted source artifact, and download bytes identical to the
workspace source. Reloaded QC and artifact deep links must recover from
SQLite-backed APIs. Existing separate controlled runs retain these semantics:

- empty: SUCCEEDED plus `qc_metrics_indexed(metric_count=0)` before a true empty
  message;
- malformed: SUCCEEDED plus a redacted QC indexing failure;
- cancel: HTTP 202 remains RUNNING until worker acknowledgement, then CANCELLED,
  with helper and child processes gone and no output markers.

The authoring path also exercises a structurally valid but scientifically
invalid config (`stage5=true`, `stage4b=false`). Backend validation must return
no snapshot and the browser must recover after the user corrects the draft.
This keeps the adapter validator authoritative without inventing client-side
scientific rules.

## Browser, history, layout, and cleanup

The existing desktop project owns the complete success/results path and checks
1440x900, 1024x768, 390x844, and 360x800 result layouts. The mobile project owns
the complete authoring-to-cancellation path. Existing workbench tests continue
to prove invalid YAML, malformed/manual sample handling, Back/Forward draft
retention, stale validation invalidation, and safe reload reset.

Screenshots must show non-empty editable controls and results without document
horizontal overflow, clipping, overlap, or inaccessible actions. The existing
owner sentinel, unique queue, process-session cleanup, and global teardown stay
authoritative. A successful test run must remove its temporary runtime and
leave no Redis, worker, horse, Snakemake, Vite, or fixture helper process.

## Local demo

`--input-authoring-demo` becomes a documented alias for the existing controlled
results demo. It retains the same sentinel-owned project, durable runtime root,
unique queue, bounded input manifest, readiness checks, and Ctrl-C cleanup.
The legacy `--results-visibility-demo` spelling remains supported. The input
manifest's config also omits `samples`; a user imports the separate TSV in the
workbench and completes the same validated snapshot path as Playwright.

No interpreter, home directory, or environment-specific path is committed.

## Roadmap and acceptance

The maintained roadmap will distinguish historical architecture phases from
delivered platform milestones and state that Milestone C is complete only on a
revision containing this gate with green required CI. PR139 remains the next
adapter-conformance foundation; no real second adapter enters PR138.

Acceptance requires focused and full Python tests, the real Redis/RQ/SIGALRM/
cancellation/tiny Snakemake gate, frontend tests/typecheck/build, both
Playwright projects with zero skips, screenshot inspection, OpenAPI/Orval zero
drift, Ruff/format/diff/hardcoded-path checks, clean worktree and process state,
green exact-HEAD GitHub checks, and one independent adversarial read-only
review with no blocking or important findings.

## Residual limits and non-goals

The fixture proves platform integration, not biological validity, throughput,
or real-data execution. Download retains the documented no-checksum risk for a
same-size same-metadata replacement. Browser drafts remain intentionally
memory-only. Authentication, multi-user isolation, HPC, object storage,
immutable bundles, FASTQ upload, omics-intake runtime coupling, another real
adapter, and Agent write actions remain outside this PR.
