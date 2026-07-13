# PR136 Schema-Driven Input Workbench Design

## Objective and boundary

PR136 adds an adapter-driven browser workbench at
`/workflows/:workflowId/new-run`. It renders the versioned PR135 authoring
contract, lets a user edit config in a schema form or advanced YAML, imports
and edits inline sample rows, edits adapter options, and reviews one stable
in-memory draft.

The workbench does not call validation, create a run, preflight, or start
execution. The UI labels its output as a client-side structural draft and
continues to treat the backend scientific validator as authoritative. PR137
will define immutable validated snapshots and submission.

## Contract loading and compatibility

The route calls the generated `getWorkflowSchema` operation directly from a
TanStack Query. It does not use the existing compatibility client because that
adapter intentionally discards version, dialect, coverage, modes, and limits.

The first implementation supports exactly schema contract `1.0.0` and JSON
Schema Draft 2020-12. It fails closed with a controlled unsupported-contract
state when:

- the response has no schema;
- the version or dialect is unsupported;
- the advertised authoring modes needed by the workbench are absent; or
- config/options are not object schemas or samples is not an array-of-object
  schema.

Loading uses a stable-sized placeholder. Network and contract failures expose
only controlled messages and an explicit retry action. No exception text,
request URL, filesystem path, or server technical detail is rendered.

## Draft and URL state

One reducer owned by the workbench route contains all config, YAML buffer,
sample rows, options, and local structural issues. Individual panels are
projections and never own authoritative input state.

The active section is URL-driven with
`?step=config|samples|options|review`. Config mode is represented by
`mode=form|yaml` when the config step is active. Tab changes push browser
history, so deep links, Back, and Forward restore the view without discarding
the reducer while the route remains mounted.

The draft is deliberately memory-only. No value is written to localStorage,
sessionStorage, a URL, or the backend. A hard reload retains the route/section
but creates a fresh draft from schema defaults; the UI says that reload clears
the current-page draft. This avoids silently retaining input paths or sample
metadata under a stale schema version.

## Canonical config and options

The canonical config value is one plain JSON object. RJSF is controlled by
that object and updates it directly. `omitExtraData` and `liveOmit` remain
disabled, so form edits preserve advanced keys not covered by the partial
schema.

Advanced YAML is a transient editing buffer:

1. Entering YAML regenerates the buffer from the latest canonical config.
2. A valid object-root YAML document replaces the canonical config.
3. Invalid YAML remains visible with a controlled line/column issue while the
   last-known-good config is retained internally.
4. Review is blocked while the YAML buffer is invalid, and switching back to
   Form is refused until the document parses. The UI never reviews a stale
   object as though the invalid buffer were accepted.
5. After a later Form edit, the next YAML entry regenerates from canonical
   state. Formatting and comments may be normalized; values and unknown keys
   are retained.

YAML parsing uses `yaml.parseDocument` with the core schema, strict duplicate
key checking, and a bounded alias expansion. Multiple documents, custom tags,
array/scalar/null roots, cyclic aliases, NaN/Infinity, `undefined`, and other
non-JSON values fail closed. User-visible errors are normalized to a short
message and optional line/column only.

Options use a second controlled RJSF form and canonical JSON object. RJSF's
Draft 2020-12 validator is constructed with
`customizeValidator({ AjvClass: Ajv2020 })`; the default Draft-07 validator is
not used. Initial values come from RJSF's default-form-state utility so visible
defaults and Review cannot diverge. Tests compile and exercise all three real
adapter schemas because RJSF documents Draft 2020-12 as non-default support.

## Sample TSV and table model

Column order, labels, required fields, enum choices, and descriptions are
derived solely from `sample_schema.items.properties`. The platform contains no
ENCODE field list and does not infer assay, replicate, control, or other
scientific values.

Papa Parse reads tab-delimited files in array mode with dynamic typing off.
Array mode is intentional: duplicate headers must be rejected instead of being
silently renamed. A successful import is atomic and replaces the current rows
only after every check passes. The importer:

- rejects a file larger than `limits.max_request_bytes` before reading;
- accepts quoted fields and CRLF while preserving exact strings and empty
  cells;
- rejects an empty/header-only file, duplicate or overlong headers, missing
  required headers, unknown headers, ragged rows, and parser errors;
- rejects embedded tab/newline/NUL values because the backend inline-row
  contract does too;
- enforces row, column, and cell limits using Unicode code-point length;
- fills absent optional columns with `""`; and
- normalizes every row into adapter schema order.

Only `File.name` is shown; browser file paths are neither available nor sent.
Failed imports leave the previous rows unchanged.

TanStack Table provides the desktop row/column model. Rows carry a client-only
stable ID that is stripped during Review. Editing, adding, and deleting update
the parent reducer; every cell remains a string. The desktop table lives in a
bounded overflow container. Below the desktop breakpoint the same rows render
as border-separated vertical records, not a squeezed wide table.

## Review and structural readiness

Review serializes exactly:

```json
{
  "config": {},
  "samples": [],
  "options": {}
}
```

Config and option keys use stable lexical ordering; sample rows retain visible
row order and adapter column order. The client rechecks JSON safety and the
UTF-8 request byte ceiling with `TextEncoder`. Structural errors, including
zero rows where the schema requires at least one, keep Review visibly not
ready. Nothing is described as scientifically valid and no submission control
is present.

## Interaction and presentation

The workbench is one quiet work surface with a persistent workflow/schema
header, compact accessible tabs, and separated editor regions. Radix Tabs owns
keyboard tab behavior; React Router owns the URL. Config Form/YAML is a compact
segmented control. Existing Button, Panel, CSS variables, typography, and
Lucide icons are reused.

Desktop uses the dense table and stable editor height. Mobile uses vertical
sample records and full-width controls. Containers are `min-width: 0`, long
values wrap, and CodeMirror receives an accessible label. Loading, errors, and
dynamic messages use stable regions with `role=status` or `aria-live` as
appropriate. The design adds no nested card stack, hero, gradients, or global
visual rewrite.

## Dependencies

- `@rjsf/core`, `@rjsf/utils`, and `@rjsf/validator-ajv8` at one pinned minor;
- direct `ajv` for the Draft 2020-12 class;
- `yaml` for safe parse/stringify;
- `papaparse` and its TypeScript declarations;
- `@tanstack/react-table` for the editable row model;
- `@uiw/react-codemirror`, `@codemirror/lang-yaml`, and the CodeMirror view
  package for the accessible YAML editor; and
- `@radix-ui/react-tabs` for accessible tab interaction.

React Hook Form and Zod are not added in PR136: there is no platform-owned
submission form, while RJSF/AJV already owns the two adapter schema forms.
Adding a second form state engine would not remove implementation logic in this
scope.

## Tests and acceptance

Pure tests cover schema compatibility, Draft 2020-12 canaries, YAML safety and
state transitions, deterministic Review serialization, TSV quoting/CRLF,
duplicate/missing/unknown headers, empty files, ragged rows, and every relevant
limit. Component/route tests cover loading, null schema, failures and retry,
unknown-field retention, invalid YAML, sample editing/add/delete, stable column
order, Review readiness, URL deep links, and Back/Forward draft retention.

Real Playwright coverage uses the actual FastAPI schema endpoint. Desktop and
mobile flows import a deterministic TSV, edit rows, edit advanced YAML, inspect
Review, navigate history, and prove CodeMirror is non-empty and interactive.
Screenshots at the existing desktop/mobile viewports are inspected and the
page asserts no document-level horizontal overflow. Hard reload is tested as a
safe fresh draft, not falsely as persistence.

Acceptance also requires all frontend tests, typecheck/build, existing real
browser gates, Python regression, applicable Redis/RQ/Snakemake gates,
OpenAPI/Orval zero drift, Ruff/format/diff checks, a clean worktree, green CI,
and an independent read-only merge-gate review.

## Non-goals and residual limits

PR136 does not validate scientifically, persist an immutable snapshot, create
or start a run, modify lifecycle/worker/artifact/QC behavior, upload FASTQ,
couple to omics-intake, add another adapter, or add Agent writes, auth,
multi-user isolation, HPC, object storage, or immutable workflow bundles.

The config schema is intentionally partial, so advanced keys are editable in
YAML but not necessarily rendered as specialized controls. Drafts are lost on
hard reload by design. Client structural validation is advisory and may still
be rejected by the backend scientific validator in PR137.
