# Workflow Platform Frontend UX / Product Flow Spec

This document defines the frontend user experience and product flow for the
validation-first workflow platform. It is a design and specification document
only. It does not include React/Next code, FastAPI route code, production Python
edits, design assets, or CI changes.

## What this spec is (and is not)

This spec is a contract between product, UX, UI, backend, and multi-agent
system perspectives. It describes the screens, behavior, API mapping, and
boundaries that a future frontend implementation should follow.

This spec is **not** an implementation plan. It does not prescribe component
libraries, state management, styling systems, or backend implementation. It
does not add endpoints, change platform models, or introduce execution/run
lifecycle concepts.

## File touch guardrail

The associated PR changes only this markdown file. No CI files, lockfiles,
Snakefile, workflow rule files, `docs/superpowers/`, `research/`, or
`encode-pipeline-architecture-report.md` are modified. Commit messages must not
contain attribution, tool, or vendor trailers.

## Product goal

The validation MVP lets a scientist or bioinformatics operator choose a
workflow, inspect its schema hints, enter `config`, `samples`, and `options`,
and validate those inputs before any execution concept exists. The frontend
must make structured validation issues actionable and help the user converge on
a valid input bundle. The frontend may copy or export only the user-submitted
`WorkflowInputs` payload; it must not copy, export, or branch on the
adapter-private validation response `value` unless a future public export model
is defined.

### MVP user goals

- Select the correct workflow for an assay (ENCODE-style ChIP-seq, CUT&Tag,
  ATAC-seq, or MNase-seq initially).
- Author valid `config`, `samples`, and `options` inputs without writing raw
  YAML/TSV blind.
- Validate inputs and resolve issues before any run lifecycle exists.
- Understand what each field means and what values are allowed.

## Primary user journey

1. **Choose** a workflow from the catalog.
2. **Inspect** its schema hints to learn required fields and sample columns.
3. **Enter** config values, provide a sample sheet path or TSV contents, and set
   options.
4. **Validate** the inputs.
5. **Review** structured `Issue`s mapped to fields by `path`, fix them, and
   revalidate.
6. **Copy/export** the user-submitted `WorkflowInputs` payload (`config`,
   `samples`, `options`) for CLI use or a future run-submission feature. The
   frontend must not copy, export, or branch on the adapter-private validation
   response `value` unless a future public export model is defined.

## Core screens

### 1. Workflow catalog

- Lists registered workflows from `GET /api/v1/workflows`.
- Each card or row shows `name`, `version`, `engines`, `tags`, and a
  "Validation" capability badge when supported.
- Clicking a workflow opens the workflow detail view.

### 2. Workflow detail / schema reference

- Reads `GET /api/v1/workflows/{workflow_id}/schema`.
- Displays workflow description, version, supported engines/tags, and adapter-owned
  schema hints for `config`, `samples`, and `options`.
- Primary call-to-action: **Validate inputs**.
- Schema hints are read-only reference, not a guaranteed JSON Schema form
  generator. The UI may render them as a collapsible reference panel.

### 3. Validation workspace

- A single-page editor where the user provides:
  - `config` (object, required)
  - `samples` (file path string, inline row list, or `null`)
  - `options` (object, defaults to `{}`)
- The workspace is organized into sections for config, samples, and options.
- Schema hints appear beside or below the inputs in a collapsible reference
  panel.
- Pre-fill input fields from schema hints **only** when an adapter-owned hint
  explicitly provides a default. Schema hints are not guaranteed JSON Schema,
  so the frontend must not run a generic form engine that infers defaults from
  absent hints.
- The primary action is **Validate**.

### 4. Sample input surface

- Supports two modes for the current ENCODE-style adapter:
  - **Path mode**: a string path to a sample sheet TSV on the compute/shared-storage
    environment.
  - **Paste/upload mode**: reserved for future adapter support; inline sample rows
    are currently unsupported and should surface `ENCODE_ADAPTER_UNSUPPORTED`.
- Path semantics must be explicit: paths refer to the compute environment, not
  the browser.
- Sample sheets should allow batch editing via TSV upload/paste plus row-level
  validation errors once adapters support inline rows.

### 5. Validation result panel

- Renders inline next to or below the editor after a validation call.
- Shows a pass/fail summary and the structured issue list.
- Empty success state: "No issues found — inputs are valid."
- On success, expose a **Copy request payload** action that copies the
  user-submitted `WorkflowInputs` (`config`, `samples`, `options`). The frontend
  must not copy, export, or branch on the adapter-private `value` in the
  validation response unless a future public export model is defined.

### 6. Issue panel

- The central artifact of the validation experience.
- Columns/fields: severity, source, path, message, hint, and a details action.
- Group issues first by `source`, then sort by severity: error, warning, info.
- Each issue row shows:
  - User-facing `message`
  - Severity badge with color + icon + label
  - `path` chip (for example `config.samples`)
  - Expandable `hint`
  - Hidden `technical_message` behind a **Details** disclosure
  - `code` for support, analytics, and future localization
- Clicking a path scrolls to and focuses the offending field.

### 7. Agent assist sidebar

- A narrow (280–320 px), read-only, collapsible panel that scrolls independently.
- Header label: **Validation Assistant — Read Only**.
- Capability statement on first open: "I can help you understand validation
  errors and schema requirements. I cannot run, modify, or delete workflows."
- The sidebar is a UI surface only. This spec does not define backend model
  logic, prompt engineering, or agent architecture.

### 8. Future placeholders

The following areas are intentionally left blank. They will be specified in
separate follow-up documents once the validation MVP is stable:

- Execution / run lifecycle screens
- Log streaming / SSE views
- Artifact, QC, and provenance pages
- User authentication and authorization
- Plugin marketplace or dynamic adapter loading
- DAG visualization

## UI behavior

### Loading, empty, and error states

| Screen | Empty | Loading | Error |
|---|---|---|---|
| Catalog | "No workflows registered" with guidance to check backend registration. | Skeleton rows. | Inline retry banner with API issue message and refresh action. |
| Workflow detail | Capability-unsupported notice or 404 "Workflow not found" message. | Spinner overlay on detail card. | Persistent banner with issue code and retry option. |
| Validation workspace | Blank input fields, pre-filled only when an adapter-owned schema hint explicitly provides a default. Schema hints are not guaranteed JSON Schema, so the frontend must not infer defaults from absent hints. | "Validating…" progress indicator and disabled Validate button. | Network/server or malformed-request banner with issue code and hint. |
| Results panel | "No issues found — inputs are valid." | Validating spinner overlay. | Validation-service or request-error banner. |

### Issue grouping and severity

- Group issues by `source` (`config`, `samples`, `options`, `adapter`, `api`,
  `runtime`).
- Within each source, sort: error → warning → info.
- Summary header shows total counts for each severity.
- Errors block validation success; warnings and info do not.

### Field highlighting

- When an issue row is hovered or clicked, highlight the corresponding form
  field or config path.
- Use a left border accent on issue rows matching severity color.
- Issue paths must be actionable: clicking a path scrolls to the field.

### Message, hint, and technical_message

- `message` is the primary user-facing text and is always visible.
- `hint` is shown inline as remediation help when present.
- `technical_message` is hidden behind a **Details** disclosure below the
  message.
- Expanded `technical_message` renders in a monospace code block with a
  copy-to-clipboard button.
- Disclosures are state-independent per issue.

### Fix-and-revalidate loop

- Keep inputs editable while results are visible. Do not navigate away.
- Provide a primary **Validate again** button.
- Preserve scroll position and cursor focus across revalidations.
- Animate new or changed issues in place and remove resolved ones.

### Config readability

- Use monospace typefaces for paths, config keys, file names, and messages.
- Provide syntax highlighting for YAML/JSON config views.
- Show line numbers and issue gutter markers in the config editor.
- Use sticky section headers in long forms.
- Add a search/filter bar for config keys and issue text.
- Prefer visible hints and doc links over tooltips.

## API mapping

| Frontend concept | API endpoint | Notes |
|---|---|---|
| Workflow catalog | `GET /api/v1/workflows` | Renders `metadata` + `capabilities`. Treat missing `validation` capability as unsupported. |
| Workflow detail / schema | `GET /api/v1/workflows/{workflow_id}/schema` | Adapter-owned hints. Do not auto-generate forms from them without adapter-specific heuristics. |
| Validate inputs | `POST /api/v1/workflows/{workflow_id}/validate` | Submit `WorkflowInputs`: `config`, `samples`, `options`. Read `ok`, `workflow_id`, and `issues`. The adapter-private `value` must not be consumed by the frontend. |

### Translating user input into `WorkflowInputs`

- Build a form model with three top-level fields:
  - `config` (object, required)
  - `samples` (string path, inline row list, or `null`)
  - `options` (object, default `{}`)
- Submit the top-level `samples` field. Do not duplicate it inside `config.samples`.
- The ENCODE-style adapter route copies top-level `samples` into
  `config["samples"]` before validation.
- If both top-level `samples` and `config.samples` are present in a legacy
  payload, top-level `samples` wins. The UI should only submit the top-level
  field.

### Adapter-private success payload

- `value` in the validation response is adapter-private.
- The frontend must treat `ok`, `issues`, and `workflow_id` as the stable API
  surface and must not branch on ENCODE-specific fields inside `value`.
- Any display of validated internals should be driven by adapter-issued `Issue`
  objects only.

### API contract ambiguities for future work

The following gaps were flagged by the backend/API review and should be
resolved in future contract revisions, not in this spec:

- Schema hints use `schema_kind: "hints"` and are not guaranteed to be JSON
  Schema. The frontend cannot safely auto-generate forms without
  adapter-specific heuristics.
- The schema endpoint does not indicate which `samples` input modes the adapter
  supports (string path, inline list, or both).
- `GET /api/v1/workflows` has no pagination, search, or filtering. A growing
  registry will need query parameters.
- There is no endpoint for listing supported engines, tags, or filtering by
  capability.

## Agent boundaries

For the validation MVP, the agent sidebar may:

- List workflows.
- Inspect workflow schemas.
- Submit validation requests on behalf of the user.
- Explain returned issues using structured, read-only context.

The agent sidebar must not:

- Run, kill, delete, or mutate workflows.
- Bypass services or adapters.
- Treat validation summaries as provenance.
- Infer execution state from validation responses.

### Context provided to the agent

- Current validation issue JSON.
- Relevant schema field definitions.
- User-submitted input values.
- API response metadata (endpoint, timestamp, workflow ID).

Do not include execution history, run IDs, log snippets, artifact metadata, or
inferred pipeline state.

### Future state-changing actions

Any future agent action beyond read-only validation (for example, "fix and
resubmit" or "save corrected inputs") must route through a blocking human-in-the-loop
approval gate. The UI must show a structured diff/summary and require explicit
confirmation. Default timeout behavior should be reject, not approve.

### Communicating read-only status

- Persistent header label: **Validation Assistant — Read Only**.
- Capability statement on first open.
- Disabled or hidden states for any action the agent cannot perform.

## Visual / product style

- **Tone**: clinical, utility-first SaaS. Neutral grays, one brand accent color,
  zero marketing fluff.
- **Density**: dense but readable. Compact line heights, tight padding inside
  data tables, generous whitespace only around primary actions.
- **Typography**: sans-serif for labels and navigation; monospace for technical
  values.
- **Theme**: light theme by default. Defer dark mode unless requested.
- **Accessibility**: WCAG AA contrast minimums (4.5:1 text, 3:1 UI components).
- **Severity encoding**:
  - Error: red accent + octagon/alert icon + label.
  - Warning: amber accent + triangle icon + label.
  - Info: blue accent + info-circle icon + label.
  - Never rely on color alone. Use left border accents on rows/panels.

## Non-goals

- No execution / run lifecycle UI.
- No log streaming / SSE views.
- No artifact, QC, or provenance pages.
- No authentication or authorization.
- No plugin marketplace or dynamic adapter loading.
- No agent backend implementation, model logic, or prompt engineering.
- No auto-generated form engine built solely from adapter-owned schema hints.
- No React/Next implementation, FastAPI code, or production Python edits.

## Recommended next PRs

1. **FastAPI validation MVP**: implement the PR84 API contract routes
   (`GET /api/v1/workflows`, `GET /api/v1/workflows/{workflow_id}/schema`,
   `POST /api/v1/workflows/{workflow_id}/validate`) using the existing default
   validation service composition.
2. **Minimal frontend app shell**: scaffold the application layout, routing,
   and global error boundaries needed to host the validation MVP screens.
3. **Validation UI implementation**: build the workflow catalog, detail/schema
   view, validation workspace, issue panel, and agent assist sidebar per this
   spec.
4. **API contract amendment for schema/sample-mode discovery**: if frontend
   form generation requires it, extend the API contract to clarify
   `schema_kind` stability and which `samples` input modes an adapter supports.
5. **Run lifecycle / execution UX spec**: define run submission, status,
   cancel, and retry screens once the platform supports execution.
6. **Logs / SSE UX spec**: define log streaming, filtering, and download
   behavior.
7. **Artifacts / QC / provenance UX spec**: define result browsing, QC tables,
   and provenance views.
8. **Agent confirmation boundary spec**: detail the human-in-the-loop gate for
   any future state-changing agent action.

Execution, logs, artifacts, QC, provenance, authentication, authorization,
plugin marketplace, dynamic adapter loading, and state-changing agent actions
remain future work outside the validation MVP.

## Design inputs

This spec synthesizes feedback from six agency-agent perspectives:

- **Product Manager**: focused MVP user goals and journey; deferred run
  lifecycle, dashboards, DAG visualization, and collaboration; flagged
  bioinformatics-specific needs such as path semantics, sample batch editing,
  control/input visibility, coercion warnings, draft save/load, and keeping
  transport failures separate from validation issues.
- **UX Architect**: proposed flat navigation with the workflow catalog as the
  only top-level destination; defined catalog → detail → workspace → results
  flow; specified empty/loading/error states per screen; recommended grouping
  issues by `source`, path chips, and fix-and-revalidate interaction patterns.
- **UI Designer**: recommended a clinical SaaS tone, issue table as the central
  artifact, severity encoding with color + icon + label, hidden
  `technical_message` details disclosure, and dense config readability through
  syntax highlighting, sticky headers, and search/filter.
- **Backend/API Architect**: confirmed `value` is adapter-private and the stable
  surface is `ok`/`issues`/`workflow_id`; recommended the three-field
  `WorkflowInputs` form model; flagged contract ambiguities around schema hints,
  sample input modes, pagination, and missing capability filtering.
- **Multi-Agent Systems Architect**: scoped the agent sidebar to validation
  assistance only; defined read-only context, future human-in-the-loop
  confirmation gate, execution-state hallucination guardrails, and UI labels
  that communicate read-only status.
- **Code Reviewer / Scope Guard**: flagged risks in future placeholders, agent
  sidebar framing, API mapping, agent boundaries, and next-PR recommendations;
  recommended explicit "what this spec is not" and file touch guardrail notes
  and a non-goals list covering execution, logs/SSE, artifacts/QC/provenance,
  auth, marketplace, and agent backend.
