# PR140 Frontend Navigation and Workflow Discoverability Design

## Scope

PR140 makes the existing React/Vite application navigable without adding a
second frontend or inventing destinations the platform does not yet support.
It adds one semantic primary navigation region to `AppShell` and moves the
existing workflow-authoring action ahead of the developer-oriented schema
payloads on workflow detail pages.

This change does not add run history, alter URL-driven workbench or result
tabs, modify API contracts, or change validation and execution behavior.

## Information architecture

The primary navigation contains only real destinations:

- the `Workflow Platform` brand links to `/workflows`;
- `Workflows` links to `/workflows` and is current on the catalog and an exact
  workflow detail route;
- `New analysis` appears only while the URL identifies a workflow and links to
  that workflow's existing `/new-run` route. It is current only on that route.

There is no `Runs` item because the platform has no run-list route. A run
detail page still exposes the brand and `Workflows` navigation, but neither
primary item claims that the run page is its current destination.

The workflow identifier is derived from the matched URL. The shell does not
fetch workflow data, hard-code the ENCODE identifier, or create a second copy
of workflow server state.

## Workflow detail disclosure

The workflow detail panel uses this document and visual order:

1. workflow name and identifier;
2. the primary `Author inputs` action;
3. a native, default-closed `Developer schema` disclosure.

The disclosure retains the complete config, sample, and options JSON Schema in
copyable `pre` elements. Native `details`/`summary` supplies keyboard and screen
reader behavior without a custom disclosure state machine. Schema code may
scroll inside its own bounded region; it must not create page-level horizontal
overflow.

## Responsive and accessibility behavior

The header uses one semantic `nav` with a stable accessible name. Its inner
container shares the application's maximum width, wraps naturally on narrow
screens, and keeps every destination as a normal link. Current-page state uses
`aria-current="page"` plus existing button styling. Lucide icons are decorative
and their text labels remain the accessible names.

The implementation must remain operable without pointer-only interactions and
must not overlap or overflow at 1440×900, 1024×768, 390×844, or 360×800.

## Verification boundary

Unit and route tests cover link destinations, active state, arbitrary workflow
IDs, detail document order, default disclosure state, expansion, and existing
unknown-route/workflow behavior. Real Playwright tests cover desktop and mobile
navigation, first-view authoring discoverability, schema expansion, keyboard
operation, and page-level overflow. OpenAPI and Orval output must remain
unchanged.

## Non-goals and residual limits

PR140 does not provide run history, draft persistence, breadcrumbs, Agent write
actions, a second adapter, authentication, multi-user features, HPC support, or
backend changes. On run detail pages the shell cannot offer a workflow-specific
authoring link because the current route does not carry a workflow identifier;
the reliable global escape remains `Workflows`.
