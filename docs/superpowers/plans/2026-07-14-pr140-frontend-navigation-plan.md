# PR140 Frontend Navigation and Workflow Discoverability Implementation Plan

**Goal:** Make the existing application globally navigable and place workflow
authoring ahead of raw developer schemas without changing validation, run, or
API behavior.

## Task 1: Lock navigation and disclosure behavior

- [x] Verify PR139 ancestry, clean isolated worktree, and protected root
  worktree state.
- [x] Define the only real global destinations and mutually exclusive active
  states.
- [x] Choose URL-derived workflow context and a native schema disclosure.
- [x] Add failing component and route tests for links, active state, arbitrary
  workflow IDs, document order, and disclosure behavior.

## Task 2: Implement the minimum shell navigation

- [x] Make the brand a normal link to `/workflows`.
- [x] Add a semantic primary navigation with `Workflows` and conditional
  `New analysis` destinations.
- [x] Apply existing button styles, Lucide icons, visible focus behavior, and
  `aria-current` without introducing another layout or client.
- [x] Keep the header wrapping and page width safe on desktop and mobile.

## Task 3: Make authoring discoverable

- [x] Separate workflow identity from developer schema presentation.
- [x] Render `Author inputs` before the schema disclosure in document order.
- [x] Keep all three schemas available in a default-closed, keyboard-operable
  disclosure with bounded code overflow.
- [x] Preserve legacy validation, run creation, Agent, and route behavior.

## Task 4: Verify the user experience

- [x] Run focused and complete frontend tests, typecheck, and production build.
- [x] Run real desktop/mobile Playwright at all four target viewports and
  inspect screenshots for order, overlap, clipping, and horizontal overflow.
- [x] Confirm OpenAPI export and Orval regeneration have zero drift.
- [x] Run `git diff --check` and confirm no unrelated or root-worktree files
  changed.

## Task 5: Review and deliver

- [ ] Perform one independent exact-HEAD read-only adversarial review.
- [ ] Fix every blocking/important finding and rerun affected checks.
- [ ] Commit one scoped change, push it, and open Draft PR140.
- [ ] Require exact local/remote/PR/CI HEAD agreement and all CI checks green.
- [ ] Keep the PR Draft and do not merge.
