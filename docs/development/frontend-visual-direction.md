# Frontend Visual Direction

## 1. Decision

**Selected direction:** Scientific SaaS Console with clinical/bioinformatics constraints.

This direction was chosen after an active-state visual exploration that produced two candidate screenshots:

- `/tmp/visual-directions-active/direction-2-saas-active.png`
- `/tmp/visual-directions-active/direction-1-clinical-active.png`

Direction 2 was selected over Direction 1 because it best balances:

- Product polish appropriate for a modern scientific platform.
- Readability and information density for validation workflows.
- Long-session comfort through a quiet light-gray canvas rather than pure white.
- A teal accent that reads as scientific/clinical without defaulting to generic SaaS purple or indigo.

## 2. Visual principles

- **Light theme only.** No dark mode in this phase.
- **Quiet light-gray canvas** with white panels/cards on top.
- **Scientific teal/blue accent.** Avoid generic purple/indigo.
- **Dense but readable validation workspace.** Inputs should be compact without feeling cramped.
- **Issue panel is the central artifact.** Validation results must be the most visually prominent outcome of a run.
- **Agent sidebar is visibly bounded and read-only.** Label it "Validation Assistant — Read Only" and keep its scope explicit.
- **No decorative gradients, orbs, or marketing landing-page elements.**
- **No fake execution state.** Only validation results, schema hints, and read-only assistant context.

## 3. Design tokens to implement later

These tokens are guidance for the follow-up implementation PR. They are not changed in this document PR.

| Token | Guidance |
|-------|----------|
| Canvas background | Quiet light gray, e.g. `#f3f4f6` or `#f8fafc` |
| Surface/panel background | White or near-white, e.g. `#ffffff` |
| Border color | Soft neutral, e.g. `#d1d5db` or `#e2e8f0` |
| Text | Near-black, e.g. `#111827` or `#0f172a` |
| Muted text | Mid gray, e.g. `#4b5563` or `#475569` |
| Accent | Scientific teal, e.g. `#0d9488` or `#0f766e` |
| Error | Strong red, e.g. `#dc2626` or `#b91c1c` |
| Warning | Amber/orange, e.g. `#b45309` or `#d97706` |
| Info | Blue, e.g. `#0369a1` or `#2563eb` |
| Radius | Small and consistent: `rounded` or `rounded-lg` |
| Shadow | Subtle, e.g. `shadow-sm` |
| Spacing density | Tight gaps between panels; compact padding inside panels |

## 4. Component polish plan

This section describes what to update in the follow-up implementation PR. The component and data-flow structure from PR87b remains unchanged.

### App shell/header
- White header with a bottom border.
- Title in teal accent, medium weight, no uppercase transformation.
- Centered max-width layout to keep content readable on wide screens.

### Panel/card treatment
- White surface, soft border, subtle shadow.
- Panel headers use muted uppercase text at small size.
- Consistent padding across all panels.

### Workflow catalog item
- Compact card with workflow name and ID.
- Selected state uses teal border and a light tinted background.
- Hover state uses a soft surface highlight.

### Schema hints blocks
- Render config/sample/options schemas in a compact multi-column grid on wide screens.
- Use a bordered preformatted block with a light background.
- Keep font small and monospaced for readability of JSON-like hints.

### Validation workspace textareas/buttons
- Three input blocks in a responsive grid.
- Compact textarea rows.
- Primary Validate button in teal with white text and a darker hover.

### Issue panel rows, badges, chips, details disclosure
- Issue panel is the visual focus.
- Each issue uses a left border accent matching severity.
- Severity badge, path chip, and code chip appear on one line when possible.
- Message in normal text; hint in muted text.
- Details disclosure remains accessible and shows technical details in a bordered pre block.

### Agent sidebar
- Visibly bounded card with its own border and background.
- Header clearly states "Validation Assistant — Read Only".
- One- or two-line scope disclaimer.
- Workflow and issue count summary below a divider.

## 5. Accessibility

- Maintain **WCAG AA contrast** for all text and interactive elements.
- Severity must **not rely on color alone**: use text labels, left border accents, and grouped structure.
- Keep **keyboard focus visible** with a clear focus ring on inputs, buttons, and disclosures.
- The **Details disclosure** must remain keyboard-operable and screen-reader friendly.

## 6. Non-goals

The following are explicitly out of scope for this visual direction and the immediate follow-up implementation PR:

- No frontend implementation in this document PR.
- No real API/fetch integration.
- No workflow execution, run lifecycle, logs, SSE, artifacts, QC, provenance, auth, or plugin marketplace.
- No real agent backend.
- No dark mode.
- No changes to the component ownership or data-flow architecture established in PR87b.
- No backend changes.

## 7. Next implementation PR

The follow-up PR will:

- Apply the tokens and spacing above to the existing frontend implementation.
- Keep the current component and data-flow structure.
- Update only presentational styling: CSS tokens, layout spacing, panel/card treatment, and component polish.
- Include no backend changes.
