# Artifact Runtime Boundary Decision

**Type:** Architecture Decision Record
**Date:** 2026-06-19
**Author:** YangZiChen-glitch (Kaslana)
**Status:** Approved
**Refines:** Stage 51 pause decision

## 1. Background

Stages 41-51 built an Artifact-oriented design layer: a YAML inventory cataloging
every pipeline output, a frozen `Artifact` dataclass with validation and query
helpers, and a suite of contract tests verifying bidirectional equivalence of
inventory, manifest, and output-contract documentation.

Stage 51 issued a pause decision: runtime artifact adoption was deferred because
no maintenance pain justified DAG-risky changes. The pause was correct, but it
left the boundary between Artifact and Snakemake imprecise. This ADR defines
that boundary explicitly.

## 2. What Remains Valuable

The following Artifact components are **proven, low-risk, and should be
maintained**:

| Component | Value |
|-----------|-------|
| `docs/architecture/artifact-inventory.yaml` | Canonical output catalog — answers "what does this pipeline produce?" without reading Snakemake rules |
| `workflow/lib/artifact.py` | Frozen `Artifact` dataclass, `validate_artifact()`, `load_artifacts()`, query helpers — zero DAG coupling |
| Contract tests (Stages 43-50) | Bidirectional equivalence validation of inventory ↔ manifest ↔ output-contract — catches drift without touching the DAG |
| `docs/architecture/artifact-roadmap.md` | Staged path documentation — updated by this ADR |

These components operate entirely at the **schema/documentation layer**. They
describe pipeline outputs without participating in how those outputs are
produced.

## 3. What Was Over-Scoped

The original roadmap's Stage 50+ vision was aspirational rather than
evidence-driven. The following directions were premature:

| Direction | Problem |
|-----------|---------|
| `artifact_path()` in `.smk` runtime | Would duplicate Snakemake's path resolution with no demonstrated benefit |
| Artifact-backed `rule all` rewrite | Replaces a working, readable target-helper architecture with a lateral abstraction |
| Global target resolver before any platform exists | Solves a problem the pipeline doesn't have today |
| Artifact as "stronger Snakemake abstraction" | Trigger conditions ("gated on demonstrated maintenance pain") were too vague — what pain, measured how? |

None of these ideas are wrong in principle. They are wrong **now** because the
current target-helper architecture works well, is independently readable and
testable, and has no demonstrated failure mode that an Artifact runtime layer
would fix.

## 4. Correct Boundary

```
Artifact layer     →  product/result language
                      What outputs exist.
                      What they mean.
                      What schemas they serve.

Snakemake layer    →  workflow execution engine
                      How outputs are produced.
                      What they depend on.
                      When they run.
```

**The boundary is: Artifact describes outputs. Snakemake produces them.**

Artifact must not:
- Generate rule bodies
- Own a dependency graph parallel to Snakemake's DAG
- Replace Snakemake wildcard, input, or output logic
- Rewrite `rule all` or target helpers in the current architecture without a
  future ADR

Artifact may:
- Serve as the canonical output catalog for documentation, manifest, and API schemas
- Be queried at the schema level (contract tests, manifest generation, export scripts)

## 5. Artifact Levels

This taxonomy defines how Artifact participates in the system at increasing
levels of integration. Levels 1-2 are current or near-term. Levels 2.5-4
require explicit triggers before adoption.

### Level 1: Catalog (Current — Maintain)

Artifact as a documentation artifact. The inventory YAML lists every output
with its path template, producing rule, scope, level, and assay gate. Contract
tests validate consistency with manifest and output-contract docs.

**DAG impact:** Zero.
**Maintenance:** Add one YAML entry per new output.

### Level 2: Schema Export (Near-term — Plan)

Artifact as a schema source. Export the inventory to JSON/TSV for downstream
consumption by documentation generators, release checklists, or static site
builders.

**DAG impact:** Zero.
**Examples:** `scripts/export_run_view.py` reading inventory plus manifest/QC
summaries to produce per-run JSON for a static SPA.

### Level 2.5: Backend Target Resolver (Future — Gated on platform need)

A thin function that maps an `artifact_id` plus concrete `sample_id` or
`experiment_id` to a filesystem path. This is a **lookup**, not target
generation.

```python
def resolve_artifact_path(artifact_id, *, sample=None, experiment=None):
    """Return the concrete filesystem path for a known artifact."""
    artifact = artifacts_by_id[artifact_id]
    template = artifact.path_template
    if sample:
        template = template.replace("<sample>", sample)
    if experiment:
        template = template.replace("<experiment>", experiment)
    return template
```

This is target resolution, not target expansion or DAG generation.

A FastAPI backend can use this to answer "where is the output for experiment X?"
and call `snakemake <target>` externally. The resolver does not run inside
Snakemake and does not influence the DAG.

**DAG impact:** Zero — runs outside Snakemake.
**Trigger:** A FastAPI/Flask backend exists and needs to serve artifact metadata.

### Level 3: Target Generation (Distant — Gated on extraordinary evidence)

Artifact entries with `rule_all: true` are expanded into `rule all` targets.
This replaces per-assay target helpers with an inventory-driven expansion.

**DAG impact:** High — changes `rule all` construction.
**Trigger:** Demonstrated, measurable maintenance pain with the current
target-helper architecture that cannot be resolved by incremental refactoring.
Requires explicit team decision.

### Level 4: Rule/Dependency Generation (Avoid)

Artifact entries define rule bodies, input/output chains, and dependencies.
This creates a parallel workflow engine alongside Snakemake.

**Status:** **Avoid.** No plausible trigger justifies this complexity.
Snakemake already handles dependencies, wildcards, and scheduling. Artifact
should not duplicate these responsibilities.

## 6. CLI-Only Recommendation

For the current CLI-only pipeline:

- Maintain `artifact-inventory.yaml` and contract tests.
- Add entries for new outputs as they are created (consensus peaks, ATAC IDR
  outputs, etc.).
- Continue the `paths.smk` pattern for assay-by-assay path centralization.
- Do NOT adopt `artifact_path()` in `.smk` runtime.
- Do NOT rewrite target helpers or `rule all` around artifacts.

## 7. FastAPI Platform Recommendation

When a FastAPI backend is built:

- `artifact-inventory.yaml` becomes the API schema seed: each entry maps to an
  API endpoint parameter or response field.
- The `Artifact` dataclass becomes a seed for Pydantic/ORM schemas, not a
  complete database table by itself.
- A thin `resolve_artifact_path()` function enables the backend to answer
  "where is output X for experiment Y?" and call `snakemake` externally.
- Snakemake remains the executor. The backend reads Artifact metadata but does
  not drive Snakemake's internal logic.

## 8. Explicit Non-Goals

The following are excluded from the current Artifact scope unless a future ADR
explicitly reverses or narrows this decision:

1. **No rule body generation.** Artifact entries describe what a rule produces,
   not how. Shell commands, container images, and thread counts belong to
   Snakemake rules.
2. **No dependency graph modeling.** Snakemake's DAG is the single source of
   truth for dependencies. Artifact does not encode `input:` chains.
3. **No replacement of Snakemake wildcards.** Path templates in the inventory
   use `<sample>` and `<experiment>` as documentation placeholders, not as a
   parallel wildcard system.
4. **No automatic `rule all` rewrite now.** Target helpers remain the
   DAG-integration mechanism. Level 3 is a distant possibility, not a plan.

## 9. Revised Roadmap

| Priority | Activity | Timeframe |
|----------|----------|-----------|
| P0 | Release hardening and scientific features | Now |
| P0 | Maintain inventory + contract tests | Continuous |
| P1 | `export_run_view.py` / static SPA for per-run output browsing | When release is stable |
| P2 | Assay-by-assay `paths.smk` adoption beyond MNase | Per-assay, as maintenance pain warrants |
| P3 | FastAPI backend with thin `resolve_artifact_path()` | When platform project exists |
| Future ADR only | Level 3 target generation | Requires explicit ADR |
| Avoid | Level 4 rule/dependency generation | No plausible trigger |
