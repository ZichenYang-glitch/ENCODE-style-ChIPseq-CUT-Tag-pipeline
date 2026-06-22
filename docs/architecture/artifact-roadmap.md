# Artifact Roadmap

This document defines the staged path from the current target-helper
architecture toward an artifact-oriented design. Each stage has a clear goal,
scope, non-goals, and trigger condition. The trigger condition is the signal
that says "this stage is now the highest-value next step" — it prevents
premature abstraction.

**Boundary (Stage 56 ADR):** Artifact is a result semantics / manifest /
frontend-backend schema layer. Snakemake is the workflow execution engine.
Artifact describes outputs; Snakemake produces them. A thin backend target
resolver (Level 2.5: `artifact_id` + sample/experiment → concrete path) is
allowed for future FastAPI consumption outside Snakemake. No `artifact_path()`
in `.smk` runtime, no `rule all` rewrite, no rule body generation, no
dependency graph ownership. See:
`docs/superpowers/specs/2026-06-19-artifact-runtime-boundary-adr.md`.

---

## Current Baseline (post-Stage 52 / pre-v0.2.0)

**Architecture:** `workflow/Snakefile` is the entry point: config loading,
sample sheet parsing, `include:` directives, and `rule all`. Derived metadata,
QC gating, genome resource helpers, and dispatch wrappers live in
`workflow/rules/metadata.smk` (Stage 42). Target-helper functions used by
`rule all` live in `workflow/rules/targets.smk` (Stage 42). MNase path helpers
live in `workflow/rules/paths.smk` (Stage 44). The `workflow/lib/artifact.py`
module owns the docs/tests-only `Artifact` dataclass, loader, and query helpers
(Stage 45/48). Executable rules live in `workflow/rules/common.smk`,
`peaks.smk`, `mnase.smk`, `replicates.smk`, `idr.smk`, `qc.smk`, `report.smk`.

**What works well:**
- Each target helper is short, independently readable, and gates on derived
  sample lists (`PEAK_SAMPLE_IDS`, `MNASE_SAMPLE_IDS`, etc.).
- Adding a new output type requires touching ~4 places (target helper, rule
  output, pipeline_done, manifest), which is manageable.
- Four assays are supported (ChIP-seq, CUT&Tag, ATAC-seq, MNase-seq) with
  clear non-peak-centric separation for MNase.
- Per-rule Conda environments provide reliable tool isolation.
- Artifact contract test infrastructure (Stages 43-50) validates bidirectional
  equivalence of inventory ↔ manifest ↔ output-contract without DAG changes.

**What shows strain:**

1. **Output path strings are still duplicated for non-MNase assays.** MNase paths
   are centralized in `paths.smk`, but ChIP-seq, CUT&Tag, and ATAC paths
   remain duplicated across target helpers, rule outputs, manifest, and docs.

2. **Artifact contract tests are mature, but runtime artifact adoption is paused**
   because no maintenance pain justifies DAG-risky changes (Stage 51 decision
   record). The contract tests catch drift without touching the DAG.

**These problems are not automatically solved by a full artifact-centric
rewrite.** An Artifact dataclass with `artifact_path()` would consolidate
paths in code, but the blast radius (target helpers, `rule all`, manifest
generator) is high and no test would validate it better than the existing
stress tests already do. The contract test infrastructure proves equivalence
without adopting runtime artifactization.

---

## Staged Path

### Stage 41: Artifact Readiness / Release Stabilization (implemented 2026-06-06)

**Goal:** Documentation-only — archive old report, create this roadmap, add
developer checklist, output-contract consistency pass.

**Scope:** `docs/archive/`, `docs/architecture/`, `docs/developer/`,
`docs/output-contract.md` consistency check.

**Non-goals:** Artifact dataclass, AssayPolicy YAML, artifact_path(),
global target resolver, DAG changes, Snakefile refactoring.

**Trigger:** Always — this is baseline readiness.

---

### Stage 42: Extract metadata.smk and targets.smk

**Goal:** Reduce `workflow/Snakefile` responsibility overload by extracting
derived metadata computation and target-helper functions into separate files,
with **zero behavior change** and **zero DAG change**.

**Scope:**
- `workflow/rules/metadata.smk` — derived sample lists (`PEAK_SAMPLE_IDS`,
  `MNASE_SAMPLE_IDS`, `MULTI_BIOREP_EXPERIMENTS`, etc.), QC config gating,
  genome resource helpers, and the three dispatch functions.
- `workflow/rules/targets.smk` — the ten `_*_targets()` functions and
  `_manifest_dependency_targets()`.
- `workflow/Snakefile` — reduced to config loading, sample sheet parsing,
  `include:` statements, and `rule all:` definition (~200 lines down from
  ~1,000).

**Non-goals:** Path consolidation, Artifact dataclass, AssayPolicy YAML,
global target resolver, target-helper logic changes, rule changes.

**Trigger:** When the next feature requires touching the Snakefile in a way
that makes its size unmanageable, or when a new contributor struggles to
navigate the single-file architecture.

---

### Stage 43: Artifact Inventory

**Goal:** Produce a single machine-readable inventory of every output type
the pipeline can produce, including assay eligibility, gating conditions,
file path template, and generating rule.

**Scope:** A YAML file under `docs/` that lists all output types with their
metadata. No code changes — the inventory is hand-maintained documentation
in machine-readable format.

**Non-goals:** Code generation from the inventory, runtime artifact resolution,
DAG behavior changes.

**Why inventory before paths.smk:** You can't centralize what you haven't
catalogued. `paths.smk` needs to know all path patterns before it can
centralize them. Without an inventory, `paths.smk` would grow organically
as new outputs are added — path centralization without path design. The
inventory provides the reference vocabulary, identifies which paths benefit
most from centralization (e.g., MNase fragment BAMs sharing a
`03_fragments/<sample>.<class>.bam` pattern), and is zero-risk (a docs
file can be wrong without breaking the DAG). `paths.smk` is then scoped
by the inventory rather than grown ad-hoc.

**Trigger:** When the output-contract.md table becomes too large to maintain
by hand, or when a manifest / output-contract inconsistency is found.

---

### Stage 44: Introduce paths.smk (Limited Pilot)

**Goal:** Create a single source of truth for output path templates, starting
with ONE assay or output category as a pilot. Not a global all-at-once path
rewrite.

**Scope:**
- `workflow/rules/paths.smk` — a thin module with functions like
  `mnase_fragment_bam(sample_id, class_name)`, `mnase_signal_bw(sample_id, kind)`.
  Each function returns a path string. Used by the MNase rules and target
  helper only.
- Existing paths in ChIP-seq / CUT&Tag / ATAC rules are NOT migrated.
- Used by exactly one target helper (`_mnase_targets()`) and its
  corresponding rules as a pilot.

**Non-goals:** Universal path resolver, Artifact dataclass, migrating all
assays, changing path formats.

**Trigger:** When a path change (e.g., renaming `03_fragments` to
`04_fragments`) requires touching 4+ files and the pilot assay is stable
enough to benefit from path centralization.

**Warning:** `paths.smk` must NOT become a hidden Artifact model. It is a
path utility module, not an abstraction layer. Do not add assay dispatch,
artifact kind enums, or `artifact_path()` semantics.

---

### Stage 45: Artifact Dataclass Spike (implemented 2026-06-06)

**Status:** Implemented as `workflow/lib/artifact.py` (frozen dataclass,
loader, validation helpers) with `test/test_stage45_artifact_model.py`
(22 checks). The Snakemake DAG and rules are NOT changed.

**Goal:** Introduce a `workflow/lib/artifact.py` with a frozen `Artifact`
dataclass used ONLY by docs/tests validation. The Snakemake DAG, rules,
`make_manifest.py`, and manifest generation are NOT changed.

**Scope:**
- Frozen 13-field `Artifact` dataclass matching `artifact-inventory.yaml`.
- `validate_artifact()` for raw dict validation, including malformed-input hardening.
- `load_artifacts()` loader for docs/tests validation.
- `test/test_stage45_artifact_model.py` with 22 checks.
- Does NOT change `make_manifest.py` or manifest generation.

**Non-goals:** Snakemake rule `input:` / `output:` using artifacts,
target-helper replacement, rule-all changes, path consolidation.

**Trigger:** When the artifact inventory (Stage 43) is stable and manually
maintaining manifest + docs + inventory in sync becomes a measurable burden,
AND Stages 42-44 (Snakefile extraction + inventory + path pilot) have proven the
extraction architecture works.

---

### Stage 46: Artifact-Backed Inventory Tests (implemented 2026-06-07)

Test-layer adoption: `test_stage43_artifact_inventory.py` now uses
`load_artifacts()` and `Artifact` from `workflow/lib/artifact.py`
as its canonical inventory read path. Schema completeness is derived
from `dataclasses.fields(Artifact)` — no duplicate `REQUIRED_FIELDS`.
No DAG/runtime changes.

---

### Stage 47: MNase Path-Helper Inventory Contract Tests (implemented 2026-06-07)

Test adoption: `test/test_stage47_mnase_path_contract.py` verifies all 7
`paths.smk` helpers produce path strings matching all 13 MNase entries in
the artifact inventory. No DAG/runtime changes.

---

### Stage 48: Artifact Query Helpers (implemented 2026-06-07)

Added `artifacts_by_id()`, `artifacts_by_manifest_output_type()`, and
`filter_artifacts()` to `workflow/lib/artifact.py`. Refactored Stage 43
and Stage 47 contract tests to use them. Stage 45 model test extended
with 7 query-helper checks (22→29). No DAG/runtime changes.

---

### Stage 49: Manifest Artifact Contract Tests (implemented 2026-06-07)

Test adoption: `test/test_stage49_manifest_artifact_contract.py` verifies
manifest output_type vocabulary matches the artifact inventory via
`artifacts_by_manifest_output_type()`. 8 checks. No DAG/runtime changes.

---

### Stage 50: Output Contract Generation Dry-Run (implemented 2026-06-07)

Dry-run verification: `test/test_stage50_output_contract_dry_run.py`
proves bidirectional equivalence between `docs/output-contract.md` and
artifact inventory, including path_template matching with `<exp>` →
`<experiment>` normalization. No docs generation adoption yet.

---

### Stage 51: Artifact Adoption Decision Record (implemented 2026-06-07)

Decision record: `docs/architecture/artifact-adoption-decision.md`.
**Pause runtime artifact adoption.** Contract test infrastructure (Stages 41-50)
is complete and self-sustaining. Return to release hardening and scientific
features. Target-helper adoption remains gated on demonstrated maintenance pain.

---

### Stage 52: Release Readiness Sweep (implemented 2026-06-07)

Documentation/test-readiness sweep after Stages 40-51. Updated README
(pre-release description, repo layout, test list), CHANGELOG (Stage 41-51
entries), KNOWN_ISSUES (Current Status), configuration.md (MNase config
expansion), assay-policy.md, and output-contract.md (pooled MNase gating).

---

### Stage 53+: Artifact-Assisted Target Helpers (gated — only if maintenance pain emerges)

**Goal:** Target-helper functions use artifact definitions to expand targets,
reducing per-output boilerplate.

**Scope:** `_build_targets_for(sample_ids, artifact_kinds)` replaces repeated
`expand(...)` calls with a loop over the artifact registry.

**Non-goals:** Global target resolver, rule-all replacement, assay dispatch
changes.

**Trigger:** When target-helper repetition is the top maintenance pain point
AND stages 41-45 documentation + extraction + inventory are complete AND
the extraction architecture has been stable for multiple releases.

---

### Stage 53+: Global Target Resolver (only if explicitly justified)

**Goal:** A single `build_run_targets()` function replaces the target-helper
functions. `rule all: input:` calls it via a lambda.

**Trigger:** ONLY if ALL of:
- A fifth assay type is added and dispatch becomes fragile.
- The artifact inventory has been stable for multiple releases.
- Stages 42-46 have proven the extraction + artifact model.
- The team explicitly decides the resolver complexity is justified.

**If the trigger is never met, the current extracted architecture is
perfectly acceptable.**

---

## Explicit Warnings

These approaches are tempting but WRONG for the current project state:

1. **Full artifact-centric rewrite now.** The pipeline has 4 assays working
   correctly. A rewrite would touch every rule, every test, and every doc
   simultaneously — high regression risk with zero user-visible benefit.
   The right path is incremental extraction, inventory, and path piloting
   (Stages 42-44) followed by targeted abstraction (Stages 45-46+), each
   stage independently verifiable.

2. **Global target resolver before extraction.** A `build_run_targets()`
   function that replaces all ten target helpers at once is a single massive
   PR. If it breaks, every assay breaks. Extract first, consolidate later.

3. **All-path migration in one PR.** Migrating every output path in the
   pipeline to a centralized `paths.smk` in a single change is high-risk
   and hard to review. Pilot with one assay (MNase), validate, then expand.

4. **Using `paths.smk` as a hidden Artifact model.** `paths.smk` is a path
   utility module. Do not add `Artifact.kind` enums, assay dispatch logic,
   or `artifact_path()` semantics to it. When the time comes for artifact
   abstractions, they go in a separate `workflow/lib/artifact.py` module
   with clear boundaries.

## Hard Gates

These must NOT be introduced before their trigger condition is met:

| Item | Earliest Stage | Trigger |
|:---|:---|:---|
| Artifact dataclass | 45 | Implemented |
| AssayPolicy YAML | 47+ | Fifth assay or dispatch fragility |
| `artifact_path()` | 47+ | Not implemented; only after artifact-backed tests/model prove useful |
| Global target resolver | 53+ | Explicit team decision |
| `rule all` rewrite | 53+ | Resolver proven in tests |
