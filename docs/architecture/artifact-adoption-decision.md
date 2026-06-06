# Artifact Adoption Decision Record

**Date:** 2026-06-07
**Stage:** 51
**Status:** Decision recorded. No runtime artifact adoption.

## Context

Stages 41 through 50 built a layered artifact contract test infrastructure:

| Stage | Description | What was proven |
|:---|:---|:---|
| 41 | Artifact readiness roadmap | Architecture plan |
| 42 | Snakefile extraction | Metadata/targets can be extracted from Snakefile |
| 43 | Artifact inventory | 62-entry YAML inventory with 13 fields per entry |
| 44 | MNase paths.smk pilot | 7 thin path helpers produce correct paths for MNase |
| 45 | Artifact dataclass | Frozen dataclass + loader + malformed-input validation |
| 46 | Inventory tests use model | `load_artifacts()` as canonical read path |
| 47 | MNase path-helper contract | All 7 helpers match all 13 MNase inventory entries |
| 48 | Query helpers | `artifacts_by_id()`, `artifacts_by_manifest_output_type()`, `filter_artifacts()` |
| 49 | Manifest contract | Manifest output_type vocabulary == inventory manifest types |
| 50 | Output-contract dry-run | Bidirectional equivalence proven (62/62 types match, paths match with `<exp>` normalization) |

**Key finding:** The artifact inventory fully describes the pipeline's output surface. Every output type, path template, and manifest entry is consistent across inventory ↔ manifest ↔ output-contract. But this consistency was achieved through manual testing, not automated generation — and the contract tests are read-only validators, not generators.

## Decision: Pause Runtime Artifact Adoption

**Do NOT proceed to target-helper or `rule all` runtime artifact adoption at this time.**

Rationale:

1. **No demonstrated maintenance pain.** The contract tests prove equivalence exists, but no one has reported drift between inventory, manifest, or output-contract. The existing Stage 43/49/50 tests would catch drift immediately if it occurred.

2. **Blast radius is high.** Artifact-assisted `_mnase_targets()` or `rule all` would touch the DAG's demand-driven resolution. A bug here would silently drop targets from the pipeline. The existing `test_stage39_mnase_stress.py` (38 checks) already validates DAG correctness — replacing it with artifact-backed helpers would replace a proven validation with an unproven one.

3. **No test would validate artifact-backed targets better than the existing stress tests already do.** Artifact-backed targets would require rewriting `_mnase_targets()` to iterate `filter_artifacts(..., assay_gate="mnase")` instead of using hardcoded `expand()` calls. The DAG output would be identical — so the only test value is "does the new code produce the same targets as the old code?", which is already validated by the stress tests.

4. **The valuable artifact layer is already built.** `load_artifacts()`, `validate_artifact()`, `filter_artifacts()`, and `artifacts_by_id()` are proven and usable. The contract tests cover all critical cross-references. This infrastructure is ready to use without adopting it into runtime.

## Candidate Next Directions

### A. Output-Contract Generation Adoption

- **Value:** `docs/output-contract.md` could be generated from the artifact inventory, eliminating the `<exp>` placeholder mismatch and ensuring path accuracy.
- **Risk:** Low — read-only generator, no DAG impact.
- **Preconditions:** Explicit decision that generated docs are worth the maintenance overhead. The current manual contract has 62 entries and changes rarely.
- **Blast radius:** `docs/output-contract.md` (one file).
- **Recommended stage if chosen:** 52 (adoption spike, not final generation).
- **Non-goals:** Writing to output-contract.md from CI; replacing manual authoring without review; generating prose sections not in inventory.

### B. Manifest Generation from Artifact Inventory

- **Value:** `scripts/make_manifest.py` could read `manifest_output_type` from the inventory instead of hardcoding `_add_row()` calls.
- **Risk:** Medium — `make_manifest.py` is standalone Python that doesn't import `.smk` files. Adding `workflow/lib/artifact.py` as a dependency introduces a new import path. The manifest stress test (18 checks) would need to be refactored.
- **Preconditions:** `make_manifest.py` would need to import from `workflow/lib/`. Currently it only imports from `scripts/`. This is a cross-directory dependency.
- **Blast radius:** `scripts/make_manifest.py` + `test/test_stage25_manifest_stress.py`.
- **Recommended stage if chosen:** 53 (requires import-path decision first).
- **Non-goals:** Changing manifest output schema; adding new columns; generating from inventory at CI time.

### C. Broader paths.smk Migration Beyond MNase

- **Value:** ChIP-seq, CUT&Tag, and ATAC paths could use helpers like the MNase ones.
- **Risk:** Medium — each assay has unique path conventions. Migrating all at once would touch 3+ `.smk` files and dozens of rule outputs/inputs.
- **Preconditions:** Path helper patterns proven stable for MNase (Stage 44) and verified by contract tests (Stage 47).
- **Blast radius:** `workflow/rules/paths.smk` + `chipseq.smk` + `cuttag.smk` + `atac.smk` + `common.smk` + `targets.smk`.
- **Recommended stage if chosen:** 52-53 (assay-by-assay, not all at once).
- **Non-goals:** Path format changes; single-helper-to-rule-them-all; migration without per-assay contract tests.

### D. Artifact-Assisted Target Helpers

- **Value:** `_mnase_targets()` could iterate `filter_artifacts(..., assay_gate="mnase")` instead of hardcoding 13 `expand()` calls.
- **Risk:** High — `targets.smk` feeds `rule all`. A missing target means missing pipeline output.
- **Preconditions:** All artifact contract tests pass; stress tests verify identical DAG output before and after; explicit decision to accept blast radius.
- **Blast radius:** `targets.smk` + `rule all` + DAG resolution.
- **Recommended stage if chosen:** 54+ (only after manifest/paths adoption proves low-risk).
- **Non-goals:** Replacing `_base_targets()` for non-MNase assays; removing existing stress tests.

### E. Global Target Resolver

- **Value:** A single function that resolves `rule all` inputs from the artifact inventory, eliminating per-assay target builders.
- **Risk:** Very high — `rule all` is the central DAG entry point. A bug here breaks the entire pipeline.
- **Preconditions:** All earlier adoption stages proven stable; artifact inventory complete for all assays; explicit team decision.
- **Blast radius:** `rule all` + all target builders + DAG resolution.
- **Recommended stage if chosen:** 55+ (explicit team decision gate).
- **Non-goals:** Rewriting `rule all` without per-assay target validation; removing stress tests.

### F. Pause Artifact Work — Return to Release Hardening / Science

- **Value:** The pipeline has working peak calling, BigWig generation, MNase fragment analysis, and QC. Scientific features (MNase caller/QC evaluation, CUT&Tag peak caller/spike-in support, histone mark IDR policy) and release hardening (packaging, CI, user docs) would directly benefit users.
- **Risk:** Low — artifact infrastructure is stable and usable. Pausing now does not degrade what exists.
- **Preconditions:** None. The artifact tests continue to run and catch drift.
- **Blast radius:** None — existing tests unchanged.
- **Recommended stage if chosen:** N/A (direction change, not a stage).
- **Non-goals:** Abandoning artifact work permanently; removing artifact tests.

## Recommendation

1. **Pause runtime artifact adoption.** Do not proceed to target helpers (D), global resolver (E), or broad paths.smk migration (C) at this time. The contract test infrastructure (Stages 43-50) is complete and self-sustaining.

2. **If continuing artifact work,** the next defensible step is **Option A — output-contract generation adoption** — but only after explicitly deciding that generating `docs/output-contract.md` from inventory is worth the maintenance trade-off. Stage 50 already proves dry-run equivalence; generation would eliminate the `<exp>` placeholder mismatch at the source.

3. **Preferred next direction: Option F — pause artifact work and return to release hardening and scientific features.** The pipeline would benefit more from:
   - MNase caller/QC evaluation (DANPOS/iNPS/SEM, periodicity, NFR/TSS profiles)
   - CUT&Tag peak caller evaluation (SEACR/GoPeaks) and spike-in normalization
   - Histone mark IDR and pseudoreplicate policy
   - Release packaging, CI smoke tests, and user documentation
   - FRiP/blacklist/library-complexity integration testing

The artifact infrastructure built in Stages 41-50 remains in place. Contract tests continue to validate consistency. If drift or maintenance pain emerges, the adoption path is documented and gated — but the project should not adopt runtime artifactization without demonstrated need.
