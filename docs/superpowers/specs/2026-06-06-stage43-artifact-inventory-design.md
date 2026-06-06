# Stage 43: Artifact Inventory — Design Spec

**Date:** 2026-06-06
**Status:** implemented
**Scope:** Create machine-readable YAML inventory of all pipeline output types; docs-only, no DAG changes
**Excluded:** Artifact dataclass, artifact_path(), AssayPolicy YAML, paths.smk, code generation, manifest generation

---

## 1. Goals

1. Create a single machine-readable YAML file listing every stable output type in the pipeline (62 entries).
2. Each entry records: id, description, scope, level, assay gate, path template (with `results/` prefix), producing rule, tool, manifest_output_type (string or null), and gating conditions.
3. Provide a hand-maintained reference vocabulary for future stages (paths.smk, artifact dataclass).
4. Add a PyYAML-based consistency test with 10 checks: YAML parses, IDs unique, no workstation paths, no future feature keywords, manifest coverage (ast extraction), manifest extras, OC cross-reference (both directions, section-based parsing), schema completeness (exactly 13 fields, booleans verified).

---

## 2. Inventory File Path

**Decision: `docs/architecture/artifact-inventory.yaml`**

Rationale:
- The artifact roadmap (`docs/architecture/artifact-roadmap.md`) defines Stage 43 as a documentation/machine-readable contract step. Placing the inventory alongside the roadmap keeps architecture documents together.
- `docs/artifact-inventory.yaml` would suggest a top-level operational file; the inventory is a reference document, not a runtime config.
- Future `paths.smk` and artifact dataclass stages will reference this file, and `docs/architecture/` is the natural source-of-truth location for architectural contracts.

---

## 3. Inventory Schema

Each artifact entry has these fields:

| Field | Type | Description |
|:---|:---|:---|
| `id` | string (required) | Stable identifier matching `output_type` in output-contract.md |
| `description` | string | Human-readable one-liner |
| `scope` | enum | `sample`, `experiment`, `project`, `reference` |
| `level` | enum | `per_sample`, `pooled_experiment`, `project`, `reference` |
| `assay_gate` | string | Which samples produce this: `all`, `peak_centric`, `mnase`, `cuttag`, `chipseq`, `atac`, `idr` |
| `path_template` | string | Relative path using output-contract style with default `results/` prefix and `<sample>`/`<experiment>`/`<genome>` placeholders |
| `producing_rule` | string | Snakemake rule name or `(common.smk)` / `(peaks.smk)` for shared rules |
| `tool` | string | Primary tool used (e.g., `bowtie2+samtools`, `macs3 callpeak`, `bamCoverage --MNase`) |
| `manifest_output_type` | string or null | Exact string used by `make_manifest.py`, or `null` if not manifest-tracked |
| `pipeline_done` | boolean | Whether this output is a named input to `pipeline_done` |
| `rule_all` | boolean | Whether this output is expanded in a `_*_targets()` function called by `rule all` |
| `config_gate` | string or null | Config condition that gates this output (e.g., `qc.signal_tracks: true`, `stage5: true`), or `null` if always produced |
| `notes` | string | Optional: deprecation status, stage introduced, deferred features |

Example entry:

```yaml
- id: mnase_mono_bam
  description: Mono-nucleosome fragment BAM
  scope: sample
  level: per_sample
  assay_gate: mnase
  path_template: results/<sample>/03_fragments/<sample>.mono.bam
  producing_rule: mnase_split_mono
  tool: alignmentSieve
  manifest_output_type: mnase_mono_bam
  pipeline_done: true
  rule_all: true
  config_gate: null
  notes: Stage 39
```

---

## 4. Initial Coverage

**Decision: Option A — all stable output-contract.md rows (62 output types).**

Scope includes:
- 36 sample-level outputs
- 21 experiment-level outputs
- 4 project-level outputs
- 1 reference output

62 total entries. BAI companion files (`final_bai`, `mnase_sub_bai`, etc.) are included. Stage 43 also fixed two documentation-only output-contract omissions: `pooled_final_bai` and `biorep_final_bai`.

Excluded:
- Internal temporary files (sorted bedGraph intermediates, `.tmp_*` files)
- Snakemake internal markers (`.snakemake/` timestamps)
- Legacy script outputs not produced by the Snakemake workflow

---

## 5. Consistency Test

**New test: `test/test_stage43_artifact_inventory.py`**

Design:
- PyYAML is required (`import yaml`). Import/load failure is FAIL.
- Manifest cross-reference uses Python `ast` to inspect `_add_row()` calls in `make_manifest.py`, extracting the `output_type` argument including f-string `JoinedStr` values.
- Output-contract cross-reference parses only the section between `## Current Output Types` and `## Manifest Field Schema`, extracting first-column table entries.
- 10 checks:
  1. YAML parses successfully; import/load failures fail the test
  2. Inventory contains entries
  3. All `id` values are unique
  4. No `path_template` contains absolute workstation paths
  5. No entry mentions future-feature keywords as implemented
  6. Every manifest `output_type` has a matching inventory entry
  7. Every non-null inventory `manifest_output_type` is emitted by the manifest
  8. All output-contract types are present in the inventory
  9. All inventory IDs are present in the output-contract (bidirectional exact match)
  10. Schema completeness: exactly 13 required fields per entry, `pipeline_done` and `rule_all` are booleans

---

## 6. Non-Goals (repeated for emphasis)

- No DAG behavior changes
- No Snakemake rule changes
- No target helper changes
- No manifest generation from inventory
- No `paths.smk`
- No Artifact dataclass
- No `artifact_path()`
- No AssayPolicy YAML
- No global target resolver

---

## 7. Files

| File | Action |
|:---|:---|
| `docs/architecture/artifact-inventory.yaml` | **Create** — 62 entries, ~300 lines |
| `test/test_stage43_artifact_inventory.py` | **Create** — 10 consistency checks |
| `docs/output-contract.md` | Modify — added `pooled_final_bai` and `biorep_final_bai` rows (docs-only contract fix) |
| `test/test_no_hardcoded_paths.py` | Modify — added Stage 43 test to skip list |
| `docs/superpowers/specs/2026-06-06-stage43-artifact-inventory-design.md` | This file |
| `docs/superpowers/plans/2026-06-06-stage43-artifact-inventory.md` | Plan file |

No modifications to `workflow/Snakefile`, workflow rules `.smk` files, `scripts/make_manifest.py`, `scripts/validate_samples.py`, or any runtime config.

---

## 8. Risks

| Risk | Mitigation |
|:---|:---|
| Inventory drifts from code | Test cross-references `make_manifest.py` output_type strings; drift detected at CI time |
| YAML hard to hand-maintain | 62 entries; comparable in scope to output-contract.md |
| PyYAML not available in test env | Test fails, because PyYAML is part of the project test/runtime environment |
