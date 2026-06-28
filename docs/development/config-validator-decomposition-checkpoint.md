# Config Validator Decomposition Checkpoint

> State as of PR66 (`a357348 Extract MNase validation helpers (#66)`).

This document is a checkpoint after the config-block extraction phase of the
validator decomposition. It records what has been extracted, what remains in
`validator.py`, the current test safety net, and recommended next steps. No
production code changes are included; this is architecture audit and planning
only.

---

## 1. Completed extraction map

| Module | Lines | Owns |
|---|---|---|
| `encode_pipeline.config.defaults` | 114 | Pure constants, allowed-value sets, static default ranges (assays, layouts, roles, peak modes, MNase ranges, tool key sets, etc.). |
| `encode_pipeline.config.coercion` | 33 | Stateless type coercion: `coerce_int`, `coerce_bool`. |
| `encode_pipeline.config.genome` | 120 | Effective genome size validation, genome resource mapping, Picard/TSS resource gates. |
| `encode_pipeline.config.tools` | 311 | Tool parameter validation for `tool_parameters` block (MACS3, IDR, samtools, etc.). |
| `encode_pipeline.config.reproducibility` | 217 | Reproducibility/IDR/consensus policy validation and warnings. |
| `encode_pipeline.config.qc` | 45 | QC config block normalization (unknown keys silently ignored). |
| `encode_pipeline.config.cuttag` | 115 | CUT&Tag config validation (`peak_caller`, `seacr` block, threshold). |
| `encode_pipeline.config.mnase` | 173 | MNase config validation (`mono_range`, `fragments`, `dyad_range`, `callers`, range helpers). |

All extracted modules depend only on `defaults` and/or `coercion`. None import
`validator.py` or reference `ValidationError` directly; they receive an
`error_cls` parameter.

---

## 2. Current `validator.py` remaining responsibilities

`src/encode_pipeline/config/validator.py` is now **1,217 lines** (down from
~2,070 at the start of the decomposition). It still contains:

### Public API and orchestration
- `ValidationError` — must stay in this module for backward compatibility.
- `validate_config()` — top-level orchestrator (~213 lines). Still owns:
  - samples path validation
  - top-level scalar coercion (`threads`, `mapq`, `binsize`)
  - `remove_dup`, `trim`, `extend_reads`, `use_control`, `multiqc`, `stage4b`,
    `stage5` normalization
  - stage-gate checks (`stage5` requires `stage4b`, IDR modes require `stage4b`)
  - delegation to extracted config-block helpers
- `load_and_validate_samples()` — TSV loading and per-row/cross-reference
  validation (~243 lines).
- `validate_replicate_groups()` — cross-replicate consistency and IDR
  eligibility (~241 lines).
- `validate_picard_reference_resources()` / `validate_tss_annotation_resources()`
  — thin wrappers around `config.genome`.
- `main()` — CLI entry point for `scripts/validate_samples.py`.

### Private helpers still in `validator.py`
- `_check_fastq_exists()` / `_check_bowtie2_index()` / `_validate_strict_inputs()`
  — strict input validation (Stage 30).
- `_coerce_int()` — thin wrapper around `coercion.coerce_int`.
- `_validate_effective_genome_size()` — thin wrapper around `genome` module.
- `_sanitize_identifier()` / `_parse_positive_int()` — sample-sheet helpers.
- `_load_yaml()` / `_parse_config_minimal()` — minimal YAML parsing.
- `_validate_genome_resources()` / `_validate_tool_params()` /
  `_validate_idr_settings()` / `_validate_reproducibility()` /
  `_validate_cuttag_config()` / `_validate_mnase_config()` /
  `_validate_qc_config()` — thin wrappers around extracted modules.

### Thin wrappers vs substantial logic
- **Thin wrappers:** all config-block validators, `_coerce_int`,
  `_validate_effective_genome_size`, Picard/TSS resource wrappers.
- **Still substantial:** `validate_config`, `load_and_validate_samples`,
  `validate_replicate_groups`, `_load_yaml` / `_parse_config_minimal`,
  `_validate_strict_inputs`, `_sanitize_identifier`, `_parse_positive_int`,
  `main()`.

---

## 3. Current test safety net

### Config-block characterization tests
- `test/config/test_coercion.py`
- `test/config/test_genome_resources.py`
- `test/config/test_tool_parameters.py`
- `test/config/test_reproducibility.py`
- `test/config/test_qc.py`
- `test/config/test_cuttag.py`
- `test/config/test_mnase.py`

### Helper unit tests for extracted modules
- `test/config/test_coercion_helpers.py`
- `test/config/test_genome_helpers.py`
- `test/config/test_tool_helpers.py`
- `test/config/test_reproducibility_helpers.py`
- `test/config/test_qc_helpers.py`
- `test/config/test_cuttag_helpers.py`
- `test/config/test_mnase_helpers.py`

### Integration / smoke / behavior guards
- `test/config/test_validation.py` — primary characterization suite.
- `test/test_stage8_smoke_profiles.py` — dry-run smoke profiles.
- `scripts/validate_samples.py --config config/config.yaml` — happy-path CLI.
- `test/test_no_hardcoded_paths.py` — environment-safety guard.
- `test/check_snakemake_lint.py` — workflow lint baseline.

---

## 4. Recommended next phase

The config-block decomposition is complete. The remaining large blocks in
`validator.py` are sample-sheet validation and the CLI/orchestrator surface.
Options for the next pin/extract pairs:

| Option | Risk | Notes |
|---|---|---|
| **A. Sample loading / row validation** | Medium | `load_and_validate_samples` is ~243 lines with clear internal passes (TSV parsing, per-row validation, control cross-references, replicate groups, strict inputs). Splitting into `samples.loader`, `samples.row`, `samples.controls`, `samples.replicates` matches the long-term architecture. Risk is medium because sample-sheet behavior is user-facing and subtle (control references, role gates, experiment/condition defaulting). |
| **B. Strict input checks** | Low | `_validate_strict_inputs` and its file-existence helpers are isolated and I/O-bound, but small (~58 lines). Could be extracted quickly, though it may make more sense as part of the sample-sheet extraction rather than a standalone PR. |
| **C. Replicate group validation** | Medium-High | `validate_replicate_groups` is ~241 lines and entangled with IDR eligibility logic for Stage 5 / Stage 55 / Stage 64 / Stage 65. It needs careful characterization before extraction. |
| **D. Minimal YAML parsing** | Medium | `_load_yaml` / `_parse_config_minimal` are ~165 lines and only needed by the CLI. Extracting to `encode_pipeline.cli.validate` or a `config.yaml_loader` module is reasonable, but it is CLI-adjacent and should wait until sample-sheet work is done. |
| **E. Validator orchestrator thinning** | Low-Medium | Once sample loading and YAML parsing are extracted, `validate_config` can be further thinned, but the top-level orchestrator should remain in `validator.py` for backward compatibility. |

### Suggested pin/extract PR pairs

- **PR68:** Pin sample loading / row normalization behavior
  - `load_and_validate_samples` Pass 1 and Pass 2 behavior
  - required column checks
  - per-row field validation and defaulting
  - control cross-reference behavior
  - duplicate sample ID detection
- **PR69:** Extract sample loading / row validation helpers
  - Target: `encode_pipeline.samples.loader`, `encode_pipeline.samples.row`,
    `encode_pipeline.samples.controls`
  - Keep `load_and_validate_samples` as a thin orchestrator in `validator.py`

- **PR70:** Pin replicate group behavior
  - Consistency checks, duplicate `(bio_rep, tech_rep)` detection
  - Control consistency rules
  - Stage 5 / Stage 55 / Stage 64 / Stage 65 IDR eligibility
- **PR71:** Extract replicate group helpers
  - Target: `encode_pipeline.samples.replicates`

- **PR72+:** Strict input checks and YAML parsing can follow as smaller,
  lower-risk cleanups after the sample-sheet phase.

---

## 5. Recommended immediate next step

**Proceed with PR68/PR69 as the preferred next step** unless inspection reveals
a blocker.

Rationale:
- `load_and_validate_samples` is the largest remaining block and has clear
  internal boundaries.
- Per-row validation is mostly stateless and easy to characterize.
- Control cross-reference logic is localized and testable.
- Replicate group validation is more entangled with IDR eligibility and should
  come after the row/control extraction stabilizes.

If PR68 discovers that per-row and control logic cannot be cleanly separated,
the fallback is to pin them together and extract a single larger
`encode_pipeline.samples` submodule first, then split later.

---

## 6. Explicit non-goals

- No FastAPI or external service integration.
- No changes to validation behavior, error messages, warning text, or defaults.
- No schema, Snakefile, rule, CI, lockfile, or CLI behavior changes.
- No migration of legacy `test/test_stage53*.py` / `test/test_stage55*.py`
  scripts in this phase.
- No scientific default or threshold changes.
- No new `encode_pipeline.config.cli` module; the CLI stays in
  `encode_pipeline.cli`.

---

## 7. Dependency direction (unchanged)

- Extracted submodules must not import `validator.py`.
- `validator.py` may import submodules.
- Public APIs (`validate_config`, `load_and_validate_samples`, etc.) remain in
  `validator.py` as thin orchestrators or re-exports.
- Tests should continue to prefer public APIs over private helpers.
