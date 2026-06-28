# Config Validator Decomposition

`src/encode_pipeline/config/validator.py` is the legacy sample-sheet and
config validator. It is behavior-critical and currently ~2,070 lines. This
document records its responsibilities, the contracts that must not change,
proposed extraction boundaries, and the recommended first extraction for PR52.

## Current responsibilities

The validator module handles all of the following:

1. **Config loading and merging**
   - Minimal YAML parsing (`_load_yaml`, `_parse_config_minimal`).
   - Default value injection for missing config keys.

2. **Top-level config validation**
   - `validate_config`: orchestrates all config-level checks.
   - Thread, mapq, and numeric coercion (`_coerce_int`, `_parse_positive_int`).
   - Boolean coercion (`_coerce_bool`).
   - `remove_dup` allowed values.

3. **Assay-specific config validation**
   - CUT&Tag-specific checks (`_validate_cuttag_config`).
   - MNase-specific checks (`_validate_mnase_config`, ranges, callers, dyad).

4. **Genome resources validation**
   - Effective genome size checks (`_validate_effective_genome_size`).
   - Chrom sizes / Picard / TSS annotation resource checks
     (`validate_picard_reference_resources`, `validate_tss_annotation_resources`,
     `_validate_genome_resources`).

5. **Tool parameter validation**
   - MACS3, IDR, and other tool-specific settings (`_validate_tool_params`).

6. **Reproducibility policy validation**
   - IDR settings (`_validate_idr_settings`).
   - Consensus and reproducibility policy (`_validate_reproducibility`).

7. **Sample sheet loading and validation**
   - TSV parsing and row normalization (`load_and_validate_samples`).
   - Column presence and type coercion.
   - Sample ID sanitization (`_sanitize_identifier`).
   - Layout, assay, target, peak_mode, genome validation.
   - Control sample referencing and role checks.
   - Replicate group validation (`validate_replicate_groups`).

8. **Strict input validation**
   - FASTQ and Bowtie2 index existence checks
     (`_validate_strict_inputs`, `_check_fastq_exists`, `_check_bowtie2_index`).

9. **CLI and reporting**
   - `main()` for `scripts/validate_samples.py`.
   - Structured output and exit codes used by CI and the Snakefile.

## Behavior contracts that must not change

- `validate_config` and `load_and_validate_samples` must remain importable
  from `encode_pipeline.config.validator` for backward compatibility.
- `ValidationError` must remain the exception type raised for invalid input.
- Error and warning messages must remain stable unless a release intentionally
  changes user-facing text.
- Default values produced by `validate_config` must match the current behavior
  exactly; downstream rules depend on them.
- Exit codes from the CLI `main()` must remain stable.
- `encode_pipeline.config.validate`, `encode_pipeline.samples.load`,
  `encode_pipeline.samples.validate`, and `encode_pipeline.cli.validate` must
  continue to re-export or wrap the legacy implementation without changing
  their public signatures.

## Proposed module boundaries

The long-term target architecture keeps `validator.py` as a thin orchestrator
and moves specialized logic into focused submodules:

| Module | Responsibility |
|---|---|
| `encode_pipeline.config.defaults` | Pure constants, allowed-value sets, and default numeric ranges. |
| `encode_pipeline.config.coercion` | Stateless type coercion helpers (`_coerce_int`, `_coerce_bool`, etc.). |
| `encode_pipeline.config.genome` | Effective genome size and genome resource validation. |
| `encode_pipeline.config.reproducibility` | IDR, consensus, and reproducibility policy validation. |
| `encode_pipeline.config.tools` | Tool-specific parameter validation (MACS3, IDR, etc.). |
| `encode_pipeline.config.mnase` | MNase-specific config validation. |
| `encode_pipeline.config.cuttag` | CUT&Tag-specific config validation. |
| `encode_pipeline.samples.loader` | TSV parsing and row normalization. |
| `encode_pipeline.samples.row` | Per-sample row validation and sanitization. |
| `encode_pipeline.samples.controls` | Control sample referencing and role checks. |
| `encode_pipeline.samples.replicates` | Replicate group validation. |
| `encode_pipeline.cli.validate` | CLI entry point and output formatting; remains outside config modules. |

The config decomposition must not introduce an `encode_pipeline.config.cli`
module. The CLI stays in the existing `encode_pipeline.cli` package.

## Dependency direction

- Submodules must not import `validator.py`.
- `validator.py` may import submodules.
- `encode_pipeline.config.validate`, `encode_pipeline.samples.load`, etc.,
  continue to expose the public API; their implementation can switch from
  re-exporting `validator.py` to importing submodules without changing signatures.
- Tests should prefer the public API (`encode_pipeline.config.validate.validate_config`,
  `encode_pipeline.samples.load.load_and_validate_samples`) over private helpers.

## Characterization test matrix

Before any extraction, the following existing test coverage characterizes the
behavior that must be preserved:

- `test/config/test_validation.py` — 41 tests covering config defaults,
  sample-sheet structure, controls, and MNase validation.
- `test/test_stage8_smoke_profiles.py` — dry-run smoke profiles that exercise
  default config and sample validation through the CLI.
- `scripts/validate_samples.py --config config/config.yaml` — happy-path CLI.
- Legacy config validation scripts in `test/test_stage53*.py` and
  `test/test_stage55*.py` are classified as `migrate-to-pytest`. They should be
  reviewed for missing assertions and migrated into native pytest if needed;
  they are not run by default.

Any PR that extracts logic must run the active native tests and CLI checks
above; legacy quarantined scripts should be reviewed for coverage gaps, not run
by default.

## Recommended PR52 first extraction

**Target:** Extract pure constants/defaults/allowed-values into
`src/encode_pipeline/config/defaults.py`.

**Why this is the safest first extraction:**

1. **Low behavior risk.** Constants have no control flow and produce the same
   values before and after extraction.
2. **Easy to characterize.** Existing tests already assert the effects of
   defaults through `validate_config` output; no new heavy tests are needed.
3. **Minimal dependency direction.** A `defaults` module has no imports from
   `validator.py`; `validator.py` simply reads from it.
4. **Public API unchanged.** `validate_config` continues to live in
   `validator.py`; only internal references move.
5. **Small review surface.** A constants module is easy to diff and reason about.

**Acceptable content for `defaults.py`:**

- Allowed assay names, layout values, and role values.
- Default numeric ranges that are static (e.g., MNase range defaults).
- Known validation allowed-value sets (e.g., `remove_dup` values).
- Static warning/default constants that are not entangled with control flow.

**Not acceptable for PR52:**

- `validate_config` control flow.
- Sample loading or row validation.
- Warning emission behavior.
- Error message formatting.
- MNase validation logic beyond static constants.
- Reproducibility/IDR/consensus validation logic.
- Any change to public imports or CLI behavior.

## PR52 fallback

If closer inspection shows that candidate constants are too intertwined with
control flow, conditional defaults, or mutable state, PR52 must **stop** and add
characterization tests instead of forcing an extraction. The decomposition spec
stays valid; the sequence simply adds more tests before the next extraction.

## Non-goals

- Do not rewrite `validate_config` in PR52.
- Do not change validation behavior, error messages, warning text, defaults,
  schema files, Snakemake rules, IDR thresholds, target lists, or config
  structure.
- Do not refactor the CLI entry point.
- Do not migrate legacy tests in PR52 unless their coverage is needed to
  characterize the extracted constants.

## Recommended PR53: pin coercion behavior

**Target:** Add direct-API pytest characterization tests for
`_coerce_int`-like behavior in `validate_config` (`threads`, `mapq`,
`binsize`) and boolean coercion patterns (`trim`, `use_control`, `multiqc`,
`stage4b`, `stage5`) as well as `_coerce_bool` reproducibility call sites
(`reproducibility.enabled`, `reproducibility.consensus.enabled`, and
`reproducibility.idr.*` flags).

**Why:** Extraction PR54 will move `_coerce_int` and `_coerce_bool` into
`encode_pipeline.config.coercion`. Without pinned behavior, subtle
normalization changes could go unnoticed.

**Acceptable content:** Tests that call
`encode_pipeline.config.validate.validate_config` directly and assert
normalized return values or `ValidationError` substrings. No changes to
`validator.py` or `defaults.py`.

## Recommended PR54: extract coercion helpers

**Target:** Move `_coerce_int` and `_coerce_bool` into
`encode_pipeline.config.coercion` and update `validate_config` (and any
reproducibility helpers) to import them, preserving the behavior pinned by
PR53.

**Why:** This is the next-lowest-risk extraction after constants. The helpers
are stateless, have no I/O, and the PR53 tests provide a safety net.

## Recommended PR55: pin genome resource validation

**Target:** Add direct-API pytest characterization tests for
`genome_resources`, `effective_genome_size`, optional genome resource path
fields, and the Picard/TSS resource gates.

**Why:** Extraction PR56 will move effective genome size checks, genome resource
mapping validation, and `validate_picard_reference_resources` /
`validate_tss_annotation_resources` into `encode_pipeline.config.genome`.
These checks combine config normalization, filesystem path validation, and
treatment-genome-only resource gates, so they need a focused safety net before
the module move.

**Acceptable content:** Tests that call
`encode_pipeline.config.validate.validate_config`,
`validate_picard_reference_resources`, and
`validate_tss_annotation_resources` directly. No changes to
`validator.py`, `defaults.py`, schema files, CLI behavior, or Snakemake rules.

## Recommended PR56: extract genome resource helpers

**Target:** Move effective genome size validation, genome resource mapping
validation, and Picard/TSS resource gate functions into
`encode_pipeline.config.genome`, preserving the behavior pinned by PR55.

**Why:** This is the next focused validator boundary after coercion. The logic
is still small enough to review, but more stateful than pure coercion because
it validates filesystem-backed resource paths.

## Recommended PR57: pin tool parameter validation

**Target:** Add direct-API pytest characterization tests for
`tool_parameters` structure validation, unknown tool/key rejection,
per-tool default expansion, integer/float/bool/string normalization, samtools
filter flag parsing, and unsupported value rejection.

**Why:** Extraction PR58 will move `_validate_tool_params` and its nested
normalization helpers into `encode_pipeline.config.tools`. The current logic
has several tool-specific branches and mixed normalization rules, so focused
tests should pin behavior before the move.

**Acceptable content:** Tests that call
`encode_pipeline.config.validate.validate_config` directly and assert
normalized `tool_parameters` values or `ValidationError` substrings. No changes
to `validator.py`, `defaults.py`, schema files, CLI behavior, or Snakemake
rules.

## Recommended PR58: extract tool parameter helpers

**Target:** Move `_validate_tool_params` and its nested helper logic into
`encode_pipeline.config.tools`, preserving the behavior pinned by PR57.

**Why:** Tool parameter validation is a coherent boundary after genome
resources. It depends on static defaults and pure normalization logic, but it
is broad enough to warrant a dedicated characterization PR first.

## Recommended PR59: pin reproducibility and IDR behavior

**Target:** Add direct-API characterization tests for
`reproducibility` structure validation, consensus numeric settings, IDR
settings, Stage 4b gates, and reproducibility warnings.

**Why:** Extraction PR60 will move `_validate_reproducibility` and
`_validate_idr_settings` into `encode_pipeline.config.reproducibility`. These
helpers contain subtle behavior that must stay unchanged, including disabled
reproducibility short-circuiting nested validation, IDR settings only being
validated when Stage 5 or a reproducibility IDR mode is enabled, and warnings
for experimental broad IDR and contradictory legacy chipseq narrow IDR config.

**Acceptable content:** Tests that call
`encode_pipeline.config.validate.validate_config` directly and assert
normalized `reproducibility` / `idr` values, `ValidationError` substrings, and
warning messages. No changes to `validator.py`, `defaults.py`, schema files,
CLI behavior, Snakemake rules, lockfiles, or CI.

## Recommended PR60: extract reproducibility helpers

**Target:** Move `_validate_reproducibility` and `_validate_idr_settings` into
`encode_pipeline.config.reproducibility`, preserving the behavior pinned by
PR59.

**Why:** Reproducibility and IDR policy is a coherent validator boundary after
tool-parameter extraction. It depends on static defaults, coercion helpers, and
an injected error class, but it should not import `validator.py` or introduce
new public APIs.

## Recommended PR61: pin QC config behavior

**Target:** Add direct-API pytest characterization tests for
`_validate_qc_config` behavior, including default expansion, boolean
normalization, invalid-value rejection, mapping validation, and silent
ignoring of unknown QC keys.

**Why:** Extraction PR62 will move `_validate_qc_config` into
`encode_pipeline.config.qc`. Without pinned behavior, default values and
unknown-key handling could drift.

**Acceptable content:** Tests that call
`encode_pipeline.config.validate.validate_config` directly and assert
normalized `validated["qc"]` values or `ValidationError` substrings. No
changes to `validator.py`, `defaults.py`, schema files, CLI behavior, or
Snakemake rules.

## Recommended PR62: extract QC helpers

**Target:** Move `_validate_qc_config` and any QC-related constants into
`encode_pipeline.config.qc`, preserving the behavior pinned by PR61.

**Why:** QC validation is a coherent boundary after reproducibility
extraction. It is mostly stateless boolean normalization, but the
unknown-key behavior and default mapping must stay stable.

## Recommended PR63: pin CUT&Tag config behavior

**Target:** Add direct-API pytest characterization tests for
`_validate_cuttag_config` behavior, including default expansion,
`peak_caller` gating, `seacr` mapping validation, unknown-key rejection for
both top-level and nested keys, boolean/string-boolean normalization for
`seacr.enabled`, allowed-value checks for `seacr.mode` and
`seacr.normalization`, and open-interval `(0, 1)` validation for
`seacr.threshold` including explicit bool rejection.

**Why:** Extraction PR64 will move `_validate_cuttag_config` into
`encode_pipeline.config.cuttag`. CUT&Tag validation differs from QC in that
unknown keys are rejected rather than silently ignored, and it contains
float-threshold and enum checks that must stay stable.

**Acceptable content:** Tests that call
`encode_pipeline.config.validate.validate_config` directly and assert
normalized `validated["cuttag"]` values or `ValidationError` substrings.
Prefer stable substrings (e.g. `cuttag: unknown key`,
`cuttag.seacr: unknown key`) over exact key lists so tests do not become
fragile when the known-key set changes. No changes to `validator.py`,
`defaults.py`, schema files, CLI behavior, or Snakemake rules.

## Recommended PR64: extract CUT&Tag helpers

**Target:** Move `_validate_cuttag_config` and any CUT&Tag-related constants
into `encode_pipeline.config.cuttag`, preserving the behavior pinned by
PR63.

**Why:** CUT&Tag validation is the next focused boundary after QC. It is
stateless, has no I/O, and the PR63 tests provide a safety net for the
move.

## Files involved

- `src/encode_pipeline/config/validator.py` — current legacy module.
- `src/encode_pipeline/config/validate.py` — public config validation API.
- `src/encode_pipeline/samples/load.py` — public sample loading API.
- `src/encode_pipeline/samples/validate.py` — row-level validation wrapper.
- `src/encode_pipeline/cli/validate.py` — CLI entry point.
- `scripts/validate_samples.py` — backward-compatible CLI wrapper.
- `test/config/test_validation.py` — primary characterization test suite.
- `test/test_stage8_smoke_profiles.py` — smoke profile characterization.
