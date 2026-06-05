# Adding a New Output Type or Assay — Developer Checklist

Use this checklist when adding a new output type (e.g., a new QC metric,
BigWig variant, or fragment class) or a new assay to the pipeline. Each
numbered item names the specific file(s) to modify and what to change.

---

## A. Adding a new output type

### 1. Sample sheet schema / validation

- [ ] `workflow/schemas/samples.schema.yaml`: If the output requires a new
  sample sheet column, add it to the schema.
- [ ] `scripts/validate_samples.py`: Add validation for any new column value
  or config key.

### 2. Config schema

- [ ] `workflow/schemas/config.schema.yaml`: If the output is gated by a new
  config key, document the key, its type, default, and allowed values.
- [ ] `scripts/validate_samples.py`: Add `_validate_*_config()` logic for the
  new key. Gate-checked keys should produce clear error messages.

### 3. Workflow rule file

- [ ] Add the Snakemake rule in the appropriate `workflow/rules/*.smk` file.
- [ ] Use `:q` quoting on all shell path expansions.
- [ ] Use the smallest appropriate Conda environment (`../envs/<name>.yml`).
- [ ] Gate the rule's inputs/outputs on assay or config presence as appropriate
  (e.g., `_is_mnase(wc)` for MNase-only rules).

### 4. Target helper wiring

- [ ] `workflow/Snakefile`: Add the new output path to the correct
  `_*_targets()` function (e.g., `_mnase_targets()` for MNase outputs,
  `_single_sample_qc_targets()` for QC outputs).
- [ ] Gate on the correct sample ID list (`MNASE_SAMPLE_IDS`,
  `PEAK_SAMPLE_IDS`, `TREATMENT_SAMPLE_IDS`, etc.).
- [ ] If the output is config-gated (e.g., `qc.signal_tracks: true`),
  add the gating condition.

### 5. `pipeline_done` wiring

- [ ] `workflow/rules/report.smk`: Add the new output as a named input to
  the `pipeline_done` rule. Use lambda gating (e.g., `_is_mnase(wc)` or
  `QC_CONFIG.get(...)`) to return `[]` when the output is not applicable.
- [ ] Verify that `snakemake <sample>/logs/<sample>.pipeline.done -n` resolves
  the new output.

### 6. Manifest rows

- [ ] `scripts/make_manifest.py`: Add `_add_row(...)` calls in the appropriate
  builder function (`_build_sample_rows` or `_build_experiment_rows`).
- [ ] Use the existing `output_type` vocabulary (snake_case, tool-derived names).
- [ ] Gate on the same conditions as the Snakemake rule (assay, config, stage).

### 7. MultiQC custom content

- [ ] `workflow/multiqc_config.yaml`: Add a `custom_data.<name>` section and a
  corresponding `sp.<name>` entry with a `fn:` glob pattern.
- [ ] Write a `headers:` block for each column you want visible in the report.
- [ ] Verify that existing custom content sections are not modified.

### 8. Output contract update

- [ ] `docs/output-contract.md`: Add a row to the appropriate table
  (single-sample, experiment-level, IDR, or project-level).
- [ ] Include `output_type`, `method`, `rule`, `path`, and `status` columns.
- [ ] Add any new gating condition to the Gating Conditions table.

### 9. Smoke / dry-run tests

- [ ] `test/test_stage8_smoke_profiles.py`: If the output is tied to a specific
  assay or config, add or update a test profile.
- [ ] `test/test_stage39_mnase_stress.py` (or equivalent stage test): Add
  assertions that the new rule name appears in dry-run job output.
- [ ] Verify that `pipeline.done` direct-target dry-run schedules the new rule.

### 10. Negative tests for assay gating

- [ ] Add assertions that the new output is NOT scheduled for assays where it
  should not apply (e.g., peak-centric outputs must not appear in MNase
  dry-run output).
- [ ] `test/test_validation_stress.py`: Add config rejection cases for invalid
  values of any new config keys.

### 11. Documentation

- [ ] `docs/assay-policy.md`: Update the relevant assay section with the new
  output, its config keys, and any constraints.
- [ ] `docs/configuration.md`: Document any new config keys with examples.
- [ ] `README.md`: Update the Features and/or Limitations sections if the new
  output changes the pipeline's capability surface.
- [ ] `CHANGELOG.md`: Add an entry under `[Unreleased]` in the appropriate
  category (`### Added`, `### Changed`, `### Fixed`).

---

## B. Adding a new assay

All of the above (items 1-11), plus:

- [ ] **Sample sheet:** `scripts/validate_samples.py` — add assay to allowlist,
  define layout constraints (PE-only, etc.), define `peak_mode` constraints.
- [ ] **Sample sheet schema:** `workflow/schemas/samples.schema.yaml` — add
  assay to the `assay` enum.
- [ ] **Config:** Add assay-specific config block if needed.
- [ ] **Policy file:** Create `workflow/rules/<assay>.smk` with
  `get_remove_dup_<assay>()`, `get_extend_reads_<assay>()`,
  `get_macs3_args_<assay>()` policy functions. For non-peak-centric assays,
  the MACS3 function should raise `ValueError`.
- [ ] **Dispatch functions:** `workflow/Snakefile` — add `<assay>` branch to
  `get_remove_dup()`, `get_macs3_args()`, `get_extend_reads()`.
- [ ] **Derived lists:** `workflow/Snakefile` — add `<ASSAY>_SAMPLE_IDS` and
  `<ASSAY>_MULTI_BIOREP_EXPERIMENTS` derived lists.
- [ ] **Target helper:** Add `_<assay>_targets()` function.
- [ ] **Rule file include:** Add `include: "rules/<assay>.smk"` to Snakefile.
- [ ] **Manifest:** Add assay-specific manifest rows with `not_applicable`
  status for outputs that don't apply (e.g., peak counts for non-peak-centric
  assays).
- [ ] **Test profile:** `test/profiles/<assay>_<layout>_noctrl/` with
  `config.yaml` and `samples.tsv`.
