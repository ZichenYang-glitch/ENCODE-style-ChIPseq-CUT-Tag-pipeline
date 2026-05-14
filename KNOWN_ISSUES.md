# Known Issues and Follow-ups

This file records non-blocking issues left after the Stage 1 Snakemake rule
modularization. Stage 1 has passed dry-run and rule-list validation, but these
items should be handled before calling the workflow production-hardened.

## Stage 1 Status

- The modular Snakemake DAG is in place under `workflow/rules/`.
- `workflow/Snakefile` is now the orchestration entry point.
- ChIP-seq and CUT&Tag policy functions are separated in assay-specific rule files.
- Control samples can be modeled as first-class sample rows through `role` and
  `control_sample`.
- Dry-run validation passes for the current test configuration.

## High Priority

1. Document the new sample model.
   - Add `role` and `control_sample` examples to `README.md`.
   - Add an example treatment/control pair to `config/samples.tsv`.
   - Explain the difference between `control_sample` and external `control_bam`.

2. Add schema/sample validation.
   - Add schemas for `config/config.yaml` and `config/samples.tsv`.
   - Validate optional columns: `role`, `control_sample`, and `control_bam`.
   - Keep backward compatibility with the current sample sheet.

3. Run a real end-to-end smoke test.
   - Dry-run checks DAG correctness, not tool behavior.
   - Test with real Bowtie2 index paths and small FASTQs.
   - Include at least one treatment sample and one FASTQ-based control sample.

## Medium Priority

4. Restore optional `plotFingerprint` QC.
   - The legacy `scripts/chipseq.sh` runs `plotFingerprint` when available.
   - The modular Snakemake workflow does not yet expose this as a rule.

5. Clean exported Conda environment metadata.
   - `chipseq.yml` may contain a machine-specific `prefix`.
   - Remove local prefixes before sharing or publishing the environment file.

6. Harden the legacy single-sample script.
   - `scripts/chipseq.sh` still uses input-derived Trim Galore output discovery
     with `ls ... | tail -n 1`.
   - This can select stale files if an output directory is reused.
   - The Snakemake path normalization avoids this issue, but the legacy script
     remains available for compatibility.

7. Reduce silent error handling.
   - Some optional QC steps intentionally continue on failure.
   - This is convenient for exploratory runs, but can hide missing QC in stricter
     production contexts.

## Low Priority

8. Split environments by responsibility.
   - Stage 1 keeps using `workflow/envs/chipseq.yml`.
   - Future stages can split this into core, ChIP-seq, CUT&Tag, and reporting
     environments if dependency resolution becomes slow or fragile.

9. Add richer reporting.
   - A later report step can summarize alignment rate, duplicate rate, peak count,
     and control usage across all treatment samples.

10. Decide how much control output to publish.
    - Stage 1 produces control bigWigs because they are useful for QC.
    - Future config may make control bigWig/report publication optional.

## Notes

- `control_bam` is checked only when `use_control: true`; when `use_control: false`,
  it remains ignored for backward compatibility.
- Control rows are preprocessed but do not run peak calling.
- SE ChIP-seq treatment samples may depend on MACS3 peak calling before bigWig
  generation so the fragment-size estimate can be reused.
