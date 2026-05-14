# Known Issues and Follow-ups

This file records non-blocking issues and planned follow-up stages after the
Stage 1 Snakemake rule modularization and Stage 1.5 validation closeout.

## Current Status

- The modular Snakemake DAG is in place under `workflow/rules/`.
- `workflow/Snakefile` is now the orchestration entry point.
- ChIP-seq and CUT&Tag policy functions are separated in assay-specific rule files.
- Control samples can be modeled as first-class sample rows through `role` and
  `control_sample`.
- Dry-run and rule-list validation pass for the current test configuration.
- Real end-to-end validation has been run on bundled test FASTQs.
- The bundled test FASTQs align much better to GRCh38 than GRCm39, so the test
  data should be treated as human data.
- Trim Galore 2.2 text and JSON reports are both preserved; MultiQC 1.35
  recognizes Trim Galore v2+ reports through the JSON files.

## Scope Gap

The current workflow is ENCODE-inspired, not a full ENCODE-compliant ChIP-seq
pipeline. It currently provides single-sample preprocessing, MACS3 peak calling,
CPM-normalized BigWig generation, and MultiQC aggregation.

Major ENCODE-aligned features still missing include:

- Replicate-aware experiment modeling.
- IDR and pseudoreplicate analysis for TF ChIP-seq.
- Pooled replicate peak sets.
- Blacklist filtering.
- FRiP calculation.
- Library complexity metrics.
- Cross-correlation metrics such as NSC/RSC.
- Fold-enrichment and p-value signal tracks.
- Reproducibility reports.
- Genome-specific resource management, including chromosome sizes, blacklist
  BED files, effective genome sizes, and optional reference annotations.

## Roadmap

### Stage 2: Schema and Genome Resources Foundation

Goal: make incorrect configs fail early and give later ENCODE-like QC rules the
genome resources they need.

Estimated effort: 0.5-1 day.

- Add `workflow/schemas/config.schema.yaml`.
- Add `workflow/schemas/samples.schema.yaml`.
- Add or formalize `scripts/validate_samples.py`.
- Validate required sample sheet columns and optional columns.
- Validate `layout`, `assay`, `peak_mode`, `role`, `trim`, `remove_dup`,
  `extend_reads`, and `use_control`.
- Validate `control_sample` references and `control_bam` paths when
  `use_control: true`.
- Keep backward compatibility with the current minimal sample sheet.
- Add genome resource config keys such as `chrom_sizes`, `blacklist`,
  `effective_genome_size`, and optional `gtf` / `reference_fasta`.
- Validate genome resource paths when downstream rules require them.

### Stage 3: ENCODE-like Single-Sample QC

Goal: add the most valuable ENCODE-like QC metrics without changing the sample
model yet.

Estimated effort: 1-2 days.

- Add blacklist filtering for BAMs and/or peak files.
- Add FRiP calculation.
- Add library complexity metrics.
- Add cross-correlation metrics such as NSC/RSC.
- Add fold-enrichment and p-value signal tracks where appropriate.
- Add sample-level summary TSV for alignment rate, duplicate rate, peak count,
  FRiP, library complexity, cross-correlation, and control usage.
- Add a MultiQC custom config if needed to improve naming and report layout.

### Stage 4: Replicate-Aware Sample Model

Goal: represent biological replicates, pooled samples, and treatment-control
relationships explicitly in the DAG.

Estimated effort: 1-2 days.

- Extend the sample sheet with replicate/group concepts such as `replicate`,
  `condition`, or `experiment`.
- Keep backward compatibility for single-sample runs.
- Define grouped outputs for replicate sets.
- Define pooled replicate BAMs and pooled peak targets.
- Decide how technical replicates should be represented or merged.

### Stage 5: TF ChIP-seq IDR and Reproducibility

Goal: implement the core ENCODE-style TF ChIP-seq reproducibility layer.

Estimated effort: 2-4 days.

- Add individual replicate peak calling targets.
- Add pooled replicate peak calling targets.
- Add self-pseudoreplicate and pooled-pseudoreplicate generation.
- Add IDR rules.
- Emit conservative and optimal peak sets.
- Add reproducibility summary outputs.

### Stage 6: Histone ChIP-seq Branch

Goal: support histone ChIP-seq expectations without forcing all histone assays
through TF-style IDR.

Estimated effort: 1-2 days.

- Clarify narrow vs broad histone behavior.
- Add pooled replicate handling for histone marks.
- Add fold-enrichment and p-value signal tracks.
- Add histone-appropriate reproducibility and QC summaries.

### Stage 7: CUT&Tag-Specific Branch

Goal: keep CUT&Tag useful while avoiding false equivalence with ENCODE ChIP-seq.

Estimated effort: 1-3 days.

- Add CUT&Tag-specific peak callers such as SEACR or GoPeaks.
- Add fragment-size QC.
- Add Tn5 insertion-aware metrics.
- Design optional spike-in normalization without breaking the current sample
  model.

### Stage 8: Test Profiles, CI, and Packaging

Goal: make the workflow easier to test, move, and run on different machines.

Estimated effort: 1-2 days.

- Add small test profiles for PE, SE, no-control, `control_bam`, and
  `control_sample`.
- Add or document tiny FASTQ subsets and a tiny Bowtie2 index suitable for CI.
- Add a reproducible smoke-test command that runs quickly on a laptop.
- Include both `chipseq` and `cuttag` dispatch coverage.
- Remove local `prefix` metadata from exported Conda environment files.
- Decide whether to split `workflow/envs/chipseq.yml` into smaller
  responsibility-specific environments.
- Document recommended execution on workstation/server environments.

## High Priority

1. Add schema/sample validation and genome resources.
   - Add schemas for `config/config.yaml` and `config/samples.tsv`.
   - Validate optional columns: `role`, `control_sample`, and `control_bam`.
   - Add and validate genome resources needed by blacklist, FRiP, and signal
     track rules.
   - Keep backward compatibility with the current sample sheet.

2. Add ENCODE-like single-sample QC metrics.
   - Add blacklist filtering.
   - Add FRiP, library complexity, and cross-correlation metrics.
   - Add a cross-sample summary TSV.

3. Add small, reproducible smoke-test data.
   - The current bundled FASTQs are useful but large.
   - Add tiny FASTQ subsets and a small index for quick validation.
   - Include at least one treatment/control sample pair.

4. Keep user-facing examples synchronized.
   - `README.md`, `config/config.yaml`, and `config/samples.tsv` should agree.
   - Avoid committing machine-specific absolute paths in shared examples.

## Medium Priority

5. Restore optional `plotFingerprint` QC.
   - The legacy `scripts/chipseq.sh` runs `plotFingerprint` when available.
   - The modular Snakemake workflow does not yet expose this as a rule.

6. Clean exported Conda environment metadata.
   - `chipseq.yml` may contain a machine-specific `prefix`.
   - Remove local prefixes before sharing or publishing the environment file.

7. Harden the legacy single-sample script.
   - `scripts/chipseq.sh` still uses input-derived Trim Galore output discovery
     with `ls ... | tail -n 1`.
   - This can select stale files if an output directory is reused.
   - The Snakemake path normalization avoids this issue, but the legacy script
     remains available for compatibility.

8. Reduce silent error handling.
   - Some optional QC steps intentionally continue on failure.
   - This is convenient for exploratory runs, but can hide missing QC in stricter
     production contexts.

## Low Priority

9. Split environments by responsibility.
   - Stage 1 keeps using `workflow/envs/chipseq.yml`.
   - Future stages can split this into core, ChIP-seq, CUT&Tag, and reporting
     environments if dependency resolution becomes slow or fragile.

10. Decide how much control output to publish.
    - Stage 1 produces control bigWigs because they are useful for QC.
    - Future config may make control bigWig/report publication optional.

## Notes

- `control_bam` is checked only when `use_control: true`; when `use_control: false`,
  it remains ignored for backward compatibility.
- Control rows are preprocessed but do not run peak calling.
- SE ChIP-seq treatment samples may depend on MACS3 peak calling before bigWig
  generation so the fragment-size estimate can be reused.
