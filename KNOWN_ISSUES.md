# Known Issues and Follow-ups

This file records non-blocking issues and planned follow-up stages after the
Stage 1 Snakemake rule modularization and Stage 1.5 validation closeout.

## Current Status

- Core workflow, validation, replicate model, TF ChIP-seq IDR,
  histone/CUT&Tag support, baseline ATAC-seq support, smoke tests, and CI are
  implemented.
- FE/ppois BigWig conversion, QC summary Python refactor, target builder
  cleanup, minimal result manifest, and assay/IDR contract audit are
  complete (Stages 20-26).
- Remaining roadmap: optional release polish.
- ✅ Public data validation plan (Stage 27a — planned; no data downloaded).
- ✅ Metadata verification and CI/CD plan (Stage 27b — designed; CI Tier 1 wired in Stage 27c).
- ✅ CI/CD wiring (Stage 27c — fast-checks expanded to validation plus 9 Python test suites).

## Scope Gap

The current workflow is ENCODE-inspired, not a full ENCODE-compliant
ChIP-seq pipeline. It now provides:

- Single-sample preprocessing, MACS3 peak calling, CPM BigWig generation
- Single-sample QC: FRiP, library complexity (Picard + NRF/PBC),
  MACS3 FE/ppois bedGraph signal tracks
- Replicate-aware experiment model (Stage 4)
- TF ChIP-seq IDR for narrowPeak experiments with exactly 2 treatment
  biological replicates (Stage 5)
- Histone pooled QC summaries (Stage 6)
- CUT&Tag fragment-size QC and optional SEACR sidecar peaks (Stage 7)
- Optional TSS enrichment-style profiles from GTF annotations (Stage 18)
- Baseline ATAC-seq dispatch for narrowPeak MACS3 runs (Stage 19)
- Test profiles, CI, and execution documentation (Stage 8)

Still missing / deferred beyond v0.2:

- ✅ BigWig conversion for FE/ppois signal tracks. (Stage 22 — completed)
- ✅ QC summary Python refactor (Stage 24 — completed)
- ✅ Target builder cleanup (Stage 23 — completed)
- ✅ Minimal result manifest (Stage 25 — completed)
- ✅ Assay policy and IDR contract documentation (Stage 26 — completed)
- GC bias metrics (Picard CollectGcBiasMetrics)
- Browser-ready genome resource management

## Roadmap

### Stage 2: Schema and Genome Resources Foundation ✅

**Completed 2026-05-15.**

Goal: make incorrect configs fail early and give later ENCODE-like QC rules the
genome resources they need.

- ✅ Add `workflow/schemas/config.schema.yaml`.
- ✅ Add `workflow/schemas/samples.schema.yaml`.
- ✅ Add or formalize `scripts/validate_samples.py`.
- ✅ Validate required sample sheet columns and optional columns.
- ✅ Validate `layout`, `assay`, `peak_mode`, `role`, `trim`, `remove_dup`,
  `extend_reads`, and `use_control`.
- ✅ Validate `control_sample` references and `control_bam` paths when
  `use_control: true`.
- ✅ Keep backward compatibility with the current minimal sample sheet.
- ✅ Add genome resource config keys such as `chrom_sizes`, `blacklist`,
  `effective_genome_size`, and optional `gtf` / `reference_fasta`.
- ✅ Validate genome resource paths when downstream rules require them.
- ✅ Update `_normalize_genome()` to prefer configured effective genome sizes;
  mm39 resolves to explicit numeric value when configured.

### Stage 3: ENCODE-like Single-Sample QC

Goal: add the most valuable ENCODE-like QC metrics without changing the sample
model yet.

**Stage 3a completed 2026-05-15** — blacklist filtering, FRiP, peak counts,
and QC summary foundation.
**Stage 3b-1 completed 2026-05-15** — library complexity from Picard
MarkDuplicates metrics.
**Stage 3b-2 completed 2026-05-15** — MACS3 fold-enrichment and p-value
bedGraph signal tracks.
**Stage 3c-1 completed 2026-05-15** — BAM-derived NRF/PBC library complexity.

Estimated effort: 1-2 days (Stage 3a + 3b + 3c-1: ~1.25 days).

- ✅ Add blacklist filtering for BAMs and/or peak files. (Stage 3a)
- ✅ Add FRiP calculation. (Stage 3a, read-record based)
- ✅ Add sample-level summary TSV for alignment rate, duplicate rate, peak count,
  FRiP, library complexity, cross-correlation, and control usage. (Stage 3a foundation)
- ✅ Add library complexity metrics. (Stage 3b-1, duplication-derived from Picard)
- ✅ Add fold-enrichment and p-value signal tracks where appropriate. (Stage 3b-2,
  bedGraph output)
- ✅ Add NRF/PBC library complexity metrics. (Stage 3c-1, BAM-derived)
- ✅ Add pooled experiment FE/ppois signal tracks. (Stage 6a, completed 2026-05-17)
- ✅ Add cross-correlation metrics such as NSC/RSC. (Stage 12 / 3c-2, completed 2026-05-21)
- ✅ Add preseq-style library complexity extrapolation. (Stage 12 / 3c-2, completed 2026-05-21)
- ✅ Add Picard CollectMultipleMetrics (alignment summary, insert size, quality
  distribution). (Stage 12 / 3c-2, completed 2026-05-21)
- ✅ Add TSS enrichment-style profiles from GTF annotations. (Stage 18,
  completed 2026-05-22)
- ✅ Add bigWig conversion for FE/ppois signal tracks when chrom sizes and
  conversion tooling are configured. (Stage 22 — completed)
- ✅ Add project-level result manifest (`results/multiqc/result_manifest.tsv`) covering
  core per-sample, experiment-level, IDR, and project outputs with existence
  status. (Stage 25 — completed)
- ✅ Add a MultiQC custom config for cross-correlation summary visibility.
  (Stages 15-16)

### Stage 4: Replicate-Aware Sample Model

Goal: represent biological replicates, pooled samples, and treatment-control
relationships explicitly in the DAG.

**Stage 4a completed 2026-05-15** — sample sheet metadata foundation with
`experiment`, `condition`, `replicate`, `biological_replicate`, and
`technical_replicate` columns. Derived metadata structures
(`EXPERIMENT_IDS`, `SAMPLES_BY_EXPERIMENT`, etc.) are available in the
Snakefile for future stages.

Estimated effort: 1-2 days (Stage 4a: ~0.25 day).

- ✅ Extend the sample sheet with replicate/group concepts such as `replicate`,
  `condition`, or `experiment`. (Stage 4a, metadata-only)
- ✅ Keep backward compatibility for single-sample runs. (Stage 4a)
- ✅ Define grouped outputs for replicate sets. (Stage 4b, completed 2026-05-16)
- ✅ Define pooled replicate BAMs and pooled peak targets. (Stage 4b, completed 2026-05-16)
- ✅ Decide how technical replicates should be represented or merged. (Stage 4b, completed 2026-05-16)
  Technical replicates merged into biological-replicate BAMs; single tech-rep uses symlink.

### Stage 4b Limitations

The following are intentionally deferred to Stage 5:

- **No IDR or pseudoreplicates.** Pooled peak calling runs on the merged BAM.
  Self-pseudoreplicates, pooled pseudoreplicates, and IDR-based conservative/
  optimal peak sets are Stage 5.
- **No cross-replicate QC.** Per-biorep peak calling, reproducibility metrics,
  and replicate-level signal track comparison are not implemented.
- **Pooled signal tracks are Stage 6a.** MACS3 FE/ppois bedGraph generation
  is available for pooled experiment peaks when `qc.signal_tracks: true`.
- **Pooled control BAM merging is separate from treatment.** Pooled control
  BAMs are produced as independent targets. If both pooled treatment and
  pooled control are needed, Snakemake schedules both.

### Stage 5: TF ChIP-seq IDR and Reproducibility

Goal: implement the core ENCODE-style TF ChIP-seq reproducibility layer.

**Stage 5a completed 2026-05-17** — true-replicate IDR foundation with
config, validation, IDR-ready MACS3 calls, and raw/thresholded IDR output.

- ✅ Add individual replicate peak calling targets. (Stage 5a: IDR-ready MACS3 per biorep)
- ✅ Add IDR rules. (Stage 5a: true-replicate IDR, raw + thresholded)
- ✅ Add self-pseudoreplicate and pooled-pseudoreplicate generation. (Stage 5b, completed 2026-05-17)
- ✅ Emit conservative and optimal peak sets. (Stage 5b, completed 2026-05-17)
- ✅ Add reproducibility summary outputs. (Stage 5b, completed 2026-05-17)

**Stage 5b completed 2026-05-17:**
- ✅ Pseudoreplicate BAM generation (deterministic hash split script)
- ✅ Self-pseudoreplicate IDR (per biorep)
- ✅ Pooled-pseudoreplicate IDR
- ✅ Final conservative/optimal peak set assembly
- ✅ Reproducibility summary TSV

### Stage 6: Histone ChIP-seq Branch

Goal: support histone ChIP-seq expectations without forcing all histone assays
through TF-style IDR.

**Stage 6a completed 2026-05-17** — pooled experiment FE/ppois bedGraph
signal tracks for multi-biorep experiments.

Estimated effort: 1-2 days.

- Clarify narrow vs broad histone behavior.
- Add pooled replicate handling for histone marks.
- ✅ Add fold-enrichment and p-value signal tracks for pooled experiments.
- ✅ Add histone-appropriate reproducibility and QC summaries. (Stage 6b, completed 2026-05-17)

### Stage 7: CUT&Tag-Specific Branch

Goal: keep CUT&Tag useful while avoiding false equivalence with ENCODE ChIP-seq.

Estimated effort: 1-3 days.

- ✅ Add CUT&Tag fragment-size QC. (Stage 7a, completed 2026-05-18)
- ✅ Add optional SEACR sidecar peak calling. (Stage 7b, completed 2026-05-18)
- Add GoPeaks support as an alternative peak caller. (Stage 7c+)
- Add Tn5 insertion-aware metrics.
- Design optional spike-in normalization without breaking the current sample
  model.

### Stage 8: Test Profiles, CI, and Packaging

Goal: make the workflow easier to test, move, and run on different machines.

Estimated effort: 1-2 days.

**Stage 8a completed 2026-05-18** — test profiles and smoke harness.
**Stage 8b completed 2026-05-18** — tiny real execution smoke test.

- ✅ Add small test profiles for PE, SE, no-control, `control_bam`, and
  `control_sample`. (7 profiles under `test/profiles/`) (Stage 8a)
- ✅ Add a reproducible smoke-test command that runs quickly on a laptop.
  (`python3 test/test_stage8_smoke_profiles.py`, all dry-run) (Stage 8a)
- ✅ Include `chipseq` and `cuttag` dispatch coverage. (Stage 8a)
- ✅ Include baseline `atac` dispatch coverage. (Stage 19)
- ✅ Refine `.gitignore` so test profile files track normally. (Stage 8a)
- ✅ Add tiny FASTQ subsets and a tiny Bowtie2 index generated at runtime
  under `/tmp` (no committed binary fixtures). (Stage 8b)
- ✅ Add real end-to-end execution with tiny synthetic data (1 ChIP-seq PE
  no-control run, preprocessing + signal path only). (Stage 8b)
**Stage 8c completed 2026-05-18** — GitHub Actions CI wiring.

- ✅ GitHub Actions CI workflow (`.github/workflows/ci.yml`) with fast
  validation + dry-run checks on PR/push and manual real-execution via
  workflow_dispatch. (Stage 8c)
**Stage 8d completed 2026-05-19** — environment cleanup and execution
documentation.

- ✅ Remove local `prefix` metadata from exported Conda environment files.
  (Already absent from chipseq.yml and ci-fast.yml; confirmed.)
- ✅ Remove `defaults` channel from `workflow/envs/chipseq.yml`.
  Environment files use `conda-forge` + `bioconda` with `nodefaults`.
- ✅ Decide whether to split environments further: deferred in Stage 8d,
  then implemented in Stage 10e after solve/runtime pain persisted.
  The workflow now uses a lightweight runner plus rule-specific tool
  environments.
- ✅ Document local execution, validation, smoke tests, and full
  workflow run in `README.md`.

**Stage 10e completed 2026-05-19** — environment reliability fix.

- ✅ Replace the single-install path with `workflow/envs/runner.yml`
  (`snakemake-minimal`, Python, PyYAML) plus rule-specific tool envs.
- ✅ Isolate IDR in `workflow/envs/idr.yml` to contain its Bioconda
  `python <3.10` constraint.
- ✅ Isolate SEACR and MultiQC in their own envs so R/reporting stacks are
  not part of the first install.
- ✅ Rewire Snakemake `conda:` directives to the smallest appropriate env.

**Stage 10f completed 2026-05-20** — environment user documentation.

- ✅ Add `docs/environments.md` explaining runner, core, CI, and
  rule-specific environments.
- ✅ Document first-run `--use-conda` behavior, Snakemake env caching,
  cleanup commands, and troubleshooting.
- ✅ Link README and quickstart docs to the environment guide.

**Engineering hardening completed 2026-05-20** — no-hardcoding guard.

- ✅ Remove author-local Snakemake/samtools fallback paths from test harnesses.
- ✅ Add `test/test_no_hardcoded_paths.py` to prevent workstation-specific
  runtime paths and stale rule-env references.
- ✅ Add `docs/developer/no-hardcoding.md` as the path-handling policy.

## High Priority

1. ✅ Add schema/sample validation and genome resources. (Stage 2, completed 2026-05-15)
   - ✅ Add schemas for `config/config.yaml` and `config/samples.tsv`.
   - ✅ Validate optional columns: `role`, `control_sample`, and `control_bam`.
   - ✅ Add and validate genome resources needed by blacklist, FRiP, and signal
     track rules.
   - ✅ Keep backward compatibility with the current sample sheet.

2. Add ENCODE-like single-sample QC metrics. (Stage 3)
   - ✅ Add blacklist filtering.
   - ✅ Add FRiP and duplication-derived library complexity metrics.
   - ✅ Add a cross-sample summary TSV.
   - ✅ Add MACS3 FE/ppois bedGraph signal tracks.
   - ✅ Add cross-correlation, NRF/PBC, and preseq-style complexity metrics.
   - ✅ Add Picard CollectMultipleMetrics (alignment summary, insert size,
     quality distribution).
   - ✅ Stage 14: project-level cross-correlation summary
     (`results/multiqc/cross_correlation_summary.tsv`) and QC interpretation
     guide (`docs/qc-interpretation.md`).

3. ✅ Add synthetic tiny FASTQ/index smoke execution (one ChIP-seq PE
   no-control real execution). (Stage 8b)
   - Stage 8a covers control paths by dry-run.
   - Pending: optional small public real-data treatment/control
     validation set.

4. Keep user-facing examples synchronized.
   - `README.md`, `config/config.yaml`, and `config/samples.tsv` should agree.
   - Avoid committing machine-specific absolute paths in shared examples.
   - ✅ Stage 17 added `docs/reference-resources.md` with Bowtie2 index,
     FASTA `.fai`/`.dict`, chrom sizes, blacklist, GTF, and mm39/GRCm39
     preparation guidance.

## Medium Priority

5. Restore optional `plotFingerprint` QC.
   - The legacy `scripts/chipseq.sh` runs `plotFingerprint` when available.
   - The modular Snakemake workflow does not yet expose this as a rule.

6. ✅ Clean exported Conda environment metadata.
   - Completed in Stage 8d: `workflow/envs/chipseq.yml` and
     `workflow/envs/ci-fast.yml` contain no machine-specific `prefix`.
   - Environment files now use `conda-forge` + `bioconda` with `nodefaults`.

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

9. ✅ Split environments by responsibility.
   - Completed in Stage 10e: `workflow/envs/runner.yml` is the first-install
     environment, while rule-specific envs cover FastQC, trimming, alignment,
     samtools/bedtools, MACS3, deepTools, MultiQC, IDR, SEACR, and Python
     helper rules.

10. Decide how much control output to publish.
    - Stage 1 produces control bigWigs because they are useful for QC.
    - Future config may make control bigWig/report publication optional.

## Stage 13 Real-Data Validation (2026-05-22)

Status: **PASS with follow-ups**

- **CUT&Tag PE 2-replicate real run:** PASS. mm39/GRCm39, SMARCA5, no control,
  SEACR enabled. All outputs generated (final.bam, CPM BigWig, MACS3 narrowPeak,
  FE/ppois bedGraph, SEACR stringent BED, FRiP, NRF/PBC, library_complexity,
  qc_summary, fragment size QC, cross_correlation, preseq, Picard metrics,
  pooled outputs, MultiQC report).
- **ChIP-seq PE 2-replicate real run:** PASS. mm39/GRCm39, H3K27ac, no control.
  All outputs generated (same categories as CUT&Tag, minus SEACR).
- **Stage 12 QC outputs** generated on both assays. MultiQC captured Picard and
  preseq outputs.
- **Picard metric fixes applied** (see below).

### Issues found and fixed

1. **Missing `.dict` for Picard reference.** Picard CollectMultipleMetrics
   requires `reference_fasta` plus matching `.fai` and `.dict`. The initial
   run failed when `.dict` was absent.
2. **Picard STRICT validation failed on `INVALID_FLAG_MATE_UNMAPPED`.** Real PE
   BAMs after MAPQ/flag filtering can contain mate-field validation warnings.
   Fixed by adding `VALIDATION_STRINGENCY=LENIENT` to the
   `picard_collect_multiple_metrics` rule.

### Follow-ups

- **MultiQC visibility for phantompeakqualtools `cc.qc`:** project-level
  cross-correlation summary (`cross_correlation_summary.tsv`) added in
  Stage 14 and exposed through a MultiQC custom section in Stage 15. A native
  MultiQC plugin is not planned unless the project needs package-level report
  extensions.
- **PE mate-field cleanup:** evaluate `samtools fixmate` after MAPQ filtering
  if mate-field warnings become problematic for downstream tools beyond Picard.
- **No real data bundled in repo:** real FASTQ/BAM/BAI/BW/BDG/HTML/MultiQC
  outputs remain external and are not committed.

See `docs/release-checks/stage13-real-data-validation.md` for the full record.
See `docs/release-checks/stage16-real-run-report-audit.md` for the Stage 15
MultiQC custom-section audit on external CUT&Tag real-run outputs.

## Stage 18-19 TSS and ATAC Baseline (2026-05-22)

Status: **Implemented; dry-run validated**

- **TSS enrichment-style QC:** opt-in via `qc.tss_enrichment: true`. Requires
  `genome_resources.<genome>.gtf`, derives `results/reference/<genome>.tss.bed`,
  and produces per-sample deepTools matrix/profile outputs under
  `results/<sample>/05_qc/tss/`.
- **ATAC-seq baseline:** `assay: atac` is accepted for `peak_mode: narrow`.
  ATAC uses Tn5-aware MACS3 parameters, duplicate removal when
  `remove_dup: auto`, and the shared preprocessing/QC/report chain.
- **Scope boundary:** ATAC broad peaks, footprinting, nucleosome positioning,
  and ATAC-specific IDR are not implemented yet.

## Notes

- `control_bam` is checked only when `use_control: true`; when `use_control: false`,
  it remains ignored for backward compatibility.
- Control rows are preprocessed but do not run peak calling.
- SE ChIP-seq treatment samples may depend on MACS3 peak calling before bigWig
  generation so the fragment-size estimate can be reused.
