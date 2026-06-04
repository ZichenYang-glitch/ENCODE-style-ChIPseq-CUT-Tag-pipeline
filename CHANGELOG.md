# Changelog

All notable changes to this project will be documented in this file.
The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- MNase-seq nucleosome positioning support (Stage 39): `assay: mnase` for PE
  MNase-seq. Produces mono-nucleosome BAM (alignmentSieve), dyad BigWig
  (bamCoverage --MNase), and mono occupancy BigWig at both sample and pooled
  experiment levels. Reuses existing preprocessing (FastQC, Trim Galore,
  Bowtie2, MAPQ filter, duplicate handling, final.bam).
- New `peak_mode: nucleosome` — required for MNase, rejected for non-MNase
  assays. PE-only constraint enforced in sample validation.
- Optional `mnase.mono_range` config key (default [140, 200]) for
  mono-nucleosome fragment length.
- MNase test profile (8th smoke profile), validation stress tests (7 new
  cases), and Stage 39 stress tests (28 cases).
- Manifest support for MNase outputs: mono BAM, dyad/mono BigWigs at sample
  and pooled level. Peak-centric outputs marked not_applicable for MNase rows.
- Optional `--strict-inputs` validation mode for FASTQ and Bowtie2 index
  file existence (Stage 30). Default is non-strict for dry-run compatibility;
  use `--strict-inputs` for pre-run or release validation.

### Changed
- `pipeline_done` rule in report.smk now conditionally gates peak-centric
  inputs on non-MNase samples.
- Target builder split into `PEAK_SAMPLE_IDS` (non-MNase treatment) and
  `MNASE_SAMPLE_IDS` (MNase treatment) with corresponding experiment-level
  subsets for pooled output gating.
- SE ChIP-seq MACS3 fragment-size fallback now emits a clear warning to stderr
  when using the 200 bp default, rather than failing silently.

### Documentation
- Documented strict input validation in README, configuration docs, release
  checklist, and assay policy.
- Public data execution report template and four per-queue stubs (Stage 32).
  All stubs state "not executed yet" — no data committed.
- `scripts/prepare_public_validation_inputs.py --report-stubs` prints
  queue names and expected report paths (stdlib-only, no downloads).
- Containerization planning (Stage 33): Apptainer/Docker strategy for v0.3-dev.
  Runner-only image first; full bioinformatics image deferred. Conda env files
  remain source of truth. No images built or committed.
- Runner container definition files (Stage 34): `containers/Dockerfile.runner`,
  `containers/Apptainer.runner.def`, `containers/README.md` with Docker and
  Apptainer usage examples. Static files only — build/smoke deferred to Stage 35.
  18 static-file contract tests. No images built.
- Docker runner build smoke report (Stage 35): build, pull, and dry-run
  verification documented. `/.cache` PermissionError root cause and
  HOME/XDG_CACHE_HOME fix documented in `containers/README.md`.
  `chipseq-runner:stage35-smoke` local image verified — no image artifacts
  committed.
- Singularity runner build smoke report (Stage 36): `singularity-container`
  installed as `singularity-ce version 4.1.1`; `containers/Apptainer.runner.def`
  built to a local `.sif`; Snakemake 8.30.0, PyYAML import, and bind-mounted
  14-job dry-run verified. Singularity HOME override warning documented; no
  image artifacts committed.
- Full 30-sample histone ChIP-seq real-run report (Stage 37): 5 conditions x
  3 marks x 2 biological replicates completed externally with sample-level and
  pooled outputs, MultiQC, result manifest, signal tracks, QC summaries, and
  documented follow-ups for bedGraph disk pressure and MultiQC naming.
- Container user experience polish (Stage 37): `docs/container-usage.md` full
  user guide with Docker/SingularityCE/Apptainer build, run, bind mount, HPC,
  and troubleshooting sections. `scripts/smoke_container_runner.sh` for
  automated container smoke tests (docker and singularity modes). Image
  publishing deferred.

## [v0.2.0-rc1] - 2026-05-24

### Added
- FE/ppois bedGraph-to-BigWig conversion for single-sample and pooled experiment
  signal tracks (gated on `genome_resources.<genome>.chrom_sizes`)
- Minimal project-level result manifest (`scripts/make_manifest.py`) recording core
  output existence with 10-column TSV. Uses `validate_samples` for DAG-consistent
  sample normalization and replicate/IDR eligibility gating
- QC summary assembly migrated from shell printf+tail/cut to Python scripts
  (`scripts/assemble_qc_summary.py`, `scripts/aggregate_qc_summary.py`), using
  header-based `csv.DictReader` parsing while preserving the 37-column output contract
- `scripts/prepare_public_validation_inputs.py` — stdlib-only dataset inventory helper
  with `--check-metadata`, `--json`, `--dry-run` modes (no downloads)
- v0.2 roadmap (`ROADMAP_v0.2.md`) and release checklist (`RELEASE_CHECKLIST.md`)
- Output contract document (`docs/output-contract.md`) with implemented manifest schema
- Conda install alternative in README

### Changed
- `rule all` target builder refactored into 9 helper functions (Stage 23) — zero DAG change
- Updated `docs/configuration.md` to document `chrom_sizes` role in BigWig gating
- Updated `docs/qc-interpretation.md` to reflect FE/ppois BigWig availability
- Updated `KNOWN_ISSUES.md` and README Limitations to reflect completed stages

### Documentation
- Assay policy contract (`docs/assay-policy.md`) documenting ChIP-seq, CUT&Tag, and
  ATAC-seq behavioral contracts from actual code
- IDR contract (`docs/idr-contract.md`) documenting Stage 5 eligibility, outputs, and
  scope boundaries
- Public data validation plan (`docs/release-checks/stage27-public-data-validation-plan.md`)
  with candidate datasets, expected outputs, and execution tiers
- Public metadata verification checklist and CI/CD expansion plan
  (`docs/release-checks/stage27b-metadata-ci-plan.md`) with acceptance criteria,
  artifact policy, and 3-tier CI/CD design

### Tests
- Stage 22 BigWig DAG gating stress tests (6 cases)
- Stage 24 QC summary unit tests (9 cases)
- Stage 25 manifest stress tests (14 cases) including replicate/IDR eligibility
- Stage 27a public validation plan tests (7 cases)
- Stage 27b metadata/CI plan tests (7 cases)
- Stage 27c CI workflow tests (7 cases)
- Stage 28 release readiness tests (11 cases)
- No-hardcoded-paths guard for all runtime and documentation files

### CI
- Fast PR checks expanded from 3 to 10 steps: validation plus 9 Python test suites
  covering all completed stages
- Manual `workflow_dispatch` tiny real execution retained
- Public data execution remains manual/external (not in CI)

## [v0.1.0-beta] - 2026-05-19

### Added
- Shared preprocessing pipeline: FastQC, Trim Galore, Bowtie2 alignment, MAPQ
  filtering, Picard duplicate handling, flagstat, idxstats, BigWig generation
- ChIP-seq / CUT&Tag assay-aware MACS3 peak calling with configurable parameters
- Optional control support via external BAM or FASTQ-based control sample rows
- Single-sample QC: FRiP, library complexity (Picard + NRF/PBC), MACS3 FE/ppois
  bedGraph signal tracks, per-sample QC summaries, and project-level aggregation
- Replicate-aware experiment model with technical replicate merging, biological
  replicate BAMs, and pooled treatment/control outputs
- TF ChIP-seq IDR: true-replicate IDR, pseudoreplicate-based self-IDR and
  pooled-IDR, reproducibility metrics, conservative/optimal peak sets
- Histone pooled QC summaries with target classification and peak-mode
  compatibility status
- CUT&Tag fragment-size QC and optional SEACR sidecar peak calling
- Configurable tool parameters for MACS3, bamCoverage, and IDR
- Schema-based config and sample sheet validation
- 7 test profiles covering major SE/PE/control/SEACR/IDR dispatch paths
- Tiny real-execution smoke test with synthetic data (no binary fixtures)
- GitHub Actions CI: fast validation + dry-run on PR/push, full execution on
  manual dispatch
- Streamlined environment setup with a full bioinformatics runtime environment
  and a lightweight CI environment
- Release checklist in README Developer Notes

### Documentation
- README with quick start, sample sheet reference, configuration guide,
  output structure, limitations, and developer notes
- KNOWN_ISSUES.md with full roadmap and current status
- Design specs and implementation plans under docs/superpowers/

### Known limitations
- Cross-correlation metrics (NSC/RSC) not yet implemented
- FE/ppois BigWig conversion not yet implemented
- IDR requires exactly 2 treatment biological replicates (no 3+ selection)
- CUT&Tag and broad-peak IDR not yet supported
- No bundled real public dataset
- ATAC-seq not supported

See README Limitations and KNOWN_ISSUES.md for full details.
