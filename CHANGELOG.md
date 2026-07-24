# Changelog

All notable changes to this project will be documented in this file.
The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

This is the curated release-candidate summary for v0.3.0. The release has not
yet been tagged or published.

### Added
- HelixWeave's workflow-neutral local platform: schema-driven input authoring,
  structured validation, immutable validated snapshots, file-backed SQLite
  lifecycle state, Redis/RQ execution, cancellation acknowledgement, run
  history, logs, indexed artifacts, machine-readable QC, and revision-bound
  downloads.
- Browser workflows for desktop and mobile covering workflow discovery,
  authoring, validation, run creation and start, run history, artifacts, QC,
  and downloads.
- The default `bulk-rnaseq` adapter, pinned to nf-core/rnaseq 3.26.0 at
  `e7ca46272c8f9d5ceee3f71759f4ba551d3217a4`, with authoring available without
  runtime assets and fail-closed execution admission.
- Protected synthetic execution acceptance for Bulk RNA-seq through the public
  product path, SQLite/RQ/Nextflow, STAR+Salmon, SortMeRNA, artifact/QC
  extraction, lifecycle, cancellation, timeout, and offline reuse contracts.
- A pinned, read-only Omics Intake Bundle 0.2 inspection boundary for the
  ENCODE adapter. It verifies the public contract and local file identities
  without creating snapshots or runs.
- ENCODE-style MNase-seq support, including fragment classes, dyad and occupancy
  tracks, sample and pooled outputs, QC summaries, and manifest coverage.
- A deterministic local ENCODE input-to-results demonstration and runtime
  doctor for workstation trial and operational checks.

### Changed
- The public product name is HelixWeave. The `encode-pipeline` distribution,
  `encode_pipeline` import namespace, and existing `encode-*` commands remain
  compatibility identities for v0.3.0.
- The default registry now contains both the ENCODE-style epigenomics adapter
  and Bulk RNA-seq. Shared API, lifecycle, artifact, QC, and frontend surfaces
  remain workflow-neutral.
- Workflow availability is explicit: authoring can remain available while
  execution is `not_configured` or `unavailable`; create and start recheck
  execution admission at backend boundaries.
- `mnase.mono_range` is deprecated in favor of `mnase.fragments.mono`;
  `fragments.mono` takes precedence when both are present.

### Fixed
- Public errors and workflow availability remain path-free and redact
  deployment coordinates, Redis credentials, environment values, command
  lines, and private adapter payloads.
- Artifact and QC generations, cursors, revisions, and download closure now
  fail closed across refresh, re-indexing, and stale requests.
- Local process ownership, timeout, and cancellation paths persist canonical
  terminal outcomes and clean up process groups and managed containers.

### Release boundary
- v0.3.0 does not publish a Docker, Apptainer, or Singularity image. The
  committed container definitions build an optional local runner for the
  ENCODE Snakemake workflow only.
- Bulk RNA-seq OCI images, JDK, Nextflow runtime, references, indexes, and
  fixture payloads remain operator-owned deployment assets and are not shipped
  in the Python distribution or source release.
- Synthetic gates prove execution and product contracts, not biological
  validity or production-scale performance.

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
