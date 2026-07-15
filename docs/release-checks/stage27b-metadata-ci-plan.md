# Stage 27b: Public Metadata Verification + CI/CD Plan

**Date:** 2026-05-24
**Status:** historical release-planning evidence

> This document records the CI layout and test inventory used for the
> 2026-05-24 release checkpoint. It is not the current CI definition. The
> stage-numbered Stage 27 test modules listed below were later retired during
> the maintenance baseline; current offline behavior coverage for the public
> validation inventory lives in
> `test/scripts/test_prepare_public_validation_inputs.py`.

## Metadata Verification

### Accession Metadata Checklist

Before downloading any data, verify the following fields for each candidate dataset against official ENCODE or GEO records:

| Accession | Fields to verify |
| :--- | :--- |
| ENCSR000DYI | assay=ChIP-seq, target=CEBPB, organism=human, genome=hg38, 2 treatment bioreps + 1 control, PE |
| ENCSR000AKB | assay=ChIP-seq, target=H3K27me3, organism=human, genome=hg38, 2 treatment bioreps + 1 control, PE, broad mark |
| ENCSR254KDA | assay=ATAC-seq, organism=human, genome=hg38, 2 replicates, PE, no control |
| GSE145187 | assay=CUT&Tag, target=H3K27me3 (or H3K4me3), 2 treatment + IgG control, PE |

Use `scripts/prepare_public_validation_inputs.py --check-metadata` to print the full checklist (no network requests).

### ENCODE Verification Steps

1. Navigate to `https://www.encodeproject.org/experiments/<accession>/`.
2. Confirm the experiment page lists the expected assay, target, biosample, and organism.
3. Open the "Files" tab and confirm FASTQ files exist (not just processed BAM/bigWig/peaks).
4. Note the exact file download URLs or ENCODE file accession IDs.
5. Record read length (e.g., 2×101) and library strategy.

### GEO Verification Steps

1. Navigate to `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=<accession>`.
2. Confirm the Series record lists the expected samples and targets.
3. Open the SRA Run Selector to confirm FASTQ availability and layout (PE/SE).
4. Note SRR run accessions for downstream `fasterq-dump` or `prefetch`.

### Dataset Acceptance/Rejection Criteria

| Criterion | Accept | Reject |
| :--- | :--- | :--- |
| Assay matches expected | Yes | — |
| Target matches expected | Yes | — |
| Genome build | hg38 (or mm10 for CUT&Tag alternative) | Other builds (unless explicit re-mapping planned) |
| Replicate count matches | Yes | If count differs, re-evaluate sample sheet |
| Publicly accessible | Yes | dbGaP-controlled or login-walled datasets |
| FASTQ available | Yes | BAM-only or processed-signal-only |

## CI/CD Tiers

### Historical Tier 1 snapshot (Stage 27c)

**Trigger:** `push`, `pull_request` to `main`, `stage*`

The following commands were wired into `.github/workflows/ci.yml` at that
release checkpoint:

**Historical command inventory:**
- `python3 scripts/validate_samples.py --config config/config.yaml`
- `python3 test/test_validation_stress.py` (15 tests)
- `python3 test/test_no_hardcoded_paths.py`
- `python3 test/test_stage8_smoke_profiles.py` (7 profiles, dry-run)
- `python3 test/test_stage22_bigwig_stress.py` (6 tests)
- `python3 test/test_stage24_qc_summary_unit.py` (9 tests)
- `python3 test/test_stage25_manifest_stress.py` (14 tests)
- `python3 test/test_stage27_public_validation_plan.py` (7 tests)
- `python3 test/test_stage27b_metadata_ci_plan.py` (new in Stage 27b)
- `python3 test/test_stage27c_ci_workflow.py` (wired in Stage 27c)

**Environment:** `workflow/envs/ci-fast.yml` (Python + PyYAML + Snakemake-minimal)
**Budget:** <5 minutes per run

### Tier 2: Manual workflow_dispatch (on request)

**Trigger:** `workflow_dispatch` in GitHub Actions

**Runs:**
- All Tier 1 checks
- `python3 test/test_stage8b_tiny_execution.py` (real execution with synthetic data)
- Stage-specific stress tests (4b, 4c, 5a, 5b, 6a, 6b, 7a, 7b, 12, 18-19)

**Environment:** `workflow/envs/chipseq.yml` (full bioinformatics tools)
**Budget:** <15 minutes per run

### Tier 3: External/Manual Public Data Run

**Trigger:** Manual, external to CI

**Runs:**
- Full pipeline on downsampled or full public data
- Audit outputs against the checklist in `docs/release-checks/stage27-public-data-validation-plan.md`
- Record QC summary values, manifest output, and IDR reproducibility metrics

**Environment:** User's workstation or compute cluster
**Budget:** 2-8 hours per queue (16 cores), 20-100 GB storage per queue

### Why Public Data Is Not Run on Every PR

- **Data size:** FASTQ files are 1-50 GB per queue. CI runners lack this storage.
- **Runtime:** A full workflow run takes hours, not minutes.
- **Cost:** GitHub Actions minutes and egress bandwidth are limited.
- **Reproducibility:** Public data URLs can change or require authentication tokens. CI environments are not suitable for long-lived data fetches.
- **Intent:** Public data validation is a release-gating activity, not a per-commit check.

## Artifact Policy

### Committed to repo (allowed)
- This plan document
- `scripts/prepare_public_validation_inputs.py`
- Per-queue QC summary audits (small TSV, <10 KB)
- IDR reproducibility summaries (small TSV, <1 KB)
- Manifest outputs from validation runs (small TSV, <5 KB)
- Validation run reports (`docs/release-checks/stage27-<queue>-validation-report.md`)

### NOT committed (enforced by `.gitignore` and tests)
- FASTQ / BAM / BAI / BW / BDG files
- MultiQC HTML reports
- Bowtie2 index files
- Full pipeline output directories

## Storage / Runtime Budget

| Tier | Storage (per queue) | Runtime (per queue) | Who Runs |
| :--- | :--- | :--- | :--- |
| Fast PR | 0 (dry-run only) | <5 min | CI |
| workflow_dispatch | ~200 MB (synthetic data) | <15 min | CI (manual trigger) |
| Downsampled public | ~1-5 GB | ~30-60 min | User |
| Full public | ~20-100 GB | 2-8 hours | User (release gate) |
