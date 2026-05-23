# v0.2 Roadmap

This roadmap defines the scope and sequence for the **v0.2 release** of the ENCODE-style ChIP-seq / CUT&Tag / ATAC-seq pipeline. The corrected implementation strategy is documented in [`deep-research-report.md`](deep-research-report.md).

## v0.2 Scope

### In-scope

| Area | Deliverable |
| :--- | :--- |
| Assay support | ChIP-seq, CUT&Tag, baseline ATAC-seq (semantic freeze) |
| Peak calling | MACS3 narrow/broad, CUT&Tag SEACR sidecar |
| Replicate model | Pooled BAM/peaks/signal, TF ChIP-seq exactly-2-biorep IDR |
| Signal tracks | CPM BigWig, FE/ppois bedGraph, FE/ppois BigWig (chrom_sizes-gated) |
| QC | Full single-sample QC, cross-correlation, preseq, Picard, TSS enrichment |
| Engineering | Schema/validation, test profiles, stress tests, CI, target builder cleanup |
| Documentation | Assay policy, IDR contract, output contract, manifest schema, configuration and QC interpretation guides |
| Deployment | Conda envs as truth source; containerization plan (post-v0.2) |

### Out-of-scope

| Area | Deferred |
| :--- | :--- |
| Assay support | RNA-seq, Hi-C, single-cell, differential analysis |
| Peak calling | Multi-peak-caller comparison framework |
| Replicate model | 3+ replicate IDR, broad-mark IDR, CUT&Tag IDR, ATAC IDR |
| QC | GC bias metrics, LIMS-style auto-pass/fail, `plotFingerprint` |
| Engineering | Frontend GUI, database backend, user permission system |
| Data | Bundled public datasets (reference validation only, no committed FASTQ/BAM/BW) |
| Deployment | Full server deployment, track hub auto-publish |

### Design decisions preserved

- Continue using existing `config/config.yaml`, `config/samples.tsv`, `genome_resources`, `tool_parameters` model.
- No new config schema (`assay: chip_tf`, `references:`, `peak_calling:`).
- Current directory layout: `workflow/rules/`, `workflow/schemas/`, `test/`, `test/profiles/`.

## Stage Sequence

```
Stage 20a → 20b → 22 → 23 → 24 → 25 → 26 → 27a → 27b → 27c
  │         │      │      │      │      │      │       │       │
  docs    contract BigWig  target  QC    manifest assay/  public  CI/CD
roadmap                    builder summary        IDR     data   wiring
                           refactor               contract plan
```

All stages 20-27c are complete:

| Stage | Name | Status |
| :--- | :--- | :--- |
| 20a | v0.2 roadmap and documentation | ✅ Completed |
| 20b | Output contract | ✅ Completed |
| 22 | FE/ppois BigWig conversion | ✅ Completed |
| 23 | Target builder cleanup | ✅ Completed |
| 24 | QC Summary Python refactor | ✅ Completed |
| 25 | Minimal result manifest | ✅ Completed |
| 26 | Assay/IDR contract audit + roadmap sync | ✅ Completed |
| 27a | Public data validation plan | ✅ Completed |
| 27b | Metadata verification + CI/CD plan | ✅ Completed |
| 27c | CI/CD wiring | ✅ Completed |

### Remaining

| Stage | Name | Description |
| :--- | :--- | :--- |
| — | Public data execution | Manual downsampled/full runs on validation datasets (external, not automated) |
| — | Release polish | Optional: version bump, release notes, CHANGELOG consolidation, final repo hygiene |

## Current State

As of 2026-05-24, the project has:
- Full ChIP-seq / CUT&Tag / ATAC-seq preprocessing and peak calling with assay-specific policies
- Single-sample QC (FRiP, NRF/PBC, library complexity, peak counts, QC summaries) with Python-based assembly
- Opt-in QC (cross-correlation, preseq, Picard metrics, TSS profiles)
- Replicate-aware pooled outputs (BAM, peaks, signal) with correct biological-replicate gating
- TF ChIP-seq exactly-2-biorep IDR with conservative/optimal peak sets and reproducibility summary
- FE/ppois BigWig conversion (gated on `genome_resources.<genome>.chrom_sizes`)
- Minimal result manifest (`result_manifest.tsv`) with 10-column schema and DAG-consistent gating
- 7 smoke test profiles, CI (GitHub Actions), validation framework, 14 stress test suites
- Assay policy and IDR contract documentation (`docs/assay-policy.md`, `docs/idr-contract.md`)

See [`KNOWN_ISSUES.md`](KNOWN_ISSUES.md) for the full feature status and planned follow-ups.
