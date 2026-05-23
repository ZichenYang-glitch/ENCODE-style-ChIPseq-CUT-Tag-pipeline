# Stage 25: Minimal Manifest — Design Spec

**Date:** 2026-05-24
**Status:** design review
**Scope:** Minimal project-level result manifest recording core outputs, paths, methods, status, and QC flags
**Prerequisite:** Stage 20-24 completed and verified
**Excluded:** Heavy provenance system, md5 by default, JSON output, QC threshold auto-judgment, LIMS integration

---

## 1. Goals

1. Create `scripts/make_manifest.py` that enumerates project outputs and records their existence status.
2. Add a `result_manifest` rule producing `results/multiqc/result_manifest.tsv`.
3. Cover core single-sample, experiment-level, IDR, and project-level outputs.
4. Use the manifest field schema from `docs/output-contract.md`.
5. Keep md5 off by default; provide `--md5` flag for opt-in hashing.
6. Deterministic output: stable row ordering across runs.

---

## 2. Design Questions

### Q1: Which output types does the minimal manifest cover?

From `docs/output-contract.md`, the manifest covers all stable and Stage 22 implemented types:

**Per-sample (always applicable):**
- `final_bam`, `final_bai`, `cpm_bigwig`, `macs3_peak`, `qc_summary`

**Per-sample (gated on `qc.signal_tracks`):**
- `macs3_fe_bdg`, `macs3_ppois_bdg`

**Per-sample (gated on `chrom_sizes` + `signal_tracks`):**
- `macs3_fe_bw`, `macs3_ppois_bw`

**Per-experiment (gated on `stage4b` + multi-biorep):**
- `pooled_final_bam`, `pooled_final_bai`, `biorep_final_bam`, `biorep_final_bai`, `pooled_macs3_peak`, `pooled_qc_summary`

**Per-experiment (gated on `stage4b` + `signal_tracks`):**
- `pooled_fe_bdg`, `pooled_ppois_bdg`

**Per-experiment (gated on `stage4b` + `chrom_sizes` + `signal_tracks`):**
- `pooled_fe_bw`, `pooled_ppois_bw`

**IDR (gated on `stage5`):**
- `idr_conservative`, `idr_optimal`, `idr_reproducibility_summary`

**Project-level:**
- `stage3_qc_summary`, `multiqc_report`

**Deferred (not in v0.2 minimal manifest):**
- `cross_correlation`, `preseq`, `picard_*`, `tss_*` (opt-in QC; manifest them as `not_applicable` when gated off, or omit)
- `cuttag_fragment_size`, `seacr_*` (CUT&Tag-specific)
- `blacklist_filtered_bam`, `peak_counts`, `frip`, `library_complexity`, `nrf_pbc` (intermediate QC; already summarized in qc_summary)

### Q2: Per-sample vs per-experiment vs project-level

See Q1 enumeration. The manifest iterates all treatment samples and all multi-biorep experiments, plus project-level targets.

### Q3: Manifest field schema

Adapted from `docs/output-contract.md` draft, simplified for v0.2:

```tsv
sample_id  experiment_id  assay  target  genome  output_type  method  path  status  qc_flag
```

- `experiment_id`: empty for per-sample and project-level rows
- `sample_id`: empty for per-experiment and project-level rows
- `path`: relative to `outdir`

### Q4: `status` definition

| Value | Meaning |
| :--- | :--- |
| `present` | File/directory exists on disk |
| `missing` | Output expected but not found |
| `not_applicable` | Feature gated off (e.g., `stage5: false` → all IDR outputs are NA) |

### Q5: `qc_flag` — first version

**First version: `NA` for all rows.** Auto-judgment from qc_summary thresholds (e.g., FRiP < 0.01 → `fail`) is deferred to a future stage. The manifest's job in v0.2 is to record existence.

### Q6: md5

**Not calculated by default.** Large BAM/BW files would make this prohibitively slow. Add a `--md5` CLI flag (off by default). When enabled, compute hex digest for all `present` files. The `method` field records whether md5 was computed (e.g., `bowtie2+samtools` vs `bowtie2+samtools+md5`).

### Q7: Hard fail vs mark missing?

- **Per-row:** mark `missing` — the manifest is informational.
- **Exit code:** `0` by default (always succeeds). `--strict` flag makes exit code non-zero if any row has `status: missing`.
- **Rationale:** The manifest runs as a post-processing rule. A missing file in the manifest should not crash Snakemake.

### Q8: Keeping Stage 22/23/24 unchanged

- `make_manifest.py` is a new script — no changes to existing rules or scripts.
- New `result_manifest` rule added to `workflow/rules/report.smk`.
- New manifest target added to `_replicate_targets()` or a lightweight new helper in Snakefile.
- File existence checks use `os.path.exists()` — no changes to pipeline DAG.

---

## 3. Architecture

### `scripts/make_manifest.py`

```
Usage:
    python3 scripts/make_manifest.py \
        --config config/config.yaml \
        --output results/multiqc/result_manifest.tsv \
        [--md5] [--strict]
```

**Data flow:**
1. Read `config.yaml` → extract `outdir`, `qc.signal_tracks`, `stage4b`, `stage5`, `multiqc`, genome_resources
2. Read `samples.tsv` → extract treatment samples, experiments, IDs, assay, target, genome
3. Build rows list: iterate samples, then experiments, then project-level
4. For each row, resolve path relative to outdir, check `os.path.exists()`
5. Write TSV via `csv.DictWriter` with `lineterminator="\n"`

**Row ordering (deterministic):**
1. Per-sample rows sorted by `sample_id`
2. Per-experiment rows sorted by `experiment_id`
3. Project-level rows at end

### New rule: `result_manifest`

Location: `workflow/rules/report.smk`

```python
rule result_manifest:
    output:
        f"{OUTDIR}/multiqc/result_manifest.tsv"
    input:
        # Depends on stage3_qc_summary (marks all per-sample QC complete)
        f"{OUTDIR}/multiqc/stage3_qc_summary.tsv"
    params:
        config = "config/config.yaml"
    conda:
        "../envs/python.yml"
    shell:
        """
        mkdir -p "$(dirname {output:q})"
        python3 scripts/make_manifest.py \
            --config {params.config:q} \
            --output {output:q}
        """
```

### Target in Snakefile

Add to `_base_targets()` or as a standalone target after `_replicate_targets()`:

```python
if QC_CONFIG.get("summary", True) and TREATMENT_SAMPLE_IDS:
    manifest = [f"{OUTDIR}/multiqc/result_manifest.tsv"]
else:
    manifest = []
```

---

## 4. Out of Scope

- JSON output (Stage 25+)
- md5 by default
- QC threshold auto-judgment from qc_summary.tsv
- Manifest for opt-in QC modules (cross_corr, preseq, Picard, TSS)
- `--strict` exit code mode (plumbing only; full implementation deferred)
- LIMS or database integration
- Manifest web frontend

---

## 5. Verification

```bash
python3 scripts/validate_samples.py --config config/config.yaml
python3 test/test_validation_stress.py
python3 test/test_stage8_smoke_profiles.py
python3 test/test_stage22_bigwig_stress.py
python3 test/test_stage24_qc_summary_unit.py
python3 test/test_stage25_manifest_stress.py  # new
```
