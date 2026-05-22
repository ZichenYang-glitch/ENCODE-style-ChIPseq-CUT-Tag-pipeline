# Configuration Reference

`config/config.yaml` controls workflow behavior, resource paths, and
feature gating. For the minimal config needed to start, see the
[README](../README.md#configuration).

## Core keys

| Key | Default | Accepted | Effect |
| :--- | :--- | :--- | :--- |
| `samples` | `"config/samples.tsv"` | path | Path to the sample sheet. |
| `outdir` | `"results"` | path | Output root directory. |
| `threads` | `8` | positive integer | Threads per job. |
| `mapq` | `30` | integer >= 0 | Minimum MAPQ for samtools filtering. |
| `binsize` | `10` | positive integer | Bin size in bp for bamCoverage BigWig. |
| `remove_dup` | `"auto"` | `"auto"`, `"yes"`, `"no"` | Duplicate removal mode. See below. |
| `trim` | `true` | boolean | Enable Trim Galore adapter/quality trimming. |
| `extend_reads` | `"auto"` | `"auto"`, `"yes"`, `"no"` | Extend SE reads to fragment size. |
| `use_control` | `false` | boolean | Enable control sample resolution. |
| `multiqc` | `true` | boolean | Enable MultiQC report aggregation. |

### `remove_dup` modes

| Value | Behavior |
| :--- | :--- |
| `"auto"` | Remove duplicates for `chipseq` narrow-peak; skip for `broad` peak mode and all `cuttag`. |
| `"yes"` | Always remove duplicates. The bundled core env uses `samtools markdup`; Picard is used only if available in a custom runtime environment. |
| `"no"` | Never remove duplicates. `final.bam` may be a symlink to the filtered BAM. |

## Genome resources

```yaml
genome_resources:
  hg38:
    effective_genome_size: 2913022398
    chrom_sizes: "/data/genomes/hg38/hg38.chrom.sizes"    # optional
    blacklist: "/data/genomes/hg38/hg38.blacklist.bed"    # optional
    gtf: ""        # optional
    reference_fasta: ""   # optional
```

- **`effective_genome_size`**: Either a MACS3 shortcut (`hs` = 2.7e9, `mm` = 1.87e9)
  or a positive integer. Used for MACS3 peak calling and bamCoverage normalization.
- **`chrom_sizes`**: Two-column file (chr, size). Used for bedGraph-to-BigWig
  conversion (not yet implemented).
- **`blacklist`**: BED file of problematic regions. Used when
  `qc.blacklist_filter: true`.
- **`gtf`, `reference_fasta`**: Reserved for future QC modules. Optional.

All optional paths, if non-empty, must exist on disk. The `genome` column in
`config/samples.tsv` must match a key in this block.

## `use_control`

**Default:** `false`

When `false`, `control_sample` and `control_bam` are ignored. Rows with
`role: control` are not scheduled as active outputs; MACS3 runs treatment
samples without a control (`-c` omitted).

When `true`, treatment rows may reference a control via `control_sample`
(another sample row) or `control_bam` (an external BAM path). Control rows
(`role: control`) are preprocessed but skip peak calling.

## `multiqc`

**Default:** `true`

When `true`, MultiQC aggregates QC artifacts from all active samples into
`results/multiqc/multiqc_report.html`.

MultiQC runs in `workflow/envs/multiqc.yml` when Snakemake is invoked with
`--use-conda`, keeping reporting dependencies out of the core runtime.

## QC block

```yaml
qc:
  blacklist_filter: true
  frip: true
  library_complexity: true
  nrf_pbc: true
  signal_tracks: true
  summary: true
  cuttag_fragment_size: true
  cross_correlation: false       # Stage 12 — opt-in
  preseq_complexity: false       # Stage 12 — opt-in
  picard_metrics: false          # Stage 12 — opt-in
```

| Switch | Default | Effect |
| :--- | :--- | :--- |
| `blacklist_filter` | `true` | Filter BAMs and peaks against the blacklist BED in genome resources. Requires a blacklist path in `genome_resources`. |
| `frip` | `true` | Compute FRiP (Fraction of Reads in Peaks) per sample. |
| `library_complexity` | `true` | Parse duplicate metrics. Picard fields are reported when Picard is available; otherwise fallback fields are emitted as `NA`. |
| `nrf_pbc` | `true` | Compute NRF and PBC1/PBC2 from the BAM file (library complexity from read counts). |
| `signal_tracks` | `true` | Produce MACS3 FE (fold-enrichment) and ppois (Poisson p-value) bedGraph tracks per treatment sample and for pooled experiments. BigWig conversion is not yet implemented. |
| `summary` | `true` | Emit per-sample QC summary TSV and a project-level aggregate at `results/multiqc/stage3_qc_summary.tsv`. |
| `cuttag_fragment_size` | `true` | Compute CUT&Tag fragment-size statistics for active samples with assay=cuttag. |
| `cross_correlation` | `false` | Run phantompeakqualtools cross-correlation QC per treatment sample. Produces NSC/RSC metrics (`.cc.qc`) and a cross-correlation plot (`.cc.plot.pdf`). |
| `preseq_complexity` | `false` | Run preseq library complexity extrapolation (`lc_extrap -B`) per treatment sample. Produces `.preseq.txt`. Complements existing NRF/PBC metrics. |
| `picard_metrics` | `false` | Run Picard CollectMultipleMetrics per treatment sample. Produces alignment summary, insert size, and quality distribution metrics. Requires `genome_resources.<genome>.reference_fasta`.

All QC switches default to `true`. Set individual switches to `false` to skip
specific metrics.

## Replicate and IDR features

```yaml
stage4b: true        # replicate-aware pooled BAMs and peaks (default on)
stage5: false        # TF ChIP-seq IDR (default off; requires stage4b: true)

idr:
  seed: 42
  threshold: 0.05    # IDR -log10(threshold) for thresholded peak set
  rank: "p.value"    # ranking measure: p.value or signal.value
```

- **`stage4b`** (default `true`): Enables technical replicate merging,
  biological-replicate BAMs, pooled treatment/control BAMs, pooled MACS3 peak
  calls, and pooled signal tracks. Applies to experiments with 2+ biological
  replicates.
- **`stage5`** (default `false`): Enables TF ChIP-seq IDR. Requires
  `stage4b: true`, `chipseq` assay, `narrow` peak_mode, and exactly 2
  treatment biological replicates per experiment.
- **`idr.*`**: The idr block is only read when `stage5: true`. `threshold` is
  the IDR threshold (0.05 = IDR -log10(0.05) ≈ 1.3). `rank` chooses between
  p-value and signal-value ranking for IDR input.
- IDR rules use `workflow/envs/idr.yml` when Snakemake is run with
  `--use-conda`, keeping IDR's older Python dependency constraints out of the
  core environment.

### IDR gating summary

| Condition | Required |
| :--- | :--- |
| `stage4b: true` | Always required for IDR |
| `assay` | `chipseq` only |
| `peak_mode` | `narrow` only |
| Biological replicates | Exactly 2 treatment bioreps per experiment |

CUT&Tag IDR and 3+ replicate IDR are not yet supported.

## CUT&Tag SEACR

```yaml
cuttag:
  seacr:
    enabled: false
```

- **`cuttag.seacr.enabled`** (default `false`): When `true`, runs SEACR peak
  calling as a sidecar alongside MACS3 for CUT&Tag samples. Outputs land under
  `results/<sample>/04_peaks_seacr/`.
- Only affects samples with `assay: cuttag`. ChIP-seq samples ignore this
  setting.
- SEACR rules use `workflow/envs/seacr.yml` when Snakemake is run with
  `--use-conda`, keeping SEACR's R runtime out of the core environment.

## Tool parameters

Optional `tool_parameters` blocks tune individual tool behavior. Absent keys
use built-in defaults.

```yaml
tool_parameters:
  macs3:
    qvalue: 0.01
  idr_macs3:
    pvalue: 0.1        # relaxed threshold for IDR-ready peak calls
  bamcoverage:
    normalize_using: "CPM"
```

| Block | Key | Default | Effect |
| :--- | :--- | :--- | :--- |
| `macs3` | `qvalue` | `0.01` | MACS3 q-value cutoff for peak calling. |
| `idr_macs3` | `pvalue` | `0.1` | Relaxed p-value for IDR-ready per-biorep peak calls. Only used when `stage5: true`. |
| `bamcoverage` | `normalize_using` | `"CPM"` | deepTools bamCoverage normalization method. |
| any block | `extra_args` | `""` | Additional CLI arguments passed through to the tool. Use sparingly. |

## Defaults and key dependencies

| Key / block | Default | Takes effect when |
| :--- | :--- | :--- |
| `use_control` | `false` | Always; enables control resolution when `true` |
| `multiqc` | `true` | Always; enables MultiQC aggregation when `true` |
| `qc.*` | all `true` | Always; individual switches gate each metric |
| `stage4b` | `true` | Always; enables pooled outputs for multi-biorep experiments |
| `stage5` | `false` | `stage4b: true` + chipseq + narrow + exactly 2 bioreps |
| `idr.*` | listed above | Only when `stage5: true` |
| `cuttag.seacr.enabled` | `false` | Always; affects CUT&Tag samples only |
| `tool_parameters.*` | tool defaults | Always; absent keys use built-in tool defaults |
