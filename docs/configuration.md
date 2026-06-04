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
| `"auto"` | Remove duplicates for narrow-peak samples; skip for `broad` peak mode. |
| `"yes"` | Always remove duplicates. The bundled core env uses `samtools markdup`; Picard is used only if available in a custom runtime environment. |
| `"no"` | Never remove duplicates. `final.bam` may be a symlink to the filtered BAM. |

## Genome resources

```yaml
genome_resources:
  hg38:
    effective_genome_size: 2913022398
    chrom_sizes: "/data/genomes/hg38/hg38.chrom.sizes"    # optional
    blacklist: "/data/genomes/hg38/hg38.blacklist.bed"    # optional
    gtf: ""        # required when qc.tss_enrichment: true
    reference_fasta: ""   # optional
```

- **`effective_genome_size`**: Either a MACS3 shortcut (`hs` = 2913022398,
  `mm` = 2652783500)
  or a positive integer. Used for MACS3 peak calling and bamCoverage normalization.
- **`chrom_sizes`**: Two-column file (chr, size). Required for FE/ppois
  BigWig conversion. When non-empty and `qc.signal_tracks: true`, FE/ppois
  bedGraph tracks are also converted to BigWig format. When empty, only
  bedGraph signal tracks are produced. See the BigWig section below.
- **`blacklist`**: BED file of problematic regions. Used when
  `qc.blacklist_filter: true`.
- **`gtf`**: Annotation file required when `qc.tss_enrichment: true`; used to
  derive transcription start sites for deepTools TSS profile QC.
- **`reference_fasta`**: Required when `qc.picard_metrics: true`; must have a
  matching `.fai` and `.dict`.

All optional paths, if non-empty, must exist on disk. The `genome` column in
`config/samples.tsv` must match a key in this block.
For preparation commands and mm39/GRCm39 examples, see
[docs/reference-resources.md](reference-resources.md).

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
  tss_enrichment: false          # Stage 18 — opt-in; requires GTF
```

| Switch | Default | Effect |
| :--- | :--- | :--- |
| `blacklist_filter` | `true` | Filter BAMs and peaks against the blacklist BED in genome resources. Requires a blacklist path in `genome_resources`. |
| `frip` | `true` | Compute FRiP (Fraction of Reads in Peaks) per sample. |
| `library_complexity` | `true` | Parse duplicate metrics. Picard fields are reported when Picard is available; otherwise fallback fields are emitted as `NA`. |
| `nrf_pbc` | `true` | Compute NRF and PBC1/PBC2 from the BAM file (library complexity from read counts). |
| `signal_tracks` | `true` | Produce MACS3 FE (fold-enrichment) and ppois (Poisson p-value) bedGraph tracks per treatment sample and for pooled experiments. FE/ppois BigWig conversion is available when `genome_resources.<genome>.chrom_sizes` is also configured. |
| `summary` | `true` | Emit per-sample QC summary TSV and a project-level aggregate at `results/multiqc/stage3_qc_summary.tsv`. |
| `cuttag_fragment_size` | `true` | Compute CUT&Tag fragment-size statistics for active samples with assay=cuttag. |
| `cross_correlation` | `false` | Run phantompeakqualtools cross-correlation QC per treatment sample. Produces NSC/RSC metrics (`.cc.qc`) and a cross-correlation plot (`.cc.plot.pdf`). When enabled, also generates a project-level summary at `results/multiqc/cross_correlation_summary.tsv` and exposes it as a MultiQC custom section. See [docs/qc-interpretation.md](qc-interpretation.md) for interpretation guidance. |
| `preseq_complexity` | `false` | Run preseq library complexity extrapolation (`lc_extrap -B`) per treatment sample. Produces `.preseq.txt`. Complements existing NRF/PBC metrics. |
| `picard_metrics` | `false` | Run Picard CollectMultipleMetrics per treatment sample. Produces alignment summary, insert size, and quality distribution metrics. Requires `genome_resources.<genome>.reference_fasta` with a matching samtools FASTA index (`.fai`) and Picard sequence dictionary (`.dict`) next to the FASTA (e.g. `GRCm39.dict` for `GRCm39.fa`). Uses `VALIDATION_STRINGENCY=LENIENT` because real PE BAMs after MAPQ/flag filtering can trigger mate-field validation warnings (e.g. `INVALID_FLAG_MATE_UNMAPPED`) that would fail the default STRICT mode.
| `tss_enrichment` | `false` | Run deepTools TSS profile QC per treatment sample. Requires `genome_resources.<genome>.gtf`. Produces `results/reference/<genome>.tss.bed`, `<sample>.tss_matrix.gz`, `<sample>.tss_profile.tsv`, and `<sample>.tss_profile.pdf`. |

The original Stage 3 QC switches default to `true`. Heavier optional modules
(`cross_correlation`, `preseq_complexity`, `picard_metrics`,
`tss_enrichment`) default to `false`. Set individual switches explicitly for
production runs.

## MNase-seq fragment ranges

```yaml
mnase:
  mono_range: [140, 200]   # alignmentSieve min/max fragment length for mono-nucleosome BAM
```

The `mnase` block is optional. If absent, `mono_range` defaults to `[140, 200]`.
Only relevant for samples with `assay: mnase`. The dyad BigWig uses `bamCoverage
--MNase` which applies its own default fragment range (130–200 bp).

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
- Only affects samples with `assay: cuttag`. ChIP-seq and ATAC-seq samples
  ignore this setting.
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

## Strict input validation

`scripts/validate_samples.py --strict-inputs` enables optional file-existence
checks for FASTQ paths and Bowtie2 index prefixes. This is off by default so
placeholder-path dry-run profiles continue to work.

| Mode | FASTQ | Bowtie2 index | Use case |
| :--- | :--- | :--- | :--- |
| Default (non-strict) | Not checked | Not checked | Dry-runs, smoke profiles, placeholder paths |
| `--strict-inputs` | Must exist | Complete `.bt2` or `.bt2l` set must exist | Pre-run validation, release checks |

When `--strict-inputs` is active:
- `fastq_1` must exist as a regular file.
- PE `fastq_2` must exist as a regular file.
- The full Bowtie2 index set (6 files) must exist at the configured prefix,
  in either standard `.bt2` or large-index `.bt2l` format.
- A clear error message lists the first missing file.

This flag does NOT re-validate `genome_resources` paths or `control_bam` —
those are already checked regardless of strict mode.

## Defaults and key dependencies

| Key / block | Default | Takes effect when |
| :--- | :--- | :--- |
| `use_control` | `false` | Always; enables control resolution when `true` |
| `multiqc` | `true` | Always; enables MultiQC aggregation when `true` |
| `qc.*` | legacy QC `true`; heavier optional QC `false` | Always; individual switches gate each metric |
| `stage4b` | `true` | Always; enables pooled outputs for multi-biorep experiments |
| `stage5` | `false` | `stage4b: true` + chipseq + narrow + exactly 2 bioreps |
| `idr.*` | listed above | Only when `stage5: true` |
| `cuttag.seacr.enabled` | `false` | Always; affects CUT&Tag samples only |
| `qc.tss_enrichment` | `false` | Requires `genome_resources.<genome>.gtf` |
| `tool_parameters.*` | tool defaults | Always; absent keys use built-in tool defaults |
