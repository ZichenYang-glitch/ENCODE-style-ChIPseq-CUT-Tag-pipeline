# Configuration Reference

`config/config.yaml` controls workflow behavior, resource paths, and feature
gating. This page is the maintained configuration reference; the
[scientific quick start](quickstart.md) provides the shortest setup path.

## Core keys

| Key | Default | Accepted | Effect |
| :--- | :--- | :--- | :--- |
| `samples` | `"config/samples.tsv"` | path | Path to the sample sheet. |
| `outdir` | `"results"` | path | Output root directory. |
| `threads` | `8` | positive integer | Threads per job. |
| `mapq` | `30` | integer >= 0 | Minimum MAPQ for samtools filtering. |
| `binsize` | `10` | positive integer | Bin size in bp for bamCoverage BigWig. |
| `remove_dup` | `"auto"` | `"auto"`, `"yes"`, `"no"` | Legacy alias for the duplicate policy. See below. |
| `duplicate_analysis` | `"auto"` | `"auto"`, `"exclude"`, `"retain"` | Duplicate policy selector; per-sample assay resolution in auto. See below. |
| `trim` | `true` | boolean | Enable Trim Galore adapter/quality trimming and trimmed-read FastQC. Raw-input FastQC always runs. |
| `extend_reads` | `"auto"` | `"auto"`, `"yes"`, `"no"` | Extend SE reads to fragment size. |
| `use_control` | `false` | boolean | Enable control sample resolution. |
| `multiqc` | `true` | boolean | Enable MultiQC report aggregation. |

### Duplicate policy (`duplicate_analysis` and `remove_dup`)

The duplicate policy is decided per sample at the BAM layer; every MACS3
call (per-sample, pooled, consensus, and IDR-ready) then renders
`--keep-dup all` so peak callers consume the chosen molecule set without a
second implicit coordinate collapse.

| `duplicate_analysis` | Behavior |
| :--- | :--- |
| `"auto"` (default) | Resolve per sample by assay: CUT&Tag **retains** duplicates (apparent duplicates at high-affinity sites are mostly real independent fragments); ChIP-seq, ATAC-seq, and MNase-seq **exclude** them. The narrow/broad peak mode no longer changes this decision — legacy broad ChIP-seq now excludes duplicates like every other exclude assay. |
| `"exclude"` | Remove duplicates for all samples. |
| `"retain"` | Keep duplicates for all samples (marked, for metrics, where the backend supports marking). |

`remove_dup` remains as a legacy alias: `"yes"` maps to `exclude`, `"no"`
maps to `retain`, and `"auto"` never overrides an explicit
`duplicate_analysis`. Setting both to conflicting values
(e.g. `remove_dup: "no"` with `duplicate_analysis: exclude`) is rejected.
The startup summary and `<outdir>/.helixweave/effective_config.json`
record the per-sample resolution (`resolved_duplicate_analysis`).

In legacy runs, exclusion uses Picard `REMOVE_DUPLICATES=true` when Picard
is available, otherwise the `samtools markdup` fallback; retained samples
keep a symlinked `final.bam`. In updated runs, Picard marks duplicates in
`final.bam` and `analysis.bam` applies the policy (`-F 0x400` for exclude
samples).

A legacy `--keep-dup` in `tool_parameters.macs3.extra_args` or
`tool_parameters.idr_macs3.extra_args` is migrated: a synonymous
`--keep-dup all` is normalized away (the workflow renders the flag
itself); any other value is rejected with migration guidance rather than
being decided by argument order.

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
samples without a control (`-c` omitted). If the sample sheet still contains a
control row or control reference, validation emits a `use_control=False`
warning; it does not change behavior.

When `true`, treatment rows may reference a control via `control_sample`
(another sample row) or `control_bam` (an external BAM path). Control rows
(`role: control`) are preprocessed but skip peak calling. Control-reference
cycles are rejected.

## `multiqc`

**Default:** `true`

When `true`, MultiQC aggregates QC artifacts from all active samples into
`results/multiqc/multiqc_report.html`.

FastQC is reported in native `FastQC (raw input)` and, when `trim: true`,
`FastQC (trimmed reads)` sections. PE samples keep R1/R2 identities and SE
samples have R1 only. `trim: false` does not rescan the alignment-input
symlink. No filtered FastQC section is created because filtering produces BAM,
not FASTQ; native `Samtools (aligned, pre-filter)` and `Samtools (final
alignment)` sections report that transition instead. These behaviors require
no additional configuration keys.

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
  cross_correlation: false       # opt-in
  preseq_complexity: false       # opt-in
  picard_metrics: false          # opt-in
  tss_enrichment: false          # opt-in; requires GTF
```

| Switch | Default | Effect |
| :--- | :--- | :--- |
| `blacklist_filter` | `true` | Create additional blacklist-filtered BAM and peak copies for eligible treatment samples with a configured blacklist BED. Raw outputs remain available; consumers are described below. |
| `frip` | `true` | Compute read-record FRiP for peak-assay treatment samples, using both filtered copies when blacklist QC applies. |
| `library_complexity` | `true` | Parse duplicate metrics. Picard fields are reported when Picard is available; otherwise fallback fields are emitted as `NA`. |
| `nrf_pbc` | `true` | Compute NRF and PBC1/PBC2 from the MAPQ-filtered, duplicate-containing BAM (`{sample}.mapq{mapq}.bam`), before duplicate handling. |
| `signal_tracks` | `true` | Produce MACS3 FE (fold-enrichment) and ppois (Poisson p-value) bedGraph tracks for eligible peak-assay treatment samples and pooled replicate experiments. FE/ppois BigWig conversion is available when `genome_resources.<genome>.chrom_sizes` is also configured. |
| `summary` | `true` | Emit per-sample and project QC summaries for peak-assay treatment samples. The project TSV is `results/multiqc/project_qc_summary.tsv`; MNase has its own summary, independent of this switch. |
| `cuttag_fragment_size` | `true` | Compute CUT&Tag fragment-size statistics for active samples with assay=cuttag. |
| `cross_correlation` | `false` | Run phantompeakqualtools cross-correlation QC per treatment sample. Produces NSC/RSC metrics (`.cc.qc`) and a cross-correlation plot (`.cc.plot.pdf`). When enabled, also generates a project-level summary at `results/multiqc/cross_correlation_summary.tsv` and makes it available to the custom section when MultiQC is enabled. See [docs/qc-interpretation.md](qc-interpretation.md) for interpretation guidance. |
| `preseq_complexity` | `false` | Run preseq on the same MAPQ-filtered, duplicate-containing BAM per treatment sample. Uses `lc_extrap -B -P` for PE and `lc_extrap -B` for SE, producing `.preseq.txt` when extrapolation succeeds. |
| `picard_metrics` | `false` | Run Picard CollectMultipleMetrics per treatment sample. Produces alignment summary, insert size, and quality distribution metrics. Requires `genome_resources.<genome>.reference_fasta` with a matching samtools FASTA index (`.fai`) and Picard sequence dictionary (`.dict`) next to the FASTA (e.g. `GRCm39.dict` for `GRCm39.fa`). Uses `VALIDATION_STRINGENCY=LENIENT` because real PE BAMs after MAPQ/flag filtering can trigger mate-field validation warnings (e.g. `INVALID_FLAG_MATE_UNMAPPED`) that would fail the default STRICT mode.
| `tss_enrichment` | `false` | Run deepTools TSS profile QC per treatment sample. Requires `genome_resources.<genome>.gtf`. Produces `results/reference/<genome>.tss.bed`, `<sample>.tss_matrix.gz`, `<sample>.tss_profile.tsv`, and `<sample>.tss_profile.pdf`. |

These 11 defaults come from [`config/qc.py`](../src/encode_pipeline/config/qc.py):
seven are `true` and four are `false`. Heavier optional modules
(`cross_correlation`, `preseq_complexity`, `picard_metrics`,
`tss_enrichment`) default to `false`. Set individual switches explicitly for
production runs. Default target selection also depends on assay, role and resources;
for example, a peak summary pulls its component metrics even if their individual
switches are false. Explicitly requesting a rule output is distinct from default
target selection.

Blacklist handling is not a global replacement of BAMs, peaks or signal tracks.
In [`qc.smk`](../workflow/rules/qc.smk), `_has_blacklist_qc` requires both the
switch and a sample in `BLACKLIST_SAMPLES` (treatment samples with a nonempty
blacklist genome resource). `_frip_inputs` then requests the filtered BAM **and**
filtered peaks; otherwise it requests `final.bam` and the original peak directory.
This is configuration-based DAG selection: missing selected outputs must be
produced, not silently replaced with raw inputs. `_peak_counts_inputs` requests
raw peaks and optional filtered peaks, never the filtered BAM.

### Analysis-filtered blacklist branch

`analysis_profile: updated` also supports a whole-analysis blacklist branch:

```yaml
analysis_profile: updated
blacklist_mode: analysis
genome_resources:
  hg38:
    effective_genome_size: "hs"
    blacklist: "/refs/hg38/hg38.blacklist.bed"
```

`blacklist_mode: analysis` requires a non-empty blacklist for every active
sample genome and is rejected under `analysis_profile: legacy`. Raw
preprocessing outputs stay under `results/<sample>/`; analysis outputs move to
`results/analysis_filtered/<sample>/` and
`results/analysis_filtered/experiments/<experiment>/`. The analysis BAM is
derived after duplicate/organelle filtering. PE templates are removed when the
fragment span (`min(POS, PNEXT)` through `min + abs(TLEN) - 1`) overlaps the
blacklist; SE reads use their alignment span. Non-proper/orphan PE records fail
closed. Per-template drop counts are written to
`05_qc/<sample>.blacklist_dropped.tsv`.

In this branch, MACS3, BigWig, FRiP, SEACR, cross-correlation, TSS, MNase, IDR
and consensus consume the filtered analysis BAM/namespace. Per-sample peaks are
also written as `04_peaks/<sample>_blacklist_filtered/`, and final IDR/consensus
peak sets are filtered before promotion. The raw `final.bam` and raw QC remain
available under `results/`.

`final.bam` may retain duplicates according to the duplicate policy.
Default CPM coverage uses that final BAM in
[`common.smk`](../workflow/rules/common.smk); MACS3 FE/ppois tracks use MACS3
pileup/background bedGraphs, not the optional blacklist copies. See
[QC interpretation](qc-interpretation.md#frip-fraction-of-reads-in-peaks) for
counting and signal-consumer details. A workspace TSV does not by itself prove
platform listing/download or complete MultiQC HTML visibility.

MultiQC and reproducibility have separate defaults in
[`validator.py`](../src/encode_pipeline/config/validator.py) and
[`reproducibility.py`](../src/encode_pipeline/config/reproducibility.py):
`multiqc=true`, `replicate_analysis=true`, `chipseq_idr=false`, and
`reproducibility.enabled=false`. The consensus child default (`true`) does not
bypass its parent or replicate/assay eligibility. See the
[reproducibility policy](reproducibility-policy.md#5-config-reference).

## Platform reference binding v2

Platform reference profiles may use `encode-reference-binding-v2`. In addition
to the v1 content hashes, v2 requires a `reference_fai` descriptor and verifies
that FASTA headers/lengths, `.fai`, `chrom_sizes`, GTF contigs, blacklist
contigs, and an optional `reference_dict` agree. It also accepts
`organelle_chroms` and `canonical_chroms`; every listed contig must exist in the
FASTA. Symbolic effective genome sizes are constrained to matching assembly
families (`hs` for GRCh/hg, `mm` for GRCm/mm); numeric sizes remain explicit.

`organelle_chroms` is injected into `genome_resources.<assembly>` for the
workflow. Binding v1 remains accepted and keeps its existing identity. The
binding identity includes the v2 contract and metadata and remains path-free.
FAI verification now checks contig names, lengths, offsets, line bases, and
line widths against the same hash-verified FASTA byte stream, so random-access
usability is part of v2 admission without a second unverified file read.
Bowtie2 index contents remain hash-bound without internal contig parsing.

## Fragment-level CPM tracks

```yaml
analysis_profile: updated
fragment_cpm: auto       # off | auto | on; legacy requires off
```

`fragment_cpm` adds an additive `.fragment.CPM.bw` track without changing the
existing `.CPM.bw` path or read-level semantics. In `auto` mode, tracks are
emitted for ATAC, MNase, and CUT&Tag targets classified as histone
(broad-like, narrow-like, or context-dependent). CUT&Tag targets the classifier
cannot recognize are skipped; explicit `on` also emits them. ChIP-seq remains
read-level in this release.

Non-MNase tracks live at
`results/<sample>/03_bigwig/<sample>.fragment.CPM.bw`; MNase mono-nucleosome
tracks live at
`results/<sample>/04_signal/<sample>.mono.fragment.CPM.bw`. In the
`analysis_filtered` branch the same paths are rooted under
`results/analysis_filtered/`.

PE tracks use `bamCoverage --extendReads --samFlagInclude 64`, counting one
first-mate fragment per pair; SE tracks use `--extendReads` without the flag.
MNase dyad tracks keep their existing `--MNase --binSize 1` behavior.

## Fragment-occupancy TSS profile

```yaml
analysis_profile: updated
fragment_cpm: auto
tss_fragment_profile: true
```

The opt-in fragment TSS profile uses the fragment CPM BigWig and requires a
GTF for every eligible ATAC/CUT&Tag sample. It writes
`05_qc/tss/<sample>.fragment.tss_matrix.gz`, `.fragment.tss_profile.tsv`, and
`.fragment.tss_profile.pdf`. This is a descriptive fragment-occupancy profile,
not an ENCODE TSS enrichment score; the default `qc.tss_enrichment` profile is
unchanged. MNase is not eligible for this profile.

## MNase-seq fragment ranges and caller surface

```yaml
mnase:
  fragments:
    sub: [1, 139]          # alignmentSieve min/max for sub-nucleosome BAM
    mono: [140, 200]       # alignmentSieve min/max for mono-nucleosome BAM
    di: [300, 400]         # alignmentSieve min/max for di-nucleosome BAM
  dyad_range: [130, 200]   # bamCoverage --MNase min/max fragment length
  callers:
    danpos3: false         # reserved; setting true raises validation error
    inps: false            # reserved; setting true raises validation error
    sem: false             # reserved; setting true raises validation error
```

The `mnase` block is optional. If absent, `fragments.mono` defaults to `[140, 200]`.
The legacy key `mono_range` is still accepted but `fragments.mono` takes precedence.
Only relevant for samples with `assay: mnase`. The dyad BigWig uses `bamCoverage
--MNase` with `dyad_range` defaulting to `[130, 200]`.

## Replicate and IDR features

```yaml
replicate_analysis: true        # replicate-aware pooled BAMs and peaks (default on)
chipseq_idr: false        # TF ChIP-seq IDR (default off; requires replicate_analysis: true)

idr:
  seed: 42
  threshold: 0.05    # IDR -log10(threshold) for thresholded peak set
  rank: "p.value"    # ranking measure: p.value or signal.value
```

- **`replicate_analysis`** (default `true`): Enables technical replicate merging,
  biological-replicate BAMs, pooled treatment/control BAMs, pooled MACS3 peak
  calls, and pooled signal tracks. Applies to experiments with 2+ biological
  replicates.
- **`chipseq_idr`** (default `false`): Enables TF ChIP-seq IDR. Requires
  `replicate_analysis: true`, `chipseq` assay, `narrow` peak_mode, and exactly 2
  treatment biological replicates per experiment.
- **`idr.*`**: The idr block is only read when `chipseq_idr: true`. `threshold` is
  the IDR threshold (0.05 = IDR -log10(0.05) ≈ 1.3). `rank` chooses between
  p-value and signal-value ranking for IDR input.
- IDR rules use `workflow/envs/idr.yml` when Snakemake is run with
  `--use-conda`, keeping IDR's older Python dependency constraints out of the
  core environment.

### ChIP-seq narrow IDR (`chipseq_idr`) gating summary

| Condition | Required |
| :--- | :--- |
| `replicate_analysis: true` | Always required for IDR |
| `assay` | `chipseq` only |
| `peak_mode` | `narrow` only |
| Biological replicates | Exactly 2 treatment bioreps per experiment |

This ChIP-seq narrow IDR path does not support CUT&Tag or 3+ replicate IDR. The separate
`reproducibility` block provides opt-in CUT&Tag narrow IDR and the maintained
consensus policies described in the
[reproducibility policy](reproducibility-policy.md).

## CUT&Tag SEACR

```yaml
cuttag:
  seacr:
    enabled: false
```

- **`cuttag.seacr.enabled`** (default `false`): When `true`, runs SEACR peak
  calling as a sidecar alongside MACS3 for CUT&Tag samples. Outputs land under
  `results/<sample>/04_peaks_seacr/`. If no PE CUT&Tag treatment sample is
  eligible, the workflow emits a warning and produces no SEACR output.
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
| `macs3` | `qvalue` | `0.01` | MACS3 q-value cutoff for peak calling. Must be a finite probability in (0, 1]. |
| `macs3` | `broad_cutoff` | `0.1` | MACS3 broad-peak q-value cutoff. Must be a finite probability in (0, 1]. |
| `idr_macs3` | `pvalue` | `0.1` | Relaxed p-value for IDR-ready per-biorep peak calls. Only used when `chipseq_idr: true`. Must be a finite probability in (0, 1]. |
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
| `replicate_analysis` | `true` | Always; enables pooled outputs for multi-biorep experiments |
| `chipseq_idr` | `false` | `replicate_analysis: true` + chipseq + narrow + exactly 2 bioreps |
| `idr.*` | listed above | Only when `chipseq_idr: true` |
| `cuttag.seacr.enabled` | `false` | Always; affects CUT&Tag samples only |
| `qc.tss_enrichment` | `false` | Requires `genome_resources.<genome>.gtf` |
| `tool_parameters.*` | tool defaults | Always; absent keys use built-in tool defaults |

## Run identity and effective configuration

Every run — both `legacy` and `updated` profiles — records a versioned
identity marker at `<outdir>/.helixweave/analysis_identity.json` before any
scientific task executes. The marker carries a marker schema version, the
analysis profile, the internal `scientific_policy_revision`, and a hash of the
merged, validated configuration plus the parsed sample sheet. The identity
covers every validated (effective) config value — including the `mnase`,
`cuttag`/SEACR, `qc`, `idr`/`reproducibility`, genome-resource, duplicate,
organelle, blacklist, and tool-parameter blocks — with only these explicitly
excluded scheduling/location keys: `outdir`, `threads`, and the `samples`
path (the parsed sample content is hashed instead).
`scientific_policy_revision` is an internal workflow identity; it cannot be
set by user configuration. Path-valued inputs (FASTQ, Bowtie2 index prefix,
external control BAM, and genome-resource file paths) enter the identity
through their canonical absolute locations, so the same relative string in a
different working directory cannot alias a different actual input; the same
actual input referenced through equivalent path forms stays equivalent.

| Outdir / marker state | Behavior |
| :--- | :--- |
| New directory, or only scheduler metadata (`.snakemake`, `.helixweave`) | The marker is created, then scientific tasks run |
| Marker matches profile, policy revision, and identity | The same run may continue |
| Scientific outputs present but no marker | Refused; choose a fresh outdir (existing content is kept) |
| Old-format marker, other profile, other policy revision, or different identity | Refused; choose a fresh outdir |
| Marker or record corrupted, not a regular file, or a symlink; `.helixweave` a symlink or not a directory | Refused; never auto-repaired or overwritten |
| Marker valid but `effective_config.json` missing, unreadable, missing required fields, or inconsistent with the run identity | Refused; the run does not continue without an intact, truthful configuration record |

The record stores scientific input paths in the canonical form resolved at
creation time, so resuming from a different working directory with
equivalent absolute paths of the same actual inputs stays valid, while a
different actual resource is refused. The record is validated beyond file
integrity: it must carry the profile, schema version (consistent with the
recorded configuration), policy revision, identity, and a `config` block
whose canonical form recomputes to the recorded identity, so a bare
identity string or a tampered inner configuration is refused like any other
mismatch. Nothing is auto-repaired from current parameters.

The marker and the effective-config record form one metadata boundary: each
is written to an exclusively created temporary regular file inside the
verified metadata directory and atomically replaced into place, record
first, so a valid marker never publishes without its record and neither
symlinked parents nor pre-placed temporary symlinks can redirect the write.

Dry-runs perform the same checks (a dry-run reports a mismatch) but never
create the marker or any run record.

### CLI overrides and the effective configuration record

Snakemake `--config KEY=VALUE` overrides are merged with `--configfile` before
validation, so an override goes through the same validator, the run identity,
and the effective configuration record. The merged value must be legal for the
active profile; conflicting or internal keys are rejected.

```bash
snakemake -s workflow/Snakefile --configfile config/config.yaml \
  --config mapq=20
```

At startup the workflow prints a short, path-free summary of the effective
configuration (profile, policy revision, mapq/trim/binsize, duplicate,
organelle, and blacklist modes, sample assay counts). An actual run also
saves the merged, normalized configuration to
`<outdir>/.helixweave/effective_config.json`. A dry-run prints the summary but
writes nothing into the run directory.

Overrides that change scientific identity require a fresh outdir — they are
not applied on top of an existing run directory. The standalone
`validate --config FILE` command validates only that file; it does not
validate `--config` values injected later. To check an overridden run, use
the workflow's own merged-config validation and dry-run.
