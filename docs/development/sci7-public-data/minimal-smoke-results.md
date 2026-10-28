# SCI-7 minimal real-data smoke results

Status: **execution smoke only**. One CUT&Tag Rep2 treatment/IgG default arm
and one isolated MACS3 keep-dup comparison. This is not biological acceptance
and no scientific default is frozen.

## Inputs and reference

- Treatment: SRR12842272 R1/R2, IgG: SRR12842282 R1/R2 (4 FASTQ).
- All four FASTQ and the GRCh38 FASTA matched the official byte counts and MD5.
  Total new download: 1,228,325,160 bytes.
- Local SHA256 recorded in the task evidence directory.
- Reference: GENCODE GRCh38 primary assembly; 194 contigs; chrM present.
- Derived FAI/chrom_sizes and six Bowtie2 index files built from that FASTA;
  194 inspected contigs match the FAI. Reference build: 23:55 wall,
  7,305,712 KiB max RSS, 7.7 GiB reference directory.
- Blacklist: ENCFF356LFX, 8,211-byte gzip MD5
  `393688b4f06c9ce26165d47433dd8c37`; decompressed to 910 BED intervals.

## Production run

- Runtime source identity recorded by `.helixweave/analysis_identity.json`:
  `c223dd9283d83c25a6f0eeb911c9e9c50a45b6a2`, schema 2, policy revision 6,
  identity `7e5599eebc879bce37fb2147b148f75996b8637bf5e88d846d84de3afd49cbe9`.
  The current
  head `b1fcdd6` differs from that commit only by this results document; no
  `src/`, `scripts/`, or `workflow/` file changed between the run and this
  record.
- Actual launch was **not** the canonical bootstrap: the run script invoked
  `/home/yangzichen/miniconda3/envs/chipseq/bin/snakemake` directly with
  `--workflow-profile none` and explicitly set
  `PYTHONPATH=/tmp/sci7-work/checkout/src`. The canonical bootstrap exposes
  `verify-checkout`, `installed-artifact`, `pytest`, `openapi`, `validate`, and
  `local-platform`; it has no `snakemake` command. The run therefore used the
  checkout source path directly, but did not obtain a bootstrap provenance
  check for the workflow process itself. A later `verify-checkout` at the
  current head passes; this is a post-run verification, not evidence that the
  original launch passed bootstrap.
- Tool identity recorded in the run log: Snakemake 8.30.0, samtools 1.23.1,
  MACS3 3.0.4. The historical log did **not** print the Bowtie2, Picard,
  FastQC, Trim Galore, cutadapt, or bamCoverage versions. Current task-local
  binaries and the explicit locks are: Bowtie2 2.5.5 (`align.lock`), Picard
  3.5.0 / OpenJDK 17.0.18 (`duplicate.lock`), FastQC 0.12.1 (`fastqc.lock`),
  Trim Galore 2.2.0 / cutadapt 5.2, and bamCoverage 3.5.6 (`deeptools.lock`).
  Those current/lock identities are not a retroactive claim about the
  historical process.
- `analysis_profile=updated`, CUT&Tag PE, treatment `T2`, control `IgG2`.
- `duplicate_analysis=auto -> retain`; `max_insert=auto -> 700`;
  `blacklist_mode=analysis`; `organelle_mode=auto -> exclude` with `chrM`.
- 48/48 Snakemake steps completed, exit 0.
- Wall time 19:11.70; user CPU 4784.74 s; system CPU 183.61 s;
  431% average CPU; `/usr/bin/time` max RSS 3,483,136 KiB. This is the maximum
  RSS reported for the timed Snakemake process, not the simultaneous total of
  all child processes.
- Outputs include final/analysis BAMs, MAPQ30 BAM, peaks, CPM and fragment-CPM
  bigWigs, FE/ppois bedGraphs, FRiP, NRF/PBC, library complexity, fragment-size
  and mito QC, policy, sample QC summary and result manifest.
- Reported totals: 8,289,326 reads; 360,552 reads in peaks; FRiP 0.043496;
  12,458 narrow peaks (raw = delivered = 12,458); organelle reads 3,000 /
  8,295,528 (0.000362); `analysis_valid_templates=4,144,663`,
  `invalid=0`; blacklist dropped 1,601 templates / 3,202 alignments;
  upstream dropped 3,101 templates / 6,202 alignments.
- Disk: final reference directory 7.7 GiB, final results directory 3.7 GiB,
  inputs about 1.2 GiB. These are final occupancy values, not measured scratch
  high-water marks. No simultaneous scratch peak was captured; the start/end
  `df` samples include other tasks on the host and are not attributable to
  this run.

## MACS3 keep-dup comparison

Same frozen analysis treatment/control BAMs, MACS3 3.0.4, same parameters
(`-f BAMPE -g hs -q 0.01 -B`), only `--keep-dup` and output directory differ.
Input BAM and peak-file SHA256 are stored with the task evidence.

| metric | value |
| --- | ---: |
| `--keep-dup all` peaks | 12,458 |
| `--keep-dup 1` peaks | 2,885 |
| all → 1, at least one overlapping bp | 2,736 |
| 1 → all, at least one overlapping bp | 2,882 |
| exact chr/start/end intersection | 23 |
| exact-coordinate union | 15,320 |
| bp intersection (sorted overlap merged) | 1,892,369 |
| bp union (both peak sets merged) | 8,985,542 |
| bp Jaccard (intersection / union) | 0.210602 |

The earlier summary used “overlap” for all → 1 and “union” for the
exact-coordinate union; those are different definitions from bp overlap and
must not be mixed into a Jaccard value. Peak splitting/merging changes peak
counts, so these numbers describe this MACS3 run, not biological superiority.

The direct command log records the applied `--keep-dup` value in both arms.
The parameter materially changes the peak set on this real library; the
production policy uses `--keep-dup all` as previously specified. This is a
mechanism check, not a decision to change the default.

## Limits

- One treatment/control pair, one default arm, no Rep1, no confirmation cohort,
  no duplicate-exclude arm, no second max_insert arm.
- No biological replicate-concordance or default-freezing conclusion.
- Raw FASTQ/BAM/reference/index are task-local and not committed to Git.
