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

- Checkout `c223dd9283d83c25a6f0eeb911c9e9c50a45b6a2`, canonical bootstrap.
- `analysis_profile=updated`, CUT&Tag PE, treatment `T2`, control `IgG2`.
- `duplicate_analysis=auto -> retain`; `max_insert=auto -> 700`;
  `blacklist_mode=analysis`; `organelle_mode=auto -> exclude` with `chrM`.
- 48/48 Snakemake steps completed, exit 0.
- Wall time 19:11.70; user CPU 4784.74 s; system CPU 183.61 s;
  431% average CPU; max RSS 3,483,136 KiB; results directory 3.7 GiB.
- Outputs include final/analysis BAMs, MAPQ30 BAM, peaks, CPM and fragment-CPM
  bigWigs, FE/ppois bedGraphs, FRiP, NRF/PBC, library complexity, fragment-size
  and mito QC, policy, sample QC summary and result manifest.
- Reported totals: 8,289,326 reads; 360,552 reads in peaks; FRiP 0.043496;
  12,458 narrow peaks (raw = delivered = 12,458); organelle reads 3,000 /
  8,295,528 (0.000362); `analysis_valid_templates=4,144,663`,
  `invalid=0`; blacklist dropped 1,601 templates / 3,202 alignments;
  upstream dropped 3,101 templates / 6,202 alignments.
- `analysis_identity.json`: schema 2, policy revision 6, identity
  `7e5599eebc879bce37fb2147b148f75996b8637bf5e88d846d84de3afd49cbe9`.

## MACS3 keep-dup comparison

Same frozen analysis treatment/control BAMs, MACS3 3.0.4, same parameters
(`-f BAMPE -g hs -q 0.01 -B`), only `--keep-dup` differs:

| arm | peaks |
| --- | ---: |
| `--keep-dup all` | 12,458 |
| `--keep-dup 1` | 2,885 |
| overlap (all as query) | 2,736 |
| union | 15,320 |

The direct command log records the applied `--keep-dup` value in both arms.
The parameter materially changes the peak set on this real library; the
production policy uses `--keep-dup all` as previously specified. This is a
mechanism check, not a decision to change the default.

## Limits

- One treatment/control pair, one default arm, no Rep1, no confirmation cohort,
  no duplicate-exclude arm, no second max_insert arm.
- No biological replicate-concordance or default-freezing conclusion.
- Raw FASTQ/BAM/reference/index are task-local and not committed to Git.
