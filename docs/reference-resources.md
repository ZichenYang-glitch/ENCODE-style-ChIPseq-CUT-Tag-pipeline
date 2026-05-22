# Reference Resources

This guide explains the reference files used by the workflow, what each file
is for, and how to prepare common genome assets without hardcoding local
machine paths into the repository.

For the config keys that point to these files, see
[docs/configuration.md](configuration.md).

## Resource Map

| Resource | Config / sample key | Required? | Used by |
| :--- | :--- | :--- | :--- |
| Bowtie2 index prefix | `samples.tsv: bowtie2_index` | Yes | Bowtie2 alignment |
| Effective genome size | `genome_resources.<genome>.effective_genome_size` | Yes | MACS3 peak calling and deepTools signal normalization |
| Reference FASTA | `genome_resources.<genome>.reference_fasta` | Required when `qc.picard_metrics: true` | Picard CollectMultipleMetrics |
| FASTA index | `<reference_fasta>.fai` | Required with Picard metrics | Picard / htsjdk reference access |
| Sequence dictionary | `<reference_fasta basename>.dict` | Required with Picard metrics | Picard / htsjdk sequence metadata |
| Chrom sizes | `genome_resources.<genome>.chrom_sizes` | Optional now | Future bedGraph-to-BigWig and genome-window operations |
| Blacklist BED | `genome_resources.<genome>.blacklist` | Required when blacklist filtering is enabled for that genome | BAM and peak blacklist filtering |
| GTF / annotation | `genome_resources.<genome>.gtf` | Optional now | Future TSS enrichment / gene annotation modules |

## Path Rules

- Repository examples should use portable placeholder paths.
- Real local runs may use absolute paths in local config files.
- Do not commit local FASTA, index, BAM, FASTQ, BigWig, or generated result
  files.
- Keep all genome resource files outside the repository or under an ignored
  local run directory.
- The `genome` column in `samples.tsv` must match a key under
  `genome_resources`.

## Bowtie2 Index Prefix

The sample sheet `bowtie2_index` value is a prefix, not a directory.

If these files exist:

```text
/refs/mm39/bowtie2/GRCm39.1.bt2
/refs/mm39/bowtie2/GRCm39.2.bt2
/refs/mm39/bowtie2/GRCm39.3.bt2
/refs/mm39/bowtie2/GRCm39.4.bt2
/refs/mm39/bowtie2/GRCm39.rev.1.bt2
/refs/mm39/bowtie2/GRCm39.rev.2.bt2
```

then `bowtie2_index` should be:

```text
/refs/mm39/bowtie2/GRCm39
```

Build an index from a FASTA:

```bash
mkdir -p /refs/mm39/bowtie2
bowtie2-build /refs/mm39/GRCm39.fa /refs/mm39/bowtie2/GRCm39
```

## Effective Genome Size

`effective_genome_size` is the mappable genome size used by MACS3 and related
signal calculations. It is not the same thing as a per-chromosome sizes file.

Accepted forms:

```yaml
effective_genome_size: "hs"       # MACS3 human shortcut
effective_genome_size: "mm"       # MACS3 mouse shortcut
effective_genome_size: 2654621783 # explicit integer, useful for mm39/GRCm39
```

Use the shortcuts only when they match the intended assembly. For mm39/GRCm39,
prefer an explicit integer rather than the `mm` shortcut.

## Reference FASTA, `.fai`, and `.dict`

Picard metrics require:

```text
GRCm39.fa
GRCm39.fa.fai
GRCm39.dict
```

The `.fai` file indexes FASTA sequence offsets and lengths. The `.dict` file
stores Picard/htsjdk sequence dictionary metadata. They overlap in information
but are not interchangeable.

Prepare both files:

```bash
samtools faidx /refs/mm39/GRCm39.fa
picard CreateSequenceDictionary \
  R=/refs/mm39/GRCm39.fa \
  O=/refs/mm39/GRCm39.dict
```

If `picard` is not on PATH but `samtools` supports dictionary creation:

```bash
samtools dict /refs/mm39/GRCm39.fa > /refs/mm39/GRCm39.dict
```

Then configure:

```yaml
genome_resources:
  mm39:
    reference_fasta: "/refs/mm39/GRCm39.fa"
```

## Chrom Sizes

A chrom sizes file is a two-column text file:

```text
chr1    195154279
chr2    181755017
...
```

Create it from a FASTA index:

```bash
samtools faidx /refs/mm39/GRCm39.fa
cut -f1,2 /refs/mm39/GRCm39.fa.fai > /refs/mm39/GRCm39.fa.sizes
```

Configure:

```yaml
genome_resources:
  mm39:
    chrom_sizes: "/refs/mm39/GRCm39.fa.sizes"
```

The workflow does not currently require `chrom_sizes` for default outputs, but
future bedGraph-to-BigWig and TSS enrichment features will use it.

## Blacklist BED

Blacklist files mark problematic regions such as repeats, assembly artifacts,
and anomalous signal regions. Configure them per genome:

```yaml
genome_resources:
  hg38:
    blacklist: "/refs/hg38/hg38.blacklist.bed"
```

When `qc.blacklist_filter: true`, the workflow uses the blacklist path if it
is configured for the sample genome. If a genome has no blacklist available,
set the value to an empty string or disable blacklist filtering for that run.

Use an assembly-matched blacklist only. Do not reuse hg38 blacklists for hg19,
or mm10 blacklists for mm39, unless the file has been lifted over and checked.

## GTF / Annotation

`gtf` is optional in the current workflow. It is reserved for annotation-driven
QC such as TSS enrichment and gene-centric summaries.

Configure when available:

```yaml
genome_resources:
  mm39:
    gtf: "/refs/mm39/gencode.annotation.gtf"
```

The annotation assembly must match the FASTA and Bowtie2 index assembly.

## mm39 / GRCm39 Minimal Example

```yaml
genome_resources:
  mm39:
    effective_genome_size: 2654621783
    chrom_sizes: "/refs/mm39/GRCm39.fa.sizes"
    blacklist: ""
    gtf: "/refs/mm39/gencode.annotation.gtf"
    reference_fasta: "/refs/mm39/GRCm39.fa"
```

Sample sheet:

```tsv
sample	fastq_1	fastq_2	layout	assay	target	peak_mode	genome	bowtie2_index
sample1	/data/sample1_R1.fq.gz	/data/sample1_R2.fq.gz	PE	chipseq	H3K27ac	narrow	mm39	/refs/mm39/bowtie2/GRCm39
```

Before running with Picard metrics:

```bash
test -f /refs/mm39/GRCm39.fa
test -f /refs/mm39/GRCm39.fa.fai
test -f /refs/mm39/GRCm39.dict
test -f /refs/mm39/GRCm39.fa.sizes
ls /refs/mm39/bowtie2/GRCm39*.bt2
```
