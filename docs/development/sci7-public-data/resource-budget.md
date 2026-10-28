# SCI-7 resource budget (from official ENCODE metadata)

## Selected official files (2026-10-28)

| Group | Files | Official bytes | GiB |
|---|---:|---:|---:|
| TF ChIP-seq ENCSR054ZMK treatment | 2 PE | 5,972,185,158 | 5.56 |
| TF matched control ENCSR142JFU | 2 PE | 5,424,758,870 | 5.05 |
| Histone ENCSR817LUF treatment | 5 SE | 1,716,109,509 | 1.60 |
| Histone control ENCSR603VTC | 2 SE | 794,586,045 | 0.74 |
| ATAC ENCSR844TIU (2 reps) | 4 PE | 3,952,560,652 | 3.68 |
| **Total compressed FASTQ** | **15** | **17,860,200,234** | **16.63 GiB** |

No MNase FASTQ or CUT&Tag file is included: ENCODE MNase is SOLiD-only and the
CUT&Tag API path returned no experiment. These are gaps, not zero-cost data.

## Reference/index (identity to confirm; not downloaded)
- TF/histone: GRCh38 primary assembly, matching GTF, blacklist, Bowtie2 index.
- ATAC: mm10 primary assembly, matching GTF, blacklist, Bowtie2 index.
Estimate 3-5 GiB for the GRCh38 Bowtie2 index plus annotation/blacklist; mm10
similar. Do not count these as FASTQ.

## Storage peak (estimate, to be measured)
- Raw FASTQ: 16.63 GiB.
- Per-branch BAM/index/intermediates: 1.5-3x compressed FASTQ depending on
  compression and duplicate handling; peak planning range 60-100 GiB for the
  full selected set, plus 10-20 GiB for peaks/logs/temporary sort files.
- Keep raw FASTQ until checksums and a minimal trial pass; results retention
  policy is separate.

## CPU/RAM/wall time (estimate, to be measured)
- Alignment/QC/peak calling: roughly 20-60 CPU-hours for the full set,
  assuming 1-3 CPU-hours per 10M PE reads and small peak-calling cost.
- RAM: plan 16 GiB per sample, 32 GiB peak for sort/index; Bowtie2 GRCh38
  index must fit in RAM.
- Wall time: 4-8 hours on 8-16 cores for the full set; the minimal trial is
  much shorter. These are estimates, not measured values.

## Run order
1. Minimal trial: histone treatment+control (2.34 GiB) end-to-end.
2. Add ATAC (3.68 GiB) for PE/max_insert behavior.
3. Add TF (10.61 GiB) for the TF/control duplicate-policy comparison.
4. Only after verified checksums, reference identity, and a written decision
   rule; no bulk download before step 1 passes.
