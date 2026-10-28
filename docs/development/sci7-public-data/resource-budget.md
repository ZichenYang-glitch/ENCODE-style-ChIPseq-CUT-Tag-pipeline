# SCI-7 resource budget

Status: planning values for the final registered candidates on 2026-10-28. No
FASTQ/BAM, large FASTA/GTF or index was downloaded. “Official” values come from
portal metadata or HTTPS HEAD; all storage multipliers, CPU/RAM and wall times
are estimates to be measured, not official figures.

## 1. Official input download volume

### Discovery

| Group | FASTQ | Official bytes | GiB |
| --- | ---: | ---: | ---: |
| TF ENCSR054ZMK treatment + ENCSR142JFU control | 4 | 11,396,944,028 | 10.613 |
| Histone ENCSR817LUF + ENCSR603VTC | 7 | 2,510,695,554 | 2.338 |
| ATAC ENCSR844TIU | 4 | 3,952,560,652 | 3.681 |
| CUT&Tag GSE159608 / PRJNA669859 treatment + IgG | 8 | 443,578,936 | 0.413 |
| MNase PRJNA928451 LP_pMSCs_FM | 6 | 4,402,061,533 | 4.100 |
| **Total** | **29** | **22,705,840,703** | **21.146** |

### Confirmation candidates

| Group | FASTQ | Official bytes | GiB |
| --- | ---: | ---: | ---: |
| TF ENCSR785RQR + ENCSR194GTV | 4 | 74,110,406,175 | 69.021 |
| Histone ENCSR213YKJ rep1 + ENCSR148TRG rep1 | 16 | 7,462,062,697 | 6.949 |
| ATAC ENCSR639PPN | 8 | 11,262,765,629 | 10.489 |
| Optional CUT&Tag GSE145187 / PRJNA606273 | 6 | 232,017,631 | 0.216 |
| Optional MNase PRJNA928451 HP_pMSCs_FM | 4 | 3,139,436,756 | 2.924 |
| **Total planned, not downloaded** | **38** | **96,206,688,888** | **89.599** |

The discovery plus all registered confirmation inputs total **67 FASTQ,
118,912,529,591 bytes (110.746 GiB)**. Optional candidates are included so the
budget is not understated; their execution still depends on the independence
and depth limitations in `confirmation-review.md`.

## 2. Reference and index

Official HTTPS HEAD sizes for compressed primary-assembly FASTA are 844,691,642
bytes (GRCh38) and 773,507,376 bytes (mm10), totalling **1,618,199,018 bytes
(1.507 GiB)**. Optional GTF HEAD sizes total 78,287,254 bytes (0.073 GiB), but
GTF is not required by the default comparison. See
`metadata/confirmation/reference-head-checks.json`.

Estimated reference-side working storage is **15–20 GiB for both species**:
compressed plus decompressed FASTA, FAI/chrom_sizes, small blacklists and two
Bowtie2 2.5.5 indexes. This is an estimate because decompressed FASTA and index
sizes have not been measured. The two ENCODE blacklist payloads (9,474 bytes)
are already verified. Reuse the same validated references/indexes across arms;
do not rebuild them per run.

## 3. Unique scientific configurations and intermediate volume

`auto` resolution checks add no duplicate scientific arm. A baseline default
configuration is shared by the duplicate and max_insert comparisons. Thus each
PE dataset needs three core full-chain configurations (default, alternative
duplicate policy, alternative max_insert); histone SE needs two (max_insert not
applicable). MNase 200 is included as a sensitivity diagnostic, not a proposed
default.

| Cohort | PE datasets x 3 | SE datasets x 2 | Core configurations |
| --- | ---: | ---: | ---: |
| Discovery | TF, ATAC, CUT&Tag, MNase | Histone | 14 |
| Confirmation | TF, ATAC, CUT&Tag, MNase | Histone | 14 |
| **Total** |  |  | **28** |

The isolated MACS3 keep-dup experiment reuses frozen BAMs and adds no alignment
BAM. One optional combined regression per dataset would add up to 10 further
configurations and is not included in the core total below.

For rough capacity planning, each full-chain configuration's BAM/index/sort
intermediates are estimated at **1.5–3x its compressed FASTQ input**. The
28-configuration weighted input is 65.607 GB for discovery and 281.155 GB for
confirmation. Therefore:

- discovery intermediates: about 98–197 GB (91–183 GiB);
- confirmation intermediates: about 422–843 GB (393–785 GiB);
- **all core intermediates if everything is retained: 520–1,040 GB
  (484–969 GiB)**.

Final peaks/QC/logs/bigWigs across the matrix are separately estimated at
**20–80 GiB**. Retaining all raw inputs, all branch intermediates, references
and final reports therefore requires roughly **1.2 TiB free**, with additional
headroom for the optional combined regression. This is a conservative planning
figure, not a measured requirement.

## 4. Retention policy and peak-storage control

Do not retain every alternative-arm intermediate by default.

1. Retain raw FASTQ and official/local checksum evidence: 110.746 GiB for all
   registered candidates.
2. Retain shared validated references/indexes: estimated 15–20 GiB.
3. Retain only the selected candidate arm's BAM/index for each dataset:
   estimated 167–334 GiB across all datasets; retain alternative-arm BAMs only
   while their metrics and attribution checks are being computed.
4. Retain final peaks, QC, manifests, commands, logs and effect-size tables for
   every arm: estimated 20–80 GiB.
5. Process one data group at a time, largest TF confirmation last. After the
   comparison evidence is complete and its checksums verified, future
   execution may clean only its own explicitly listed alternative-arm scratch
   outputs; no cleanup is authorized or performed by this preparation step.

The resulting long-term retention estimate is about **313–545 GiB**, plus one
active group's alternative-arm scratch. Processing TF confirmation alone can
temporarily add about **207–414 GiB** for its two non-retained core arms. Keep
at least 1 TiB free for a conservative full plan, or reduce peak usage by
running one group at a time and applying the retention policy. Available space
must be rechecked before each download; this document does not authorize a
download.

## 5. CPU, RAM and wall-time estimates

The previous planning assumption was 20–60 CPU-hours per 17.86 GB of compressed
input for one alignment/QC/peak arm. Scaling that rough assumption by the
28-configuration weighted input gives:

- discovery core matrix: approximately **70–220 CPU-hours**;
- confirmation core matrix: approximately **310–950 CPU-hours**.

These estimates are highly uncertain because read length, PE/SE mix, duplicate
rate and reference differ. The optional combined regressions and repeated
failures are extra. On 16 dedicated cores, the core matrices correspond to a
rough sequential wall-time range of **6–20 hours for discovery** and **20–60
hours for confirmation**, plus reference setup and download time. Do not treat
larger `-X` cost as linear; measure it in the max_insert arms.

Plan **16 GiB RAM per ordinary sample and 32 GiB for sort/index peaks**. For
concurrent samples or reference construction, plan 32–64 GiB total RAM and keep
the Bowtie2 index resident. Record CPU time, wall time and maximum RSS for every
arm and replace these estimates with measurements after the minimal trial.

## 6. Minimal trial and expansion order

1. **Minimal execution smoke test:** CUT&Tag Rep2 treatment + IgG Rep2
   (SRR12842272/SRR12842282), 383,633,518 compressed bytes. It exercises PE,
   matched-control, duplicate and max_insert paths at low cost. It is not
   sufficient for a biological default decision.
2. **Reference/QC expansion:** full CUT&Tag discovery (0.413 GiB), then MNase
   Rep1 (about 1.43 GiB) to check mono/di windows and PE accounting. Resolve
   shallow CUT&Tag Rep1 separately; do not silently drop it from a biological
   comparison.
3. **Remaining discovery:** histone treatment/control (2.338 GiB), ATAC (3.681
   GiB), then TF treatment/control (10.613 GiB). Run all single-variable arms
   before the combined regression.
4. **Confirmation only after discovery rules are fixed:** optional
   CUT&Tag/MNase first, then histone, ATAC, and the 69.021 GiB TF confirmation
   last. Do not tune on confirmation data.

Each expansion requires verified FASTQ bytes/MD5, reference admission, adequate
free storage and review of the preceding stage. Failure or inconclusive depth
stops biological claims for that group, not the documentation of the failure.
