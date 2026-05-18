# Stage 8b: Tiny Real Execution — Design Spec

**Date:** 2026-05-18
**Status:** design review
**Scope:** one real ChIP-seq PE no-control execution through preprocessing + signal
**Excluded:** MACS3 peak calling, CUT&Tag execution, CI wiring, Conda env cleanup

---

## 1. Goals

1. Prove the workflow runs end-to-end with real tool invocations on a laptop,
   not just a dry-run DAG check.
2. Keep all generated fixtures and outputs under `/tmp` — zero binary files
   committed to the repository.
3. Intentionally skip MACS3 peak calling. Stage 8a already verifies MACS3,
   IDR, and SEACR DAG connectivity via dry-run. Ultra-tiny synthetic data
   cannot reliably produce peaks, and designing a peak-calling fixture is a
   separate research task deferred to a later slice.
4. Provide a single reproducible command:
   `python3 test/test_stage8b_tiny_execution.py`

---

## 2. Execution scope

Stage 8b runs **one** real execution path:

```
FASTQ → FastQC → Trim Galore → Bowtie2 → sort → index →
MAPQ filter → dedup → final.bam → bamCoverage → CPM.bw
```

Explicit Snakemake targets (no MACS3):

```
results/T1/logs/T1.fastqc.done
results/T1/logs/T1.trim.done
results/T1/02_align/T1.final.bam
results/T1/02_align/T1.final.bam.bai
results/T1/03_bigwig/T1.CPM.bw
```

### Explicit design decision: MACS3 skipped

- Ultra-tiny synthetic data does not produce ChIP-seq enrichment peaks.
  MACS3 `callpeak` would find zero peaks and fail because the rule asserts
  a `.narrowPeak` file must exist.
- Stage 8a dry-run smoke already verifies full MACS3/IDR/SEACR DAG
  connectivity for all seven profiles.
- Real MACS3 execution is deferred to a later slice (small real data or
  a carefully designed peak-calling fixture).
- This decision is documented here so future readers understand the
  intentional gap.

---

## 3. Fixture generation (all under `/tmp`)

Nothing is committed to the repository except the harness script itself.
All fixtures are generated deterministically (seed=42) under
`/tmp/stage8b_XXXX/` at runtime.

### 3.1 Reference FASTA

- 1 contig: `chr1`
- Length: 20 000 bp (20 kb)
- Generated with `random.Random(42)`, uniform ATGC distribution
- No repeated patterns, no low-complexity regions
- 20 kb ensures 1 000 PE reads with ~200 bp fragments map uniquely
  with high MAPQ, surviving `samtools_filter -q 30`

```
>chr1
<20 kb pseudo-random DNA>
```

### 3.2 PE FASTQ reads (compressed `.fq.gz`)

Deterministic, generated from the 20 kb reference with seed=42.

| Parameter | Value |
|-----------|-------|
| Number of pairs | 1 000 |
| Read length | 50 bp |
| Insert size | 200 bp (fixed) |
| Quality string | `I` repeated 50× (Phred+33, constant high quality) |
| Format | **Compressed `.fq.gz`** (written via Python `gzip.open`) |

**Fragment layout:**

```
|-------- insert_size=200 --------|
[start]                   [start+insert_size]
|-- R1=50bp --|                |-- R2=50bp --|
              [start+insert_size-read_len]
```

For fragment start position `s` (0-based):

- R1 forward: `reference[s : s + 50]`
- R2 reverse-complement: `reference[s + 200 - 50 : s + 200]`
  i.e. `reference[s + 150 : s + 200]`

**Start position spacing:**

```
max_start = ref_len - insert_size          # 20000 - 200 = 19800
start_i   = round(i * max_start / (n_pairs - 1))   # for i in 0..999
```

This spreads 1 000 read pairs evenly across the valid range of the
reference, ensuring every pair lands within `chr1`.

**Read name format:**

R1 reads: `@read000001/1` through `@read001000/1`
R2 reads: `@read000001/2` through `@read001000/2`

Zero-padded 6-digit numbering ensures correct lexicographic pairing
by downstream tools. The `/1` and `/2` suffixes are the standard
convention for paired-end read identifiers.

**Why compressed `.fq.gz`:** Trim Galore's output format mirrors the input
format.  Compressed `.fq.gz` input causes Trim Galore to produce
`*_val_1.fq.gz` / `*_val_2.fq.gz` output, which matches the workflow
rule's hardcoded output globs.  Uncompressed `.fq` input would produce
`.fq` output that the rule's `mv "$TMPD"/*_val_1.fq.gz"` glob cannot
find, causing a shell error.  Files are under `/tmp` — the repo's
`*.gz` gitignore rule is irrelevant.

### 3.3 Bowtie2 index

Built at runtime by the harness with `bowtie2-build` from the active
Conda environment:

```bash
bowtie2-build /tmp/stage8b_XXXX/ref.fa /tmp/stage8b_XXXX/bt2_idx
```

The index prefix path is written into the config so the Bowtie2 rule
can find it. Not committed.

### 3.4 chrom.sizes

Not generated. `bamCoverage` in CPM mode does not require
`--chromSizes`. `genome_resources.hs.chrom_sizes` stays `""`.

### 3.5 Profile config

Based on Stage 8a `chipseq_pe_noctrl` baseline with these overrides:

```yaml
samples: "/tmp/stage8b_XXXX/samples.tsv"
outdir: "results"      # resolves under cwd (/tmp/stage8b_XXXX)
threads: 1
trim: true             # Trim Galore runs; .fq.gz outputs under /tmp are safe
remove_dup: "no"       # exercises final.bam contract, skips actual dedup
use_control: false
multiqc: false
stage4b: true
stage5: false
qc:
  blacklist_filter: false
  frip: false
  library_complexity: false
  nrf_pbc: false
  signal_tracks: false
  summary: false
  cuttag_fragment_size: false
```

`trim: true` means Trim Galore actually executes (not just the no-trim
symlink branch). Trim Galore's `.fq.gz` outputs land under `/tmp` and
are cleaned up with the workdir — the repo's `*.gz` gitignore is
irrelevant.

### 3.6 Profile samples.tsv

Single treatment row, PE, ChIP-seq narrow:

```
sample	fastq_1	fastq_2	layout	assay	target	peak_mode	genome	bowtie2_index	experiment	condition	replicate	biological_replicate	technical_replicate	role	control_sample	control_bam
T1	/tmp/stage8b_XXXX/R1.fq.gz	/tmp/stage8b_XXXX/R2.fq.gz	PE	chipseq	H3K27ac	narrow	hs	/tmp/stage8b_XXXX/bt2_idx	exp1	H3K27ac	1	1	1	treatment
```

All paths are absolute temp-directory paths. Nothing relative that
could accidentally resolve into the repo.

---

## 4. Harness: `test/test_stage8b_tiny_execution.py`

### 4.1 Form

Single Python script, stdlib only (`subprocess`, `tempfile`, `shutil`,
`os`, `sys`, `random`). Independent of the Stage 8a harness. Follows
the same patterns: temp workdir, PASS/FAIL reporting, exit code.

### 4.2 Tool availability check (before anything else)

The harness checks for 8 required executables:

```python
REQUIRED_TOOLS = [
    "snakemake",         # via SNAKEMAKE resolution
    "bowtie2-build",
    "bowtie2",
    "samtools",
    "fastqc",
    "trim_galore",
    "cutadapt",
    "bamCoverage",
]
```

Resolution order for snakemake: `SNAKEMAKE` env var → known Conda path
(`/home/irenadler/miniconda3/envs/chipseq/bin/snakemake`) → bare
`"snakemake"` on PATH. For all others: `shutil.which()`.

If `SNAKEMAKE` or the resolved snakemake path is an absolute path to a
known Conda `bin/` directory, the harness prepends that `bin/` directory
to `PATH` for all tool-discovery and subprocess calls. This ensures
`bowtie2-build`, `bowtie2`, `samtools`, `fastqc`, `trim_galore`,
`cutadapt`, and `bamCoverage` are found from the same Conda environment
even when the env is not fully activated (e.g. `conda activate` not run).

If any tool is missing or not executable, print:
`SKIP: <tool> not found` and **exit with code 2**. Exit code 2
distinguishes "environment not ready" from test failure (exit code 1).

### 4.3 Execution flow

```
1. SKIP if any required tool is missing → exit 2

2. WORKDIR = mktemp -d /tmp/stage8b_XXXX/

3. Generate ref.fa (20 kb pseudo-random DNA, seed=42) → WORKDIR/ref.fa

4. Generate R1.fq.gz, R2.fq.gz (1 000 PE pairs, seed=42) → WORKDIR/

5. bowtie2-build WORKDIR/ref.fa WORKDIR/bt2_idx

6. Write WORKDIR/samples.tsv (absolute paths to R1.fq.gz, R2.fq.gz, bt2_idx)

7. Write WORKDIR/config.yaml (trim: true, remove_dup: no, all QC off,
   samples: WORKDIR/samples.tsv)

8. Run Snakemake (cwd=WORKDIR):
   snakemake -s <REPO>/workflow/Snakefile \
     --configfile WORKDIR/config.yaml \
     --cores 1 \
     results/T1/logs/T1.fastqc.done \
     results/T1/logs/T1.trim.done \
     results/T1/02_align/T1.final.bam \
     results/T1/02_align/T1.final.bam.bai \
     results/T1/03_bigwig/T1.CPM.bw
   All output lands under WORKDIR/results/ — repo is untouched.

9. Assert outputs:
   - fastqc.done      exists (sentinel via touch — may be zero bytes)
   - trim.done        exists (sentinel via touch — may be zero bytes)
   - final.bam        exists + nonzero size
   - final.bam.bai    exists + nonzero size
   - CPM.bw           exists + nonzero size
   - samtools view -c final.bam >= 100  (mapped reads)

10. Assert no new generated artifacts leaked into the repo.
    Capture `git status --short --untracked-files=all` before and after
    the Snakemake run.  Compare the two snapshots: any new paths matching
    forbidden suffixes (results/, .snakemake/, *.fq, *.fq.gz, *.bam,
    *.bai, *.bw) under the repo root are FAIL.  The working tree may
    already contain intentional untracked files (specs, new harness code,
    profile fixtures); the check only flags newly appeared artifacts.

11. rmtree WORKDIR
```

### 4.4 Output assertions detail

| Target | Check |
|--------|-------|
| `fastqc.done` | `os.path.isfile` (sentinel, may be 0 bytes) |
| `trim.done` | `os.path.isfile` (sentinel, may be 0 bytes) |
| `final.bam` | `os.path.isfile` + `os.path.getsize > 0` |
| `final.bam.bai` | `os.path.isfile` + `os.path.getsize > 0` |
| `CPM.bw` | `os.path.isfile` + `os.path.getsize > 0` |
| BAM has reads | `samtools view -c <bam>` ≥ 100 |
| Repo clean | before/after `git status` diff, no new generated artifacts |

### 4.5 Error reporting and exit codes

| Exit code | Meaning |
|-----------|---------|
| 0 | PASS — all steps succeeded, all assertions passed |
| 1 | FAIL — a Snakemake rule failed or an output assertion failed |
| 2 | SKIP — required tool not found, environment not ready |

### 4.6 Environment hygiene

- `PYTHONDONTWRITEBYTECODE=1` set for all subprocess calls.
- `scripts/__pycache__` removed if created.
- Snakemake's `.snakemake/` metadata directory is created under
  `WORKDIR` (CWD during execution), not in the repo.

---

## 5. Rule coverage table

| Rule | Real execution | Verified by |
|------|---------------|-------------|
| fastqc | ✅ | sentinel exists |
| trim_galore | ✅ (trim: true) | sentinel exists |
| bowtie2 | ✅ | BAM has mapped reads |
| samtools_sort | ✅ | BAM is coordinate-sorted |
| samtools_index | ✅ | `.bai` exists, nonzero |
| samtools_filter | ✅ (-q 30) | mapped reads survive |
| duplicate_handling | ✅ | final.bam exists |
| final_bam | ✅ | file exists, nonzero |
| bamCoverage | ✅ | CPM.bw exists, nonzero |
| macs3_callpeak | ❌ deferred | Stage 8a dry-run |
| All IDR / SEACR | ❌ deferred | Stage 8a dry-run |

---

## 6. Runtime estimate

| Step | Time |
|------|------|
| Fixture generation | < 1 s |
| bowtie2-build (20 kb) | < 2 s |
| Snakemake (8 rules, 1 core) | 15–40 s |
| Output assertions | < 2 s |
| **Total** | **20–45 s** |

Laptop-appropriate.

---

## 7. Stage boundaries

| Stage | Content | Excluded |
|-------|---------|----------|
| **8a** | 7 dry-run smoke profiles | Real execution |
| **8b** (this spec) | 1 real ChIP-seq PE no-control run, preprocessing + signal only | MACS3, CUT&Tag, CI |
| **8b2** (future) | MACS3 real execution, CUT&Tag real execution, trim variations | CI |
| **8c** | GitHub Actions workflow, CI triggers, matrix | Container images |
| **8d** | Conda env cleanup, env splitting, docs | New features |

---

## 8. Files created / modified

| File | Action |
|------|--------|
| `test/test_stage8b_tiny_execution.py` | **New** — harness script |
| `KNOWN_ISSUES.md` | **Edit** — mark Stage 8b complete |
| `README.md` | **Edit** — add tiny-execution smoke command |

No existing workflow, config, rule, or test files are modified.

---

## 9. Self-review notes

- No TODOs or TBDs.
- All fixtures generated under `/tmp` at runtime — zero binary files
  committed.
- `trim: true` means Trim Galore actually runs. Its `.fq.gz` outputs
  are under `/tmp` and cleaned up — the repo's `*.gz` gitignore is not
  a constraint on `/tmp` files.
- PE FASTQ reads are written as compressed `.fq.gz` (via Python
  `gzip.open`) so Trim Galore produces `.fq.gz` output matching the
  rule's hardcoded `*_val_1.fq.gz` glob. Files are under `/tmp` — the
  repo's `*.gz` gitignore rule is irrelevant for temp-directory files.
- Compressed fixtures are implemented in Stage 8b (not deferred).
  "Compressed fixtures" was removed from the 8b2 deferred list.
- Sentinel files (`.done`) are asserted for existence only — they are
  created via `touch` and may be zero bytes.
- Missing tools → exit code 2, distinct from PASS (0) and FAIL (1).
- Tool PATH robustness: if the resolved snakemake path is absolute to a
  Conda `bin/` directory, that directory is prepended to `PATH` so all 8
  required tools (including `cutadapt`) are discovered from the same
  environment even without `conda activate`.
- Repo-clean check uses before/after `git status` diff, filtering only
  for forbidden generated suffixes (`results/`, `.snakemake/`, `*.fq`,
  `*.fq.gz`, `*.bam`, `*.bai`, `*.bw`). Pre-existing untracked files
  are not flagged.
- PE read names use `@read000001/1` / `@read000001/2` format with
  zero-padded 6-digit numbering, ensuring correct paired-end pairing.
- Rule coverage table uses `duplicate_handling` (the rule name) rather
  than `picard_markduplicates` (the tool). With `remove_dup: "no"` the
  rule exercises the final.bam contract but does not validate real
  duplicate removal.
- No dependency on machine-specific paths except the optional
  `SNAKEMAKE` env var.
- MACS3 is explicitly documented as intentionally skipped with rationale.
- No CI wiring — that is Stage 8c.
- Runtime under 60 seconds on a laptop.
