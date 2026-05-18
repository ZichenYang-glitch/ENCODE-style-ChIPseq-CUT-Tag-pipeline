# Stage 8a: Test Profiles and Smoke-Test Harness — Design Spec

**Date:** 2026-05-18
**Status:** design review
**Scope:** test profiles, smoke harness, .gitignore refinement
**Excluded:** tiny FASTQ/index, real execution, CI wiring, Conda env cleanup

## 1. Goals

1. Define a repo-local directory structure for smoke-test profiles.
2. Cover the most important ChIP-seq and CUT&Tag dispatch paths with
   minimal, text-only profile fixtures.
3. Provide a single reproducible smoke command that validates every
   profile and dry-runs the Snakemake DAG.
4. Never write generated outputs into the repository. All temporary
   files live under `/tmp`.
5. Keep full backward compatibility: no DAG, config, sample-sheet, or
   rule-file changes for ordinary users.

---

## 2. Profile directory structure

```
test/profiles/
├── chipseq_se_noctrl/
│   ├── config.yaml
│   └── samples.tsv
├── chipseq_pe_noctrl/
│   ├── config.yaml
│   └── samples.tsv
├── chipseq_pe_ctrlsample/
│   ├── config.yaml
│   └── samples.tsv
├── cuttag_pe_noctrl/
│   ├── config.yaml
│   └── samples.tsv
├── cuttag_pe_seacr/
│   ├── config.yaml
│   └── samples.tsv
├── chipseq_idr_dryrun/
│   ├── config.yaml
│   └── samples.tsv
└── chipseq_pe_external_ctrlbam/
    ├── config.yaml
    └── samples.tsv
```

Every profile directory contains exactly two tracked text files:
`config.yaml` and `samples.tsv`. No placeholder FASTQs, BAMs, or
index files are committed. The smoke harness creates any required
empty placeholder files under `/tmp` at runtime.

---

## 3. Profile definitions

All profiles share a common configuration baseline:

```yaml
outdir: "results"
threads: 1
mapq: 30
binsize: 10
remove_dup: "auto"
trim: true
extend_reads: "auto"
use_control: false
multiqc: false
stage4b: true
stage5: false
genome_resources:
  hs:
    effective_genome_size: "hs"
    chrom_sizes: ""
    blacklist: ""
    gtf: ""
    reference_fasta: ""
qc:
  blacklist_filter: false
  frip: false
  library_complexity: false
  nrf_pbc: false
  signal_tracks: false
  summary: false
  cuttag_fragment_size: false
```

Per-profile overrides are listed below. The `samples` key in each
profile's `config.yaml` points to `samples.tsv` (relative to the
profile directory). The harness rewrites this to an absolute `/tmp`
path at runtime.

Every profile's `samples.tsv` MUST use the full canonical header
with all 17 columns, even when some are left empty for a given
profile. This avoids relying on column-defaulting behavior and
makes control_sample / control_bam profiles explicit:

```
sample	fastq_1	fastq_2	layout	assay	target	peak_mode	genome	bowtie2_index	experiment	condition	replicate	biological_replicate	technical_replicate	role	control_sample	control_bam
```

### 3.1 `chipseq_se_noctrl` — ChIP-seq SE, no control

- Config: baseline (no overrides).
- Samples: 1 treatment row, SE layout, assay=chipseq, peak_mode=narrow,
  target=H3K27ac, genome=hs, bowtie2_index=bt2_idx.
- FASTQ paths: `R1.fq` (single-end).

### 3.2 `chipseq_pe_noctrl` — ChIP-seq PE, no control

- Config: baseline.
- Samples: 1 treatment row, PE layout, assay=chipseq, peak_mode=narrow,
  target=H3K27ac, genome=hs, bowtie2_index=bt2_idx.
- FASTQ paths: `R1.fq`, `R2.fq` (paired-end).

### 3.3 `chipseq_pe_ctrlsample` — ChIP-seq PE, FASTQ-based control (primary control profile)

- Config: baseline + `use_control: true`.
- Samples: 2 rows.
  - Row 1: treatment, PE, chipseq, narrow, H3K27ac, hs, bt2_idx,
    experiment=exp1, biological_replicate=1,
    control_sample=Ctrl1. FASTQs: `R1.fq`, `R2.fq`.
  - Row 2: control, PE, chipseq, narrow, Input, hs, bt2_idx,
    experiment=exp1, biological_replicate=1, role=control.
    FASTQs: `R3.fq`, `R4.fq`.
- This is the **preferred** smoke-test coverage for control handling.
  Both treatment and control enter through the full FASTQ preprocessing
  pipeline.

### 3.4 `cuttag_pe_noctrl` — CUT&Tag PE, no control

- Config: baseline + `qc.cuttag_fragment_size: true`.
- Samples: 1 treatment row, PE, assay=cuttag, peak_mode=narrow,
  target=H3K4me3, genome=hs, bowtie2_index=bt2_idx.
- FASTQ paths: `R1.fq`, `R2.fq`.
- Verifies the CUT&Tag fragment-size QC rule appears in the DAG.

### 3.5 `cuttag_pe_seacr` — CUT&Tag PE, SEACR enabled

- Config: same as `cuttag_pe_noctrl` plus:
  ```yaml
  cuttag:
    seacr:
      enabled: true
      mode: stringent
      normalization: non
      threshold: 0.01
  ```
- Samples: identical to `cuttag_pe_noctrl`.
- Verifies `seacr_bedgraph` and `seacr_call` rules appear in the DAG.

### 3.6 `chipseq_idr_dryrun` — Stage 5 IDR, dry-run only

- Config: baseline + `stage5: true` + `idr` block:
  ```yaml
  idr:
    seed: 42
    threshold: 0.05
    rank: "p.value"
  ```
- Samples: 2 treatment rows, PE, assay=chipseq, peak_mode=narrow,
  target=CTCF, genome=hs, bowtie2_index=bt2_idx,
  experiment=exp1, biological_replicate=1 and 2.
  - Row 1: `R1.fq`, `R2.fq`, biological_replicate=1.
  - Row 2: `R3.fq`, `R4.fq`, biological_replicate=2.
- Verifies the full Stage 5 IDR DAG resolves: IDR-ready MACS3 peaks,
  true-replicate IDR, pseudoreplicate BAM splitting, self-IDR,
  pooled-IDR, and final conservative/optimal peak sets.

### 3.7 `chipseq_pe_external_ctrlbam` — External control_bam (backward compatibility)

- Config: baseline + `use_control: true`.
- Samples: 1 treatment row, PE, assay=chipseq, peak_mode=narrow,
  target=H3K27ac, genome=hs, bowtie2_index=bt2_idx,
  control_bam=ctrl.bam.
- FASTQ paths: `R1.fq`, `R2.fq`.
- This profile covers the external `control_bam` path for backward
  compatibility. FASTQ-based `control_sample` (profile 3.3) is the
  preferred workflow for normal use. The profile is committed and
  run by the harness alongside the other six profiles.

---

## 4. Smoke harness: `test/test_stage8_smoke_profiles.py`

### 4.1 Form

A single Python script using stdlib only (`subprocess`, `tempfile`,
`shutil`, `os`, `sys`, `csv`). No package dependencies. Follows the
existing stress-test pattern established by `test_stage4b_stress.py`,
`test_stage7a_stress.py`, and `test_stage7b_stress.py`.

### 4.2 Execution flow (per profile)

```
Given profile name (e.g. "chipseq_se_noctrl"):

1. WORKDIR = mktemp -d /tmp/smoke_profile_XXXX/

2. Read test/profiles/<name>/samples.tsv to discover paths that
   need placeholder files:
   - fastq_1 value
   - fastq_2 value (only when non-empty)
   - control_bam value (only for external_ctrlbam profile)
   bowtie2_index is NOT used to create a placeholder — Bowtie2
   indexes are multi-file prefixes; a fake single file is misleading.
   Real tiny indexes belong in Stage 8b.

3. Create each placeholder as an empty file under WORKDIR.

4. Copy samples.tsv to WORKDIR/samples.tsv.
   Paths inside the TSV stay unchanged; they resolve correctly when
   CWD is WORKDIR.

5. Copy config.yaml to WORKDIR/config.yaml, rewriting the samples
   key to the absolute path WORKDIR/samples.tsv.

6. Validate (cwd=WORKDIR, absolute script/config paths):
   python3 <REPO_ROOT>/scripts/validate_samples.py \
     --config <WORKDIR>/config.yaml
   Using cwd=WORKDIR ensures relative paths such as control_bam
   (checked via os.path.isfile) resolve against the temp directory
   where placeholders were created.

7. Dry-run (cwd=WORKDIR, absolute paths):
   snakemake -s <REPO_ROOT>/workflow/Snakefile \
     --configfile <WORKDIR>/config.yaml \
     -n --quiet

8. Profile-specific rule scheduling checks (cwd=WORKDIR, absolute paths).
   --list-rules only proves a rule is defined, not that it is scheduled
   for a profile. To verify scheduling, re-run the dry-run without
   --quiet and grep the job listing for expected rule names:
   snakemake -s <REPO_ROOT>/workflow/Snakefile \
     --configfile <WORKDIR>/config.yaml -n
   The output (stdout) lists each scheduled job as "rule <name>:" lines.

   Profile                  | Expected rules in dry-run job listing
   -------------------------|----------------------------------------
   cuttag_pe_noctrl         | cuttag_fragment_size
   cuttag_pe_seacr          | seacr_bedgraph, seacr_call
   chipseq_idr_dryrun       | macs3_idr_biorep, idr_true_replicates,
                            | split_pseudoreps, idr_self_pseudoreps,
                            | idr_pooled_pseudoreps, stage5b_summary

   --list-rules is still useful as a debug aid (verify a rule exists
   in the Snakefile at all), but not as a scheduling gate.

9. rmtree WORKDIR
```

### 4.3 Snakemake executable

Discovered via the `SNAKEMAKE` environment variable, falling back to
`"snakemake"`. The Conda environment `chipseq` provides the binary at
a known path; the harness documents how to set the variable:

```bash
SNAKEMAKE=/home/irenadler/miniconda3/envs/chipseq/bin/snakemake \
  python3 test/test_stage8_smoke_profiles.py
```

### 4.4 Environment hygiene

- `PYTHONDONTWRITEBYTECODE=1` is set in the subprocess environment
  for every `python3` and `snakemake` invocation.
- After all profiles complete, the harness removes `scripts/__pycache__`
  if it exists.
- Subprocess `cwd` is set to `WORKDIR` for **all three steps**
  (validate, dry-run, --list-rules) so that relative paths in
  `samples.tsv` (FASTQs, bowtie2_index prefix, control_bam) resolve
  against the temp directory where placeholders were created.

### 4.5 Error reporting and exit code

- Each profile prints a single line: `PASS: <name>` or
  `FAIL: <name> — <reason>`.
- After all profiles: `X/Y profiles passed`.
- Exit code 0 when all profiles pass; exit code 1 when any profile
  fails.
- On failure, the harness prints the captured stderr (last 500 chars)
  so the user can diagnose without re-running manually.

### 4.6 Profile discovery

The harness defines an explicit list of 7 profiles to run. All seven
are committed and executed by default:

```python
PROFILES = [
    "chipseq_se_noctrl",
    "chipseq_pe_noctrl",
    "chipseq_pe_ctrlsample",
    "cuttag_pe_noctrl",
    "cuttag_pe_seacr",
    "chipseq_idr_dryrun",
    "chipseq_pe_external_ctrlbam",
]
```

`chipseq_pe_external_ctrlbam` is marked as backward-compatibility
coverage — it is not the preferred workflow — but it is a normal
committed profile that runs every time. If a profile directory is
missing, the harness fails with a clear error rather than silently
skipping.

---

## 5. `.gitignore` changes

### 5.1 Current rules

```
test/
reference/
miniconda3/
.snakemake/
results/
*.sh
*.gz
*.bam
```

### 5.2 Proposed rules

```
# --- Test data files (generated at runtime, never committed) ---
test/**/*.fq
test/**/*.fastq
test/**/*.fastq.gz
test/**/*.fq.gz
test/**/*.bam
test/**/*.bai
test/**/results/
test/**/.snakemake/
test/**/__pycache__/
test/**/.pytest_cache/

# --- Repo-level generated directories ---
reference/
miniconda3/
.snakemake/
results/

# --- Legacy ignores ---
*.sh
*.gz
```

### 5.3 Rationale

- The broad `test/` directory ignore is removed. Without a parent
  ignore, everything under `test/` is tracked by default — no `!`
  negation patterns are needed. Negation patterns are only effective
  when a parent rule ignores the directory first; removing the parent
  ignore makes them unnecessary.
- Data-file patterns (`*.fq`, `*.fastq`, `*.bam`, etc.) are scoped
  to `test/**/` so they do not affect real workflow outputs outside
  `test/`. These prevent accidental commits of generated test
  artifacts while allowing `test/*.py` and `test/profiles/**` to be
  tracked normally without `git add -f`.
- The global `*.bam` ignore is removed. Real BAM outputs under
  `results/` are already covered by the `results/` directory ignore.
  Test BAM placeholders are covered by `test/**/*.bam`.
- The legacy `*.sh` and `*.gz` global ignores are preserved for
  backward compatibility.
- **Stage 8b note:** The global `*.gz` ignore will conflict with
  committing tiny compressed FASTQ fixtures (e.g. `*.fastq.gz`).
  Stage 8b must either scope the `*.gz` ignore to specific
  directories or replace it with more precise patterns before
  committing any `.gz` test fixtures.

---

## 6. Stage boundaries

| Stage | Content | Explicitly excluded |
|-------|---------|---------------------|
| **8a** (this spec) | 7 profiles, smoke harness script, `.gitignore` refinement, dry-run + validation only | Real FASTQ data, Bowtie2 index construction, real execution |
| **8b** | Tiny FASTQ subsets, tiny Bowtie2 index, at least one profile run end-to-end | CI integration |
| **8c** | GitHub Actions workflow, CI triggers, matrix strategy | Container images |
| **8d** | Conda env cleanup (remove `prefix`), possible env splitting, workstation/server documentation | New pipeline features |

---

## 7. Backward compatibility

- `config/config.yaml` — unchanged.
- `config/samples.tsv` — unchanged.
- `workflow/Snakefile` — unchanged.
- All existing rule files under `workflow/rules/` — unchanged.
- `scripts/validate_samples.py` — unchanged.
- Existing stress tests (`test/test_validation_stress.py`,
  `test/test_stage4b_stress.py`, `test/test_stage7a_stress.py`,
  `test/test_stage7b_stress.py`) — unchanged.
- Ordinary users running `snakemake -s workflow/Snakefile
  --configfile config/config.yaml` see zero behavioral difference.

---

## 8. Implementation notes

- New profile files committed under `test/profiles/` must not require
  `git add -f` after the `.gitignore` update.
- No profile file contains absolute paths; all paths are relative
  simple names (`R1.fq`, `bt2_idx`, `ctrl.bam`).
- The harness is the only component that creates files outside the
  repository tree (under `/tmp`), and it cleans them up.
- Snakemake's `--quiet` flag suppresses per-job logging during dry-run
  while still reporting DAG errors on stderr.

---

## 9. Self-review notes

- No TODOs or TBDs remain in this spec.
- All seven profile names are consistent between the directory
  structure (Section 2), profile definitions (Section 3), and harness
  profile list (Section 4.6).
- Profile 3.7 (`chipseq_pe_external_ctrlbam`) is committed and runs
  by default alongside the other six. It is "legacy" in workflow
  preference, not in file presence or harness execution.
- Placeholder strategy matches validation behavior:
  `validate_samples.py` does **not** check FASTQ or bowtie2_index
  file existence; it **does** check `control_bam` existence when
  `use_control: true`. The harness creates an empty `ctrl.bam` only
  for the external_ctrlbam profile. All three subprocess steps
  (validate, dry-run, scheduling check) use `cwd=WORKDIR` so
  relative paths resolve correctly.
- Rule scheduling is verified via non-quiet dry-run (`-n`) job
  listings, not `--list-rules`. The latter only proves a rule is
  defined, not that it is scheduled for a given profile.
- Stage 8a produces zero generated files in the repository. The
  harness writes exclusively to `/tmp` and cleans up.
- Runtime estimate: ~2-5 seconds per profile on a laptop (fork +
  Python import + Snakemake DAG resolution). Total smoke suite:
  under 35 seconds for 7 profiles.
- `.gitignore` changes do not affect any existing tracked files.
- Sample-sheet header has 17 columns. The harness must tolerate
  empty trailing fields (tab-separated, no content) for rows that
  don't use all optional columns.
