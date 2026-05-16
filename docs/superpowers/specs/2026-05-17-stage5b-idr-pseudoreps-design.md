# Stage 5b: Pseudoreplicate IDR and Final Reproducibility Outputs

## Scope

Stage 5b adds pseudoreplicate-based IDR and final reproducibility outputs on top
of the Stage 5a true-replicate IDR foundation.

**Stage 5b includes:**
- Deterministic pseudoreplicate BAM splitting via `scripts/split_pseudoreps.py`
- IDR-ready MACS3 calls on all pseudorep BAMs (reusing `_idr_macs3_args`)
- Self-pseudoreplicate IDR (per biological replicate)
- Pooled-pseudoreplicate IDR
- Final conservative/optimal peak set assembly
- Reproducibility QC summary (`reproducibility_summary.tsv`)

**Deferred to Stage 6/7:**
- Histone broad-peak IDR
- CUT&Tag IDR
- 3+ biological replicate support
- MultiQC integration of IDR results

## Configuration

### config.yaml

```yaml
stage5: false

idr:
  seed: 42               # deterministic pseudorep split seed (positive int)
  threshold: 0.05        # IDR threshold (Stage 5a)
  rank: "p.value"        # IDR ranking measure (Stage 5a)

tool_parameters:
  idr_macs3:
    pvalue: 0.1
    extra_args: ""
```

`idr.seed` feeds the deterministic pseudoreplicate split script. No separate
`pseudorep.seed` key.

### Validated Config Keys

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `stage5` | boolean/string boolean | `false` | Master IDR gate |
| `idr.seed` | positive int | 42 | Seed for deterministic pseudoreplicate split |
| `idr.threshold` | float in (0,1) | 0.05 | IDR threshold |
| `idr.rank` | `"p.value"` or `"signal.value"` | `"p.value"` | IDR ranking column |

### Validation Additions

**`idr.seed` in `_validate_idr_settings()`:**

```python
seed = idr.get("seed", 42)
if isinstance(seed, bool):
    raise ValidationError("idr.seed must be a positive integer")
if isinstance(seed, int):
    if seed <= 0:
        raise ValidationError("idr.seed must be positive, got {}".format(seed))
elif isinstance(seed, str):
    text = seed.strip()
    if not text.isdigit() or int(text) <= 0:
        raise ValidationError(
            "idr.seed must be a positive integer, got {!r}".format(seed)
        )
    seed = int(text)
else:
    raise ValidationError(
        "idr.seed must be a positive integer, got {!r}".format(seed)
    )
```

**Unknown-key validation in `idr` block:**

`_validate_idr_settings()` must reject unknown keys in the `idr` mapping.
Known keys: `seed`, `threshold`, `rank`. Any other key raises:

```
"idr: unknown key {key!r}. Known: ['rank', 'seed', 'threshold']"
```

No new eligibility checks — Stage 5b reuses Stage 5a constraints:
chipseq + narrow + exactly 2 bio-reps.

## Output Namespace

```
results/experiments/<experiment>/
├── 05_pseudorep/                              # Stage 5b NEW
│   ├── <exp>_biorep<bio_rep1>.pr1.bam         # (and .bai for each)
│   ├── <exp>_biorep<bio_rep1>.pr2.bam
│   ├── <exp>_biorep<bio_rep2>.pr1.bam
│   ├── <exp>_biorep<bio_rep2>.pr2.bam
│   ├── <exp>_pooled.pr1.bam
│   └── <exp>_pooled.pr2.bam
├── 04_peaks/idr/                              # Stage 5a + 5b
│   ├── <exp>_biorep<bio_rep1>_idr_peaks.narrowPeak      # 5a
│   ├── <exp>_biorep<bio_rep2>_idr_peaks.narrowPeak      # 5a
│   ├── <exp>_biorep<bio_rep1>_pr1_idr_peaks.narrowPeak  # 5b
│   ├── <exp>_biorep<bio_rep1>_pr2_idr_peaks.narrowPeak  # 5b
│   ├── <exp>_biorep<bio_rep2>_pr1_idr_peaks.narrowPeak  # 5b
│   ├── <exp>_biorep<bio_rep2>_pr2_idr_peaks.narrowPeak  # 5b
│   ├── <exp>_pooled_pr1_idr_peaks.narrowPeak            # 5b
│   └── <exp>_pooled_pr2_idr_peaks.narrowPeak            # 5b
├── 06_idr/
│   ├── true_replicates/                       # Stage 5a
│   │   ├── idr.txt
│   │   └── idr.thresholded.narrowPeak
│   ├── self_pseudoreplicates/                 # Stage 5b NEW
│   │   ├── biorep<bio_rep1>.idr.txt
│   │   ├── biorep<bio_rep1>.idr.thresholded.narrowPeak
│   │   ├── biorep<bio_rep2>.idr.txt
│   │   └── biorep<bio_rep2>.idr.thresholded.narrowPeak
│   ├── pooled_pseudoreplicates/               # Stage 5b NEW
│   │   ├── idr.txt
│   │   └── idr.thresholded.narrowPeak
│   └── final/                                 # Stage 5b NEW
│       ├── conservative.narrowPeak
│       ├── optimal.narrowPeak
│       └── reproducibility_summary.tsv
└── logs/
    ├── ... (separate raw/thresholded IDR logs)
    ├── <experiment>_<source>.pr{1,2}.idr.macs3.log
    ├── <exp>_<source>.split.log
    └── <experiment>.stage5b.summary.log
```

## New Script: `scripts/split_pseudoreps.py`

### Interface

```
python3 scripts/split_pseudoreps.py \
    --input in.bam \
    --out1 <exp>_biorepN.pr1.bam \
    --out2 <exp>_biorepN.pr2.bam \
    --seed 42 \
    --threads 4
```

### Behavior

1. Stream reads via `samtools view -h <input>`.
2. For each read, extract the query name. Compute the canonical template name
   by stripping trailing `/1` or `/2` suffixes from the query name. Do NOT
   modify the read record itself — only use the canonical name for hashing.
3. Compute a deterministic hash of `seed + canonical_name` using **hashlib.sha256**.
   Do NOT use Python's built-in `hash()` (its seed changes across interpreter
   restarts and is not reproducible). Use `hashlib.sha256(f"{seed}:{canonical_name}".encode()).hexdigest()`
   or equivalent.
4. Convert the hex digest to an integer. Assign to pr1 if the integer is even,
   pr2 if odd. Both mates of a PE pair always hash to the same canonical name,
   so they always go to the same pseudoreplicate.
5. Write header lines to both output streams. Pipe each stream through
   `samtools sort -@ threads -o <output>`.
6. Run `samtools index` on both output BAMs.
7. If any subprocess fails, exit non-zero with a clear error message.
8. Ensure complementarity by counting with `samtools view -c`:
   `samtools view -c pr1.bam + samtools view -c pr2.bam == samtools view -c input.bam`.
   Do NOT use `wc -l` on SAM output — BAM is compressed, and SAM line count
   includes headers.
9. Do not drop reads — no intentional read loss.
10. Check that `samtools` is callable before starting. Fail clearly if not found.

### Dependencies

- Python 3 stdlib only: `subprocess`, `hashlib` (sha256), `argparse`, `os`, `sys`
- `samtools` on PATH (view, sort, index)

## Derived Structures (Snakefile)

### Precomputed Expansion Lists

For each experiment in `IDR_EXPERIMENTS`, derive:

```python
# Source labels for pseudorep splitting: biorep<label_a>, biorep<label_b>, pooled
IDR_SPLIT_SOURCE_EXP = []   # experiment column
IDR_SPLIT_SOURCE_NAME = []  # source name column (zip-compatible)

# Pseudorep MACS3 targets: each source × pr{1,2}
IDR_PR_PEAK_EXP = []
IDR_PR_PEAK_SRC = []   # source name
IDR_PR_PEAK_PR = []    # "1" or "2"

# Self-IDR targets: one per bio_rep
IDR_SELF_EXP = []
IDR_SELF_BR = []       # actual bio_rep label

for exp in IDR_EXPERIMENTS:
    bioreps = sorted(_bioreps_for(exp, "treatment"))
    br_a = bioreps[0]
    br_b = bioreps[1]

    # Split sources: 2 bioreps + pooled
    for br in (br_a, br_b):
        src = f"biorep{br}"
        IDR_SPLIT_SOURCE_EXP.append(exp)
        IDR_SPLIT_SOURCE_NAME.append(src)
        # Pseudorep MACS3: pr1, pr2 for each source
        for pr in ("1", "2"):
            IDR_PR_PEAK_EXP.append(exp)
            IDR_PR_PEAK_SRC.append(src)
            IDR_PR_PEAK_PR.append(pr)
        # Self-IDR: one per bio_rep
        IDR_SELF_EXP.append(exp)
        IDR_SELF_BR.append(br)

    # Pooled split + pooled pseudorep MACS3
    IDR_SPLIT_SOURCE_EXP.append(exp)
    IDR_SPLIT_SOURCE_NAME.append("pooled")
    for pr in ("1", "2"):
        IDR_PR_PEAK_EXP.append(exp)
        IDR_PR_PEAK_SRC.append("pooled")
        IDR_PR_PEAK_PR.append(pr)
```

## New Rules (in `workflow/rules/idr.smk`)

### `split_pseudoreps` — deterministic BAM splitting

```python
rule split_pseudoreps:
    output:
        pr1 = f"{OUTDIR}/experiments/{{experiment}}/05_pseudorep/{{experiment}}_{{source}}.pr1.bam",
        p1i = f"{OUTDIR}/experiments/{{experiment}}/05_pseudorep/{{experiment}}_{{source}}.pr1.bam.bai",
        pr2 = f"{OUTDIR}/experiments/{{experiment}}/05_pseudorep/{{experiment}}_{{source}}.pr2.bam",
        p2i = f"{OUTDIR}/experiments/{{experiment}}/05_pseudorep/{{experiment}}_{{source}}.pr2.bam.bai",
    input:
        lambda wc: _split_input(wc),
    params:
        seed    = IDR_SEED,
    wildcard_constraints:
        source  = r"pooled|biorep\d+",
    log:
        f"{OUTDIR}/experiments/{{experiment}}/logs/{{experiment}}_{{source}}.split.log",
    threads: THREADS,
    conda:
        "../envs/chipseq.yml",
    shell:
        """
        mkdir -p "$(dirname {output.pr1})" "$(dirname {log})"
        python3 scripts/split_pseudoreps.py \
            --input {input:q} \
            --out1 {output.pr1:q} \
            --out2 {output.pr2:q} \
            --seed {params.seed} \
            --threads {threads} \
            2>&1 | tee {log:q}
        """
```

**`_split_input()`:** if `source == "pooled"`, return `{experiment}.pooled.final.bam`.
If `source` starts with `"biorep"`, parse the bio_rep label and return that
biorep BAM. Never hardcode labels 1/2.

### `macs3_idr_pseudorep` — IDR-ready MACS3 on pseudorep BAM

```python
rule macs3_idr_pseudorep:
    output:
        f"{OUTDIR}/experiments/{{experiment}}/04_peaks/idr/{{experiment}}_{{source}}_pr{{pr}}_idr_peaks.narrowPeak",
    input:
        lambda wc: _idr_pseudorep_inputs(wc),
    params:
        macs3_args = lambda wc: _idr_macs3_args(wc),
        sample     = lambda wc: f"{wc.experiment}_{wc.source}_pr{wc.pr}_idr",
    wildcard_constraints:
        source = r"pooled|biorep\d+",
        pr     = r"[12]",
    log:
        f"{OUTDIR}/experiments/{{experiment}}/logs/{{experiment}}_{{source}}_pr{{pr}}.idr.macs3.log",
    threads: THREADS,
    conda:
        "../envs/chipseq.yml",
    shell:
        """
        set -e -o pipefail
        mkdir -p "$(dirname {output})" "$(dirname {log})"

        set -- {input:q}
        TREATMENT="$1"
        if [[ $# -ge 3 ]]; then
            macs3 callpeak -t "$TREATMENT" -c "$3" \
                -n {params.sample:q} \
                --outdir "$(dirname {output:q})" \
                {params.macs3_args} 2>&1 | tee {log:q}
        else
            macs3 callpeak -t "$TREATMENT" \
                -n {params.sample:q} \
                --outdir "$(dirname {output:q})" \
                {params.macs3_args} 2>&1 | tee {log:q}
        fi

        [[ -f "{output}" ]] || {{ echo "ERROR: Expected {output}" >&2; exit 1; }}
        """
```

**`_idr_pseudorep_inputs()`**: returns `[pseudorep_bam, pseudorep_bam.bai,
...optional_pooled_control_bam]`. Pooled control BAM dependency mirrors
`_idr_biorep_peaks_inputs`. Reuses `_idr_macs3_args()` from Stage 5a —
same `-p` relaxed threshold, same `-f`/`-g` format/genome resolution.

### `idr_self_pseudoreps` — self-IDR per biorep

Two-invocation pattern matching Stage 5a `idr_true_replicates`:

```python
rule idr_self_pseudoreps:
    output:
        idr_out = f"{OUTDIR}/experiments/{{experiment}}/06_idr/self_pseudoreplicates/biorep{{bio_rep}}.idr.txt",
        thr_out = f"{OUTDIR}/experiments/{{experiment}}/06_idr/self_pseudoreplicates/biorep{{bio_rep}}.idr.thresholded.narrowPeak",
    input:
        peaks1 = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/04_peaks/idr/"
            f"{wc.experiment}_biorep{wc.bio_rep}_pr1_idr_peaks.narrowPeak"
        ),
        peaks2 = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/04_peaks/idr/"
            f"{wc.experiment}_biorep{wc.bio_rep}_pr2_idr_peaks.narrowPeak"
        ),
    params:
        threshold     = IDR_THRESHOLD,
        rank          = IDR_RANK,
        raw_log_file  = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/logs/"
            f"{wc.experiment}.idr.self.biorep{wc.bio_rep}.raw.log"
        ),
        thr_log_file  = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/logs/"
            f"{wc.experiment}.idr.self.biorep{wc.bio_rep}.thr.log"
        ),
    ...
```

Two invocations (same pattern as Stage 5a):
1. Raw: `idr --samples ... --output-file {output.idr_out} --log-output-file {params.raw_log_file}`
2. Thresholded: `idr --samples ... --idr-threshold {params.threshold} --output-file {output.thr_out} --log-output-file {params.thr_log_file}`

`{bio_rep}` wildcard constraint: `r"\d+"`. Labels come from precomputed
`IDR_SELF_EXP`/`IDR_SELF_BR` expansion lists.

### `idr_pooled_pseudoreps` — pooled-IDR

Same two-invocation pattern. Inputs are `{exp}_pooled_pr1_idr_peaks` and
`{exp}_pooled_pr2_idr_peaks`. Outputs under `06_idr/pooled_pseudoreplicates/`.

### `stage5b_summary` — reproducibility QC and final assembly

Prefer a small Python helper script (`scripts/stage5b_summary.py`) to compute
ratios, handle zero denominators, and emit the reproducibility summary TSV.
No fragile shell arithmetic or heredocs.

```python
rule stage5b_summary:
    output:
        summary      = f"{OUTDIR}/experiments/{{experiment}}/06_idr/final/reproducibility_summary.tsv",
        conservative = f"{OUTDIR}/experiments/{{experiment}}/06_idr/final/conservative.narrowPeak",
        optimal      = f"{OUTDIR}/experiments/{{experiment}}/06_idr/final/optimal.narrowPeak",
    input:
        true_thresh   = f"{OUTDIR}/experiments/{{experiment}}/06_idr/true_replicates/idr.thresholded.narrowPeak",
        pool_thresh   = f"{OUTDIR}/experiments/{{experiment}}/06_idr/pooled_pseudoreplicates/idr.thresholded.narrowPeak",
        self1_thresh  = lambda wc: _self_thresh_path(wc.experiment, 0),
        self2_thresh  = lambda wc: _self_thresh_path(wc.experiment, 1),
    params:
        experiment    = lambda wc: wc.experiment,
        bio_rep_a     = lambda wc: str(sorted(_bioreps_for(wc.experiment, "treatment"))[0]),
        bio_rep_b     = lambda wc: str(sorted(_bioreps_for(wc.experiment, "treatment"))[1]),
    log:
        f"{OUTDIR}/experiments/{{experiment}}/logs/{{experiment}}.stage5b.summary.log",
    conda:
        "../envs/chipseq.yml",
    shell:
        """
        mkdir -p "$(dirname {output.summary})" "$(dirname {log})"
        python3 scripts/stage5b_summary.py \
            --true-peaks {input.true_thresh:q} \
            --pooled-peaks {input.pool_thresh:q} \
            --self1-peaks {input.self1_thresh:q} \
            --self2-peaks {input.self2_thresh:q} \
            --experiment {params.experiment:q} \
            --bio-rep-a {params.bio_rep_a} \
            --bio-rep-b {params.bio_rep_b} \
            --output-tsv {output.summary:q} \
            --output-cons {output.conservative:q} \
            --output-opt {output.optimal:q} \
            2>&1 | tee {log:q}
        """
```

**`scripts/stage5b_summary.py`**: stdlib Python script. Reads peak files with
`open(f).readlines()`, counts peaks (minus header), computes:

```
Nt = line count of true_replicates/thresholded
Np = line count of pooled_pseudoreps/thresholded
N1 = line count of self_pseudoreps/biorep_a/thresholded
N2 = line count of self_pseudoreps/biorep_b/thresholded

rescue_ratio = max(Np, Nt) / min(Np, Nt)   (NA if both 0, inf if one 0)
self_consistency_ratio = max(N1, N2) / min(N1, N2)

reproducibility_status = "pass" if both ratios < 2 and both finite, else "fail"
```

Copies conservative (true_replicates thresholded) and optimal (pooled
pseudoreps thresholded) to `06_idr/final/`. Writes `reproducibility_summary.tsv`
with columns: `experiment`, `true_peaks(Nt)`, `pooled_peaks(Np)`,
`self1_peaks(N1)`, `self2_peaks(N2)`, `rescue_ratio`, `self_consistency_ratio`,
`reproducibility_status`.

**`_self_thresh_path(experiment, index)`:** returns the thresholded narrowPeak
path for self-IDR of the bio_rep at 0-based index.

## Rule All Targets (Snakefile additions)

```python
# Stage 5b: pseudorep BAMs (pr1 + pr2 for each source)
stage5b_split_pr1 = (
    expand(
        "{outdir}/experiments/{experiment}/05_pseudorep/{experiment}_{source}.pr1.bam",
        zip, outdir=OUTDIR,
        experiment=IDR_SPLIT_SOURCE_EXP, source=IDR_SPLIT_SOURCE_NAME,
    )
    if STAGE5 and IDR_EXPERIMENTS else []
),
stage5b_split_pr2 = (
    expand(
        "{outdir}/experiments/{experiment}/05_pseudorep/{experiment}_{source}.pr2.bam",
        zip, outdir=OUTDIR,
        experiment=IDR_SPLIT_SOURCE_EXP, source=IDR_SPLIT_SOURCE_NAME,
    )
    if STAGE5 and IDR_EXPERIMENTS else []
),
# BAIs are produced by the same rule; Snakemake tracks them automatically.

# Stage 5b: pseudorep MACS3 peaks (6 targets per experiment)
...

# Stage 5b: self-IDR + pooled-IDR + final summary
...
```

Including both pr1.bam and pr2.bam in `rule all` ensures Snakemake tracks
both outputs of the split rule. If downstream targets already require these
BAMs (e.g., pseudorep MACS3 input), the pr1/pr2 split targets can be omitted
from `rule all` since they'll be scheduled via DAG resolution. Choose one
approach and be consistent.

## Test Plan

### `test/test_stage5b_stress.py` — DAG and validation tests (13 tests)

Reuses Stage 4b/4c/5a harness: `SNAKEMAKE` env var with conda fallback,
temp file cleanup, `scripts/__pycache__` removal.

1. `idr.seed` positive int → validates
2. `idr.seed` negative → rejected
3. `idr.seed` zero → rejected
4. `idr.seed` non-integer string → rejected
5. `idr.seed = "42"` (string) → accepted and normalized
6. Stage 5 config + 2 bio-reps → `split_pseudoreps` in DAG
7. Pseudorep paths use actual bio_rep labels (fields 2 and 4 in sample sheet)
8. `macs3_idr_pseudorep`, `idr_self_pseudoreps`, `idr_pooled_pseudoreps` in DAG
9. `stage5b_summary` in DAG
10. `split_pseudoreps.py --seed 42` visible in dry-run
11. Self-IDR `--idr-threshold` visible
12. Pooled-IDR `--idr-threshold` visible
13. Unknown key in `idr` block → rejected with clear message

### `test/test_split_pseudoreps.py` — split script unit tests (6 tests)

Creates a tiny SAM fixture and runs `split_pseudoreps.py` against it.

14. **pr1 + pr2 count equals original:** uses `samtools view -c` to verify
    `count(pr1) + count(pr2) == count(input)`.
15. **No read appears in both outputs:** extracts read names from pr1 and pr2,
    verifies intersection is empty.
16. **PE mates stay together:** fixture includes one PE pair; both mates must
    be in the same output BAM.
17. **Both BAM indexes produced:** verify `.bam.bai` files exist and are
    non-empty.
18. **Same seed deterministic:** run split twice with seed=42, verify output
    BAMs are identical (same read counts per pseudorep).
19. **Different seed may differ:** run with seed=42 and seed=43, verify output
    read counts differ (choose fixture and seeds that guarantee difference;
    skip this test if no such fixture can be guaranteed).

## Files Changed

| File | Change |
|------|--------|
| `config/config.yaml` | Add `seed: 42` to `idr` block |
| `scripts/validate_samples.py` | Add `seed` validation, unknown-key check in `_validate_idr_settings()` |
| `scripts/split_pseudoreps.py` | **New** — deterministic BAM-splitting (hashlib.sha256, samtools view -c) |
| `scripts/stage5b_summary.py` | **New** — reproducibility QC summary helper |
| `workflow/Snakefile` | Add `IDR_SEED`, precompute expansion lists, rule-all targets |
| `workflow/rules/idr.smk` | Add `split_pseudoreps`, `macs3_idr_pseudorep`, `idr_self_pseudoreps`, `idr_pooled_pseudoreps`, `stage5b_summary`, helpers |
| `workflow/schemas/config.schema.yaml` | Document `idr.seed` |
| `README.md` | Update Stage 5 section |
| `KNOWN_ISSUES.md` | Mark 5b items as completed |
| `test/test_stage5b_stress.py` | **New** — 13 DAG/validation tests |
| `test/test_split_pseudoreps.py` | **New** — 6 unit tests |

## Verification

```bash
python3 scripts/validate_samples.py --config config/config.yaml
conda run -n chipseq snakemake -s workflow/Snakefile --configfile config/config.yaml -n --quiet
python3 test/test_validation_stress.py
python3 test/test_stage4b_stress.py
python3 test/test_stage4c_stress.py
python3 test/test_stage5a_stress.py
python3 test/test_stage5b_stress.py
python3 test/test_split_pseudoreps.py
```
