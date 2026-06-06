# Stage 44: MNase paths.smk Limited Pilot — Design Spec

**Date:** 2026-06-06
**Branch:** stage44-mnase-paths-pilot
**Status:** design approved, plan written, awaiting implementation

## Purpose

Create a zero-behavior-change MNase-only path helper file (`paths.smk`). This is a thin pilot to centralize MNase path construction, not artifactization.

## Scope

### In scope

1. `workflow/rules/paths.smk` — 7 thin helper functions, zero `rule` declarations.
2. Replace inline MNase path f-strings in three consumer files:
   - `workflow/rules/mnase.smk` — rule outputs and inputs (9 rules touched)
   - `workflow/rules/targets.smk` — `_mnase_targets()` `expand()` calls
   - `workflow/rules/report.smk` — `pipeline_done` MNase conditional inputs
3. Include `paths.smk` in `workflow/Snakefile` before its consumers, preserving
   `rule all` as the first rule-bearing declaration.

### Out of scope

- Path format changes of any kind.
- ChIP-seq / CUT&Tag / ATAC path migration.
- `make_manifest.py` changes (standalone Python, cannot import `.smk`).
- `output-contract.md` path changes.
- Artifact inventory generation.
- `Artifact` dataclass, `artifact_path()`, `AssayPolicy` YAML, global target resolver.
- `rule all` rewrite.
- Validation or enum enforcement of `fragment_class` or `kind` parameters.
- Log file path helpers.

## Design

### `workflow/rules/paths.smk`

Seven helpers, no imports, no side effects, no `rule` declarations. References global `OUTDIR` (defined in `workflow/Snakefile` before this file is included).

```python
def mnase_fragment_bam(sample, fragment_class):
    return f"{OUTDIR}/{sample}/03_fragments/{sample}.{fragment_class}.bam"

def mnase_fragment_bai(sample, fragment_class):
    return f"{OUTDIR}/{sample}/03_fragments/{sample}.{fragment_class}.bam.bai"

def mnase_signal_bw(sample, kind):
    return f"{OUTDIR}/{sample}/04_signal/{sample}.{kind}.CPM.bw"

def mnase_qc_summary_tsv(sample):
    return f"{OUTDIR}/{sample}/01_qc/{sample}.mnase_qc_summary.tsv"

def mnase_pooled_fragment_bam(experiment, fragment_class):
    return f"{OUTDIR}/experiments/{experiment}/03_fragments/{experiment}.pooled.{fragment_class}.bam"

def mnase_pooled_fragment_bai(experiment, fragment_class):
    return f"{OUTDIR}/experiments/{experiment}/03_fragments/{experiment}.pooled.{fragment_class}.bam.bai"

def mnase_pooled_signal_bw(experiment, kind):
    return f"{OUTDIR}/experiments/{experiment}/04_signal/{experiment}.pooled.{kind}.CPM.bw"
```

**Single `sample`/`experiment` parameter** works in all three consumer contexts:

| Context | Call site | Return value |
|:---|:---|:---|
| Rule output/input | `mnase_fragment_bam("{sample}", "mono")` | `results/{sample}/03_fragments/{sample}.mono.bam` (Snakemake wildcard) |
| `expand()` template | `expand(mnase_fragment_bam("{sample}", "mono"), sample=ids)` | same wildcard template, filled by `expand()` |
| Lambda input function | `mnase_fragment_bam(wc.sample, "mono")` | `results/AC_1/03_fragments/AC_1.mono.bam` (concrete) |

**Not covered by helpers** (intentionally):

- `final.bam` / `final.bai` — shared preprocessing, not MNase-specific
- `pipeline_done` sentinel — lives in `common.smk` pattern
- Picard `insert_size_metrics` param path — not MNase-specific
- Log file paths — peripheral, not structural outputs

### Snakefile Include Order

`paths.smk` is inserted after assay policy includes and before `metadata.smk`:

```
include: "rules/chipseq.smk"
include: "rules/cuttag.smk"
include: "rules/atac.smk"
include: "rules/paths.smk"          <-- NEW
include: "rules/metadata.smk"
include: "rules/targets.smk"        <-- consumer
...
rule all:                             <-- still first rule-bearing decl
...
include: "rules/common.smk"
include: "rules/peaks.smk"
include: "rules/mnase.smk"          <-- consumer
include: "rules/replicates.smk"
include: "rules/idr.smk"
include: "rules/qc.smk"
include: "rules/report.smk"         <-- consumer
```

Rationale:

- `OUTDIR` is defined before any includes — `paths.smk` sees it.
- `paths.smk` loads before `targets.smk`, `mnase.smk`, and `report.smk`.
- Zero `rule` declarations in `paths.smk` means `rule all` is unchanged as first rule-bearing declaration.

### Consumer Changes

#### `workflow/rules/mnase.smk` — 9 rules touched

| Rule | Paths replaced with helpers |
|:---|:---|
| `mnase_split_mono` | output bam → `mnase_fragment_bam`, bai → `mnase_fragment_bai` |
| `mnase_split_sub` | output bam → `mnase_fragment_bam`, bai → `mnase_fragment_bai` |
| `mnase_split_di` | output bam → `mnase_fragment_bam`, bai → `mnase_fragment_bai` |
| `mnase_dyad_bigwig` | output bw → `mnase_signal_bw` |
| `mnase_mono_bigwig` | output bw → `mnase_signal_bw` |
| `mnase_qc_summary` | output tsv → `mnase_qc_summary_tsv`; fragment bam inputs → `mnase_fragment_bam`; signal bw inputs → `mnase_signal_bw` |
| `mnase_pooled_mono` | output bam → `mnase_pooled_fragment_bam`, bai → `mnase_pooled_fragment_bai` |
| `mnase_pooled_dyad_bigwig` | output bw → `mnase_pooled_signal_bw` |
| `mnase_pooled_mono_bigwig` | output bw → `mnase_pooled_signal_bw` |

Not changed: `final.bam`/`final.bai` inputs, log file paths, `insert_size_metrics` param path.

#### `workflow/rules/targets.smk` — `_mnase_targets()`

Replace `expand()` template strings with helper calls. Remove unused `outdir=OUTDIR` keyword arguments (the helper already embeds global `OUTDIR`).

Before:
```python
expand(
    f"{OUTDIR}/{{sample}}/03_fragments/{{sample}}.mono.bam",
    sample=MNASE_SAMPLE_IDS
)
```

After:
```python
expand(
    mnase_fragment_bam("{sample}", "mono"),
    sample=MNASE_SAMPLE_IDS
)
```

Nine sample-level and four experiment-level `expand()` calls follow the same pattern.

#### `workflow/rules/report.smk` — `pipeline_done` MNase inputs

Replace inline f-strings in lambda input functions. Preserve the current lambda return shape: `helper(wc.sample, ...) if _is_mnase(wc) else []` (no list wrapping on the helper result).

Before:
```python
mnase_mono_bam = lambda wc: (
    [f"{OUTDIR}/{wc.sample}/03_fragments/{wc.sample}.mono.bam"]
    if _is_mnase(wc) else []
),
```

After:
```python
mnase_mono_bam = lambda wc: (
    mnase_fragment_bam(wc.sample, "mono") if _is_mnase(wc) else []
),
```

### `make_manifest.py` — intentionally untouched

`make_manifest.py` is standalone Python. It imports from `validate_samples.py` but cannot and should not import `.smk` files. It continues to use its own inline path construction. Duplicating MNase path logic in `make_manifest.py` is an explicit non-goal for Stage 44; that would be part of a future artifactization stage where paths are stored in a language-agnostic format (e.g., YAML).

## Verification

### Hard acceptance tests (all must pass)

```bash
# 1. paths.smk contains no rule declarations
rg -n '^rule ' workflow/rules/paths.smk
# Expected: no output

# 2. rule all remains first rule-bearing declaration
rg -n '^rule ' workflow/Snakefile workflow/rules/paths.smk \
  workflow/rules/metadata.smk workflow/rules/targets.smk
# Expected: first rule match is in Snakefile at the rule all line

# 3. git diff clean
git diff --check

# 4-10. Existing test suite unchanged
python3 test/test_stage43_artifact_inventory.py
python3 test/test_validation_stress.py
python3 test/test_stage28_release_readiness.py
python3 test/test_no_hardcoded_paths.py
SNAKEMAKE=/home/irenadler/miniconda3/envs/chipseq/bin/snakemake \
  python3 test/test_stage39_mnase_stress.py
SNAKEMAKE=/home/irenadler/miniconda3/envs/chipseq/bin/snakemake \
  python3 test/test_stage8_smoke_profiles.py
SNAKEMAKE=/home/irenadler/miniconda3/envs/chipseq/bin/snakemake \
  python3 test/test_stage25_manifest_stress.py
```

### Optional manual check

A before/after dry-run DAG diff can confirm byte-identical path output, but is optional because it requires capturing the "before" snapshot in the same branch before implementation:

```bash
XDG_CACHE_HOME=/tmp/chipseq-snakemake-cache \
  /home/irenadler/miniconda3/envs/chipseq/bin/snakemake \
  -s workflow/Snakefile --configfile config/config.yaml -n --quiet \
  2>&1 | sort > /tmp/dag_before.txt

# ... implement changes ...

XDG_CACHE_HOME=/tmp/chipseq-snakemake-cache \
  /home/irenadler/miniconda3/envs/chipseq/bin/snakemake \
  -s workflow/Snakefile --configfile config/config.yaml -n --quiet \
  2>&1 | sort > /tmp/dag_after.txt

diff /tmp/dag_before.txt /tmp/dag_after.txt
# Expected: no differences
```

### Zero-behavior-change rationale

Every replacement is a mechanical 1:1 substitution. Each helper is a single f-string that produces the identical string as the inline version it replaces. The helpers use the same `OUTDIR` global, same string interpolation, and same directory structure. A Snakemake dry-run before and after produces identical DAG edge sets.

| Risk | Mitigation |
|:---|:---|
| Helper produces wrong path | 1-line f-string, mechanically identical to inline version |
| Snakemake wildcard handling differs | `helper("{sample}", cls)` returns same `{sample}` wildcard string |
| `expand()` placeholder mismatch | `{sample}` in helper return is the same `{sample}` that `expand()` fills |
| Lambda input behavior changes | `helper(wc.sample, cls)` resolves same way as `f"{OUTDIR}/{wc.sample}/..."` |
| `OUTDIR` not accessible | Included after `OUTDIR` defined, shared Snakemake namespace |
