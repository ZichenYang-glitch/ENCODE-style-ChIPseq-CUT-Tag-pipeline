# Stage 54: Consensus Engine — Design Spec

**Date:** 2026-06-18
**Status:** Design approved; implementation pending
**Author:** YangZiChen-glitch (Kaslana)
**Depends on:** Stage 53 (reproducibility policy)
**Blocks:** Stage 55 (ATAC narrow IDR)

## 1. Overview

Stage 54 implements the consensus/overlap engine: a standalone Python script
that reads per-biorep peak files and produces replicate-validated consensus
peak outputs. No DAG integration, no Snakemake rules, no target expansion.

### 1.1 Motivation

Stage 53 defined the reproducibility policy: consensus is N-of-M replicate
support baseline for all six peak-calling modes. Stage 54 implements the
algorithm that computes consensus peaks from per-replicate peak calls.

### 1.2 Scope

| In scope | Out of scope |
|----------|-------------|
| `scripts/compute_consensus.py` | Snakemake rules |
| `test/test_stage54_consensus.py` | DAG target expansion |
| Pure Python stdlib | BEDTools dependency |
| narrowPeak, broadPeak, BED input | SEACR IDR |
| Synthetic unit tests (31) | Manifest/output-contract integration |

### 1.3 Dependencies

Pure Python stdlib only: `sys`, `argparse`, `os`, `collections`, `math`,
`itertools`. No BEDTools, no Conda environment, no external packages.

Future optimization: a BEDTools backend may be added in a later stage if
performance becomes a real issue. That is not Stage 54 scope.

## 2. Algorithm

### 2.1 Overview

```
INPUT: N peak files + N biorep labels
       ↓
PARSE:  source peaks → (chrom, start, end, biorep, original_fields)
       ↓
GRAPH:  nodes = source peaks
        edges = pairwise reciprocal overlap ≥ threshold, same chrom
       ↓
COMPONENTS: connected components of overlap graph
       ↓
FILTER:  support_count = |distinct bioreps in component|
         keep if support_count ≥ min_replicates
       ↓
OUTPUT:  consensus peaks (merged intervals) + summary TSV
```

### 2.2 Step 1: Parse Peak Files

For each `(peak_file, biorep_label)` pair, parse the file according to
`--format`:

- **narrowPeak:** 10 tab-separated columns.
  `chrom, start, end, name, score, strand, signalValue, pValue, qValue, peak`
- **broadPeak:** 9 tab-separated columns.
  `chrom, start, end, name, score, strand, signalValue, pValue, qValue`
- **bed:** BED3+ columns. Minimum: `chrom, start, end`.
  If name/score/strand absent, synthesize defaults.

Per-line validation:
- `start` and `end` must be non-negative integers, `end > start`.
- `signalValue` parsed as float; if unparseable, treat as 0.0.
- `score` parsed as int; if unparseable, treat as 0.
- `peak` (narrowPeak only) parsed as int; if unparseable or out of bounds,
  treat as invalid (triggers summit fallback).
- Skip header lines starting with `#`, `track`, or `browser`.
  Physical file line numbers are preserved in error messages (headers/blank
  lines count toward the line number).
- Blank/whitespace-only lines are skipped.
- Malformed lines (fewer than required columns, non-integer coordinates,
  `end <= start`) → non-zero exit with physical line number.

Output: list of `(chrom, start, end, biorep_label, fields_dict)` tuples.

### 2.3 Step 2: Build Overlap Graph

Each source peak is a node. An undirected edge exists between peaks A and B
iff all of:

1. `A.chrom == B.chrom`
2. `overlap_len = min(A.end, B.end) - max(A.start, B.start) > 0`
3. `overlap_len / (A.end - A.start) >= reciprocal_overlap`
4. `overlap_len / (B.end - B.start) >= reciprocal_overlap`

Implementation: sort peaks by `(chrom, start)`. For each peak, scan forward
within the same chromosome until `other.start > self.end`.

### 2.4 Step 3: Find Connected Components

Run DFS/BFS over the undirected graph. Each connected component is a candidate
consensus cluster.

**Chained overlaps are intentional.** If peak A overlaps B, and B overlaps C,
but A does not directly overlap C — all three form one connected component.
Support count includes all distinct bioreps in the chain. This behavior is
documented and tested (T9). Future stages may add stricter clique or
anchor-based modes.

Chromosomes are strict separation boundaries — no edges cross chromosomes.

### 2.5 Step 4: Filter by Support Count

For each component:
- Collect all distinct `biorep_label` values.
- `support_count = len(distinct_bioreps)`.
- Keep if `support_count >= min_replicates`.

Multiple peaks from the same biological replicate within a component do NOT
increase `support_count`. A biorep contributes at most 1 to support.

Validation: `min_replicates > n_bioreps` → non-zero exit (config/user error).

### 2.6 Step 5: Produce Consensus Intervals

For each retained component, output a single consensus peak.

| Column | Strategy |
|--------|---------|
| `chrom` | Component chromosome |
| `start` | Minimum start among all peaks in component |
| `end` | Maximum end among all peaks in component |
| `name` | `consensus_peak_N` (1-based, in sorted component order) |
| `score` | `int(1000 * support_count / n_bioreps + 0.5)` (round half up) |
| `strand` | `.` always |
| `signalValue` | Maximum signalValue among peaks in component |
| `pValue` | `-1` (placeholder, not new MACS3/IDR statistic) |
| `qValue` | `-1` (placeholder) |
| `peak` | narrowPeak only: summit derived from max-signalValue source peak (§2.7) |

BED output: 6 columns (chrom, start, end, name, score, strand).

Output sorted by `(chrom, start, end, name)` deterministically. Tie-breaking
on identical coordinates uses `(sorted_biorep_labels, first_source_order)`.

### 2.7 narrowPeak Summit Derivation

1. From the component, find the source peak with maximum `signalValue`.
2. If tie: pick smallest `(chrom, start, end)`, then biorep label, then input
   file order / line number.
3. If selected source peak has a valid summit column (non-negative integer,
   `source_start + source_summit < source_end`):
   - `consensus_summit = source_start + source_summit - consensus_start`.
4. Else (fallback):
   - `consensus_summit = (consensus_end - consensus_start) // 2`.
5. Clamp: `0 <= consensus_summit < consensus_end - consensus_start`.

Note: consensus summit is inherited from the max-signal source peak. It is not
newly computed from read pileup.

### 2.8 Step 6: Write Summary TSV

Columns as specified in Stage 53 §4.5:

```
experiment   assay   peak_mode   caller   n_bioreps   min_replicates
reciprocal_overlap   consensus_peak_count   support_distribution
biorep_labels   source_peak_files   final_method   final_output
```

- `support_distribution`: JSON object, e.g. `{"2": 15000, "3": 8000}`,
  `{}` when no consensus peaks.
- `biorep_labels`: comma-separated sorted list of all `--bioreps` labels.
- `source_peak_files`: JSON array of `--peaks` paths in order.
- `final_method` / `final_output`: from CLI flags (empty string if not provided).

Summary TSV is always written — even when `consensus_peak_count = 0`.

## 3. CLI Interface

```
python3 scripts/compute_consensus.py \
  --peaks <file1> <file2> [...] \
  --bioreps <label1> <label2> [...] \
  --format narrowPeak|broadPeak|bed \
  --min-replicates <int>       (default 2)
  --reciprocal-overlap <float> (default 0.5)
  --output <path> \
  --summary <path> \
  --experiment <str> \
  --assay <str> \
  --caller <str> \
  --peak-mode <str> \
  [--final-method <str>]       (default "")
  [--final-output <str>]       (default "")
```

### Required flags

`--peaks`, `--bioreps`, `--format`, `--output`, `--summary`, `--experiment`,
`--assay`, `--caller`, `--peak-mode`.

### Optional flags with defaults

`--min-replicates` (2), `--reciprocal-overlap` (0.5), `--final-method` (""),
`--final-output` ("").

### Validation

| Check | Behavior on failure |
|-------|-------------------|
| `len(--peaks) < 2` | Exit: "at least two peak files / biological replicates are required" |
| `len(--peaks) != len(--bioreps)` | Exit: "N peak files but M biorep labels" |
| `--bioreps` has duplicates | Exit: "duplicate biorep label: <label>" |
| `--min-replicates < 2` | Exit: "min-replicates must be >= 2" |
| `--min-replicates > len(--bioreps)` | Exit: "min-replicates exceeds n_bioreps" |
| `--reciprocal-overlap <= 0 or > 1` | Exit: "reciprocal-overlap must be in (0, 1]" |
| `--format` not in {narrowPeak, broadPeak, bed} | Exit: "unknown format" |
| Peak file not readable | Exit: "cannot read: <path>" |
| Malformed peak line | Exit: "<file>: line N: <reason>" |

### Exit codes

- 0: success (even when consensus_peak_count = 0)
- non-zero: error (validation failure, malformed input)

### Edge cases

| Scenario | Behavior |
|----------|---------|
| No consensus peaks found | Empty output file (0 bytes), summary with count=0 and `{}` |
| Unsorted input | Output sorted deterministically |
| Header lines (`#`, `track`, `browser`) | Skipped as records, count toward physical line numbers |
| All peaks from same biorep | support_count=1, nothing retained if min_replicates >= 2 |
| Parent directories don't exist | Script creates them for `--output` and `--summary` |

## 4. Output Semantics

### Consensus peak file (`--output`)

- Tab-separated, no header.
- narrowPeak: 10 columns. broadPeak: 9 columns. BED: 6 columns.
- Sorted by `(chrom, start, end, name)`.
- Empty file (0 bytes) when no components retained.
- `strand` is always `.`.
- `score` encodes reproducibility: `int(1000 * support_count / n_bioreps + 0.5)`.
- `pValue` and `qValue` are `-1` (documented placeholders, not computed statistics).
- Summit (narrowPeak column 10): `0 <= summit < consensus_end - consensus_start`.

### Summary TSV (`--summary`)

- Always written with header row + one data row.
- `support_distribution` uses JSON formatting with double quotes.
- `source_peak_files` uses JSON array with double quotes.
- `biorep_labels`: comma-separated, sorted order.
- Always includes all 13 columns from the Stage 53 schema.

## 5. Test Plan

**File:** `test/test_stage54_consensus.py`
**Style:** Standalone `main()`, subprocess invocation, synthetic peak files in
temporary directories. No pytest.

### Core algorithm (T1–T8)

| # | Test | Expected |
|---|------|---------|
| T1 | 2 bioreps, perfect overlap | 1 consensus, support=2, score=1000 |
| T2 | 2 bioreps, no overlap | 0 consensus, empty output |
| T3 | 3 bioreps, N-of-M=2 (2 of 3 overlap) | 1 consensus, support=2, score=667 |
| T4 | 3 bioreps, all overlap | 1 consensus, support=3, score=1000 |
| T5 | reciprocal_overlap fails (500bp peak, 50bp overlap) | 0 consensus |
| T6 | reciprocal_overlap passes (500bp peak, 400bp overlap) | 1 consensus |
| T7 | reciprocal_overlap at exact threshold (0.5) | 1 consensus |
| T8 | Same biorep multiple peaks (within one input file) | support=1, not 2 |

### Graph topology (T9–T12)

| # | Test | Expected |
|---|------|---------|
| T9 | Chained overlap (A-B-C, A not overlapped with C) | 1 component, support=3 |
| T10 | Adjacent intervals (100-200, 200-300, overlap=0) | no consensus |
| T11 | Cross-chromosome separation | never cluster |
| T12 | Unsorted input | output sorted (chrom, start, end, name) |

### CLI validation (T13–T17)

| # | Test | Expected |
|---|------|---------|
| T13 | Mismatched peaks/bioreps count | non-zero exit |
| T14 | Duplicate biorep labels | non-zero exit |
| T15 | `--min-replicates 1` | non-zero exit |
| T16 | `--min-replicates > n_bioreps` | non-zero exit |
| T17 | Only 1 peak file (`< 2`) | non-zero exit, "at least two required" |

### Output format (T18–T21)

| # | Test | Expected |
|---|------|---------|
| T18 | narrowPeak output shape | 10 columns, valid summit |
| T19 | broadPeak output shape | 9 columns, no summit column |
| T20 | BED output shape | 6 columns (BED6) |
| T21 | Score calculation | 2/2→1000, 2/3→667, 3/4→750 |

### Edge cases (T22–T25)

| # | Test | Expected |
|---|------|---------|
| T22 | Empty consensus | output empty, summary count=0, distribution={} |
| T23 | Malformed line (non-integer start) | non-zero exit, physical line number |
| T24 | end <= start | non-zero exit |
| T25 | Summit fallback (invalid source summit) | midpoint used, within [0, end-start) |

### Summary TSV (T26–T27)

| # | Test | Expected |
|---|------|---------|
| T26 | All 13 columns present | header verified |
| T27 | Values correct | experiment, assay, support_distribution, biorep_labels, source_peak_files |

### Additional validation (T28–T31)

| # | Test | Expected |
|---|------|---------|
| T28 | Invalid reciprocal_overlap (0, 1.5) | non-zero exit |
| T29 | Invalid format | non-zero exit |
| T30 | BED3 input accepted | synthesized name/score/strand in BED6 output |
| T31 | Parent directories created | script creates dirs for --output and --summary |

**Total: 31 tests.** All synthetic, no real data, no Snakemake/DAG tests.

## 6. Non-goals

- No DAG integration (no `workflow/rules/*.smk` or `workflow/Snakefile` changes)
- No target expansion, no `rule all` changes
- No manifest/output-contract runtime integration
- No new IDR rules
- No artifact runtime adoption
- No BEDTools dependency (future optimization only)
- No stdin/YAML/JSON config input
- No Snakemake rule wrappers for consensus
- No pipeline-level `06_reproducibility/` DAG outputs
- No SEACR IDR

## 7. Deferred to Stage 55+

- ATAC narrow IDR implementation
- Snakemake rules wrapping `compute_consensus.py`
- Target expansion in `rule all`
- Manifest entries for `06_reproducibility/` outputs
- MultiQC integration

## 8. Verification

```bash
python3 test/test_stage54_consensus.py
python3 test/test_stage53_config_validation.py
python3 test/test_stage53_stage5_invariant.py
python3 test/test_stage53_pooled_not_validated.py
python3 test/test_stage53_experimental_warnings.py
python3 test/test_stage53_reproducibility_policy_contract.py
python3 test/test_stage53_output_path_templates.py
python3 test/test_stage28_release_readiness.py
python3 test/test_no_hardcoded_paths.py
git diff --check
```
