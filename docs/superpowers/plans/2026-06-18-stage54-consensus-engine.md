# Stage 54 Implementation Plan: Consensus Engine

**Date:** 2026-06-18
**Status:** Ready for implementation
**Depends on:** Stage 53 (reproducibility policy)
**Blocks:** Stage 55 (ATAC narrow IDR)
**Design spec:** `docs/superpowers/specs/2026-06-18-stage54-consensus-engine-design.md`

## Context

Stage 53 defined the reproducibility policy. Stage 54 implements the consensus
engine: a standalone Python (stdlib) script that computes consensus peaks from
per-biorep peak calls using N-of-M replicate support with pairwise reciprocal
overlap.

No Snakemake integration. No DAG changes. No BEDTools.

## Files in scope

### To create

| File | Purpose |
|------|---------|
| `scripts/compute_consensus.py` | Consensus engine (~300-400 lines, pure stdlib) |
| `test/test_stage54_consensus.py` | 31 synthetic unit tests |

### Explicitly NOT modified

- `workflow/Snakefile`
- `workflow/rules/*.smk`
- `config/config.yaml`
- `scripts/validate_samples.py`
- `docs/reproducibility-policy.md`
- Any Stage 53 files

## Implementation steps

### Step 1: Implement `scripts/compute_consensus.py`

Core structure:

```
argparse → validate → parse_peak_files → build_graph → find_components →
filter_by_support → write_peaks → write_summary
```

Module organization:

1. `parse_args()` — argparse with all required/optional flags
2. `validate_inputs(args)` — enforce all validation rules from design §3
3. `parse_peak_file(path, biorep, format)` → `list[SourcePeak]`
   - narrowPeak: 10 columns
   - broadPeak: 9 columns
   - bed: BED3+ columns
   - Returns `(chrom, start, end, biorep, fields_dict)`
4. `build_overlap_graph(peaks, reciprocal_overlap)` → adjacency dict
   - Pairwise reciprocal overlap check
   - Sorted scan within chromosome
5. `find_components(graph)` → `list[list[int]]` (node index lists)
   - DFS/BFS
6. `filter_by_support(components, min_replicates)` → filtered components
7. `compute_consensus_peak(component, n_bioreps, format)` → output row
   - score = `int(1000 * support_count / n_bioreps + 0.5)`
   - strand = `.`
   - signalValue = max source signal
   - pValue/qValue = `-1`
   - summit from max-signalValue source peak with midpoint fallback
8. `write_summary_tsv(...)` — 13-column TSV
9. `main()` — orchestrator, exit codes

### Step 2: Implement `test/test_stage54_consensus.py`

Standalone `main()` style. All tests use `subprocess.run()` to invoke
`scripts/compute_consensus.py` with synthetic peak files in `tempfile`
directories.

Test structure:
```
main()
  → check("T1: perfect overlap", t1)
  → check("T2: no overlap", t2)
  ...
  → print summary
  → sys.exit(0 if all_pass else 1)
```

Helper functions:
- `make_narrowpeak(path, rows)` — write synthetic narrowPeak file
- `make_broadpeak(path, rows)` — write synthetic broadPeak file
- `make_bed(path, rows, ncols)` — write synthetic BED file
- `run_consensus(args)` — subprocess wrapper
- `assert_file_has_n_columns(path, n)` — count columns
- `assert_file_has_n_lines(path, n)` — count data lines

### Step 3: Run tests and verify

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

## Key design decisions (from spec)

1. **Pairwise source-peak reciprocal overlap**, not candidate-region overlap
2. **Connected components** via graph DFS, not merge-first or anchor iteration
3. **Chained overlaps intentional** (A-B-C with A-/-C forms one component)
4. **Score = reproducibility fraction**: `int(1000 * support_count / n_bioreps + 0.5)`
5. **Strand always `.`**
6. **pValue/qValue always `-1`** (placeholders)
7. **Summit from max-signalValue** source peak with midpoint fallback
8. **BED3+ input accepted**, output always BED6
9. **Duplicate biorep labels rejected** at CLI level
10. **min_replicates > n_bioreps rejected** at CLI level

## Non-goals (Stage 54)

- No DAG integration
- No Snakemake rules
- No target expansion
- No manifest/output-contract changes
- No new IDR rules
- No artifact runtime adoption
- No BEDTools dependency
- No SEACR IDR

## Verification

```bash
cd /home/irenadler/workflow/chipseq
python3 test/test_stage54_consensus.py
```
