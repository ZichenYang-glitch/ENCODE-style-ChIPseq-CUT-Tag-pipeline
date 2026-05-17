# Stage 7a: CUT&Tag Fragment-Size QC Foundation

## Scope

Stage 7a adds CUT&Tag-specific fragment-size QC as the lowest-risk first slice
of the Stage 7 CUT&Tag branch. It computes fragment-length distributions from
`final.bam` for samples with `assay: cuttag`.

**Stage 7a includes:**
- `scripts/calc_cuttag_fragment_size.py` — stdlib-only QC helper
- `cuttag_fragment_size` rule in `qc.smk`
- `qc.cuttag_fragment_size: true` config switch (default on)
- Output: `results/<sample>/01_qc/<sample>.cuttag_fragment_size.tsv`

**Scheduling:** The rule is scheduled for all active samples with
`assay == "cuttag"` (using `ACTIVE_SAMPLE_IDS` filtered to cuttag). When
`use_control: true`, this includes FASTQ-based CUT&Tag control rows so
their fragment-size distributions are also available for QC.

**Deferred to Stage 7b+:**
- SEACR/GoPeaks peak callers
- Tn5 insertion-aware metrics
- Spike-in normalization
- Integration with per-sample `qc_summary.tsv`

## Configuration

```yaml
qc:
  cuttag_fragment_size: true   # default on; only schedules for assay=cuttag
```

New key in the existing `qc` block. Default `true`. Accepts Python `True`/`False`
or string `"true"`/`"false"`.

## Validation Changes

### `_validate_qc_config()` in `validate_samples.py`

Add `cuttag_fragment_size` to the fixed-key set so it is validated (not silently
dropped) and appears in the normalized output:

```python
def _validate_qc_config(qc: dict) -> dict:
    ...
    return {
        "blacklist_filter": _normalize_bool("blacklist_filter"),
        "frip": _normalize_bool("frip"),
        "library_complexity": _normalize_bool("library_complexity"),
        "nrf_pbc": _normalize_bool("nrf_pbc"),
        "signal_tracks": _normalize_bool("signal_tracks"),
        "summary": _normalize_bool("summary"),
        "cuttag_fragment_size": _normalize_bool("cuttag_fragment_size"),  # NEW
    }
```

If `cuttag_fragment_size` is absent from the `qc` block, `qc.get("cuttag_fragment_size", True)`
defaults it to `True`.

### `config.schema.yaml`

Add `cuttag_fragment_size` to the `qc` schema section (alongside other boolean switches):

```yaml
cuttag_fragment_size:
  type: boolean
  default: true
  description: >
    Enable CUT&Tag fragment-size QC. Only schedules for samples with
    assay=cuttag.
```

## `scripts/calc_cuttag_fragment_size.py`

### Interface

```
python3 scripts/calc_cuttag_fragment_size.py \
    --sample <sample_id> \
    --bam <final.bam> \
    --layout PE|SE \
    --output <output.tsv>
```

### Behavior — PE samples

1. Stream `samtools view <bam>` (no `-h`).
2. Parse each SAM record.
3. Include only records where:
   - Flag `0x1` (paired) is set
   - Flag `0x2` (properly paired) is set
   - Flag `0x4` (unmapped) is NOT set
   - Flag `0x100` (secondary) is NOT set
   - Flag `0x800` (supplementary) is NOT set
   - Flag `0x40` (first in pair / read1) IS set — count only read1 to avoid
     double counting each fragment
4. Fragment length = absolute value of TLEN (column 9).
5. Compute:
   - `fragment_count` — total number of read1 properly paired primary alignments
   - `fragment_mean` — mean fragment length (float, 1 decimal place)
   - `fragment_median` — median fragment length (float, 1 decimal place; .5 for even counts)
   - `fragment_mode` — mode fragment length (int, smallest value in case of tie)
   - `fragment_min` — minimum fragment length (int)
   - `fragment_max` — maximum fragment length (int)
   - Four **mutually exclusive** size bins that partition all fragments:
     - `fraction_lt_150` — length < 150
     - `fraction_150_300` — 150 ≤ length < 300
     - `fraction_300_500` — 300 ≤ length < 500
     - `fraction_ge_500` — length ≥ 500
   - One extra **diagnostic** bin (NOT part of the sum-to-1 group):
     - `fraction_lt_120` — length < 120 (sub-nucleosomal / TF binding range)

   All fractions are computed as: `count_in_bin / fragment_count`
   and formatted to 4 decimal places. The four mutually exclusive bins
   must sum to approximately 1.0 (tolerance ≤ 0.0002, due to rounding).
   `fraction_lt_120` overlaps `fraction_lt_150`; it is a diagnostic subset
   reported independently.

6. Write TSV with header + single data row.

### Behavior — SE samples

1. Write the same TSV header.
2. Write a single data row with sample ID, `layout_not_supported` for the
   status field, and `NA` for all numeric columns.
3. Do NOT fail or exit non-zero.

### Behavior — PE with zero fragments

If after filtering there are zero fragments:
- `fragment_count` = 0
- All numeric/stat columns = `NA`
- `status` = `no_fragments`

### TSV Columns (14)

| Column | PE value | SE value | Zero-fragment value |
|--------|----------|----------|---------------------|
| `sample` | sample ID | sample ID | sample ID |
| `layout` | `PE` | `SE` | `PE` |
| `fragment_count` | integer | `NA` | `0` |
| `fragment_mean` | float (1 dp) | `NA` | `NA` |
| `fragment_median` | float (1 dp) | `NA` | `NA` |
| `fragment_mode` | integer | `NA` | `NA` |
| `fragment_min` | integer | `NA` | `NA` |
| `fragment_max` | integer | `NA` | `NA` |
| `fraction_lt_150` | float (4 dp) | `NA` | `NA` |
| `fraction_150_300` | float (4 dp) | `NA` | `NA` |
| `fraction_300_500` | float (4 dp) | `NA` | `NA` |
| `fraction_ge_500` | float (4 dp) | `NA` | `NA` |
| `fraction_lt_120` | float (4 dp) | `NA` | `NA` |
| `status` | `ok` | `layout_not_supported` | `no_fragments` |

14 columns total.

### Dependencies

The Python script uses stdlib only (`argparse`, `statistics`, `subprocess`).
It requires `samtools` on `PATH` for `samtools view` streaming, which is provided
by the workflow Conda environment (`workflow/envs/chipseq.yml`). Mode is computed
manually (frequency table, smallest value on tie).

### Error handling

- `samtools` not found → print error, exit 1
- `samtools view` fails → print error, exit 1
- Output file cannot be written → let Python raise; Snakemake captures it

## Snakemake Rule

### Placement

`workflow/rules/qc.smk` — added after the existing `nrf_pbc` rule.

### Rule definition

```python
rule cuttag_fragment_size:
    output:
        f"{OUTDIR}/{{sample}}/01_qc/{{sample}}.cuttag_fragment_size.tsv",
    input:
        f"{OUTDIR}/{{sample}}/02_align/{{sample}}.final.bam",
    params:
        layout = lambda wc: SAMPLE_MAP[wc.sample]["layout"],
    log:
        f"{OUTDIR}/{{sample}}/logs/{{sample}}.cuttag_fragment_size.log",
    conda:
        "../envs/chipseq.yml",
    shell:
        """
        set -e -o pipefail
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        python3 scripts/calc_cuttag_fragment_size.py \
            --sample {wildcards.sample:q} \
            --bam {input:q} \
            --layout {params.layout:q} \
            --output {output:q} \
            2>&1 | tee {log:q}
        """
```

## Config and Scheduling

### QC_CONFIG addition (Snakefile)

```python
QC_CONFIG = VALIDATED_CONFIG.get("qc", {
    "blacklist_filter": True,
    "frip": True,
    "library_complexity": True,
    "nrf_pbc": True,
    "signal_tracks": True,
    "summary": True,
    "cuttag_fragment_size": True,   # ← NEW
})
```

### config.yaml addition

```yaml
qc:
  cuttag_fragment_size: true
```

### Precomputed sample list (Snakefile)

```python
CUTTAG_ACTIVE_SAMPLE_IDS = [
    sid for sid in ACTIVE_SAMPLE_IDS
    if SAMPLE_MAP[sid]["assay"] == "cuttag"
]
```

`ACTIVE_SAMPLE_IDS` already includes control rows when `use_control: true`,
so FASTQ-based CUT&Tag control samples also get fragment-size QC.

### Rule-all target

```python
cuttag_frag_size = (
    expand(
        "{outdir}/{sample}/01_qc/{sample}.cuttag_fragment_size.tsv",
        outdir=OUTDIR,
        sample=CUTTAG_ACTIVE_SAMPLE_IDS,
    )
    if QC_CONFIG.get("cuttag_fragment_size", True)
    and CUTTAG_ACTIVE_SAMPLE_IDS
    else []
),
```

### Backward compatibility

- `qc.cuttag_fragment_size: false` → validated and respected; rule not scheduled
- No CUT&Tag samples → `CUTTAG_ACTIVE_SAMPLE_IDS` is empty, rule not scheduled
- No `qc:` block → defaults to `true`, but only CUT&Tag samples trigger it
- ChIP-seq behavior unchanged — chipseq samples are excluded from `CUTTAG_ACTIVE_SAMPLE_IDS`
- `_validate_qc_config` now validates the key; `qc: {cuttag_fragment_size: false}` is
  correctly normalized to `False`, not silently dropped

## Test Plan

### `test/test_calc_cuttag_fragment_size.py` — helper unit tests (6 tests)

1. **PE fixture** — creates a tiny SAM with 5 properly paired read1/read2
   pairs (known TLENs). Runs script with `--layout PE`, verifies:
   - fragment_count == 5 (only read1 counted)
   - mean, median (float, .5 OK), mode, min, max match expected
   - the four mutually exclusive bin fractions sum to ~1.0 (tolerance ≤ 0.0002)
   - fraction_lt_120 is NOT included in the sum (it is a diagnostic subset)

2. **SE fixture** — runs script with `--layout SE`, verifies:
   - `status == layout_not_supported`
   - all numeric columns are `NA`

3. **Zero fragments** — SAM with only unmapped reads, verifies:
   - fragment_count == 0
   - status == `no_fragments`
   - all numeric columns are `NA`

4. **Mixed records** — SAM with mapped + unmapped + secondary reads;
   verifies only properly paired read1 records are counted

5. **Read2 exclusion** — read1 and read2 have different TLEN signs;
   verifies read2 records are excluded, fragment count correct

6. **Exclusive bin sum** — larger fixture with fragments across all bins;
   verifies the four mutually exclusive bin fractions sum to ~1.0
   (tolerance ≤ 0.0002), and fraction_lt_120 ≤ fraction_lt_150 (subset
   relationship)

### `test/test_stage7a_stress.py` — Snakemake scheduling tests (10 tests)

7. **Default config + ChIP-seq sample** — `cuttag_fragment_size` NOT in DAG
8. **Default config + CUT&Tag PE sample** — rule IS in DAG
9. **Default config + CUT&Tag SE sample** — rule IS in DAG (produces NA output)
10. **`qc.cuttag_fragment_size: false`** — validated and rule NOT in DAG
11. **`qc: {}` (empty block)** — defaults to true, rule IS in DAG for CUT&Tag
12. **`--list-rules`** — shows `cuttag_fragment_size` rule
13. **`-n -p` dry-run** — `calc_cuttag_fragment_size.py` visible with
    `--layout PE` and `--sample` args
14. **Multi-assay sheet** — ChIP-seq excluded from scheduling, CUT&Tag included
15. **`use_control: true` + CUT&Tag control row** — control row also scheduled
    for fragment-size QC (via `CUTTAG_ACTIVE_SAMPLE_IDS`)
16. **`qc.cuttag_fragment_size: maybe`** — rejected by `_normalize_bool`

## Files Changed

| File | Change |
|------|--------|
| `scripts/calc_cuttag_fragment_size.py` | **New** — stdlib QC helper |
| `scripts/validate_samples.py` | Add `cuttag_fragment_size` to `_validate_qc_config` return dict |
| `workflow/rules/qc.smk` | Add `cuttag_fragment_size` rule |
| `workflow/Snakefile` | Add to `QC_CONFIG` defaults, precompute `CUTTAG_ACTIVE_SAMPLE_IDS`, add rule-all target |
| `config/config.yaml` | Add `cuttag_fragment_size: true` to `qc` block |
| `workflow/schemas/config.schema.yaml` | Document `cuttag_fragment_size` in qc schema |
| `test/test_calc_cuttag_fragment_size.py` | **New** — 6 unit tests |
| `test/test_stage7a_stress.py` | **New** — 10 DAG tests |

## Verification

```bash
python3 scripts/validate_samples.py --config config/config.yaml
conda run -n chipseq snakemake -s workflow/Snakefile --configfile config/config.yaml -n --quiet
python3 test/test_validation_stress.py
python3 test/test_calc_cuttag_fragment_size.py
python3 test/test_stage7a_stress.py
```
