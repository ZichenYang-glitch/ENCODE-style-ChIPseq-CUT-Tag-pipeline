# Stage 30: Strict Input Validation / Runtime Robustness — Design Spec

**Date:** 2026-05-24
**Status:** implemented (Stage 30 completed; all 8 tests pass; backward-compatible)
**Scope:** Planning for stricter runtime validation — optional strict mode, FASTQ/Bowtie2 index existence checks, MACS3 fragment fallback warnings
**Excluded:** Mandatory strict mode, DAG changes, new config schema keys, large data files

---

## 1. Goals

1. Add an optional `--strict-inputs` flag to `scripts/validate_samples.py` that checks real file existence for FASTQ and Bowtie2 index paths.
2. Keep the default non-strict behavior unchanged — all 7 smoke profiles must still pass.
3. Add a clear warning when the MACS3 fragment-size fallback (200bp) is used in SE ChIP-seq.
4. Add tests covering strict-mode rejections and non-strict backward compatibility.
5. No DAG changes, no new config keys, no mandatory behavior changes.

---

## 2. Design Questions

### Q1: Strict input validation interface

**Decision: CLI flag only — `--strict-inputs`**

| Option | Pros | Cons |
| :--- | :--- | :--- |
| `--strict-inputs` CLI flag | Simple, no config schema changes, explicit per-invocation | User must remember to use it |
| `validation.strict_inputs: true` config key | Persistent, checked automatically | Adds config schema complexity, config-killing dry-run compatibility |
| Both, with CLI override | Most flexible | Over-engineered for v0.2 scope |

**Recommendation: `--strict-inputs` CLI flag only.** Default remains non-strict so dry-run smoke profiles continue to work with placeholder paths. Users add `--strict-inputs` when doing real-data validation or release checks.

### Q2: FASTQ validation

**Strict mode behavior:**
- `fastq_1` must exist as a regular file
- For PE samples, `fastq_2` must exist as a regular file
- Support `.fq`, `.fastq`, `.fq.gz`, `.fastq.gz` extensions (no additional check beyond `os.path.isfile`)
- Non-strict mode: no change from current behavior (paths are accepted as-is)

**Non-goal:** The validator does NOT check for truncated/corrupt FASTQ files or validate FASTQ format. File existence only.

### Q3: Bowtie2 index validation

**Strict mode behavior:**
- Accept standard index files: `{prefix}.1.bt2`, `{prefix}.2.bt2`, `{prefix}.3.bt2`, `{prefix}.4.bt2`, `{prefix}.rev.1.bt2`, `{prefix}.rev.2.bt2`
- Accept large-index format: `{prefix}.1.bt2l`, etc.
- If none of the expected `.1.bt2` / `.1.bt2l` files exist, produce a clear error message listing which files were expected
- Non-strict mode: no change from current behavior

**Error message example:**
```
ERROR: Bowtie2 index not found for sample S1: /path/to/index
  Expected .1.bt2 or .1.bt2l at: /path/to/index.1.bt2
```

### Q4: Existing resource validation — do not duplicate

Current `validate_samples.py` already validates:
- `control_bam` existence when `use_control: true` (line ~1085)
- `genome_resources.<genome>.chrom_sizes/blacklist/gtf/reference_fasta` existence when non-empty (lines 460-466)

**Decision:** `--strict-inputs` does NOT re-validate these. They are already checked regardless of strict mode. The `--strict-inputs` flag only adds FASTQ + Bowtie2 index checks.

### Q5: MACS3 fragment-size fallback policy

**Current behavior** (lines 71-90 of `chipseq.smk`):
- SE ChIP-seq `extend_reads=auto/yes` calls `_extract_macs3_fragment_size()`
- Tries to read MACS3 log; if unavailable or invalid, returns `None`
- Caller falls back to `200` (default SE extension)
- Exception is silently caught — user never knows fallback happened

**Design improvement:**
- Keep the fallback — do NOT break runs where MACS3 log is genuinely unavailable
- Emit a clear warning to stderr when fallback is used: `"Warning: MACS3 fragment size unavailable for <sample>, using --extendReads 200"`
- Warning should NOT fail the run
- Future (Stage 30+): optionally surface fallback status in QC summary as an informational column

**Not in scope:** Making the fallback configurable, or adding a config key for default extension.

### Q6: Test plan

**New tests in `test/test_stage30_strict_inputs_stress.py`:**

| # | Test | Expected |
|---|------|----------|
| 1 | `--strict-inputs` with missing `fastq_1` | Validation fails, "file not found" for fastq_1 |
| 2 | `--strict-inputs` with PE missing `fastq_2` | Validation fails, "file not found" for fastq_2 |
| 3 | `--strict-inputs` with missing Bowtie2 `.1.bt2` | Validation fails, "index not found" |
| 4 | `--strict-inputs` with `.bt2` set present | Validation passes |
| 5 | `--strict-inputs` with `.bt2l` set present | Validation passes |
| 6 | Non-strict mode with placeholder paths | Smoke profile behavior unchanged (7/7) |
| 7 | MACS3 fragment fallback emits warning to stderr | Warning contains "extendReads 200" |
| 8 | `--strict-inputs` accepts `.fq.gz` extension | Validation passes if file exists |

**Existing tests must still pass:**
- `test_validation_stress.py` (15/15) — non-strict, should be unchanged
- `test_stage8_smoke_profiles.py` (7/7) — smoke profiles unchanged
- `test_no_hardcoded_paths.py` — unchanged

### Q7: Documentation plan

| File | Updates |
| :--- | :--- |
| `README.md` | Add `--strict-inputs` to validation section |
| `docs/configuration.md` | Document strict mode as optional pre-run check |
| `RELEASE_CHECKLIST.md` | Add `--strict-inputs` as optional real-data validation step |
| `docs/assay-policy.md` | Note MACS3 fragment fallback warning in ChIP-seq SE extension section |

---

## 3. Architecture

### CLI change

```
python3 scripts/validate_samples.py --config config/config.yaml           # non-strict (default)
python3 scripts/validate_samples.py --config config/config.yaml --strict-inputs  # strict mode
```

### validate_samples.py additions

```python
def _check_fastq_exists(path, sample_id, label):
    """In strict mode, verify a FASTQ path exists."""

def _check_bowtie2_index(prefix, sample_id):
    """In strict mode, verify Bowtie2 index files exist."""

def validate_input_strictness(samples, config, strict_inputs):
    """If strict_inputs, validate FASTQ and Bowtie2 index existence."""
```

### MACS3 fallback warning

In `chipseq.smk` `_extract_macs3_fragment_size()`, replace silent `except Exception: pass` with:

```python
except Exception:
    return None

# In get_extend_reads_chipseq, when fallback is used:
print(f"Warning: MACS3 fragment size unavailable for {sample}, "
      f"using --extendReads 200", file=sys.stderr)
```

---

## 4. Out of Scope

- Making strict mode the default
- Config key for strict validation
- FASTQ format/content validation (beyond existence)
- Checking Bowtie2 index internal integrity
- GC bias, adaptor content, or sequence quality validation
- Network-accessible path validation
- Pipeline aborts on MACS3 fallback
- New conda environment dependencies
