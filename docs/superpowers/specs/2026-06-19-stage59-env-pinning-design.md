# Stage 59: Conda Env Version Pinning — Design Spec

**Date:** 2026-06-19
**Status:** Implemented
**Author:** YangZiChen-glitch (Kaslana)
**Scope:** P0 hardening — pin core tool versions in Conda env files

## 1. Purpose

Pin version constraints on all top-level core CLI tools in workflow/envs/*.yml
to reduce reproducibility drift. Without pins, `conda create` / `micromamba
create` can resolve different tool versions on different dates, producing
subtly different results.

## 2. Strategy

- `>=major.minor,<next-major` ranges for mature tools
- Pin only top-level CLI tools, not transitive dependencies
- Do not generate lockfiles
- Do not change Python version constraints

## 3. Pins applied

| Tool | Constraint |
|------|-----------|
| fastqc | `>=0.12,<1` |
| trim-galore | `>=0.6,<1` |
| cutadapt | `>=4.0,<5` |
| bowtie2 | `>=2.5,<3` |
| samtools | `>=1.16,<2` |
| bedtools | `>=2.31,<3` |
| macs3 | `>=3.0,<4` |
| deeptools | `>=3.5,<4` |
| idr | `>=2.0,<3` |
| multiqc | `>=1.20,<2` |
| picard | `>=3.0,<4` |
| phantompeakqualtools | `>=1.2,<2` |
| preseq | `>=3.1,<4` |
| seacr | `>=1.3,<2` |
| ucsc-bedgraphtobigwig | `>=377,<378` |

Also fixed: `openjdk >=17` → `openjdk >=17,<18` in picard.yml for consistency.

## 4. Files

| File | Change |
|------|--------|
| `workflow/envs/chipseq.yml` | 8 pins |
| `workflow/envs/align.yml` | 2 pins |
| `workflow/envs/samtools.yml` | 2 pins |
| `workflow/envs/deeptools.yml` | 2 pins |
| `workflow/envs/trim.yml` | 2 pins |
| `workflow/envs/seacr.yml` | 2 pins |
| `workflow/envs/phantompeakqualtools.yml` | 2 pins |
| `workflow/envs/fastqc.yml` | 1 pin |
| `workflow/envs/idr.yml` | 1 pin |
| `workflow/envs/macs3.yml` | 1 pin |
| `workflow/envs/multiqc.yml` | 1 pin |
| `workflow/envs/picard.yml` | 2 pins + openjdk cap |
| `workflow/envs/preseq.yml` | 1 pin |
| `workflow/envs/ucsc.yml` | 1 pin |
| `test/test_stage59_env_pinning.py` | Enforcement test |

## 5. Non-goals

- No rule, target, config, manifest, or Artifact changes
- No lockfiles
- No transitive dependency pinning
- No Python version changes
- No Co-Authored-By

## 6. Verification

```bash
python3 test/test_stage59_env_pinning.py
git diff --check
```
