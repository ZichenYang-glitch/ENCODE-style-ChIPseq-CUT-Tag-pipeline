# Stage 10b: User Documentation Polish — Design Spec

**Date:** 2026-05-19
**Status:** design review
**Scope:** documentation only — three new user docs under `docs/`, conservative README trim with cross-links
**Excluded:** workflow/script/schema/config/test/env changes, large README rewrites, release tagging

---

## 1. Goals

1. Create three focused user-facing docs under `docs/`: quickstart, sample sheet reference, and configuration reference.
2. Conservatively trim README by moving full reference tables and long explanations into the new docs, keeping README self-sufficient for a first run.
3. Add short cross-links from README sections to the corresponding detailed docs.
4. No workflow, rule, script, schema, config, test, or env changes.

---

## 2. Files

| File | Action | Purpose |
|------|--------|---------|
| `docs/quickstart.md` | Create | Extended install/run guide + smoke tests + troubleshooting |
| `docs/sample-sheet.md` | Create | Full column reference, role semantics, examples, pitfalls |
| `docs/configuration.md` | Create | Complete configuration key reference with defaults and dependencies |
| `README.md` | Edit | Conservative trim; add cross-links; keep self-sufficient for first run |

No other files created or modified.

---

## 3. README changes

### 3.1 Quick Start section

**Remove:** explanatory paragraphs between the numbered steps (~12 lines of prose that restate what the commands do).

**Keep:**
- Step 1: clone + conda env create + activate
- Step 2: "Edit `config/samples.tsv` (see Sample Sheet below)."
- Step 3: core config YAML block (threads, mapq, binsize, remove_dup, trim, extend_reads, use_control, multiqc)
- Step 4: validate + dry-run commands
- Step 5: dry-run + run + resume commands

**Add** after step 5:

```markdown
> For troubleshooting, smoke-test instructions, and detailed run guidance, see
> [docs/quickstart.md](docs/quickstart.md).
```

### 3.2 Sample Sheet section

**Remove:**
- Large optional-columns table for replicate metadata (columns `experiment`, `condition`, `replicate`, `biological_replicate`, `technical_replicate`) — move to `docs/sample-sheet.md`.
- Large optional control-columns table (`role`, `control_sample`, `control_bam`) — move to `docs/sample-sheet.md`.
- Technical-replicate merging paragraph — move to `docs/sample-sheet.md`.
- Minimal-sample-sheet backward-compatibility sentence — move to `docs/sample-sheet.md`.

**Keep:**
- Required-columns table (unchanged).
- Compact PE two-replicate example table (unchanged).
- One sentence noting that optional replicate/control columns are documented in the schema and the sample-sheet doc.

**Add** after the example:

```markdown
> Full column reference, role/control semantics, additional examples
> (treatment+control, biological replicates), and common pitfalls:
> [docs/sample-sheet.md](docs/sample-sheet.md)
```

### 3.3 Configuration section

**Remove:**
- QC switches YAML block and explanation (~12 lines) — move to `docs/configuration.md`.
- Replicate and IDR features YAML block and explanation (~10 lines) — move to `docs/configuration.md`.
- Tool parameters YAML block and explanation (~10 lines) — move to `docs/configuration.md`.

**Keep:**
- Core config keys YAML block (already shown in Quick Start; keep as reference).
- Genome resources YAML block and explanation (unchanged).

**Add** after the genome resources block:

```markdown
> Full configuration reference covering core keys, genome resources, QC
> switches, replicates/IDR, CUT&Tag SEACR, tool parameters, `use_control`,
> `multiqc`, and key dependencies:
> [docs/configuration.md](docs/configuration.md)
```

### 3.4 Workflow Capabilities / Output Structure

Mostly unchanged. Add short cross-links where useful:

- After the Controls subsection: link to `docs/sample-sheet.md` for role/control semantics.
- After the Replicates subsection: link to `docs/configuration.md` for `stage4b`/`stage5` gating.

### 3.5 Developer Notes — Smoke-test profiles subsection

Add one sentence pointing to `docs/quickstart.md` for troubleshooting and smoke-test details.

### 3.6 What stays untouched

- H1 title and badge row.
- Overview + beta status sentence.
- Key Features section.
- Quick Start commands.
- Required-columns table + PE example.
- Core config YAML + genome resources.
- All Workflow Capabilities subsections.
- Limitations section (full).
- Output Structure section.
- Developer Notes (except the smoke-test cross-link).
- License section.

### 3.7 Expected line-count impact

Target net reduction: roughly 40–80 lines. Final README roughly 430–480 lines (down from ~520). Line-count reduction is a guideline, not a hard acceptance criterion — clarity and self-sufficiency take priority.

---

## 4. `docs/quickstart.md`

### Outline

```
# Quick Start

## 1. Install
## 2. Configure samples
## 3. Adjust workflow options
## 4. Validate
## 5. Dry-run
## 6. Run the workflow
## Smoke and test execution
## Troubleshooting
```

### Section content

**1. Install** — clone, `conda env create -f workflow/envs/chipseq.yml`, `conda activate chipseq`. Mention micromamba as alternative.

**2. Configure samples** — edit `config/samples.tsv`, point to `docs/sample-sheet.md` for full reference.

**3. Adjust workflow options** — edit `config/config.yaml`, point to `docs/configuration.md` for full reference.

**4. Validate** — `python3 scripts/validate_samples.py --config config/config.yaml` and `snakemake -s workflow/Snakefile --configfile config/config.yaml -n`. Explain that validation exit code != 0 means fix before running.

**5. Dry-run** — `snakemake -s workflow/Snakefile --configfile config/config.yaml -n`. Explain what a dry-run checks (DAG resolution, rule connectivity, file path existence warnings). Add a note that `--use-conda` is optional and only needed if the user wants Snakemake to manage rule environments.

**6. Run the workflow** — `--cores N` vs `--jobs` distinction, `--rerun-incomplete` for resuming interrupted runs, `--latency-wait` for NFS filesystems.

**Smoke and test execution** — three bundled test commands:

```bash
# Validation stress tests (15 checks, <5 s)
python3 test/test_validation_stress.py

# Dry-run smoke profiles (7 profiles, <30 s, requires snakemake on PATH)
python3 test/test_stage8_smoke_profiles.py

# Tiny real execution (preprocessing + signal, <60 s, requires full chipseq env)
python3 test/test_stage8b_tiny_execution.py
```

Each listed with its purpose, expected duration, and prerequisites. Distinguish which commands need user data and which run on bundled synthetic fixtures.

**Troubleshooting** — common beginner issues:
- Slow micromamba/conda solve (channel priority, using libmamba solver)
- `snakemake: command not found` (conda activate not run, wrong environment)
- FASTQ path not found (typo in sample sheet, absolute vs relative paths, WSL path translation)
- Bowtie2 index prefix issues (pointing to directory instead of prefix, missing `.1.bt2` files)
- `control_sample` or `control_bam` path not found (file doesn't exist, use_control still false)
- Missing external tools for tiny execution (snakemake not found, bowtie2 not on PATH)

---

## 5. `docs/sample-sheet.md`

### Outline

```
# Sample Sheet Reference

## Required columns
## Optional columns
### Replicate metadata
### Control columns
## Role semantics
## Examples
### Minimal single-sample
### Treatment with control sample
### Biological replicates
## Common pitfalls
```

### Section content

**Required columns** — full table with column name, type, description, and constraints (e.g., `sample` must match `[A-Za-z0-9_.-]+`). Same information as current README table but with additional constraint detail.

**Optional columns — Replicate metadata** — table: `experiment`, `condition`, `replicate`, `biological_replicate`, `technical_replicate` with defaults and descriptions. Explanation of how technical replicates merge into biological-replicate BAMs.

**Optional columns — Control columns** — table: `role`, `control_sample`, `control_bam` with defaults and descriptions. Note that these are ignored when `use_control: false`.

**Role semantics** — explains:
- `role: treatment` (default): preprocessing + peak calling.
- `role: control`: preprocessing only, no peak calling.
- `control_sample` references another row's `sample` ID; that row must have `role: control`.
- `control_bam` is an external path; not processed by the pipeline.

**Examples:**
- *Minimal single-sample*: required columns only, SE ChIP-seq, no controls, no replicates.
- *Treatment with control sample*: two-row sheet, treatment row referencing a control row.
- *Biological replicates*: two treatment rows with same `experiment`, different `biological_replicate`.

**Common pitfalls:**
- PE sample missing `fastq_2`.
- Sample ID contains spaces or special characters.
- `control_sample` references a row with `role: treatment`.
- `control_sample` references itself or forms a cycle.
- `control_bam` path set while `use_control` is `false`; controls are disabled unless `use_control: true`.
- Forgetting to set `genome` to a key that exists in `genome_resources`.

---

## 6. `docs/configuration.md`

### Outline

```
# Configuration Reference

## Core keys
## Genome resources
## `use_control`
## `multiqc`
## QC block
## Replicate and IDR features
## CUT&Tag SEACR
## Tool parameters
## Defaults and key dependencies
```

### Section content

**Core keys** — `samples`, `outdir`, `threads`, `mapq`, `binsize`, `remove_dup`, `trim`, `extend_reads`. Each with default value, accepted values, and effect. `remove_dup` options explained: `auto` / `yes` / `no` with assay-dependent behavior.

**Genome resources** — `effective_genome_size` dual syntax (MACS3 shortcut string like `hs`/`mm` or positive integer). Optional `chrom_sizes`, `blacklist`, `gtf`, `reference_fasta`. Path validation behavior.

**`use_control`** — boolean, default `false`. When `true`, enables `control_sample` and `control_bam` resolution. When `false`, both are ignored.

**`multiqc`** — boolean, default `true`. When `true`, aggregates QC artifacts into `results/multiqc/multiqc_report.html`.

**QC block** — each switch (`blacklist_filter`, `frip`, `library_complexity`, `nrf_pbc`, `signal_tracks`, `summary`) with default and effect. Notes: `blacklist_filter` requires a blacklist BED in genome resources; `signal_tracks` produces bedGraph files (BigWig conversion not yet implemented).

**Replicate and IDR features** — `stage4b` (default `true`, replicate-aware pooled outputs), `stage5` (default `false`, requires `stage4b: true`). `idr` block keys: `seed`, `threshold`, `rank`. IDR gating: only runs when `stage5: true`, `chipseq` assay, `narrow` peak_mode, exactly 2 treatment biological replicates.

**CUT&Tag SEACR** — `cuttag.seacr.enabled` (boolean, default `false`). When `true`, produces SEACR peak calls alongside MACS3 peaks for CUT&Tag samples. Output under `results/<sample>/04_peaks_seacr/`.

**Tool parameters** — `macs3.qvalue`, `idr_macs3.pvalue`, `bamcoverage.normalize_using`, optional `extra_args` escape hatches. Each with default and effect.

**Defaults and key dependencies** — summary table:

| Key | Default | Takes effect when |
|-----|---------|-------------------|
| `use_control` | `false` | Always; enables control resolution |
| `multiqc` | `true` | Always; enables MultiQC aggregation |
| `qc.*` | all `true` | Always; individual switches gate each metric |
| `stage4b` | `true` | Always; enables pooled outputs for multi-biorep experiments |
| `stage5` | `false` | `stage4b: true` + `chipseq` + `narrow` + exactly 2 bioreps |
| `idr.*` | listed | Only when `stage5: true` |
| `cuttag.seacr.enabled` | `false` | Always; affects CUT&Tag samples only |
| `tool_parameters.*` | tool defaults | Always; absent keys use built-in tool defaults |

---

## 7. Verification plan

1. **File scope:** `git diff --stat` shows only `README.md` + 3 new files under `docs/`. No workflow, script, schema, config, test, or env changes.
2. **No pipeline regressions:**
   - `python3 test/test_validation_stress.py` — 15/15 PASS.
   - `python3 test/test_stage8_smoke_profiles.py` — 7/7 PASS.
3. **README self-sufficiency check:** A new user reading only README can install, configure a basic sample sheet, validate, dry-run, and run without clicking any docs link.
4. **No duplicate/conflicting commands:** Verify that the same command appearing in README and a doc file uses identical syntax. Example: README and quickstart.md both show `snakemake -s workflow/Snakefile --configfile config/config.yaml -n` — wording must match.
5. **Command labeling:** Commands that require user data (`/data/ac1_R1.fq.gz`) are clearly distinguishable from bundled smoke-test commands that run on synthetic fixtures.
6. **Cross-link validity:** All three `docs/*.md` links in README point to existing files. Anchor links like `#troubleshooting` resolve correctly.
7. **Line-count check:** README line count before and after; confirm reduction is moderate (guideline, not gate).
8. **Manual read-through:** Each doc file read independently — coherent, no dangling references, no TBDs.

---

## 8. Out-of-scope

- No workflow, rule, script, schema, config, test, or env changes.
- Do not add large new README sections; short cross-links or a compact documentation link block are acceptable.
- No `docs/running.md` or fourth user doc file.
- No changes under `docs/superpowers/`.
- No release tagging.
- No new features.

---

## 9. Self-review

- **No README bloat:** Three links added, ~40–80 lines of reference material moved out. README remains self-sufficient.
- **No duplicated or conflicting commands:** Each command appears in full in one authoritative location; cross-links avoid re-stating commands at length in two places. Quick start commands stay in README; docs provide expansion and troubleshooting, not divergent copies.
- **No workflow/config/schema/script/test/env changes:** All changes are under `docs/*.md` and `README.md` only.
- **Doc outlines are concrete:** Every section has a described purpose. No TBDs or placeholders.
- **Cross-links are bidirectional in spirit:** README links to docs; each doc's intro references the README for context.
