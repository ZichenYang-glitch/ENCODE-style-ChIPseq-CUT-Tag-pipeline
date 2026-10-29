---
name: helixweave-qc
description: >-
  Explain an existing HelixWeave ENCODE-style ChIP-seq, CUT&Tag, ATAC-seq,
  or MNase-seq run: which products exist, which stage and denominator a
  metric uses, and how duplicate, organelle, and blacklist policy affect
  the numbers. Use when the user asks what finished, how to read QC, or
  how to trace an artifact. Do not use to start a workflow, download data,
  or troubleshoot a crash.
---

# HelixWeave QC

Explain products and QC for a run that already exists. Command inventory and
the short status distinctions are in
[../_shared/entries.md](../_shared/entries.md). Leave that inventory there.
If it and the checkout disagree, follow the checkout.

A crash, validator failure, or source that will not parse belongs to
[helixweave-troubleshoot](../helixweave-troubleshoot/SKILL.md). This skill
explains values that can be read.

## Read-only

Answering a QC question does not start, resume, cancel, or dry-run a
workflow. Do not download sequencing data.

A platform artifact download waits until the user asks for that download.
Use only `GET /api/v1/runs/{run_id}/artifacts/{artifact_id}/download` for an
artifact on the generation currently listed. On
`RUN_ARTIFACT_GENERATION_CHANGED` or `RUN_QC_METRIC_GENERATION_CHANGED`,
reload from the first page before naming an artifact or a metric.

## Workflow id

ENCODE sample columns, assays, paths, and the tables below apply only to
workflow id `encode-style-chipseq-cuttag-atac-mnase`
(`ENCODE_WORKFLOW_ID` in `src/encode_pipeline/platform/registry.py`). For
any other id, including `bulk-rnaseq` and `hitrac-preprocess`, read that
adapter's QC extractor and schema. Do not reuse ENCODE columns.

## Where to read

Platform runs: list artifacts and QC with the routes in
[../_shared/entries.md](../_shared/entries.md), including
`GET /api/v1/runs/{run_id}/artifacts` and
`GET /api/v1/runs/{run_id}/qc-metrics` and their generation cursors.
SQLite is lifecycle truth. Run status `succeeded` means execution finished.

Direct runs: read files under the configured outdir. `workflow/Snakefile`
sets `RAW_OUTDIR` to that directory and sets `OUTDIR` to
`<outdir>/analysis_filtered` only when the analysis blacklist branch is
active. The manifest is `<OUTDIR>/multiqc/result_manifest.tsv` (`status` is
`present`, `missing`, or `not_applicable`; `qc_flag` is always `NA`). For a
treatment sample, follow that manifest to `analysis_policy`, `qc_summary`,
and `mito_fraction`. Do not assume the outdir is named `results`.

Raw `final.bam` stays at
`<RAW_OUTDIR>/<sample>/02_align/<sample>.final.bam`. Cite the stored
`output_type` and `path`. Some filtered rows use an `analysis_filtered_*`
type (`src/encode_pipeline/manifest/make.py`); the output-contract table
lists those types, while its "Analysis Namespace Boundary" sentence says
the type stays canonical. The row wins. `analysis_branch` on a platform
artifact is separate metadata.

A workspace TSV does not prove platform listing, download, or complete
MultiQC HTML visibility (`docs/configuration.md`).

## Schema before columns

Read `schema_version` before using any column. TSV versions are the string
constants in `src/encode_pipeline/qc/contract.py` (current writers use
`"2"`). Policy JSON `schema_version` is the integer `1` or `2`
(`ANALYSIS_POLICY_VERSION` in `src/encode_pipeline/qc/policy.py`).

Version 1 stays readable and is not rewritten during extraction. Version 2
puts `schema_version` first, preserves the version 1 column order, and
appends analysis-stage fields. Headers are the `*_V1_COLUMNS` and
`*_V2_COLUMNS` tuples in `contract.py`.

`src/encode_pipeline/adapters/encode_qc.py` (`EncodeQcSummaryError`) accepts
an exact v1 header or an exact v2 header. Unknown versions, extra or
duplicate columns, loose numbers, and undefined statuses fail the source.
Policy documents fail as `AnalysisPolicyError`. That failure is not a zero.

## Status and policy

Allowed sets, what `NA` means, and manifest statuses are in
[references/status-values.md](references/status-values.md). A numeric `0` is
a measurement only when the status for that column allows a defined value.

How duplicate, organelle, blacklist, namespace, and input stage affect the
numbers is in [references/policy-effects.md](references/policy-effects.md).
Read the sample `analysis_policy` JSON. Do not recompute it.

## What an answer can claim

Name the file or generation, the schema version, the status beside each
number, and the recorded input stage. Then say what is still unknown:

- an optional QC switch left off, or a product the gating table does not
  schedule (`docs/qc-interpretation.md`, "Availability and defaults";
  `docs/output-contract.md`, "Gating Conditions")
- whether MultiQC ingested that module; a workspace TSV can exist without
  HTML visibility
- a producer that did not run, or a manifest row that is `missing`

Checklist numbers in `docs/qc-interpretation.md` ("Practical Checklists")
are guidelines. "Known Caveats" says the pipeline does not abort on low QC
metrics. Execution success, a green CI result, or a larger peak count is
not biological acceptance. Manifest `qc_flag` does not implement threshold
conclusions.

When the user asks about an IDR summary, read the rescue and
self-consistency rules in `docs/qc-interpretation.md` and
`docs/idr-contract.md`. Do not apply those rules to other metrics.

## Answer shape

1. Workflow id. For a platform run, add run id and the artifact or QC generation.
2. Products that are present, missing, or not applicable.
3. Each cited metric: schema version, status, recorded stage and denominator, and the value only when that status defines it.
4. Policy fields actually read, including the CUT&Tag `auto` retain candidate when the recorded selector is `auto`.
5. What those records still do not show.
