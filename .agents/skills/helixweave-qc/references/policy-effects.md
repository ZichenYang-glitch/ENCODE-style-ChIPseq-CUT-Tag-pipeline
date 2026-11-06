# Policy effects on QC numbers

Read the per-sample `analysis_policy` JSON. It records the resolved
duplicate policy, PE insert bound, organelle mode, status, and matched
contigs, blacklist mode and application, the scientific input stage for
each QC product, stage counts, and the QC schema versions
(`docs/output-contract.md`, "QC Contract (v2)"). It has no absolute
workspace paths or environment values. Do not recompute those fields from
assay name, peak mode, or a checklist.

JSON `schema_version` `1` and `2` are both readable. Version 2 adds
`analysis_template_status` and `input_stages.analysis_template_metrics`.
Keys are `ANALYSIS_POLICY_V1_KEYS` and `ANALYSIS_POLICY_KEYS` in
`src/encode_pipeline/qc/contract.py`. `parse_analysis_policy_bytes` in
`src/encode_pipeline/qc/policy.py` rejects any other shape as
`AnalysisPolicyError`. Nested `qc_schema_versions` is checked as the string
`"2"` for sample, mnase, fragment, and mito. That field does not replace
`schema_version` on each TSV.

The sample v2 table is a view of the same policy plus metric cells. Its
`duplicate_analysis` column is `resolved_duplicate_analysis`, not the
config selector `auto` (`scripts/assemble_qc_summary.py`). One sample row
mixes stages: `analysis_input_stage` is the analysis-template stage, while
organelle cells follow the mito table. Use `input_stages` per product.

## Duplicate

`POLICY_DUPLICATE_POLICIES` is `exclude` and `retain`. The run selector is
`duplicate_analysis`: `auto`, `exclude`, or `retain`. The legacy alias
`remove_dup` maps `yes` to `exclude` and `no` to `retain`; conflicting
explicit values are rejected (`docs/configuration.md`, "Duplicate policy").

Read the selector from the validated config. On a direct run, use
`.helixweave/effective_config.json` (`config.duplicate_analysis` and the
per-sample `resolved_duplicate_analysis` map). Platform artifact runtime
metadata copies the config selector (`auto`, `exclude`, or `retain`) from
`_artifact_runtime_metadata` in `src/encode_pipeline/adapters/encode.py`.
That field is not the per-sample resolution. The sample policy and the
sample QC column hold the resolved value only.

Read `profile` on the policy (`legacy` or `updated`) before describing the
BAM. Legacy exclusion uses Picard `REMOVE_DUPLICATES=true` when Picard is
available, otherwise the `samtools markdup` fallback; a retained legacy
sample keeps a symlinked `final.bam`. Updated runs mark duplicates in
`final.bam`, and `analysis.bam` applies the policy (`-F 0x400` for
exclude). Every MACS3 call then renders `--keep-dup all`, so the caller
does not collapse coordinates again. The profile split is
`docs/configuration.md` ("Duplicate policy").

`cuttag_fragment_size` counts the canonical analysis BAM and does not itself
remove `0x400` records (`docs/output-contract.md`). NRF/PBC stays on the
pre-duplicate MAPQ BAM (`input_stages.nrf_pbc`).

## Organelle

`POLICY_ORGANELLE_MODES` is `auto`, `exclude`, `report`, and `keep`. The
value that affected the numbers is `organelle_mode` on the sample policy
(and `resolved_organelle_mode` on `effective_config.json`). Do not apply
`assay_organelle_mode` over that recorded mode.

`mito_fraction` uses the pre-duplicate MAPQ BAM. Its `input_stage` is
`pre_duplicate_mapq_bam`. Effects of the recorded mode are the invariants in
`scripts/calc_organelle_metrics.py`:

- `exclude` with no matched candidate contig fails closed before analysis
  outputs
- `report` with no match records `unresolved_contigs`, fraction `NA`, and
  organelle reads `NA`
- matched contigs and zero mapped reads record `no_mapped_reads`, fraction
  `NA`, and organelle reads `NA`
- matched contigs, mapped reads greater than zero, and zero organelle reads
  record `ok` and a real zero fraction
- `keep` does not enforce exclusion and does not fail on missing contigs

`policy.py` requires `stage_counts.organelle_reads` to be null when
`organelle_status` is not `ok`. A null count is `NA`. Do not turn
`unresolved_contigs` or `no_mapped_reads` into a zero fraction.

## Blacklist and analysis namespace

`POLICY_BLACKLIST_MODES` is `off`, `copies`, and `analysis`. Read
`blacklist_mode` and `blacklist_applied` from the policy JSON. The decision
table is `docs/configuration.md` ("Effective blacklist mode decision table").

`analysis` filters the canonical analysis BAM and final peaks. Those
analysis products are written under `<outdir>/analysis_filtered/`. Raw
`final.bam` stays under the configured outdir in every mode. Cite the
`output_type` and `path` on the row you read.
`src/encode_pipeline/manifest/make.py` names several of those rows
`analysis_filtered_*` when the analysis root differs from the raw outdir.
The output-contract type table lists the same names. A later sentence in
"Analysis Namespace Boundary" says the type stays canonical; do not rename
the stored type to match that sentence. Platform `analysis_branch` is
metadata on the artifact, not a substitute for the row.

`copies` leaves the canonical analysis BAM unfiltered by blacklist and adds
sidecar filtered BAM and peak copies. `off` does not blacklist-filter.

Read `raw_peak_count` and `delivered_peak_count` on the v2 sample summary.
Do not recompute the delivered count. The assembler uses the filtered peak
count as the delivered count when the recorded policy applied a blacklist.

## Input stage

Allowed stage tokens are `_ALLOWED_INPUT_STAGES` in `policy.py`:
`analysis_bam`, `pre_duplicate_markdup_bam`, `pre_duplicate_mapq_bam`,
`final_bam_to_analysis_bam`, and `not_applicable`. The recorded
`input_stages` object is the denominator map. Stage-count keys are
`POLICY_STAGE_COUNT_KEYS`. JSON null is `NA`, not zero.

`resolved_max_insert` is the PE bound the run used, or the string
`not_applicable` for SE (`POLICY_MAX_INSERT_NOT_APPLICABLE`). It is the
Bowtie2 `-X` alignment parameter, not a post-alignment TLEN filter
(`docs/assay-policy.md`, "Paired-end alignment insert limit"). Do not
replace the recorded integer with the assay table.

`analysis_template_status` explains null template counts. `ok` and
`no_valid_templates` are counts from `analysis_bam`. `not_applicable` and
`not_collected` use stage `not_applicable` and null counts. `parse_failed`
also has null counts. Version 1 policy has no `analysis_template_status`;
do not invent one.

## Candidate CUT&Tag default

`duplicate_analysis: auto` resolves per sample by assay. CUT&Tag retains
duplicates; ChIP-seq, ATAC-seq, and MNase-seq exclude them. Peak mode does
not make that choice. An explicit `exclude` or `retain` (including the
`remove_dup` alias) overrides the assay matrix for the run. The matrix and
the candidate statement are `docs/assay-policy.md` ("CUT&Tag", "Duplicate
policy") and `docs/configuration.md` ("Duplicate policy").

Say the retain default is still a candidate pending the SCI-7 real-library
data freeze only when the recorded selector is `auto`, the sample assay is
CUT&Tag, and the resolved value is `retain`. An explicit `retain` is that
configured choice. If the selector is not in the record you read, say the
selector is unknown. Do not rewrite a historical inconclusive result as a
frozen default.
